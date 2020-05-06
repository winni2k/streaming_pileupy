"""streaming_pileupy"""
__author__ = """Warren W. Kretzschmar"""
__email__ = "winni@warrenwk.com"
__version__ = "0.5.2"

import logging
import sys
from collections import defaultdict
from contextlib import closing
from dataclasses import dataclass
from typing import IO, DefaultDict, Dict, List, Tuple

import click
import pysam

# def get_regions(bed_file):
#     with open(bed_file, 'r') as fh:
#         for line in bed_file:
#             yield line.rstrip('\n').split(maxsplit=3)[0:3]

logger = logging.getLogger(__name__)


def set_logging(verbose):
    if verbose == 1:
        logging.basicConfig(level=logging.INFO)
    elif verbose > 1:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARNING)


def log_record_info(match_base, rec, sequence):
    logger.debug(f"match_base: {match_base}")
    ref_positions = rec.get_reference_positions()
    ref_sequence = rec.get_reference_sequence()
    logger.debug(str(rec))
    logger.debug(f"ref_positions: len={len(ref_positions)}, positions={ref_positions}")
    logger.debug(f"ref_sequence: len={len(ref_sequence)}, seq={ref_sequence}")
    logger.debug(f"sequence: len={len(sequence)}, seq={sequence}")


def get_rg_lookup_table(header):
    return {rg["ID"]: rg["SM"] for rg in header["RG"]}


@dataclass
class MpileupWriter:
    out_fh: IO[str]
    samples: List[str]
    buffer: DefaultDict[int, DefaultDict[str, List[Tuple[str, int]]]]
    ref_bases: Dict[int, str]
    _chrom: str
    min_bq: int

    @classmethod
    def from_filehandle_and_samples(cls, fh, samples):
        return cls(
            out_fh=fh,
            samples=samples,
            buffer=defaultdict(lambda: defaultdict(list)),
            ref_bases={},
            _chrom="",
            min_bq=0,
        )

    def add_base(self, sample: str, pos: int, base: str, qual: int, ref: str) -> None:
        self.buffer[pos][sample].append((base, qual))
        self.ref_bases[pos] = ref
        if qual < self.min_bq:
            self.buffer[pos][sample].pop()

    @property
    def chrom(self):
        return self._chrom

    @chrom.setter
    def chrom(self, chrom):
        if self._chrom != chrom:
            self.flush_all()
            self._chrom = chrom

    def flush_to_position(self, pos):
        """Print the stored pileups for all records up to pos"""
        for flush_pos in sorted(p for p in self.buffer.keys() if p < pos):
            self._flush_position(flush_pos)

    def flush_all(self):
        """Print the stored pileups for all records"""
        for flush_pos in sorted(self.buffer.keys()):
            self._flush_position(flush_pos)

    def _flush_position(self, pos):
        ref_base = self.ref_bases.pop(pos)
        pos_bases = self.buffer.pop(pos)
        outline = [f"{self.chrom}\t{pos + 1}\t{ref_base.upper()}"]
        outline += list(self._get_sample_strings(pos_bases))
        outline.append("\n")
        self.out_fh.write("".join(outline))

    def _get_sample_strings(self, pos_bases):
        for sample in self.samples:
            sample_data = pos_bases.get(sample, None)
            if not sample_data:
                yield "\t0\t*\t*"
            else:
                bases, quals = list(zip(*sample_data))
                quals = [chr(q + 33) for q in quals]
                assert len(bases) == len(quals)
                n_bases = len(bases)
                bases = "".join(bases)
                quals = "".join(quals)
                yield f"\t{n_bases}\t{bases}\t{quals}"


@dataclass
class RecordReader:
    current_pos: int
    rg_table: Dict[str, str]
    writer: MpileupWriter

    def move_pos_to(self, pos):
        assert pos >= self.current_pos
        self.writer.flush_to_position(pos)
        self.current_pos = pos

    def ingest(self, contig, rec):
        """Ingest a pysam record"""
        self.writer.chrom = contig
        self.move_pos_to(rec.pos)
        rg_id = rec.get_tag("RG")
        if rec.is_unmapped:
            return
        if rec.is_reverse:
            match_base = ","
        else:
            match_base = "."

        sequence = rec.query_sequence
        qualities = rec.query_qualities
        if logger.isEnabledFor(logging.DEBUG):
            log_record_info(match_base, rec, sequence)

        for read_pos, ref_pos, ref_base in rec.get_aligned_pairs(
            matches_only=True, with_seq=True
        ):
            logger.debug("read: %s, ref: %s, ref_base: %s", read_pos, ref_pos, ref_base)

            base = sequence[read_pos].upper()
            if base == ref_base.upper():
                base = match_base
            if rec.is_reverse:
                base = base.lower()

            self.writer.add_base(
                sample=self.rg_table[rg_id],
                pos=ref_pos,
                base=base,
                qual=qualities[read_pos],
                ref=ref_base,
            )

    def close(self):
        self.writer.flush_all()


@click.command()
@click.version_option()
@click.option("-v", "--verbose", count=True)
@click.option(
    "-Q",
    "--min-BQ",
    help="skip bases with baseQ/BAQ smaller than INT",
    type=int,
    default=0,
)
@click.argument("input", type=click.Path(allow_dash=False))
@click.argument(
    "sample_file", type=click.File(),
)
def main(input, sample_file, verbose, min_bq):
    """Create a read-group-aware pileup from a single file.

    INPUT: Stream of a SAM/BAM file (including header), '-' reads from stdin.
    SAMPLE_FILE: File containing each sample to pileup on a separate line.

    The read group's SM tag is used to infer samples.
    Pileup columns correspond to samples in SAMPLE_FILE in the order given in the file.
    """
    set_logging(verbose)
    infile = pysam.AlignmentFile(input, "r")
    samples = [s.rstrip("\n") for s in sample_file]
    rg_table = get_rg_lookup_table(infile.header)

    writer = MpileupWriter.from_filehandle_and_samples(sys.stdout, samples)
    writer.min_bq = min_bq
    piler = RecordReader(current_pos=0, rg_table=rg_table, writer=writer,)
    with closing(piler):
        for rec in infile:
            piler.ingest(infile.get_reference_name(rec.reference_id), rec)

    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover