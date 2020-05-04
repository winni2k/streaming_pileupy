"""streaming_pileupy"""
__author__ = """Warren W. Kretzschmar"""
__email__ = "winni@warrenwk.com"
__version__ = "0.1.0"

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

QUERY_CONSUMERS = (
    pysam.CMATCH,
    pysam.CINS,
    pysam.CSOFT_CLIP,
    pysam.CEQUAL,
    pysam.CDIFF,
)


def get_rg_lookup_table(header):
    return {rg["ID"]: rg["SM"] for rg in header["RG"]}


@dataclass
class MpileupWriter:
    out_fh: IO[str]
    samples: List[str]
    buffer: DefaultDict[int, DefaultDict[str, List[Tuple[str, int]]]]
    ref_bases: Dict[int, str]
    _chrom: str

    @classmethod
    def from_filehandle_and_samples(cls, fh, samples):
        return cls(
            out_fh=fh,
            samples=samples,
            buffer=defaultdict(lambda: defaultdict(list)),
            ref_bases={},
            _chrom="",
        )

    def add_base(self, sample: str, pos: int, base: str, qual: int, ref: str) -> None:
        self.buffer[pos][sample].append((base, qual))
        self.ref_bases[pos] = ref

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
        self.out_fh.write(f"{self.chrom}\t{pos+1}\t{ref_base.upper()}")
        for sample in self.samples:
            try:
                bases, quals = list(zip(*pos_bases[sample]))
            except KeyError:
                n_bases, bases, quals = 0, "*", "*"
            else:
                quals = [chr(q + 33) for q in quals]
                assert len(bases) == len(quals)
                n_bases = len(bases)
                bases = "".join(bases)
                quals = "".join(quals)
            self.out_fh.write(f"\t{n_bases}\t{bases}\t{quals}")


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

        current_read_pos = 0
        ref_positions = rec.get_reference_positions()
        ref_sequence = rec.get_reference_sequence()
        sequence = rec.get_forward_sequence()
        qualities = rec.get_forward_qualities()
        for operation, cigar_len in rec.cigartuples:
            if operation in (pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF):
                for read_pos in range(current_read_pos, cigar_len + current_read_pos):
                    base = sequence[read_pos]
                    if base == ref_positions[read_pos]:
                        base = match_base
                    self.writer.add_base(
                        sample=self.rg_table[rg_id],
                        pos=ref_positions[read_pos],
                        base=base,
                        qual=qualities[read_pos],
                        ref=ref_sequence[read_pos],
                    )
            if operation in QUERY_CONSUMERS:
                current_read_pos += cigar_len

    def close(self):
        self.writer.flush_all()


@click.command()
@click.argument("input", type=click.Path(allow_dash=False))
@click.argument(
    "sample_file", type=click.File(),
)
def main(input, sample_file):
    """Create a read-group-aware pileup from a single file.

    INPUT: Stream of a SAM/BAM file (including header).
    SAMPLE_FILE: File containing each sample to pileup on a separate line.

    The read group's SM tag is used to infer samples.
    """
    infile = pysam.AlignmentFile(input, "r")
    samples = [s.rstrip("\n") for s in sample_file]
    rg_table = get_rg_lookup_table(infile.header)

    writer = RecordReader(
        current_pos=0,
        rg_table=rg_table,
        writer=MpileupWriter.from_filehandle_and_samples(sys.stdout, samples),
    )
    with closing(writer):
        for rec in infile:
            writer.ingest(infile.get_reference_name(rec.reference_id), rec)

    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
