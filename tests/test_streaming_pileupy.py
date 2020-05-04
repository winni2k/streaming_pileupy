#!/usr/bin/env python

"""Tests for `streaming_pileupy` package."""
import re
from dataclasses import dataclass
from pathlib import Path
from subprocess import PIPE, run
from typing import List

from click.testing import CliRunner

from streaming_pileupy import main


def test_command_line_interface():
    """Test the CLI."""
    runner = CliRunner()
    result = runner.invoke(main)
    assert result.exit_code != 0
    assert "Usage:" in result.output
    help_result = runner.invoke(main, ["--help"])
    assert help_result.exit_code == 0
    assert "--help  Show this message and exit." in help_result.output


@dataclass
class TestFileBuilder:
    tmpdir: Path
    sam: str
    samples: List

    def with_sam(self, sam, samples):
        self.sam = re.sub(r" +", "\t", sam)
        self.samples = samples

    def build(self):
        input_sam = self.tmpdir / "input.sam"
        input_samples = self.tmpdir / "input.samples"
        with open(input_sam, "w") as fh:
            fh.write(self.sam)
        with open(input_samples, "w") as fh:
            fh.write("\n".join(self.samples) + "\n")
        return input_sam, input_samples


def test_two_record_sam(tmpdir):
    # given
    builder = TestFileBuilder(tmpdir, "", [])
    builder.with_sam(
        "@HD VN:1.6  SO:unknown\n"
        "@SQ SN:1 LN:249250621\n"
        "@RG ID:0    SM:sample_0\n"
        "@RG ID:1    SM:sample_1\n"
        "r0  0   1    24  0   1M  *   0   0   G   I   RG:Z:0	NM:i:1	MD:Z:0N0\n"
        "r1  0   1    24  0   1M  *   0   0   G   I   RG:Z:1	NM:i:1	MD:Z:0N0\n",
        ["sample_0", "sample_1"],
    )

    # when
    sam, samples = builder.build()
    result = run(f"spileup {sam} {samples}", shell=True, stdout=PIPE, text=True)

    # then
    expected = "1 24 N 1 G I 1 G I".replace(" ", "\t")
    assert expected == result.stdout
