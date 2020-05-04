#!/usr/bin/env python

"""Tests for `streaming_pileupy` package."""

from click.testing import CliRunner

from streaming_pileupy import streaming_pileupy
from streaming_pileupy import cli


def test_command_line_interface():
    """Test the CLI."""
    runner = CliRunner()
    result = runner.invoke(cli.main)
    assert result.exit_code == 0
    assert "streaming_pileupy.cli.main" in result.output
    help_result = runner.invoke(cli.main, ["--help"])
    assert help_result.exit_code == 0
    assert "--help  Show this message and exit." in help_result.output
