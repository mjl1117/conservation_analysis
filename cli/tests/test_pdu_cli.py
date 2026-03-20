# cli/tests/test_pdu_cli.py
from click.testing import CliRunner
from unittest.mock import patch, MagicMock
import json
import pytest

# Import pdu from cli/pdu.py — sys.path must include cli/
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from pdu import cli

runner = CliRunner()

def test_cli_help():
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "fetch" in result.output
    assert "align" in result.output
    assert "metrics" in result.output
    assert "compare" in result.output
    assert "colorpdb" in result.output

def test_cli_fetch_calls_fetch_homologs(tmp_path):
    with patch("pdu.fetch_reference_sequence", return_value=("header", "MKKLL")) as mock_ref, \
         patch("pdu.fetch_homologs", return_value=[("h2", "MRKLL")]) as mock_hom, \
         patch("pdu.write_fasta") as mock_write:
        result = runner.invoke(cli, ["fetch", "P0A1C7", "--cluster", "50",
                                     "--out", str(tmp_path / "out.fasta")])
        assert result.exit_code == 0
        mock_hom.assert_called_once()

def test_cli_metrics_outputs_json(tmp_path):
    mock_fasta = tmp_path / "aligned.fasta"
    mock_fasta.write_text(">P0A1C7\nMKKLL\n>seq2\nMRKLL\n")
    out = tmp_path / "conservation.json"
    result = runner.invoke(cli, [
        "metrics", str(mock_fasta), "--ref", "P0A1C7", "--out", str(out)
    ])
    assert result.exit_code == 0, result.output
    assert out.exists()
    data = json.loads(out.read_text())
    assert data["schema_version"] == "1"

def test_cli_run_executes_full_pipeline(tmp_path):
    """Full pipeline smoke test with mocked network calls."""
    with patch("pdu.fetch_reference_sequence", return_value=("P0A1C7", "MKKLL")), \
         patch("pdu.fetch_homologs", return_value=[("h2", "MRKLL"), ("h3", "MKRLL")]), \
         patch("pdu.run_clustalo", return_value=str(tmp_path / "aligned.fasta")) as mock_align, \
         patch("pdu.check_clustalo_available", return_value=True):
        # Write a mock aligned FASTA that run_clustalo would produce
        aligned = tmp_path / "aligned.fasta"
        aligned.write_text(">P0A1C7\nMKKLL\n>h2\nMRKLL\n>h3\nMKRLL\n")
        mock_align.return_value = str(aligned)
        result = runner.invoke(cli, [
            "run", "P0A1C7", "--cluster", "50", "--out", str(tmp_path)
        ])
        assert result.exit_code == 0, result.output
        assert (tmp_path / "conservation.json").exists()
