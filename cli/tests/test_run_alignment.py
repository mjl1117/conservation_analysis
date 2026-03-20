import subprocess
from unittest.mock import patch, MagicMock
import pytest
from scripts.run_alignment import run_clustalo, check_clustalo_available

def test_check_clustalo_available_returns_true_when_found():
    with patch("shutil.which", return_value="/usr/local/bin/clustalo"):
        assert check_clustalo_available() is True

def test_check_clustalo_available_returns_false_when_missing():
    with patch("shutil.which", return_value=None):
        assert check_clustalo_available() is False

def test_run_clustalo_calls_correct_command(tmp_path):
    input_fasta = tmp_path / "input.fasta"
    output_fasta = tmp_path / "output.fasta"
    input_fasta.write_text(">seq1\nMKKLL\n")
    expected_output = ">seq1\nMKKLL\n"

    with patch("scripts.run_alignment.shutil.which", return_value="/usr/local/bin/clustalo"):
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)
            output_fasta.write_text(expected_output)  # simulate clustalo writing output
            result = run_clustalo(str(input_fasta), str(output_fasta))

    call_args = mock_run.call_args[0][0]
    assert "clustalo" in call_args[0]
    assert str(input_fasta) in call_args
    assert str(output_fasta) in call_args
    assert "--force" in call_args

def test_run_clustalo_raises_on_nonzero_exit(tmp_path):
    input_fasta = tmp_path / "input.fasta"
    input_fasta.write_text(">seq1\nMKKLL\n")
    with patch("scripts.run_alignment.shutil.which", return_value="/usr/local/bin/clustalo"):
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=1, stderr="error")
            with pytest.raises(RuntimeError, match="ClustalOmega failed"):
                run_clustalo(str(input_fasta), str(tmp_path / "out.fasta"))

def test_run_clustalo_raises_if_not_installed(tmp_path):
    input_fasta = tmp_path / "input.fasta"
    input_fasta.write_text(">seq1\nMKKLL\n")
    with patch("scripts.run_alignment.shutil.which", return_value=None):
        with pytest.raises(RuntimeError, match="clustalo not found"):
            run_clustalo(str(input_fasta), str(tmp_path / "out.fasta"))
