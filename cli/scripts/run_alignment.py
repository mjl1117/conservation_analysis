"""Run ClustalOmega alignment on a FASTA file."""
from __future__ import annotations
import shutil
import subprocess
from pathlib import Path


def check_clustalo_available() -> bool:
    return shutil.which("clustalo") is not None


def run_clustalo(
    input_fasta: str,
    output_fasta: str,
    threads: int = 4,
) -> str:
    """
    Run ClustalOmega on input_fasta, write aligned FASTA to output_fasta.
    Returns the output_fasta path. Raises RuntimeError on failure.
    """
    if not check_clustalo_available():
        raise RuntimeError(
            "clustalo not found. Install with: conda install -c bioconda clustalo"
        )
    cmd = [
        "clustalo",
        "-i", str(input_fasta),
        "-o", str(output_fasta),
        f"--threads={threads}",
        "--force",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"ClustalOmega failed: {result.stderr}")
    return str(output_fasta)


def read_aligned_fasta(fasta_path: str) -> str:
    """Return the contents of an aligned FASTA file as a string."""
    return Path(fasta_path).read_text()
