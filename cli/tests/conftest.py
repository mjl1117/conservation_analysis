"""Shared fixtures for all CLI tests."""
import json
from pathlib import Path
import pytest

FIXTURES_DIR = Path(__file__).parent / "fixtures"

@pytest.fixture
def tiny_fasta():
    """A minimal 3-sequence aligned FASTA for fast tests."""
    return (
        ">seq1\nMKKLLVIGGSGSIGKSTLLQRLAQKSGI\n"
        ">seq2\nMKKLLVIGGAGSIGKSTLLQRLAEKSGI\n"
        ">seq3\nMKKLLVIGGAGSIGKSTLLQRLAEKTGV\n"
    )

@pytest.fixture
def tiny_unaligned_fasta():
    return (
        ">P0A1C7\nMQQEALGMVETKGLTAAIEAADAMVKSANVMLVGYEKIGSGLVTVIVRGDVGAVKAATDAGAAAARNVGEVKAVHVIPRPHTDVEKILPKGISQ\n"
        ">Q9X1X2\nMQQEALGMVETKGLTAAIEAADAMVKSANVMLVGYEKIGSGLVTVIVRGDVGAVKAATDAGAAAARNVGEVKAVHVIPRPHTDVEKILP\n"
    )

@pytest.fixture
def sample_conservation_json():
    return {
        "schema_version": "1",
        "accession": "P0A1C7",
        "cluster": "UniRef50",
        "n_sequences": 10,
        "n_columns": 28,
        "gap_zscore_k": 1.5,
        "gap_threshold_zscore": 1.5,
        "gap_threshold_fraction": 0.3,
        "residues": [
            {
                "position": 1,
                "amino_acid": "M",
                "shannon_entropy": 0.0,
                "rel_entropy": 4.32,
                "mutability": 0.0,
                "conservation": 1.0,
                "gap_fraction": 0.0,
                "gap_zscore": -1.0,
                "aa_types_count": 1,
                "aa_types_observed": ["M"],
                "mutability_class": "invariant",
            }
        ],
    }

@pytest.fixture
def minimal_pdb_lines():
    """Minimal PDB ATOM records for testing."""
    return [
        "ATOM      1  CA  MET A   1      10.000  10.000  10.000  1.00 50.00           C",
        "ATOM      2  CA  GLN A   2      11.000  11.000  11.000  1.00 50.00           C",
        "HETATM    3  O   HOH A 101      20.000  20.000  20.000  1.00 30.00           O",
    ]
