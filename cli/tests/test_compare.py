import json
import pytest
from scripts.compare import (
    compute_sequence_comparison,
    compute_phylogenetic_comparison,
    build_comparison_json,
    pairwise_align,
)

TINY_ALIGNED_FOR_COMPARE = (
    ">P0A1C7\nMKKLLVIGGSGSIGKSTLL\n"
    ">Q9X1X2\nMKKLLVIGGAGSIGKSTLL\n"
    ">Q8ZRE4\nMKKLLVIGGAGSIGKSTLL\n"
    ">P12345\nMRRLLVIGGAGSIGKATLL\n"
)

ALIGNED_A = "MKKLLVIGGSGSIGKSTLL"
ALIGNED_B = "MKKLLVIGGSGSIGKSTLL"
ALIGNED_B_DIFF = "MKKLLVIGGSGSIGKATLL"

def test_pairwise_align_returns_equal_length_strings():
    seq_a = "MKKLLVIGGSGSIG"
    seq_b = "MKKLLVIGGSGSIGKSTLL"  # longer
    aligned_a, aligned_b = pairwise_align(seq_a, seq_b)
    assert len(aligned_a) == len(aligned_b)
    # Gaps inserted in shorter sequence
    assert "-" in aligned_a or "-" in aligned_b

def test_pairwise_align_identical_sequences_no_gaps():
    seq = "MKKLLVIGGSGSIG"
    aligned_a, aligned_b = pairwise_align(seq, seq)
    assert aligned_a == seq
    assert aligned_b == seq
    assert "-" not in aligned_a

def test_sequence_comparison_identical():
    result = compute_sequence_comparison(ALIGNED_A, ALIGNED_A)
    assert result["pairwise_identity"] == pytest.approx(1.0)

def test_sequence_comparison_one_diff():
    result = compute_sequence_comparison(ALIGNED_A, ALIGNED_B_DIFF)
    n = len(ALIGNED_A)
    expected_identity = (n - 1) / n
    assert result["pairwise_identity"] == pytest.approx(expected_identity)

def test_sequence_comparison_has_aligned_strings():
    result = compute_sequence_comparison(ALIGNED_A, ALIGNED_B_DIFF)
    assert result["aligned_a"] == ALIGNED_A
    assert result["aligned_b"] == ALIGNED_B_DIFF

def test_sequence_comparison_raises_on_unequal_length():
    with pytest.raises(ValueError, match="equal length"):
        compute_sequence_comparison("MKK", "MKKLL")

def test_build_comparison_json_valid_schema():
    import jsonschema
    from pathlib import Path
    result = build_comparison_json(
        protein_a="P0A1C7",
        protein_b="Q9X1X2",
        cluster="UniRef50",
        sequence=compute_sequence_comparison(ALIGNED_A, ALIGNED_B),
        phylogenetic={"tree_distance": 0.05, "lca_node": "abc12345",
                      "shared_cluster_fraction": 0.8},
        structural={"available": False},
    )
    schema_path = Path(__file__).parent.parent.parent / "schema/v1/comparison.schema.json"
    schema = json.loads(schema_path.read_text())
    jsonschema.validate(result, schema)

def test_compute_phylogenetic_comparison_returns_valid_output():
    from scripts.phylogenetic_tree import build_nj_tree
    newick = build_nj_tree(TINY_ALIGNED_FOR_COMPARE)
    result = compute_phylogenetic_comparison(newick, "P0A1C7", "P12345")
    assert "tree_distance" in result
    assert result["tree_distance"] >= 0
    assert len(result["lca_node"]) == 8  # SHA-1 first 8 hex chars
    assert 0.0 <= result["shared_cluster_fraction"] <= 1.0

def test_compute_phylogenetic_comparison_same_protein():
    from scripts.phylogenetic_tree import build_nj_tree
    newick = build_nj_tree(TINY_ALIGNED_FOR_COMPARE)
    result = compute_phylogenetic_comparison(newick, "P0A1C7", "P0A1C7")
    assert result["tree_distance"] == pytest.approx(0.0, abs=1e-6)
