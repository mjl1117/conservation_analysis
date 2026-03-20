# cli/tests/test_conservation_metrics.py
import math
import pytest
from scripts.conservation_metrics import (
    GAP_CHARS,
    shannon_entropy,
    rel_entropy,
    mutability,
    conservation_score,
    gap_fraction,
    compute_gap_zscores,
    mutability_class,
    compute_metrics,
)

# --- Unit tests for individual metrics ---

def test_shannon_entropy_all_same():
    col = ["M"] * 10
    assert shannon_entropy(col) == pytest.approx(0.0)

def test_shannon_entropy_two_equal_types():
    col = ["M"] * 5 + ["K"] * 5
    assert shannon_entropy(col) == pytest.approx(1.0)

def test_shannon_entropy_excludes_gap_chars():
    col = ["-"] * 5 + ["M"] * 5
    # gaps excluded → 5/5 M → entropy = 0
    assert shannon_entropy(col) == pytest.approx(0.0)

def test_shannon_entropy_all_gaps_returns_zero():
    assert shannon_entropy(["-", "-", "X", "."]) == pytest.approx(0.0)

def test_shannon_entropy_excludes_X_and_dot():
    # X and . are gap chars
    col = ["X"] * 5 + ["M"] * 5
    assert shannon_entropy(col) == pytest.approx(0.0)

def test_rel_entropy_uniform_returns_zero():
    # 20 equal AA types → uniform → KL divergence from uniform = 0
    aas = list("ACDEFGHIKLMNPQRSTVWY")
    col = aas * 5  # 100 residues, 20 types
    assert rel_entropy(col) == pytest.approx(0.0, abs=1e-6)

def test_mutability_is_normalized_entropy():
    col = ["M"] * 5 + ["K"] * 5
    h = shannon_entropy(col)
    assert mutability(h) == pytest.approx(h / math.log2(20))

def test_conservation_is_one_minus_mutability():
    col = ["M"] * 10
    h = shannon_entropy(col)
    m = mutability(h)
    assert conservation_score(m) == pytest.approx(1.0 - m)

def test_gap_fraction_counts_dash_X_dot():
    col = ["-", "X", ".", "M", "K"]
    assert gap_fraction(col) == pytest.approx(3 / 5)

def test_gap_fraction_no_gaps():
    col = ["M", "K", "L"]
    assert gap_fraction(col) == pytest.approx(0.0)

def test_compute_gap_zscores_known_values():
    # fracs=[0.0, 0.5, 1.0]: mean=0.5, population std=sqrt(1/6)≈0.4082
    # z[0] = (0.0 - 0.5) / 0.4082 ≈ -1.2247
    fracs = [0.0, 0.5, 1.0]
    zscores = compute_gap_zscores(fracs)
    assert zscores[0] == pytest.approx(-1.2247, rel=1e-3)
    assert zscores[0] < zscores[1] < zscores[2]

def test_compute_gap_zscores_all_same_returns_zeros():
    fracs = [0.3, 0.3, 0.3]
    zscores = compute_gap_zscores(fracs)
    assert all(z == pytest.approx(0.0) for z in zscores)

def test_mutability_class_invariant():
    assert mutability_class(0.0) == "invariant"
    assert mutability_class(0.20) == "invariant"

def test_mutability_class_conserved():
    assert mutability_class(0.21) == "conserved"
    assert mutability_class(0.45) == "conserved"

def test_mutability_class_variable():
    assert mutability_class(0.46) == "variable"
    assert mutability_class(0.70) == "variable"

def test_mutability_class_hypervariable():
    assert mutability_class(0.71) == "hypervariable"
    assert mutability_class(1.0) == "hypervariable"

# --- Integration test for compute_metrics ---

TINY_ALIGNED_FASTA = (
    ">P0A1C7\nMKK\n"
    ">seq2\nMKK\n"
    ">seq3\nMKK\n"
    ">seq4\nMK-\n"
    ">seq5\nMRL\n"
)

def test_compute_metrics_returns_valid_schema(sample_conservation_json):
    import jsonschema, json
    from pathlib import Path
    result = compute_metrics(TINY_ALIGNED_FASTA, ref_accession="P0A1C7", gap_zscore_k=1.5)
    schema_path = Path(__file__).parent.parent.parent / "schema/v1/conservation.schema.json"
    schema = json.loads(schema_path.read_text())
    jsonschema.validate(result, schema)  # should not raise

def test_compute_metrics_residue_count():
    result = compute_metrics(TINY_ALIGNED_FASTA, ref_accession="P0A1C7", gap_zscore_k=1.5)
    # 3 reference positions (no gap in ref), but col 3 has a gap in seq4
    assert result["n_sequences"] == 5
    ref_residues = [r for r in result["residues"]]
    assert len(ref_residues) == 3

def test_compute_metrics_conserved_position():
    result = compute_metrics(TINY_ALIGNED_FASTA, ref_accession="P0A1C7", gap_zscore_k=1.5)
    pos1 = result["residues"][0]
    assert pos1["position"] == 1
    assert pos1["amino_acid"] == "M"
    assert pos1["shannon_entropy"] == pytest.approx(0.0)
    assert pos1["mutability_class"] == "invariant"
