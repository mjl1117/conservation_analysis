import math
import json
import pytest
import numpy as np
from scripts.coevolution import (
    mutual_information,
    apply_apc,
    compute_coevolution,
    GAP_CHARS,
)


def test_mi_perfectly_correlated_columns():
    # col_i and col_j are identical → max MI for 2 types = 1.0 bit
    col_i = ["A"] * 5 + ["B"] * 5
    col_j = ["A"] * 5 + ["B"] * 5
    assert mutual_information(col_i, col_j) == pytest.approx(1.0)


def test_mi_independent_columns_near_zero():
    # col_i alternates A/B, col_j is constant C → independent (knowing col_i
    # tells you nothing about col_j, which is always C).
    # p(A,C)=0.5, p(B,C)=0.5, p(A)*p(C)=0.5*1=0.5 → log2(p_ij/p_i*p_j)=log2(1)=0
    col_i = ["A", "B"] * 5
    col_j = ["C"] * 10
    mi = mutual_information(col_i, col_j)
    assert mi == pytest.approx(0.0, abs=1e-9)


def test_mi_excludes_gaps():
    col_i = ["-"] * 3 + ["A"] * 5 + ["B"] * 5
    col_j = ["-"] * 3 + ["A"] * 5 + ["B"] * 5
    # gaps excluded → perfectly correlated
    assert mutual_information(col_i, col_j) == pytest.approx(1.0)


def test_mi_all_gaps_returns_zero():
    assert mutual_information(["-"] * 5, ["-"] * 5) == pytest.approx(0.0)


def test_apply_apc_shape_preserved():
    # Hand-crafted symmetric MI matrix with zero diagonal
    mi = np.array([
        [0.0, 0.8, 0.2, 0.1],
        [0.8, 0.0, 0.3, 0.1],
        [0.2, 0.3, 0.0, 0.05],
        [0.1, 0.1, 0.05, 0.0],
    ])
    apc_mi = apply_apc(mi)
    assert apc_mi.shape == (4, 4)
    # Diagonal of output is zero (self-MI undefined, clipped to 0)
    assert np.diag(apc_mi) == pytest.approx([0.0, 0.0, 0.0, 0.0], abs=1e-6)
    # All values are >= 0 (APC clipped at 0)
    assert (apc_mi >= 0).all()
    # The strongly-correlated pair (0,1) should remain the highest after APC
    assert apc_mi[0, 1] == apc_mi.max()


def test_compute_coevolution_valid_schema():
    import jsonschema
    from pathlib import Path
    aligned = (
        ">seq1\nMKKLL\n"
        ">seq2\nMKRLL\n"
        ">seq3\nMRKLL\n"
        ">seq4\nMKKLL\n"
        ">seq5\nMRRLL\n"
    )
    result = compute_coevolution(aligned, ref_accession="seq1", gap_zscore_k=1.5)
    schema_path = Path(__file__).parent.parent.parent / "schema/v1/coevolution.schema.json"
    schema = json.loads(schema_path.read_text())
    jsonschema.validate(result, schema)


def test_compute_coevolution_positions_are_sparse():
    # Column 2 (0-indexed) is all gaps in seq2 and seq3, but seq1 has G there.
    # gap_fraction for col 2 = 2/3, while other columns have gap_fraction = 0.
    # With k=0.0 and non-zero std, col 2 gets a positive z-score (> 0.0) and
    # is excluded from positions. Reference position 3 (col_idx=2, ref_aa=G)
    # therefore should NOT appear in result["positions"].
    aligned = (
        ">seq1\nMKGLL\n"
        ">seq2\nMK-LL\n"
        ">seq3\nMK-LL\n"
    )
    result = compute_coevolution(aligned, ref_accession="seq1", gap_zscore_k=0.0)
    # col 2: gap_fraction=2/3 → z-score > 0 → excluded when k=0.0
    # ref_pos 3 corresponds to col_idx=2 (G in seq1) → should be excluded
    assert 3 not in result["positions"]
