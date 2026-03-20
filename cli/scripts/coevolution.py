"""Compute mutual information co-evolution analysis from aligned FASTA."""
from __future__ import annotations
import io
import math
from collections import Counter

import numpy as np
from Bio import AlignIO

from scripts.conservation_metrics import GAP_CHARS, compute_gap_zscores, gap_fraction


def _consensus(col: list[str]) -> str:
    """Most common non-gap residue in a column, or '?' if all gaps."""
    residues = [a for a in col if a not in GAP_CHARS]
    return Counter(residues).most_common(1)[0][0] if residues else "?"


def mutual_information(col_i: list[str], col_j: list[str]) -> float:
    """Mutual information (bits) between two alignment columns, gaps excluded."""
    pairs = [
        (a, b) for a, b in zip(col_i, col_j)
        if a not in GAP_CHARS and b not in GAP_CHARS
    ]
    if not pairs:
        return 0.0
    n = len(pairs)
    freq_i = Counter(a for a, _ in pairs)
    freq_j = Counter(b for _, b in pairs)
    freq_ij = Counter(pairs)
    mi = 0.0
    for (a, b), cnt in freq_ij.items():
        p_ij = cnt / n
        p_i = freq_i[a] / n
        p_j = freq_j[b] / n
        mi += p_ij * math.log2(p_ij / (p_i * p_j))
    return float(max(0.0, mi))


def apply_apc(mi: np.ndarray) -> np.ndarray:
    """
    Apply Average Product Correction (Dunn et al. 2008) to a symmetric MI matrix.
    Row means and global mean are computed excluding the diagonal (self-MI is undefined).
    Uses numpy masked arrays to safely exclude diagonal without reshape assumptions.
    """
    n = mi.shape[0]
    # Mask diagonal: True = excluded from calculations
    masked = np.ma.array(mi, mask=np.eye(n, dtype=bool))
    mean_all = float(masked.mean())
    if mean_all == 0.0:
        return mi.copy()
    # Per-row off-diagonal mean: masked.mean(axis=1) ignores the masked diagonal
    mean_i = np.array(masked.mean(axis=1))  # shape (n,)
    apc = np.outer(mean_i, mean_i) / mean_all
    return np.maximum(0.0, mi - apc)


def compute_coevolution(
    aligned_fasta: str,
    ref_accession: str,
    gap_zscore_k: float = 1.5,
    top_n: int = 50,
    cluster: str = "UniRef50",  # reserved for future use; not currently included in output
) -> dict:
    """
    Compute MI co-evolution matrix with APC correction.
    Returns dict conforming to coevolution.schema.json v1.
    Positions with gap_zscore > gap_zscore_k are excluded (sparse positions array).
    """
    alignment = AlignIO.read(io.StringIO(aligned_fasta), "fasta")
    n_seqs = len(alignment)
    n_cols = alignment.get_alignment_length()

    # Find reference row to get residue numbering
    ref_row = next(
        (rec for rec in alignment if ref_accession in rec.id), alignment[0]
    )
    ref_seq = str(ref_row.seq)

    all_columns = [[str(alignment[i, col]) for i in range(n_seqs)] for col in range(n_cols)]
    all_gap_fracs = [gap_fraction(col) for col in all_columns]
    gap_zscores = compute_gap_zscores(all_gap_fracs)

    # Select non-gap reference positions that pass gap filter
    included: list[tuple[int, int]] = []  # (ref_residue_number, col_idx)
    ref_pos = 0
    for col_idx, ref_aa in enumerate(ref_seq):
        if ref_aa in GAP_CHARS:
            continue
        ref_pos += 1
        if gap_zscores[col_idx] <= gap_zscore_k:
            included.append((ref_pos, col_idx))

    positions = [rp for rp, _ in included]
    n_pos = len(positions)
    if n_pos == 0:
        raise ValueError(
            f"No positions passed the gap filter (gap_zscore_k={gap_zscore_k}). "
            "Try increasing gap_zscore_k."
        )

    # Build MI matrix over included columns
    included_cols = [all_columns[ci] for _, ci in included]
    mi_matrix = np.zeros((n_pos, n_pos))
    for i in range(n_pos):
        for j in range(i + 1, n_pos):
            val = mutual_information(included_cols[i], included_cols[j])
            mi_matrix[i, j] = mi_matrix[j, i] = val

    apc_matrix = apply_apc(mi_matrix)

    # Build top pairs
    pairs = []
    for i in range(n_pos):
        for j in range(i + 1, n_pos):
            if mi_matrix[i, j] > 0:
                pairs.append({
                    "position_i": positions[i],
                    "position_j": positions[j],
                    "aa_i": _consensus(included_cols[i]),
                    "aa_j": _consensus(included_cols[j]),
                    "mi_score": round(float(mi_matrix[i, j]), 6),
                    "mi_score_apc": round(float(apc_matrix[i, j]), 6),
                })
    pairs.sort(key=lambda p: p["mi_score_apc"], reverse=True)
    top_pairs = pairs[:top_n]

    return {
        "schema_version": "1",
        "accession": ref_accession,
        "method": "mutual_information_apc",
        "n_positions": n_pos,
        "positions": positions,
        "matrix": apc_matrix.tolist(),
        "top_pairs": top_pairs,
    }
