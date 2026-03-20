"""Compute per-residue conservation metrics from an aligned FASTA."""
from __future__ import annotations
import io
import json
import math
from collections import Counter
from pathlib import Path
from typing import Optional

import numpy as np
from Bio import AlignIO

GAP_CHARS = frozenset({"-", "X", "."})
LOG2_20 = math.log2(20)


def shannon_entropy(column: list[str]) -> float:
    residues = [aa for aa in column if aa not in GAP_CHARS]
    if not residues:
        return 0.0
    n = len(residues)
    counts = Counter(residues)
    return float(-sum((c / n) * math.log2(c / n) for c in counts.values()))


def rel_entropy(column: list[str]) -> float:
    """KL divergence from uniform distribution over 20 amino acids."""
    residues = [aa for aa in column if aa not in GAP_CHARS]
    if not residues:
        return 0.0
    n = len(residues)
    counts = Counter(residues)
    q = 1 / 20
    kl = sum((c / n) * math.log2((c / n) / q) for c in counts.values())
    return float(max(0.0, kl))


def mutability(h: float) -> float:
    return h / LOG2_20


def conservation_score(m: float) -> float:
    return 1.0 - m


def gap_fraction(column: list[str]) -> float:
    return sum(1 for aa in column if aa in GAP_CHARS) / len(column)


def compute_gap_zscores(gap_fractions: list[float]) -> list[float]:
    arr = np.array(gap_fractions, dtype=float)
    std = float(arr.std())
    if std == 0.0:
        return [0.0] * len(gap_fractions)
    return list((arr - arr.mean()) / std)


def mutability_class(m: float) -> str:
    if m <= 0.20:
        return "invariant"
    elif m <= 0.45:
        return "conserved"
    elif m <= 0.70:
        return "variable"
    return "hypervariable"


def compute_metrics(
    aligned_fasta: str,
    ref_accession: str,
    gap_zscore_k: float = 1.5,
    cluster: str = "UniRef50",
) -> dict:
    """
    Compute per-residue conservation metrics from an aligned FASTA string.
    Returns a dict conforming to conservation.schema.json v1.
    """
    alignment = AlignIO.read(io.StringIO(aligned_fasta), "fasta")
    n_seqs = len(alignment)
    n_cols = alignment.get_alignment_length()

    # Find reference row
    ref_row = next(
        (rec for rec in alignment if ref_accession in rec.id), alignment[0]
    )
    ref_seq = str(ref_row.seq)

    # Compute per-column metrics (only for non-gap reference positions)
    all_columns = [[str(alignment[i, col]) for i in range(n_seqs)] for col in range(n_cols)]
    all_gap_fracs = [gap_fraction(col) for col in all_columns]
    gap_zscores = compute_gap_zscores(all_gap_fracs)

    mean_gap = float(np.mean(all_gap_fracs))
    std_gap = float(np.std(all_gap_fracs))
    gap_threshold_fraction = min(1.0, mean_gap + gap_zscore_k * std_gap)

    residues = []
    ref_pos = 0
    for col_idx, ref_aa in enumerate(ref_seq):
        if ref_aa in GAP_CHARS:
            continue
        ref_pos += 1
        col = all_columns[col_idx]
        h = shannon_entropy(col)
        re = rel_entropy(col)
        m = mutability(h)
        non_gap_aas = [aa for aa in col if aa not in GAP_CHARS]
        observed = sorted(set(non_gap_aas))

        residues.append({
            "position": ref_pos,
            "amino_acid": ref_aa,
            "shannon_entropy": round(h, 6),
            "rel_entropy": round(re, 6),
            "mutability": round(m, 6),
            "conservation": round(conservation_score(m), 6),
            "gap_fraction": round(all_gap_fracs[col_idx], 6),
            "gap_zscore": round(gap_zscores[col_idx], 6),
            "aa_types_count": len(observed),
            "aa_types_observed": observed,
            "mutability_class": mutability_class(m),
        })

    return {
        "schema_version": "1",
        "accession": ref_accession,
        "cluster": cluster,
        "n_sequences": n_seqs,
        "n_columns": n_cols,
        "gap_zscore_k": gap_zscore_k,
        "gap_threshold_zscore": gap_zscore_k,
        "gap_threshold_fraction": round(gap_threshold_fraction, 6),
        "residues": residues,
    }
