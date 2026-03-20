"""Build comparison.json between two proteins."""
from __future__ import annotations
import hashlib

from Bio import Align


def pairwise_align(seq_a: str, seq_b: str) -> tuple[str, str]:
    """
    Run global pairwise alignment on two raw sequences.
    Returns (aligned_a, aligned_b) as equal-length strings with '-' gap characters.
    Uses Biopython PairwiseAligner with BLOSUM62 scoring.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    alignment = next(iter(aligner.align(seq_a, seq_b)))
    aligned_a, aligned_b = _extract_aligned_strings(alignment, seq_a, seq_b)
    return aligned_a, aligned_b


def _extract_aligned_strings(alignment, seq_a: str, seq_b: str) -> tuple[str, str]:
    """
    Extract gap-containing aligned strings from a Biopython PairwiseAlignment object.
    alignment[0] and alignment[1] return the aligned target/query as strings with '-'
    gap characters — this is the documented PairwiseAlignment.__getitem__ API.
    """
    return str(alignment[0]), str(alignment[1])


def compute_sequence_comparison(aligned_a: str, aligned_b: str) -> dict:
    """
    Compute pairwise identity and similarity from pre-aligned sequences of equal length.
    For raw (unaligned) sequences, call pairwise_align() first.
    """
    if len(aligned_a) != len(aligned_b):
        raise ValueError(
            f"Sequences must be equal length for comparison. "
            f"Got {len(aligned_a)} and {len(aligned_b)}. "
            "Call pairwise_align(seq_a, seq_b) first."
        )
    n = len(aligned_a)
    matches = sum(a == b for a, b in zip(aligned_a, aligned_b))
    # Blosum62-based similarity: positive score = similar
    positive = _count_positive_subs(aligned_a, aligned_b)
    return {
        "pairwise_identity": matches / n,
        "pairwise_similarity": positive / n,
        "aligned_a": aligned_a,
        "aligned_b": aligned_b,
    }


def _count_positive_subs(seq_a: str, seq_b: str) -> int:
    """Count positions where amino acids are identical or biochemically similar."""
    SIMILAR_GROUPS = [
        set("IVLM"), set("FYW"), set("KRH"), set("DE"), set("NQ"), set("ST"),
    ]
    count = 0
    for a, b in zip(seq_a, seq_b):
        if a == b:
            count += 1
        elif any(a in g and b in g for g in SIMILAR_GROUPS):
            count += 1
    return count


def compute_phylogenetic_comparison(
    newick: str, accession_a: str, accession_b: str
) -> dict:
    """Compute tree distance and LCA node between two accessions in a Newick tree."""
    from Bio import Phylo
    import io
    tree = Phylo.read(io.StringIO(newick), "newick")

    # Initialize with safe defaults in case of parse failure
    tree_distance = -1.0
    lca_id = "unknown"
    leaves: list[str] = []

    try:
        path_a = tree.root.get_path(accession_a) or []
        path_b = tree.root.get_path(accession_b) or []
        # LCA: last common element in the paths
        lca = tree.root
        for a, b in zip(path_a, path_b):
            if a == b:
                lca = a
            else:
                break
        tree_distance = float(tree.distance(accession_a, accession_b))
        # LCA node ID: SHA-1 of sorted leaf accession IDs (first 8 hex chars)
        leaves = sorted(t.name for t in lca.get_terminals() if t.name)
        lca_id = hashlib.sha1("|".join(leaves).encode()).hexdigest()[:8]
    except Exception:
        pass

    total_leaves = len(tree.get_terminals())
    shared_frac = round(len(leaves) / total_leaves, 4) if total_leaves > 0 else 0.0

    return {
        "tree_distance": round(tree_distance, 6),
        "lca_node": lca_id,
        "shared_cluster_fraction": shared_frac,
    }


def build_comparison_json(
    protein_a: str,
    protein_b: str,
    cluster: str,
    sequence: dict,
    phylogenetic: dict,
    structural: dict,
) -> dict:
    return {
        "schema_version": "1",
        "protein_a": protein_a,
        "protein_b": protein_b,
        "cluster": cluster,
        "sequence": sequence,
        "phylogenetic": phylogenetic,
        "structural": structural,
    }
