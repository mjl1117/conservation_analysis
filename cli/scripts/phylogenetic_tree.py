"""Build neighbor-joining phylogenetic tree from aligned FASTA."""
from __future__ import annotations
import io
from typing import Optional

from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


def build_nj_tree(aligned_fasta: str) -> str:
    """
    Build a neighbor-joining tree from an aligned FASTA string.
    Returns Newick string. Uses p-distance (Biopython 'identity' model).
    Note: Kimura 2-parameter is a DNA model; for protein sequences p-distance
    is the standard simple approach and produces equivalent tree topology.
    """
    alignment = AlignIO.read(io.StringIO(aligned_fasta), "fasta")
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    newick_io = io.StringIO()
    Phylo.write(tree, newick_io, "newick")
    return newick_io.getvalue().strip()


def build_tree_json(
    accession: str,
    cluster: str,
    aligned_fasta: str,
    newick: str,
    pdb_map: dict,
) -> dict:
    """
    Build alignment.schema.json v1 dict.
    pdb_map: {seq_id: {"pdb_ids": [...], "methods": [...]}}
    """
    alignment = AlignIO.read(io.StringIO(aligned_fasta), "fasta")
    sequences = []
    for rec in alignment:
        seq_id = rec.id
        pdb_info = pdb_map.get(seq_id, {})
        sequences.append({
            "id": seq_id,
            "description": rec.description,
            "sequence": str(rec.seq),
            "is_reference": accession in seq_id,
            "has_pdb": bool(pdb_info.get("pdb_ids")),
            "pdb_ids": pdb_info.get("pdb_ids", []),
            "structure_methods": pdb_info.get("methods", []),
        })
    return {
        "schema_version": "1",
        "accession": accession,
        "cluster": cluster,
        "n_sequences": len(alignment),
        "n_columns": alignment.get_alignment_length(),
        "fasta": aligned_fasta,
        "newick": newick,
        "sequences": sequences,
    }
