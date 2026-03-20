"""Replace PDB B-factors with per-residue mutability scores."""
from __future__ import annotations
from pathlib import Path
from typing import Optional


def _replace_bfactor(line: str, value: float) -> str:
    """Replace the B-factor field (cols 60-65) in a PDB ATOM line."""
    return line[:60] + f"{value:6.2f}" + line[66:]


def color_pdb_by_mutability(
    pdb_path: str,
    conservation_data: dict,
    chain: str = "A",
    output_path: Optional[str] = None,
) -> str:
    """
    Replace B-factors in pdb_path with mutability scores from conservation_data.
    Only ATOM records for the specified chain are modified.
    HETATM records are retained unchanged.
    Returns the modified PDB as a string.
    """
    mut_map = {r["position"]: r["mutability"] for r in conservation_data["residues"]}
    accession = conservation_data.get("accession", "unknown")
    n_seqs = conservation_data.get("n_sequences", "?")
    gap_thresh = conservation_data.get("gap_threshold_fraction", "?")

    # Any existing REMARK  99 lines in the input are filtered out and replaced with this header.
    output_lines = [
        "REMARK  99 COLORED BY MUTABILITY",
        f"REMARK  99 ACCESSION: {accession}",
        f"REMARK  99 N_SEQUENCES: {n_seqs}",
        f"REMARK  99 GAP_THRESHOLD_FRACTION: {gap_thresh}",
        f"REMARK  99 CHAIN: {chain}",
        "REMARK  99 COLOR SCALE: 0.00=blue(conserved) 1.00=red(hypervariable)",
    ]

    with open(pdb_path) as f:
        for line in f:
            stripped = line.rstrip("\n")
            if stripped.startswith("ATOM") and len(stripped) >= 66 and stripped[21] == chain:
                try:
                    resnum = int(stripped[22:26])
                except ValueError:
                    output_lines.append(stripped)
                    continue
                val = mut_map.get(resnum, 0.0)
                # ljust(80) normalises non-standard short ATOM lines to canonical PDB width
                output_lines.append(_replace_bfactor(stripped.ljust(80), val))
            elif not stripped.startswith("REMARK  99"):
                output_lines.append(stripped)

    result = "\n".join(output_lines) + "\n"
    if output_path:
        Path(output_path).write_text(result)
    return result
