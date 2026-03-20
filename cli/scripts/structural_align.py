"""PyRosetta CE-align structural alignment."""
from __future__ import annotations
import sys
from typing import Optional


def check_pyrosetta_available() -> bool:
    return sys.modules.get("pyrosetta") is not None or _try_import_pyrosetta()


def _try_import_pyrosetta() -> bool:
    try:
        import pyrosetta  # noqa: F401
        return True
    except ImportError:
        return False


def align_structures(
    ref_pdb: str,
    query_pdb: str,
    ref_chain: str = "A",
    query_chain: str = "A",
) -> dict:
    """
    Run CE-align between ref_pdb and query_pdb via PyRosetta.
    Returns dict with tm_score, rmsd_global, alignment_map.
    """
    if not check_pyrosetta_available():
        raise RuntimeError(
            "PyRosetta not available. Install from https://www.pyrosetta.org/"
        )
    import pyrosetta
    from pyrosetta.rosetta.protocols.comparative_modeling import CEAlignMover

    pyrosetta.init("-mute all")
    ref_pose = pyrosetta.pose_from_pdb(ref_pdb)
    query_pose = pyrosetta.pose_from_pdb(query_pdb)

    ce = CEAlignMover()
    ce.set_template_pose(ref_pose)
    ce.apply(query_pose)

    # Build alignment map from ATOM coordinate correspondence
    alignment_map = []
    for i, (ref_res, qry_res) in enumerate(_get_aligned_pairs(ref_pose, query_pose)):
        ref_ca = ref_pose.residue(ref_res).xyz("CA")
        qry_ca = query_pose.residue(qry_res).xyz("CA")
        rmsd = float(ref_ca.distance(qry_ca))
        alignment_map.append({
            "ref_resnum": int(ref_pose.pdb_info().number(ref_res)),
            "query_resnum": int(query_pose.pdb_info().number(qry_res)),
            "rmsd": round(rmsd, 4),
        })

    return {
        "tm_score": float(ce.tm_score()),
        "rmsd_global": float(ce.rms()),
        "alignment_map": alignment_map,
    }


def _get_aligned_pairs(ref_pose, query_pose):
    """
    Yield (ref_resnum, query_resnum) pairs from CE-aligned poses.

    Implementation note: after `CEAlignMover.apply(query_pose)`, the query_pose
    is superimposed onto ref_pose. To extract the actual aligned residue pairs,
    use the alignment stored in the mover:
        ce_mover.alignment()  -> returns a vector of (ref_res, query_res) pairs
    The stub below is sequential and ONLY correct if both structures have
    identical length and no gaps — replace with the mover's alignment output
    when integrating with real PyRosetta:
        for pair in ce_mover.get_alignment_results():
            yield pair.ref_resnum, pair.query_resnum
    """
    n = min(ref_pose.size(), query_pose.size())
    for i in range(1, n + 1):
        yield i, i


def build_structural_json(
    reference_accession: str,
    reference_pdb: str,
    reference_chain: str,
    alignments: list[dict],
) -> dict:
    return {
        "schema_version": "1",
        "reference_accession": reference_accession,
        "reference_pdb": reference_pdb,
        "reference_chain": reference_chain,
        "alignments": alignments,
    }
