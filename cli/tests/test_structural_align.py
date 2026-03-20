import json
import pytest
from unittest.mock import patch, MagicMock
from scripts.structural_align import (
    check_pyrosetta_available,
    align_structures,
    build_structural_json,
)

def test_check_pyrosetta_available_true():
    with patch.dict("sys.modules", {"pyrosetta": MagicMock()}):
        assert check_pyrosetta_available() is True

def test_check_pyrosetta_available_false():
    import sys
    with patch.dict("sys.modules", {"pyrosetta": None}):
        assert check_pyrosetta_available() is False

def test_align_structures_raises_without_pyrosetta(tmp_path):
    pdb_a = tmp_path / "a.pdb"
    pdb_b = tmp_path / "b.pdb"
    pdb_a.write_text("ATOM...")
    pdb_b.write_text("ATOM...")
    with patch("scripts.structural_align.check_pyrosetta_available", return_value=False):
        with pytest.raises(RuntimeError, match="PyRosetta not available"):
            align_structures(str(pdb_a), str(pdb_b))

def test_build_structural_json_valid_schema():
    import jsonschema
    from pathlib import Path
    mock_results = [
        {
            "query_id": "Q8ZRE4",
            "query_pdb": "4RBT",
            "query_chain": "A",
            "tm_score": 0.94,
            "rmsd_global": 0.82,
            "method": "CE-align",
            "alignment_map": [
                {"ref_resnum": 1, "query_resnum": 1, "rmsd": 0.10},
                {"ref_resnum": 2, "query_resnum": 2, "rmsd": 0.30},
            ],
        }
    ]
    result = build_structural_json(
        reference_accession="P0A1C7",
        reference_pdb="3NGK",
        reference_chain="A",
        alignments=mock_results,
    )
    schema_path = Path(__file__).parent.parent.parent / "schema/v1/structural.schema.json"
    schema = json.loads(schema_path.read_text())
    jsonschema.validate(result, schema)
