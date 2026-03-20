import pytest
from scripts.color_pdb import color_pdb_by_mutability, _replace_bfactor

PDB_LINES = [
    "ATOM      1  CA  MET A   1      10.000  10.000  10.000  1.00 50.00           C\n",
    "ATOM      2  CA  GLN A   2      11.000  11.000  11.000  1.00 50.00           C\n",
    "ATOM      3  CA  MET B   1      12.000  12.000  12.000  1.00 50.00           C\n",
    "HETATM    4  O   HOH A 101      20.000  20.000  20.000  1.00 30.00           O\n",
]

CONSERVATION = {
    "schema_version": "1",
    "accession": "P0A1C7",
    "n_sequences": 100,
    "gap_threshold_fraction": 0.3,
    "residues": [
        {"position": 1, "mutability": 0.0},
        {"position": 2, "mutability": 0.75},
    ],
}

def test_replace_bfactor_modifies_correct_columns():
    line = "ATOM      1  CA  MET A   1      10.000  10.000  10.000  1.00 50.00           C\n"
    result = _replace_bfactor(line, 0.42)
    assert "  0.42" in result
    # B-factor is columns 60-65 (0-indexed)
    assert float(result[60:66]) == pytest.approx(0.42)

def test_color_pdb_chain_A_only(tmp_path):
    pdb = tmp_path / "test.pdb"
    pdb.write_text("".join(PDB_LINES))
    result = color_pdb_by_mutability(str(pdb), CONSERVATION, chain="A")
    lines = result.splitlines()
    # Chain A residue 1 → mutability 0.00
    atom_lines = [l for l in lines if l.startswith("ATOM") and l[21] == "A"]
    assert float(atom_lines[0][60:66]) == pytest.approx(0.0)
    assert float(atom_lines[1][60:66]) == pytest.approx(0.75)

def test_color_pdb_chain_B_not_modified(tmp_path):
    pdb = tmp_path / "test.pdb"
    pdb.write_text("".join(PDB_LINES))
    result = color_pdb_by_mutability(str(pdb), CONSERVATION, chain="A")
    chain_b = [l for l in result.splitlines() if l.startswith("ATOM") and l[21] == "B"]
    # Chain B should still have original B-factor 50.00
    assert float(chain_b[0][60:66]) == pytest.approx(50.0)

def test_color_pdb_hetatm_retained(tmp_path):
    pdb = tmp_path / "test.pdb"
    pdb.write_text("".join(PDB_LINES))
    result = color_pdb_by_mutability(str(pdb), CONSERVATION, chain="A")
    hetatm_lines = [l for l in result.splitlines() if l.startswith("HETATM")]
    assert hetatm_lines  # HETATM is present
    assert float(hetatm_lines[0][60:66]) == pytest.approx(30.0)  # B-factor unchanged

def test_color_pdb_remark_header_written(tmp_path):
    pdb = tmp_path / "test.pdb"
    pdb.write_text("".join(PDB_LINES))
    result = color_pdb_by_mutability(str(pdb), CONSERVATION, chain="A")
    assert "REMARK  99 COLORED BY MUTABILITY" in result
    assert "P0A1C7" in result

def test_color_pdb_writes_output_file(tmp_path):
    pdb = tmp_path / "test.pdb"
    out = tmp_path / "colored.pdb"
    pdb.write_text("".join(PDB_LINES))
    color_pdb_by_mutability(str(pdb), CONSERVATION, chain="A", output_path=str(out))
    assert out.exists()
    assert "REMARK" in out.read_text()
