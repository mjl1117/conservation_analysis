import json
import pytest
from scripts.phylogenetic_tree import build_nj_tree, build_tree_json

TINY_ALIGNED = (
    ">P0A1C7\nMKKLLVIGGSGSIGKSTLL\n"
    ">Q9X1X2\nMKKLLVIGGAGSIGKSTLL\n"
    ">Q8ZRE4\nMKKLLVIGGAGSIGKSTLL\n"
    ">P12345\nMRRLLVIGGAGSIGKATLL\n"
)

def test_build_nj_tree_returns_newick_string():
    newick = build_nj_tree(TINY_ALIGNED)
    assert isinstance(newick, str)
    assert newick.strip().endswith(";")
    # Newick has parentheses
    assert "(" in newick and ")" in newick

def test_build_nj_tree_contains_all_sequences():
    newick = build_nj_tree(TINY_ALIGNED)
    for seq_id in ["P0A1C7", "Q9X1X2", "Q8ZRE4", "P12345"]:
        assert seq_id in newick

def test_build_tree_json_valid_schema():
    import jsonschema
    from pathlib import Path
    newick = build_nj_tree(TINY_ALIGNED)
    result = build_tree_json(
        accession="P0A1C7",
        cluster="UniRef50",
        aligned_fasta=TINY_ALIGNED,
        newick=newick,
        pdb_map={"P0A1C7": {"pdb_ids": ["3NGK"], "methods": ["X-RAY DIFFRACTION"]}},
    )
    schema_path = Path(__file__).parent.parent.parent / "schema/v1/alignment.schema.json"
    schema = json.loads(schema_path.read_text())
    jsonschema.validate(result, schema)

def test_build_tree_json_marks_reference():
    newick = build_nj_tree(TINY_ALIGNED)
    result = build_tree_json(
        accession="P0A1C7",
        cluster="UniRef50",
        aligned_fasta=TINY_ALIGNED,
        newick=newick,
        pdb_map={},
    )
    ref_seq = next(s for s in result["sequences"] if s["is_reference"])
    assert "P0A1C7" in ref_seq["id"]
