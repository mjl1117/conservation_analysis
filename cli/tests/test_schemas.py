"""Tests for JSON schema validation."""
import json
from pathlib import Path
import jsonschema
import pytest

SCHEMA_DIR = Path(__file__).parent.parent.parent / "schema" / "v1"


def load_schema(name: str) -> dict:
    return json.loads((SCHEMA_DIR / name).read_text())


# ---------------------------------------------------------------------------
# conservation.schema.json
# ---------------------------------------------------------------------------

def test_conservation_schema_validates_good_data(sample_conservation_json):
    schema = load_schema("conservation.schema.json")
    jsonschema.validate(sample_conservation_json, schema)  # should not raise


def test_conservation_schema_rejects_missing_field(sample_conservation_json):
    schema = load_schema("conservation.schema.json")
    bad = dict(sample_conservation_json)
    del bad["accession"]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, schema)


def test_conservation_schema_rejects_bad_cluster(sample_conservation_json):
    schema = load_schema("conservation.schema.json")
    bad = dict(sample_conservation_json)
    bad["cluster"] = "UniRef42"  # not in enum
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, schema)


def test_conservation_schema_rejects_bad_mutability_class(sample_conservation_json):
    schema = load_schema("conservation.schema.json")
    bad = {**sample_conservation_json}
    bad["residues"] = [
        {**sample_conservation_json["residues"][0], "mutability_class": "unknown"}
    ]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, schema)


# ---------------------------------------------------------------------------
# alignment.schema.json
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_alignment_json():
    return {
        "schema_version": "1",
        "accession": "P0A1C7",
        "cluster": "UniRef50",
        "n_sequences": 2,
        "n_columns": 28,
        "fasta": ">seq1\nMKKLLVIGG\n>seq2\nMKKLLVIGG\n",
        "newick": "((seq1:0.1,seq2:0.1):0.0);",
        "sequences": [
            {
                "id": "seq1",
                "description": "Test sequence 1",
                "sequence": "MKKLLVIGG",
                "is_reference": True,
                "has_pdb": True,
                "pdb_ids": ["1PF9"],
                "structure_methods": ["X-ray"],
            },
            {
                "id": "seq2",
                "description": "Test sequence 2",
                "sequence": "MKKLLVIGG",
                "is_reference": False,
                "has_pdb": False,
                "pdb_ids": [],
                "structure_methods": [],
            },
        ],
    }


def test_alignment_schema_validates_good_data(sample_alignment_json):
    schema = load_schema("alignment.schema.json")
    jsonschema.validate(sample_alignment_json, schema)


def test_alignment_schema_rejects_missing_newick(sample_alignment_json):
    schema = load_schema("alignment.schema.json")
    bad = dict(sample_alignment_json)
    del bad["newick"]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, schema)


# ---------------------------------------------------------------------------
# coevolution.schema.json
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_coevolution_json():
    return {
        "schema_version": "1",
        "accession": "P0A1C7",
        "method": "mutual_information_apc",
        "n_positions": 2,
        "positions": [1, 2],
        "matrix": [[0.0, 0.5], [0.5, 0.0]],
        "top_pairs": [
            {
                "position_i": 1,
                "position_j": 2,
                "aa_i": "M",
                "aa_j": "K",
                "mi_score": 0.5,
                "mi_score_apc": 0.3,
            }
        ],
    }


def test_coevolution_schema_validates_good_data(sample_coevolution_json):
    schema = load_schema("coevolution.schema.json")
    jsonschema.validate(sample_coevolution_json, schema)


def test_coevolution_schema_rejects_bad_method(sample_coevolution_json):
    schema = load_schema("coevolution.schema.json")
    bad = {**sample_coevolution_json, "method": "pearson"}
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, schema)


# ---------------------------------------------------------------------------
# structural.schema.json
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_structural_json():
    return {
        "schema_version": "1",
        "reference_accession": "P0A1C7",
        "reference_pdb": "1PF9",
        "reference_chain": "A",
        "alignments": [
            {
                "query_id": "Q9X1X2",
                "query_pdb": "2A1P",
                "query_chain": "A",
                "tm_score": 0.95,
                "rmsd_global": 0.8,
                "method": "TM-align",
                "alignment_map": [
                    {"ref_resnum": 1, "query_resnum": 1, "rmsd": 0.2},
                    {"ref_resnum": 2, "query_resnum": 2, "rmsd": 0.3},
                ],
            }
        ],
    }


def test_structural_schema_validates_good_data(sample_structural_json):
    schema = load_schema("structural.schema.json")
    jsonschema.validate(sample_structural_json, schema)


def test_structural_schema_rejects_tm_score_out_of_range(sample_structural_json):
    schema = load_schema("structural.schema.json")
    bad = {**sample_structural_json}
    bad["alignments"] = [{**sample_structural_json["alignments"][0], "tm_score": 1.5}]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, schema)


# ---------------------------------------------------------------------------
# comparison.schema.json
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_comparison_json():
    return {
        "schema_version": "1",
        "protein_a": "P0A1C7",
        "protein_b": "Q9X1X2",
        "cluster": "UniRef50",
        "sequence": {
            "pairwise_identity": 0.72,
            "pairwise_similarity": 0.85,
            "aligned_a": "MKKLL-VIGG",
            "aligned_b": "MKKLLRVIGG",
        },
        "phylogenetic": {
            "tree_distance": 0.15,
            "lca_node": "node_42",
            "shared_cluster_fraction": 0.6,
        },
        "structural": {
            "available": True,
            "tm_score": 0.95,
            "rmsd_global": 0.8,
            "pdb_a": "1PF9",
            "pdb_b": "2A1P",
            "alignment_map": [],
        },
    }


def test_comparison_schema_validates_good_data(sample_comparison_json):
    schema = load_schema("comparison.schema.json")
    jsonschema.validate(sample_comparison_json, schema)


def test_comparison_schema_validates_structural_unavailable(sample_comparison_json):
    schema = load_schema("comparison.schema.json")
    # structural with only "available: false" should be valid
    minimal = {**sample_comparison_json, "structural": {"available": False}}
    jsonschema.validate(minimal, schema)


def test_comparison_schema_rejects_missing_protein_b(sample_comparison_json):
    schema = load_schema("comparison.schema.json")
    bad = dict(sample_comparison_json)
    del bad["protein_b"]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, schema)
