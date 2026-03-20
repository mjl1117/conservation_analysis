# cli/tests/test_fetch_homologs.py
import responses as resp_lib
import pytest
from scripts.fetch_homologs import (
    fetch_reference_sequence,
    fetch_homologs,
    write_fasta,
    UNIPROT_API,
    UNIREF_API,
)

MOCK_REF_FASTA = (
    ">sp|P0A1C7|PDUA_SALTY PduA OS=Salmonella typhimurium LT2\n"
    "MQQEALGMVETKGLTAAIEAADAMVKSANVMLVGYEKIGSGLVTVIVRGDVGAVKAATDAGAAAARNVGEVKAVHVIPRPHTDVEKILPKGISQ\n"
)

MOCK_MEMBER_FASTA = (
    ">UniRef50_P0A1C7 Propanediol utilization protein PduA n=100\n"
    "MQQEALGMVETKGLTAAIEAADAMVKSANVMLVGYEKIGSGLVTVIVRGDVGAVKAATDAGAAAARNVGEVKAVHVIPRPHTDVEKILPKGISQ\n"
    ">UniRef50_Q8ZRE4 PduA homolog n=50\n"
    "MQQEALGMVETKGLTAAIEAADAMVKSANVMLVGYEKIGSGLVTVIVRGDVGAVKAATDAGAAAARNVGEVKAVHVIPRPHTDVEKILP\n"
)

@resp_lib.activate
def test_fetch_reference_sequence_returns_header_and_seq():
    resp_lib.add(resp_lib.GET, f"{UNIPROT_API}/P0A1C7.fasta",
                 body=MOCK_REF_FASTA, status=200)
    header, seq = fetch_reference_sequence("P0A1C7")
    assert "P0A1C7" in header
    assert seq.startswith("MQQEALGMVE")
    assert len(seq) == 94

@resp_lib.activate
def test_fetch_homologs_returns_filtered_sequences():
    resp_lib.add(
        resp_lib.GET,
        f"{UNIREF_API}/UniRef50_P0A1C7/members",
        body=MOCK_MEMBER_FASTA,
        status=200,
        headers={"Link": ""},  # no next page
    )
    seqs = fetch_homologs("P0A1C7", cluster_identity=50, min_length=56, max_length=150)
    assert len(seqs) == 2
    for header, seq in seqs:
        assert 56 <= len(seq) <= 150

@resp_lib.activate
def test_fetch_homologs_deduplicates():
    # Two sequences with same sequence content → only one retained
    duplicate_fasta = MOCK_MEMBER_FASTA + (
        ">UniRef50_Q8ZRE4_dup duplicate\n"
        "MQQEALGMVETKGLTAAIEAADAMVKSANVMLVGYEKIGSGLVTVIVRGDVGAVKAATDAGAAAARNVGEVKAVHVIPRPHTDVEKILP\n"
    )
    resp_lib.add(resp_lib.GET, f"{UNIREF_API}/UniRef50_P0A1C7/members",
                 body=duplicate_fasta, status=200, headers={"Link": ""})
    seqs = fetch_homologs("P0A1C7", cluster_identity=50)
    seq_bodies = [s for _, s in seqs]
    assert len(seq_bodies) == len(set(seq_bodies))

def test_write_fasta_formats_correctly(tmp_path):
    seqs = [("header1", "MKKLLV"), ("header2", "MKKLLX")]
    out = tmp_path / "out.fasta"
    write_fasta(seqs, str(out))
    content = out.read_text()
    assert ">header1\nMKKLLV\n" in content
    assert ">header2\nMKKLLX\n" in content
