#!/usr/bin/env python3
"""PDU Conservation Analysis — Unified CLI."""
import json
import sys
from pathlib import Path

import click

# Add scripts dir to path for relative imports
sys.path.insert(0, str(Path(__file__).parent))

from scripts.fetch_homologs import fetch_reference_sequence, fetch_homologs, write_fasta
from scripts.run_alignment import run_clustalo, check_clustalo_available
from scripts.conservation_metrics import compute_metrics
from scripts.phylogenetic_tree import build_nj_tree, build_tree_json
from scripts.coevolution import compute_coevolution
from scripts.color_pdb import color_pdb_by_mutability
from scripts.compare import (
    compute_sequence_comparison, compute_phylogenetic_comparison,
    build_comparison_json, pairwise_align,
)

CLUSTER_MAP = {"50": "UniRef50", "90": "UniRef90", "100": "UniRef100"}


@click.group()
def cli():
    """PDU Conservation Analysis CLI.\n\nExample: python pdu.py run P0A1C7 --cluster 50 --out ./results/"""
    pass


@cli.command()
@click.argument("accession")
@click.option("--cluster", type=click.Choice(["50", "90", "100"]), default="50", show_default=True)
@click.option("--out", type=Path, default=Path("."), show_default=True, help="Output directory")
@click.option("--gap-zscore-k", default=1.5, show_default=True)
def run(accession, cluster, out, gap_zscore_k):
    """Run full pipeline: fetch → align → metrics → phylo → coevo."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    click.echo(f"[1/5] Fetching sequences for {accession} (UniRef{cluster})...")
    ref_header, ref_seq = fetch_reference_sequence(accession)
    homologs = fetch_homologs(accession, cluster_identity=int(cluster))
    all_seqs = [(ref_header, ref_seq)] + homologs
    fasta_path = out / "homologs.fasta"
    write_fasta(all_seqs, str(fasta_path))
    click.echo(f"    → {len(all_seqs)} sequences written to {fasta_path}")

    click.echo("[2/5] Running ClustalOmega alignment...")
    aligned_path = out / "aligned.fasta"
    run_clustalo(str(fasta_path), str(aligned_path))
    aligned_fasta = aligned_path.read_text()

    click.echo("[3/5] Computing conservation metrics...")
    conservation = compute_metrics(aligned_fasta, accession, gap_zscore_k=gap_zscore_k,
                                   cluster=CLUSTER_MAP[cluster])
    (out / "conservation.json").write_text(json.dumps(conservation, indent=2))

    click.echo("[4/5] Building phylogenetic tree...")
    newick = build_nj_tree(aligned_fasta)
    tree_json = build_tree_json(accession, CLUSTER_MAP[cluster], aligned_fasta, newick, {})
    (out / "tree.json").write_text(json.dumps(tree_json, indent=2))

    click.echo("[5/5] Computing co-evolution...")
    coevo = compute_coevolution(aligned_fasta, accession, gap_zscore_k=gap_zscore_k,
                                cluster=CLUSTER_MAP[cluster])
    (out / "coevolution.json").write_text(json.dumps(coevo, indent=2))

    click.echo(f"\n✓ Done. Results in {out}/")


@cli.command()
@click.argument("accession")
@click.option("--cluster", type=click.Choice(["50", "90", "100"]), default="50")
@click.option("--out", default="homologs.fasta")
def fetch(accession, cluster, out):
    """Fetch homolog sequences from UniRef cluster to a FASTA file."""
    ref_header, ref_seq = fetch_reference_sequence(accession)
    homologs = fetch_homologs(accession, cluster_identity=int(cluster))
    all_seqs = [(ref_header, ref_seq)] + homologs
    write_fasta(all_seqs, out)
    click.echo(f"Wrote {len(all_seqs)} sequences to {out}")


@cli.command()
@click.argument("input_fasta")
@click.option("--out", default="aligned.fasta")
@click.option("--threads", default=4)
def align(input_fasta, out, threads):
    """Run ClustalOmega alignment on a FASTA file."""
    run_clustalo(input_fasta, out, threads=threads)
    click.echo(f"Aligned FASTA written to {out}")


@cli.command()
@click.argument("aligned_fasta")
@click.option("--ref", required=True, help="Reference UniProt accession")
@click.option("--out", default="conservation.json")
@click.option("--gap-zscore-k", default=1.5)
@click.option("--cluster", default="UniRef50")
def metrics(aligned_fasta, ref, out, gap_zscore_k, cluster):
    """Compute conservation metrics from aligned FASTA."""
    fasta_text = Path(aligned_fasta).read_text()
    result = compute_metrics(fasta_text, ref, gap_zscore_k=gap_zscore_k, cluster=cluster)
    Path(out).write_text(json.dumps(result, indent=2))
    click.echo(f"Conservation metrics written to {out}")


@cli.command()
@click.argument("aligned_fasta")
@click.option("--ref", required=True)
@click.option("--out", default="tree.json")
def phylo(aligned_fasta, ref, out):
    """Build NJ phylogenetic tree from aligned FASTA."""
    fasta_text = Path(aligned_fasta).read_text()
    newick = build_nj_tree(fasta_text)
    tree_json = build_tree_json(ref, "UniRef50", fasta_text, newick, {})
    Path(out).write_text(json.dumps(tree_json, indent=2))
    click.echo(f"Tree JSON written to {out}")


@cli.command()
@click.argument("aligned_fasta")
@click.option("--ref", required=True)
@click.option("--out", default="coevolution.json")
@click.option("--gap-zscore-k", default=1.5)
@click.option("--top-n", default=50)
def coevo(aligned_fasta, ref, out, gap_zscore_k, top_n):
    """Compute mutual information co-evolution analysis."""
    fasta_text = Path(aligned_fasta).read_text()
    result = compute_coevolution(fasta_text, ref, gap_zscore_k=gap_zscore_k, top_n=top_n)
    Path(out).write_text(json.dumps(result, indent=2))
    click.echo(f"Co-evolution JSON written to {out}")


@cli.command()
@click.argument("conservation_json")
@click.argument("pdb_file")
@click.option("--out", default="colored.pdb")
@click.option("--chain", default="A")
def colorpdb(conservation_json, pdb_file, out, chain):
    """Color PDB B-factors by mutability from conservation.json."""
    conservation = json.loads(Path(conservation_json).read_text())
    color_pdb_by_mutability(pdb_file, conservation, chain=chain, output_path=out)
    click.echo(f"Colored PDB written to {out}")


@cli.command()
@click.argument("accession_a")
@click.argument("accession_b")
@click.option("--cluster", type=click.Choice(["50", "90", "100"]), default="50")
@click.option("--out", default="comparison.json")
@click.option("--conservation-a", default=None, help="Path to conservation.json for protein A")
@click.option("--conservation-b", default=None, help="Path to conservation.json for protein B")
@click.option("--tree-json", default=None, help="Path to tree.json containing both proteins")
def compare(accession_a, accession_b, cluster, out, conservation_a, conservation_b, tree_json):
    """Compare two proteins: sequence, phylogenetic, and structural summary."""
    click.echo(f"Comparing {accession_a} vs {accession_b}...")
    # Fetch and align both proteins together
    ref_a_header, ref_a_seq = fetch_reference_sequence(accession_a)
    ref_b_header, ref_b_seq = fetch_reference_sequence(accession_b)
    # Pairwise global alignment before position-wise comparison
    aligned_a, aligned_b = pairwise_align(ref_a_seq, ref_b_seq)
    seq_result = compute_sequence_comparison(aligned_a, aligned_b)
    # Phylogenetic comparison from tree.json if provided
    phylo_result = {"tree_distance": -1.0, "lca_node": "unknown", "shared_cluster_fraction": 0.0}
    if tree_json:
        tree_data = json.loads(Path(tree_json).read_text())
        phylo_result = compute_phylogenetic_comparison(
            tree_data["newick"], accession_a, accession_b
        )
    structural_result = {"available": False}
    result = build_comparison_json(
        protein_a=accession_a,
        protein_b=accession_b,
        cluster=CLUSTER_MAP[cluster],
        sequence=seq_result,
        phylogenetic=phylo_result,
        structural=structural_result,
    )
    Path(out).write_text(json.dumps(result, indent=2))
    click.echo(f"Comparison JSON written to {out}")


if __name__ == "__main__":
    cli()
