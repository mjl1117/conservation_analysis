"""Fetch reference and homologous sequences from UniProtKB/UniRef."""
from __future__ import annotations
import re
import requests
from pathlib import Path
from typing import Optional

UNIPROT_API = "https://rest.uniprot.org/uniprotkb"
UNIREF_API  = "https://rest.uniprot.org/uniref"
CLUSTER_MAP = {100: "UniRef100", 90: "UniRef90", 50: "UniRef50"}


def fetch_reference_sequence(accession: str) -> tuple[str, str]:
    """Return (header, sequence) for the reference UniProt entry."""
    url = f"{UNIPROT_API}/{accession}.fasta"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    lines = r.text.strip().splitlines()
    header = lines[0][1:]  # strip leading >
    seq = "".join(lines[1:])
    return header, seq


def _parse_fasta_text(text: str) -> list[tuple[str, str]]:
    """Parse a FASTA string into list of (header, sequence) tuples."""
    records = []
    current_header = None
    seq_parts: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header is not None:
                records.append((current_header, "".join(seq_parts)))
            current_header = line[1:]
            seq_parts = []
        else:
            seq_parts.append(line)
    if current_header is not None:
        records.append((current_header, "".join(seq_parts)))
    return records


def _next_page_url(response: requests.Response) -> Optional[str]:
    """Extract next page URL from Link header, or None."""
    link = response.headers.get("Link", "")
    match = re.search(r'<([^>]+)>;\s*rel="next"', link)
    return match.group(1) if match else None


def fetch_homologs(
    accession: str,
    cluster_identity: int = 50,
    min_length: int = 56,
    max_length: int = 150,
    max_seqs: Optional[int] = None,
) -> list[tuple[str, str]]:
    """
    Fetch all members of the UniRef cluster for accession.
    Returns list of (header, sequence) after length filtering and deduplication.
    """
    cluster_id = f"{CLUSTER_MAP[cluster_identity]}_{accession}"
    url = f"{UNIREF_API}/{cluster_id}/members"
    params: Optional[dict] = {"format": "fasta", "size": 500}

    raw: list[tuple[str, str]] = []
    while url:
        r = requests.get(url, params=params, timeout=60)
        r.raise_for_status()
        raw.extend(_parse_fasta_text(r.text))
        url = _next_page_url(r)
        params = None  # pagination URL already has params encoded

    # Length filter
    filtered = [(h, s) for h, s in raw if min_length <= len(s) <= max_length]

    # Deduplication by sequence
    seen: set[str] = set()
    deduped: list[tuple[str, str]] = []
    for h, s in filtered:
        if s not in seen:
            seen.add(s)
            deduped.append((h, s))

    if max_seqs is not None:
        deduped = deduped[:max_seqs]

    return deduped


def write_fasta(sequences: list[tuple[str, str]], output_path: str) -> None:
    """Write list of (header, sequence) tuples to a FASTA file."""
    lines = []
    for header, seq in sequences:
        lines.append(f">{header}")
        lines.append(seq)
    Path(output_path).write_text("\n".join(lines) + "\n")
