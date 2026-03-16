# Schema

This directory contains JSON Schema definitions (Draft-07) for all JSON files
produced by the `pdu-conservation` CLI. These schemas are the **source of
truth** for the data contract between the Python CLI and the React web app.

## Directory layout

```
schema/
└── v1/                         # Schema version 1
    ├── conservation.schema.json
    ├── alignment.schema.json
    ├── coevolution.schema.json
    ├── structural.schema.json
    └── comparison.schema.json
```

## Schemas

### `conservation.schema.json`

Output of `conservation_metrics.py`. Contains per-residue conservation
statistics computed from a multiple sequence alignment (MSA).

Key fields:
- `accession` — UniProt accession of the query/reference protein
- `cluster` — MSA redundancy level (`UniRef100`, `UniRef90`, or `UniRef50`)
- `residues` — array of per-position metrics (Shannon entropy, relative
  entropy, mutability, conservation, gap statistics, observed amino acid types,
  and a qualitative `mutability_class`)

### `alignment.schema.json`

Output of `run_alignment.py`. Contains the full multiple sequence alignment
together with a phylogenetic tree.

Key fields:
- `fasta` — aligned FASTA string (gap characters included)
- `newick` — Newick-format tree string
- `sequences` — per-sequence metadata (UniProt id, description, reference
  flag, associated PDB structures)

### `coevolution.schema.json`

Output of `coevolution.py`. Contains pairwise co-evolution scores computed via
mutual information with average-product correction (MI-APC).

Key fields:
- `method` — always `"mutual_information_apc"` in v1
- `positions` — list of alignment column indices included in the analysis
- `matrix` — symmetric N×N matrix of MI-APC scores
- `top_pairs` — ranked list of strongly co-evolving residue pairs

### `structural.schema.json`

Output of `structural_align.py`. Contains pairwise structural alignment
results relative to a single reference structure.

Key fields:
- `reference_pdb` / `reference_chain` — RCSB PDB id and chain of the reference
- `alignments` — per-query TM-score, global RMSD, alignment method, and a
  residue-level `alignment_map`

### `comparison.schema.json`

Output of `compare.py`. Aggregates sequence, phylogenetic, and structural
similarity between two proteins.

Key fields:
- `sequence` — pairwise identity/similarity and gapped alignment strings
- `phylogenetic` — tree distance, lowest common ancestor node,
  shared-cluster fraction
- `structural` — optional block (present only when PDB structures are
  available for both proteins); mirrors fields from `structural.schema.json`

## Versioning

The `schema_version` field (always the string `"1"` in v1) is present in
every output file. Future breaking changes will increment this version and
live in a new `schema/v2/` directory. Non-breaking additions (new optional
fields) may be made within a version.

## Validation

Schemas are validated in `cli/tests/test_schemas.py` using the `jsonschema`
Python package. Run with:

```bash
pytest cli/tests/test_schemas.py -v
```
