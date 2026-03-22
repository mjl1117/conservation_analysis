const UNIPROT_API = "https://rest.uniprot.org/uniprotkb";
const UNIREF_API = "https://rest.uniprot.org/uniref";
const CLUSTER_MAP = { 100: "UniRef100", 90: "UniRef90", 50: "UniRef50" };

function parseFasta(text) {
  const records = [];
  let header = null, seqParts = [];
  for (const line of text.split("\n")) {
    const l = line.trim();
    if (!l) continue;
    if (l.startsWith(">")) {
      if (header !== null) records.push({ header, seq: seqParts.join("") });
      header = l.slice(1);
      seqParts = [];
    } else {
      seqParts.push(l);
    }
  }
  if (header !== null) records.push({ header, seq: seqParts.join("") });
  return records;
}

function nextPageUrl(headers) {
  const link = headers.get("Link") || "";
  const m = link.match(/<([^>]+)>;\s*rel="next"/);
  return m ? m[1] : null;
}

export async function fetchReferenceSequence(accession) {
  const url = `${UNIPROT_API}/${accession}.fasta`;
  const res = await fetch(url);
  if (!res.ok) throw new Error(`UniProt fetch failed: ${res.status}`);
  const text = await res.text();
  const records = parseFasta(text);
  return records[0];
}

export async function fetchHomologs(accession, clusterIdentity = 50, { minLength = 56, maxLength = 150, maxSeqs = 500, seed = 42 } = {}) {
  const clusterId = `${CLUSTER_MAP[clusterIdentity]}_${accession}`;
  let url = `${UNIREF_API}/${clusterId}/members?format=fasta&size=500`;
  const raw = [];

  while (url) {
    const res = await fetch(url);
    if (!res.ok) throw new Error(`UniRef fetch failed: ${res.status}`);
    raw.push(...parseFasta(await res.text()));
    url = nextPageUrl(res.headers);
  }

  // Filter and deduplicate
  const seen = new Set();
  const filtered = raw.filter(({ seq }) => {
    if (seq.length < minLength || seq.length > maxLength) return false;
    if (seen.has(seq)) return false;
    seen.add(seq);
    return true;
  });

  const total = filtered.length;
  // Reproducible Fisher-Yates shuffle with seeded PRNG, then take cap
  if (filtered.length > maxSeqs) {
    const rng = mulberry32(seed);
    for (let i = filtered.length - 1; i > 0; i--) {
      const j = Math.floor(rng() * (i + 1));
      [filtered[i], filtered[j]] = [filtered[j], filtered[i]];
    }
    return { sequences: filtered.slice(0, maxSeqs), total };
  }
  return { sequences: filtered, total };
}

// Seeded PRNG (mulberry32) for reproducible subsampling
function mulberry32(seed) {
  let s = seed;
  return () => {
    s |= 0; s = s + 0x6D2B79F5 | 0;
    let t = Math.imul(s ^ (s >>> 15), 1 | s);
    t = t + Math.imul(t ^ (t >>> 7), 61 | t) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

export async function fetchPdbIdsForAccession(accession) {
  const url = `${UNIPROT_API}/${accession}.json`;
  const res = await fetch(url);
  if (!res.ok) return { has_pdb: false, pdb_ids: [], structure_methods: [] };
  const data = await res.json();
  const pdbRefs = (data.uniProtKBCrossReferences || []).filter(r => r.database === "PDB");
  const pdb_ids = pdbRefs.map(r => r.id);
  const structure_methods = pdbRefs.flatMap(r =>
    (r.properties || []).filter(p => p.key === "Method").map(p => p.value)
  );
  return { has_pdb: pdb_ids.length > 0, pdb_ids, structure_methods };
}
