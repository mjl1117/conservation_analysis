const RCSB_DATA_API = "https://data.rcsb.org/rest/v1/core/entry";

const EXPERIMENTAL_METHODS = new Set([
  "X-RAY DIFFRACTION",
  "SOLUTION NMR",
  "ELECTRON MICROSCOPY",
  "ELECTRON CRYSTALLOGRAPHY",
]);
const EXCLUDED_METHODS = new Set(["THEORETICAL MODEL", "COMPUTATIONAL"]);

export function filterExperimentalStructures(structures) {
  return structures.filter(s => {
    const method = (s.method || "").toUpperCase();
    if (EXCLUDED_METHODS.has(method)) return false;
    return EXPERIMENTAL_METHODS.has(method);
  });
}

export async function fetchPdbStructure(pdbId) {
  const url = `https://files.rcsb.org/download/${pdbId}.pdb`;
  const res = await fetch(url);
  if (!res.ok) throw new Error(`PDB fetch failed for ${pdbId}: ${res.status}`);
  return res.text();
}

export async function fetchStructureMethod(pdbId) {
  const res = await fetch(`${RCSB_DATA_API}/${pdbId.toLowerCase()}`);
  if (!res.ok) return null;
  const data = await res.json();
  return data?.exptl?.[0]?.method || null;
}
