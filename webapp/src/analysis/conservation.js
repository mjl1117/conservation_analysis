// Conservation metrics — mirrors cli/scripts/conservation_metrics.py exactly

export const GAP_CHARS = new Set(["-", "X", "."]);
const LOG2_20 = Math.log2(20);

export function shannonEntropy(column) {
  const residues = column.filter(aa => !GAP_CHARS.has(aa));
  if (residues.length === 0) return 0.0;
  const n = residues.length;
  const counts = {};
  for (const aa of residues) counts[aa] = (counts[aa] || 0) + 1;
  return -Object.values(counts).reduce((sum, c) => sum + (c/n) * Math.log2(c/n), 0);
}

export function relEntropy(column) {
  const residues = column.filter(aa => !GAP_CHARS.has(aa));
  if (residues.length === 0) return 0.0;
  const n = residues.length;
  const counts = {};
  for (const aa of residues) counts[aa] = (counts[aa] || 0) + 1;
  const q = 1 / 20;
  const kl = Object.values(counts).reduce((sum, c) => sum + (c/n) * Math.log2((c/n) / q), 0);
  return Math.max(0, kl);
}

export function mutability(h) { return h / LOG2_20; }
export function conservationScore(m) { return 1.0 - m; }

export function gapFraction(column) {
  return column.filter(aa => GAP_CHARS.has(aa)).length / column.length;
}

export function computeGapZscores(gapFractions) {
  const n = gapFractions.length;
  const mean = gapFractions.reduce((a, b) => a + b, 0) / n;
  // Population std (ddof=0) — matches numpy default
  const variance = gapFractions.reduce((s, x) => s + (x - mean) ** 2, 0) / n;
  const std = Math.sqrt(variance);
  if (std === 0) return gapFractions.map(() => 0);
  return gapFractions.map(x => (x - mean) / std);
}

export function mutabilityClass(m) {
  if (m <= 0.20) return "invariant";
  if (m <= 0.45) return "conserved";
  if (m <= 0.70) return "variable";
  return "hypervariable";
}

/**
 * Compute per-residue conservation metrics from an alignment array.
 * @param {Array<{id: string, seq: string}>} alignment - aligned sequences
 * @param {string} refAccession - reference sequence ID
 * @param {number} gapZscoreK - threshold multiplier (default 1.5)
 * @param {string} cluster - UniRef cluster name
 * @returns {Object} conforming to conservation.schema.json v1
 */
export function computeMetrics(alignment, refAccession, gapZscoreK = 1.5, cluster = "UniRef50") {
  const nSeqs = alignment.length;
  const nCols = alignment[0].seq.length;
  const refRow = alignment.find(r => r.id.includes(refAccession)) || alignment[0];

  const columns = Array.from({ length: nCols }, (_, ci) =>
    alignment.map(r => r.seq[ci])
  );
  const allGapFracs = columns.map(gapFraction);
  const gapZscores = computeGapZscores(allGapFracs);
  const meanGap = allGapFracs.reduce((a, b) => a + b, 0) / allGapFracs.length;
  const stdGap = Math.sqrt(allGapFracs.reduce((s, x) => s + (x - meanGap) ** 2, 0) / allGapFracs.length);
  const gapThresholdFraction = Math.min(1.0, meanGap + gapZscoreK * stdGap);

  const residues = [];
  let refPos = 0;
  for (let ci = 0; ci < nCols; ci++) {
    const refAa = refRow.seq[ci];
    if (GAP_CHARS.has(refAa)) continue;
    refPos++;
    const col = columns[ci];
    const h = shannonEntropy(col);
    const re = relEntropy(col);
    const m = mutability(h);
    const nonGap = col.filter(aa => !GAP_CHARS.has(aa));
    const observed = [...new Set(nonGap)].sort();
    residues.push({
      position: refPos,
      amino_acid: refAa,
      shannon_entropy: +h.toFixed(6),
      rel_entropy: +re.toFixed(6),
      mutability: +m.toFixed(6),
      conservation: +conservationScore(m).toFixed(6),
      gap_fraction: +allGapFracs[ci].toFixed(6),
      gap_zscore: +gapZscores[ci].toFixed(6),
      aa_types_count: observed.length,
      aa_types_observed: observed,
      mutability_class: mutabilityClass(m),
    });
  }

  return {
    schema_version: "1",
    accession: refAccession,
    cluster,
    n_sequences: nSeqs,
    n_columns: nCols,
    gap_zscore_k: gapZscoreK,
    gap_threshold_zscore: gapZscoreK,
    gap_threshold_fraction: +gapThresholdFraction.toFixed(6),
    residues,
  };
}
