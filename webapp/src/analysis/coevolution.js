/**
 * Co-evolution: mutual information + APC correction.
 * Mirrors cli/scripts/coevolution.py exactly.
 */

import { GAP_CHARS, computeGapZscores, gapFraction } from "./conservation.js";

export function mutualInformation(colI, colJ) {
  const pairs = colI.map((a, k) => [a, colJ[k]])
    .filter(([a, b]) => !GAP_CHARS.has(a) && !GAP_CHARS.has(b));
  if (pairs.length === 0) return 0;
  const n = pairs.length;
  const freqI = {}, freqJ = {}, freqIJ = {};
  for (const [a, b] of pairs) {
    freqI[a] = (freqI[a] || 0) + 1;
    freqJ[b] = (freqJ[b] || 0) + 1;
    const key = `${a}|${b}`;
    freqIJ[key] = (freqIJ[key] || 0) + 1;
  }
  let mi = 0;
  for (const [key, cnt] of Object.entries(freqIJ)) {
    const [a, b] = key.split("|");
    const pIJ = cnt / n;
    const pI = freqI[a] / n;
    const pJ = freqJ[b] / n;
    mi += pIJ * Math.log2(pIJ / (pI * pJ));
  }
  return Math.max(0, mi);
}

export function applyApc(mi) {
  // mi is a 2D array (n×n), symmetric, zero diagonal
  const n = mi.length;
  // Compute row means and global mean excluding diagonal
  const rowMeans = mi.map((row, i) => {
    const offDiag = row.filter((_, j) => j !== i);
    return offDiag.reduce((a, b) => a + b, 0) / offDiag.length;
  });
  const allOffDiag = mi.flatMap((row, i) => row.filter((_, j) => j !== i));
  const globalMean = allOffDiag.reduce((a, b) => a + b, 0) / allOffDiag.length;
  if (globalMean === 0) return mi.map(row => [...row]);

  return mi.map((row, i) =>
    row.map((val, j) => {
      const apc = (rowMeans[i] * rowMeans[j]) / globalMean;
      return Math.max(0, val - apc);
    })
  );
}

/**
 * Compute co-evolution matrix from alignment.
 * @param {Array<{id: string, seq: string}>} alignment
 * @param {string} refAccession
 * @param {number} gapZscoreK
 * @param {number} topN
 * @returns {Object} conforming to coevolution.schema.json v1
 */
export function computeCoevolution(alignment, refAccession, gapZscoreK = 1.5, topN = 50) {
  const nSeqs = alignment.length;
  const nCols = alignment[0].seq.length;
  const refRow = alignment.find(r => r.id.includes(refAccession)) || alignment[0];

  const columns = Array.from({ length: nCols }, (_, ci) => alignment.map(r => r.seq[ci]));
  const allGapFracs = columns.map(gapFraction);
  const gapZscores = computeGapZscores(allGapFracs);

  // Included positions: non-gap in ref + gap_zscore <= k
  const included = [];
  let refPos = 0;
  for (let ci = 0; ci < nCols; ci++) {
    if (GAP_CHARS.has(refRow.seq[ci])) continue;
    refPos++;
    if (gapZscores[ci] <= gapZscoreK) included.push({ refPos, ci });
  }

  const positions = included.map(p => p.refPos);
  const nPos = positions.length;
  const includedCols = included.map(p => columns[p.ci]);

  // Build MI matrix
  const miMatrix = Array.from({ length: nPos }, () => Array(nPos).fill(0));
  for (let i = 0; i < nPos; i++) {
    for (let j = i + 1; j < nPos; j++) {
      const val = mutualInformation(includedCols[i], includedCols[j]);
      miMatrix[i][j] = miMatrix[j][i] = val;
    }
  }
  const apcMatrix = applyApc(miMatrix);

  // Top pairs
  const pairs = [];
  for (let i = 0; i < nPos; i++) {
    for (let j = i + 1; j < nPos; j++) {
      if (miMatrix[i][j] > 0) {
        const consensus = col => {
          const r = col.filter(a => !GAP_CHARS.has(a));
          if (!r.length) return "?";
          const f = {};
          r.forEach(a => f[a] = (f[a] || 0) + 1);
          return Object.entries(f).sort((a, b) => b[1] - a[1])[0][0];
        };
        pairs.push({
          position_i: positions[i],
          position_j: positions[j],
          aa_i: consensus(includedCols[i]),
          aa_j: consensus(includedCols[j]),
          mi_score: +miMatrix[i][j].toFixed(6),
          mi_score_apc: +apcMatrix[i][j].toFixed(6),
        });
      }
    }
  }
  pairs.sort((a, b) => b.mi_score_apc - a.mi_score_apc);

  return {
    schema_version: "1",
    accession: refAccession,
    method: "mutual_information_apc",
    n_positions: nPos,
    positions,
    matrix: apcMatrix,
    top_pairs: pairs.slice(0, topN),
  };
}
