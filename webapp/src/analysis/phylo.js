/**
 * Phylogenetic tree utilities: Newick parser, p-distance matrix, neighbor-joining.
 * Mirrors cli/scripts/phylogenetic_tree.py (using p-distance, not Kimura 2p).
 */

import { GAP_CHARS } from "./conservation.js";

// --- Newick parser ---
// Returns a tree node: { name, length, children }
export function parseNewick(newick) {
  let pos = 0;
  const str = newick.trim().replace(/;$/, "");

  function parseNode() {
    const node = { name: "", length: 0, children: [] };
    if (str[pos] === "(") {
      pos++; // consume '('
      node.children.push(parseNode());
      while (str[pos] === ",") { pos++; node.children.push(parseNode()); }
      pos++; // consume ')'
    }
    // Parse label
    const labelStart = pos;
    while (pos < str.length && !"(),;:".includes(str[pos])) pos++;
    node.name = str.slice(labelStart, pos);
    // Parse branch length
    if (str[pos] === ":") {
      pos++;
      const lenStart = pos;
      while (pos < str.length && !"(),;".includes(str[pos])) pos++;
      node.length = parseFloat(str.slice(lenStart, pos)) || 0;
    }
    return node;
  }
  return parseNode();
}

// --- p-distance matrix ---
export function buildDistanceMatrix(sequences) {
  const n = sequences.length;
  const dm = Array.from({ length: n }, () => Array(n).fill(0));
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const dist = pDistance(sequences[i].seq, sequences[j].seq);
      dm[i][j] = dm[j][i] = dist;
    }
  }
  return dm;
}

function pDistance(seqA, seqB) {
  const len = Math.min(seqA.length, seqB.length);
  if (len === 0) return 0;
  let diff = 0;
  let compared = 0;
  for (let k = 0; k < len; k++) {
    if (GAP_CHARS.has(seqA[k]) || GAP_CHARS.has(seqB[k])) continue;
    compared++;
    if (seqA[k] !== seqB[k]) diff++;
  }
  return compared === 0 ? 0 : diff / compared;
}

// --- Neighbor-joining (Saitou & Nei 1987) ---
export function neighborJoining(names, dm) {
  let ids = [...names];
  let D = dm.map(row => [...row]);
  let labels = ids.map(name => `${name}`);

  while (ids.length > 2) {
    const size = ids.length;
    // Compute Q matrix
    const rowSums = D.map(row => row.reduce((a, b) => a + b, 0));
    let minQ = Infinity, minI = 0, minJ = 1;
    for (let i = 0; i < size; i++) {
      for (let j = i + 1; j < size; j++) {
        const q = (size - 2) * D[i][j] - rowSums[i] - rowSums[j];
        if (q < minQ) { minQ = q; minI = i; minJ = j; }
      }
    }
    // Branch lengths to new internal node
    const dij = D[minI][minJ];
    const liI = dij / 2 + (rowSums[minI] - rowSums[minJ]) / (2 * (size - 2));
    const liJ = dij - liI;
    // New node Newick label
    const newLabel = `(${labels[minI]}:${liI.toFixed(6)},${labels[minJ]}:${liJ.toFixed(6)})`;
    // Indices of nodes to keep (all except minI and minJ)
    const keepIndices = [];
    for (let k = 0; k < size; k++) {
      if (k !== minI && k !== minJ) keepIndices.push(k);
    }
    // Distances from each kept node to the new internal node
    const dist2new = keepIndices.map(k => (D[minI][k] + D[minJ][k] - dij) / 2);
    // Build replacement distance matrix
    const newSize = keepIndices.length;
    const newD = Array.from({ length: newSize + 1 }, () => Array(newSize + 1).fill(0));
    for (let ni = 0; ni < newSize; ni++) {
      for (let nj = 0; nj < newSize; nj++) {
        newD[ni][nj] = D[keepIndices[ni]][keepIndices[nj]];
      }
      newD[ni][newSize] = newD[newSize][ni] = dist2new[ni];
    }
    D = newD;
    labels = [...keepIndices.map(k => labels[k]), newLabel];
    ids = [...keepIndices.map(k => ids[k]), "__new__"];
  }
  // Final two nodes
  if (ids.length === 2) {
    return `(${labels[0]}:${(D[0][1]/2).toFixed(6)},${labels[1]}:${(D[0][1]/2).toFixed(6)});`;
  }
  return `${labels[0]};`;
}
