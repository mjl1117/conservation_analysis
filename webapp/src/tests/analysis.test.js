// webapp/src/tests/analysis.test.js
import { describe, it, expect } from "vitest";
import { parseNewick, buildDistanceMatrix, neighborJoining } from "../analysis/phylo.js";
import { mutualInformation, applyApc, computeCoevolution } from "../analysis/coevolution.js";
import {
  shannonEntropy,
  relEntropy,
  mutability,
  conservationScore,
  gapFraction,
  computeGapZscores,
  mutabilityClass,
  computeMetrics,
  GAP_CHARS,
} from "../analysis/conservation.js";

describe("shannonEntropy", () => {
  it("returns 0 for all-identical column", () => {
    expect(shannonEntropy(["M","M","M","M"])).toBeCloseTo(0.0);
  });
  it("returns 1.0 for two equal-frequency types", () => {
    expect(shannonEntropy(["M","M","K","K"])).toBeCloseTo(1.0);
  });
  it("excludes gap chars (-, X, .)", () => {
    // 5 gaps + 5 M → only M present → entropy = 0
    expect(shannonEntropy(["-","-","X",".",".",  "M","M","M","M","M"])).toBeCloseTo(0.0);
  });
  it("returns 0 for all-gap column", () => {
    expect(shannonEntropy(["-","X","."])).toBeCloseTo(0.0);
  });
});

describe("relEntropy", () => {
  it("returns ~0 for uniform distribution over 20 AAs", () => {
    const aas = "ACDEFGHIKLMNPQRSTVWY".split("");
    const col = [...aas, ...aas, ...aas, ...aas, ...aas]; // 100 residues
    expect(relEntropy(col)).toBeCloseTo(0.0, 3);
  });
  it("returns log2(20) ≈ 4.32 for invariant column", () => {
    expect(relEntropy(["M","M","M","M","M"])).toBeCloseTo(Math.log2(20), 2);
  });
});

describe("mutability", () => {
  it("normalizes entropy to [0,1]", () => {
    const h = Math.log2(20); // maximum entropy
    expect(mutability(h)).toBeCloseTo(1.0);
    expect(mutability(0)).toBeCloseTo(0.0);
  });
});

describe("gapFraction", () => {
  it("counts -, X, and . as gaps", () => {
    expect(gapFraction(["-","X",".","M","K"])).toBeCloseTo(3/5);
  });
  it("returns 0 for no gaps", () => {
    expect(gapFraction(["M","K","L"])).toBeCloseTo(0.0);
  });
});

describe("computeGapZscores", () => {
  it("returns ordered zscores", () => {
    const zscores = computeGapZscores([0.0, 0.5, 1.0]);
    expect(zscores[0]).toBeLessThan(zscores[1]);
    expect(zscores[1]).toBeLessThan(zscores[2]);
    expect(zscores[0]).toBeCloseTo(-1.2247, 3);
  });
  it("returns zeros for all-same input", () => {
    const zscores = computeGapZscores([0.3, 0.3, 0.3]);
    zscores.forEach(z => expect(z).toBeCloseTo(0.0));
  });
});

describe("mutabilityClass", () => {
  it("classifies invariant [0, 0.20]", () => {
    expect(mutabilityClass(0.0)).toBe("invariant");
    expect(mutabilityClass(0.20)).toBe("invariant");
  });
  it("classifies conserved (0.20, 0.45]", () => {
    expect(mutabilityClass(0.21)).toBe("conserved");
    expect(mutabilityClass(0.45)).toBe("conserved");
  });
  it("classifies variable (0.45, 0.70]", () => {
    expect(mutabilityClass(0.46)).toBe("variable");
    expect(mutabilityClass(0.70)).toBe("variable");
  });
  it("classifies hypervariable (0.70, 1.0]", () => {
    expect(mutabilityClass(0.71)).toBe("hypervariable");
    expect(mutabilityClass(1.0)).toBe("hypervariable");
  });
});

describe("computeMetrics", () => {
  const tinyAlignment = [
    { id: "P0A1C7", seq: "MKK" },
    { id: "seq2",   seq: "MKK" },
    { id: "seq3",   seq: "MKK" },
    { id: "seq4",   seq: "MK-" },
    { id: "seq5",   seq: "MRL" },
  ];

  it("returns correct number of residues (reference non-gap positions)", () => {
    const result = computeMetrics(tinyAlignment, "P0A1C7");
    expect(result.residues).toHaveLength(3);
  });
  it("position 1 is invariant M", () => {
    const result = computeMetrics(tinyAlignment, "P0A1C7");
    const pos1 = result.residues[0];
    expect(pos1.position).toBe(1);
    expect(pos1.amino_acid).toBe("M");
    expect(pos1.mutability_class).toBe("invariant");
  });
  it("schema_version is '1'", () => {
    const result = computeMetrics(tinyAlignment, "P0A1C7");
    expect(result.schema_version).toBe("1");
  });
});

// --- phylo.js ---

describe("parseNewick", () => {
  it("parses a simple Newick string into a tree object", () => {
    const tree = parseNewick("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);");
    expect(tree).toBeDefined();
    expect(tree.children).toHaveLength(3);
  });
  it("leaf nodes have name and no children", () => {
    const tree = parseNewick("(A:0.1,B:0.2);");
    const leaves = [];
    const walk = (node) => { if (!node.children?.length) leaves.push(node.name); else node.children.forEach(walk); };
    walk(tree);
    expect(leaves).toContain("A");
    expect(leaves).toContain("B");
  });
});

describe("buildDistanceMatrix", () => {
  it("returns symmetric matrix with zero diagonal", () => {
    const seqs = [
      { id: "a", seq: "MKKLL" },
      { id: "b", seq: "MRKLL" },
      { id: "c", seq: "MRRLL" },
    ];
    const dm = buildDistanceMatrix(seqs);
    expect(dm.length).toBe(3);
    expect(dm[0][0]).toBeCloseTo(0);
    expect(dm[0][1]).toBeCloseTo(dm[1][0]); // symmetric
  });
  it("distance between identical sequences is 0", () => {
    const seqs = [{ id: "a", seq: "MKKLL" }, { id: "b", seq: "MKKLL" }];
    const dm = buildDistanceMatrix(seqs);
    expect(dm[0][1]).toBeCloseTo(0);
  });
});

describe("neighborJoining", () => {
  it("returns a Newick string ending with semicolon", () => {
    const seqs = [
      { id: "P0A1C7", seq: "MKKLLVIGG" },
      { id: "Q9X1X2", seq: "MKKLLVIGG" },
      { id: "Q8ZRE4", seq: "MKKLLVIAG" },
      { id: "P12345", seq: "MRRLLVIAG" },
    ];
    const dm = buildDistanceMatrix(seqs);
    const newick = neighborJoining(seqs.map(s => s.id), dm);
    expect(newick.endsWith(";")).toBe(true);
    expect(newick).toContain("P0A1C7");
    expect(newick).toContain("P12345");
  });
});

// --- coevolution.js ---

describe("mutualInformation", () => {
  it("returns 1.0 for perfectly correlated columns", () => {
    const colI = ["A","A","A","A","A","B","B","B","B","B"];
    const colJ = ["A","A","A","A","A","B","B","B","B","B"];
    expect(mutualInformation(colI, colJ)).toBeCloseTo(1.0);
  });
  it("returns 0 for all-gap column pair", () => {
    expect(mutualInformation(["-","-"],["-","-"])).toBeCloseTo(0.0);
  });
});

describe("applyApc", () => {
  it("returns matrix of same shape", () => {
    const mi = [
      [0, 0.8, 0.2, 0.1],
      [0.8, 0, 0.3, 0.1],
      [0.2, 0.3, 0, 0.05],
      [0.1, 0.1, 0.05, 0],
    ];
    const result = applyApc(mi);
    expect(result.length).toBe(4);
    expect(result[0].length).toBe(4);
  });
  it("diagonal is zero after APC", () => {
    const mi = [[0,0.8],[0.8,0]];
    const result = applyApc(mi);
    expect(result[0][0]).toBeCloseTo(0);
    expect(result[1][1]).toBeCloseTo(0);
  });
  it("all values are >= 0 (clipped)", () => {
    const mi = [[0,0.1],[0.1,0]];
    const result = applyApc(mi);
    result.forEach(row => row.forEach(v => expect(v).toBeGreaterThanOrEqual(0)));
  });
});

describe("computeCoevolution", () => {
  const TINY = [
    { id: "P0A1C7", seq: "MKKLL" },
    { id: "seq2",   seq: "MKKLL" },
    { id: "seq3",   seq: "MRKLL" },
    { id: "seq4",   seq: "MRRLL" },
    { id: "seq5",   seq: "MKK-L" },
  ];
  it("returns schema_version '1'", () => {
    const result = computeCoevolution(TINY, "P0A1C7");
    expect(result.schema_version).toBe("1");
  });
  it("positions array contains only non-gap reference positions", () => {
    const result = computeCoevolution(TINY, "P0A1C7");
    expect(Array.isArray(result.positions)).toBe(true);
    result.positions.forEach(p => expect(typeof p).toBe("number"));
  });
  it("matrix dimensions match positions length", () => {
    const result = computeCoevolution(TINY, "P0A1C7");
    expect(result.matrix.length).toBe(result.n_positions);
    expect(result.matrix[0].length).toBe(result.n_positions);
  });
  it("top_pairs entries have required fields", () => {
    const result = computeCoevolution(TINY, "P0A1C7");
    if (result.top_pairs.length > 0) {
      const pair = result.top_pairs[0];
      expect(pair).toHaveProperty("position_i");
      expect(pair).toHaveProperty("position_j");
      expect(pair).toHaveProperty("mi_score");
      expect(pair).toHaveProperty("mi_score_apc");
    }
  });
});
