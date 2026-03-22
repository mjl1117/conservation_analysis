// webapp/src/tests/analysis.test.js
import { describe, it, expect } from "vitest";
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
