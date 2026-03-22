// webapp/src/tests/api.test.js
import { describe, it, expect, vi, beforeEach } from "vitest";

// Mock global fetch
const mockFetch = vi.fn();
global.fetch = mockFetch;

import { fetchReferenceSequence, fetchHomologs, fetchPdbIdsForAccession } from "../api/uniprot.js";
import { filterExperimentalStructures } from "../api/rcsb.js";

beforeEach(() => { mockFetch.mockReset(); });

// --- uniprot.js ---
describe("fetchReferenceSequence", () => {
  it("parses FASTA response into {header, seq}", async () => {
    mockFetch.mockResolvedValueOnce({
      ok: true,
      text: async () => ">sp|P0A1C7|PDUA_SALTY PduA\nMQQEALGMVE\n",
    });
    const result = await fetchReferenceSequence("P0A1C7");
    expect(result.seq).toBe("MQQEALGMVE");
    expect(result.header).toContain("P0A1C7");
  });
  it("throws on HTTP error", async () => {
    mockFetch.mockResolvedValueOnce({ ok: false, status: 404, text: async () => "" });
    await expect(fetchReferenceSequence("BAD")).rejects.toThrow("404");
  });
});

describe("fetchHomologs", () => {
  it("returns {sequences, total} with sequence objects", async () => {
    mockFetch.mockResolvedValueOnce({
      ok: true,
      headers: { get: () => null },
      text: async () => ">UniRef50_P0A1C7 cluster\nMQQEALGMVE\n>Q8ZRE4 homolog\nMQQEALGMVD\n",
    });
    const result = await fetchHomologs("P0A1C7", 50, { minLength: 1, maxLength: 9999 });
    expect(result.sequences).toHaveLength(2);
    expect(result.sequences[0].seq).toBe("MQQEALGMVE");
    expect(result.total).toBe(2);
  });
  it("deduplicates identical sequences", async () => {
    mockFetch.mockResolvedValueOnce({
      ok: true,
      headers: { get: () => null },
      text: async () => ">seq1\nMKKLL\n>seq2\nMKKLL\n",
    });
    const result = await fetchHomologs("P0A1C7", 50, { minLength: 1, maxLength: 9999 });
    expect(result.sequences).toHaveLength(1);
    expect(typeof result.total).toBe("number");
  });
  it("returns sequences array and total count", async () => {
    mockFetch.mockResolvedValueOnce({
      ok: true,
      headers: { get: () => null },
      text: async () => ">seq1\nMKKLL\n>seq2\nMRKLL\n",
    });
    const result = await fetchHomologs("P0A1C7", 50, { minLength: 1, maxLength: 9999 });
    expect(result.sequences).toHaveLength(2);
    expect(result.total).toBe(2);
  });
});

describe("fetchPdbIdsForAccession", () => {
  it("extracts PDB cross-references", async () => {
    mockFetch.mockResolvedValueOnce({
      ok: true,
      json: async () => ({
        uniProtKBCrossReferences: [
          { database: "PDB", id: "3NGK", properties: [{ key: "Method", value: "X-RAY DIFFRACTION" }] },
          { database: "GO", id: "GO:001" },
        ]
      }),
    });
    const result = await fetchPdbIdsForAccession("P0A1C7");
    expect(result.pdb_ids).toEqual(["3NGK"]);
    expect(result.structure_methods).toContain("X-RAY DIFFRACTION");
  });
  it("returns empty arrays for no PDB references", async () => {
    mockFetch.mockResolvedValueOnce({
      ok: true,
      json: async () => ({ uniProtKBCrossReferences: [] }),
    });
    const result = await fetchPdbIdsForAccession("P0A1C7");
    expect(result.pdb_ids).toHaveLength(0);
    expect(result.has_pdb).toBe(false);
  });
});

// --- ebi_clustalo.js ---
import { submitAlignmentJob, fetchAlignmentResult } from "../api/ebi_clustalo.js";

describe("submitAlignmentJob", () => {
  it("returns job ID from response text", async () => {
    mockFetch.mockResolvedValueOnce({ ok: true, text: async () => "clustalo-R20250316-abc123\n" });
    const jobId = await submitAlignmentJob(">seq1\nMKK\n");
    expect(jobId).toBe("clustalo-R20250316-abc123");
  });
  it("throws on HTTP error", async () => {
    mockFetch.mockResolvedValueOnce({ ok: false, status: 503, text: async () => "" });
    await expect(submitAlignmentJob(">seq1\nMKK\n")).rejects.toThrow("503");
  });
});

describe("fetchAlignmentResult", () => {
  it("returns FASTA text from result endpoint", async () => {
    mockFetch.mockResolvedValueOnce({ ok: true, text: async () => ">seq1\nM-KK\n" });
    const fasta = await fetchAlignmentResult("job-123");
    expect(fasta).toContain(">seq1");
  });
});

// --- rcsb.js ---
import { fetchPdbStructure, fetchStructureMethod } from "../api/rcsb.js";

describe("fetchPdbStructure", () => {
  it("returns PDB text on success", async () => {
    mockFetch.mockResolvedValueOnce({ ok: true, text: async () => "ATOM      1  CA  MET\n" });
    const text = await fetchPdbStructure("3NGK");
    expect(text).toContain("ATOM");
  });
  it("throws on HTTP error", async () => {
    mockFetch.mockResolvedValueOnce({ ok: false, status: 404, text: async () => "" });
    await expect(fetchPdbStructure("XXXX")).rejects.toThrow("404");
  });
});

describe("fetchStructureMethod", () => {
  it("returns experimental method string", async () => {
    mockFetch.mockResolvedValueOnce({ ok: true, json: async () => ({ exptl: [{ method: "X-RAY DIFFRACTION" }] }) });
    const method = await fetchStructureMethod("3NGK");
    expect(method).toBe("X-RAY DIFFRACTION");
  });
  it("returns null on HTTP error", async () => {
    mockFetch.mockResolvedValueOnce({ ok: false });
    const method = await fetchStructureMethod("XXXX");
    expect(method).toBeNull();
  });
});

describe("filterExperimentalStructures", () => {
  it("keeps X-ray, NMR, cryo-EM, electron crystallography", () => {
    const structs = [
      { method: "X-RAY DIFFRACTION" },
      { method: "SOLUTION NMR" },
      { method: "ELECTRON MICROSCOPY" },
      { method: "ELECTRON CRYSTALLOGRAPHY" },
      { method: "THEORETICAL MODEL" },
      { method: "COMPUTATIONAL" },
    ];
    const filtered = filterExperimentalStructures(structs);
    expect(filtered).toHaveLength(4);
  });
  it("excludes computational methods", () => {
    const filtered = filterExperimentalStructures([{ method: "COMPUTATIONAL" }]);
    expect(filtered).toHaveLength(0);
  });
});
