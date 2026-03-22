// webapp/src/components/StructureViewer.jsx
import React, { useEffect, useRef, useState } from "react";
import { fetchPdbStructure } from "../api/rcsb.js";

const COLOR_SCHEMES = [
  { label: "Mutability", value: "mutability" },
  { label: "Spectrum", value: "spectrum" },
  { label: "Chain", value: "chain" },
];

export default function StructureViewer({ conservation, alignment }) {
  const viewerRef = useRef(null);
  const viewerObj = useRef(null);
  const [colorScheme, setColorScheme] = useState("mutability");
  const [pdbId, setPdbId] = useState("");
  const [pdbText, setPdbText] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  // Find the first experimental PDB from alignment data
  const defaultPdb = alignment?.sequences?.find(s => s.has_pdb)?.pdb_ids?.[0];

  useEffect(() => {
    if (defaultPdb && !pdbId) setPdbId(defaultPdb);
  }, [defaultPdb]);

  useEffect(() => {
    if (!viewerRef.current) return;
    // Initialize 3Dmol viewer — window.$3Dmol injected by CDN script
    if (!window.$3Dmol) return;
    viewerObj.current = window.$3Dmol.createViewer(viewerRef.current, {
      backgroundColor: "#080c14",
    });
    return () => viewerObj.current?.clear();
  }, []);

  const loadStructure = async (id) => {
    if (!id || !viewerObj.current || !window.$3Dmol) return;
    setLoading(true);
    setError(null);
    try {
      const text = await fetchPdbStructure(id);
      setPdbText(text);
      const viewer = viewerObj.current;
      viewer.clear();
      viewer.addModel(text, "pdb");
      applyColoring(viewer, colorScheme, conservation);
      viewer.zoomTo();
      viewer.render();
    } catch (err) {
      setError(err.message);
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    if (pdbId) loadStructure(pdbId);
  }, [pdbId, colorScheme, conservation]);

  return (
    <div className="card">
      <div className="card-title">Structure Viewer</div>
      <div style={{ display: "flex", gap: "0.75rem", alignItems: "center", marginBottom: "1rem", flexWrap: "wrap" }}>
        <input
          type="text" value={pdbId} onChange={e => setPdbId(e.target.value.toUpperCase())}
          placeholder="PDB ID (e.g. 3NGK)"
          style={{ background: "var(--surface)", border: "1px solid var(--border)",
            borderRadius: "6px", color: "var(--text)", padding: "0.35rem 0.7rem",
            fontSize: "0.85rem", width: "120px" }}
          onKeyDown={e => e.key === "Enter" && loadStructure(pdbId)}
        />
        <button className="btn btn-secondary" onClick={() => loadStructure(pdbId)} disabled={loading}>
          {loading ? "Loading…" : "Load"}
        </button>
        <div className="tab-bar" style={{ margin: 0 }}>
          {COLOR_SCHEMES.map(cs => (
            <button key={cs.value} className={`tab ${colorScheme === cs.value ? "active" : ""}`}
              onClick={() => setColorScheme(cs.value)}>{cs.label}</button>
          ))}
        </div>
        <button className="btn btn-secondary"
          onClick={() => downloadColoredPdb(pdbText, conservation, pdbId)}
          disabled={!pdbText || !conservation}>
          Download Colored PDB
        </button>
      </div>
      {error && <div className="error-banner">{error}</div>}
      <div ref={viewerRef} style={{ width: "100%", height: "450px", borderRadius: "8px", overflow: "hidden", background: "#080c14" }} />
      <div style={{ fontSize: "0.75rem", color: "var(--text-2)", marginTop: "0.5rem" }}>
        Color scale: <span style={{ color: "#4a9eff" }}>blue = conserved</span> · <span style={{ color: "#f87171" }}>red = hypervariable</span>
      </div>
    </div>
  );
}

function applyColoring(viewer, scheme, conservation) {
  const mutMap = {};
  if (conservation?.residues) {
    conservation.residues.forEach(r => { mutMap[r.position] = r.mutability; });
  }
  if (scheme === "mutability" && Object.keys(mutMap).length > 0) {
    viewer.setStyle({}, { cartoon: {
      colorfunc: (atom) => {
        const m = mutMap[atom.resi] ?? 0.5;
        return blueWhiteRed(m);
      }
    }});
  } else if (scheme === "spectrum") {
    viewer.setStyle({}, { cartoon: { color: "spectrum" } });
  } else {
    viewer.setStyle({}, { cartoon: { color: "chain" } });
  }
}

function blueWhiteRed(t) {
  // t=0 → blue (#4a9eff), t=0.5 → white, t=1 → red (#f87171)
  const r = t < 0.5 ? Math.round(74 + (255 - 74) * (t * 2)) : 255;
  const g = t < 0.5 ? Math.round(158 + (255 - 158) * (t * 2)) : Math.round(255 - (255 - 113) * ((t - 0.5) * 2));
  const b = t < 0.5 ? 255 : Math.round(255 - (255 - 113) * ((t - 0.5) * 2));
  return `rgb(${r},${g},${b})`;
}

function downloadColoredPdb(pdbText, conservation, pdbId) {
  if (!pdbText || !conservation || !pdbId) return;
  const mutMap = {};
  conservation.residues.forEach(r => { mutMap[r.position] = r.mutability * 100; });

  const remark = [
    `REMARK   1 PDU CONSERVATION COLORING`,
    `REMARK   1 ACCESSION: ${conservation.accession}`,
    `REMARK   1 CLUSTER: ${conservation.cluster}`,
    `REMARK   1 N_SEQUENCES: ${conservation.n_sequences}`,
    `REMARK   1 SCHEME: mutability (0=conserved, 100=hypervariable) in B-factor column`,
    `REMARK   1 GAP_THRESHOLD_FRACTION: ${conservation.gap_threshold_fraction}`,
  ].join("\n") + "\n";

  const colored = pdbText.split("\n").map(line => {
    const rec = line.slice(0, 6).trim();
    if (rec !== "ATOM") return line;
    const resSeq = parseInt(line.slice(22, 26), 10);
    const bfac = mutMap[resSeq] ?? 50;
    return line.slice(0, 60) + String(bfac.toFixed(2)).padStart(6) + line.slice(66);
  }).join("\n");

  const blob = new Blob([remark + colored], { type: "text/plain" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = `${pdbId}_colored.pdb`;
  a.click(); URL.revokeObjectURL(url);
}
