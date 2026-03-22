// webapp/src/components/AccessionForm.jsx
import React, { useState } from "react";

const CLUSTERS = [
  { label: "100% identity (UniRef100)", value: 100 },
  { label: "90% identity (UniRef90)", value: 90 },
  { label: "50% identity (UniRef50)", value: 50 },
];

export default function AccessionForm({ onSubmit, disabled }) {
  const [accession, setAccession] = useState("");
  const [cluster, setCluster] = useState(50);

  const handleSubmit = (e) => {
    e.preventDefault();
    const trimmed = accession.trim().toUpperCase();
    if (!trimmed) return;
    onSubmit({ accession: trimmed, cluster });
  };

  return (
    <div className="card">
      <div className="card-title">Analyze Protein</div>
      <form onSubmit={handleSubmit} style={{ display: "flex", flexDirection: "column", gap: "1rem" }}>
        <div>
          <label style={{ display: "block", color: "var(--text-2)", fontSize: "0.8rem", marginBottom: "0.4rem" }}>
            UniProt Accession
          </label>
          <input
            type="text"
            value={accession}
            onChange={e => setAccession(e.target.value)}
            placeholder="e.g. P0A1C7"
            style={{
              background: "var(--surface)", border: "1px solid var(--border)",
              borderRadius: "8px", color: "var(--text)", padding: "0.6rem 0.9rem",
              fontSize: "1rem", width: "100%", outline: "none",
            }}
            disabled={disabled}
          />
        </div>
        <div>
          <label style={{ display: "block", color: "var(--text-2)", fontSize: "0.8rem", marginBottom: "0.4rem" }}>
            Sequence Cluster
          </label>
          <div style={{ display: "flex", gap: "0.75rem", flexWrap: "wrap" }}>
            {CLUSTERS.map(c => (
              <label key={c.value} style={{ display: "flex", alignItems: "center", gap: "0.4rem", cursor: "pointer", color: "var(--text-2)", fontSize: "0.85rem" }}>
                <input type="radio" name="cluster" value={c.value}
                  checked={cluster === c.value}
                  onChange={() => setCluster(c.value)}
                  disabled={disabled}
                />
                {c.label}
              </label>
            ))}
          </div>
        </div>
        <button type="submit" className="btn" disabled={disabled || !accession.trim()}>
          {disabled ? "Analyzing…" : "Analyze"}
        </button>
      </form>
    </div>
  );
}
