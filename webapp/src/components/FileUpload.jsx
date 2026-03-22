// webapp/src/components/FileUpload.jsx
import React, { useCallback, useState } from "react";

const SCHEMA_FILES = {
  conservation: "conservation.schema.json",
  alignment: "alignment.schema.json (tree.json)",
  coevolution: "coevolution.schema.json",
  structural: "structural.schema.json",
  comparison: "comparison.schema.json",
};

export default function FileUpload({ onLoad }) {
  const [dragging, setDragging] = useState(false);
  const [errors, setErrors] = useState([]);

  const processFile = useCallback(async (file) => {
    try {
      const text = await file.text();
      const data = JSON.parse(text);
      if (data.schema_version !== "1") {
        throw new Error(`Unsupported schema_version: ${data.schema_version}. Expected "1".`);
      }
      // Detect schema type by required fields
      let schemaType = null;
      if ("residues" in data && "mutability" in (data.residues?.[0] ?? {})) schemaType = "conservation";
      else if ("newick" in data && "sequences" in data) schemaType = "alignment";
      else if ("method" in data && data.method === "mutual_information_apc") schemaType = "coevolution";
      else if ("alignments" in data && "reference_pdb" in data) schemaType = "structural";
      else if ("protein_a" in data && "protein_b" in data) schemaType = "comparison";
      if (!schemaType) throw new Error("Could not detect schema type from file content.");
      onLoad({ type: schemaType, data });
    } catch (err) {
      setErrors(prev => [...prev, `${file.name}: ${err.message}`]);
    }
  }, [onLoad]);

  const handleDrop = useCallback((e) => {
    e.preventDefault();
    setDragging(false);
    setErrors([]);
    Array.from(e.dataTransfer.files).forEach(processFile);
  }, [processFile]);

  return (
    <div className="card">
      <div className="card-title">Load CLI Output</div>
      <div
        onDrop={handleDrop}
        onDragOver={e => { e.preventDefault(); setDragging(true); }}
        onDragLeave={() => setDragging(false)}
        style={{
          border: `2px dashed ${dragging ? "var(--accent)" : "var(--border)"}`,
          borderRadius: "8px", padding: "2rem", textAlign: "center",
          color: "var(--text-2)", fontSize: "0.85rem", cursor: "pointer",
          background: dragging ? "rgba(99,202,183,0.05)" : "transparent",
          transition: "all 0.2s",
        }}
        onClick={() => document.getElementById("file-input").click()}
      >
        <div>Drop JSON files here or click to browse</div>
        <div style={{ marginTop: "0.5rem", fontSize: "0.75rem", opacity: 0.7 }}>
          Accepts: {Object.values(SCHEMA_FILES).join(", ")}
        </div>
        <input id="file-input" type="file" accept=".json" multiple hidden
          onChange={e => Array.from(e.target.files).forEach(processFile)} />
      </div>
      {errors.map((err, i) => (
        <div key={i} className="error-banner" style={{ marginTop: "0.5rem" }}>{err}</div>
      ))}
    </div>
  );
}
