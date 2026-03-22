// webapp/src/App.jsx
import React, { useState, useCallback } from "react";
import AccessionForm from "./components/AccessionForm.jsx";
import FileUpload from "./components/FileUpload.jsx";
import ConservationPanel from "./components/ConservationPanel.jsx";
import PhyloTree from "./components/PhyloTree.jsx";
import StructureViewer from "./components/StructureViewer.jsx";
import CoevolutionPanel from "./components/CoevolutionPanel.jsx";
import ComparePanel from "./components/ComparePanel.jsx";

import { fetchReferenceSequence, fetchHomologs, fetchPdbIdsForAccession } from "./api/uniprot.js";
import { alignSequences } from "./api/ebi_clustalo.js";
import { computeMetrics } from "./analysis/conservation.js";
import { buildDistanceMatrix, neighborJoining } from "./analysis/phylo.js";
import { computeCoevolution } from "./analysis/coevolution.js";

function parseFastaToAlignment(fastaText) {
  const records = [];
  let header = null, seqParts = [];
  for (const line of fastaText.split("\n")) {
    const l = line.trim();
    if (!l) continue;
    if (l.startsWith(">")) {
      if (header !== null) records.push({ id: header.split(" ")[0], seq: seqParts.join("") });
      header = l.slice(1);
      seqParts = [];
    } else seqParts.push(l);
  }
  if (header !== null) records.push({ id: header.split(" ")[0], seq: seqParts.join("") });
  return records;
}

function downloadJson(data, filename) {
  const blob = new Blob([JSON.stringify(data, null, 2)], { type: "application/json" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a"); a.href = url; a.download = filename;
  a.click(); URL.revokeObjectURL(url);
}

function downloadText(text, filename) {
  const blob = new Blob([text], { type: "text/plain" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a"); a.href = url; a.download = filename;
  a.click(); URL.revokeObjectURL(url);
}

export default function App() {
  const [session, setSession] = useState({
    accession: null, cluster: "UniRef50",
    conservation: null, alignment: null,
    coevolution: null, structural: null, comparison: null,
  });
  const [loading, setLoading] = useState(false);
  const [loadingStep, setLoadingStep] = useState("");
  const [error, setError] = useState(null);
  const [samplingNotice, setSamplingNotice] = useState(null);

  const runPipeline = useCallback(async ({ accession, cluster }) => {
    setLoading(true);
    setError(null);
    setSamplingNotice(null);
    const clusterName = `UniRef${cluster}`;

    try {
      setLoadingStep("Fetching reference sequence…");
      const ref = await fetchReferenceSequence(accession);

      setLoadingStep("Fetching homologs…");
      const { sequences: homologs, total: homologTotal } = await fetchHomologs(accession, cluster);
      const allSeqs = [ref, ...homologs];
      if (homologs.length < homologTotal) {
        setSamplingNotice(`Subsampled ${homologs.length} of ${homologTotal} homologs for alignment.`);
      }

      setLoadingStep("Submitting alignment to ClustalOmega (this may take a few minutes)…");
      const alignedFasta = await alignSequences(allSeqs.slice(0, 500), ({ elapsed }) => {
        const mins = Math.floor(elapsed / 60000);
        setLoadingStep(`Waiting for ClustalOmega alignment… (${mins}m elapsed, max 10m)`);
      });

      setLoadingStep("Computing conservation metrics…");
      const alignedSeqs = parseFastaToAlignment(alignedFasta);
      const conservation = computeMetrics(alignedSeqs, accession, 1.5, clusterName);

      setLoadingStep("Building phylogenetic tree…");
      const dm = buildDistanceMatrix(alignedSeqs);
      const ids = alignedSeqs.map(s => s.id);
      const newick = neighborJoining(ids, dm);

      setLoadingStep("Resolving PDB IDs…");
      const pdbMap = {};
      for (const seq of alignedSeqs.slice(0, 50)) {
        pdbMap[seq.id] = await fetchPdbIdsForAccession(seq.id.split("|")[1] || seq.id);
      }

      const alignment = {
        schema_version: "1",
        accession, cluster: clusterName,
        n_sequences: alignedSeqs.length,
        n_columns: alignedSeqs[0]?.seq.length || 0,
        fasta: alignedFasta, newick,
        sequences: alignedSeqs.map(s => ({
          id: s.id, description: s.id, sequence: s.seq,
          is_reference: s.id.includes(accession),
          ...(pdbMap[s.id] || { has_pdb: false, pdb_ids: [], structure_methods: [] }),
        })),
      };

      setLoadingStep("Computing co-evolution…");
      const coevolution = computeCoevolution(alignedSeqs, accession);

      setSession({ accession, cluster: clusterName, conservation, alignment, coevolution, structural: null, comparison: null });
    } catch (err) {
      setError(err.message);
    } finally {
      setLoading(false);
      setLoadingStep("");
    }
  }, []);

  const handleFileLoad = useCallback(({ type, data }) => {
    setSession(prev => ({ ...prev, [type]: data }));
  }, []);

  return (
    <div className="app">
      <header className="app-header"><h1>PDU Conservation Analysis</h1></header>
      <main className="app-main">
        {loading && <div className="loading-bar" />}
        {loadingStep && <div style={{ color: "var(--text-2)", fontSize: "0.85rem", marginBottom: "0.75rem" }}>{loadingStep}</div>}
        {samplingNotice && <div style={{ color: "var(--accent-am)", fontSize: "0.8rem", marginBottom: "0.75rem" }}>⚠ {samplingNotice}</div>}
        {error && <div className="error-banner">{error}</div>}
        <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: "1.5rem", marginBottom: "2rem" }}>
          <AccessionForm onSubmit={runPipeline} disabled={loading} />
          <FileUpload onLoad={handleFileLoad} />
        </div>

        {session.conservation && <ConservationPanel data={session.conservation} />}
        {session.alignment && <PhyloTree data={session.alignment} />}
        {(session.conservation || session.alignment) && (
          <StructureViewer conservation={session.conservation} alignment={session.alignment} />
        )}
        {session.coevolution && <CoevolutionPanel data={session.coevolution} />}
        <ComparePanel
          data={session.comparison}
          sessionA={{ conservation: session.conservation, alignment: session.alignment }}
          sessionB={null}
        />

        {session.accession && (
          <div className="card">
            <div className="card-title">Downloads</div>
            <div style={{ display: "flex", gap: "0.75rem", flexWrap: "wrap" }}>
              {session.conservation && (
                <button className="btn btn-secondary"
                  onClick={() => downloadJson(session.conservation, `${session.accession}_conservation.json`)}>
                  conservation.json
                </button>
              )}
              {session.alignment && (
                <>
                  <button className="btn btn-secondary"
                    onClick={() => downloadJson(session.alignment, `${session.accession}_tree.json`)}>
                    tree.json
                  </button>
                  <button className="btn btn-secondary"
                    onClick={() => downloadText(session.alignment.fasta, `${session.accession}_aligned.fasta`)}>
                    aligned.fasta
                  </button>
                  <button className="btn btn-secondary"
                    onClick={() => downloadText(session.alignment.newick, `${session.accession}_tree.nwk`)}>
                    tree.nwk
                  </button>
                </>
              )}
              {session.coevolution && (
                <button className="btn btn-secondary"
                  onClick={() => downloadJson(session.coevolution, `${session.accession}_coevolution.json`)}>
                  coevolution.json
                </button>
              )}
              {session.conservation && (
                <button className="btn btn-secondary"
                  onClick={() => {
                    const header = "position,amino_acid,shannon_entropy,rel_entropy,mutability,conservation,gap_fraction,gap_zscore,mutability_class";
                    const rows = session.conservation.residues.map(r =>
                      [r.position, r.amino_acid, r.shannon_entropy, r.rel_entropy, r.mutability,
                       r.conservation, r.gap_fraction, r.gap_zscore, r.mutability_class].join(",")
                    );
                    downloadText([header, ...rows].join("\n"), `${session.accession}_conservation.csv`);
                  }}>
                  conservation.csv
                </button>
              )}
            </div>
          </div>
        )}
      </main>
    </div>
  );
}
