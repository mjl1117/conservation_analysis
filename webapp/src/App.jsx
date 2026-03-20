import React, { useState } from "react";

export default function App() {
  // Session state — all results keyed by accession
  const [session, setSession] = useState({
    accession: null,
    cluster: "UniRef50",
    conservation: null,   // conservation.schema.json v1
    alignment: null,      // alignment.schema.json v1
    coevolution: null,    // coevolution.schema.json v1
    structural: null,     // structural.schema.json v1 (upload only)
    comparison: null,     // comparison.schema.json v1
  });
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  return (
    <div className="app">
      <header className="app-header">
        <h1>PDU Conservation Analysis</h1>
      </header>
      <main className="app-main">
        <p>Scaffold — components will be added in subsequent tasks.</p>
      </main>
    </div>
  );
}
