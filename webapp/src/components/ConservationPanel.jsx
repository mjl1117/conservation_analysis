// webapp/src/components/ConservationPanel.jsx
import React, { useState } from "react";
import {
  BarChart, Bar, XAxis, YAxis, Tooltip, ResponsiveContainer, ReferenceLine, Cell,
} from "recharts";

const CHARTS = [
  { key: "shannon_entropy", label: "Shannon Entropy", color: "#63cab7", domain: [0, Math.log2(20)] },
  { key: "rel_entropy",     label: "Relative Entropy", color: "#4a9eff", domain: [0, Math.log2(20)] },
  { key: "mutability",      label: "Mutability", color: "#f87171", domain: [0, 1] },
  { key: "conservation",    label: "Conservation", color: "#34d399", domain: [0, 1] },
  { key: "gap_fraction",    label: "Gap Fraction", color: "#fbbf24", domain: [0, 1] },
];

const CLASS_COLORS = {
  invariant: "#34d399", conserved: "#63cab7", variable: "#fbbf24", hypervariable: "#f87171",
};

export default function ConservationPanel({ data }) {
  const [activeChart, setActiveChart] = useState(CHARTS[0].key);

  if (!data?.residues) return null;
  const residues = data.residues;
  const gapThresholdZscore = data.gap_threshold_zscore ?? 1.5;
  const gapThresholdFraction = data.gap_threshold_fraction ?? 0.5;

  const chart = CHARTS.find(c => c.key === activeChart);
  const chartData = residues.map(r => ({
    position: r.position,
    value: r[activeChart],
    amino_acid: r.amino_acid,
    mutability_class: r.mutability_class,
  }));

  return (
    <div className="card">
      <div className="card-title">Conservation Analysis</div>
      <div className="tab-bar">
        {CHARTS.map(c => (
          <button key={c.key} className={`tab ${activeChart === c.key ? "active" : ""}`}
            onClick={() => setActiveChart(c.key)}>{c.label}</button>
        ))}
      </div>
      <div style={{ fontSize: "0.75rem", color: "var(--text-2)", marginBottom: "0.5rem" }}>
        {data.n_sequences} sequences · {residues.length} residues · cluster: {data.cluster}
        {activeChart === "gap_fraction" && (
          <span style={{ color: "var(--accent-am)", marginLeft: "1rem" }}>
            — threshold: {gapThresholdFraction.toFixed(3)} (Z={gapThresholdZscore})
          </span>
        )}
      </div>
      <ResponsiveContainer width="100%" height={220}>
        <BarChart data={chartData} margin={{ top: 5, right: 10, left: -10, bottom: 5 }}>
          <XAxis dataKey="position" tick={{ fontSize: 10, fill: "var(--text-2)" }}
            label={{ value: "Residue Position", position: "insideBottom", offset: -2, fill: "var(--text-2)", fontSize: 11 }} />
          <YAxis domain={chart.domain} tick={{ fontSize: 10, fill: "var(--text-2)" }} />
          <Tooltip
            contentStyle={{ background: "var(--card)", border: "1px solid var(--border)", borderRadius: "8px", fontSize: "0.8rem" }}
            formatter={(value, _, props) => [
              `${value.toFixed(4)}`,
              `${props.payload.amino_acid}${props.payload.position} (${props.payload.mutability_class})`
            ]}
          />
          {activeChart === "gap_fraction" && (
            <ReferenceLine y={gapThresholdFraction} stroke="var(--accent-am)"
              strokeDasharray="4 2" label={{ value: "threshold", fill: "var(--accent-am)", fontSize: 10 }} />
          )}
          <Bar dataKey="value" radius={[2, 2, 0, 0]}>
            {chartData.map((entry, i) => (
              <Cell key={i}
                fill={activeChart === "mutability" || activeChart === "conservation"
                  ? CLASS_COLORS[entry.mutability_class]
                  : chart.color}
                opacity={0.85}
              />
            ))}
          </Bar>
        </BarChart>
      </ResponsiveContainer>
    </div>
  );
}
