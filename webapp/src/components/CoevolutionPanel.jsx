// webapp/src/components/CoevolutionPanel.jsx
import React, { useEffect, useRef, useState } from "react";
import * as d3 from "d3";

export default function CoevolutionPanel({ data }) {
  const [view, setView] = useState("heatmap"); // "heatmap" | "network"
  const [threshold, setThreshold] = useState(0.1);
  const heatmapRef = useRef(null);
  const networkRef = useRef(null);

  useEffect(() => {
    if (view === "heatmap" && data?.matrix && heatmapRef.current)
      renderHeatmap(heatmapRef.current, data);
  }, [view, data]);

  useEffect(() => {
    if (view === "network" && data?.top_pairs && networkRef.current)
      renderNetwork(networkRef.current, data, threshold);
  }, [view, data, threshold]);

  if (!data?.matrix) return null;

  return (
    <div className="card">
      <div className="card-title">Co-evolution (Mutual Information + APC)</div>
      <div style={{ display: "flex", gap: "0.75rem", alignItems: "center", marginBottom: "1rem" }}>
        <div className="tab-bar" style={{ margin: 0 }}>
          <button className={`tab ${view === "heatmap" ? "active" : ""}`} onClick={() => setView("heatmap")}>Heatmap</button>
          <button className={`tab ${view === "network" ? "active" : ""}`} onClick={() => setView("network")}>Network</button>
        </div>
        {view === "network" && (
          <label style={{ display: "flex", alignItems: "center", gap: "0.5rem", color: "var(--text-2)", fontSize: "0.8rem" }}>
            MI threshold:
            <input type="range" min="0" max="1" step="0.01" value={threshold}
              onChange={e => setThreshold(+e.target.value)}
              style={{ width: "100px" }} />
            {threshold.toFixed(2)}
          </label>
        )}
      </div>
      {view === "heatmap"
        ? <svg ref={heatmapRef} style={{ width: "100%", height: "500px" }} />
        : <svg ref={networkRef} style={{ width: "100%", height: "500px" }} />}
    </div>
  );
}

function renderHeatmap(svgEl, data) {
  const svg = d3.select(svgEl);
  svg.selectAll("*").remove();
  const W = svgEl.clientWidth || 600;
  const H = 480;
  const margin = { top: 30, right: 30, bottom: 60, left: 60 };
  const w = W - margin.left - margin.right;
  const h = H - margin.top - margin.bottom;
  const n = data.positions.length;
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

  const matrix = data.matrix;
  const allVals = matrix.flat().filter(v => v > 0);
  const maxVal = d3.max(allVals) || 1;

  const colorScale = d3.scaleSequential(d3.interpolateYlOrRd).domain([0, maxVal]);
  const cellW = w / n;
  const cellH = h / n;

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      g.append("rect")
        .attr("x", j * cellW).attr("y", i * cellH)
        .attr("width", cellW).attr("height", cellH)
        .attr("fill", matrix[i][j] > 0 ? colorScale(matrix[i][j]) : "var(--surface)");
    }
  }

  // Axis ticks (every 10 positions)
  const tickPositions = data.positions.filter((_, i) => i % Math.max(1, Math.floor(n / 10)) === 0);
  g.append("g").attr("transform", `translate(0,${h})`)
    .call(d3.axisBottom(d3.scaleLinear().domain([0, n]).range([0, w]))
      .tickValues(tickPositions.map((_, i) => i * Math.max(1, Math.floor(n / 10))))
      .tickFormat((_, i) => tickPositions[i] || ""))
    .selectAll("text").style("fill", "var(--text-2)").style("font-size", "9px");
}

function renderNetwork(svgEl, data, threshold) {
  const svg = d3.select(svgEl);
  svg.selectAll("*").remove();
  const W = svgEl.clientWidth || 600;
  const H = 480;

  const filteredPairs = data.top_pairs.filter(p => p.mi_score_apc >= threshold);
  if (filteredPairs.length === 0) {
    svg.append("text").attr("x", W / 2).attr("y", H / 2)
      .attr("text-anchor", "middle").attr("fill", "var(--text-2)")
      .attr("font-size", "14px")
      .text("No pairs above threshold — lower the MI threshold slider.");
    return;
  }
  const nodeIds = [...new Set(filteredPairs.flatMap(p => [p.position_i, p.position_j]))];
  const nodes = nodeIds.map(id => ({ id, position: id }));
  const links = filteredPairs.map(p => ({ source: p.position_i, target: p.position_j, value: p.mi_score_apc }));

  const sim = d3.forceSimulation(nodes)
    .force("link", d3.forceLink(links).id(d => d.id).strength(d => d.value))
    .force("charge", d3.forceManyBody().strength(-60))
    .force("center", d3.forceCenter(W / 2, H / 2));

  const g = svg.append("g");

  const link = g.selectAll(".link")
    .data(links).enter().append("line").attr("class", "link")
    .attr("stroke", "var(--border)").attr("stroke-width", d => Math.max(1, d.value * 4));

  const node = g.selectAll(".node")
    .data(nodes).enter().append("g").attr("class", "node")
    .call(d3.drag()
      .on("start", (e, d) => { if (!e.active) sim.alphaTarget(0.3).restart(); d.fx = d.x; d.fy = d.y; })
      .on("drag", (e, d) => { d.fx = e.x; d.fy = e.y; })
      .on("end", (e, d) => { if (!e.active) sim.alphaTarget(0); d.fx = null; d.fy = null; }));

  node.append("circle").attr("r", 6).attr("fill", "var(--accent)").attr("opacity", 0.9);
  node.append("text").text(d => d.id)
    .attr("x", 9).attr("dy", "0.35em")
    .style("font-size", "10px").style("fill", "var(--text-2)");

  sim.on("tick", () => {
    link.attr("x1", d => d.source.x).attr("y1", d => d.source.y)
        .attr("x2", d => d.target.x).attr("y2", d => d.target.y);
    node.attr("transform", d => `translate(${d.x},${d.y})`);
  });
}
