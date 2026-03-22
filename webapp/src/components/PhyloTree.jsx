// webapp/src/components/PhyloTree.jsx
import React, { useEffect, useRef, useState } from "react";
import * as d3 from "d3";
import { parseNewick } from "../analysis/phylo.js";

export default function PhyloTree({ data }) {
  const svgRef = useRef(null);
  const [layout, setLayout] = useState("rectangular"); // "rectangular" | "radial"
  const [search, setSearch] = useState("");

  useEffect(() => {
    if (!data?.newick || !svgRef.current) return;
    renderTree(svgRef.current, data, layout, search);
  }, [data, layout, search]);

  if (!data?.newick) return null;

  return (
    <div className="card">
      <div className="card-title">Phylogenetic Tree</div>
      <div style={{ display: "flex", gap: "0.75rem", alignItems: "center", marginBottom: "1rem", flexWrap: "wrap" }}>
        <div className="tab-bar" style={{ margin: 0 }}>
          <button className={`tab ${layout === "rectangular" ? "active" : ""}`} onClick={() => setLayout("rectangular")}>Rectangular</button>
          <button className={`tab ${layout === "radial" ? "active" : ""}`} onClick={() => setLayout("radial")}>Radial</button>
        </div>
        <input
          type="text" placeholder="Search sequence…" value={search}
          onChange={e => setSearch(e.target.value)}
          style={{ background: "var(--surface)", border: "1px solid var(--border)", borderRadius: "6px",
            color: "var(--text)", padding: "0.3rem 0.7rem", fontSize: "0.8rem", width: "180px" }}
        />
        <span style={{ fontSize: "0.75rem", color: "var(--text-2)" }}>
          {data.n_sequences} sequences · orange leaf = experimental PDB
        </span>
      </div>
      <svg ref={svgRef} style={{ width: "100%", height: layout === "radial" ? "600px" : "500px", overflow: "visible" }} />
    </div>
  );
}

function renderTree(svgEl, data, layout, search) {
  const root = d3.hierarchy(parseNewick(data.newick), d => d.children?.length ? d.children : null);
  const pdbAccessions = new Set(
    (data.sequences || []).filter(s => s.has_pdb).map(s => s.id)
  );

  const svg = d3.select(svgEl);
  svg.selectAll("*").remove();

  const W = svgEl.clientWidth || 800;
  const H = parseInt(svgEl.style.height) || 500;

  if (layout === "rectangular") {
    renderRectangular(svg, root, W, H, pdbAccessions, search);
  } else {
    renderRadial(svg, root, W, H, pdbAccessions, search);
  }
}

function renderRectangular(svg, root, W, H, pdbAccessions, search) {
  const margin = { top: 20, right: 180, bottom: 20, left: 20 };
  const w = W - margin.left - margin.right;
  const h = H - margin.top - margin.bottom;

  const cluster = d3.cluster().size([h, w]);
  cluster(root);

  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

  g.selectAll(".link")
    .data(root.links()).enter().append("path").attr("class", "link")
    .attr("fill", "none").attr("stroke", "var(--border)").attr("stroke-width", 1)
    .attr("d", d3.linkHorizontal().x(d => d.y).y(d => d.x));

  const leaf = g.selectAll(".leaf")
    .data(root.leaves()).enter().append("g").attr("class", "leaf")
    .attr("transform", d => `translate(${d.y},${d.x})`);

  leaf.append("circle").attr("r", 3)
    .attr("fill", d => pdbAccessions.has(d.data.name) ? "var(--accent-am)" : "var(--accent)");

  leaf.append("text")
    .attr("x", 6).attr("dy", "0.31em")
    .style("font-size", "10px").style("fill", d => {
      if (search && d.data.name.toLowerCase().includes(search.toLowerCase()))
        return "var(--accent-am)";
      return "var(--text-2)";
    })
    .text(d => d.data.name.slice(0, 20));
}

function renderRadial(svg, root, W, H, pdbAccessions, search) {
  const radius = Math.min(W, H) / 2 - 80;
  const cluster = d3.cluster().size([2 * Math.PI, radius]);
  cluster(root);

  const g = svg.append("g").attr("transform", `translate(${W/2},${H/2})`);

  g.selectAll(".link")
    .data(root.links()).enter().append("path").attr("class", "link")
    .attr("fill", "none").attr("stroke", "var(--border)").attr("stroke-width", 1)
    .attr("d", d3.linkRadial().angle(d => d.x).radius(d => d.y));

  const leaf = g.selectAll(".leaf")
    .data(root.leaves()).enter().append("g").attr("class", "leaf")
    .attr("transform", d => `rotate(${d.x * 180 / Math.PI - 90}) translate(${d.y},0)`);

  leaf.append("circle").attr("r", 3)
    .attr("fill", d => pdbAccessions.has(d.data.name) ? "var(--accent-am)" : "var(--accent)");

  leaf.append("text")
    .attr("x", d => d.x < Math.PI ? 6 : -6)
    .attr("text-anchor", d => d.x < Math.PI ? "start" : "end")
    .attr("transform", d => d.x >= Math.PI ? "rotate(180)" : null)
    .style("font-size", "9px").style("fill", "var(--text-2)")
    .text(d => d.data.name.slice(0, 16));
}
