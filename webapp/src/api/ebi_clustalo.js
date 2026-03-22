const PROXY_URL = import.meta.env.VITE_CLUSTALO_PROXY_URL || "http://localhost:8787";
const POLL_INTERVAL_MS = 10_000;
const MAX_WAIT_MS = 10 * 60 * 1000; // 10 minutes

export async function submitAlignmentJob(fastaText) {
  const body = new URLSearchParams({
    sequence: fastaText,
    outfmt: "fa",
    stype: "protein",
    email: "pdu-conservation@example.com",
  });
  const res = await fetch(`${PROXY_URL}/clustalo/run`, {
    method: "POST",
    headers: { "Content-Type": "application/x-www-form-urlencoded" },
    body: body.toString(),
  });
  if (!res.ok) throw new Error(`ClustalOmega job submission failed: ${res.status}`);
  const jobId = (await res.text()).trim();
  return jobId;
}

export async function pollUntilDone(jobId, onProgress) {
  const start = Date.now();
  while (Date.now() - start < MAX_WAIT_MS) {
    await sleep(POLL_INTERVAL_MS);
    const res = await fetch(`${PROXY_URL}/clustalo/status/${jobId}`);
    if (!res.ok) throw new Error(`Status poll failed: ${res.status}`);
    const status = (await res.text()).trim();
    const elapsed = Date.now() - start;
    onProgress?.({ status, elapsed, maxWait: MAX_WAIT_MS });
    if (status === "FINISHED") return;
    if (status === "ERROR" || status === "FAILURE") {
      throw new Error(`ClustalOmega job failed with status: ${status}`);
    }
  }
  throw new Error(
    "ClustalOmega job timed out after 10 minutes. For large datasets, use the CLI tool instead."
  );
}

export async function fetchAlignmentResult(jobId) {
  const res = await fetch(`${PROXY_URL}/clustalo/result/${jobId}/aln-fasta`);
  if (!res.ok) throw new Error(`Result fetch failed: ${res.status}`);
  return res.text();
}

export async function alignSequences(sequences, onProgress) {
  const fastaText = sequences.map(({ header, seq }) => `>${header}\n${seq}`).join("\n");
  const jobId = await submitAlignmentJob(fastaText);
  await pollUntilDone(jobId, onProgress);
  return fetchAlignmentResult(jobId);
}

function sleep(ms) { return new Promise(resolve => setTimeout(resolve, ms)); }
