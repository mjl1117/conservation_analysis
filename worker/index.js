// Cloudflare Worker: CORS proxy for EBI ClustalOmega REST API
const EBI_BASE = "https://www.ebi.ac.uk/Tools/services/rest/clustalo";

const CORS_HEADERS = {
  "Access-Control-Allow-Origin": "*",
  "Access-Control-Allow-Methods": "GET, POST, OPTIONS",
  "Access-Control-Allow-Headers": "Content-Type, Accept",
};

export default {
  async fetch(request) {
    // Handle CORS preflight
    if (request.method === "OPTIONS") {
      return new Response(null, { status: 204, headers: CORS_HEADERS });
    }

    const url = new URL(request.url);
    // Strip the worker's own path prefix, proxy the rest to EBI
    const ebiPath = url.pathname.replace(/^\/clustalo/, "");
    const ebiUrl = `${EBI_BASE}${ebiPath}${url.search}`;

    const ebiRequest = new Request(ebiUrl, {
      method: request.method,
      headers: { "Content-Type": request.headers.get("Content-Type") || "text/plain" },
      body: request.method !== "GET" ? request.body : undefined,
    });

    const ebiResponse = await fetch(ebiRequest);
    const body = await ebiResponse.text();

    return new Response(body, {
      status: ebiResponse.status,
      headers: {
        ...CORS_HEADERS,
        "Content-Type": ebiResponse.headers.get("Content-Type") || "text/plain",
      },
    });
  },
};
