import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [react()],
  base: "./",  // relative paths for GitHub Pages subdirectory
  optimizeDeps: {
    exclude: ["3Dmol"],  // loaded via CDN script tag
  },
  build: {
    rollupOptions: {
      external: ["3Dmol"],  // prevent Rollup from trying to bundle it
    },
  },
  test: {
    environment: "jsdom",
    setupFiles: ["./src/tests/setup.js"],
  },
});
