// Shared by the node tests of the page's modules: where the modules are, and the engine's schema, which the runner
// (tests/pyazr/ui/node_test.py) writes to a file named by AZURE2_UI_SCHEMA before it starts node.
import { readFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

export const STATIC = join(dirname(fileURLToPath(import.meta.url)), "..", "..", "..", "..", "pyazr", "ui", "static");
export const mod = (name) => import(join(STATIC, name));

/** The engine's schema (pyazr.schema()) when the runner gave it, else null (the tests that need it skip). */
export function engineSchema() {
  const path = process.env.AZURE2_UI_SCHEMA;
  return path ? JSON.parse(readFileSync(path, "utf8")) : null;
}
