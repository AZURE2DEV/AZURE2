"""The node tests of the page's modules (tests/pyazr/ui/js/*.test.mjs): the widget table and the schema-to-form
description of forms.js, sorting and filtering of tables.js, the revision and error handling of api.js, the plot
table of plots.js, the Fit tab's rows, job form and progress views of fitting.js, the Physics tab's rules of
physicsform.js.  The engine's schema, the server's job table and the function tables (the physics rows, the
session's rows, the bound of a request) are written to files for them (AZURE2_UI_SCHEMA, AZURE2_UI_JOBS,
AZURE2_UI_FUNCTIONS).

Needs ``node`` (version 18 or later: ``node --test``) on PATH or named by the environment variable NODE; without it,
or without the engine module, prints ``skip: <reason>`` and exits 77.  Nothing is installed: the tests use node's own
test runner and the page's modules as they are served.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

import uitest


def main():
    node = os.environ.get("NODE") or shutil.which("node")
    if not node:
        print("skip: node is not on PATH (set NODE=/path/to/node to run the page's module tests)")
        return 77
    version = subprocess.run([node, "--version"], capture_output=True, text=True).stdout.strip()
    if not version.startswith("v") or int(version[1:].split(".")[0]) < 18:
        print(f"skip: node {version} has no test runner (node --test needs 18 or later)")
        return 77
    reason = uitest.missing_requirements()
    if reason:
        print(f"skip: {reason}")
        return 77
    import pyazr
    from pyazr.ui.jobs import job_table
    with tempfile.TemporaryDirectory() as tmp:
        schema = os.path.join(tmp, "schema.json")
        with open(schema, "w") as f:
            json.dump(pyazr.schema(), f)
        jobs = os.path.join(tmp, "jobs.json")
        with open(jobs, "w") as f:
            json.dump(job_table(), f)
        tests = sorted(os.path.join(uitest.HERE, "js", n) for n in os.listdir(os.path.join(uitest.HERE, "js"))
                       if n.endswith(".test.mjs"))
        functions = os.path.join(tmp, "functions.json")
        from pyazr.ui.physics import MAX_ELEMENTS, session_rows
        with open(functions, "w") as f:
            json.dump({"tables": [list(pyazr.kinds()["physics"].values()), list(session_rows().values())],
                       "max_elements": MAX_ELEMENTS}, f)
        env = dict(os.environ, AZURE2_UI_SCHEMA=schema, AZURE2_UI_JOBS=jobs, AZURE2_UI_FUNCTIONS=functions)
        proc = subprocess.run([node, "--test", *tests], env=env, capture_output=True, text=True, timeout=300)
    out = proc.stdout + proc.stderr
    failed = [line for line in out.splitlines() if line.startswith("not ok")]
    summary = [line for line in out.splitlines() if line.startswith(("# tests", "# pass", "# fail"))]
    print("\n".join(summary))
    if proc.returncode != 0 or failed:
        print(out[-6000:])
        print(f"FAILED (node {version}: {len(failed)} failed)")
        return 1
    print(f"ok (node {version})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
