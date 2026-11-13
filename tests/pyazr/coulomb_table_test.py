#!/usr/bin/env python3
"""The Coulomb table's fetch and build (pyazr/coulomb_table.py; docs/upgrade/pyazr3.md section 2.9).

  1. build() writes a small table (l = 0, 5 nodes per decade) that the engine loads (coulomb_table_info);
  2. fetch() from a file:// URL and from a local HTTP server (http.server on 127.0.0.1) with the right SHA-256: the
     same bytes, the same content identity in the engine, no .part left;
  3. a corrupted copy (one byte changed) with the published SHA-256, by file:// and by HTTP: ChecksumMismatch, and
     neither the file nor its .part exists afterwards;
  4. an existing file is kept without overwrite and replaced with it; no URL or no hash known: ValueError naming
     build(); where() says where and whether.
Nothing is downloaded from outside this machine.  Scratch in the working directory (ctest: the build's tests/pyazr).
Exits 77 (skipped) without the compiled engine module.
"""
import functools
import http.server
import os
import shutil
import sys
import tempfile
import threading

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, ROOT)
try:
    from pyazr import _azure2  # noqa: F401
    from pyazr import coulomb_table
    from pyazr.described import coulomb_table as table_in_use
except ImportError as e:
    print(f"skip: {e}")
    sys.exit(77)

failures = []
# The local server must not be reached through a proxy of the environment.
os.environ["no_proxy"] = os.environ["NO_PROXY"] = "127.0.0.1,localhost"


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {str(detail)[:600]}"))
    if not ok:
        failures.append(name)


def refused(call, kind):
    try:
        call()
    except kind as e:
        return str(e)
    except Exception as e:  # a wrong class is a failure, shown with its text
        return f"WRONG CLASS {type(e).__name__}: {e}"
    return None


def main(scratch):
    built = os.path.join(scratch, "served", "table.bin")
    info = coulomb_table.build(built, threads=1, l_max=0, per_decade=5)
    check("build: the file is there", info["exists"] and info["path"] == built and info["bytes"] > 0, info)
    engine = table_in_use(coulomb_table=built)
    check("build: the engine loads it", engine["loaded"] and engine["identity"], engine)
    digest = coulomb_table.file_sha256(built)
    bad = os.path.join(scratch, "served", "corrupted.bin")
    shutil.copyfile(built, bad)
    with open(bad, "r+b") as f:
        f.seek(os.path.getsize(bad) // 2)
        byte = f.read(1)
        f.seek(-1, os.SEEK_CUR)
        f.write(bytes([byte[0] ^ 0xFF]))

    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=os.path.dirname(built))
    handler.log_message = lambda *a: None
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    http_base = f"http://127.0.0.1:{server.server_address[1]}/"
    try:
        for how, url, bad_url in (("file://", "file://" + built, "file://" + bad),
                                  ("HTTP", http_base + "table.bin", http_base + "corrupted.bin")):
            dest = os.path.join(scratch, how.strip(":/").lower(), "table.bin")
            got = coulomb_table.fetch(dest, url=url, sha256=digest)
            check(f"fetch by {how}: fetched, the same bytes, no .part", got["fetched"] and got["exists"]
                  and coulomb_table.file_sha256(dest) == digest and not os.path.exists(dest + ".part"), got)
            loaded = table_in_use(coulomb_table=dest)
            check(f"fetch by {how}: the engine loads it, the same content identity",
                  loaded["loaded"] and loaded["identity"] == engine["identity"], loaded)
            dest_bad = os.path.join(scratch, how.strip(":/").lower(), "refused.bin")
            text = refused(lambda: coulomb_table.fetch(dest_bad, url=bad_url, sha256=digest),
                           coulomb_table.ChecksumMismatch)
            check(f"fetch by {how} of a corrupted copy: refused (ChecksumMismatch, an OSError)",
                  text is not None and not text.startswith("WRONG") and digest in text, text)
            check(f"fetch by {how} of a corrupted copy: removed (neither the file nor .part)",
                  not os.path.exists(dest_bad) and not os.path.exists(dest_bad + ".part"),
                  os.listdir(os.path.dirname(dest_bad)))
    finally:
        server.shutdown()
        server.server_close()

    kept = os.path.join(scratch, "kept.bin")
    with open(kept, "wb") as f:
        f.write(b"old")
    got = coulomb_table.fetch(kept, url="file://" + built, sha256=digest)
    check("fetch onto an existing file: kept (fetched False)", not got["fetched"] and open(kept, "rb").read() == b"old",
          got)
    got = coulomb_table.fetch(kept, url="file://" + built, sha256=digest, overwrite=True)
    check("fetch with overwrite: replaced", got["fetched"] and coulomb_table.file_sha256(kept) == digest, got)
    text = refused(lambda: coulomb_table.fetch(os.path.join(scratch, "x.bin"), url="", sha256=digest), ValueError)
    check("fetch without a URL: ValueError naming build()", text is not None and "build" in text, text)
    text = refused(lambda: coulomb_table.fetch(os.path.join(scratch, "x.bin"), url="file://" + built, sha256=""),
                   ValueError)
    check("fetch without a hash: ValueError", text is not None and not text.startswith("WRONG"), text)
    check("TABLE_URL and TABLE_SHA256 are set together", bool(coulomb_table.TABLE_URL) ==
          bool(coulomb_table.TABLE_SHA256), (coulomb_table.TABLE_URL, coulomb_table.TABLE_SHA256))

    w = coulomb_table.where(os.path.join(scratch, "absent.bin"))
    check("where(FILE): the path, absent, source 'argument'", w["path"].endswith("absent.bin") and not w["exists"]
          and w["source"] == "argument", w)
    w = coulomb_table.where()
    check("where(): the engine's default path", w["path"] == coulomb_table.default_path()
          and w["source"] in ("AZURE2_COULOMB_TABLE", "default"), w)


if __name__ == "__main__":
    with tempfile.TemporaryDirectory(prefix="coulomb_table_test.", dir=os.getcwd()) as scratch:
        main(scratch)
    if failures:
        print(f"{len(failures)} failure(s)")
        sys.exit(1)
    print("all checks passed")
