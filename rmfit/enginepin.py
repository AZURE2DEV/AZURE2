"""Engine pin: which AZURE2 engine a campaign was started on, and a refusal to run on another.

Why (retrospective 2026-11-07, DESIGN.md): the shared ``~/AZURE2/pyazr`` .so was
replaced under running campaigns many times in October 2026 ("replaced three times on 10-21"), one fit
needed eight frozen builds in eight days, and campaigns compared objectives across engine changes without
noticing.  The practice that grew up by hand -- a frozen clone named by commit, ``setenv AZURE2_ROOT`` and an
md5 guard in every job script (an archived campaign, an archived campaign
run_crc_merged``) -- is made a mechanism here:

* ``fingerprint()`` describes the engine THIS process imports (``rmfit.AZURE2_ROOT`` -> ``pyazr``): the
  ``_azure2`` .so md5, the pyazr version, the git commit of the root when it is a clone, the newest
  ``.keep_*_pre-<hash>`` install marker otherwise, and the ``AZURE2`` CLI on the PATH (md5).
* ``rmfit init`` records it as ``campaign.json["engine"]``.
* ``Campaign(...)`` compares the recorded .so md5 with the imported one and raises ``EnginePinMismatch`` for
  every command that evaluates the model.  Read-only commands (status, report, log) only warn.  Set
  ``RMFIT_ALLOW_ENGINE_CHANGE=1`` to run anyway; the mismatch is then written to the campaign readme, because
  objectives before and after are not the same number.
* The job template (``templates/run_crc_rmfit.csh``) sets ``AZURE2_ROOT`` from the pin and checks the md5
  before Python starts.

A campaign without an ``engine`` entry (made before 2026-11-07) is never refused.  ``pyazr`` exports no build id
yet (``pyazr.__version__`` only); the md5 of the compiled module is the honest identity until it does.
"""
import hashlib
import glob
import os
import shutil
import subprocess
import sys


class EnginePinMismatch(RuntimeError):
    pass


def _md5(path, limit=None):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(1 << 20)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def _git_head(root):
    try:
        out = subprocess.run(["git", "-C", root, "rev-parse", "HEAD"], capture_output=True, text=True, timeout=20)
    except Exception:
        return None
    return out.stdout.strip() if out.returncode == 0 and out.stdout.strip() else None


def _keep_marker_hash(pyazr_dir):
    """The newest ``*.keep_<date>_pre-<hash>`` marker names the install that replaced it, i.e. the
    installed build when the root is the shared, un-committed checkout (``~/AZURE2``)."""
    marks = glob.glob(os.path.join(pyazr_dir, "*.keep_*_pre-*"))
    if not marks:
        return None
    newest = max(marks, key=os.path.getmtime)
    return os.path.basename(newest).rsplit("_pre-", 1)[-1]


def so_path(pyazr_dir):
    # the engine module: _azure2*.so, on Windows _azure2*.pyd
    cands = sorted(glob.glob(os.path.join(pyazr_dir, "_azure2*.so"))) or sorted(glob.glob(os.path.join(pyazr_dir, "_azure2*.pyd")))
    return cands[0] if cands else None


def engine_so(azure2_root):
    """The engine module the pyazr of ``azure2_root`` imports: ``pyazr/_azure2*`` first, else the one in
    ``$PYAZR_ENGINE_DIR`` (a build configured with -DPYAZR_OUTPUT_DIR writes it there and pyazr appends that directory
    to its package path, feature/upgrade 2026-11-17)."""
    so = so_path(os.path.join(os.path.abspath(azure2_root), "pyazr"))
    if not so and os.environ.get("PYAZR_ENGINE_DIR"):
        so = so_path(os.environ["PYAZR_ENGINE_DIR"])
    return so


def keeps_close_levels(azure2_root=None):
    """Does this engine keep two same-J^pi levels with different level IDs apart when they are within 1 keV?

    True since dev e239a12 (2026-10-24, CNuc::Fill: a line joins a level only with the same ID within 1 keV, or at the
    identical energy).  Older builds merge any pair within 1 keV on read, so only there a fit holding such a pair cannot
    be written exactly.  A campaign can be pinned to either, and a commit comparison needs a clone, so the test is whether
    the .so carries that branch's warning text.  No .so found: False, the old rule."""
    if azure2_root is None:
        from . import AZURE2_ROOT
        azure2_root = AZURE2_ROOT
    so = engine_so(azure2_root)
    if not so:
        return False
    key = (so, os.path.getmtime(so))
    if key not in _keeps_close_levels:
        with open(so, "rb") as fh:
            _keeps_close_levels[key] = b"they are kept as two levels" in fh.read()
    return _keeps_close_levels[key]


_keeps_close_levels = {}


def fingerprint(azure2_root=None, cli=None):
    """Describe the engine this process would import.  ``azure2_root`` defaults to ``rmfit.AZURE2_ROOT``."""
    if azure2_root is None:
        from . import AZURE2_ROOT
        azure2_root = AZURE2_ROOT
    root = os.path.abspath(azure2_root)
    pyazr_dir = os.path.join(root, "pyazr")
    so = engine_so(root)
    fp = {
        "azure2_root": root,
        "so": os.path.basename(so) if so else None,
        "so_md5": _md5(so) if so else None,
        "pyazr_version": None,
        "git_head": _git_head(root),
        "install_marker": _keep_marker_hash(pyazr_dir),
        "cli": None,
        "cli_md5": None,
        "python": sys.version.split()[0],
    }
    if so and os.path.dirname(os.path.abspath(so)) != pyazr_dir:
        fp["engine_dir"] = os.path.dirname(os.path.abspath(so))   # recorded only when outside pyazr/ (PYAZR_ENGINE_DIR)
    init_py = os.path.join(pyazr_dir, "__init__.py")
    if os.path.exists(init_py):
        with open(init_py) as fh:
            for line in fh:
                if line.startswith("__version__"):
                    fp["pyazr_version"] = line.split("=", 1)[1].strip().strip("'\"")
                    break
    # the engine's own command-line program before whatever AZURE2 is first on PATH (2026-11-26: a frozen engine of
    # AZURE2_engines/ has bin/AZURE2, a build tree build/src/AZURE2; PATH in a job usually still finds the old shared
    # ~/bin/AZURE2, which recorded the wrong binary).  Informational: check() compares the .so md5 only.
    if not cli:
        cli = next((p for p in (os.path.join(root, "bin", "AZURE2"), os.path.join(root, "build", "src", "AZURE2"))
                    if os.path.isfile(p)), None) or shutil.which("AZURE2")
    if cli and os.path.exists(cli):
        fp["cli"] = os.path.abspath(cli)
        fp["cli_md5"] = _md5(cli)
    # a frozen engine made by AZURE2_engines/freeze_engine.sh has no git; its MANIFEST names the commit
    man = os.path.join(root, "MANIFEST")
    if os.path.isfile(man):
        try:
            with open(man) as fh:
                for line in fh:
                    if line.startswith("commit "):
                        fp["manifest_commit"] = line.split()[1]
                        break
        except OSError:
            pass
    return fp


def short(fp):
    """One line for logs: root, .so md5 (8), commit/marker."""
    tag = (fp.get("git_head") or fp.get("manifest_commit") or fp.get("install_marker") or "?")[:7]
    return f"{fp.get('azure2_root')} .so {str(fp.get('so_md5'))[:8]} ({tag}, pyazr {fp.get('pyazr_version')})"


def check(recorded, strict=True, log=None):
    """Compare the recorded pin with the engine this process imports.

    Returns ``(current fingerprint, mismatch message or None)``.  On a mismatch of the .so md5: raise
    ``EnginePinMismatch`` when ``strict`` and ``RMFIT_ALLOW_ENGINE_CHANGE`` is unset; otherwise call ``log`` with
    the warning and return the message so the caller can record it."""
    cur = fingerprint()
    if not recorded or not recorded.get("so_md5"):
        return cur, None
    if recorded["so_md5"] == cur["so_md5"]:
        return cur, None
    msg = (f"engine pin mismatch: campaign pinned {short(recorded)}; this process imports {short(cur)}. "
           "Objectives are not comparable across engines.  Run from the pinned root (the job script sets "
           "AZURE2_ROOT from campaign.json), or set RMFIT_ALLOW_ENGINE_CHANGE=1 to continue and record the change.")
    allowed = os.environ.get("RMFIT_ALLOW_ENGINE_CHANGE", "") not in ("", "0")
    if strict and not allowed:
        raise EnginePinMismatch(msg)
    if log:
        log("WARNING " + msg + ("  (allowed by RMFIT_ALLOW_ENGINE_CHANGE)" if allowed else ""))
    return cur, msg


def guard_lines_sh(recorded):
    """The same guard for a bash/sh job script (the Slurm templates)."""
    if not recorded or not recorded.get("so_md5"):
        return ("# no engine pin recorded in campaign.json: the shared checkout is used\n"
                'export AZURE2_ROOT="${AZURE2_ROOT:-$HOME/AZURE2}"\n')
    root = recorded["azure2_root"]
    so = recorded.get("so") or "_azure2.cpython-312-x86_64-linux-gnu.so"
    md5 = recorded["so_md5"]   # (no quotes inside the f-strings: Python < 3.12 rejects them)
    edir = recorded.get("engine_dir")
    return (f'export AZURE2_ROOT="{root}"\n'
            f'export PYAZR_ROOT="$AZURE2_ROOT"\n'
            + (f'export PYAZR_ENGINE_DIR="{edir}"\n'
               f'm=$(md5sum "$PYAZR_ENGINE_DIR/{so}" | cut -c1-32)\n' if edir else
               f'm=$(md5sum "$AZURE2_ROOT/pyazr/{so}" | cut -c1-32)\n') +
            f'if [ "$m" != "{md5}" ]; then\n'
            f'  echo "$(date) ENGINE PIN FAILED: $m != {md5} ({short(recorded)})"\n'
            f'  exit 1\n'
            f'fi\n')


def guard_lines_csh(recorded):
    """csh lines for a job script: export the pinned root, verify the .so md5 before Python runs."""
    if not recorded or not recorded.get("so_md5"):
        return ("# no engine pin recorded in campaign.json (campaign made before 2026-11-07): the shared checkout is used\n"
                "if (! $?AZURE2_ROOT) setenv AZURE2_ROOT ${HOME}/AZURE2\n")
    root = recorded["azure2_root"]
    so = recorded.get("so") or "_azure2.cpython-312-x86_64-linux-gnu.so"
    edir = recorded.get("engine_dir")
    return (f"setenv AZURE2_ROOT {root}\n"
            f"setenv PYAZR_ROOT $AZURE2_ROOT\n"
            + (f"setenv PYAZR_ENGINE_DIR {edir}\n"
               f"set m = `md5sum ${{PYAZR_ENGINE_DIR}}/{so} | cut -c1-32`\n" if edir else
               f"set m = `md5sum ${{AZURE2_ROOT}}/pyazr/{so} | cut -c1-32`\n") +
            f"if (\"$m\" != \"{recorded['so_md5']}\") then\n"
            f"  echo \"`date` ENGINE PIN FAILED: $m != {recorded['so_md5']} ({short(recorded)})\"\n"
            f"  exit 1\n"
            f"endif\n")


def repin_job_scripts(directory, old, new, log=None):
    """After a re-pin (``init`` with ``RMFIT_ALLOW_ENGINE_CHANGE=1``), carry the new pin into the campaign's job scripts.

    Every file in ``directory`` that holds the old engine guard exactly as ``guard_lines_csh`` / ``guard_lines_sh`` wrote it
    gets the new guard in its place; the rest of the script (a drive's wall time, hand edits) is kept.  A file that names
    the old .so md5 without that block (a hand-written script) is listed with a WARNING, not edited.  Returns
    ``(rewritten, warned)``, lists of file names.  2026-11-16: before this, a re-pin left every job script checking the
    old md5, so each job stopped with ENGINE PIN FAILED."""
    if not old or not old.get("so_md5") or not new or old.get("so_md5") == new.get("so_md5"):
        return [], []
    pairs = [(guard_lines_csh(old), guard_lines_csh(new)), (guard_lines_sh(old), guard_lines_sh(new))]
    rewritten, warned = [], []
    for name in sorted(os.listdir(directory)):
        path = os.path.join(directory, name)
        if not os.path.isfile(path) or os.path.getsize(path) > 1_000_000:
            continue
        try:
            with open(path) as f:
                text = f.read()
        except (UnicodeDecodeError, OSError):
            continue
        if old["so_md5"] not in text:
            continue
        new_text = text
        for a, b in pairs:
            new_text = new_text.replace(a, b)
        if new_text != text and old["so_md5"] not in new_text:
            mode = os.stat(path).st_mode
            with open(path, "w") as f:
                f.write(new_text)
            os.chmod(path, mode)
            rewritten.append(name)
        elif name not in ("campaign.json", "campaign.log", "init.log", "LOG.jsonl", "readme", "STATUS.md") \
                and not name.endswith((".log", ".jsonl", ".md", ".json", ".sqlite")):
            warned.append(name)
    if log:
        if rewritten:
            log(f"re-pin: engine guard updated in {', '.join(rewritten)}")
        if warned:
            log(f"WARNING re-pin: {', '.join(warned)} name the old .so md5 {old['so_md5'][:8]} outside the generated "
                f"guard block; not edited -- update them by hand")
    return rewritten, warned
