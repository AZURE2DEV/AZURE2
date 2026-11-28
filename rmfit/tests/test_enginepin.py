"""Engine pin (enginepin.py, 2026-11-07): fingerprint, mismatch refusal, job-script guard.  Seconds, no engine session."""
import json
import os
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from rmfit import enginepin  # noqa: E402
from rmfit.campaign import Campaign  # noqa: E402


def test_fingerprint():
    fp = enginepin.fingerprint()
    assert fp["so_md5"] and len(fp["so_md5"]) == 32, fp
    assert fp["so"].startswith("_azure2"), fp
    # the commit is read with git: a runner without git on PATH (the Windows MSYS2 shell) has neither
    assert fp["git_head"] or fp["install_marker"] or fp.get("manifest_commit") or not shutil.which("git"), fp
    print("fingerprint:", enginepin.short(fp))
    return fp


def test_check(fp):
    same = dict(fp)
    cur, msg = enginepin.check(same, strict=True)
    assert msg is None and cur["so_md5"] == fp["so_md5"]
    cur, msg = enginepin.check(None, strict=True)
    assert msg is None
    other = dict(fp, so_md5="0" * 32, azure2_root="/nowhere/AZURE2_frozen_deadbeef")
    os.environ.pop("RMFIT_ALLOW_ENGINE_CHANGE", None)
    try:
        enginepin.check(other, strict=True)
        raise AssertionError("strict mismatch did not raise")
    except enginepin.EnginePinMismatch as e:
        assert "deadbeef" in str(e) or "/nowhere" in str(e)
    logged = []
    cur, msg = enginepin.check(other, strict=False, log=logged.append)
    assert msg and logged and logged[0].startswith("WARNING")
    os.environ["RMFIT_ALLOW_ENGINE_CHANGE"] = "1"
    cur, msg = enginepin.check(other, strict=True, log=logged.append)
    assert msg and "allowed" in logged[-1]
    os.environ.pop("RMFIT_ALLOW_ENGINE_CHANGE")
    print("check: ok")


def test_campaign_refuses(fp):
    d = tempfile.mkdtemp(prefix="rmfit_pin_")
    try:
        cfg = {"base_azr": "x.azr", "flags": {}, "nshards": 2, "policy": {}, "search": {}, "objective": {},
               "engine": dict(fp, so_md5="f" * 32, azure2_root="/nowhere/AZURE2_frozen_deadbeef")}
        with open(os.path.join(d, "campaign.json"), "w") as f:
            json.dump(cfg, f)
        try:
            Campaign(d)
            raise AssertionError("Campaign opened on a mismatched engine")
        except enginepin.EnginePinMismatch:
            pass
        c = Campaign(d, enforce_engine=False)          # read-only commands: warn only
        assert c.engine["so_md5"] == fp["so_md5"]
        assert "engine pin mismatch" in open(os.path.join(d, "campaign.log")).read()
        c.close()
        cfg["engine"] = fp                                # matching pin opens
        with open(os.path.join(d, "campaign.json"), "w") as f:
            json.dump(cfg, f)
        c = Campaign(d)
        js = c.write_job_script(overwrite=True)
        txt = open(js).read()
        assert "{engine_guard}" not in txt and "{pythonpath}" not in txt
        assert f"setenv AZURE2_ROOT {fp['azure2_root']}" in txt and fp["so_md5"] in txt
        assert "ENGINE PIN FAILED" in txt and "PYTHONUNBUFFERED" in txt
        c.close()
        cfg.pop("engine")                                  # pre-2026-11-07 campaign: never refused
        with open(os.path.join(d, "campaign.json"), "w") as f:
            json.dump(cfg, f)
        c = Campaign(d)
        txt = open(c.write_job_script(overwrite=True)).read()
        assert "no engine pin recorded" in txt
        c.close()
        print("campaign refusal / job script: ok")
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_repin_job_scripts(fp):
    """A re-pin rewrites the generated guard in every job script and keeps the rest; a hand-written mention is only listed."""
    d = tempfile.mkdtemp(prefix="rmfit_repin_")
    try:
        old = dict(fp, so_md5="a" * 32, azure2_root="/nowhere/AZURE2_frozen_old")
        new = dict(fp)
        csh = "#!/bin/csh\n#$ -N x\n" + enginepin.guard_lines_csh(old) + "python3 -m rmfit.cli . drive --max-minutes 123\n"
        sh = "#!/bin/bash\n" + enginepin.guard_lines_sh(old) + "python3 -m rmfit.cli . polish\n"
        hand = "#!/bin/csh\nif (`md5sum x` != " + old["so_md5"] + ") exit 1\n"
        for name, text in (("run_crc_rmfit", csh), ("run_slurm_drive", sh), ("run_hand", hand),
                           ("campaign.log", "pinned " + old["so_md5"] + "\n")):
            with open(os.path.join(d, name), "w") as f:
                f.write(text)
        os.chmod(os.path.join(d, "run_crc_rmfit"), 0o755)
        logged = []
        rw, warned = enginepin.repin_job_scripts(d, old, new, log=logged.append)
        assert rw == ["run_crc_rmfit", "run_slurm_drive"] and warned == ["run_hand"], (rw, warned)
        t = open(os.path.join(d, "run_crc_rmfit")).read()
        assert enginepin.guard_lines_csh(new) in t and old["so_md5"] not in t and "--max-minutes 123" in t
        # the executable bit: a Windows file system has none (st_mode there comes from the file's extension)
        assert os.name == "nt" or os.stat(os.path.join(d, "run_crc_rmfit")).st_mode & 0o111
        t = open(os.path.join(d, "run_slurm_drive")).read()
        assert enginepin.guard_lines_sh(new) in t and old["so_md5"] not in t
        assert old["so_md5"] in open(os.path.join(d, "run_hand")).read()
        assert any("WARNING" in m and "run_hand" in m for m in logged), logged
        assert enginepin.repin_job_scripts(d, new, new) == ([], [])          # same pin: nothing to do
        print("re-pin of job scripts: ok")
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_engine_dir():
    """A build with -DPYAZR_OUTPUT_DIR puts _azure2 outside pyazr/ and the tests set PYAZR_ENGINE_DIR (2026-11-19):
    the pin finds the module there, records engine_dir, and the job-script guards check that file.  Engine-free."""
    d = tempfile.mkdtemp(prefix="rmfit_enginedir_")
    old_env = os.environ.get("PYAZR_ENGINE_DIR")
    try:
        root, edir = os.path.join(d, "root"), os.path.join(d, "engine")
        os.makedirs(os.path.join(root, "pyazr")); os.makedirs(edir)
        with open(os.path.join(root, "pyazr", "__init__.py"), "w") as fh:
            fh.write('__version__ = "9.9.9"\n')
        with open(os.path.join(edir, "_azure2.cpython-312-x86_64-linux-gnu.so"), "wb") as fh:
            fh.write(b"not an engine")
        os.environ.pop("PYAZR_ENGINE_DIR", None)
        fp0 = enginepin.fingerprint(root)
        assert fp0["so"] is None and "engine_dir" not in fp0, fp0          # nothing in pyazr/, no engine dir: no module
        os.environ["PYAZR_ENGINE_DIR"] = edir
        fp1 = enginepin.fingerprint(root)
        assert fp1["so"] == "_azure2.cpython-312-x86_64-linux-gnu.so" and fp1["engine_dir"] == edir, fp1
        assert fp1["so_md5"] == enginepin._md5(os.path.join(edir, fp1["so"])), fp1
        assert f'export PYAZR_ENGINE_DIR="{edir}"' in enginepin.guard_lines_sh(fp1)
        assert f'md5sum "$PYAZR_ENGINE_DIR/{fp1["so"]}"' in enginepin.guard_lines_sh(fp1)
        assert f"setenv PYAZR_ENGINE_DIR {edir}" in enginepin.guard_lines_csh(fp1)
        # a module in pyazr/ wins over the engine dir, and then nothing extra is recorded or written
        shutil.copy(os.path.join(edir, fp1["so"]), os.path.join(root, "pyazr", fp1["so"]))
        fp2 = enginepin.fingerprint(root)
        assert "engine_dir" not in fp2 and "PYAZR_ENGINE_DIR" not in enginepin.guard_lines_sh(fp2), fp2
        print("engine module in PYAZR_ENGINE_DIR: ok")
    finally:
        if old_env is None:
            os.environ.pop("PYAZR_ENGINE_DIR", None)
        else:
            os.environ["PYAZR_ENGINE_DIR"] = old_env
        shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    fp = test_fingerprint()
    test_check(fp)
    test_campaign_refuses(fp)
    test_repin_job_scripts(fp)
    test_engine_dir()
    print("test_enginepin: PASS")
