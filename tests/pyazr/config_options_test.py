#!/usr/bin/env python3
"""The options of the calculation as ``<config>`` keys, through pyazr (docs/upgrade/project_format3.md section 3.12).

1. key = keyword: for every option and both of its values, a session of tests/upgrade/T1 with the key in ``<config>``
   (written by ``Project.set_config``) computes what a session with the keyword computes: the start vector, every
   residual and the objective, exactly.  A value that moves T1's numbers is shown to.
2. Precedence: a keyword given wins over the key (even with the default's value), the key over the default;
   ``option_sources`` and ``describe()["option_sources"]`` say "keyword", "project" or "default"; ``options`` holds the
   values in effect, and passed back as keywords they rebuild the same session.
3. Refusals: ``set_config`` of a bad value (ValueError) and of each contradiction (the engine's code, the project
   unchanged); a keyword against a key (AZURE2Error naming both and their sources); contradicting keywords without a
   key (``use_rmc=True`` with Brune at its default, ``use_park=True, use_brune=False``, ``ec_convection_current=False``
   with the long-wavelength form at its default), refused as the program refuses the flags.
4. ``set_config`` round trip (bool -> "1"/"0", ``None`` removes, the text back byte for byte); ``Project.new`` writes no
   option key T1 does not have (it keeps T1's ``<config>``); ``save_fit`` writes none it did not find.
5. ``pyazr.kinds()["options"]``: every option of the calculation has its ``<config>`` key (the same word) and every
   other has a reason (the one table).

Needs the compiled engine; without it exits 77 (ctest: Skipped).  Threads: as the caller's environment sets them.

Run from anywhere:  python3 tests/pyazr/config_options_test.py
"""
import os
import shutil
import sys
import tempfile
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(("  ok    " if ok else "  FAIL  ") + name + ("" if ok else f"  ({detail})"))
    if not ok:
        failures.append(name)


try:
    import numpy as np
    import pyazr
    from pyazr import _azure2  # noqa: F401
except Exception as err:  # pragma: no cover
    print(f"skip: the engine module is not available ({type(err).__name__}: {err})")
    sys.exit(77)

warnings.simplefilter("ignore", pyazr.EngineWarning)

# Every option, both values; options that need a partner carry it in both sessions (as the program's check).
CASES = [("use_brune", True, None, False), ("use_brune", False, None, True),
         ("use_park", True, None, True), ("use_park", False, None, False),
         ("transform", True, None, False), ("transform", False, None, True),
         ("ignore_externals", True, None, False), ("ignore_externals", False, None, False),
         ("use_long_wavelength", True, None, False), ("use_long_wavelength", False, None, True),
         ("ec_convection_current", True, None, False),
         ("ec_convection_current", False, ("use_long_wavelength", False), True),
         ("use_rmc", False, None, False), ("use_rmc", True, ("use_brune", False), True),
         ("use_gsl_coul", False, None, False), ("use_gsl_coul", True, None, False),
         ("use_coulomb_table", True, None, False), ("use_coulomb_table", False, None, False),
         ("use_ec_table", True, None, False), ("use_ec_table", False, None, False),
         ("replace_grids", True, None, False), ("replace_grids", False, None, False)]


def numbers(s):
    """What a session computes at its start: the start vector, the residuals, the objective."""
    x = np.asarray(s.params_rwa, float)
    return x, np.asarray(s.residuals(x), float), np.asarray([s.objective(x)])


def equal(a, b):
    return all(np.array_equal(p, q) for p, q in zip(a, b))


tmp = tempfile.mkdtemp(prefix="pyazr_config_options_")
try:
    t1 = os.path.join(tmp, "T1")
    shutil.copytree(os.path.join(ROOT, "tests", "upgrade", "T1"), t1, ignore=shutil.ignore_patterns("output", "checks"))
    os.makedirs(os.path.join(t1, "output"), exist_ok=True)
    base_path = os.path.join(t1, "T1.azr")
    base_text = open(base_path, "rb").read()

    def variant(name, **keys):
        """T1 with these <config> keys (Project.set_config), saved beside T1; the text must change."""
        p = pyazr.Project.open(base_path)
        p.set_config(**keys)
        path = os.path.join(t1, name + ".azr")
        p.save(path)
        if open(path, "rb").read() == base_text:
            raise AssertionError(f"{name}: set_config({keys}) changed nothing")
        return path

    with pyazr.azure2(base_path) as s:
        base = numbers(s)
        check("T1 alone: every option from its default", set(s.option_sources.values()) == {"default"},
              s.option_sources)

    # 1. key = keyword
    for i, (key, value, partner, moves) in enumerate(CASES):
        keys = {key: value}
        if partner:
            keys[partner[0]] = partner[1]
        path = variant(f"case{i}", **keys)
        label = f"{key}={int(value)}" + (f" ({partner[0]}={int(partner[1])})" if partner else "")
        outcome = {}
        for how, call in (("key", lambda: pyazr.azure2(path)), ("keyword", lambda: pyazr.azure2(base_path, **keys))):
            try:
                with call() as s:
                    outcome[how] = (numbers(s), all(s.option_sources[k] == ("project" if how == "key" else "keyword")
                                                    and s.options[k] == v for k, v in keys.items()))
            except pyazr.AZURE2Error as err:  # use_rmc on T1: the session refuses the capture (as before, both ways)
                outcome[how] = ("refused", str(err).replace(path, "").replace(base_path, ""))
        if outcome["key"][0] == "refused" or outcome["keyword"][0] == "refused":
            check(f"{label}: key = keyword (both refused alike)", outcome["key"] == outcome["keyword"], outcome)
            continue
        check(f"{label}: key = keyword (start vector, residuals, objective identical)",
              equal(outcome["key"][0], outcome["keyword"][0]))
        check(f"{label}: sources 'project' and 'keyword', values in effect", outcome["key"][1] and outcome["keyword"][1])
        if moves:
            check(f"{label}: T1's numbers move", not equal(outcome["key"][0], base))

    # 2. precedence
    brune0 = variant("brune0", use_brune=False)
    with pyazr.azure2(brune0, use_brune=True) as s:
        check("keyword use_brune=True (the default's value) over the key use_brune=0: T1 alone's numbers",
              equal(numbers(s), base) and s.option_sources["use_brune"] == "keyword")
    with pyazr.azure2(brune0) as s:
        options = dict(s.options)
        d = s.describe()
        check("describe(): option_sources", d["option_sources"]["use_brune"] == "project"
              and d["option_sources"]["transform"] == "default" and d["options"]["use_brune"] is False)
        mine = numbers(s)
    with pyazr.azure2(base_path, **options) as s:
        check("options passed back as keywords rebuild the same session (on T1 without the key)", equal(numbers(s), mine))
    check("options holds data_mode and coulomb_table besides the options of the calculation",
          {"data_mode", "coulomb_table", "use_brune", "replace_grids"} <= set(options))

    # 3. refusals
    p = pyazr.Project.open(base_path)
    try:
        p.set_config(use_brune="yes")
        check("set_config(use_brune='yes') refused", False)
    except ValueError as err:
        check("set_config(use_brune='yes') refused: ValueError with the values taken", "0, 1, true or false" in str(err),
              err)
    for keys, words in [({"use_park": True, "use_brune": False}, ["use_park=1", "use_brune=0"]),
                        ({"use_rmc": True}, ["use_rmc=1", "use_brune=1 (the default)"]),
                        ({"ec_convection_current": False}, ["ec_convection_current=0", "use_long_wavelength=1"])]:
        p = pyazr.Project.open(base_path)
        try:
            p.set_config(**keys)
            check(f"set_config({keys}) refused", False)
        except Exception as err:  # noqa: BLE001 - the class is pyazr's validation error
            text = str(err) + " ".join(str(getattr(d, "message", d)) for d in getattr(err, "diagnostics", []) or [])
            check(f"set_config({keys}) refused with E_CONFIG_OPTIONS_CONTRADICT, both keys named, project unchanged",
                  "E_CONFIG_OPTIONS_CONTRADICT" in text and all(w in text for w in words)
                  and p.text.encode("utf-8", "surrogateescape") == base_text, text[:300])
    brune1 = variant("brune1", use_brune=True)
    try:
        pyazr.azure2(brune1, use_rmc=True).close()
        check("use_rmc=True against the key use_brune=1 refused", False)
    except pyazr.AZURE2Error as err:
        check("use_rmc=True against the key use_brune=1: AZURE2Error naming both and their sources",
              "E_CONFIG_OPTIONS_CONTRADICT" in str(err) and "use_rmc=1 (keyword)" in str(err)
              and "use_brune=1 (this file)" in str(err), err)
    # Contradicting keywords, without a key: refused as the program refuses the flags (the same code, the same three
    # parts: what, allowed, fix), the fix in keywords.
    for keys, words in [({"use_rmc": True}, ["use_rmc=1 (keyword) with use_brune=1 (the default) cannot hold together",
                                             "; allowed: use_rmc=1 with use_brune=0, or use_rmc=0",
                                             "; fix: pass use_brune=False (or use_rmc=False in place of use_rmc=True)"]),
                        ({"use_park": True, "use_brune": False},
                         ["use_park=1 (keyword) with use_brune=0 (keyword) cannot hold together",
                          "; allowed: use_park=1 with use_brune=1, or use_park=0",
                          "; fix: use_brune=True in place of use_brune=False (or use_park=False in place of "
                          "use_park=True)"]),
                        ({"ec_convection_current": False},
                         ["ec_convection_current=0 (keyword) with use_long_wavelength=1 (the default)",
                          "; fix: pass use_long_wavelength=False"])]:
        try:
            pyazr.azure2(base_path, **keys).close()
            check(f"{keys} on T1 without keys refused", False)
        except pyazr.AZURE2Error as err:
            check(f"{keys} on T1 without keys: AZURE2Error E_CONFIG_OPTIONS_CONTRADICT, what, allowed, fix",
                  "E_CONFIG_OPTIONS_CONTRADICT" in str(err) and all(w in str(err) for w in words), err)
    for keys in ({"use_park": True}, {"use_brune": False}):  # an implied value counts as given: no contradiction
        try:
            with pyazr.azure2(base_path, **keys) as s:
                check(f"{keys} alone: accepted (the implied option given)", s.options["use_park"] == s.options["use_brune"])
        except pyazr.AZURE2Error as err:
            check(f"{keys} alone: accepted", False, err)

    # 4. round trip and writers
    p = pyazr.Project.open(base_path)
    p.set_config(use_brune=False, transform=True)
    check("set_config writes use_brune=0 and transform=1", p.config().get("use_brune") == "0"
          and p.config().get("transform") == "1")
    p.set_config(use_brune=None, transform=None)
    check("None removes the keys: the text back byte for byte", p.text.encode("utf-8", "surrogateescape") == base_text)
    option_keys = [k for k, o in pyazr.kinds()["options"].items() if o["key"]]
    new = pyazr.Project.new(os.path.join(t1, "new.azr"), nuclear_part=base_path)
    check("Project.new writes no option key", not set(option_keys) & set(new.config()), new.config())
    with pyazr.azure2(base_path) as s:
        out = os.path.join(t1, "fit.azr")
        s.save_fit(out)
    check("save_fit writes no option key it did not find",
          not set(option_keys) & set(pyazr.Project.open(out).config()))
    with pyazr.azure2(brune0) as s:
        out = os.path.join(t1, "fit_brune0.azr")
        # verify=False: T1's 15 MeV background levels do not read back in the formal basis (the same with the keyword
        # use_brune=False before the keys existed); what is checked here is the <config> written
        s.save_fit(out, verify=False)
    check("save_fit keeps the key it found", pyazr.Project.open(out).config().get("use_brune") == "0")

    # 5. the one table
    rows = pyazr.kinds()["options"]
    check("eleven options of the calculation, each its own key", sorted(k for k, o in rows.items() if o["key"]) ==
          sorted(o["key"] for o in rows.values() if o["key"]) and len(option_keys) == 11)
    check("each other option (not project state) has its reason",
          all(o["reason"] for o in rows.values() if not o["key"]) and
          {k for k, o in rows.items() if not o["key"]} == {"data_mode", "coulomb_table", "cwd"})
    positional, keys = _azure2.config_keys()
    check("the option keys are keys of <config> (bool)", all((k, "bool") in [(n, t) for n, t, _ in keys]
                                                               for k in option_keys))
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print(f"{len(failures)} failure(s)" if failures else "all checks passed")
sys.exit(1 if failures else 0)
