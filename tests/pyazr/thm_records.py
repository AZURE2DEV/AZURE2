"""The THM part of a format-3 project text edited on its records (project_format3.md section 3.9), for the THM tests.

A data set joins an experiment by its own field `thm <experiment>` (`thm global`: none); the <thm> section holds the
global options and each experiment's keys, without a member list.  The same edits as tests/lib/thm_records.sh:

    thm_set(text, block, members)   text with its <thm> section replaced by BLOCK (appended when there is none; an
                                    empty BLOCK removes it) and the data sets MEMBERS (names) in the one experiment
                                    BLOCK defines
    thm_join(text, experiment, members)   the data sets MEMBERS name EXPERIMENT
    datasets(text)                  the names of the data sets, in the order of <data>
    no_fold(text)                   text without its folds: the records of <targets> and the beam and quadrature
                                    fields of its data sets and extrapolations (the THM projects' only target
                                    records are their beams)

On a pyazr.Project (the same edits through the engine, each validated):

    thm_datasets(project)           the names of the THM data sets (those with a `thm` field), in record order
    set_experiment(project, name, members, **keys)
                                    experiment[NAME] with exactly KEYS (keys it had and KEYS lacks are removed; a
                                    value may be a number, a bool (on/off), a (lo, hi) range or a list, written as
                                    the engine reads it), and the data sets MEMBERS (names, or 1-based positions
                                    among thm_datasets) naming it; data sets that named it and are not MEMBERS go
                                    back to `thm global`
    session_refusal(project, path)  PROJECT saved to PATH and a pyazr session made on it in a child process (the
                                    engine prints on the process's stdout): (refused, the engine's output).  The
                                    rules that need the pairs (a ps window with spectatorEnergy, a global optical
                                    potential's projectile and mass range) and entranceL=coherent with a theta
                                    window are refused here, not by the edit
    refusal(make, path)             the refusal of a variant: MAKE() returns the edited project, or raises
                                    ValueError (the edit refused: its message); else the engine's output when a
                                    session on it (saved to PATH) is refused; None when both take it

A name that is not a THM data set of the text raises ValueError: a misspelt name must not leave a data set out of
its experiment unnoticed.
"""
import os
import re
import subprocess
import sys


def datasets(text):
    out, inside = [], False
    for line in text.split("\n"):
        if line.startswith("<data>"):
            inside = True
        elif line.startswith("</data>"):
            inside = False
        elif inside and line.split()[:1] == ["dataset"]:
            out.append(line.split()[1])
    return out


def thm_join(text, experiment, members):
    lines, current, joined, inside = text.split("\n"), None, set(), False
    for i, line in enumerate(lines):
        t = line.split()
        if line.startswith("<data>"):
            inside = True
        elif line.startswith("</data>"):
            inside = False
        elif inside and t[:1] == ["dataset"]:
            current = t[1]
        elif inside and current in members and t[:1] == ["thm"]:
            lines[i] = f"    {'thm':<14} {experiment}"
            joined.add(current)
    missing = set(members) - joined
    if missing:
        raise ValueError(f"no THM data set {sorted(missing)}")
    return "\n".join(lines)


def thm_section(text, block):
    new = f"<thm>\n{block}\n</thm>" if block else None
    m = re.search(r"^<thm>\n.*?^</thm>\n?", text, re.S | re.M)
    if m:
        return text[:m.start()] + (new + "\n" if new else "") + text[m.end():]
    return text.rstrip("\n") + "\n" + new + "\n" if new else text


def thm_set(text, block, members=()):
    text = thm_section(text, block)
    names = sorted(set(re.findall(r"^\s*experiment\[([^\]]*)\]", block or "", re.M)))
    if not members or not names:
        return text
    if len(names) > 1:
        raise ValueError("the block defines several experiments; join them one by one")
    return thm_join(text, names[0], members)


def no_fold(text):
    out, skip = [], False
    for line in text.split("\n"):
        if line.startswith("<targets>"):
            out.append(line)
            skip = True
            continue
        if line.startswith("</targets>"):
            skip = False
        if skip or line.split()[:1] in (["beam"], ["quadrature"]):
            continue
        out.append(line)
    return "\n".join(out)


def thm_datasets(project):
    return [n for n, r in project.records().items() if r.record_kind == "dataset" and "thm" in r]


def _key_text(key, value):
    if isinstance(value, bool):
        return "on" if value else "off"
    if isinstance(value, (tuple, list)):
        if key in ("theta", "spectatorAngles"):
            return "-".join(repr(v) if not isinstance(v, str) else v for v in value)
        joiner = ";" if key == "cbackground" else ","
        return joiner.join(v if isinstance(v, str) else repr(v) for v in value)
    return value


def set_experiment(project, name, members, **keys):
    names = thm_datasets(project)
    members = [names[m - 1] if isinstance(m, int) else m for m in members]
    unknown = sorted(set(members) - set(names))
    if unknown:
        raise ValueError(f"no THM data set {unknown}")
    th = project.thm() or {"experiments": {}}
    had = th["experiments"].get(name, {})
    changes = {k: None for k in had if k not in keys}
    changes.update({k: _key_text(k, v) for k, v in keys.items() if v is not None})
    project.set_thm_experiment(name, datasets=[n for n in members if project[n]["thm"] != name], **changes)
    for n in names:
        if n not in members and project[n]["thm"] == name:
            project.set(n, thm="global")
    return project


_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))


def session_refusal(project, path):
    project.save(path)
    p = subprocess.run([sys.executable, "-c", "import sys; sys.path.insert(0, sys.argv[1]); from pyazr import azure2;"
                        " azure2(sys.argv[2], cwd=sys.argv[3])", _ROOT, path, os.path.dirname(path)],
                       capture_output=True, text=True, timeout=600)
    return p.returncode != 0, p.stdout + p.stderr


def refusal(make, path):
    try:
        project = make()
    except ValueError as err:
        return str(err)
    refused, out = session_refusal(project, path)
    return out if refused else None
