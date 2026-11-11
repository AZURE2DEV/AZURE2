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

A name that is not a THM data set of the text raises ValueError: a misspelt name must not leave a data set out of
its experiment unnoticed.
"""
import re


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
