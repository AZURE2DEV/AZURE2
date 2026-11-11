"""Compound nuclei of a project: which nucleus a pair forms, and merging several
single-nucleus projects into one file (natural-abundance targets).

A project can hold the level schemes of several compound nuclei.  No column
says which nucleus a level belongs to: the pair columns of its channel lines
do, by the rule the engine applies (``PPair``):

* particle pair (pType 0):  Z = Z1 + Z2,  A = round(M1 + M2)
* gamma pair (pType 10):    Z = Z2,       A = round(M2)
* beta pair (pType 20):     Z = Z2 - Z1,  A = round(M2)   (Z1 = -1 for beta-minus)

Pairs with the same (Z, A) form one nucleus; the engine keeps their J^pi
groups apart, so the level matrix is block-diagonal.  
The identity in the file is the COMPOUND nucleus: n + 20Ne forms 21Ne, n + 22Ne
forms 23Ne, n + 21Ne forms 22Ne.  The abundance belongs to the target isotope,
so a label shows both -- :meth:`Nucleus.label` gives the compound, and
:attr:`pyazr.ProjectLevel.nucleus` the compound nucleus of a level.  In a format-3 project a natural or enriched
target is a ``target`` record with ``isotope`` fields (the target nuclei); a data set on it sums one term per
isotope (``12C(p,p)12C + 13C(p,p)13C``), each weighted by its abundance.
"""
import os
from typing import NamedTuple

__all__ = ["SYMBOLS", "Nucleus", "nucleus_of", "nucleus_of_pair", "merge_projects", "MergeReport"]


def _element_symbols():
    """The engine's table of element symbols (physics/Elements.h), through the binding; also when this file is
    loaded by path outside its package."""
    import importlib
    import sys
    if __package__:
        return list(importlib.import_module(f"{__package__}._azure2").element_symbols())
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if root not in sys.path:
        sys.path.insert(0, root)
    return list(importlib.import_module("pyazr._azure2").element_symbols())


#: Element symbols by charge number (0 = the neutron, to 118), the engine's one table.
SYMBOLS = _element_symbols()


class Nucleus(NamedTuple):
    """A compound nucleus, identified by its charge and mass numbers."""
    Z: int
    A: int

    @property
    def symbol(self) -> str:
        return SYMBOLS[self.Z] if 0 <= self.Z < len(SYMBOLS) else f"Z{self.Z}"

    @property
    def label(self) -> str:
        """``"23Ne"`` style label of the compound nucleus."""
        return f"{self.A}{self.symbol}"

    @classmethod
    def parse(cls, spec):
        """``Nucleus``, ``(Z, A)`` or a label such as ``"23Ne"`` / ``"Ne23"``."""
        if isinstance(spec, Nucleus):
            return spec
        if isinstance(spec, (tuple, list)) and len(spec) == 2:
            return cls(int(spec[0]), int(spec[1]))
        s = str(spec).strip()
        digits = "".join(ch for ch in s if ch.isdigit())
        letters = "".join(ch for ch in s if ch.isalpha())
        if not digits or not letters:
            raise ValueError(f"cannot read a nucleus from {spec!r}; give (Z, A) or e.g. '23Ne'")
        sym = letters[0].upper() + letters[1:].lower()
        if sym not in SYMBOLS:
            raise ValueError(f"unknown element symbol {letters!r} in {spec!r}")
        return cls(SYMBOLS.index(sym), int(digits))

    def __str__(self):
        return self.label


def nucleus_of(ptype, Z1, Z2, M1, M2) -> Nucleus:
    """The engine's rule, from a pair's type, charges and masses."""
    if ptype == 20:
        return Nucleus(int(Z2) - int(Z1), int(round(M2)))
    if ptype == 10:
        return Nucleus(int(Z2), int(round(M2)))
    return Nucleus(int(Z1) + int(Z2), int(round(M1 + M2)))


def nucleus_of_pair(pair) -> Nucleus:
    """The compound nucleus of a runtime :class:`pyazr.parameters.Pair` (the
    engine reports it; this recomputes it from the pair's columns)."""
    return nucleus_of(pair.ptype, pair.Z1, pair.Z2, pair.M1, pair.M2)


class MergeReport(NamedTuple):
    """What :func:`merge_projects` did: per member (member index -> {old: new}) the pair-key map of its ``<levels>``
    rows and the record-name map (only the records renamed because an earlier member had the name), and warnings."""
    pair_map: list
    record_map: list
    warnings: list


def merge_projects(projects, out, *, isotopes=None, target="natural", isotope_sum="none", shared=()):
    """Merge several format-3 projects of single compound nuclei into one, written to ``out``.

    ``projects`` are paths or :class:`pyazr.Project` objects (none is modified).  The ``<levels>`` rows of all are
    put in one section by the engine (``_azure2.merge_level_rows``, experiment::MergeLevelRows): each member's pair
    keys move past those of the members before it, its level IDs likewise.  Every record of every member is kept,
    a name an earlier member uses renamed ``<name>_<k>`` (k the member's index from 1) with the references to it
    (:meth:`pyazr.Project.rename`); data-file paths are rewritten relative to ``out``.  Member 0's ``<config>``,
    ``<thm>``, ``<parameterSettings>`` and ``<mcmc>`` are kept; the other members' are not (listed in the warnings
    when they have any).

    ``isotopes``: ``{target nucleus: abundance}`` (``"12C": 0.6``, or the text of the field after the nucleus,
    ``"12C": "0.6 fit prior 0.6 0.02"``; a bare number is fixed) adds a ``target`` record named ``target`` (kind none) with one ``isotope``
    field per entry and ``isotope_sum`` as given.  ``shared`` names data sets or extrapolations the members have in
    common, one dict per set mapping member index -> that member's record name: the first member's record carries
    the set, its reaction the sum of every member's reaction and its target ``target``; the other members' records
    are dropped.  The project is validated as written: an error raises :class:`pyazr.FormatError`.
    """
    from .project import Project, FormatError, _bytes, _text, _engine_module as engine
    out = os.path.abspath(out)
    out_dir = os.path.dirname(out)
    members = []
    for m in projects:
        src = m if isinstance(m, Project) else Project.open(os.fspath(m))
        members.append((src, Project(src.to_bytes(), out)))  # a private copy, its data paths read from out's directory
    warnings = []
    rows, pair_map = engine().merge_level_rows([c._doc for _, c in members])

    # record names: a name an earlier member has gets the member's number
    used, record_map = set(), []
    for k, (_, c) in enumerate(members):
        renamed = {}
        for name in list(c.records()):
            if name in used:
                new = f"{name}_{k + 1}"
                c.rename(name, new)
                renamed[name] = new
            used.add(renamed.get(name, name))
        record_map.append(renamed)
    # data files relative to the merged project
    for (src, c) in members:
        for name, r in c.records().items():
            f = r.get("file")
            if isinstance(f, str) and f and not os.path.isabs(f.strip('"')):
                c.set(name, file=os.path.relpath(os.path.join(src.directory, f.strip('"')), out_dir))
    # the shared sets: the carrier's reaction is the sum of the members', the others dropped
    drop = set()
    carriers = {}
    for group in shared:
        group = {int(k): v for k, v in group.items()}
        first = min(group)
        names = {k: record_map[k].get(v, v) for k, v in group.items()}
        reactions = [members[k][1][names[k]]["reaction"] for k in sorted(group)]
        carriers[(first, names[first])] = " + ".join(reactions)
        drop.update((k, n) for k, n in names.items() if k != first)

    sections = {}
    for k, (_, c) in enumerate(members):
        for s in engine().record_sections(c._doc):
            body = sections.setdefault(_text(s["section"]), [])
            for r in s["records"]:
                name = _text(r["name"])
                if (k, name) in drop:
                    continue
                if (k, name) in carriers:
                    c.set(name, reaction=carriers[(k, name)], target=target)
                body.append(c[name].text)
        if k > 0:
            for tag in ("thm", "parameterSettings", "mcmc"):
                if any(_text(t[0]) == tag for t in c._doc.sections()):
                    warnings.append(f"project {k + 1}: <{tag}> not merged (project 1's is kept)")
    if isotopes:
        fields = [("kind", "none")] + [("isotope", f"{n} {v!r} fixed" if isinstance(v, (int, float)) else f"{n} {v}")
                                       for n, v in isotopes.items()] + \
                 [("isotope_sum", str(isotope_sum))]
        lines, problems = engine().record_lines("target", target, fields)
        if problems:
            raise FormatError(f"target {target!r}:", [{"severity": "error", "code": c, "record": target, "message": m}
                                                      for c, m in problems])
        sections.setdefault("targets", []).insert(0, "\n".join(_text(x) for x in lines))

    merged = Project(members[0][0].to_bytes(), out)
    doc = merged._doc
    for s in sorted(engine().record_sections(doc), key=lambda s: -s["open"]):
        name = _text(s["section"])
        if s["close"] > 0:
            body = "\n\n".join(sections.pop(name, [])).split("\n") if sections.get(name) else []
            doc.replace_lines(s["open"], s["close"] - 1, [_bytes(x) for x in body])
    levels = [t for t in doc.sections() if _text(t[0]) == "levels"][0]
    doc.replace_lines(levels[1] + 1, levels[2], [_bytes(r) for r in rows] + [b""])
    for name, records in sections.items():  # sections member 0 does not have
        if records:
            merged._replace_record_section(name, "\n\n".join(records).split("\n"))
    errors = [d for d in merged._diagnostics()[0] if d.get("severity") == "error"]
    if errors:
        raise FormatError(f"{out}: the merged project has {len(errors)} error(s)", errors)
    merged.save(out)
    return MergeReport([dict(m) for m in pair_map], record_map, warnings)
