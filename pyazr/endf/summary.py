"""What an ENDF-6 material holds, in words an evaluator recognises (``Material.summary()``).

From MF=1/MT=451 (pp. 38-42): MAT, the target (ZSYMAM), the projectile (NSUB, Table 3 p. 7), the library (HSUB),
lab, date and authors, LRP (pp. 38-39), EMAX; per section its file and reaction (Table 4 p. 13, Appendix B) and its
lines; for MF=3 the Q values and the energy range (pp. 103-104), for MF=6 the frame and each product's law (pp.
134-147), for MF=2 the isotopes and ranges with their format and, for LRF=7, the R-matrix flags and counts (pp.
63-64, 71-74).  The descriptive text follows.  Nothing here interprets physics.
"""

import dataclasses
import os

from ..errors import Error
from . import names


@dataclasses.dataclass
class Summary:
    """The summary as data (``to_dict()``, for ``pyazr endf summary --json``) that prints itself."""
    path: str
    mat: int
    target: str
    projectile: str
    sublibrary: str
    library: str
    lab: str
    date: str
    authors: str
    lrp: str
    emax_mev: float
    awr: float
    awi: float
    sections: list
    text: list

    def to_dict(self):
        return dataclasses.asdict(self)

    def __str__(self):
        out = [f"{os.path.basename(self.path)}: MAT {self.mat}, {self.target} + {self.projectile} "
               f"({self.sublibrary}); {self.library}",
               f"  lab {self.lab}, {self.date}, {self.authors}",
               f"  AWR {self.awr:g}, AWI {self.awi:g} (neutron masses); EMAX {self.emax_mev:g} MeV; LRP: {self.lrp}",
               "  sections:"]
        for s in self.sections:
            out.append(f"    MF={s['mf']:<3} MT={s['mt']:<4} {s['reaction']:<28} {s['lines']:>6} lines  {s['what']}")
        out.append("  text:")
        out += [f"    {line.rstrip()}" for line in self.text]
        return "\n".join(out)


def _mev(e_ev):
    return f"{e_ev * 1e-6:.6g}"


def _describe_section(material, mf, mt):
    try:
        if mf == 3:
            xs = material.cross_section(mt)
            lo, hi = xs.energy_range
            return (f"sigma(E) at {len(xs.table.x)} energies, {_mev(lo)}-{_mev(hi)} MeV lab; QM {_mev(xs.qm)} MeV, "
                    f"QI {_mev(xs.qi)} MeV")
        if mf == 6:
            d = material.distribution(mt)
            frame = {1: "lab", 2: "c.m.", 3: "c.m. light, lab heavy", 4: "c.m., breakup lab"}.get(d.lct, d.lct)
            parts = []
            for p in d.products:
                extra = ""
                if p.law in (2, 5) and p.panels:
                    flags = sorted({q.flag for q in p.panels})
                    extra = (f" {'LANG' if p.law == 2 else 'LTP'}={','.join(map(str, flags))}, "
                             f"{len(p.panels)} energies {_mev(p.panels[0].energy)}-{_mev(p.panels[-1].energy)} MeV")
                    if p.law == 5:
                        extra += f", LIDP={p.lidp}"
                parts.append(f"ZAP={p.zap} LAW={p.law}{extra}")
            return f"LCT={d.lct} ({frame}); " + "; ".join(parts)
        if mf == 2 and mt == 151:
            return _describe_mf2(material)
        if mf == 1 and mt == 451:
            return "header, text and directory"
    except Error as e:
        return f"not read: {e}"
    return ""


def _describe_mf2(material):
    c = material.section(2, 151).cursor(material.path)
    head = c.cont()
    parts = []
    for _ in range(max(head.n1, 0)):
        iso = c.cont()
        for _ in range(max(iso.n1, 0)):
            rng = c.cont()
            fmt = names.RESOLVED_FORMATS.get(rng.l2, f"LRF={rng.l2}") if rng.l1 == 1 else f"LRU={rng.l1} LRF={rng.l2}"
            text = f"ZAI {iso.c1:g}: {_mev(rng.c1)}-{_mev(rng.c2)} MeV, {fmt}"
            if rng.l1 == 1 and rng.l2 == 7 and rng.n1 == 0:
                rml = c.cont()
                pairs = c.list()
                text += (f" (KRM={rml.l2} {names.RML_FORMULAS.get(rml.l2, '')}, IFG={rml.l1}, KRL={rml.n2}, "
                         f"{pairs.l1} pairs, {rml.n1} spin groups)")
            parts.append(text)
            return "; ".join(parts)  # the records after the first range are read by Material.resonances()
    return "; ".join(parts)


def summarize(material):
    h = material.require_header()
    std = h.standard_lines
    library = std["HSUB"] or f"NLIB={h.nlib} NVER={h.nver}"
    sections = []
    for (mf, mt), s in sorted(material.sections.items()):
        sections.append({"mf": mf, "mt": mt, "file": names.FILES.get(mf, f"MF={mf}"),
                         "reaction": names.reaction_name(mt), "lines": len(s.lines),
                         "what": _describe_section(material, mf, mt)})
    return Summary(path=material.path, mat=material.mat, target=std["ZSYMAM"],
                   projectile=names.PARTICLES.get(h.projectile_za, f"ZA={h.projectile_za}"),
                   sublibrary=names.SUBLIBRARIES.get(h.nsub, f"NSUB={h.nsub}"), library=library, lab=std["ALAB"],
                   date=std["EDATE"], authors=std["AUTH"], lrp=f"{h.lrp} ({names.LRP.get(h.lrp, 'not defined')})",
                   emax_mev=h.emax * 1e-6, awr=h.awr, awi=h.awi, sections=sections, text=[t for t in h.text])
