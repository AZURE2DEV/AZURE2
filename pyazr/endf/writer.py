"""A session's model written as an ENDF-6 material (ENDF-102, 2023 edition, cited by printed page;
docs/upgrade/endf.md sections 4 and 11).

What is written:
- MF=1/MT=451 (pp. 38-42): the header (LRP=2, NFOR=6, the projectile's AWI, EMAX, NSUB), the three standard text
  lines and HSUB, the text (what wrote the file, the project, the sentence on MF=2 and MF=3, the photon widths, the
  losses), the directory made from the sections (NC without SEND, p. 41).
- MF=2/MT=151 (pp. 64-65, 71-74): the particle channels in the formal basis, LRF=7 KRM=4 IFG=1
  (:func:`pyazr.endf.from_session`, the engine's door).
- MF=3 per reaction (pp. 103-104): every particle exit (MT 600 + k, 650 + k, ..., 800 + k, 50 + k) and the total
  capture (MT=102) from the engine on one grid (:mod:`pyazr.endf.grid`: the engine's refinement of an adaptive
  extrapolation grid, started at the level energies and thresholds, checked at quarter points) to the tolerance of
  linear interpolation, lab energies of the projectile in eV,
  barns, INT=2.  Before the file is written, MF=2 is read back into a project and its
  reaction cross sections compared with MF=3 at every grid point (the PI's decision: LRP=2 with the statement that
  MF=2 reproduces MF=3); a difference above ``MF2_TOLERANCE`` refuses the write.

What is not written is listed in the MF=1 text and in the report (never dropped silently): MF=3/MT=2 and MF=6
(the elastic channel needs MF=6 LAW=5, the exits MF=6 LAW=2, the primary photons MF=6/MT=102), the photon channels
(MF=1 text), and what :class:`pyazr.endf.Losses` names.
"""

import dataclasses
import datetime
import os
import tempfile

import numpy as np

from ..errors import FormatRefused
from . import grid, names
from .material import M_NEUTRON_U
from .model import EV_PER_MEV, from_session, to_project
from .records import RecordWriter, NUMBER_FORMS

#: LRP of MF=1 (pp. 38-39): 2, "parameters given in File 2 but not to be used in calculating cross sections, which
#: are assumed to be represented completely in File 3" (the PI's decision; p. 42 reserves it for derived files, which
#: CHECKR reports).  One named value, so that another choice is one line.
LRP = 2
#: MF=2 must reproduce MF=3 to this (relative, at every grid point where MF=3 is above 1e-12 of its largest value).
MF2_TOLERANCE = 1e-5
ELECTRON_MASS_U = 0.000548579909  # CODATA 2018
#: Why a section a project's reactions need is not written: (MF, role) -> the reason.  The writer lists every section
#: that a reaction of the project needs and the file does not hold (``Written.not_written`` and the MF=1 text).
NOT_BUILT = {(3, "elastic"): "charged-particle elastic scattering: MF=3/MT=2 goes with MF=6 LAW=5 (pp. 144-147, "
                             "150), which this writer does not write",
             (6, "elastic"): "LAW=5 of charged-particle elastic scattering (pp. 144-147) is not written",
             (6, "exit"): "LAW=2 of a two-body exit (p. 143) is not written",
             (6, "capture"): "the primary photons to each final state (LAW=2, p. 143) and the recoils (p. 136) are "
                             "not written"}


@dataclasses.dataclass(frozen=True)
class Written:
    """What :func:`write` wrote: the file, its sections (MF, MT, lines), the grid of each MF=3 section (points, lab
    eV range), the largest relative difference of MF=2's reconstruction from MF=3 per MT, the losses and the
    sections not written.  ``print()`` lists it."""
    path: str
    sections: tuple
    grids: dict
    mf2_check: dict
    losses: tuple
    not_written: tuple

    def __str__(self):
        out = [f"{self.path}"]
        for mf, mt, nc in self.sections:
            g = self.grids.get(mt) if mf == 3 else None
            extra = f", {g[0]} points {g[1]:.6g}-{g[2]:.6g} eV" if g else ""
            check = f", MF=2 reproduces it to {self.mf2_check[mt]:.1e}" if mf == 3 and mt in self.mf2_check else ""
            out.append(f"  MF={mf}/MT={mt} {names.reaction_name(mt)}: {nc} lines{extra}{check}")
        out += [f"  lost: {t}" for t in self.losses]
        out += [f"  not written: {t}" for t in self.not_written]
        return "\n".join(out)


def _symbol(z):
    from .. import _azure2
    return _azure2.element_symbols()[int(z)]


def _zsymam(z, a):
    """ZSYMAM (p. 40): Z right-justified in columns 1-3, '-', the symbol in 5-6, '-', A in 8-10."""
    return f"{int(z):3d}-{_symbol(z):<2s}-{int(a):3d}"


def _engine_grid(session, x, reactions, energy_range_MeV, tolerance, max_points, workdir):
    """The engine's values on one grid for every reaction: :func:`pyazr.endf.grid.refine` (the engine's refinement of
    an adaptive extrapolation grid with the quarter-point check) from :func:`pyazr.endf.grid.start_grid` (even points,
    the level energies and the thresholds).  Each sweep is one session of a project made of the session's nuclear
    part, at the fit vector x, with one extrapolation per reaction at the sweep's lab energies.  Returns (lab energies
    in MeV, {MT: values in b}, the refinement's report)."""
    from ..project import Project
    from ..session import Session
    path = os.path.join(workdir, "grid.azr")
    p = Project.new(path, nuclear_part=session.file)
    lo, hi = energy_range_MeV
    grid_file = os.path.join(workdir, "grid.dat")  # absolute: the project is also read without its directory
    np.savetxt(grid_file, [lo, hi], fmt="%.17g")  # read when the record is checked
    for k, (label, _, pair_keys) in enumerate(reactions):
        p.add_extrapolation(f"endf{k}", reaction=label, quantity="sigma", unit="b",
                            energies=f"file {grid_file} MeV lab", **pair_keys)
    for d in ("output", "checks"):
        os.makedirs(os.path.join(workdir, d), exist_ok=True)
    p.save()
    pf = _parameter_file(session, x, workdir)
    opts = {k: v for k, v in session.options.items() if k in ("use_brune", "use_park", "use_rmc", "use_gsl_coul",
                                                              "use_coulomb_table", "ignore_externals",
                                                              "use_long_wavelength", "ec_convection_current")}

    def model(xs):
        np.savetxt(grid_file, np.asarray(xs, float), fmt="%.17g")
        s = Session(path, extrapolations=True, parameter_file=pf, **opts)
        try:
            data = s.evaluate(s.parameters.x0).data
            return np.column_stack([np.asarray(d.model, float) for d in data])
        finally:
            s.close()

    refined = grid.refine(grid.start_grid(session, x, (lo, hi)), model, tolerance, max_points)
    e = np.asarray(refined["x"], float)
    v = np.asarray(refined["values"], float)
    return e, {mt: v[:, k] for k, (_, mt, _) in enumerate(reactions)}, refined


def _pair_keys(project, entrance_key, exit_key):
    """The pair keys an extrapolation of entrance -> exit names: a key only where its label names two pairs of the
    same nuclei (the engine's rule, E_PAIR_KEY / E_AMBIGUOUS_PAIR)."""
    pairs = {p.key: (p.light, p.heavy, p.light_excitation, p.heavy_excitation) for p in project.pairs()}
    out = {}
    for field, key in (("entrance_pair", entrance_key), ("exit_pair", exit_key)):
        if key > 0 and sum(v == pairs[key] for v in pairs.values()) > 1:
            out[field] = key
    return out


def _parameter_file(session, x, workdir):
    """Every parameter of the session at x, full precision, in the session's basis (a param.sav as the engine reads
    it by name)."""
    P = session.parameters
    x = P.x0 if x is None else np.asarray(x, float)
    path = os.path.join(workdir, "x.par")
    with open(path, "w") as fh:
        fh.write(f"#parametrization {session.basis}\n")
        for p in P:
            v = p.start if p.fixed else x[p.free_index]
            fh.write(f"{p.name} {float(v)!r}\n")
    return path


def _reconstruct(sections, s_entrance, energies_MeV, workdir):
    """MF=2 and MF=3 as written (MF=3 gives the residuals' excitations, QM - QI), read back into a project
    (:func:`pyazr.endf.to_project`) and its reaction cross sections at
    the given lab energies: {MT: values in b}."""
    from . import read
    path = os.path.join(workdir, "mf2.endf")
    _write_material(path, 9999, sections, header_rows=None, tpid="MF=2 check")
    m = read(path)
    project, conversion = to_project(m.resonances(), os.path.join(workdir, "mf2.azr"), separation_energy=s_entrance,
                                     material=m)
    key_of = {p["mt"]: p["key"] for p in conversion.pairs}
    out = {}
    for mt, e in energies_MeV.items():
        np.savetxt(os.path.join(workdir, f"e{mt}.dat"), e, fmt="%.17g")
        exit_key = key_of[mt]
        project.add_extrapolation(f"mt{mt}", reaction=project.reaction(1, exit_key), quantity="sigma", unit="b",
                                  energies=f"file {os.path.join(workdir, f'e{mt}.dat')} MeV lab",
                                  **_pair_keys(project, 1, exit_key))
    project.save()
    s = project.session(extrapolations=True)
    try:
        ev = s.evaluate(s.parameters.x0)
        for mt, d in zip(energies_MeV, ev.data):
            out[mt] = np.asarray(d.model, float)
    finally:
        s.close()
    return out


def _write_material(path, mat, sections, header_rows, tpid):
    """The lines of a tape (p. 32): TPID, the sections in (MF, MT) order with sequence numbers, SEND (NS 99999), FEND,
    MEND, TEND (p. 34).  ``header_rows`` are MF=1/MT=451's rows (None: no MF=1)."""
    allsec = dict(sections)
    if header_rows is not None:
        allsec[(1, 451)] = header_rows
    zero = NUMBER_FORMS["classic"](0.0) * 2 + f"{0:11d}" * 4

    def line(content, m, mf, mt, ns):
        return content.ljust(66) + f"{m:4d}{mf:2d}{mt:3d}{ns:5d}"

    out = [line(tpid[:66], 1, 0, 0, 0)]
    last_mf = None
    for mf, mt in sorted(allsec):
        if last_mf is not None and mf != last_mf:
            out.append(line(zero, mat, 0, 0, 0))
        for k, row in enumerate(allsec[(mf, mt)], start=1):
            out.append(line(row, mat, mf, mt, k))
        out.append(line(zero, mat, mf, 0, 99999))
        last_mf = mf
    out.append(line(zero, mat, 0, 0, 0))
    out.append(line(zero, 0, 0, 0, 0))
    out.append(line(zero, -1, 0, 0, 0))
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")


def _header_rows(record, header, sections, text, number_form):
    """MF=1/MT=451 (p. 42): HEAD, three CONT, NWD text lines, NXC directory lines (NC of every section and of this
    one; MOD = NMOD)."""
    ent = record.pairs[0]
    nxc = len(sections) + 1
    nwd = len(text)
    n451 = 4 + nwd + nxc
    w = RecordWriter(number_form)
    za_projectile = 1000 * int(round(ent.za)) + int(round(ent.ma * M_NEUTRON_U))
    # AWI: the projectile's mass in neutron masses, nuclear (p. 24): its AME2020 atomic mass minus Z electron masses
    # (electron binding neglected, < 1e-7 u for Z <= 2), not the project's mass, which a project may give nuclear
    from ..sources import masses
    z_projectile, a_projectile = int(round(ent.za)), int(round(ent.ma * M_NEUTRON_U))
    awi = (masses.atomic_mass((z_projectile, a_projectile)) - z_projectile * ELECTRON_MASS_U) / M_NEUTRON_U
    w.cont(record.za, record.awr, LRP, 0, header["NLIB"], header["NMOD"])
    w.cont(0.0, 0.0, 0, 0, 0, 6)
    w.cont(awi, record.eh, header["LREL"], 0, 10 * za_projectile, header["NVER"])
    w.cont(0.0, 0.0, 0, 0, nwd, nxc)
    for t in text:
        w.text(t)
    for mf, mt in [(1, 451)] + sorted(sections):
        nc = n451 if (mf, mt) == (1, 451) else len(sections[(mf, mt)])
        w.rows.append(" " * 22 + "".join(f"{v:11d}" for v in (mf, mt, nc, header["NMOD"])))
    return w.rows


def _wrap(text, width=66):
    """Text wrapped into TEXT records of at most 66 characters (p. 33)."""
    out, line = [], ""
    for word in text.split():
        if line and len(line) + 1 + len(word) > width:
            out.append(line)
            line = word
        else:
            line = f"{line} {word}" if line else word
    if line:
        out.append(line)
    return out


def write(session, path, x=None, *, header=None, energy_range_MeV=None, tolerance=1e-3, max_points=4000,
          number_form="widest", overwrite=False):
    """Write ``session``'s model at the fit vector ``x`` (default its start) as an ENDF-6 material at ``path``.

    ``header``: the evaluator's fields of MF=1 (p. 42), a dict with any of ``MAT`` (default 9999), ``NLIB``
    (0), ``NVER`` (1), ``LREL`` (0), ``NMOD`` (1), ``ALAB``, ``EDATE``, ``AUTH``, ``REF``, ``DDATE``, ``RDATE``,
    ``LIBRARY`` (HSUB's first line), ``TEXT`` (lines of free text).  ``energy_range_MeV``: (low, high) lab energy
    of the projectile (default: the data's, as the session has it); EL, EH of MF=2 and EMAX.  ``tolerance`` and
    ``max_points``: of the engine's adaptive grid (relative, for linear interpolation between its points).
    ``number_form``: ``"widest"`` (default), ``"classic"`` or ``"e"`` (:data:`pyazr.endf.NUMBER_FORMS`).

    Returns a :class:`Written`.  Refused (:class:`pyazr.FormatRefused`): an existing ``path`` without ``overwrite``,
    what :func:`pyazr.endf.from_session` refuses, a session without exit reactions, MF=2 that does not reproduce
    MF=3 to ``MF2_TOLERANCE``."""
    path = os.path.abspath(path)
    if os.path.exists(path) and not overwrite:
        raise FormatRefused(what=f"{path} exists", allowed="a new file (an existing file is replaced only when that is "
                            "asked for)", fix="write to a file name that does not exist yet, or ask for the "
                            "existing file to be replaced (overwrite)")
    if number_form not in NUMBER_FORMS:
        raise FormatRefused(what=f"number_form {number_form!r}", allowed=", ".join(NUMBER_FORMS),
                            fix="choose the widest form (the default): fixed point where it fits, 8-10 digits")
    h = dict(MAT=9999, NLIB=0, NVER=1, LREL=0, NMOD=1, ALAB="", EDATE="", AUTH="", REF="", DDATE="", RDATE="",
             LIBRARY="AZURE2", TEXT=())
    unknown = set(header or {}) - set(h)
    if unknown:
        raise FormatRefused(what=f"header fields {sorted(unknown)}", allowed=", ".join(h),
                            fix="name the MF=1 fields as listed")
    h.update(header or {})
    if energy_range_MeV is None:
        lab = [e for d in session.evaluate(session.parameters.x0 if x is None else x).data
               for e in np.asarray(d.extra.get("energy_lab", []), float)]
        if not lab:
            raise FormatRefused(what="the session has no data to take an energy range from",
                                allowed="an energy range given with the request (lab MeV of the projectile)",
                                fix="give the lowest and highest lab energy of the file")
        energy_range_MeV = (min(lab), max(lab))
    lo, hi = (float(v) for v in energy_range_MeV)
    record, losses = from_session(session, x, energy_range_eV=(lo * EV_PER_MEV, hi * EV_PER_MEV))
    from ..project import Project
    project = Project.open(session.file)
    pairs = {p.key: p.to_dict() for p in session.pairs}
    entrance = next(p for p in pairs.values() if p["is_entrance"])
    s_ent = entrance["sep_energy"] + entrance["excitation"]
    reactions = []  # (label, MT, the pair keys an extrapolation names when the label is ambiguous)
    keys = [k for k, p in pairs.items() if p["pair_kind"] == "particle" and k != entrance["key"]]
    mt_of = {}
    for k, p in zip([entrance["key"]] + keys, record.pairs):
        mt_of[k] = p.mt
    for k in keys:
        reactions.append((project.reaction(entrance["key"], k), mt_of[k],
                          _pair_keys(project, entrance["key"], k)))
    if any(p["pair_kind"] == "photon" for p in pairs.values()):
        reactions.append((project.reaction(entrance["key"], -1), 102, {}))
    if not reactions:
        raise FormatRefused(what="the session has no exit pair (elastic only)", allowed="a project with a reaction "
                            "or capture exit", fix="the elastic channel needs MF=6 LAW=5, which this writer does not "
                            "write yet")
    work = tempfile.mkdtemp(prefix="endf_write_")
    e_grid, grid, refined = _engine_grid(session, x, reactions, (lo, hi), tolerance, max_points, work)
    if refined.get("budget_reached"):
        raise FormatRefused(what=f"the grid reached max_points={max_points} before the tolerance {tolerance:g} "
                            f"(largest deviation {refined.get('largest_deviation'):.2e} at "
                            f"{refined.get('largest_deviation_at')} MeV)", allowed="a grid that meets the tolerance",
                            fix="allow more grid points, a looser tolerance or a narrower energy range")
    mf2 = record.rows(number_form)
    sections = {(2, 151): mf2}
    grids, values = {}, {}
    for _, mt, _ in reactions:
        order = np.argsort(e_grid)
        e_MeV, sigma = e_grid[order], grid[mt][order]
        if mt == 102:
            qm = qi = s_ent * EV_PER_MEV
        else:
            pk = next(k for k, m in mt_of.items() if m == mt)
            qi = record.pairs[[p.mt for p in record.pairs].index(mt)].q
            qm = qi + pairs[pk]["excitation"] * EV_PER_MEV
        w = RecordWriter(number_form)
        w.cont(record.za, record.awr, 0, 0, 0, 0)
        w.tab1(qm, qi, 0, 0, [len(e_MeV)], [2], list(e_MeV * EV_PER_MEV), list(sigma))
        sections[(3, mt)] = w.rows
        grids[mt] = (len(e_MeV), e_MeV[0] * EV_PER_MEV, e_MeV[-1] * EV_PER_MEV)
        values[mt] = (e_MeV, sigma)
    # MF=2 must reproduce MF=3 for the particle reactions (the photon channels are not in MF=2)
    particle = {mt: values[mt][0] for _, mt, _ in reactions if mt != 102}
    mf2_check = {}
    if particle:
        back = _reconstruct(sections, s_ent, particle, work)
        for mt, e in particle.items():
            sigma = values[mt][1]
            keep = np.abs(sigma) > 1e-12 * np.max(np.abs(sigma)) if np.any(sigma) else np.zeros(len(sigma), bool)
            mf2_check[mt] = float(np.max(np.abs(back[mt][keep] - sigma[keep]) / np.abs(sigma[keep]))) \
                if keep.any() else 0.0
            if mf2_check[mt] > MF2_TOLERANCE:
                raise FormatRefused(what=f"MF=2 as written reproduces MF=3/MT={mt} only to {mf2_check[mt]:.2e}",
                                    allowed=f"agreement to {MF2_TOLERANCE:g} (LRP=2: MF=2 reproduces MF=3)",
                                    fix=f"report the project to the AZURE2 developers (the files are kept in {work})")
    # every section a reaction of the project needs, and why it is not written where it is not
    needed = [(3, 2, "elastic"), (6, 2, "elastic")] + [(mf, mt, "capture" if mt == 102 else "exit")
                                                       for _, mt, _ in reactions for mf in (3, 6)]
    not_written = tuple(f"MF={mf}/MT={mt} ({names.reaction_name(mt)}): {NOT_BUILT[(mf, role)]}"
                        for mf, mt, role in needed if (mf, mt) not in sections)
    # MF=1
    ent = record.pairs[0]
    target_z, target_a = int(round(ent.zb)), int(round(ent.mb * M_NEUTRON_U))
    projectile = names.SUBLIBRARIES.get(10 * (1000 * int(round(ent.za)) + int(round(ent.ma * M_NEUTRON_U))),
                                        "INCIDENT PARTICLE")
    today = datetime.date.today()
    std = [f"{_zsymam(target_z, target_a):<11s}{h['ALAB'][:11]:<11s}{h['EDATE'][:10]:<10s} {h['AUTH'][:33]}",
           f" {h['REF'][:21]:<21s}{h['DDATE'][:10]:<10s} {h['RDATE'][:10]:<10s}"
           f"{'':12s}{today.strftime('%Y%m%d')}",
           f"---- {h['LIBRARY'][:16]:<16s}  MATERIAL {int(h['MAT']):4d}",
           f"----- {projectile.upper()} DATA",
           "------ ENDF-6"]
    body = [f"Written by pyazr.endf.write (AZURE2) on {today.isoformat()} from the project "
            f"{os.path.basename(session.file)} ({'Park' if session.basis == 2 else 'Brune'} fit vector); the "
            "particle channels converted to the formal basis with the engine's boundary values B_c = S_c(E_1) of "
            "each spin group.",
            f"LRP={LRP}: MF=3 holds the cross sections computed by AZURE2 on an adaptive grid (relative tolerance "
            f"{tolerance:g} of linear interpolation, not below 1e-3 of the largest value); MF=2 (LRF=7, KRM=4, "
            "IFG=1) reproduces the particle reactions "
            f"of MF=3 (checked before writing to {MF2_TOLERANCE:g}). MF=2 holds no photon channel: AZURE2 keeps "
            "radiative widths out of the total width; capture is MF=3/MT=102 (internal and external capture).",
            "Separation energy of the incident pair: "
            f"{s_ent:.9f} MeV (the project's); masses as the project's atomic masses."]
    text = std
    for paragraph in body:
        text = text + [""] + _wrap(paragraph)
    if losses.photon_widths:
        text += ["", "Photon widths of the levels (eV, at the level energy):"]
        for jpi, e, k, r, l, wv in losses.photon_widths:
            text.append(f"  {jpi:5s} E_x {e:11.6f} MeV to pair {k} ({r}{l}): {wv:.6g}")
    text += ["", "Not carried by this file:"]
    for t in losses.items + not_written:
        text += _wrap("- " + t)
    for t in h["TEXT"]:
        text += _wrap(t) or [""]
    header_rows = _header_rows(record, h, sections, text, number_form)
    _write_material(path, int(h["MAT"]), sections, header_rows,
                    tpid=f" pyazr.endf.write, {os.path.basename(session.file)}")
    listed_sections = [(1, 451, len(header_rows))] + [(mf, mt, len(rows))
                                                       for (mf, mt), rows in sorted(sections.items())]
    return Written(path, tuple(listed_sections), grids, mf2_check, tuple(losses.items), not_written)
