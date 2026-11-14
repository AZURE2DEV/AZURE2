"""An evaluation's resonance parameters as an AZURE2 model, and an AZURE2 model as resonance parameters (ENDF-102,
2023 edition, cited by printed page; docs/upgrade/endf.md sections 3 and 7.1).

- :func:`to_project`: an :class:`~pyazr.endf.resonances.LRF7` record into a new format-3 project in the Brune basis:
  pairs, levels, channels, every level and width fixed.  ER and GAM go from the lab scale of the incident pair to
  E_x and the c.m. scale (pp. 66, 355-356), the file's boundary values B_c through the engine's door
  ``observed_levels_of_block`` (formal levels with a caller's B_c -> Brune energies, widths and ANCs).
- :func:`from_session`: the particle channels of a session's nucleus at a fit vector as an LRF=7 record in the
  formal basis, through the engine's door ``formal_levels_of_block`` (the engine's B_c = S_c(E_1)), and the list of
  what the record cannot carry.

The frame and unit factors (E_lab / E_cm of the incident pair, eV <-> MeV, fm <-> 1e-12 cm, u <-> neutron masses)
are the format's and live here and in :mod:`pyazr.endf.resonances`; every basis change is the engine's.
"""

import dataclasses
import math
import os
import warnings

import numpy as np

from ..errors import FormatRefused, listed
from . import names
from .material import M_NEUTRON_U
from .resonances import LRF7, Channel, Level, ParticlePair, SpinGroup, amplitude_scale, lab_factor

EV_PER_MEV = 1e6
#: The Brune iteration of to_project stops when a level moves by at most this (MeV): tighter than the engine's
#: 1e-6 MeV of a fit's report, so that the levels written reproduce the file's R-matrix to the digits it holds.
BRUNE_TOLERANCE_MEV = 1e-11
FM_PER_ENDF_LENGTH = 10.0  # radii in 1e-12 cm (p. 73) = 10 fm


def _engine():
    from .. import _azure2
    return _azure2


#: The smallest format-3 project the engine reads: the <config> lines of every format-3 project (an A-matrix flag,
#: the output and checks directories, eight unused paths, the format line) and the empty sections a project has.
#: to_project starts from it and adds every pair and level with the engine's edits (Project.add_pair, add_level).
EMPTY_PROJECT = ("<config>\ntrue\noutput/\nchecks/\n" + "none\n" * 8 + "format=3\n</config>\n<levels>\n</levels>\n"
                 "<data>\n</data>\n<extrapolations>\n</extrapolations>\n<targets>\n</targets>\n")


# --------------------------------------------------------------------------------------------- to a project
@dataclasses.dataclass(frozen=True)
class ConvertedLevel:
    """One level of the new project: its spin group, the file's ER (eV, lab) and B_c, the formal E_x (MeV) the
    conversion made of ER, the Brune E_x (MeV) and the observed values written (width in eV or ANC in fm^-1/2 per
    channel), the engine's status of the Brune iteration, and whether ER lies above EH (a background level)."""
    jpi: str
    er_eV: float
    formal_MeV: float
    energy_MeV: float
    values: tuple
    kinds: tuple
    status: str
    background: bool


@dataclasses.dataclass(frozen=True)
class Conversion:
    """What :func:`to_project` did: the project's path, the separation energy of the incident pair and where it came
    from, each particle's mass in the file against AME2020 (keV), the pairs written (key, MT, Q from the file,
    separation energy, excitation), the levels, the spin groups not written and the notes.  ``print()`` lists it."""
    path: str
    separation_MeV: float
    separation_source: str
    mass_checks: tuple
    pairs: tuple
    levels: tuple
    skipped: tuple
    notes: tuple

    def __str__(self):
        out = [f"project {self.path}",
               f"  separation energy of the incident pair: {self.separation_MeV:.9f} MeV ({self.separation_source})"]
        for name, file_u, ame_u in self.mass_checks:
            diff = "not in AME2020" if ame_u is None else f"{(file_u - ame_u) * 931494.10242:+.3f} keV against AME2020"
            out.append(f"  mass of {name} in the file: {file_u:.9f} u ({diff})")
        for p in self.pairs:
            out.append(f"  pair {p['key']}: MT={p['mt']} {names.reaction_name(p['mt'])}, Q = {p['q_eV']:.10g} eV, "
                       f"S = {p['separation_MeV']:.9f} MeV, excitation {p['excitation_MeV']:g} MeV")
        for lv in self.levels:
            vals = ", ".join(f"{v:.8g} {'eV' if k == 'width' else 'fm^-1/2'}" for v, k in zip(lv.values, lv.kinds))
            out.append(f"  {lv.jpi}: ER {lv.er_eV:.10g} eV -> formal E_x {lv.formal_MeV:.7f} MeV -> E_x "
                       f"{lv.energy_MeV:.7f} MeV ({lv.status}"
                       f"{', above EH: background' if lv.background else ''}); {vals}")
        out += [f"  not written: {s}" for s in self.skipped]
        out += [f"  note: {n}" for n in self.notes]
        return "\n".join(out)


def _nuclide(z, mass_u):
    from ..sources import masses
    entry = masses.identify(int(round(z)), mass_u)
    return entry


def _pair_name(z, mass_u):
    e = _nuclide(z, mass_u)
    return f"{int(round(mass_u))}{_engine().element_symbols()[int(round(z))]}" if e is None else \
        f"{e.A}{_engine().element_symbols()[e.Z]}"


def _excitation_of(material, mt):
    """The residual's excitation of the pair of reaction MT, from MF=3/MT's QM - QI (p. 105), or None."""
    if material is None or (3, mt) not in material.sections:
        return None
    xs = material.cross_section(mt)
    return (xs.qm - xs.qi) / EV_PER_MEV


def to_project(lrf7, path, *, separation_energy="ame2020", material=None, overwrite=False):
    """A new format-3 project at ``path`` (never an open one: an existing file is refused unless ``overwrite``) made
    of the resonance parameters ``lrf7`` (:meth:`pyazr.endf.Material.resonances`), in the Brune basis, every level
    and width fixed.  Returns ``(Project, Conversion)``; the project is saved.

    ``separation_energy``: of the incident pair from the compound nucleus, ``"ame2020"`` (default; the file holds no
    mass of the compound nucleus) or a number in MeV.  The other pairs take the file's Q values (S_pair = S_incident
    - Q - excitation); the residual's excitation is MF=3's QM - QI of the pair's MT when ``material`` (the
    :class:`~pyazr.endf.Material` the record came from) has that section, else 0.  Masses are the file's (MA, MB
    times the neutron mass, p. 394); each is compared with AME2020 in the report.  Levels keep their order; spin
    groups without levels are not written (listed)."""
    from ..project import Project
    from ..sources import masses
    path = os.path.abspath(path)
    if os.path.exists(path) and not overwrite:
        raise FormatRefused(what=f"{path} exists", allowed="a new file (to_project never writes over a project)",
                            fix="give another path, or overwrite=True to replace that file")
    entrance = lrf7.entrance
    ent = lrf7.pairs[entrance - 1]
    f_lab = lab_factor(ent)
    m_a, m_b = ent.ma * M_NEUTRON_U, ent.mb * M_NEUTRON_U
    if separation_energy == "ame2020":
        try:
            s_ent = masses.separation_energy((int(round(ent.za)), int(round(m_a))), (int(round(ent.zb)),
                                                                                      int(round(m_b))))
        except KeyError as e:
            raise FormatRefused(what=f"{lrf7.source}: the incident pair is not in AME2020 ({e})",
                                allowed="a pair of nuclides AME2020 lists, or separation_energy= in MeV",
                                fix="to_project(..., separation_energy=<S in MeV>)") from None
        source = "AME2020"
    else:
        s_ent = float(separation_energy)
        source = "given"
    notes, skipped, mass_checks, pair_rows = [], [], [], []
    seen = set()
    for p in lrf7.pairs:
        for z, m in ((p.za, p.ma * M_NEUTRON_U), (p.zb, p.mb * M_NEUTRON_U)):
            name = _pair_name(z, m)
            if name in seen:
                continue
            seen.add(name)
            e = _nuclide(z, m)
            mass_checks.append((name, m, None if e is None else e.atomic_mass))
    # the pairs, the incident one first (key 1)
    order = [entrance] + [k for k in range(1, len(lrf7.pairs) + 1) if k != entrance]
    used = {c.pair for g in lrf7.groups if g.levels for c in g.channels}
    radius = {c.pair: c.radius * FM_PER_ENDF_LENGTH for g in lrf7.groups for c in g.channels}
    key_of, records = {}, []
    project = Project(EMPTY_PROJECT, path)
    for k in order:
        p = lrf7.pairs[k - 1]
        if k not in used:
            skipped.append(f"pair {k} (MT={p.mt}): no level has a channel of it")
            continue
        excitation = 0.0
        if k != entrance:
            e = _excitation_of(material, p.mt)
            if e is None and p.mt not in (2, 600, 650, 700, 750, 800):
                notes.append(f"pair {k} (MT={p.mt}): the residual's excitation is not in MF=2 and MF=3/MT={p.mt} "
                             "is not given: written as a ground-state pair at the file's threshold")
            excitation = e or 0.0
        sep = s_ent - p.q / EV_PER_MEV - excitation
        light = dict(J=p.spin_a, parity=p.parity_a, Z=int(round(p.za)), mass=p.ma * M_NEUTRON_U)
        heavy = dict(J=p.spin_b, parity=p.parity_b, Z=int(round(p.zb)), mass=p.mb * M_NEUTRON_U)
        key = project.add_pair(light=light, heavy=heavy, separation_energy=sep, radius=radius[k],
                               heavy_excitation=excitation)
        key_of[k] = key
        pair_rows.append(dict(key=key, mt=p.mt, q_eV=p.q, separation_MeV=sep, excitation_MeV=excitation))
        records.append(dict(key=key, pair_kind="particle", Z1=light["Z"], Z2=heavy["Z"], M1=light["mass"],
                            M2=heavy["mass"], J1=light["J"], J2=heavy["J"], parity1=light["parity"],
                            parity2=heavy["parity"], excitation=excitation, light_excitation=0.0,
                            sep_energy=sep, channel_radius=radius[k]))
    scale = amplitude_scale(ent)
    levels = []
    for g in lrf7.groups:
        if not g.levels:
            skipped.append(f"spin group {g.jpi}: no level (its hard-sphere phases are AZURE2's own)")
            continue
        chans = [(key_of[c.pair], c.l, int(round(2 * c.s))) for c in g.channels]
        E = [lv.energy / f_lab / EV_PER_MEV + s_ent for lv in g.levels]
        gamma = [v / scale for lv in g.levels for v in lv.gam]
        B = [c.boundary for c in g.channels]
        out = _engine().observed_levels_of_block(records, int(round(2 * g.J)), g.parity, chans, E, gamma, B,
                                                 energy_tolerance_MeV=BRUNE_TOLERANCE_MEV)
        nc = len(chans)
        for i, lv in enumerate(g.levels):
            vals = out["value"][i * nc:(i + 1) * nc]
            kinds = out["kind"][i * nc:(i + 1) * nc]
            channels = [dict(pair=c[0], l=c[1], s=c[2] / 2, width=v, fixed=True) for c, v in zip(chans, vals)]
            project.add_level(g.J, g.parity, out["E"][i], channels, fixed=True)
            levels.append(ConvertedLevel(g.jpi, lv.energy, E[i], out["E"][i], tuple(vals), tuple(kinds),
                                         out["status"][i], lv.energy > lrf7.eh))
            if out["status"][i] not in ("converged", "bracketed"):
                notes.append(f"{g.jpi} level at ER = {lv.energy:.10g} eV: the Brune iteration ended "
                             f"'{out['status'][i]}'; its formal values were written")
        notes.extend(out["warnings"])
    notes.append("the particles' g-factors are not in an ENDF file: 0 written")
    project.save()
    return project, Conversion(path, s_ent, source, tuple(mass_checks), tuple(pair_rows), tuple(levels),
                               tuple(skipped), tuple(notes))


# ---------------------------------------------------------------------------------------------- from a session
@dataclasses.dataclass(frozen=True)
class Losses:
    """What an LRF=7 record of a session cannot carry (docs/upgrade/endf.md section 3.1), one sentence each, and the
    photon widths of the session's levels (Gamma_gamma in eV at the level, as ``Session.transform`` gives them), which
    the writer lists in the MF=1 text."""
    items: tuple
    photon_widths: tuple

    def __str__(self):
        return "\n".join([f"  lost: {t}" for t in self.items] +
                         [f"  photon width: {jpi} E_x {e:.6f} MeV, to pair {k} (E{l if r == 'E' else ''}"
                          f"{'M' + str(l) if r == 'M' else ''}): {w:.6g} eV" for jpi, e, k, r, l, w in
                          self.photon_widths])


def _refuse_session(what, allowed, fix):
    return FormatRefused(what=what, allowed=allowed, fix=fix)


def _mt_of_pairs(particle_pairs, entrance_key):
    """MT of each particle pair (Appendix B, pp. 339-341): 2 for the incident pair; 600 + k (protons), 650 + k (d),
    700 + k (t), 750 + k (3He), 800 + k (alphas), 50 + k (neutrons), k the order of the pair's excitation among the
    pairs of the same emitted particle (0 the lowest)."""
    by_family = {}
    out = {}
    for p in particle_pairs:
        if p["key"] == entrance_key:
            out[p["key"]] = 2
            continue
        za = 1000 * int(p["Z1"]) + int(round(p["M1"]))
        first = names.family_of(za)
        if first is None:
            raise _refuse_session(f"pair {p['key']}: the emitted particle (Z={p['Z1']}, A={round(p['M1'])}) has no "
                                  "discrete-level reaction family in ENDF-6 (Appendix B, pp. 339-341)",
                                  "exit pairs emitting n, p, d, t, 3He or an alpha",
                                  "remove the pair from the session's project before writing")
        by_family.setdefault(first, []).append((p["excitation"], p["key"]))
    for first, rows in by_family.items():
        for k, (_, key) in enumerate(sorted(rows)):
            out[key] = first + k
    return out


def _full_values(session, x):
    """Every parameter's value at the fit vector x (fixed ones their start), in the order of session.parameters."""
    P = session.parameters
    x = P.x0 if x is None else np.asarray(x, float)
    return [float(x[p.free_index]) if not p.fixed else float(p.start) for p in P]


def from_session(session, x=None, *, energy_range_eV=(1e-5, 2e7)):
    """The particle channels of ``session``'s nucleus at the fit vector ``x`` (default its start) as an
    :class:`~pyazr.endf.resonances.LRF7` record in the formal basis (KRM=4, IFG=1, KRL=0; the engine's
    B_c = S_c(E_1) written as BND), and the :class:`Losses`.  ``energy_range_eV``: EL, EH (lab energy of the incident
    particle, eV).  ER and GAM are on the lab scale of the incident pair (pp. 66, 355-356); MA, MB in neutron masses
    (p. 72); Q from the separation energies and excitations (p. 105 for its sign).

    Refused (:class:`pyazr.FormatRefused`): a session in the formal basis (``use_brune=False``: its B_c is not the
    one of its present energies), relativistic kinematics (KRL=1 not decided), a nuclear potential (hybrid model),
    Reich-Moore (``use_rmc``: its eliminated widths act on the particle channels, which MF=2 without photon channels
    would not reproduce), several compound nuclei, an exit pair outside ENDF's reaction families, a block that has no
    formal parameters."""
    opts = session.options
    if session.basis == 0:
        raise _refuse_session("the session is in the formal basis (use_brune=False)", "a session in the Brune or "
                              "Park basis (its formal levels and B_c are computed at its own energies)",
                              "open the project with pyazr.Session(path) (Brune is the default)")
    if opts.get("use_rmc"):
        raise _refuse_session("the session uses Reich-Moore (use_rmc): its eliminated photon widths change the "
                              "particle channels", "a session without use_rmc (photon channels stay out of MF=2, "
                              "decision 1 of the design)", "open the project without use_rmc")
    from ..project import Project
    project = Project.open(session.file)
    config = project.config()
    if str(config.get("relativistic", "0")).strip() not in ("0", "false", ""):
        raise _refuse_session("the project uses relativistic kinematics (relativistic=1)", "relativistic=0 (KRL=0)",
                              "KRL=1 has no equations in ENDF-102 (p. 71): write with relativistic=0, or ask for "
                              "the KRL=1 convention to be decided")
    if any(p.enabled for p in session.nuclear_potentials().values()):
        raise _refuse_session("a pair of the session has a nuclear potential (the hybrid model)", "plain Coulomb "
                              "functions (LRF=7 has no place for a potential: KPS would carry phases, not P and S)",
                              "remove the <potential> record before writing")
    pairs = [p.to_dict() for p in session.pairs]
    nuclei = {(p["nucleus_Z"], p["nucleus_A"]) for p in pairs}
    if len(nuclei) > 1:
        raise _refuse_session(f"the session has {len(nuclei)} compound nuclei", "one compound nucleus (one isotope of "
                              "MF=2)", "write each compound nucleus from a project of its own")
    particle = [p for p in pairs if p["pair_kind"] == "particle"]
    entrance = next(p for p in particle if p["is_entrance"])
    mts = _mt_of_pairs(particle, entrance["key"])
    order = [entrance["key"]] + [p["key"] for p in particle if p["key"] != entrance["key"]]
    by_key = {p["key"]: p for p in pairs}
    s_ent = entrance["sep_energy"] + entrance["excitation"]
    endf_pairs = []
    for key in order:
        p = by_key[key]
        q = (s_ent - p["sep_energy"] - p["excitation"]) * EV_PER_MEV
        ia = p["J1"] * (p["parity1"] if p["J1"] != 0 else 1)
        ib = p["J2"] * (p["parity2"] if p["J2"] != 0 else 1)
        endf_pairs.append(ParticlePair(p["M1"] / M_NEUTRON_U, p["M2"] / M_NEUTRON_U, float(p["Z1"]), float(p["Z2"]),
                                       ia, ib, q, 1, 1, mts[key], float(p["parity1"]), float(p["parity2"])))
    ppi = {key: i + 1 for i, key in enumerate(order)}
    lost = ["external capture, the photon channels and their energy dependence: not in MF=2 (MF=3/MT=102 and "
            "MF=6/MT=102 hold the engine's capture)",
            "the data, normalizations, energy shifts, targets, beams and detectors of the project",
            "the particles' g-factors and the light particle's excitation are not LRF=7 fields"]
    for p in pairs:
        if p["pair_kind"] == "beta":
            lost.append(f"pair {p['key']}: a beta-decay pair has no place in ENDF-6")
        if p["pair_kind"] == "particle" and p["light_excitation"]:
            lost.append(f"pair {p['key']}: the light particle's excitation ({p['light_excitation']} MeV) is in Q only")
        if p["pair_kind"] == "particle" and p["Z1"] == p["Z2"] and abs(p["M1"] - p["M2"]) < 1e-9:
            lost.append(f"pair {p['key']}: identical particles; MF=2 has no field for it (the writer sets LIDP in "
                        "MF=6 LAW=5)")
    values = _full_values(session, x)
    physical = session.transform(session.parameters.x0 if x is None else np.asarray(x, float))
    rows = [p.to_dict() for p in session.parameters]
    groups = {}
    for i, r in enumerate(rows):
        if r["kind"] == "energy":
            g = groups.setdefault(r["jgroup"], {"J": r["J"], "parity": r["parity"], "jpi": r["jpi"], "levels": {},
                                                "channels": {}})
            g["levels"][r["level"]] = {"E": values[i], "gamma": {}, "label": r["level_label"]}
    photon_widths = []
    for i, r in enumerate(rows):
        if r["kind"] != "width":
            continue
        g = groups[r["jgroup"]]
        if r["radiation_type"] == "P":
            g["channels"][r["channel"]] = (r["pair_key"], r["L"], r["S"])
            g["levels"][r["level"]]["gamma"][r["channel"]] = values[i]
        elif r["radiation_type"] in ("E", "M") and physical[i] != 0.0:
            lev = g["levels"][r["level"]]
            photon_widths.append((g["jpi"], lev["E"], r["pair_key"], r["radiation_type"], r["L"], float(physical[i])))
    spin_groups = []
    basis = "park" if session.basis == 2 else "brune"
    ent_pair = endf_pairs[0]
    scale = amplitude_scale(ent_pair)
    f_lab = lab_factor(ent_pair)
    for gk in sorted(groups):
        g = groups[gk]
        chan_ids = sorted(g["channels"])
        if not chan_ids:
            lost.append(f"spin group {g['jpi']}: no particle channel (photon channels only)")
            continue
        lev_ids = sorted(g["levels"])
        E = [g["levels"][l]["E"] for l in lev_ids]
        gam = [g["levels"][l]["gamma"].get(c, 0.0) for l in lev_ids for c in chan_ids]
        chans = [(g["channels"][c][0], g["channels"][c][1], int(round(2 * g["channels"][c][2]))) for c in chan_ids]
        try:
            out = _engine().formal_levels_of_block(pairs, int(round(2 * g["J"])), g["parity"], chans, E, gam, [],
                                                   "amplitudes", basis, bool(opts.get("use_coulomb_table", True)),
                                                   bool(opts.get("use_gsl_coul", False)))
        except ValueError as e:
            raise _refuse_session(f"spin group {g['jpi']}: no formal parameters correspond ({e})",
                                  "levels whose Brune parameters have a formal equivalent (M positive definite)",
                                  "the engine refuses the same set in the formal basis; change the levels") from None
        nc = len(chan_ids)
        channels = tuple(Channel(ppi[c[0]], c[1], c[2] / 2, out["B"][k], by_key[c[0]]["channel_radius"] /
                                 FM_PER_ENDF_LENGTH, by_key[c[0]]["channel_radius"] / FM_PER_ENDF_LENGTH)
                         for k, c in enumerate(chans))
        levels = tuple(Level((out["E"][i] - s_ent) * f_lab * EV_PER_MEV,
                             tuple(v * scale for v in out["gamma"][i * nc:(i + 1) * nc])) for i in range(len(E)))
        spin_groups.append(SpinGroup(g["J"], g["parity"], channels, levels))
        for w in out["warnings"]:
            warnings.warn(f"spin group {g['jpi']}: {w}", stacklevel=2)
    target = entrance
    za = 1000.0 * target["Z2"] + round(target["M2"])
    awr = target["M2"] / M_NEUTRON_U
    record = LRF7(za, awr, za, 1.0, float(energy_range_eV[0]), float(energy_range_eV[1]), 1, 1, 4, 0,
                  tuple(endf_pairs), tuple(spin_groups), "")
    if not spin_groups:
        raise _refuse_session("the session has no level with a particle channel", "at least one level",
                              "add levels to the project")
    return record, Losses(tuple(lost), tuple(photon_widths))


def particle_pairs_text(record):
    """The pairs of a record as a short text (for messages)."""
    return listed(f"MT={p.mt}" for p in record.pairs)
