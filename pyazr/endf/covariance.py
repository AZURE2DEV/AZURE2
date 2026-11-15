"""Covariances of the resonance parameters, MF=32/MT=151 for LRF=7 in the general form LCOMP=1 (ENDF-102, 2023 edition,
cited by printed page; docs/upgrade/endf.md section 13.13).

What a file holds (pp. 253-255, 258-259): HEAD ``ZA, AWR, 0, 0, NIS, 0``; per isotope a CONT ``ZAI, ABN, 0, LFW,
NER, 0``; per range a CONT ``EL, EH, LRU=1, LRF=7, NRO, NAPS``; then ``0, 0, 0, LCOMP=1, 0, ISR`` (p. 258), with
ISR=0 no record of radius uncertainties; ``AWRI, 0, 0, 0, NSRS, NLRS`` (NLRS = 0: "No long-range subsections are
allowed for LRF=7", p. 258); per short-range block a CONT ``0, 0, NJSX, 0, 0, 0``, per included spin group a LIST
``0, 0, NCH, NRB, 6*NX, NX`` of ER and the NCH values GAM of each included resonance (each resonance starting a new
line, as in MF=2, p. 74), and after the last group the upper triangle of the covariance of the NPARB = sum (NCH+1)
NRB parameters, row by row, a LIST ``0, 0, 0, 0, N, NPARB`` (p. 259).  The parameters are those of File 2 as listed
there (p. 253): ER in eV (lab), GAM the reduced-width amplitude in eV^1/2 (IFG=1, lab scale) with its sign.  A
spin group carries no J in MF=32: its resonances are identified by their ER and GAM against MF=2 ("for positive
identification", p. 256).

The covariance written (the PI's decision, endf.md section 11 decision 10): HESSE in a session of the formal basis
at the transformed best fit -- the session's fit vector written as a project (``Session.save_fit``) and opened with
``use_brune=False``, so that its parameters are the formal E and gamma of MF=2 and its objective is the fit's; the
second derivatives of that objective by iminuit's HESSE (the recipe of :mod:`pyazr.objective`), errordef 1 (a
chi-squared), the covariance the inverse of half the Hessian; the rows of MF=2's parameters taken out (the norms and
other nuisances of the fit, and the photon channels, are marginalised: their rows are dropped) and scaled to the
file's units (ER: the lab factor of the incident pair x 1e6; GAM: the amplitude scale of :mod:`pyazr.endf.
resonances`).  A resonance whose every parameter is fixed in the fit is not listed; a fixed parameter of a listed
resonance has a zero row and column (the format lists every parameter of a resonance it includes), and MF=1 says
which.
"""

import dataclasses
import os
import tempfile

import numpy as np

from ..errors import FormatRefused
from .records import RecordWriter
from .resonances import amplitude_scale, lab_factor

EV_PER_MEV = 1e6
#: The formal session's objective at its start must equal the fit's at x to this (relative): the same model, up to the
#: snapshot's round trip through the observed widths (measured 1.9e-7 on the fitted T2).
SAME_MODEL_TOLERANCE = 1e-5
#: The formal session's level values must equal MF=2's (door 1) to this (relative to the largest of each kind).
SAME_LEVELS_TOLERANCE = 1e-7


@dataclasses.dataclass(frozen=True)
class FormalCovariance:
    """The covariance of MF=2's parameters in the file's units: ``parameters`` [(group, level, column)] in the order
    of MF=32 (group and level positions in the LRF7 record, column 0 for ER, k >= 1 for GAM of channel k), their
    ``values`` (ER in eV, GAM in eV^1/2, as MF=2 writes them), ``matrix`` (eV^2, eV^3/2, eV), ``fixed`` (True where
    the fit did not vary the parameter: a zero row), ``groups`` (the spin groups of the record that have a listed
    resonance) and ``levels`` {group: [levels listed]}, the two objectives (``objective_fit`` of the session at x,
    ``objective_formal`` of the formal session at its start) and ``source`` (in words)."""
    parameters: tuple
    values: np.ndarray
    matrix: np.ndarray
    fixed: tuple
    groups: tuple
    levels: dict
    objective_fit: float = float("nan")
    objective_formal: float = float("nan")
    source: str = ""

    @property
    def sigma(self):
        return np.sqrt(np.clip(np.diag(self.matrix), 0.0, None))


def formal_covariance(session, x=None, record=None):
    """The covariance of the formal parameters of MF=2 at the fit vector ``x`` (default the session's start) from
    HESSE in a formal-basis session (the module's text).  ``record``: the LRF7 of ``from_session(session, x)``
    (computed when None).  Refused (:class:`pyazr.FormatRefused`): a session of the extrapolations, no iminuit, a
    formal session whose objective or levels are not the fit's, a Hessian that is not positive definite."""
    from .model import from_session
    if session.mode != "data":
        raise FormatRefused(what="a session of the extrapolations has no objective", allowed="a session of the data "
                            "(the fit's)", fix="open the project without extrapolations=True")
    try:
        import iminuit
    except ImportError:
        raise FormatRefused(what="iminuit is not installed (HESSE of the formal-basis session)",
                            allowed="pyazr with iminuit", fix="install iminuit, or write without a covariance") from None
    if record is None:
        record, _ = from_session(session, x)
    x = session.parameters.x0 if x is None else np.asarray(x, float)
    objective_fit = float(session.evaluate(x).objective)
    # the snapshot lies beside the project (its data files and covariance files are named relative to it) under a
    # name of its own, and is removed after the formal session has read it
    handle, snapshot = tempfile.mkstemp(prefix=".endf_covariance_", suffix=".azr",
                                        dir=os.path.dirname(os.path.abspath(session.file)))
    os.close(handle)
    try:
        session.save_fit(snapshot, x, param_sav=False, verify=False)
        return _hesse(session, record, snapshot, objective_fit)
    finally:
        os.remove(snapshot)


def _hesse(session, record, snapshot, objective_fit):
    """HESSE of the objective of the formal-basis session of ``snapshot`` (the module's text)."""
    import iminuit
    from ..session import Session
    from .engine_values import SWEEP_OPTIONS
    options = {k: v for k, v in session.options.items() if k in SWEEP_OPTIONS and k not in ("use_brune", "use_park")}
    formal = Session(snapshot, use_brune=False, **options)
    try:
        f = formal.objective()
        x0 = np.asarray(f.x0, float)
        objective_formal = float(f(x0))
        if abs(objective_formal - objective_fit) > SAME_MODEL_TOLERANCE * max(abs(objective_fit), 1.0):
            raise FormatRefused(what=f"the formal-basis session's objective {objective_formal!r} is not the fit's "
                                f"{objective_fit!r}", allowed="the same model in both bases",
                                fix="report the project to the AZURE2 developers")
        rows = [p.to_dict() for p in formal.parameters]
        free_of = {}
        k = 0
        for i, p in enumerate(formal.parameters):
            if not p.fixed:
                free_of[i] = k
                k += 1
        m = iminuit.Minuit(f, x0, grad=f.grad, name=list(f.names))
        m.hesse()
        cov_free = np.array(m.covariance, float)
        cov_free = 0.5 * (cov_free + cov_free.T)  # Minuit's matrix is symmetric to its rounding
        full_values = [float(x0[free_of[i]]) if i in free_of else float(p.start)
                       for i, p in enumerate(formal.parameters)]
    finally:
        formal.close()
    if not np.all(np.isfinite(cov_free)) or np.any(np.linalg.eigvalsh(0.5 * (cov_free + cov_free.T)) <= 0):
        raise FormatRefused(what="HESSE of the formal-basis session gives no positive definite covariance",
                            allowed="a fit at a minimum with every free parameter determined",
                            fix="fit the project to its minimum first; fix the parameters the data do not determine")
    entrance = next(p.to_dict() for p in session.pairs if p.is_entrance)
    s_ent = entrance["sep_energy"] + entrance["excitation"]
    return _assemble(record, s_ent, rows, full_values, free_of, cov_free, objective_fit, objective_formal)


def _assemble(record, s_ent, rows, values, free_of, cov_free, objective_fit, objective_formal):
    """The parameters of MF=2 in MF=32's order, their covariance in the file's units (the module's text)."""
    ent = record.pairs[0]
    e_scale = lab_factor(ent) * EV_PER_MEV
    g_scale = amplitude_scale(ent)
    # the formal session's rows by (jgroup, level) and (jgroup, level, channel), as from_session builds the record:
    # groups in jgroup order that have a particle channel, levels and channels in their id order
    groups = {}
    for i, r in enumerate(rows):
        if r["kind"] == "energy":
            groups.setdefault(r["jgroup"], {"levels": {}, "channels": {}})["levels"][r["level"]] = i
    for i, r in enumerate(rows):
        if r["kind"] == "width" and r.get("radiation_type") == "P":
            g = groups[r["jgroup"]]
            g["channels"][r["channel"]] = True
            g.setdefault("widths", {})[(r["level"], r["channel"])] = i
    order = [gk for gk in sorted(groups) if groups[gk]["channels"]]
    if len(order) != len(record.groups):
        raise FormatRefused(what="the formal-basis session's spin groups are not MF=2's", allowed="the same project",
                            fix="report the project to the AZURE2 developers")
    params, index, fixed, file_values, levels_listed, groups_listed = [], [], [], [], {}, []
    worst = 0.0
    for gi, gk in enumerate(order):
        g = groups[gk]
        lev_ids = sorted(g["levels"])
        ch_ids = sorted(g["channels"])
        spin_group = record.groups[gi]
        for li, lid in enumerate(lev_ids):
            cols = [g["levels"][lid]] + [g.get("widths", {}).get((lid, c)) for c in ch_ids]
            if not any(c is not None and c in free_of for c in cols):
                continue
            levels_listed.setdefault(gi, []).append(li)
            if gi not in groups_listed:
                groups_listed.append(gi)
            lv = spin_group.levels[li]
            written = [lv.energy] + list(lv.gam)
            for col, c in enumerate(cols):
                params.append((gi, li, col))
                index.append(None if c is None else free_of.get(c))
                fixed.append(c is None or c not in free_of)
                v = (values[c] - s_ent) * e_scale if col == 0 else (0.0 if c is None else values[c] * g_scale)
                file_values.append(v)
                # ER relative to 1 eV at least, GAM relative to its value
                worst = max(worst, abs(v - written[col]) / max(abs(written[col]), 1.0 if col == 0 else 1e-300))
    if worst > SAME_LEVELS_TOLERANCE:
        raise FormatRefused(what=f"the formal-basis session's levels differ from MF=2's by {worst:.2e}",
                            allowed="the levels of from_session (door 1)", fix="report the project to the AZURE2 "
                            "developers")
    n = len(params)
    cov = np.zeros((n, n))
    scale = np.array([e_scale if p[2] == 0 else g_scale for p in params])
    for a in range(n):
        for b in range(a, n):
            if index[a] is not None and index[b] is not None:
                cov[a, b] = cov[b, a] = cov_free[index[a], index[b]] * scale[a] * scale[b]
    return FormalCovariance(tuple(params), np.array(file_values), cov, tuple(fixed), tuple(groups_listed),
                            levels_listed, objective_fit, objective_formal,
                            "HESSE (iminuit, errordef 1) of the objective in a formal-basis session at the transformed "
                            "best fit")


def rows(record, covariance, number_form="widest"):
    """MF=32/MT=151 (LCOMP=1, one short-range block) of ``covariance`` (a :class:`FormalCovariance` of ``record``)
    as rows of 66 columns (the module's text, pp. 253-259)."""
    w = RecordWriter(number_form)
    w.cont(record.za, record.awr, 0, 0, 1, 0)
    w.cont(record.zai, record.abn, 0, 0, 1, 0)
    w.cont(record.el, record.eh, 1, 7, 0, record.naps)
    w.cont(0.0, 0.0, 0, 1, 0, 0)            # p. 258: 0.0, 0.0, 0, LCOMP=1, 0, ISR=0
    w.cont(record.awr, 0.0, 0, 0, 1, 0)     # AWRI, 0.0, 0, 0, NSRS=1, NLRS=0
    w.cont(0.0, 0.0, len(covariance.groups), 0, 0, 0)  # NJSX
    by_level = {}
    for k, (gi, li, col) in enumerate(covariance.parameters):
        by_level.setdefault((gi, li), [None] * (len(record.groups[gi].channels) + 1))[col] = covariance.values[k]
    for gi in covariance.groups:
        nch = len(record.groups[gi].channels)
        per_level = 6 * (nch // 6 + 1)  # each resonance starts a new line (as MF=2, p. 74)
        values = []
        for li in covariance.levels[gi]:
            row = list(by_level[(gi, li)])
            values.extend(row + [0.0] * (per_level - len(row)))
        w.list(0.0, 0.0, nch, len(covariance.levels[gi]), len(values) // 6, values)
    n = len(covariance.parameters)
    upper = [covariance.matrix[i, j] for i in range(n) for j in range(i, n)]
    w.list(0.0, 0.0, 0, 0, n, upper)
    return w.rows


@dataclasses.dataclass(frozen=True)
class ResonanceCovariance:
    """MF=32/MT=151 of a file as read: ``names`` of the parameters (``"<J^pi>#<k>.ER"`` and
    ``"<J^pi>#<k>.GAM<c>"``, k the resonance's position in its spin group of MF=2 counted from 1, c the channel's;
    the spin group found by the resonances' ER and GAM in MF=2), their ``values`` (eV, eV^1/2, as MF=2 lists them) and
    ``matrix`` (the file's units), ``groups`` (the MF=2 spin group of each LIST, by position)."""
    names: tuple
    values: np.ndarray
    matrix: np.ndarray
    groups: tuple

    @property
    def sigma(self):
        return np.sqrt(np.clip(np.diag(self.matrix), 0.0, None))

    @property
    def correlation(self):
        s = self.sigma
        with np.errstate(divide="ignore", invalid="ignore"):
            c = self.matrix / np.outer(s, s)
        return np.where(np.outer(s, s) > 0, c, 0.0)


def read_covariance(material, resonances=None):
    """MF=32/MT=151 of ``material`` (one isotope, one LRF=7 range, LCOMP=1; others refused with what, allowed, fix)
    as a :class:`ResonanceCovariance`; ``resonances``: the LRF7 of MF=2 (read when None) that names the parameters."""
    source = f"{os.path.basename(material.path)} MAT={material.mat} MF=32/MT=151"
    if (32, 151) not in material.sections:
        raise FormatRefused(what=f"{os.path.basename(material.path)} has no MF=32/MT=151",
                            allowed="a file with covariances of its resonance parameters",
                            fix="Material.summary() lists the sections")
    c = material.section(32, 151).cursor(material.path)
    head = c.cont()
    iso = c.cont()
    rng = c.cont()
    if head.n1 != 1 or iso.n1 != 1:
        raise FormatRefused(what=f"{source}: NIS={head.n1}, NER={iso.n1}", allowed="one isotope and one range",
                            fix="not read here")
    if rng.l1 != 1 or rng.l2 != 7 or rng.n1 != 0:
        raise FormatRefused(what=f"{source}: LRU={rng.l1} LRF={rng.l2} NRO={rng.n1}",
                            allowed="LRU=1, LRF=7, NRO=0 (p. 258)", fix="not read here")
    flags = c.cont()
    if flags.l2 != 1:
        raise FormatRefused(what=f"{source}: LCOMP={flags.l2}", allowed="LCOMP=1 (pp. 255, 258)",
                            fix="LCOMP=0 and 2 are not read here")
    if flags.n2 > 0:
        c.list()  # ISR > 0: the radius uncertainties (p. 258), not kept
    blocks = c.cont()
    if blocks.n2 != 0:
        raise FormatRefused(what=f"{source}: NLRS={blocks.n2}", allowed="NLRS=0 (p. 258: no long-range subsections "
                            "for LRF=7)", fix="the file breaks the format")
    if resonances is None:
        resonances = material.resonances()
    names, values, groups = [], [], []
    matrix = None
    for _ in range(blocks.n1):
        if matrix is not None:
            raise FormatRefused(what=f"{source}: {blocks.n1} short-range blocks", allowed="one block (NSRS=1)",
                                fix="not read here")
        njsx = c.cont().l1
        for _ in range(njsx):
            rec = c.list()
            nch, nrb = rec.l1, rec.l2
            per = len(rec.values) // max(nrb, 1)
            lists = [rec.values[k * per:k * per + nch + 1] for k in range(nrb)]
            gi, levels = _identify(resonances, nch, lists, source)
            groups.append(gi)
            g = resonances.groups[gi]
            jpi = f"{g.J:g}{'+' if g.parity > 0 else '-'}"
            for li, row in zip(levels, lists):
                names.append(f"{jpi}#{li + 1}.ER")
                names += [f"{jpi}#{li + 1}.GAM{k + 1}" for k in range(nch)]
                values.extend(row)
        tri = c.list()
        n = tri.n2
        if n != len(names) or len(tri.values) < n * (n + 1) // 2:
            raise FormatRefused(what=f"{source}: NPARB={n} where the spin groups list {len(names)} parameters",
                                allowed="NPARB = sum (NCH+1) NRB (p. 259)", fix="the file breaks the format")
        matrix = np.zeros((n, n))
        k = 0
        for i in range(n):
            for j in range(i, n):
                matrix[i, j] = matrix[j, i] = tri.values[k]
                k += 1
    return ResonanceCovariance(tuple(names), np.array(values, float),
                               matrix if matrix is not None else np.zeros((0, 0)), tuple(groups))


def _identify(resonances, nch, lists, source):
    """The MF=2 spin group and the positions of the listed resonances in it: the group with NCH channels whose
    resonances hold every listed ER (relative 1e-9)."""
    for gi, g in enumerate(resonances.groups):
        if len(g.channels) != nch:
            continue
        positions = []
        for row in lists:
            hit = [li for li, lv in enumerate(g.levels)
                   if abs(lv.energy - row[0]) <= 1e-9 * max(abs(row[0]), 1.0) and
                   all(abs(a - b) <= 1e-7 * max(abs(a), 1e-12) for a, b in zip(lv.gam, row[1:]))]
            if not hit:
                break
            positions.append(hit[0])
        else:
            return gi, positions
    raise FormatRefused(what=f"{source}: a spin group of {nch} channels whose resonances are not those of MF=2",
                        allowed="resonances listed as in MF=2 (p. 256: for positive identification)",
                        fix="the file's MF=32 does not belong to its MF=2")
