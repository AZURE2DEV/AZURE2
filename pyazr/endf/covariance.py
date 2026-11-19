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

The covariance written (the owner's answer of 2026-11-18 to decision 10 of endf.md section 11, which replaced HESSE
in a session of the formal basis): the fit's covariance C of the session's free parameters, in the fit's own basis
(Brune or Park) with the fit's own fixed parameters, carried through the transformation to MF=2's parameters,
C_32 = J C J^T.  C is HESSE (iminuit, errordef 1: a chi-squared; the inverse of half the Hessian) of the session's
objective at x, or the matrix a caller passes (``EngineFitResult.covariance`` of the engine's fit).  J is the
derivative of MF=2's ER and GAM in the session's free parameters at the boundary values B_c the file writes (BND):
a reader of the file varies ER and GAM with BND held, so J must hold B_c too.  Per J^pi group, J = D2^-1 D_b: D2 the
Jacobian of door 2 (``observed_levels_of_block``: formal levels at B_c -> Brune levels) at MF=2's levels, D_b that of
the group's Brune parameters in the session's (the identity in the Brune basis), both by central differences
(steps ``ENERGY_STEP``, ``AMPLITUDE_STEP``), then the file's units (ER: the lab factor of the incident pair x 1e6,
GAM: the amplitude scale of :mod:`pyazr.endf.resonances`).  So J C J^T taken back through door 2 is the fit's
covariance of its level parameters.  (Door 1's own derivative, with B_c = S_c(E_1) following the levels, is another
matrix: on T1 its J C J^T differs from this one by up to 0.93 sigma_i sigma_j, and read back with BND fixed it would
not reproduce the fit's covariance.)  Every other free parameter of the fit (norms, shifts, target parameters, the
photon widths) has a zero column of J: it is marginalised.  The transformation mixes the levels of a group, so a
formal parameter moves with every free parameter of its group: a resonance is listed when a parameter of it has a
non-zero row of J, and a parameter with a zero row (in a group none of whose parameters the fit varied, or the formal
energy of the lowest level, which equals its fixed Brune energy at B_c = S_c(E_1)) has a zero row and column, which
MF=1 names.
"""

import dataclasses
import os

import numpy as np

from ..errors import FormatRefused
from .records import RecordWriter
from .resonances import amplitude_scale, lab_factor

EV_PER_MEV = 1e6
#: Central-difference steps of the Jacobians of the doors: a level energy (MeV) by ENERGY_STEP x (|E| + 1), an
#: amplitude (MeV^1/2) by AMPLITUDE_STEP x max(|gamma|, 1e-3).  Door 2's iteration is taken to DOOR2_TOLERANCE_MEV,
#: so that its rounding (tolerance / step, 1e-7) stays below the step's truncation error.
ENERGY_STEP = 1e-5
AMPLITUDE_STEP = 1e-5
DOOR2_TOLERANCE_MEV = 1e-12


@dataclasses.dataclass(frozen=True)
class FormalCovariance:
    """The covariance of MF=2's parameters in the file's units: ``parameters`` [(group, level, column)] in the order
    of MF=32 (group and level positions in the LRF7 record, column 0 for ER, k >= 1 for GAM of channel k), their
    ``values`` (ER in eV, GAM in eV^1/2, as MF=2 writes them), ``matrix`` (eV^2, eV^3/2, eV) = J C J^T, ``fixed``
    (True where the row of J is zero: a zero row), ``groups`` (the spin groups of the record that have a listed
    resonance) and ``levels`` {group: [levels listed]}, ``objective_fit`` (the session's objective at x; NaN when the
    caller gave C), ``jacobian`` J (rows: ``parameters``, file units per the fit's unit; columns: the session's free
    parameters), ``fit_covariance`` C and ``source`` (in words)."""
    parameters: tuple
    values: np.ndarray
    matrix: np.ndarray
    fixed: tuple
    groups: tuple
    levels: dict
    objective_fit: float = float("nan")
    jacobian: np.ndarray = None
    fit_covariance: np.ndarray = None
    source: str = ""

    @property
    def sigma(self):
        return np.sqrt(np.clip(np.diag(self.matrix), 0.0, None))


def formal_covariance(session, x=None, record=None, covariance=None):
    """The covariance of the formal parameters of MF=2 at the fit vector ``x`` (default the session's start): the
    fit's covariance carried through the transformation, J C J^T (the module's text).  ``covariance``: C, the
    covariance of the session's free parameters in their order (``EngineFitResult.covariance``); None: HESSE
    (iminuit) of the session's objective at x.  ``record``: the LRF7 of ``from_session(session, x)`` (computed when
    None).  Refused (:class:`pyazr.FormatRefused`): HESSE asked of a session of the extrapolations or without
    iminuit, a C of the wrong shape or not finite, a HESSE that is not positive definite."""
    from .model import from_session
    if record is None:
        record, _ = from_session(session, x)
    x = session.parameters.x0 if x is None else np.asarray(x, float)
    n_free = int(np.size(x))
    objective_fit = float("nan")
    if covariance is None:
        C, objective_fit = _hesse(session, x)
        source = ("HESSE (iminuit, errordef 1) of the fit's objective in the session's basis at x, carried through the "
                  "transformation to MF=2's parameters (J C J^T)")
    else:
        C = np.asarray(covariance, float)
        if C.shape != (n_free, n_free) or not np.all(np.isfinite(C)):
            raise FormatRefused(what=f"a covariance of shape {C.shape} (finite: {bool(np.all(np.isfinite(C)))})",
                                allowed=f"the finite covariance of the session's {n_free} free parameters in their "
                                "order (EngineFitResult.covariance)",
                                fix="pass the fit's covariance, or covariance=None for HESSE at x")
        C = 0.5 * (C + C.T)
        source = ("the fit's covariance in the session's basis, carried through the transformation to MF=2's "
                  "parameters (J C J^T)")
    J, params, values = jacobian(session, x, record)
    return _assemble(record, J, params, values, C, objective_fit, source)


def _hesse(session, x):
    """C = HESSE (iminuit, errordef 1) of the session's objective at x, symmetrised, and the objective at x."""
    if session.mode != "data":
        raise FormatRefused(what="a session of the extrapolations has no objective", allowed="a session of the data "
                            "(the fit's), or the fit's covariance as covariance=",
                            fix="open the project without extrapolations=True")
    try:
        import iminuit
    except ImportError:
        raise FormatRefused(what="iminuit is not installed (HESSE of the fit's objective)",
                            allowed="pyazr with iminuit, or the fit's covariance as covariance=",
                            fix="install iminuit, pass EngineFitResult.covariance, or write without a covariance"
                            ) from None
    f = session.objective()
    objective_fit = float(f(x))
    m = iminuit.Minuit(f, np.asarray(x, float), grad=f.grad, name=list(f.names))
    m.hesse()
    C = np.array(m.covariance, float)
    C = 0.5 * (C + C.T)  # Minuit's matrix is symmetric to its rounding
    if not np.all(np.isfinite(C)) or np.any(np.linalg.eigvalsh(C) <= 0):
        raise FormatRefused(what="HESSE of the fit's objective gives no positive definite covariance",
                            allowed="a fit at a minimum with every free parameter determined",
                            fix="fit the project to its minimum first; fix the parameters the data do not determine")
    return C, objective_fit


def _door2(pairs, g, chans, E, gamma, B):
    """Door 2 (``observed_levels_of_block``) on one J^pi group: the Brune E and gamma of formal levels E, gamma with
    the boundary values B, the iteration taken to ``DOOR2_TOLERANCE_MEV``."""
    from .model import _engine
    out = _engine().observed_levels_of_block(pairs, int(round(2 * g["J"])), g["parity"], chans, list(E), list(gamma),
                                             list(B), DOOR2_TOLERANCE_MEV)
    return np.concatenate([np.asarray(out["E"], float), np.asarray(out["gamma"], float)])


def _central(f, v0, steps):
    """The Jacobian of the vector function f at v0 by central differences with the given steps, one column each."""
    cols = []
    for k, h in enumerate(steps):
        vp, vm = np.array(v0, float), np.array(v0, float)
        vp[k] += h
        vm[k] -= h
        cols.append((f(vp) - f(vm)) / (2.0 * h))
    return np.array(cols).T


def _steps(v, n_energies):
    return [ENERGY_STEP * (abs(t) + 1.0) if k < n_energies else AMPLITUDE_STEP * max(abs(t), 1e-3)
            for k, t in enumerate(v)]


def jacobian(session, x, record):
    """J = d(MF=2's ER and GAM, file units) / d(the session's free parameters) at x (the module's text): per J^pi
    group, the formal levels at the boundary values B_c the file writes, J = D2^-1 D_b, D2 the Jacobian of door 2 at
    MF=2's levels and B_c, D_b that of the group's Brune parameters in the session's (the identity for a Brune
    session; door 2 after door 1 for a Park session), both by central differences.  Returns J, the parameters
    [(group, level, column)] of every resonance of ``record`` in MF=32's order and their values (MF=2's).  Columns
    of free parameters that are no level energy or particle amplitude are zero."""
    from .model import _block_of, _full_values, _level_groups, formal_levels
    rows = [p.to_dict() for p in session.parameters]
    free_of = {i: p.free_index for i, p in enumerate(session.parameters) if not p.fixed}
    values = _full_values(session, x)
    pairs = [p.to_dict() for p in session.pairs]
    groups = _level_groups(rows, values)
    order = [gk for gk in sorted(groups) if groups[gk]["channels"]]
    if len(order) != len(record.groups):
        raise FormatRefused(what="the session's spin groups are not MF=2's", allowed="the record of this session",
                            fix="pass record=None (from_session of this session at x)")
    ent = record.pairs[0]
    e_scale = lab_factor(ent) * EV_PER_MEV
    g_scale = amplitude_scale(ent)
    params, file_values = [], []
    J = np.zeros((sum(len(sg.levels) * (len(sg.channels) + 1) for sg in record.groups), np.size(x)))
    for gi, gk in enumerate(order):
        g = groups[gk]
        block = _block_of(g)
        chans = block["chans"]
        n_lev, nch = len(block["E"]), len(block["chan_ids"])
        spin_group = record.groups[gi]
        first = len(params)
        for li, lv in enumerate(spin_group.levels):
            for col, v in enumerate([lv.energy] + list(lv.gam)):
                params.append((gi, li, col))
                file_values.append(v)
        inputs = list(block["E"]) + list(block["gam"])
        free_slots = [(k, free_of[r]) for k, r in enumerate(block["rows"]) if r is not None and r in free_of]
        if not free_slots:
            continue
        base = formal_levels(session, pairs, g, chans, block["E"], block["gam"])
        B = list(base["B"])
        formal = np.concatenate([np.asarray(base["E"], float), np.asarray(base["gamma"], float)])
        D2 = _central(lambda f: _door2(pairs, g, chans, f[:n_lev], f[n_lev:], B), formal, _steps(formal, n_lev))
        if session.basis == 2:  # Park: the group's Brune parameters as functions of the Park inputs
            def brune_of(v):
                out = formal_levels(session, pairs, g, chans, v[:n_lev], v[n_lev:])
                return _door2(pairs, g, chans, out["E"], out["gamma"], out["B"])
            Db = _central(brune_of, inputs, _steps(inputs, n_lev))
        else:
            Db = np.eye(len(inputs))
        Jg = np.linalg.solve(D2, Db)  # rows: formal E (levels), then formal gamma (level-major); columns: inputs
        for li in range(n_lev):
            a = first + li * (nch + 1)
            for k, col in free_slots:
                J[a, col] = Jg[li, k] * e_scale
                for c in range(nch):
                    J[a + 1 + c, col] = Jg[n_lev + li * nch + c, k] * g_scale
    return J, params, np.array(file_values, float)


def _assemble(record, J, params, values, C, objective_fit, source):
    """The listed resonances (a parameter with a non-zero row of J), their covariance J C J^T, symmetrised."""
    moved = np.any(J != 0.0, axis=1)
    listed_levels = sorted({(gi, li) for (gi, li, _), m in zip(params, moved) if m})
    keep = [k for k, (gi, li, _) in enumerate(params) if (gi, li) in set(listed_levels)]
    Jk = J[keep]
    cov = Jk @ C @ Jk.T
    cov = 0.5 * (cov + cov.T)
    levels_listed, groups_listed = {}, []
    for gi, li in listed_levels:
        levels_listed.setdefault(gi, []).append(li)
        if gi not in groups_listed:
            groups_listed.append(gi)
    return FormalCovariance(tuple(params[k] for k in keep), values[keep], cov, tuple(not moved[k] for k in keep),
                            tuple(groups_listed), levels_listed, objective_fit, Jk, C, source)


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
