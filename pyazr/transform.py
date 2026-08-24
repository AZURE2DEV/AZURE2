"""Reduced-width amplitudes -> observed partial widths: a port of ``TransformOut``.

AZURE2 fits *reduced-width amplitudes* ``gamma`` (MeV^1/2) and reports *observed
partial widths* ``Gamma`` (eV).  The two are not related by ``Gamma = 2 gamma^2
P_c``: ``CNuc::TransformOut`` (``src/CNuc.cpp``) divides by the level-shift
normalisation,

    Gamma_c = 2 gamma_c^2 P_c(E_lambda) / (1 + sum_c' gamma_c'^2 dS_c'/dE|E_lambda)

and, when the Brune parameterisation is *not* in use, first rotates the level
matrix so the boundary condition sits at each level's own energy.  Anyone
post-processing an MCMC chain, a covariance sample, or a hand-edited parameter
vector has to redo that transformation, and the denominator is not a detail: for
a strongly coupled level such as the 3/2+ state in d+t it is ~34, so dropping it
overstates the width by a factor of 34.  Push ``gamma`` high enough and
``Gamma`` *saturates* at ``2 P_c / (dS_c/dE)`` -- arbitrarily large reduced
widths map onto a finite physical width, which is exactly why an unconstrained
sampler can wander off to ``gamma ~ 1e59`` without the fit noticing.

This module reproduces that transformation in Python, at AZURE2's own constants,
channel radii and boundary conditions, with no running instance required::

    from pyazr import AzrModel
    from pyazr.transform import levels_from_azr, transform_out

    levels = levels_from_azr(AzrModel.from_file("3H+d.azr"),
                             gammas={(1, 1): [0.0, 5.3032255, 0.0, 0.780586]})
    for lv in transform_out(levels, brune=True):
        print(lv.table())

For a whole MCMC chain, :class:`LevelWidthTransformer` tabulates ``P_c(E)`` and
``dS_c/dE(E)`` once on a spline and then maps hundreds of thousands of samples
in one vectorised call.

What is reproduced, and what is not
-----------------------------------
* Coulomb penetrability and shift function for open particle channels, and the
  Whittaker-function forms below threshold, at AZURE2's ``Constants.h`` values.
* The iterative level-matrix diagonalisation of the non-Brune branch
  (``brune=False``), including its ascending-eigenvalue level ordering.
* Photon channels: ``P = (E_gamma/hbarc)^(2L+1)``, plus the bound-state M1/E2
  moment special cases.
* *External* (channel-capture) widths are **not** computed here.  AZURE2 adds
  ``g_ext`` from ``CalcExternalWidth`` before squaring; pass it yourself via
  ``Channel.external_gamma`` if you need it.  With ``--ignore-externals``, or
  for any model without external capture, it is zero and this is exact.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Dict, List, Optional, Sequence, Tuple

import mpmath as _mp

# Verbatim from include/Constants.h -- do not "improve" these to CODATA values,
# the point is to agree with AZURE2 digit for digit.
HBARC = 197.32696310
UCONV = 931.4940880
FSTRUC = 1.00 / 137.0359996790
NUCLEAR_MAGNETON = 0.105155

# Energy step for dS/dE.  AZURE2 hands gsl_deriv_central h = 1e-6 MeV, which is
# small enough that the difference of two shift-function values is dominated by
# their last few digits: probing a level of 7Be where dS/dE = 1.128e-3, AZURE2's
# derivative comes out 1.7% high (the penetrability at the same point agrees to
# 2.5e-10, so it is the differencing, not the Coulomb functions).  1e-3 MeV is
# converged to ~9 digits on both the open and sub-threshold branches.  Set this
# to 1e-6 to reproduce AZURE2's numbers bit for bit rather than correctly.
SHIFT_DERIVATIVE_STEP = 1e-3

__all__ = [
    "Channel", "Level", "TransformedChannel", "TransformedLevel",
    "transform_out", "partial_widths", "LevelWidthTransformer",
    "levels_from_azr", "levels_from_scheme",
    "sommerfeld_eta", "coulomb_rho", "coulomb_waves",
    "penetrability", "shift", "shift_derivative", "whittaker",
]


# ---------------------------------------------------------------------------
# Coulomb / Whittaker functions -- CoulFunc.cpp, ShftFunc.cpp, WhitFunc.h
# ---------------------------------------------------------------------------

def sommerfeld_eta(z1: int, z2: int, red_mass: float, energy: float) -> float:
    """Sommerfeld parameter, ``CoulFunc::operator()``."""
    return math.sqrt(UCONV / 2.0) * FSTRUC * z1 * z2 * math.sqrt(red_mass / energy)


def coulomb_rho(red_mass: float, radius: float, energy: float) -> float:
    """Dimensionless radius ``k*a``, ``CoulFunc::operator()``."""
    return math.sqrt(2.0 * UCONV) / HBARC * radius * math.sqrt(red_mass * energy)


@lru_cache(maxsize=200_000)
def coulomb_waves(L: int, eta: float, rho: float) -> Tuple[float, float, float, float]:
    """``(F, dF/drho, G, dG/drho)`` -- what GSL's ``coulomb_wave_FG_e`` returns.

    mpmath has no derivative for these, so the derivatives come from the
    standard recurrence (DLMF 33.4.1), ``X'_L = S_{L+1} X_L - R_{L+1} X_{L+1}``,
    which needs only one extra order and is exact rather than differenced.
    """
    F = float(_mp.coulombf(L, eta, rho))
    G = float(_mp.coulombg(L, eta, rho))
    F1 = float(_mp.coulombf(L + 1, eta, rho))
    G1 = float(_mp.coulombg(L + 1, eta, rho))
    lp = L + 1.0
    R = math.sqrt(1.0 + (eta / lp) ** 2)
    S = lp / rho + eta / lp
    return F, S * F - R * F1, G, S * G - R * G1


def penetrability(L: int, z1: int, z2: int, red_mass: float, radius: float,
                  energy: float) -> float:
    """``P_L = rho / (F^2 + G^2)`` -- ``CoulFunc::Penetrability``."""
    eta = sommerfeld_eta(z1, z2, red_mass, energy)
    rho = coulomb_rho(red_mass, radius, energy)
    F, _, G, _ = coulomb_waves(L, eta, rho)
    return rho / (F * F + G * G)


def shift(L: int, z1: int, z2: int, red_mass: float, radius: float,
          energy: float) -> float:
    """``S_L = rho (F F' + G G') / (F^2 + G^2)`` -- ``CoulFunc::PEShift``."""
    eta = sommerfeld_eta(z1, z2, red_mass, energy)
    rho = coulomb_rho(red_mass, radius, energy)
    F, dF, G, dG = coulomb_waves(L, eta, rho)
    if F * F == 0.0 and F * dF == 0.0:
        return rho * (dG / G)
    return rho / (F * F + G * G) * (F * dF + G * dG)


def whittaker(L: int, z1: int, z2: int, red_mass: float, radius: float,
              binding_energy: float) -> float:
    """``W_{k,m}(z)`` for a sub-threshold channel -- ``WhitFunc::operator()``.

    ``binding_energy`` is positive (AZURE2 passes ``|E - S - Ex|``).  Evaluated
    as ``exp(-z/2) z^(m+1/2) U(m-k+1/2, 1+2m, z)``, the same expression AZURE2
    builds on GSL's ``hyperg_U``.  Just below a threshold ``U`` underflows a
    double, so mpmath takes over there.
    """
    k = -math.sqrt(UCONV / 2.0) * FSTRUC * z1 * z2 * math.sqrt(red_mass / binding_energy)
    m = L + 0.5
    z = 2.0 * math.sqrt(2.0 * UCONV) / HBARC * radius * math.sqrt(red_mass * binding_energy)

    try:
        from scipy.special import hyperu
        u = float(hyperu(m - k + 0.5, 1.0 + 2.0 * m, z))
        w = math.exp(-z / 2.0) * z ** (m + 0.5) * u
        if w != 0.0 and math.isfinite(w):
            return w
    except ImportError:
        pass
    return float(_mp.whitw(k, m, z))


def _central_diff(f, x, h):
    """5-point central derivative, standing in for ``gsl_deriv_central``."""
    return (f(x - 2 * h) - 8 * f(x - h) + 8 * f(x + h) - f(x + 2 * h)) / (12 * h)


def shift_derivative(L: int, z1: int, z2: int, red_mass: float, radius: float,
                     energy: float, h: float = SHIFT_DERIVATIVE_STEP) -> float:
    """``dS_L/dE`` above threshold -- ``CoulFunc::PEShift_dE``.

    See :data:`SHIFT_DERIVATIVE_STEP` for why ``h`` is not AZURE2's 1e-6.
    """
    return _central_diff(
        lambda e: shift(L, z1, z2, red_mass, radius, e), energy, h)


# ---------------------------------------------------------------------------
# Model objects
# ---------------------------------------------------------------------------

@dataclass
class Channel:
    """One channel of one level: geometry, kinematics and a reduced width.

    The fields mirror what ``AChannel`` + its ``PPair`` supply, so a channel can
    be built from a ``.azr`` file (:func:`levels_from_azr`), from a running
    instance (:func:`levels_from_scheme`), or by hand.
    """

    L: int
    gamma: float = 0.0                    # reduced-width amplitude, MeV^1/2
    radiation_type: str = "P"             # 'P' particle, 'E'/'M' photon, 'F'/'G' beta
    S: Optional[float] = None             # channel spin, label only
    pair: Optional[int] = None            # pair number, label only

    # pair physics
    Z1: int = 0
    Z2: int = 0
    M1: float = 0.0                       # amu
    M2: float = 0.0
    channel_radius: float = 0.0           # fm
    sep_energy: float = 0.0               # MeV
    excitation: float = 0.0               # MeV, of the pair's residual nucleus

    # only needed for the bound-state M1 / E2 moment special case
    pair_J2: Optional[float] = None
    pair_parity2: Optional[int] = None

    external_gamma: complex = 0j          # g_ext; AZURE2 adds it before squaring
    boundary_condition: Optional[float] = None   # None -> S_L at the first level

    # P_c, S_c and dS_c/dE depend only on the level energy, never on gamma, and
    # each costs several Coulomb-wave evaluations -- worth remembering.
    _cache: dict = field(default_factory=dict, repr=False, compare=False)

    @property
    def red_mass(self) -> float:
        """Reduced mass of the channel's pair, in u."""
        return self.M1 * self.M2 / (self.M1 + self.M2)

    @property
    def threshold(self) -> float:
        """Excitation energy at which this channel opens (MeV)."""
        return self.sep_energy + self.excitation

    @property
    def is_particle(self) -> bool:
        """Is this a particle channel?"""
        return self.radiation_type == "P"

    @property
    def is_photon(self) -> bool:
        """Is this a photon channel?"""
        return self.radiation_type in ("E", "M")

    def local_energy(self, level_energy: float) -> float:
        """Channel energy in the c.m. system: ``E_level - S - Ex``."""
        return level_energy - self.threshold

    def is_open(self, level_energy: float) -> bool:
        """Is the channel open at the level's energy, i.e. above its threshold?"""
        return self.local_energy(level_energy) >= 0.0

    # -- the three functions TransformOut needs -----------------------------

    def _memo(self, tag: str, level_energy: float, compute):
        key = (tag, level_energy)
        if key not in self._cache:
            self._cache[key] = compute(level_energy)
        return self._cache[key]

    def penetrability(self, level_energy: float) -> float:
        """``P_c`` at this level energy.

        Below threshold AZURE2 substitutes ``mu a / (hbar^2 W^2)``, which turns
        the width into a squared ANC rather than a width -- see
        :attr:`TransformedChannel.unit`.
        """
        return self._memo("P", level_energy, self._penetrability)

    def _penetrability(self, level_energy: float) -> float:
        e = self.local_energy(level_energy)
        if self.is_particle:
            if e < 0.0:
                w = whittaker(self.L, self.Z1, self.Z2, self.red_mass,
                              self.channel_radius, abs(e))
                return (self.red_mass * self.channel_radius * UCONV
                        / HBARC ** 2 / w ** 2)
            return penetrability(self.L, self.Z1, self.Z2, self.red_mass,
                                 self.channel_radius, e)
        if self.is_photon:
            return (abs(e) / HBARC) ** (2 * self.L + 1)
        return 1.0

    def shift(self, level_energy: float) -> float:
        """``S_c`` at this level energy (particle channels only)."""
        return self._memo("S", level_energy, self._shift)

    def _shift(self, level_energy: float) -> float:
        e = self.local_energy(level_energy)
        if not self.is_particle:
            return 0.0
        if e < 0.0:
            return self._whittaker_shift(level_energy)
        return shift(self.L, self.Z1, self.Z2, self.red_mass,
                     self.channel_radius, e)

    def shift_derivative(self, level_energy: float) -> float:
        """``dS_c/dE`` -- the quantity in the level-shift normalisation."""
        return self._memo("dS", level_energy, self._shift_derivative)

    def _shift_derivative(self, level_energy: float) -> float:
        if not self.is_particle:
            return 0.0
        h = SHIFT_DERIVATIVE_STEP
        if self.local_energy(level_energy) < 0.0:
            # ShftFunc::EnergyDerivative
            return _central_diff(self._whittaker_shift, level_energy, h)
        return shift_derivative(self.L, self.Z1, self.Z2, self.red_mass,
                                self.channel_radius,
                                self.local_energy(level_energy), h)

    def _whittaker_shift(self, level_energy: float) -> float:
        """``r W'(r)/W(r)`` -- ``ShftFunc::operator()`` (GSL radius step 1e-4)."""
        eb = abs(level_energy - self.threshold)
        r = self.channel_radius

        def w(x):
            """Whittaker function used to normalize a closed channel."""
            return whittaker(self.L, self.Z1, self.Z2, self.red_mass, x, eb)

        w0 = w(r)
        if w0 == 0.0:
            return 0.0
        return r * _central_diff(w, r, 1e-4) / w0


@dataclass
class Level:
    """One R-matrix level: an energy plus the channels of its J-group.

    ``energy`` is the level's excitation energy in the compound system, the same
    number the ``.azr`` ``<levels>`` block carries.  Levels sharing a
    ``jgroup`` mix in the non-Brune transformation and must list the *same*
    channels in the same order, as they do inside AZURE2.
    """

    energy: float
    channels: List[Channel]
    J: Optional[float] = None
    parity: Optional[int] = None
    jgroup: int = 1
    level: int = 1

    @property
    def gammas(self) -> List[float]:
        """The channel reduced-width amplitudes of the level."""
        return [c.gamma for c in self.channels]

    @property
    def jpi(self) -> str:
        """J^pi as text."""
        if self.J is None:
            return "?"
        j = int(self.J) if float(self.J).is_integer() else f"{int(round(2 * self.J))}/2"
        return f"{j}{'+' if (self.parity or 1) > 0 else '-'}"


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

@dataclass
class TransformedChannel:
    """One channel after the transformation: the number AZURE2 prints."""

    channel: Channel
    gamma: float                # transformed reduced-width amplitude, MeV^1/2
    big_gamma: float            # AZURE2's "bigGamma", MeV for an open channel
    penetrability: float
    is_open: bool

    @property
    def width_eV(self) -> Optional[float]:
        """Partial width in eV, or ``None`` for a closed / non-width channel."""
        if not self.is_open or self.channel.radiation_type in ("F", "G"):
            return None
        return abs(self.big_gamma) * 1e6

    @property
    def width_keV(self) -> Optional[float]:
        """Partial width in keV."""
        w = self.width_eV
        return None if w is None else w * 1e-3

    @property
    def anc(self) -> Optional[float]:
        """Asymptotic normalisation coefficient (fm^-1/2) for a closed channel."""
        if self.is_open or not self.channel.is_particle:
            return None
        return math.sqrt(abs(self.big_gamma))

    @property
    def unit(self) -> str:
        """Unit of the transformed value: eV for a width, fm^-1/2 for an ANC."""
        if self.channel.radiation_type in ("F", "G"):
            return ""
        return "eV" if self.is_open else "fm^-1/2"

    @property
    def value(self) -> float:
        """Whichever of width / ANC / beta strength applies, in :attr:`unit`."""
        if self.channel.radiation_type in ("F", "G"):
            return self.big_gamma
        return self.width_eV if self.is_open else self.anc

    @property
    def label(self) -> str:
        """The channel as text, with its pair, L and S."""
        c = self.channel
        if c.is_photon:
            return f"{c.radiation_type}{c.L}->pair{c.pair}"
        s = "" if c.S is None else f" s={c.S:g}"
        return f"pair{c.pair} l={c.L}{s}"


@dataclass
class TransformedLevel:
    """A level after the transformation, with the shared normalisation exposed."""

    energy: float               # transformed level energy (MeV)
    channels: List[TransformedChannel]
    norm_sum: float             # sum_c gamma_c^2 dS_c/dE  (the "1 +" is added below)
    iterations: int = 0
    J: Optional[float] = None
    parity: Optional[int] = None
    jgroup: int = 1
    level: int = 1

    @property
    def normalization(self) -> float:
        """``1 + sum_c gamma_c^2 dS_c/dE`` -- the divisor.  1 means "no effect"."""
        return 1.0 + self.norm_sum

    @property
    def total_width_eV(self) -> float:
        """Sum of the level's partial widths, in eV."""
        return sum(c.width_eV or 0.0 for c in self.channels)

    def table(self) -> str:
        """The same information ``output/parameters.out`` prints, as text."""
        jpi = "?" if self.J is None else (
            f"{int(self.J) if float(self.J).is_integer() else f'{int(round(2*self.J))}/2'}"
            f"{'+' if (self.parity or 1) > 0 else '-'}")
        head = (f"J = {jpi}  E_level = {self.energy:9.4f} MeV  "
                f"ITERATIONS = {self.iterations:3d}  "
                f"1+sum(g^2 dS/dE) = {self.normalization:.6g}")
        rows = []
        for tc in self.channels:
            val = tc.value
            rows.append(f"  {tc.label:<20}  "
                        f"{'G' if tc.is_open else 'C'} = {val:14.6g} {tc.unit:<8}  "
                        f"g_int = {tc.gamma:14.6g} MeV^(1/2)")
        return "\n".join([head] + rows)

    def __str__(self):
        return self.table()


# ---------------------------------------------------------------------------
# The transformation itself -- CNuc::TransformOut
# ---------------------------------------------------------------------------

def _default_boundaries(group: Sequence[Level]) -> List[float]:
    """``CNuc::CalcBoundaryConditions``: ``B_c = S_c(E)`` of the *first* level."""
    first = group[0]
    boundaries = []
    for c in first.channels:
        boundaries.append(c.shift(first.energy) if c.is_particle else 0.0)
    # Photon channels inherit channel 1's boundary condition.
    for i, c in enumerate(first.channels):
        if not c.is_particle:
            boundaries[i] = boundaries[0]
    # An explicit per-channel override wins, as AChannel's stored value would.
    for i, c in enumerate(first.channels):
        if c.boundary_condition is not None:
            boundaries[i] = c.boundary_condition
    return boundaries


def _rotate_group(group: Sequence[Level], this_level: int,
                  tolerance: float, max_iterations: int):
    """Non-Brune branch: iterate ``B_c -> S_c(E_lambda)`` for one level.

    Mirrors the ``while(iteration<=maxIterations&&!done)`` loop of
    ``TransformOut``, including GSL's ascending-eigenvalue ordering, which is
    what keeps ``this_level`` pointing at the same level across iterations.
    """
    import numpy as np

    channels = group[0].channels
    temp_e = [lv.energy for lv in group]
    temp_gamma = [[c.gamma for c in lv.channels] for lv in group]
    temp_boundary = _default_boundaries(group)

    n = len(group)
    iteration = 1
    done = False
    while iteration <= max_iterations and not done:
        boundary_diff = []
        for ch, c in enumerate(channels):
            if c.is_particle:
                new_b = c.shift(temp_e[this_level])
                boundary_diff.append(new_b - temp_boundary[ch])
                temp_boundary[ch] = new_b
            else:
                boundary_diff.append(boundary_diff[0] if boundary_diff else 0.0)

        cmat = np.zeros((n, n))
        for mu in range(n):
            for mup in range(n):
                chan_sum = sum(
                    boundary_diff[ch] * temp_gamma[mu][ch] * temp_gamma[mup][ch]
                    for ch, c in enumerate(channels) if c.is_particle)
                cmat[mu][mup] = (temp_e[mu] - chan_sum) if mu == mup else -chan_sum

        evals, evecs = np.linalg.eigh(cmat)          # ascending, like GSL's sort
        if abs(evals[this_level] - temp_e[this_level]) <= tolerance:
            done = True

        new_gamma = [[sum(evecs[mup][mu] * temp_gamma[mup][ch] for mup in range(n))
                      for ch in range(len(channels))] for mu in range(n)]
        temp_e = list(evals)
        temp_gamma = new_gamma

        if not done:
            if iteration == max_iterations:
                # AZURE2 warns and falls back to the untransformed parameters.
                temp_e[this_level] = group[this_level].energy
                temp_gamma[this_level] = [c.gamma for c in group[this_level].channels]
                break
            iteration += 1

    return temp_e[this_level], temp_gamma[this_level], iteration


def _normalize(level_energy: float, gammas: Sequence[float],
               channels: Sequence[Channel], input_energy: float,
               J: Optional[float], parity: Optional[int]):
    """The block every branch of ``TransformOut`` ends in.

    Builds ``normSum = sum_c gamma_c^2 dS_c/dE`` over particle channels, the
    per-channel penetrability, and ``bigGamma``.
    """
    norm_sum = 0.0
    penes = []
    for c, g in zip(channels, gammas):
        if c.is_particle:
            norm_sum += c.shift_derivative(level_energy) * g ** 2
            penes.append(c.penetrability(level_energy))
        elif c.is_photon:
            penes.append(_photon_penetrability(c, level_energy, input_energy,
                                               J, parity))
        else:
            penes.append(1.0)

    out = []
    for c, g, p in zip(channels, gammas, penes):
        total = complex(g) + c.external_gamma
        sign = -1.0 if total.real < 0.0 else 1.0
        if c.radiation_type in ("F", "G"):
            big = total.real
        else:
            big = sign * 2.0 * abs(total) ** 2 * p / (1.0 + norm_sum)
        out.append(TransformedChannel(channel=c, gamma=g, big_gamma=big,
                                      penetrability=p,
                                      is_open=c.is_open(level_energy)))
    return out, norm_sum


def _photon_penetrability(c: Channel, level_energy: float, input_energy: float,
                          J: Optional[float], parity: Optional[int]) -> float:
    """Photon channel ``P``; the M1/E2 branches are static moments, not widths."""
    e = c.local_energy(level_energy)
    is_ground_transition = (
        c.pair_J2 is not None and J is not None
        and abs(input_energy - c.excitation) < 1e-3
        and J == c.pair_J2 and parity == c.pair_parity2)
    if is_ground_transition:
        pene = 1e-10
        if c.radiation_type == "M" and c.L == 1:
            pene = 3.0 * J / 4.0 / (J + 1.0) / NUCLEAR_MAGNETON ** 2
        elif c.radiation_type == "E" and c.L == 2:
            pene = 60.0 * J * (2.0 * J - 1.0) / (J + 1.0) / (2.0 * J + 3.0)
        if int(round(2 * J)) % 2 != 0:
            pene *= -1.0
        return pene
    return (abs(e) / HBARC) ** (2 * c.L + 1)


def transform_out(levels: Sequence[Level], *, brune: bool = True,
                  tolerance: float = 1e-6,
                  max_iterations: int = 1000) -> List[TransformedLevel]:
    """Transform reduced-width amplitudes into AZURE2's physical parameters.

    ``levels`` may span several J-groups; levels sharing a ``jgroup`` are
    transformed together (they mix when ``brune=False``).  Set ``brune=True``
    when the fit was run with ``--use-brune`` / the GUI's Brune checkbox, in
    which case the level energies and reduced widths pass through untouched and
    only the level-shift normalisation is applied -- exactly the ``else`` branch
    of ``CNuc::TransformOut``.

    Returns one :class:`TransformedLevel` per input level, in input order.
    """
    groups: Dict[int, List[Level]] = {}
    for lv in levels:
        groups.setdefault(lv.jgroup, []).append(lv)

    results: Dict[int, TransformedLevel] = {}
    for jg, group in groups.items():
        nchan = len(group[0].channels)
        if any(len(lv.channels) != nchan for lv in group):
            raise ValueError(
                f"J-group {jg}: levels of a J-group share one channel list, "
                f"got {[len(lv.channels) for lv in group]} channels")

        for i, lv in enumerate(group):
            if brune:
                # The Brune parameters already sit at B_c = S_c(E_lambda); only
                # the normalisation below applies.
                energy, gammas, iters = lv.energy, lv.gammas, 0
            else:
                energy, gammas, iters = _rotate_group(group, i, tolerance,
                                                      max_iterations)

            chans, norm_sum = _normalize(energy, gammas, lv.channels,
                                         lv.energy, lv.J, lv.parity)
            results[id(lv)] = TransformedLevel(
                energy=energy, channels=chans, norm_sum=norm_sum,
                iterations=iters, J=lv.J, parity=lv.parity,
                jgroup=lv.jgroup, level=lv.level)

    return [results[id(lv)] for lv in levels]


def partial_widths(level: Level, gammas: Optional[Sequence[float]] = None, *,
                   brune: bool = True, unit: str = "eV") -> List[Optional[float]]:
    """Convenience: the partial widths of one level, one number per channel.

    ``gammas`` overrides the channels' own values, so a chain sample can be fed
    straight in.  Closed channels come back as ``None`` (they carry an ANC, not
    a width).  ``unit`` is ``'eV'``, ``'keV'`` or ``'MeV'``.
    """
    if gammas is not None:
        if len(gammas) != len(level.channels):
            raise ValueError(f"expected {len(level.channels)} gammas, "
                             f"got {len(gammas)}")
        for c, g in zip(level.channels, gammas):
            c.gamma = float(g)
    scale = {"eV": 1.0, "keV": 1e-3, "MeV": 1e-6}[unit]
    out = transform_out([level], brune=brune)[0]
    return [None if c.width_eV is None else c.width_eV * scale for c in out.channels]


# ---------------------------------------------------------------------------
# Vectorised path for chains
# ---------------------------------------------------------------------------

class LevelWidthTransformer:
    """Map many ``(E_level, gammas)`` samples onto partial widths, fast.

    ``P_c`` and ``dS_c/dE`` depend only on the level energy, and a Coulomb-wave
    evaluation costs milliseconds -- far too slow to repeat per MCMC sample.
    This tabulates both on a cubic spline across the energy range you hand it
    (or a window around a single energy) and then evaluates the transformation
    with numpy, so a 500k-sample chain takes well under a second::

        tr = LevelWidthTransformer(level, energies=chain["param0"])
        widths = tr.widths(chain["param0"], chain[["param1", "param2", ...]])

    Brune only: with ``brune=False`` levels of a J-group mix and the
    transformation is no longer a per-sample closed form -- loop over
    :func:`transform_out` instead.
    """

    def __init__(self, level: Level, energies=None, *, npoints: int = 33,
                 pad: float = 1e-3):
        import numpy as np
        from scipy.interpolate import CubicSpline

        self.level = level
        self.channels = level.channels

        if energies is None:
            lo, hi = level.energy - pad, level.energy + pad
        else:
            e = np.asarray(energies, float)
            lo, hi = float(e.min()), float(e.max())
            if hi - lo < 2 * pad:
                mid = 0.5 * (lo + hi)
                lo, hi = mid - pad, mid + pad
        self.energy_range = (lo, hi)

        grid = np.linspace(lo, hi, npoints)
        pene = np.empty((npoints, len(self.channels)))
        dsde = np.empty((npoints, len(self.channels)))
        for i, e in enumerate(grid):
            for k, c in enumerate(self.channels):
                pene[i, k] = c.penetrability(e)
                dsde[i, k] = c.shift_derivative(e)

        self._pene = CubicSpline(grid, pene, axis=0)
        self._dsde = CubicSpline(grid, dsde, axis=0)
        self._open = np.array([c.is_open(0.5 * (lo + hi)) for c in self.channels])
        self._particle = np.array([c.is_particle for c in self.channels])

    def normalization(self, energy, gammas):
        """``1 + sum_c gamma_c^2 dS_c/dE``, one value per sample."""
        import numpy as np
        energy = np.asarray(energy, float)
        gammas = np.asarray(gammas, float)
        dsde = self._dsde(energy) * self._particle
        return 1.0 + np.sum(dsde * gammas ** 2, axis=-1)

    def widths(self, energy, gammas, unit: str = "eV"):
        """Partial widths, shape ``(nsamples, nchannels)``.

        Closed channels come back as ``nan``; use :meth:`ancs` for those.
        """
        import numpy as np
        energy = np.asarray(energy, float)
        gammas = np.asarray(gammas, float)
        scale = {"eV": 1e6, "keV": 1e3, "MeV": 1.0}[unit]
        big = (2.0 * gammas ** 2 * self._pene(energy)
               / self.normalization(energy, gammas)[..., None])
        big = np.where(self._open, big, np.nan)
        return big * scale

    def ancs(self, energy, gammas):
        """ANCs (fm^-1/2) for the closed channels; ``nan`` for open ones."""
        import numpy as np
        energy = np.asarray(energy, float)
        gammas = np.asarray(gammas, float)
        big = (2.0 * gammas ** 2 * self._pene(energy)
               / self.normalization(energy, gammas)[..., None])
        return np.where(self._open, np.nan, np.sqrt(np.abs(big)))


# ---------------------------------------------------------------------------
# Adapters
# ---------------------------------------------------------------------------

def _radiation_type(azr_channel) -> str:
    """``AChannel::AChannel(NucLine, int)``: pType -> 'P' / 'E' / 'M' / 'F' / 'G'."""
    ptype = azr_channel.ptype
    if ptype == 0:
        return "P"
    if ptype == 10:
        return ("E" if azr_channel.levelPi * azr_channel.parity2
                == (-1) ** azr_channel.L else "M")
    if ptype == 20:
        return "F" if azr_channel.L == 0 else "G"
    raise ValueError(f"unknown pType {ptype}")


def levels_from_azr(model, gammas: Optional[Dict[Tuple[int, int], Sequence[float]]] = None,
                    *, active_only: bool = True) -> List[Level]:
    """Build :class:`Level` objects from a :class:`pyazr.AzrModel` (or a path).

    The ``.azr`` ``<levels>`` block carries the *physical* widths, not reduced
    ones, so supply the reduced-width amplitudes yourself -- from
    ``output/param.par``, an MCMC sample, whatever.  ``gammas`` is keyed by
    ``(jgroup, level)`` where the J-group is numbered in first-seen J^pi order,
    matching AZURE2; omit it to get channels with ``gamma = 0`` that you fill in
    later.
    """
    if isinstance(model, str):
        from .azrfile import AzrModel
        model = AzrModel.from_file(model)

    jgroup_of: Dict[Tuple[float, int], int] = {}
    counter: Dict[int, int] = {}
    levels: List[Level] = []

    for azr_level in model.levels:
        if active_only and not any(c.active for c in azr_level.channels):
            continue
        key = (azr_level.J, azr_level.parity)
        if key not in jgroup_of:
            jgroup_of[key] = len(jgroup_of) + 1
        jg = jgroup_of[key]
        counter[jg] = counter.get(jg, 0) + 1
        idx = counter[jg]

        chans = []
        for c in azr_level.channels:
            chans.append(Channel(
                L=c.L, S=c.S, pair=c.pair,
                radiation_type=_radiation_type(c),
                Z1=c.Z1, Z2=c.Z2, M1=c.M1, M2=c.M2,
                channel_radius=c.channel_radius,
                sep_energy=c.sep_energy, excitation=c.excitation,
                pair_J2=c.J2, pair_parity2=c.parity2))

        lv = Level(energy=azr_level.energy, channels=chans,
                   J=azr_level.J, parity=azr_level.parity,
                   jgroup=jg, level=idx)
        if gammas and (jg, idx) in gammas:
            g = gammas[(jg, idx)]
            if len(g) != len(chans):
                raise ValueError(f"level ({jg},{idx}) has {len(chans)} channels, "
                                 f"got {len(g)} gammas")
            for c, v in zip(chans, g):
                c.gamma = float(v)
        levels.append(lv)

    return levels


def levels_from_scheme(scheme, gammas=None) -> List[Level]:
    """Build :class:`Level` objects from a :class:`pyazr.LevelScheme`.

    Uses the running instance's pair table for masses, charges and radii.  As
    with :func:`levels_from_azr`, ``SchemeChannel.width`` is a *physical* width,
    so pass reduced-width amplitudes in ``gammas`` keyed by ``(jgroup, level)``.
    """
    pairs = {p.number: p for p in scheme.pairs}
    out = []
    for lv in scheme.levels:
        chans = []
        for c in lv.channels:
            p = pairs[c.pair]
            chans.append(Channel(
                L=c.L, S=c.S, pair=c.pair,
                radiation_type=c.radiation_type or "P",
                Z1=p.Z1, Z2=p.Z2, M1=p.M1, M2=p.M2,
                channel_radius=p.channel_radius,
                sep_energy=p.sep_energy, excitation=p.excitation,
                pair_J2=p.J2, pair_parity2=p.parity2))
        level = Level(energy=lv.energy, channels=chans, J=lv.J,
                      parity=lv.parity, jgroup=lv.jgroup, level=lv.level)
        if gammas and (lv.jgroup, lv.level) in gammas:
            for c, v in zip(chans, gammas[(lv.jgroup, lv.level)]):
                c.gamma = float(v)
        out.append(level)
    return out
