"""Dimensionless widths: Wigner limits for particle channels, Weisskopf units
for gamma channels.

A fitted partial width in eV says little on its own -- whether it is physically
sensible depends on how it compares with the largest width a single-particle
picture allows.  This module builds that comparison for a whole level scheme at
once, at whatever parameter vector you hand it:

* **particle channels** -- the Wigner (single-particle) limit
  ``gamma^2_W = hbar^2 / (mu a^2)``, and the dimensionless reduced width
  ``theta^2 = Gamma_c / Gamma_W`` with ``Gamma_W = 2 P_l(E_r) gamma^2_W``.
  ``theta^2 <= 1`` is the sum-rule expectation.
* **photon channels** -- the Weisskopf single-particle estimate for the
  multipolarity, and the strength in Weisskopf units, ``W.u. = Gamma_gamma /
  Gamma_W(Weisskopf)``.  A few W.u. is normal; hundreds is not.

Two Wigner limits, and why they differ
--------------------------------------
AZURE2's API reports ``Parameter.wigner_limit`` = ``gamma^2_W`` in **MeV**:
energy-independent, already squared.  The AZURE2 **GUI** shows a different
number next to each channel, ``Gamma_W = 2 P_l(E_r) gamma^2_W`` in **eV**,
which is energy-dependent.  They give different dimensionless widths:

    theta^2_formal = gamma^2_formal / gamma^2_W        (from the API limit)
    theta^2        = Gamma_c        / Gamma_W          (from the GUI limit)

and the two differ by the level-shift derivative, because AZURE2 reports
*observed* partial widths ``Gamma_c = 2 P_c gamma_c^2 / (1 + sum_c' gamma_c'^2
dS_c'/dE)``.  For narrow levels the denominator is ~1 and the two agree; for
broad ones they can differ by tens of percent.  Prefer ``theta^2`` -- it is the
number the GUI shows, one channel at a time.  Never divide by ``wigner_limit``
twice: it is already a squared quantity.

Getting the penetrability without a Coulomb-function implementation
-------------------------------------------------------------------
``Gamma_W`` needs ``P_l(E_r)``, which the API does not expose.  Rather than
re-implement Coulomb wave functions, this module asks AZURE2 itself: it
transforms a *probe* vector in which every reduced-width amplitude is set to a
tiny ``eps`` (so the level-shift denominator is 1 to O(eps^2), and the levels of
a J-group no longer mix), and reads

    2 P_c = Gamma_c(eps) / eps^2

off the result.  One extra ``TRANSFORM_RWA`` call per model, at AZURE2's own
channel radii, boundary conditions and Coulomb functions.

Closed channels
---------------
Below its threshold a channel has ``P_l = 0`` and AZURE2 reports an ANC
(fm^-1/2) instead of a width, so ``Gamma_W`` and ``theta^2`` do not exist there.
Those rows carry ``theta2_formal`` only, and ``is_open`` is False.
"""

from dataclasses import dataclass, fields
from typing import Optional


# Weisskopf single-particle estimates, Gamma in eV with E_gamma in MeV.
# Bohr & Mottelson / Krane; A is the mass number of the emitting nucleus.
_WEISSKOPF = {
    ("E", 1): lambda A, e: 6.8e-2 * A ** (2 / 3) * e ** 3,
    ("E", 2): lambda A, e: 4.9e-8 * A ** (4 / 3) * e ** 5,
    ("E", 3): lambda A, e: 2.3e-14 * A ** 2 * e ** 7,
    ("E", 4): lambda A, e: 6.8e-21 * A ** (8 / 3) * e ** 9,
    ("M", 1): lambda A, e: 2.1e-2 * e ** 3,
    ("M", 2): lambda A, e: 1.5e-8 * A ** (2 / 3) * e ** 5,
    ("M", 3): lambda A, e: 6.8e-15 * A ** (4 / 3) * e ** 7,
    ("M", 4): lambda A, e: 2.0e-21 * A ** 2 * e ** 9,
}


def weisskopf_width(radiation_type, L, e_gamma, mass_number):
    """Weisskopf single-particle width in **eV**.

    ``radiation_type`` is ``'E'`` or ``'M'``, ``L`` the multipolarity, and
    ``e_gamma`` the transition energy in MeV.  Returns ``None`` for a
    multipolarity outside the tabulated E1-E4 / M1-M4 range or for a
    non-positive transition energy.
    """
    f = _WEISSKOPF.get((radiation_type, int(L)))
    if f is None or e_gamma is None or e_gamma <= 0:
        return None
    return f(float(mass_number), float(e_gamma))


def teichmann_wigner(wigner_gamma2):
    """The Teichmann-Wigner sum-rule unit ``3 hbar^2 / (2 mu a^2)``.

    Some references quote dimensionless widths against this instead of
    ``hbar^2/(mu a^2)``; it is 1.5x larger, so their theta^2 is 1.5x smaller.
    """
    return None if wigner_gamma2 is None else 1.5 * wigner_gamma2


@dataclass
class ChannelWidth:
    """One channel of one level, with its width made dimensionless.

    Particle channels carry :attr:`wigner_width` / :attr:`theta2`; photon
    channels carry :attr:`e_gamma`, :attr:`weisskopf` and :attr:`wu`.  Both
    carry :attr:`theta2_formal`, which needs no penetrability.
    """

    # identity
    name: str                       # AZURE2 parameter name, e.g. 'width_3_6'
    index: int                      # position among all parameters
    free_index: Optional[int]       # position in the free-parameter vector
    fixed: bool
    jgroup: int
    level: int
    jpi: str
    level_energy: float             # excitation energy of the level (MeV)
    pair: int
    L: int
    S: Optional[float]
    radiation_type: str             # 'P' for particle, 'E'/'M' for photon

    # values
    gamma: Optional[float]          # formal reduced-width amplitude (MeV^1/2)
    value: Optional[float]          # AZURE2's physical value: Gamma in eV for an
                                    # open channel, ANC in fm^-1/2 for a closed one

    # particle channels
    is_photon: bool = False
    is_open: bool = True            # above threshold at this level energy?
    threshold: Optional[float] = None       # channel threshold, in Ex (MeV)
    wigner_gamma2: Optional[float] = None   # hbar^2/(mu a^2), MeV
    wigner_width: Optional[float] = None    # 2 P_l(E_r) gamma^2_W, eV
    theta2: Optional[float] = None          # Gamma / Gamma_W
    theta2_formal: Optional[float] = None   # gamma^2 / gamma^2_W

    # photon channels
    e_gamma: Optional[float] = None         # transition energy, MeV
    weisskopf: Optional[float] = None       # single-particle estimate, eV
    wu: Optional[float] = None              # Gamma_gamma / weisskopf

    @property
    def channel(self) -> str:
        """Compact channel label, e.g. ``'pair1 L=3 S=1/2'`` or ``'E1->pair4'``."""
        if self.is_photon:
            return f"{self.radiation_type}{self.L}->pair{self.pair}"
        s = "" if self.S is None else f" S={self.S:g}"
        return f"pair{self.pair} L={self.L}{s}"

    @property
    def dimensionless(self) -> Optional[float]:
        """The headline number: W.u. for a photon, theta^2 for a particle."""
        return self.wu if self.is_photon else self.theta2


class WidthTable(list):
    """A list of :class:`ChannelWidth`, one per channel of the level scheme."""

    @property
    def photons(self):
        """Only the photon channels."""
        return WidthTable(c for c in self if c.is_photon)

    @property
    def particles(self):
        """Only the particle channels."""
        return WidthTable(c for c in self if not c.is_photon)

    @property
    def free(self):
        """Only the channels the fit varies."""
        return WidthTable(c for c in self if not c.fixed)

    @property
    def nonzero(self):
        """Only the channels with a nonzero width."""
        return WidthTable(c for c in self if c.value)

    def by_level(self):
        """``{(jgroup, level): WidthTable}`` -- the channels of each level."""
        out = {}
        for c in self:
            out.setdefault((c.jgroup, c.level), WidthTable()).append(c)
        return out

    def to_records(self):
        """The table as a list of plain dicts."""
        keys = [f.name for f in fields(ChannelWidth)]
        return [{k: getattr(c, k) for k in keys} for c in self]

    def table(self, nonzero_only=False):
        """Aligned text table; ``nonzero_only`` hides the inert channels."""
        src = self.nonzero if nonzero_only else self
        rows = [("J^pi", "Ex(MeV)", "channel", "free", "gamma", "value",
                 "unit", "Gamma_W/W.u.est", "theta2 | W.u.", "theta2_formal")]
        for c in src:
            if c.is_photon:
                unit, limit = "eV", c.weisskopf
                dim = c.wu
            elif c.is_open:
                unit, limit = "eV", c.wigner_width
                dim = c.theta2
            else:
                unit, limit = "fm^-1/2", None
                dim = None
            rows.append((
                c.jpi,
                "" if c.level_energy is None else f"{c.level_energy:.4f}",
                c.channel,
                "" if c.fixed else "*",
                "" if c.gamma is None else f"{c.gamma:+.4g}",
                "" if c.value is None else f"{c.value:.6g}",
                unit,
                "" if limit is None else f"{limit:.5g}",
                "" if dim is None else f"{dim:.4g}",
                "" if c.theta2_formal is None else f"{c.theta2_formal:.4g}",
            ))
        w = [max(len(r[i]) for r in rows) for i in range(len(rows[0]))]
        return "\n".join("  ".join(cell.ljust(w[i]) for i, cell in enumerate(r))
                         for r in rows)

    def __repr__(self):
        return (f"WidthTable({len(self)} channels, "
                f"{len(self.particles)} particle, {len(self.photons)} photon)")
