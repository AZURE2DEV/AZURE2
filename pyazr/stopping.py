"""Stopping and straggling of ions in matter: AZURE2's ``stopping/`` library (docs/upgrade/stopping.md).

SRIM-2013's stopping powers as implemented by LunaSRIM (J. Skowronski, github.com/skowrons94/LunaSRIM), ported to
C++; this module calls the C++ functions through ``_azure2.stopping`` element-wise (no second implementation).  Every
result that uses these stopping powers must cite SRIM: :data:`CITATION`.

Units (stopping/'s ``Material.h``): energies in MeV, kinetic, laboratory frame, of the projectile; stopping cross
sections in eV/(10^15 atoms/cm^2) **per atom of the material** (all atoms counted); areal densities in atoms/cm^2 of
the material; straggling widths (1 sigma) in MeV.  An energy loss is Delta E [MeV] = eps N 1e-21.  A stopping per atom
of one element of a compound is eps / (that element's atom fraction), as the engine's ``stopping evaluation`` does for
the target's active nucleus.

Example::

    >>> import numpy as np
    >>> from pyazr import stopping
    >>> ta2o5 = stopping.Material("Ta2O5")                      # solid, compound correction 1
    >>> s = stopping.stopping(ta2o5, "p", np.array([0.2, 0.5, 1.0]))
    >>> s["total"]                                              # eV/(1e15 atoms/cm2) per atom of Ta2O5
    >>> stopping.energy_after(ta2o5, "p", [1.0], [1e18])        # MeV after 1e18 atoms/cm2
    >>> stopping.straggling(ta2o5, "p", [1.0], [1e18], beam_sigma=[0.0005])["total"]   # MeV
    >>> print(stopping.provenance(ta2o5, "p"))                  # SRIM-2013, the citation, per element
"""
import numpy as np

from . import _azure2

__all__ = ["CITATION", "IMPLEMENTATION", "Material", "projectile", "stopping", "energy_loss", "energy_after",
           "areal_density_between", "range_of", "straggling", "provenance", "table", "evaluate"]

_st = _azure2.stopping

#: The citation every use of these stopping powers owes (SRIM-2013).
CITATION = _st.citation
#: What the code is: SRIM-2013's stopping, ported to C++ from LunaSRIM.
IMPLEMENTATION = _st.implementation


class Material:
    """A target material by its chemical formula ("Ta2O5", "TiN", "C", "CH4"; counts may be decimal).

    Parameters
    ----------
    formula : str
    phase : "solid", "gas" or None.  None is SRIM's default: gas for a material of one element that is a gas in its
        natural state, solid for everything else (a gaseous compound such as CH4 or CO2 must be said).
    compound_correction : SRIM's factor C of Bragg's rule (1 = Bragg's rule), applied once.
    """

    def __init__(self, formula, phase=None, compound_correction=1.0):
        self.formula = formula
        self.phase = phase or ""
        self.compound_correction = float(compound_correction)
        info = _st.material(formula, self.phase, self.compound_correction)
        #: [(Z, atom fraction)], fractions summing to 1
        self.constituents = [tuple(c) for c in info["constituents"]]
        #: the phase used ("solid" or "gas")
        self.phase_used = info["phase"]
        #: atom-fraction-weighted atomic weight (u, SRIM's isotope-averaged weights)
        self.mean_atomic_weight_u = info["mean_atomic_weight_u"]

    def atom_fraction(self, z):
        """The atom fraction of element ``z`` (0 if absent)."""
        return sum(f for zz, f in self.constituents if zz == z)

    def _args(self):
        return dict(phase=self.phase, compound_correction=self.compound_correction)

    def __repr__(self):
        return f"Material({self.formula!r}, phase={self.phase_used!r}, compound_correction={self.compound_correction})"


def projectile(name_or_z, mass_u=None):
    """(Z1, mass in u) of a projectile: "p", "d", "t", "3He", "a" (AME2020 atomic masses), or Z1 with its mass."""
    if isinstance(name_or_z, str):
        return tuple(_st.light_projectile(name_or_z))
    if isinstance(name_or_z, tuple):
        return name_or_z
    if mass_u is None:
        raise ValueError("projectile: give the mass (u) with the charge number")
    return int(name_or_z), float(mass_u)


def _material(m):
    return m if isinstance(m, Material) else Material(m)


def _pair(m, p):
    m = _material(m)
    z1, mass = projectile(p) if not isinstance(p, tuple) else p
    return m, z1, mass


def _arrays(*xs):
    b = np.broadcast_arrays(*[np.asarray(x, dtype=float) for x in xs])
    return [np.ascontiguousarray(x).ravel() for x in b], b[0].shape


def stopping(material, proj, energy):
    """Stopping cross section of ``material`` for ``proj`` at lab energies ``energy`` (MeV).

    Returns a dict of arrays: ``electronic``, ``nuclear``, ``total`` (eV/(1e15 atoms/cm^2) per atom of the material)
    and ``rel_uncertainty`` (relative 1 sigma: SRIM's published accuracy for the ion class, stopping.md section 7).
    """
    m, z1, mass = _pair(material, proj)
    (e,), shape = _arrays(energy)
    out = _st.stopping(m.formula, z1, mass, e, **m._args())
    return {k: np.asarray(v).reshape(shape) for k, v in out.items()}


def areal_density_between(material, proj, energy_in, energy_out):
    """N(E_in, E_out) = int_{E_out}^{E_in} dE/eps, atoms/cm^2 of the material (energies in MeV, lab)."""
    m, z1, mass = _pair(material, proj)
    (a, b), shape = _arrays(energy_in, energy_out)
    return np.asarray(_st.areal_density_between(m.formula, z1, mass, a, b, **m._args())).reshape(shape)


def energy_after(material, proj, energy_in, areal_density):
    """E_out (MeV) after an areal density (atoms/cm^2 of the material) from E_in (MeV); refused if N >= the range."""
    m, z1, mass = _pair(material, proj)
    (a, n), shape = _arrays(energy_in, areal_density)
    return np.asarray(_st.energy_after(m.formula, z1, mass, a, n, **m._args())).reshape(shape)


def energy_loss(material, proj, energy_in, areal_density):
    """E_in - E_out (MeV) of a layer of ``areal_density`` atoms/cm^2."""
    return np.asarray(energy_in, dtype=float) - energy_after(material, proj, energy_in, areal_density)


def range_of(material, proj, energy):
    """The range R(E) = N(E, 0), atoms/cm^2 of the material."""
    m, z1, mass = _pair(material, proj)
    (e,), shape = _arrays(energy)
    return np.asarray(_st.range(m.formula, z1, mass, e, **m._args())).reshape(shape)


def straggling(material, proj, energy_in, areal_density, beam_sigma=0.0, model="trim"):
    """Energy after a layer and its spread (MeV): TRIM's ("trim", default) or Bohr's ("bohr") straggling with
    Tschalaer's propagation, the beam's entrance spread ``beam_sigma`` (MeV) carried through the layer.

    Returns a dict of arrays: ``energy_out``, ``straggling`` (the layer's alone), ``beam`` (the entrance spread at the
    exit), ``total`` (both in quadrature).
    """
    m, z1, mass = _pair(material, proj)
    (a, n, s0), shape = _arrays(energy_in, areal_density, beam_sigma)
    out = _st.straggling(m.formula, z1, mass, a, n, s0, model, **m._args())
    return {k: np.asarray(v).reshape(shape) for k, v in out.items()}


def provenance(material, proj):
    """Where the stopping of ``material`` for ``proj`` comes from: SRIM-2013 and its citation, LunaSRIM, the data
    files, phase and compound correction, per element the data row, I, M, phase factor, and the uncertainties."""
    m, z1, mass = _pair(material, proj)
    return _st.provenance(m.formula, z1, mass, **m._args())


def table():
    """The rows of the binding's stopping table as plain data, by name (``pyazr.kinds()["stopping"]`` is the same):
    ``summary``, ``doc``, ``library``, ``axis``, ``group``, ``citation`` (SRIM's: every use cites it),
    ``arguments`` (``formula``, a text, then the numeric ones with unit, kind, domain, default, example, sweep) and
    ``returns`` (name, unit, doc).  The rows are the functions above on arrays with a status per element."""
    return {row["name"]: row for row in _st.functions()}


def evaluate(name, formula, **arguments):
    """Row ``name`` of :func:`table` for the material ``formula`` on numeric arguments that broadcast (numpy rules;
    a choice by its name, or its index): ``(result, status)``, ``result`` a dict of arrays of the broadcast shape (NaN
    where the status is not 0), ``status`` 0 computed, k > 0 numeric argument k (in the row's order after
    ``formula``) outside its domain, -2 the arguments together outside it (a layer at least the range), -1 failed.
    The same C++ functions as :func:`stopping`, :func:`range_of`, :func:`energy_after`,
    :func:`areal_density_between` and :func:`straggling`, bit for bit.  Raises ``KeyError`` for an unknown row or
    argument, ``ValueError`` for a missing argument or a formula the library refuses (its message)."""
    rows = table()
    if name not in rows:
        raise KeyError(f"stopping: no function {name!r}; pyazr.stopping.table() lists {sorted(rows)}")
    numeric = [a for a in rows[name]["arguments"] if a["kind"] != "text"]
    unknown = sorted(set(arguments) - {a["name"] for a in numeric})
    if unknown:
        raise KeyError(f"stopping {name}: unknown argument(s) {unknown}; allowed: {[a['name'] for a in numeric]}")
    values = []
    for a in numeric:
        if a["name"] not in arguments:
            if a["default"] is None:
                raise ValueError(f"stopping {name}: argument {a['name']} ({a['unit'] or 'no unit'}, {a['domain']}) "
                                 f"is missing")
            value = a["default"]
        else:
            value = arguments[a["name"]]
        if a["kind"] == "choice" and isinstance(value, str):
            if value not in a["choices"]:
                raise ValueError(f"stopping {name}: {a['name']} = {value!r}; allowed: {a['domain']}")
            value = a["choices"].index(value)
        values.append(np.asarray(value, dtype=float))
    (columns), shape = _arrays(*values)
    out, status = _st.evaluate(name, formula, columns)
    result = {r["name"]: np.asarray(out[j]).reshape(shape) for j, r in enumerate(rows[name]["returns"])}
    return result, np.asarray(status).reshape(shape)
