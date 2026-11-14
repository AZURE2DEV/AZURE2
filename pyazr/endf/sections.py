"""The rows of the MF=3 and MF=6 sections the writer fills from the engine's values (ENDF-102, 2023 edition, cited by
printed page).  Format work only: records, the masses of the products in the manual's units, the normalisation
LTP=12 asks for.

- MF=3 (pp. 103-104): HEAD ``ZA, AWR``; TAB1 ``QM, QI, 0, LR`` of sigma(E), lab eV and b, INT=2.  For
  charged-particle elastic scattering MT=2 holds sigma_NI(E) of eq. 6.19 when MF=6 uses LTP=12 (p. 146; p. 150:
  "With LPT=12 and LPT=14 representation, MT=2 in File 3 contains sigma_NI(E) in barns").
- MF=6 (pp. 134-136): HEAD ``ZA, AWR, JP=0, LCT=2, NK``; per product a TAB1 ``ZAP, AWP, LIP, LAW`` of its yield and
  the law's records.  LCT=2: LAW=2 is in the c.m. (p. 143), LAW=5's mu is the c.m. cosine (p. 144).
- LAW=2 (p. 143): TAB2 ``0, 0, 0, 0, NR, NE``; per incident energy a LIST ``0, E, LANG=0, 0, NW=NL, NL`` of
  a_1 ... a_NL (a_0 = 1 is implied by eq. 6.8).
- LAW=4 (p. 144): no records; the recoil of a two-body reaction and of radiative capture.
- LAW=5 (pp. 146-147): TAB2 ``SPI, 0, LIDP, 0, NR, NE``; per incident energy a LIST ``0, E, LTP=12, 0, NW=2 NL,
  NL`` of (mu, P_NI) pairs, each energy on its own cosines.  P_NI = (sigma_e - sigma_c) / sigma_NI with sigma_NI
  the integral over mu of sigma_e - sigma_c (eqs. 6.19-6.20): the trapezoid rule on the tabulated cosines, exact
  for P_NI linear in mu (LTP=12), so that P_NI integrates to 1 over the table.  Eq. 6.19 integrates over mu, not
  over the solid angle, though p. 146 calls sigma_NI "in barns": the equations are taken as written, so that
  sigma_e = sigma_c + sigma_NI P_NI (eq. 6.20), as the module's reader computes it.
- MF=6/MT=102 (p. 136): per final state, in order of excitation, the primary photon (ZAP=0, LAW=2) and its recoil
  (LAW=4), each with the yield y_f = sigma_f / sigma of MF=3/MT=102 (eq. 12.1, p. 200, the multiplicity of a discrete
  photon).  AWP of the photon "can contain the energy of a primary photon" (p. 135); the manual does not say at which
  incident energy: taken as EG of a primary photon in MF=12 (LP=2, p. 201), the photon energy at zero incident
  energy, the binding energy minus the final level's excitation.  The cascade below a final state (LAW=1 delta
  functions, p. 136) is not AZURE2's: it is not written.

Masses (p. 24, section 0.5.1.1): AWR of the material and of heavy products from atomic masses, the projectile and the
light products (A <= 4) from nuclear masses, all in neutron masses (p. 394).  Taken from AME2020
(:mod:`pyazr.sources.masses`), the nuclear mass as the atomic mass minus Z electron masses (electron binding below
1e-7 u for Z <= 2).
"""

import numpy as np

from .records import RecordWriter
from .material import M_NEUTRON_U

ELECTRON_MASS_U = 0.000548579909  # CODATA 2018
#: A product of mass number up to this is a light product (nuclear mass), above it a heavy one (atomic), p. 24.
LIGHT_PRODUCT_A = 4


def product_mass(z, a):
    """AWP of a product (neutron masses, p. 24): nuclear for A <= 4, atomic above; AME2020."""
    from ..sources import masses
    atomic = masses.atomic_mass((int(z), int(a)))
    return (atomic - (int(z) * ELECTRON_MASS_U if int(a) <= LIGHT_PRODUCT_A else 0.0)) / M_NEUTRON_U


def _yield(w, zap, awp, lip, law, energies_eV, values):
    """The TAB1 of a product (p. 136): ZAP, AWP, LIP, LAW and y(E), INT=2."""
    w.tab1(float(zap), float(awp), int(lip), int(law), [len(energies_eV)], [2], list(energies_eV), list(values))


def cross_section(za, awr, qm, qi, energies_eV, sigma, number_form):
    """MF=3 (pp. 103-104): HEAD and the TAB1 of sigma(E), INT=2."""
    w = RecordWriter(number_form)
    w.cont(za, awr, 0, 0, 0, 0)
    w.tab1(qm, qi, 0, 0, [len(energies_eV)], [2], list(energies_eV), list(sigma))
    return w.rows


def two_body(za, awr, emitted, residual, energies_eV, coefficients, number_form):
    """MF=6 of a two-body exit (p. 143): the emitted particle (LAW=2, yield 1) and the residual (LAW=4, yield 1).
    ``emitted`` and ``residual``: (ZAP, AWP); ``coefficients``: (n_E, L+1) a_l with a_0 = 1."""
    w = RecordWriter(number_form)
    w.cont(za, awr, 0, 2, 2, 0)
    ends = [energies_eV[0], energies_eV[-1]]
    _yield(w, emitted[0], emitted[1], 0, 2, ends, [1.0, 1.0])
    _law2(w, energies_eV, coefficients)
    _yield(w, residual[0], residual[1], 0, 4, ends, [1.0, 1.0])
    return w.rows


def _law2(w, energies_eV, coefficients):
    """The records of LAW=2 with LANG=0 (p. 143): NL the highest order with a coefficient other than 0."""
    w.tab2(0.0, 0.0, 0, 0, [len(energies_eV)], [2], len(energies_eV))
    for e, a in zip(energies_eV, coefficients):
        a = np.asarray(a, float)[1:]
        nonzero = np.nonzero(a)[0]
        nl = int(nonzero[-1]) + 1 if nonzero.size else 0
        w.list(0.0, float(e), 0, 0, nl, list(a[:nl]))


def capture(za, awr, compound, final_states, energies_eV, number_form):
    """MF=6/MT=102 (p. 136): per final state (``(photon energy EG in eV, yields, coefficients)``, in order of
    excitation) the primary photon (ZAP=0, AWP = EG, LAW=2) and the recoil (``compound``: (ZAP, AWP); LAW=4)."""
    w = RecordWriter(number_form)
    w.cont(za, awr, 0, 2, 2 * len(final_states), 0)
    for eg, yields, coefficients in final_states:
        _yield(w, 0, eg, 0, 2, energies_eV, yields)
        _law2(w, energies_eV, coefficients)
        _yield(w, compound[0], compound[1], 0, 4, energies_eV, yields)
    return w.rows


def elastic(za, awr, projectile, target, identical, spin, energies_eV, panels, number_form):
    """MF=3/MT=2 and MF=6/MT=2 of charged-particle elastic scattering, LTP=12 (pp. 144-147): ``panels`` per energy
    (cosines, sigma_e - sigma_c in b/sr at them).  Returns (MF=3 rows, MF=6 rows, sigma_NI per energy)."""
    trapezoid = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    sigma_ni = np.array([trapezoid(v, m) for m, v in panels])
    mf3 = cross_section(za, awr, 0.0, 0.0, energies_eV, sigma_ni, number_form)
    w = RecordWriter(number_form)
    w.cont(za, awr, 0, 2, 1 if identical else 2, 0)
    ends = [energies_eV[0], energies_eV[-1]]
    _yield(w, projectile[0], projectile[1], 0, 5, ends, [1.0, 1.0])
    w.tab2(float(spin), 0.0, 1 if identical else 0, 0, [len(energies_eV)], [2], len(energies_eV))
    for e, (mu, v), s in zip(energies_eV, panels, sigma_ni):
        p = np.asarray(v, float) / s if s != 0.0 else np.zeros(len(v))
        w.list(0.0, float(e), 12, 0, len(mu), [x for pair in zip(mu, p) for x in pair])
    if not identical:  # the recoil of distinguishable particles; identical ones are one product (LIDP=1)
        _yield(w, target[0], target[1], 0, 4, ends, [1.0, 1.0])
    return mf3, w.rows, sigma_ni
