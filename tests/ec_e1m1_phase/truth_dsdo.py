#!/usr/bin/env python3
"""Independent differential cross sections for p + 16O external capture with E1 and M1 (writes expected/dsdo_s.dat,
expected/dsdo_p.dat).  Companion of audit 2026-10-25, 04_external_capture_photons S5 (checklist t3-26): the relative
phase of the E and M external-capture amplitudes.

Nothing here uses multipole phase conventions.  The emission amplitude is the matrix element of the dipole
interaction H = -d.E - mu.B for the emitted photon (plane wave A = eps e^{i(k.r - w t)}: E = (i w/c) A and
B = i k x A carry the SAME factor, so H_int is proportional to d.eps* + mu.(khat x eps*)), with d = e e_1 r the
long-wavelength E1 dipole of the relative motion (e_1 = Z1 m2/M - Z2 m1/M) and mu = mu_N (g_l L + g_p S_p) the
magnetic moment of the relative motion and the proton (g_l = m_p mu (Z1/m1^2 + Z2/m2^2)), between the final bound
state u_f(r)/r Y_{lf m'} chi_{m_f} (ANC-normalized Whittaker tail) and the Coulomb + hard-sphere scattering state of
the proton with the beam along z,

  psi(+)_{m_s} = sum_l i^l sqrt(4 pi (2l+1)) e^{i(sigma_l + delta_l)} [u_l(r)/r] Y_l0(rhat) chi_{m_s},
  u_l = sqrt(2 mu/(pi hbar^2 k)) [F_l cos delta_l + G_l sin delta_l],  delta_l = -atan(F_l(ka)/G_l(ka)),

sigma_l = arg Gamma(l + 1 + i eta) the Coulomb phase; the radial functions are those of tests/ec_m1_exact/truth.py.
The photon goes out along khat = (sin th, 0, cos th) with helicity lambda = +-1, eps*_lambda = (cos th, -i lambda,
-sin th)/sqrt2 up to a phase, so khat x eps*_lambda = (i lambda cos th, 1, -i lambda sin th)/sqrt2.  Angular matrix
elements <Y_{lf m'}| rhat_i |Y_{l0}> are quadratures of scipy's spherical harmonics (Condon-Shortley), and L_i are the
matching angular-momentum matrices; the final m' and m_f are summed, so the basis convention cancels.

  A(lambda, m', m_f, m_s) = sum_{l: E1} i^l sqrt(4 pi (2l+1)) e^{i(sigma_l + delta_l)} e e_1 I_E1(l) <lf m'|rhat.eps*|l 0> delta_{m_f m_s}
                          + i^{lf} sqrt(4 pi (2 lf+1)) e^{i(sigma_lf + delta_lf)} mu_N I_M1
                                [g_l <lf m'|L.(khat x eps*)|lf 0> delta_{m_f m_s} + g_p delta_{m' 0} <m_f|S.(khat x eps*)|m_s>],
  I_E1(l) = int u_f r u_l dr,  I_M1 = int u_f u_{lf} dr,  mu_N = e hbar c/(2 m_p c^2) = e x 0.10515 fm,

(M1 connects l = lf only; E1 connects l = lf +- 1) and dsigma/dOmega = N (1/2) sum_{m_s, m', m_f, lambda} |A|^2.
The constant N is the golden-rule factor common to E1 and M1.  It is fixed by the E1-only angle integral, which must
equal the (l s)J reduced-matrix-element cross section of truth.py (sigma = pi/k^2 sum_J g_J 2 pi Gamma_J with
Gamma_J from <f||O_1||i>); the M1-only integral is then CHECKED against the same formula with truth.py's
red_spin_mscheme (orbital + spin).  Both agree to 1e-9, so the operators' relative normalization is right and the
interference carries no free factor.

Case s (17F 1/2+, lf = 0): E1 from l = 1, M1 (spin term only) from l = 0: the interference is identically zero
(E1 spin-diagonal, S traceless), so dsigma/dOmega(E1 + M1) = E1 + M1.
Case p (fictitious 1/2- at E_x 0.4996, lf = 1, ANC 20 fm^-1/2: an ANC of 84 there would put 130 % of the state outside
5 fm): E1 from l = 0, 2, M1 (orbital + spin) from l = 1: the orbital M1 term interferes with E1 in the odd Legendre orders.

expected/dsdo_<case>.dat: E_cm (MeV), theta_cm (deg), dsigma/dOmega E1 only, M1 only, E1 + M1 (barn/sr).
Run from this directory: python3 truth_dsdo.py   (needs mpmath, numpy, scipy, sympy; a few minutes).
"""
import os
import sys

import mpmath as mp
import numpy as np
from scipy import special
from sympy.physics.wigner import wigner_6j

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "ec_m1_exact"))
import truth as T          # noqa: E402  (System, red_Y, orbital, gamma_pref, red_spin_mscheme, constants)

import importlib.util      # noqa: E402

# this directory's make_models (truth.py already imported ec_m1_exact's under the same module name)
_spec = importlib.util.spec_from_file_location("ec_e1m1_models", os.path.join(HERE, "make_models.py"))
mm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mm)

GP = 5.5857
SX = np.array([[0, 1], [1, 0]], dtype=complex) / 2
SY = np.array([[0, -1j], [1j, 0]], dtype=complex) / 2
SZ = np.array([[1, 0], [0, -1]], dtype=complex) / 2


def coulomb_phase(l, eta):
    return float(mp.im(mp.loggamma(l + 1 + 1j * eta)))


def ylm(m, l, phi, theta):
    """Y_lm(theta, phi) with the Condon-Shortley phase (scipy's lpmv carries it), for any scipy version."""
    from math import factorial
    am = abs(m)
    y = np.sqrt((2 * l + 1) / (4 * np.pi) * factorial(l - am) / factorial(l + am)) * special.lpmv(am, l, np.cos(theta)) * np.exp(1j * am * phi)
    return y if m >= 0 else (-1) ** am * np.conj(y)


def angular_r(lf, li):
    """<Y_{lf m'}| rhat_i |Y_{li 0}> for i = x, y, z: array [m' index, i], m' = -lf .. lf, by quadrature."""
    x, w = np.polynomial.legendre.leggauss(48)
    phi = np.linspace(0, 2 * np.pi, 96, endpoint=False)
    th = np.arccos(x)
    TH, PH = np.meshgrid(th, phi, indexing="ij")
    WT = np.outer(w, np.full(phi.size, 2 * np.pi / phi.size))
    rhat = [np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)]
    yi = ylm(0, li, PH, TH)
    out = np.zeros((2 * lf + 1, 3), dtype=complex)
    for k, m in enumerate(range(-lf, lf + 1)):
        yf = ylm(m, lf, PH, TH)
        for i in range(3):
            out[k, i] = np.sum(WT * np.conj(yf) * rhat[i] * yi)
    return out


def l_matrices(l):
    """L_x, L_y, L_z in the |l m> basis, m = -l .. l (Condon-Shortley)."""
    n = 2 * l + 1
    ms = np.arange(-l, l + 1)
    Lp = np.zeros((n, n), dtype=complex)
    for k, m in enumerate(ms[:-1]):
        Lp[k + 1, k] = np.sqrt(l * (l + 1) - m * (m + 1))
    Lm = Lp.conj().T
    Lx = (Lp + Lm) / 2
    Ly = (Lp - Lm) / (2j)
    Lz = np.diag(ms.astype(complex))
    return Lx, Ly, Lz


class Case:
    def __init__(self, final):
        J, pi, ex, lf, anc = mm.FINAL[final]
        self.final, self.lf, self.Jf = final, lf, J
        self.B = mm.SEP - ex
        self.sys = T.System(1, 8, mm.MP, mm.M16, 5.0, self.B, lf, anc)
        self.e1 = self.sys.z[0] * self.sys.c[0] - self.sys.z[1] * self.sys.c[1]
        self.gl = mm.MP * (mm.MP * mm.M16 / (mm.MP + mm.M16)) * (1 / mm.MP ** 2 + 8 / mm.M16 ** 2)
        self.l_e1 = [l for l in (lf - 1, lf + 1) if l >= 0]
        self.ang = {l: angular_r(lf, l) for l in self.l_e1}
        self.L = l_matrices(lf)
        from sympy import S as Sy
        from sympy.physics.wigner import clebsch_gordan
        self.cg = {(mp, mf, Mf): float(clebsch_gordan(Sy(lf), Sy(1) / 2, Sy(int(2 * J)) / 2, Sy(mp), Sy(int(2 * mf)) / 2, Sy(int(2 * Mf)) / 2) if abs(mp + mf - Mf) < 1e-9 else 0.0)
                   for mp in range(-lf, lf + 1) for mf in (0.5, -0.5) for Mf in (J - np.arange(int(round(2 * J)) + 1))}

    def radial(self, E):
        """Radial integrals and phases of one energy."""
        r = {}
        for l in set(self.l_e1 + [self.lf]):
            k, eta, u = self.sys.ui(l, E)
            d = -np.arctan(float(mp.coulombf(l, eta, k * self.sys.a)) / float(mp.coulombg(l, eta, k * self.sys.a)))
            r[l] = (np.exp(1j * (coulomb_phase(l, eta) + d)), np.sum(self.sys.w * self.sys.uf * self.sys.r * u),
                    np.sum(self.sys.w * self.sys.uf * u), k)
        return r

    def amplitudes(self, rad, th_deg, with_e1=True, with_m1=True):
        """A[lambda index, M_f index, m_s] to the final level |(lf 1/2) Jf Mf> (arbitrary common constant): the
        product-basis amplitudes to |lf m'>|1/2 m_f> contracted with <lf m' 1/2 m_f | Jf Mf>.  (Summing the product
        basis instead would include the other Jf of that orbital, which is not a level of the model.)"""
        Ap = self._product_amplitudes(rad, th_deg, with_e1, with_m1)
        nM = int(round(2 * self.Jf)) + 1
        A = np.zeros((2, nM, 2), dtype=complex)
        for il in range(2):
            for iM, Mf in enumerate(self.Jf - np.arange(nM)):
                for k, mp in enumerate(range(-self.lf, self.lf + 1)):
                    for imf, mf in enumerate((0.5, -0.5)):
                        c = self.cg[(mp, mf, Mf)]
                        if c != 0.0:
                            A[il, iM] += c * Ap[il, k, imf]
        return A

    def _product_amplitudes(self, rad, th_deg, with_e1=True, with_m1=True):
        """A[lambda index, m' index, m_f, m_s] to the product states |lf m'>|1/2 m_f> (arbitrary common constant)."""
        e = np.sqrt(T.ALPHA * T.HBARC)
        muN = e * T.HBARC / (2 * T.MPROTON)
        th = np.radians(th_deg)
        nf = 2 * self.lf + 1
        A = np.zeros((2, nf, 2, 2), dtype=complex)
        for il, lam in enumerate((1, -1)):
            eps = np.array([np.cos(th), -1j * lam, -np.sin(th)]) / np.sqrt(2)                 # eps*
            v = np.array([1j * lam * np.cos(th), 1.0, -1j * lam * np.sin(th)]) / np.sqrt(2)   # khat x eps*
            if with_e1:
                for l in self.l_e1:
                    phase, I_E1, _, _ = rad[l]
                    pref = 1j ** l * np.sqrt(4 * np.pi * (2 * l + 1)) * phase * e * self.e1 * I_E1
                    for k in range(nf):
                        A[il, k] += pref * np.dot(self.ang[l][k], eps) * np.eye(2)
            if with_m1:
                phase, _, I_M1, _ = rad[self.lf]
                pref = 1j ** self.lf * np.sqrt(4 * np.pi * (2 * self.lf + 1)) * phase * muN * I_M1
                Lv = self.L[0] * v[0] + self.L[1] * v[1] + self.L[2] * v[2]
                Sv = SX * v[0] + SY * v[1] + SZ * v[2]
                m0 = self.lf                                       # index of m' = 0 (the initial m_l)
                for k in range(nf):
                    A[il, k] += pref * (self.gl * Lv[k, m0] * np.eye(2) + (GP * Sv if k == m0 else 0.0))
        return A

    def unnormalized(self, rad, th, with_e1=True, with_m1=True):
        return 0.5 * float(np.sum(np.abs(self.amplitudes(rad, th, with_e1, with_m1)) ** 2))

    def angle_integral(self, rad, with_e1, with_m1, n=48):
        x, w = np.polynomial.legendre.leggauss(n)
        return 2 * np.pi * sum(wi * self.unnormalized(rad, np.degrees(np.arccos(xi)), with_e1, with_m1) for xi, wi in zip(x, w))

    def sigma_reduced(self, rad):
        """Angle-integrated E1 and M1 cross sections (barn) from (l s)J reduced matrix elements, as truth.py does."""
        s = 0.5
        e1 = 0.0
        k = rad[self.lf][3]
        E = (k * T.HBARC) ** 2 / (2 * self.sys.mu)
        kg = (E + self.B) / T.HBARC
        for l in self.l_e1:
            for J in sorted({abs(l - s), l + s}):
                if abs(J - self.Jf) > 1 or J < 0:
                    continue
                red = T.orbital(self.lf, l, s, self.Jf, J, 1) * T.red_Y(self.lf, 1, l) * self.e1 * rad[l][1]
                e1 += (2 * J + 1) / 2 * T.gamma_pref(1, kg) * T.ALPHA * T.HBARC * red ** 2 / (2 * J + 1)
        e1 *= np.pi / k ** 2 * 0.01
        muN2 = T.ALPHA * T.HBARC * (T.HBARC / (2 * T.MPROTON)) ** 2
        m1 = 0.0
        for J in sorted({abs(self.lf - s), self.lf + s}):
            if abs(J - self.Jf) > 1:
                continue
            red = T.red_spin_mscheme(0.5, 0.0, self.lf, s, s, J, self.Jf, GP, 0.0, gl=self.gl) * rad[self.lf][2]
            Bm1 = 3 / (4 * np.pi) * muN2 * red ** 2 / (2 * J + 1)
            m1 += np.pi / k ** 2 * (2 * J + 1) / 2 * T.gamma_pref(1, kg) * Bm1 * 0.01
        return e1, m1


def main(cases=("s", "p")):
    os.makedirs(os.path.join(HERE, "expected"), exist_ok=True)
    for final in cases:
        case = Case(final)
        rows = []
        for E in mm.ENERGIES:
            rad = case.radial(E)
            sig_e1, sig_m1 = case.sigma_reduced(rad)
            int_e1 = case.angle_integral(rad, True, False)
            int_m1 = case.angle_integral(rad, False, True)
            N = sig_e1 / int_e1
            ratio = N * int_m1 / sig_m1
            print(f"case {final}, E = {E}: sigma_E1 {sig_e1:.6e} sigma_M1 {sig_m1:.6e} b; m-scheme M1 / reduced M1 = {ratio:.12f}", flush=True)
            assert abs(ratio - 1.0) < 1e-9, ratio
            int_both = case.angle_integral(rad, True, True)
            assert abs(N * int_both / (sig_e1 + sig_m1) - 1.0) < 1e-9, "the interference does not integrate to zero"
            for th in mm.ANGLES:
                rows.append((E, th, N * case.unnormalized(rad, th, True, False), N * case.unnormalized(rad, th, False, True),
                             N * case.unnormalized(rad, th, True, True)))
                print(f"  theta = {th:5.1f}: E1 {rows[-1][2]:.6e}  M1 {rows[-1][3]:.6e}  E1+M1 {rows[-1][4]:.6e}"
                      f"  interference {rows[-1][4] - rows[-1][2] - rows[-1][3]:+.3e}", flush=True)
        with open(os.path.join(HERE, "expected", f"dsdo_{final}.dat"), "w") as f:
            f.write(f"# case {final}: E_cm(MeV)  theta_cm(deg)  dsigma/dOmega E1 only   M1 only   E1 + M1   (barn/sr)"
                    "  -- truth_dsdo.py, independent of AZURE2\n")
            for E, th, a, b, c in rows:
                f.write(f"{E:.4f}  {th:6.1f}  {a:.10e}  {b:.10e}  {c:.10e}\n")


if __name__ == "__main__":
    main(tuple(sys.argv[1:]) or ("s", "p"))
