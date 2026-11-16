// ObservablesAdjoint.h against central differences of the forward functions of Observables.h in the real and the
// imaginary part of every T-matrix element.  Six elements with fixed values stand for the pathways of a reaction;
// no project file, nothing of the engine.
//  - the angle-integrated cross section of three (J, l, l') groups with a symmetry factor;
//  - the total cross section by the optical theorem; the beta-delayed sum;
//  - the phase shift of U summed over two diagonal pathways;
//  - the differential cross section sum_L B_L P_L from interference terms: by pairs (interference_pair_cotangent)
//    and through the coefficients q_L (add_interference_cotangents, legendre_contraction); a Legendre coefficient
//    itself, absolute and normalised (the normalised one with the cotangent of sigma_int);
//  - an unobserved primary with an observed secondary (upos_term);
//  - the Coulomb-nuclear interference of elastic scattering, with and without the spin-dependent exchange;
//  - the analysing power of capture, a ratio of two sums.
// Differences: Richardson's (4 D(h/2) - D(h)) / 3 of the central difference D, h = 1e-3; the agreement reached is
// printed with each check.
//   physics_observables_adjoint_test
#include <cmath>
#include <complex>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "azure2/physics/ObservablesAdjoint.h"

using namespace azure2::physics;
typedef std::complex<double> complex;
typedef std::vector<complex> Ts;

namespace {
int failures = 0;
const int nT = 6;
const double geom = 37.5, i1i2 = 0.25;

Ts start() {
  return {complex(0.31, -0.22), complex(-0.12, 0.45), complex(0.08, 0.17),
          complex(-0.4, -0.09), complex(0.27, 0.33), complex(0.05, -0.51)};
}

// The largest |Tbar - difference| over Re and Im of every element, in units of the largest difference.
void check(const std::string &what, const std::function<double(const Ts &)> &f, const Ts &bar, double tol = 1.0e-9) {
  const Ts T0 = start();
  const double h = 1.0e-3;
  double scale = 0.0, diff = 0.0;
  for (int i = 0; i < nT; i++)
    for (int part = 0; part < 2; part++) {
      const complex unit = part ? complex(0.0, 1.0) : complex(1.0, 0.0);
      auto D = [&](double step) {
        Ts a = T0, b = T0;
        a[i] += step * unit;
        b[i] -= step * unit;
        return (f(a) - f(b)) / (2.0 * step);
      };
      const double fd = (4.0 * D(0.5 * h) - D(h)) / 3.0;
      const double an = part ? std::imag(bar[i]) : std::real(bar[i]);
      scale = std::max(scale, std::fabs(fd));
      diff = std::max(diff, std::fabs(an - fd));
    }
  const bool ok = scale > 0.0 && diff <= tol * scale;
  std::printf("  %s  %s  (largest difference %.2e of the largest derivative; limit %.0e)\n", ok ? "ok  " : "FAIL",
              what.c_str(), scale > 0.0 ? diff / scale : -1.0, tol);
  if (!ok) failures++;
}

double legendre(int l, double x) {
  if (l == 0) return 1.0;
  double p0 = 1.0, p1 = x;
  for (int k = 2; k <= l; k++) {
    const double p2 = ((2 * k - 1) * x * p1 - (k - 1) * p0) / k;
    p0 = p1;
    p1 = p2;
  }
  return p1;
}
}  // namespace

int main() {
  const double bar = 0.83;  // dL/d(observable)
  const Ts T0 = start();
  const double x = std::cos(0.9);

  std::printf("cross sections\n");
  {
    // three (J, l, l') groups: elements {0, 1}, {2}, {3, 4, 5}
    const double J[nT] = {0.5, 0.5, 1.5, 2.5, 2.5, 2.5};
    const int l[nT] = {0, 0, 1, 2, 2, 2}, lp[nT] = {1, 1, 1, 0, 0, 0};
    const double symmetry = 2.0;
    auto sigma = [&](const Ts &T) {
      PartialWaveSums sums;
      for (int i = 0; i < nT; i++) sums.add(J[i], l[i], lp[i], T[i]);
      complex sum(0.0, 0.0);
      add_partial_wave_cross_section(sums, geom, i1i2, sum);
      return bar * symmetry * cross_section_barns(sum);
    };
    PartialWaveSums sums;
    for (int i = 0; i < nT; i++) sums.add(J[i], l[i], lp[i], T0[i]);
    Ts tb(nT);
    for (int i = 0; i < nT; i++)
      for (const PartialWaveSum &r : sums.rows)
        if (r.J == J[i] && r.l == l[i] && r.lp == lp[i]) tb[i] = partial_wave_cotangent(bar, geom, J[i], i1i2, symmetry, r.T);
    check("angle-integrated: partial_wave_cotangent", sigma, tb);

    auto total = [&](const Ts &T) {
      double sum = 0.0;
      for (int i = 0; i < nT; i++) add_optical_theorem_term(J[i], T[i], sum);
      return bar * total_cross_section(geom, i1i2, sum);
    };
    for (int i = 0; i < nT; i++) tb[i] = optical_theorem_cotangent(bar, geom, i1i2, J[i]);
    check("total (optical theorem): optical_theorem_cotangent", total, tb);

    auto beta = [&](const Ts &T) {
      complex sum(0.0, 0.0);
      for (int i = 0; i < nT; i++) add_beta_delayed_term(T[i], sum);
      return bar * cross_section_barns(sum);
    };
    for (int i = 0; i < nT; i++) tb[i] = beta_delayed_cotangent(bar, T0[i]);
    check("beta-delayed: beta_delayed_cotangent", beta, tb);
  }

  std::printf("phase shift\n");
  {
    const complex e1 = std::polar(1.0, 0.7), e2 = std::polar(1.0, -0.3);
    auto phase = [&](const Ts &T) {
      const complex U = collision_from_t_matrix(e1, T[1]) + collision_from_t_matrix(e2, T[4]);
      return bar * phase_shift_degrees(U, false);
    };
    const complex U = collision_from_t_matrix(e1, T0[1]) + collision_from_t_matrix(e2, T0[4]);
    const complex Ubar = phase_shift_cotangent(bar, U);
    Ts tb(nT, complex(0.0, 0.0));
    tb[1] = collision_cotangent_to_t_matrix(e1, Ubar);
    tb[4] = collision_cotangent_to_t_matrix(e2, Ubar);
    check("phase_shift_cotangent, collision_cotangent_to_t_matrix", phase, tb);
  }

  std::printf("differential cross section\n");
  {
    // interference terms (order L, elements, z); nL orders
    struct Term {
      int L, a, b;
      double z;
    };
    const std::vector<Term> terms = {{0, 0, 0, 1.1}, {0, 1, 1, 0.9},  {1, 0, 1, -0.6}, {1, 1, 0, -0.6}, {2, 2, 3, 0.45},
                                     {2, 3, 2, 0.45}, {2, 4, 4, 0.7}, {3, 4, 5, 0.3},  {3, 5, 4, 0.3},  {4, 5, 5, -0.2}};
    const int nL = 5;
    const double symmetry = 4.0;
    auto coefficients = [&](const Ts &T) {
      std::vector<complex> b(nL, complex(0.0, 0.0));
      for (const Term &t : terms) add_interference_product(t.z, T[t.a], T[t.b], b[t.L]);
      std::vector<double> c(nL);
      for (int L = 0; L < nL; L++) c[L] = std::real(b[L]);
      scale_legendre_coefficients(c, geom, i1i2, symmetry);
      return c;
    };
    auto dsigma = [&](const Ts &T) {
      const std::vector<double> c = coefficients(T);
      double s = 0.0;
      for (int L = 0; L < nL; L++) s += c[L] * legendre(L, x);
      return bar * s;
    };
    const double cR = legendre_sum_factor(geom, i1i2, symmetry);
    Ts pairs(nT, complex(0.0, 0.0));
    for (const Term &t : terms) {
      const double W = interference_weight(cR, t.z, legendre(t.L, x));
      pairs[t.a] += interference_pair_cotangent(bar, W, T0[t.b]);
      pairs[t.b] += interference_pair_cotangent(bar, W, T0[t.a]);
    }
    check("by pairs: interference_weight, interference_pair_cotangent", dsigma, pairs);

    std::vector<complex> q((std::size_t)nT * nL, complex(0.0, 0.0));
    for (const Term &t : terms) add_interference_cotangents(t.z, T0[t.a], T0[t.b], q[t.a * nL + t.L], q[t.b * nL + t.L]);
    std::vector<double> P(nL);
    for (int L = 0; L < nL; L++) P[L] = legendre(L, x);
    Ts byOrder(nT);
    for (int i = 0; i < nT; i++) byOrder[i] = complex(bar * cR) * legendre_contraction(&q[i * nL], P.data(), nL);
    check("by order: add_interference_cotangents, legendre_contraction", dsigma, byOrder);

    // The Legendre coefficient of order 2 itself: absolute (A_2 = B_2 / symmetry) and normalised
    // (a_2 = A_2 4 pi / sigma_int, sigma_int = 4 pi A_0), the latter with the cotangent of sigma_int.
    const int k = 2;
    auto absolute = [&](const Ts &T) { return bar * coefficients(T)[k] / symmetry; };
    auto normalised = [&](const Ts &T) {
      const std::vector<double> c = coefficients(T);
      return bar * c[k] / c[0];
    };
    for (int form = 0; form < 2; form++) {
      const std::vector<double> c0 = coefficients(T0);
      const double sigmaInt = 4.0 * kObservablePi * c0[0] / symmetry;
      const double factor = legendre_coefficient_factor(form == 0, symmetry, sigmaInt);
      Ts tb(nT, complex(0.0, 0.0));
      for (int i = 0; i < nT; i++) tb[i] = complex(bar * cR) * (q[i * nL + k] * factor);
      if (form == 1) {
        // sigma_int = (4 pi / symmetry) c_R Re sum_{L=0} z T1 conj(T2): its cotangent is -bar a_k / sigma_int
        const double sbar = -bar * (c0[k] / c0[0]) / sigmaInt;
        for (int i = 0; i < nT; i++) tb[i] += complex(sbar * 4.0 * kObservablePi / symmetry * cR) * q[i * nL + 0];
      }
      check(form == 0 ? "the coefficient A_2: legendre_coefficient_factor (absolute)"
                      : "the coefficient a_2 = A_2 / A_0: legendre_coefficient_factor (normalised)",
            form == 0 ? std::function<double(const Ts &)>(absolute) : std::function<double(const Ts &)>(normalised), tb);
    }

    // An unobserved primary with an observed secondary: order-2 terms weighted by R_2.
    const double R2 = upos_mixed_rk(0.6, -0.2, 0.35, 0.4);
    auto upos = [&](const Ts &T) {
      std::vector<double> c(nL, 0.0);
      for (const Term &t : terms)
        if (t.L == 2) c[2] += upos_term(t.z, T[t.a], T[t.b], 2, R2);
      scale_legendre_coefficients(c, geom, i1i2, 1.0);
      return bar * c[2] * legendre(2, x);
    };
    Ts tb(nT, complex(0.0, 0.0));
    const double c1 = legendre_sum_factor(geom, i1i2, 1.0);
    for (const Term &t : terms) {
      if (t.L != 2) continue;
      const double W = upos_interference_weight(c1, t.z, 2, R2, legendre(2, x));
      tb[t.a] += interference_pair_cotangent(bar, W, T0[t.b]);
      tb[t.b] += interference_pair_cotangent(bar, W, T0[t.a]);
    }
    check("unobserved primary: upos_interference_weight", upos, tb);
  }

  std::printf("elastic scattering: the Coulomb-nuclear interference\n");
  for (int spinDependent = 0; spinDependent < 2; spinDependent++) {
    // diagonal pathways: channel spin, l, statistical factor
    const double s[nT] = {0.0, 0.0, 1.0, 1.0, 2.0, 2.0};
    const int l[nT] = {0, 2, 1, 3, 0, 2};
    const double stat[nT] = {1.0, 5.0, 9.0, 21.0, 5.0, 25.0};
    const complex fd(0.8, -1.7), fx = spinDependent ? complex(-0.3, 0.6) : complex(0.0, 0.0);
    const bool identical = spinDependent;
    auto value = [&](const Ts &T) {
      CoulombInterference terms;
      for (int i = 0; i < nT; i++) terms.add(s[i], l[i], stat[i], T[i]);
      ElasticInterference e;
      elastic_interference(terms, geom, identical, spinDependent, 1.0, 10, e);
      complex sum(0.0, 0.0);
      for (std::size_t L = 0; L < e.direct.size(); L++) sum += fd * e.direct[L] * legendre((int)L, x);
      for (std::size_t L = 0; L < e.exchange.size(); L++) sum += fx * e.exchange[L] * legendre((int)L, x);
      return bar * std::real(sum);
    };
    const complex cI = coulomb_interference_factor(geom, identical ? 2.0 : 1.0);
    Ts tb(nT);
    for (int i = 0; i < nT; i++)
      tb[i] = coulomb_interference_cotangent(bar, cI, stat[i], channel_spin_coulomb_amplitude(s[i], fd, fx), legendre(l[i], x));
    check(spinDependent ? "identical particles with spin: f_d + (-1)^s f_x" : "coulomb_interference_cotangent", value, tb);
  }

  std::printf("analysing power of capture\n");
  {
    struct Term {
      int k, a, b;
      double ak, bk;
    };
    const std::vector<Term> terms = {{0, 0, 0, 1.3, 0.0}, {0, 1, 1, 0.8, 0.0}, {1, 0, 1, 0.5, 0.9},
                                     {1, 2, 3, -0.4, 0.6}, {2, 4, 5, 0.3, -0.7}, {2, 1, 4, 0.0, 0.35}};
    const double pk[3] = {1.0, x, legendre(2, x)};
    const double sn = std::sqrt(1.0 - x * x);
    const double pk1[3] = {0.0, -sn, -3.0 * x * sn};  // P_k^1 (the sign convention does not enter the check)
    auto ay = [&](const Ts &T) {
      double den[3] = {0.0, 0.0, 0.0}, num[3] = {0.0, 0.0, 0.0};
      for (const Term &t : terms) add_capture_ay_term(t.ak, t.bk, T[t.a], T[t.b], den[t.k], num[t.k]);
      double d = 0.0, n = 0.0;
      for (int k = 0; k < 3; k++) {
        d += den[k] * pk[k];
        n += num[k] * pk1[k];
      }
      return bar * n / d;
    };
    double den = 0.0, num = 0.0;
    for (const Term &t : terms) add_capture_ay_sums(t.ak, t.bk, pk[t.k], pk1[t.k], T0[t.a], T0[t.b], den, num);
    Ts tb(nT, complex(0.0, 0.0));
    for (const Term &t : terms) {
      complex d1, d2;
      capture_ay_cotangents(t.ak, t.bk, pk[t.k], pk1[t.k], num / den, den, T0[t.a], T0[t.b], d1, d2);
      tb[t.a] += bar * d1;
      tb[t.b] += bar * d2;
    }
    check("add_capture_ay_sums, capture_ay_cotangents", ay, tb);
  }

  std::printf("%s\n", failures ? "FAILED" : "all checks passed");
  return failures ? 1 : 0;
}
