// Observables.h against closed forms written out here; no project file, no engine class:
//  - one level, two neutral channels through collision_matrix (CollisionMatrix.h): the reaction cross section is the
//    one-level Breit-Wigner form, and the total by the optical theorem is the elastic plus the reaction one;
//  - pathways of one (J, l, l') add before the modulus, others after;
//  - the phase shift is half the argument of U, wrapped to [0, 180) for identical particles;
//  - an s-wave pathway alone: 4 pi B_0 is the angle-integrated cross section (with particle_interference of
//    ObservableCoupling.h as the coefficient);
//  - the Coulomb-nuclear interference: direct and exchange coefficients, the cross term's sign, the order limit;
//  - the capture analysing power's sums; the mixed R_k and the term of an unobserved primary; the R-matrix-capture
//    and beta-delayed terms.
//   physics_observables_test
#include <cmath>
#include <complex>
#include <cstdio>
#include <string>
#include <vector>

#include "azure2/physics/CollisionMatrix.h"
#include "azure2/physics/ObservableCoupling.h"
#include "azure2/physics/Observables.h"

using namespace azure2::physics;
typedef std::complex<double> complex;

namespace {
int failures = 0;
void check(const std::string &what, bool ok, double got = 0.0, double want = 0.0) {
  std::printf("  %s  %s", ok ? "ok  " : "FAIL", what.c_str());
  if (!ok) std::printf("  (got %.17g, want %.17g)", got, want);
  std::printf("\n");
  if (!ok) failures++;
}
bool close(double a, double b, double rel) { return std::fabs(a - b) <= rel * std::fabs(b); }
bool close(complex a, complex b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }
}  // namespace

int main() {
  const complex i(0.0, 1.0);
  const double pi = kObservablePi;

  // One level at 2 MeV, J = 3/2, two neutral channels (no Coulomb or hard-sphere phase, S - B = 0).
  const double El = 2.0, E = 1.9, J = 1.5, geom = 37.0, i1i2 = 0.25;
  const double g[2] = {0.6, 0.35}, P[2] = {0.2, 0.7};
  const double sqrtP[2] = {std::sqrt(P[0]), std::sqrt(P[1])};
  const double shift[2] = {0.0, 0.0}, zero2[2] = {0.0, 0.0}, park[1] = {1.0};
  const complex lo[2] = {complex(0.0, P[0]), complex(0.0, P[1])};
  const complex one[2] = {complex(1.0, 0.0), complex(1.0, 0.0)};
  const char rad[2] = {'P', 'P'};
  LevelMatrixTables t;
  t.n = 1;
  t.nc = 2;
  t.E = &El;
  t.gamma = g;
  t.shift = shift;
  t.shift_d1 = zero2;
  t.shift_d2 = zero2;
  t.park_norm = park;
  t.lo = lo;
  t.boundary = zero2;
  t.rad_type = rad;
  ChannelTablesAtE ch;
  ch.coulomb_phase = one;
  ch.hard_sphere_phase = one;
  ch.sqrt_penetrability = sqrtP;
  const CollisionMatrix cm = collision_matrix(t, E, ch);
  check("the block is solved", !cm.singular);

  const double width[2] = {2.0 * P[0] * g[0] * g[0], 2.0 * P[1] * g[1] * g[1]};
  const double total_width = width[0] + width[1];
  const double bw = geom * (2.0 * J + 1.0) * i1i2 * width[0] * width[1] /
                    ((E - El) * (E - El) + 0.25 * total_width * total_width) / 100.0;
  PartialWaveSums sums;
  complex sum(0.0, 0.0);
  sums.add(J, 1, 2, cm.t(0, 1));
  add_partial_wave_cross_section(sums, geom, i1i2, sum);
  const double reaction = cross_section_barns(sum);
  check("reaction cross section = one-level Breit-Wigner", close(reaction, bw, 1e-13), reaction, bw);

  sums.clear();
  complex sum_el(0.0, 0.0);
  sums.add(J, 1, 1, cm.t(0, 0));
  add_partial_wave_cross_section(sums, geom, i1i2, sum_el);
  const double elastic = cross_section_barns(sum_el);
  double optical = 0.0;
  add_optical_theorem_term(J, cm.t(0, 0), optical);
  const double total = total_cross_section(geom, i1i2, optical);
  check("optical theorem: total = elastic + reaction", close(total, elastic + reaction, 1e-13), total,
        elastic + reaction);

  // Coherent by (J, l, l'): two pathways of one row add before the modulus; another l' is a row of its own.
  const complex a(0.3, -0.2), b(-0.1, 0.45), c(0.05, 0.6);
  sums.clear();
  sums.add(0.5, 0, 1, a);
  sums.add(0.5, 0, 1, b);
  sums.add(0.5, 0, 2, c);
  sums.add(1.5, 0, 1, c);
  check("three rows, in order of first appearance",
        sums.rows.size() == 3 && sums.rows[0].lp == 1 && sums.rows[1].lp == 2 && sums.rows[2].J == 1.5,
        (double)sums.rows.size(), 3.0);
  check("a row holds the sum of its pathways", sums.rows[0].T == a + b);
  complex s3(0.0, 0.0);
  add_partial_wave_cross_section(sums, 2.0, 0.5, s3);
  const double want3 = 2.0 * 0.5 * (2.0 * std::norm(a + b) + 2.0 * std::norm(c) + 4.0 * std::norm(c));
  check("sum of (2J+1) |sum T|^2", close(std::real(s3), want3, 1e-14), std::real(s3), want3);
  check("the sum is real", std::fabs(std::imag(s3)) < 1e-16, std::imag(s3), 0.0);
  subtract_partial_wave_cross_section(sums, 2.0, 0.5, s3);
  check("subtracting the same rows gives zero", std::abs(s3) < 1e-15, std::abs(s3), 0.0);

  // R-matrix capture and beta-delayed terms.
  check("R-matrix-capture elastic term 2 geom (2J+1) g Im T",
        close(rmc_elastic_term(3.0, 1.5, 0.25, complex(0.7, 0.2)), 2.0 * 3.0 * 4.0 * 0.25 * 0.2, 1e-15));
  complex sb(1.0, 0.0);
  add_beta_delayed_term(a, sb);
  check("beta-delayed term 25 |T|^2", close(std::real(sb), 1.0 + 25.0 * std::norm(a), 1e-15), std::real(sb),
        1.0 + 25.0 * std::norm(a));

  // Phase shift.
  const double deg = pi / 180.0;
  const complex eC = std::exp(i * 0.37);
  const complex U30 = std::exp(2.0 * i * 30.0 * deg);
  const complex T30 = eC * eC * (1.0 - U30);
  check("U from T of a diagonal pathway", close(collision_from_t_matrix(eC, T30), U30, 1e-14));
  check("phase shift of exp(2 i 30 deg)", close(phase_shift_degrees(U30, false), 30.0, 1e-14),
        phase_shift_degrees(U30, false), 30.0);
  const complex Um20 = std::exp(-2.0 * i * 20.0 * deg);
  check("a negative phase stays for distinguishable particles", close(phase_shift_degrees(Um20, false), -20.0, 1e-14),
        phase_shift_degrees(Um20, false), -20.0);
  check("and is wrapped to [0, 180) for identical ones", close(phase_shift_degrees(Um20, true), 160.0, 1e-14),
        phase_shift_degrees(Um20, true), 160.0);
  check("a positive phase is not wrapped", close(phase_shift_degrees(U30, true), 30.0, 1e-14));

  // Legendre coefficients: one s-wave pathway (l = l' = 0, J = s = s' = 1/2) is isotropic, 4 pi B_0 = sigma.
  {
    const complex T(0.4, -0.3);
    const double z = particle_interference(0, 0, 1, 0, 0, 1, 1, 1, 0);
    std::vector<complex> wl(1, complex(0.0, 0.0));
    add_interference_product(z, T, T, wl[0]);
    std::vector<double> B(1, std::real(wl[0]));
    scale_legendre_coefficients(B, geom, i1i2, 1.0);
    PartialWaveSums one_row;
    one_row.add(0.5, 0, 0, T);
    complex ssum(0.0, 0.0);
    add_partial_wave_cross_section(one_row, geom, i1i2, ssum);
    check("s wave: 4 pi B_0 = sigma", close(4.0 * pi * B[0], cross_section_barns(ssum), 1e-14), 4.0 * pi * B[0],
          cross_section_barns(ssum));
    std::vector<double> B2(2, 1.0);
    scale_legendre_coefficients(B2, 2.0, 0.5, 2.0);
    check("the identical-particle factor multiplies every order", close(B2[0], 2.0 / pi / 100.0, 1e-15) &&
          B2[1] == B2[0], B2[0], 2.0 / pi / 100.0);
    complex w(0.0, 0.0);
    add_interference_product(0.5, a, b, w);
    check("interference product z T1 conj(T2)", close(w, 0.5 * a * std::conj(b), 1e-15));
  }

  // Unobserved primary, observed secondary.
  check("mixed R_k without mixing", upos_mixed_rk(0.3, 0.7, -0.2, 0.0) == 0.3);
  check("mixed R_k with delta = 1", close(upos_mixed_rk(0.3, 0.7, -0.2, 1.0), 0.5 * (0.3 + 1.4 - 0.2), 1e-15));
  {
    const double want = std::real(0.8 * a * std::conj(b)) * std::sqrt(5.0) / 4.0 * 0.6;
    check("term of an unobserved primary, order 2", close(upos_term(0.8, a, b, 2, 0.6), want, 1e-15),
          upos_term(0.8, a, b, 2, 0.6), want);
  }

  // Coulomb-nuclear interference of elastic scattering.
  {
    CoulombInterference ci;
    ci.add(0.0, 0, 0.25, a);
    ci.add(1.0, 0, 0.75, b);
    ci.add(1.0, 2, 0.75, c);
    ci.add(1.0, 0, 0.75, c);
    check("terms by (s, l), in order of first appearance",
          ci.terms.size() == 3 && ci.terms[1].s == 1.0 && ci.terms[1].l == 0 && ci.terms[2].l == 2,
          (double)ci.terms.size(), 3.0);
    check("a term holds statSpin conj(T) summed",
          close(ci.terms[1].c, 0.75 * std::conj(b) + 0.75 * std::conj(c), 1e-15));

    ElasticInterference e;
    int above = elastic_interference(ci, geom, false, false, 0.0, 10, e);
    const complex cI = i / std::sqrt(pi) * geom / 100.0;
    check("no order above the limit", above == -1, (double)above, -1.0);
    check("coulomb = geom / 100", e.coulomb == geom / 100.0, e.coulomb, geom / 100.0);
    check("no cross term or exchange without spin-dependent exchange", e.cross == 0.0 && e.exchange.empty());
    check("D_0 = (i geom / (100 sqrt(pi))) sum_s c_{s,0}",
          e.direct.size() == 3 && close(e.direct[0], cI * (ci.terms[0].c + ci.terms[1].c), 1e-14));
    check("D_1 is zero, D_2 one term", e.direct[1] == complex(0.0, 0.0) && close(e.direct[2], cI * ci.terms[2].c,
                                                                                    1e-14));

    ElasticInterference x;
    above = elastic_interference(ci, geom, true, true, 0.5, 10, x);
    check("identical particles: D_l doubled", close(x.direct[0], 2.0 * e.direct[0], 1e-14));
    check("X_0 = sum_s (-1)^s c_{s,0} with the same factor",
          x.exchange.size() == 3 && close(x.exchange[0], 2.0 * cI * (ci.terms[0].c - ci.terms[1].c), 1e-14));
    check("cross term of spin 1/2: -2 c_0 / (2j + 1)", close(x.cross, -x.coulomb, 1e-15), x.cross, -x.coulomb);
    ElasticInterference x1;
    elastic_interference(ci, geom, true, true, 1.0, 10, x1);
    check("cross term of spin 1: +2 c_0 / 3", close(x1.cross, 2.0 * x1.coulomb / 3.0, 1e-15), x1.cross,
          2.0 * x1.coulomb / 3.0);
    ElasticInterference lim;
    above = elastic_interference(ci, geom, false, false, 0.0, 1, lim);
    check("an l above the limit is returned", above == 2, (double)above, 2.0);
  }

  // Capture analysing power: Re T T'* into the denominator, Re (i T T'*) into the numerator.
  {
    double den = 0.0, num = 0.0;
    add_capture_ay_term(2.0, 3.0, complex(1.0, 0.0), i, den, num);  // T T'* = -i: Re = 0, Re (i T T'*) = 1
    check("T = 1, T' = i: only the numerator", den == 0.0 && num == 3.0, num, 3.0);
    add_capture_ay_term(2.0, 3.0, a, a, den, num);  // T T'* = |T|^2
    check("T = T': only the denominator", close(den, 2.0 * std::norm(a), 1e-15) && num == 3.0, den,
          2.0 * std::norm(a));
    double den0 = 1.0, num0 = 1.0;
    add_capture_ay_term(0.0, 0.0, a, b, den0, num0);
    check("zero coefficients add nothing", den0 == 1.0 && num0 == 1.0);
  }

  std::printf("%s: %d failure(s)\n", failures ? "FAILED" : "passed", failures);
  return failures ? 1 : 0;
}
