// Coupling coefficients of the observables (ObservableCoupling.h; docs/upgrade/cnuc_split.md section 9).  The
// expressions and the order of their factors are those of CNuc::CalcAngularDists and of the capture-A_y build (now
// BuildPathwayTable) before step 5, with AngCoeff's Clebsch-Gordan and Racah coefficients replaced by AngularCoupling's.
#include "azure2/physics/ObservableCoupling.h"

#include <cmath>
#include <cstdlib>

#include "azure2/physics/AngularCoupling.h"

namespace azure2 {
namespace physics {

namespace {

// (-1)^(twice/2) for an even argument (an integer exponent passed doubled).
double PhaseOfTwice(int twice) { return ((std::abs(twice) / 2) % 2 == 0) ? 1.0 : -1.0; }

// x^ = (2x + 1)^(1/2) for x passed doubled.
double HatOfTwice(int two_x) { return std::sqrt((double)two_x + 1.0); }

}  // namespace

double z_bar(int l1, int two_j1, int l2, int two_j2, int two_s, int k) {
  return std::sqrt(2. * l1 + 1.) * std::sqrt(2. * l2 + 1.) * HatOfTwice(two_j1) * HatOfTwice(two_j2) *
         clebsch_gordan(2 * l1, 0, 2 * l2, 0, 2 * k, 0) * racah_w(2 * l1, two_j1, 2 * l2, two_j2, two_s, 2 * k);
}

double z_bar_photon(int L1, int two_j1, int L2, int two_j2, int two_jf, int k) {
  return std::sqrt(2. * L1 + 1.) * std::sqrt(2. * L2 + 1.) * HatOfTwice(two_j1) * HatOfTwice(two_j2) *
         clebsch_gordan(2 * L1, 2, 2 * L2, -2, 2 * k, 0) * racah_w(2 * L1, two_j1, 2 * L2, two_j2, two_jf, 2 * k);
}

double particle_interference(int l1, int l1p, int two_j1, int l2, int l2p, int two_j2, int two_s, int two_sp, int k) {
  const double z1 = z_bar(l1, two_j1, l2, two_j2, two_s, k);
  const double z2 = z_bar(l1p, two_j1, l2p, two_j2, two_sp, k);
  return PhaseOfTwice(two_sp - two_s) / 4. * z1 * z2;
}

double photon_interference(int l1, int L1, int two_j1, int l2, int L2, int two_j2, int two_s, int two_jf, int k) {
  const double z1 = z_bar(l1, two_j1, l2, two_j2, two_s, k);
  const double z2 = z_bar_photon(L1, two_j1, L2, two_j2, two_jf, k);
  return PhaseOfTwice(2 + two_s - two_jf) / 4. * z1 * z2;
}

double upos_interference(int l1, int l1p, int two_j1, int l2, int two_j2, int two_s, int two_sp, int two_sp2,
                         int two_jb, int two_jB, int k) {
  return PhaseOfTwice(two_sp - two_sp2) * PhaseOfTwice(2 * k + two_sp2 - two_sp) * (two_j1 + 1.) * (two_j2 + 1.) *
         std::pow((2. * l1 + 1.) * (two_jB + 1.) * (two_sp + 1.) * (two_sp2 + 1.), 0.5) *
         clebsch_gordan(2 * k, 0, 2 * l1, 0, 2 * l2, 0) * racah_w(2 * k, two_jB, two_sp2, two_jb, two_jB, two_sp) *
         racah_w(2 * k, two_sp, two_j2, 2 * l1p, two_sp2, two_j1) * racah_w(2 * k, two_j1, 2 * l2, two_s, two_j2, 2 * l1);
}

CaptureAyCoefficients capture_ay_coefficients(const CapturePathway &t1, const CapturePathway &t2, int k, int two_x,
                                              int two_a, int two_c) {
  CaptureAyCoefficients out;
  // Eq. (15): the empty bracket, [ ] = 1/2 [1 + (-1)^(L + p + L' + p' + k)].
  if ((t1.L + t1.p + t2.L + t2.p + k) % 2 != 0) return out;
  // Both coefficients carry (l 0, l' 0 | k 0) and (L 1, L' -1 | k 0).
  const double cgL = clebsch_gordan(2 * t1.l, 0, 2 * t2.l, 0, 2 * k, 0);
  if (std::fabs(cgL) < 1.e-12) return out;
  const double cgG = clebsch_gordan(2 * t1.L, 2, 2 * t2.L, -2, 2 * k, 0);
  if (std::fabs(cgG) < 1.e-12) return out;
  const double wG = racah_w(2 * t1.L, t1.two_b, 2 * t2.L, t2.two_b, two_c, 2 * k);
  if (std::fabs(wG) < 1.e-12) return out;

  const double common = cgL * cgG * wG * std::sqrt(2. * t1.l + 1.) * std::sqrt(2. * t2.l + 1.) *
                        std::sqrt(2. * t1.L + 1.) * std::sqrt(2. * t2.L + 1.) * (t1.two_b + 1.) * (t2.two_b + 1.);

  // Eq. (20). Only the channel-spin diagonal contributes.
  if (t1.two_s == t2.two_s) {
    const double wL = racah_w(2 * t1.l, t1.two_b, 2 * t2.l, t2.two_b, t1.two_s, 2 * k);
    out.a = PhaseOfTwice(t1.two_s - two_c + 2) * common * wL;
  }

  // Eq. (21). Channel-spin off-diagonal terms are kept; the prefactor is 3 sqrt(x) x^ k^ / [(x+1) k (k+1)]^(1/2).
  if (k >= 1) {
    const double wS = racah_w(two_x, t1.two_s, two_x, t2.two_s, two_a, 2);
    if (std::fabs(wS) >= 1.e-12) {
      const double x9j = wigner_9j(2 * t1.l, t1.two_s, t1.two_b, 2 * t2.l, t2.two_s, t2.two_b, 2 * k, 2, 2 * k);
      if (std::fabs(x9j) >= 1.e-14) {
        const double x = two_x / 2.0;
        const double prefactor =
            3.0 * std::sqrt(x) * HatOfTwice(two_x) * std::sqrt(2. * k + 1.) / std::sqrt((x + 1.0) * (double)k * ((double)k + 1.0));
        out.b = prefactor * common * HatOfTwice(t1.two_s) * HatOfTwice(t2.two_s) * wS * x9j *
                PhaseOfTwice(two_a - two_x + two_c - t1.two_b - t1.two_s + 2 * t1.l) *
                PhaseOfTwice(t1.two_s - t2.two_s);
      }
    }
  }
  return out;
}

}  // namespace physics
}  // namespace azure2
