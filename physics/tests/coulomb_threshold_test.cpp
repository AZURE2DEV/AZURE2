// physics/coulomb just above a charged-particle threshold (docs/upgrade/coulomb.md section 5): the shift function
// S = rho (F F' + G G') / (F^2 + G^2) of the engine's method (coulsm) against its exact limit at E -> 0,
// S_l(0) = -l - (x/2) K_2l(x) / K_2l+1(x), x = (8 eta rho)^1/2 (zero_energy_shift), for 12C+12C, 16O+16O and
// alpha+12C, l = 0-4, 20 energies from 1 eV to 2.5 keV (c.m.) above threshold.
//
// Where Michel's continued fractions give up (eta >~ 230: below 3.6 keV for 12C+12C, 210 eV for alpha+12C) coulsm
// falls back to the zero-energy limit G -> 2 (rho/pi)^1/2 K_1(z) e^(pi eta), z = x, then the recursion in l; eta rho
// does not depend on E, so for l = 0 S is S_l(0) itself there.  Where the continued fractions work S(E) - S_l(0) =
// S'(0) E + O(E^2), the slope the limit leaves out.  The bound is therefore
//   |S(E) - S_l(0)| <= 1e-10 |S_l(0)| + 1.5 s E,
// with s = (S(E_ref) - S_l(0)) / E_ref the chord slope to E_ref = 50 keV, where the continued fractions work for all
// three channels (S is nearly linear in E there: for 12C+12C, l = 0, the chord slope is 0.648 MeV^-1 at 3.6 keV and
// 0.650 at 50 keV; the factor 1.5 covers the curvature), and 1e-10 the precision of coulsm.  Fails on K_0, K_1 from
// their power series (S of 12C+12C, l = 0, between 1 eV and 2.5 keV: +4.6 to +11.2 instead of -9.594).  Also: P =
// rho / (F^2 + G^2) finite, >= 0, and below 1e-30 everywhere on the scan.
#include <cmath>
#include <cstdio>

#include "azure2/physics/ChannelFunctions.h"
#include "azure2/physics/coulomb/Coulomb.h"

namespace ph = azure2::physics;
namespace pc = azure2::physics::coulomb;

int main() {
  struct Channel {
    const char *name;
    int z1, z2;
    double mu_u, a_fm;
  };
  const Channel channels[] = {{"12C+12C", 6, 6, 6.0, 6.5},
                              {"16O+16O", 8, 8, 15.9949 / 2.0, 7.0},
                              {"alpha+12C", 2, 6, 4.0026 * 12.0 / 16.0026, 5.5}};
  const double lowest_MeV = 1e-6, highest_MeV = 2.5e-3, referenceMeV = 0.05;
  int failures = 0, cases = 0;
  for (const Channel &c : channels) {
    ph::ChannelKinematics kin;
    kin.z1 = c.z1;
    kin.z2 = c.z2;
    kin.mu_u = c.mu_u;
    auto shiftAt = [&](int l, double E) {
      const double r = ph::rho(kin, c.a_fm, E);
      return ph::shift(pc::coulomb(l, ph::eta(kin, E), r), r);
    };
    for (int l = 0; l <= 4; l++) {
      cases++;
      const double cFm = ph::eta(kin, highest_MeV) * ph::rho(kin, c.a_fm, highest_MeV) / c.a_fm;
      const double S0 = ph::zero_energy_shift(cFm, c.a_fm, l);
      const double slope = (shiftAt(l, referenceMeV) - S0) / referenceMeV;
      bool ok = std::isfinite(S0) && std::isfinite(slope);
      double worst = 0.0, worstE = 0.0, worstS = 0.0;
      for (int i = 0; i < 20; i++) {
        const double E = lowest_MeV * std::pow(highest_MeV / lowest_MeV, i / 19.0);
        const double r = ph::rho(kin, c.a_fm, E);
        const pc::CoulombWaves w = pc::coulomb(l, ph::eta(kin, E), r);
        const double S = ph::shift(w, r), P = ph::penetrability(w, r);
        const double bound = 1e-10 * std::fabs(S0) + 1.5 * std::fabs(slope) * E;
        const double d = std::fabs(S - S0);
        if (w.status != pc::Status::ok || !std::isfinite(S) || !(d <= bound) || !std::isfinite(P) || P < 0.0 ||
            P > 1e-30)
          ok = false;
        if (!(d / bound <= worst)) worst = d / bound, worstE = E, worstS = S;
      }
      if (!ok) failures++;
      std::printf("%s %-9s l=%d  S(0) %.10f  worst |S - S(0)| / bound %.3g at %.3g keV (S %.10f)\n",
                  ok ? "ok  " : "FAIL", c.name, l, S0, worst, worstE * 1e3, worstS);
    }
  }
  std::printf("%d of %d cases fail\n", failures, cases);
  return failures ? 1 : 0;
}
