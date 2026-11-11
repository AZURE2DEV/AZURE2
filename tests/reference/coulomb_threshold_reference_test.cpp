// Coulomb functions just above a charged-particle threshold, against the exact zero-energy limit (2026-11-06).
//
// As E -> 0 at a fixed channel radius, eta rho stays constant and the Coulomb equation goes over into a Bessel
// equation; the irregular solution becomes G -> 2 sqrt(rho/pi) K_{2l+1}(z) e^{pi eta}, z = sqrt(8 eta rho), so
//     S(E -> 0) = rho G'/G = -l - (z/2) K_{2l}(z) / K_{2l+1}(z).
// S(E) moves away from it linearly, by ~0.6 E[MeV] for these channels, so within a few keV of the threshold
// AZURE2's shift function must equal this limit to ~1e-3.  Where Michel's continued fraction gives up
// (eta >~ 250), cwfcomp_accurate falls back to that very limit; it used to evaluate K_0 and K_1 from their power
// series, which cancel to nothing at z ~ 20: for 12C+12C the shift came out as +4.6 .. +11.2 instead of -9.594
// within ~3 keV of the threshold.  Checked here: 12C+12C and 16O+16O (z ~ 20, ~31) and alpha+12C (z ~ 7.4), l = 0..4,
// at 20 energies from 1 eV to 2.5 keV above the threshold, through CoulFunc::PEShift (the default route).
#include "CoulFunc.h"
#include "GSLException.h"
#include "NucLine.h"
#include "PPair.h"
#include "Kinematics.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <cmath>
#include <cstdio>
#include <sstream>

class Config;
Config *g_config = nullptr;  // defined by the executable; the library reads it

static int failures = 0;

int main() {
  gsl_set_error_handler(&GSLException::GSLErrorHandler);
  struct Chan { const char *name; int z1, z2; double m1, m2, a; };
  const Chan chans[] = {{"12C+12C", 6, 6, 12.0, 12.0, 6.5}, {"16O+16O", 8, 8, 15.9949, 15.9949, 7.0},
                        {"alpha+12C", 2, 6, 4.0026, 12.0, 5.5}};
  for (const Chan &c : chans) {
    std::ostringstream ln;
    ln << "0.5 1 1 1 1 1 0.5 0 1 1 0 1 0 1 0 1 0 " << c.m1 << " " << c.m2 << " " << c.z1 << " " << c.z2
       << " 0 0 0 0 0 0 " << c.a << " 0 0 0";
    std::istringstream is(ln.str());
    NucLine line(is);
    PPair pair(line);
    const Kinematics::Channel kin(&pair);
    for (int l = 0; l <= 4; l++) {
      CoulFunc coul(&pair, false);
      double worst = 0.0, worstE = 0.0, worstS = 0.0, worstS0 = 0.0;
      bool threw = false;
      for (int i = 0; i < 20; i++) {
        const double E = 1e-6 * std::pow(2.5e-3 / 1e-6, i / 19.0);
        const double eta = kin.Eta(E), rho = kin.Rho(c.a, E);
        const double z = std::sqrt(8.0 * eta * rho);
        const double S0 = -l - 0.5 * z * gsl_sf_bessel_Kn(2 * l, z) / gsl_sf_bessel_Kn(2 * l + 1, z);
        double S;
        try {
          S = coul.PEShift(l, c.a, E);
        } catch (...) {
          threw = true;
          continue;
        }
        const double d = std::isfinite(S) ? std::fabs(S - S0) / std::fabs(S0) : INFINITY;
        if (d > worst) worst = d, worstE = E, worstS = S, worstS0 = S0;
      }
      const bool ok = !threw && worst < 1e-3;
      std::printf("  %s  %-9s l=%d: worst |S - S(0)|/|S(0)| = %.2e (E %.3g keV: S %.6f, limit %.6f)%s\n", ok ? "ok  " : "FAIL",
                  c.name, l, worst, worstE * 1e3, worstS, worstS0, threw ? "  [threw]" : "");
      if (!ok) failures++;
    }
  }
  std::printf(failures ? "FAILED: %d\n" : "PASS\n", failures);
  return failures ? 1 : 0;
}
