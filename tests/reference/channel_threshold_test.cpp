/*!
 * The shared threshold helper ChannelFunc (include/ChannelFunc.h): shift
 * function S, penetrability P and dS/dE of a particle channel at a channel
 * energy at, just below and just above threshold, as a level energy sitting
 * there asks for them (CNuc boundary conditions, shift functions and parameter
 * transformations, the THM vertex, the external-capture entrance phase).
 *
 * For p+7Li, 6Li+d, 12C+12C (charged, eta ~ 1e3-1e4 at 1e-9 MeV) and 17O+n
 * (neutral), l = 0-3:
 *   - e = -1e-9, 0, +1e-9 MeV: S, P, dS/dE finite, P >= 0 (and P = 0 for the
 *     charged pairs), each call below 50 ms;
 *   - S continuous: |S(+-1e-9) - S(0)| <= 1e-6 (1 + |S(0)|);
 *   - dS/dE continuous and right: the three values agree to 1e-4 relative and
 *     match a one-sided difference of S from below with steps 1e-4/1e-5 MeV
 *     (Richardson), to 1e-3 relative -- the old dS/dE at a level exactly at
 *     threshold was 0 (ShftFunc mirrored the probe above threshold), and for
 *     p+7Li the slope of ShftFunc's threshold interpolation was 60 % off
 *     (l = 1).  Not for a neutral l <= 1 channel, where S - S(0) ~ sqrt(B)
 *     below threshold (dS/dE -> infinity from below): finite only;
 *   - the eta = 100 line, where the reflection S = 2 S(0) - S(-e) hands over
 *     to the Coulomb functions: S continuous to 1e-5 (1 + |S|);
 *   - away from threshold the helper is the routine it replaces, bit for bit:
 *     CoulFunc::PEShift/Penetrability/PEShift_dE at +1 MeV, ShftFunc and
 *     ShftFunc::EnergyDerivative at -1 MeV.
 *
 * Run:  tests/reference/channel_threshold_test [-v]   (ctest: channel_threshold)
 */
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

#include "ChannelFunc.h"
#include "Config.h"
#include "Constants.h"
#include "CoulFunc.h"
#include "NucLine.h"
#include "PPair.h"
#include "ShftFunc.h"

Config *g_config = nullptr;

namespace {

struct PairSpec {
  const char *name;
  int z1, z2;
  double m1, m2, a, sepE;
};

const PairSpec kPairs[] = {
    {"7Li+p", 3, 1, 7.016003, 1.007276, 4.15, 17.25},
    {"6Li+d", 3, 1, 6.015123, 2.014102, 4.5, 22.28},
    {"12C+12C", 6, 6, 12.0, 12.0, 6.0, 13.933},
    {"17O+n", 8, 0, 16.999131, 1.008665, 5.0, 4.143},
};

PPair makePair(const PairSpec &s) {
  std::ostringstream line;
  line.precision(15);
  line << "1 1 0 0 1 1 1 0 1 1 0 0 0 1 0 1 0 " << s.m1 << " " << s.m2 << " " << s.z1 << " "
       << s.z2 << " " << s.sepE << " " << s.sepE << " 0 1 0 0 " << s.a << " 0 0 0";
  std::istringstream in(line.str());
  NucLine nl(in);
  return PPair(nl);
}

int failures = 0;
bool verbose = false;

void check(bool ok, const std::string &what) {
  if (!ok) {
    failures++;
    std::printf("FAIL  %s\n", what.c_str());
  } else if (verbose) {
    std::printf("ok    %s\n", what.c_str());
  }
}

double relDiff(double a, double b) {
  return std::fabs(a - b) / std::max(std::fabs(b), 1e-300);
}

}  // namespace

int main(int argc, char **argv) {
  verbose = (argc > 1 && std::strcmp(argv[1], "-v") == 0);
  double maxMs = 0.0;
  for (const PairSpec &spec : kPairs) {
    PPair pair = makePair(spec);
    const bool charged = spec.z1 * spec.z2 != 0;
    for (int l = 0; l <= 3; l++) {
      char tag[96];
      std::snprintf(tag, sizeof tag, "%s l=%d", spec.name, l);
      const std::string t(tag);
      ChannelFunc cf(&pair, false);

      const double es[3] = {-1e-9, 0.0, 1e-9};
      double S[3], P[3], D[3];
      for (int i = 0; i < 3; i++) {
        auto t0 = std::chrono::steady_clock::now();
        S[i] = cf.Shift(l, es[i]);
        P[i] = cf.Penetrability(l, es[i]);
        D[i] = cf.ShiftDerivative(l, es[i]);
        double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        maxMs = std::max(maxMs, ms);
        char e[64];
        std::snprintf(e, sizeof e, " e=%g", es[i]);
        if (verbose)
          std::printf("      %s%s: S=%.12g P=%.3g dS/dE=%.10g (%.2f ms)\n", tag, e, S[i], P[i], D[i], ms);
        check(std::isfinite(S[i]) && std::isfinite(P[i]) && std::isfinite(D[i]), t + e + " finite");
        check(P[i] >= 0.0, t + e + " P >= 0");
        if (charged) check(P[i] == 0.0, t + e + " P = 0 (closed)");
        check(ms < 50.0, t + e + " fast (" + std::to_string(ms) + " ms)");
      }
      // Neutral l = 0: S = -kappa a below, 0 above; |S(-1e-9)| ~ 1e-5.
      const double sTol = (charged || l >= 1) ? 1e-6 : 1e-4;
      check(std::fabs(S[0] - S[1]) <= sTol * (1 + std::fabs(S[1])) &&
                std::fabs(S[2] - S[1]) <= sTol * (1 + std::fabs(S[1])),
            t + " S continuous through threshold");
      if (charged || l >= 2)
        check(relDiff(D[0], D[1]) < 1e-4 && relDiff(D[2], D[1]) < 1e-4,
              t + " dS/dE continuous through threshold");
      // One-sided reference from below, Richardson-extrapolated.
      const double s0 = cf.Shift(l, 0.0);
      const double d1 = (s0 - cf.Shift(l, -1e-4)) / 1e-4;
      const double d2 = (s0 - cf.Shift(l, -1e-5)) / 1e-5;
      const double ref = (10.0 * d2 - d1) / 9.0;
      if (verbose) std::printf("      %s dS/dE(0) reference %.10g\n", tag, ref);
      // A neutral l <= 1 channel has S - S(0) ~ sqrt(B) below threshold: dS/dE
      // is infinite from below there, and only finiteness is asked of it.
      if (charged || l >= 2) {
        check(relDiff(D[1], ref) < 1e-3, t + " dS/dE(0) matches the difference of S from below");
        check(D[1] != 0.0, t + " dS/dE(0) nonzero");
      }

      if (charged) {
        // The eta = 100 line.
        const double eta1 = std::sqrt(uconv / 2.) * fstruc * spec.z1 * spec.z2 *
                            std::sqrt(pair.GetRedMass());
        const double e100 = std::pow(eta1 / ChannelFunc::kEtaClosed, 2.0);
        const double sb = cf.Shift(l, e100 * (1 - 1e-9));
        const double sa = cf.Shift(l, e100 * (1 + 1e-9));
        if (verbose) std::printf("      %s eta=100 at e=%.4g MeV: S %.12g | %.12g\n", tag, e100, sb, sa);
        check(std::isfinite(sb) && std::isfinite(sa) &&
                  std::fabs(sa - sb) <= 1e-5 * (1 + std::fabs(sa)),
              t + " S continuous at eta = 100");
        check(cf.Penetrability(l, e100 * (1 - 1e-9)) == 0.0 &&
                  cf.Penetrability(l, e100 * (1 + 1e-9)) >= 0.0,
              t + " P at eta = 100");
      }

      // Bit for bit the old routines away from threshold.
      CoulFunc coul(&pair, false);
      ShftFunc shft(&pair);
      check(cf.Shift(l, 1.0) == coul.PEShift(l, spec.a, 1.0), t + " S(+1) = CoulFunc::PEShift");
      check(cf.Penetrability(l, 1.0) == coul.Penetrability(l, spec.a, 1.0),
            t + " P(+1) = CoulFunc::Penetrability");
      check(cf.ShiftDerivative(l, 1.0) == coul.PEShift_dE(l, spec.a, 1.0),
            t + " dS/dE(+1) = CoulFunc::PEShift_dE");
      check(cf.Shift(l, -1.0) == shft(l, spec.sepE - 1.0), t + " S(-1) = ShftFunc");
      check(cf.ShiftDerivative(l, -1.0) == shft.EnergyDerivative(l, spec.sepE - 1.0),
            t + " dS/dE(-1) = ShftFunc::EnergyDerivative");
      check(cf.Penetrability(l, -1.0) == 0.0, t + " P(-1) = 0");
    }
  }
  std::printf("slowest threshold call %.2f ms\n", maxMs);
  if (failures) {
    std::printf("%d check(s) failed\n", failures);
    return 1;
  }
  std::printf("all channel-threshold checks passed\n");
  return 0;
}
