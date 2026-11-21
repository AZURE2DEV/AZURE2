// The capture from the interior of a hybrid potential (docs/upgrade/interior_potential.md section 6).
//
// On the pair of tests/upgrade/interior_bound (alpha + 12C, the Woods-Saxon of Koonin, Tombrello and Fox in every l),
// E1 from l = 1 into an l = 0 final state, ECIntegral at the channel radii 4.0, 5.35 and 6.5 fm with the boundary
// constants 0 and -1.5, against tests/upgrade/interior_capture/reference.csv (make_reference.py: scipy and mpmath, not
// the engine's code):
//  (a) each term with the potential's interior over the same term from outside the radius alone (the final channel's
//      potential without its boundary constant): the non-resonant term, a level's channel term above threshold
//      (complex) and below it, for the potential's own 0+ bound state and for a final state that is none;
//  (b) the non-resonant term at 4.0 and 6.5 fm over the one at 5.35 fm: for the potential's bound state it does not
//      depend on the radius (the potential model's direct capture), in the long-wavelength form and in the exact one;
//  (c) M1 between two l = 1 states of one potential, no level: the amplitude from outside the radius is cancelled by
//      the interior (the wave functions are orthogonal);
//  (d) the interior charge q: the amplitude is linear in q, and q equal to the pair's own E1 charge number
//      Z1 m2/M - Z2 m1/M doubles the interior part;
//  (e) an unbound final state (a narrow resonance at E_f = 0.3, 1.0 and 2.5 MeV): the interior part of the
//      non-resonant and of the channel term over the interior part of the same term for the bound final state of
//      7.16 MeV (complex); and the interior charge is linear there as well.
//
//   interior_capture_test <interior_bound.azr> <reference.csv>      (ctest: upgrade_interior_capture)
#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>

#include <gsl/gsl_errno.h>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "AChannel.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "Constants.h"
#include "ECAmplitudeCache.h"
#include "ECIntegral.h"
#include "JGroup.h"
#include "PPair.h"
#include "SetUpForTests.h"
#include "azure2/physics/HybridPotential.h"

namespace phys = azure2::physics;
typedef std::complex<double> cplx;

namespace {
int failures = 0;
void check(const std::string &what, bool ok) {
  std::printf("  %-4s %s\n", ok ? "ok" : "FAIL", what.c_str());
  if (!ok) failures++;
}
std::string num(double v) {
  char buf[32];
  std::snprintf(buf, sizeof buf, "%.2g", v);
  return buf;
}
}  // namespace

int main(int argc, char **argv) {
  if (argc < 3) return 2;
  gsl_set_error_handler_off();  // a quadrature that stops at rounding keeps its value (tolerance 1e-6 below)
  std::ostringstream sink;
  Config cfg(sink);
  cfg.configfile = argv[1];
  if (cfg.ReadConfigFile() < 0) throw std::runtime_error(std::string("cannot read ") + argv[1]);
  cfg.ResolveMatrixRoute();
  InitializeECAmplitudeCache(cfg);
  CompoundSetUp setUp = test_setup::WithoutData(cfg);
  if (setUp.refusal != CompoundSetUp::Refusal::kNone) throw std::runtime_error("the set-up failed");
  const CompoundSystem &cn = *setUp.compound;
  const PPair *projectPair = nullptr;
  for (int j = 1; j <= cn.NumJGroups() && !projectPair; j++)
    for (int ch = 1; ch <= cn.GetJGroup(j)->NumChannels() && !projectPair; ch++)
      if (cn.GetJGroup(j)->GetChannel(ch)->GetRadType() == 'P')
        projectPair = cn.GetPair(cn.GetJGroup(j)->GetChannel(ch)->GetPairNum());
  if (!projectPair || !projectPair->GetHybridPotential(1)) throw std::runtime_error("no pair with a potential");
  const double sep = projectPair->GetSepE() + projectPair->GetExE();
  phys::HybridPotential base = *projectPair->GetHybridPotential(1);
  if (argc > 3) base.step_fm = std::stod(argv[3]);

  // The pair with boundary constant b in the entrance l and, when `finalInterior`, in the final l as well.
  auto makePair = [&](double b, int li, int lf, bool finalInterior, double q = 0.0) {
    PPair pair(*projectPair);
    phys::HybridPotential in = base, out = base;
    in.interior_boundary = b;
    in.interior_charge_e1 = q;
    if (finalInterior) out.interior_boundary = b; else out.interior_boundary.reset();
    pair.SetHybridPotential(base);
    if (li == lf) pair.SetHybridPotentialsByL({{li, finalInterior ? in : out}});
    else pair.SetHybridPotentialsByL({{li, in}, {lf, out}});
    return pair;
  };
  Config lw(cfg), exact(cfg);
  lw.paramMask |= Config::USE_LONGWAVELENGTH_APPROX;
  exact.paramMask &= ~Config::USE_LONGWAVELENGTH_APPROX;
  auto e1 = [&](const PPair &pair, const Config &c, double a, double E, double B, bool channel) {
    ECIntegral ec(&pair, c, a);
    ec.SetRelativeTolerance(1e-6);
    return ec(1, 0, 0.0, 0.0, 1.0, 0.0, 1, 'E', sep + E, sep - B, channel);
  };

  // (a), (b)
  std::ifstream in(argv[2]);
  if (!in) throw std::runtime_error(std::string("cannot read ") + argv[2]);
  std::string line;
  std::getline(in, line);
  std::map<std::string, double> worst;
  std::map<std::string, int> count;
  double B0 = 0.0, worstExact = 0.0;
  while (std::getline(in, line)) {
    for (char &c : line) if (c == ',') c = ' ';
    std::istringstream row(line);
    std::string kind;
    double a, b, B, E, re, im;
    row >> kind >> a >> b >> B >> E >> re >> im;
    if (kind == "bound") { B0 = B; continue; }
    const PPair with = makePair(b, 1, 0, true), without = makePair(b, 1, 0, false);
    double d;
    if (kind == "unbound") {  // B = -E_f; both terms against the one reference number
      d = 0.0;
      for (bool channel : {false, true}) {
        const cplx unbound = e1(with, lw, a, E, B, channel) - e1(without, lw, a, E, B, channel);
        const cplx bound = e1(with, lw, a, E, 7.16, channel) - e1(without, lw, a, E, 7.16, channel);
        d = std::fmax(d, std::abs(unbound / bound - cplx(re, im)) / std::abs(cplx(re, im)));
      }
    } else if (kind == "total") {
      d = std::fabs(std::abs(e1(with, lw, a, E, B, false) / e1(with, lw, 5.35, E, B, false)) / std::fabs(re) - 1.0);
      if (B == B0)
        worstExact = std::fmax(worstExact, std::fabs(std::abs(e1(with, exact, a, E, B, false) / e1(with, exact, 5.35, E, B, false)) - 1.0));
    } else {
      const bool channel = kind != "direct";
      const cplx ratio = e1(with, lw, a, E, B, channel) / e1(without, lw, a, E, B, channel);
      d = std::abs(ratio - cplx(re, im)) / std::abs(cplx(re, im));
    }
    if (!(d < 1e-4)) std::printf("    %s a=%g b=%g B=%g E=%g: off by %.3g\n", kind.c_str(), a, b, B, E, d);
    worst[kind] = std::fmax(worst[kind], std::isfinite(d) ? d : 1.0);
    count[kind]++;
  }
  check("(a) the non-resonant term with the interior over the exterior one: worst " + num(worst["direct"]) + " of " +
            std::to_string(count["direct"]) + " (limit 1e-4)", count["direct"] == 72 && worst["direct"] < 1e-4);
  check("(a) a level's channel term above threshold: worst " + num(worst["channel"]) + " of " +
            std::to_string(count["channel"]) + " (limit 1e-4)", count["channel"] == 72 && worst["channel"] < 1e-4);
  check("(a) the channel term below threshold: worst " + num(worst["closed"]) + " of " + std::to_string(count["closed"]) +
            " (limit 1e-4)", count["closed"] == 18 && worst["closed"] < 1e-4);
  check("(e) unbound final state, the interior part over that of the bound one, both terms: worst " + num(worst["unbound"]) +
            " of " + std::to_string(count["unbound"]) + " (limit 1e-4)", count["unbound"] == 48 && worst["unbound"] < 1e-4);
  check("(b) the non-resonant term at 4.0 and 6.5 fm over 5.35 fm: worst " + num(worst["total"]) + " of " +
            std::to_string(count["total"]) + " (limit 1e-4)", count["total"] == 48 && worst["total"] < 1e-4);
  check("(b) exact form, the potential's bound state: independent of the radius to " + num(worstExact) + " (limit 1e-4)",
        worstExact < 1e-4);

  // (c) M1, l = 1 -> the potential's own l = 1 bound state (near 9.61 MeV), J 1 -> 1, no level.
  {
    const double a = 5.35, mu = projectPair->GetRedMass() * uconv;
    auto fMinusS = [&](double B) {
      return phys::RegularSolution(base, projectPair->GetZ(1), projectPair->GetZ(2), mu, 1, 0.0, a, -B).log_derivative -
             phys::HybridBoundShift(base, projectPair->GetZ(1), projectPair->GetZ(2), mu, 1, 0.0, a, B);
    };
    double lo = 9.55, hi = 9.66;
    for (int i = 0; i < 60; i++) {
      const double mid = 0.5 * (lo + hi);
      ((fMinusS(lo) < 0.0) == (fMinusS(mid) < 0.0) ? lo : hi) = mid;
    }
    const double B1 = 0.5 * (lo + hi);
    const PPair with = makePair(0.0, 1, 1, true), without = makePair(0.0, 1, 1, false);
    double worstM1 = 0.0;
    for (double E : {0.5, 1.0, 2.0, 3.0, 5.0}) {
      auto m1 = [&](const PPair &pair) {
        ECIntegral ec(&pair, lw, a);
        ec.SetRelativeTolerance(1e-6);
        return ec(1, 1, 0.0, 0.0, 1.0, 1.0, 1, 'M', sep + E, sep - B1, false);
      };
      worstM1 = std::fmax(worstM1, std::abs(m1(with)) / std::abs(m1(without)));
    }
    check("(c) M1 between two l = 1 states of the potential (B = " + std::to_string(B1) +
              " MeV): the whole amplitude is " + num(worstM1) + " of its exterior part (limit 1e-3)",
          std::fabs(B1 - 9.6065) < 1e-3 && worstM1 < 1e-3);
  }

  // (d) the interior charge
  {
    const double a = 5.35, B = 7.16, m1 = projectPair->GetM(1), m2 = projectPair->GetM(2);
    const double own = (projectPair->GetZ(1) * m2 - projectPair->GetZ(2) * m1) / (m1 + m2);
    double worstLinear = 0.0, worstOwn = 0.0;
    for (bool channel : {false, true})
      for (double E : {1.0, 2.4, -1.0}) {
        if (E < 0.0 && !channel) continue;
        const cplx outside = e1(makePair(0.0, 1, 0, false), lw, a, E, B, channel);
        const cplx t0 = e1(makePair(0.0, 1, 0, true), lw, a, E, B, channel);
        const cplx t1 = e1(makePair(0.0, 1, 0, true, 0.01), lw, a, E, B, channel);
        const cplx t2 = e1(makePair(0.0, 1, 0, true, 0.02), lw, a, E, B, channel);
        const cplx tOwn = e1(makePair(0.0, 1, 0, true, own), lw, a, E, B, channel);
        worstLinear = std::fmax(worstLinear, std::abs((t2 - t0) - 2.0 * (t1 - t0)) / std::abs(t1 - t0));
        worstOwn = std::fmax(worstOwn, std::abs((tOwn - t0) - (t0 - outside)) / std::abs(t0 - outside));
        if (std::abs(t1 - t0) == 0.0) worstLinear = 1.0;
      }
    for (bool channel : {false, true}) {  // (e): the unbound final state at E_f = 1 MeV
      const cplx t0 = e1(makePair(0.0, 1, 0, true), lw, a, 3.0, -1.0, channel);
      const cplx t1 = e1(makePair(0.0, 1, 0, true, 0.01), lw, a, 3.0, -1.0, channel);
      const cplx t2 = e1(makePair(0.0, 1, 0, true, 0.02), lw, a, 3.0, -1.0, channel);
      worstLinear = std::fmax(worstLinear, std::abs((t2 - t0) - 2.0 * (t1 - t0)) / std::abs(t1 - t0));
      if (std::abs(t1 - t0) == 0.0) worstLinear = 1.0;
    }
    check("(d) the amplitude is linear in the interior charge to " + num(worstLinear) + " (limit 1e-9)", worstLinear < 1e-9);
    check("(d) the pair's own charge number (" + num(own) + ") as interior charge doubles the interior part, to " +
              num(worstOwn) + " (limit 1e-9)", worstOwn < 1e-9);
  }

  std::printf(failures ? "interior_capture: %d check(s) FAILED\n" : "interior_capture: all checks passed\n", failures);
  return failures ? 1 : 0;
}
