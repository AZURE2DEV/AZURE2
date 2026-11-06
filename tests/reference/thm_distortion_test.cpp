/*!
 * Reference check of the THM distortion factor (ThmDistortion,
 * src/ThmDistortion.cpp): the zero-range prior-form DWBA amplitude
 *
 *   M = Int d^3r chi(-)*_{k_sF}(r) phi_sx(r) chi(+)_{k_aA}(beta r)
 *
 * and its plane-wave limit M_PW (Mukhamedzhanov & Pang, PRC 99 (2019) 064618,
 * eqs. 20-24; Mukhamedzhanov, arXiv:2609.04498, eqs. 22-30).
 *
 *  (c) |M|^2 and M_PW against thm_distortion_reference.py: mpmath's Coulomb
 *      functions and adaptive quadrature for point Coulomb (12C(14N,d) at 30
 *      MeV, forward and at 90 deg, Whittaker and Yukawa tails, a 3 fm cutoff;
 *      a made-up 18O(3He,d) at 115 MeV with the Trojan horse as the target,
 *      forward and at 60 deg), scipy's DOP853 for a complex Woods-Saxon
 *      potential in both channels.  Nothing is shared with the Numerov /
 *      COUL / Simpson numerics of the engine.
 *  (b) No distortion (both channels plane waves): the partial-wave sum equals
 *      the Fourier transform of the Yukawa tail, 4 pi/(kappa^2 + q^2), and R
 *      is 1.
 *  (e) An optical potential whose nuclear part is off (depths 0, point
 *      charge) gives the point-Coulomb amplitude, and a uniform sphere of
 *      0.02 fm nearly so.  (The a + A wave is needed at beta r, deep inside:
 *      0.3 fm already moves 12C(14N,d) by 2 %.)
 *  Also: the interpolated weight against R evaluated directly between the
 *  grid nodes, and R(E_ref) = 1.
 *
 * Run:  tests/reference/thm_distortion_test      (ctest: thm_distortion)
 *       tests/reference/thm_distortion_test -v   (print every value)
 */
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "Config.h"
#include "Constants.h"
#include "ThmDistortion.h"
#include "ThmExperiment.h"

// Each consumer of the engine defines this itself (see thm_coulomb_term_test).
Config *g_config = nullptr;

namespace {

bool verbose = false;
int failures = 0;

void Check(const std::string &what, double got, double want, double tol) {
  double rel = std::fabs(got - want) / std::max(std::fabs(want), 1.0e-300);
  bool ok = rel <= tol;
  if (!ok) failures++;
  if (!ok || verbose)
    std::printf("  %s  %-44s got %.12e want %.12e rel %.2e (tol %.0e)\n", ok ? "ok  " : "FAIL", what.c_str(), got,
                want, rel, tol);
}

const double kAmu = 931.49410242;  // EData::BuildThmGroups
const double kN14 = 13.9992355671, kC12 = 11.9967096429, kD = 2.0135532134, kHe3 = 3.0149322434,
             kO18 = 17.9947732059, kP = 1.0072764675;

ThmDistortion::Kinematics C12() {
  ThmDistortion::Kinematics k;
  k.Za = 7, k.ZA = 6, k.Zs = 1, k.Zx = 6;
  k.ma = kN14, k.mA = kC12, k.ms = kD, k.mx = kC12;
  k.horseIsBeam = true;
  k.mBeam = kN14, k.mTarget = kC12, k.beamEnergy = 30.0;
  k.bind = (k.mx + k.ms - k.ma) * kAmu;
  return k;
}
ThmDistortion::Kinematics He3() {
  ThmDistortion::Kinematics k;
  k.Za = 2, k.ZA = 8, k.Zs = 1, k.Zx = 1;
  k.ma = kHe3, k.mA = kO18, k.ms = kD, k.mx = kP;
  k.horseIsBeam = false;
  k.mBeam = kO18, k.mTarget = kHe3, k.beamEnergy = 115.0;
  k.bind = (k.mx + k.ms - k.ma) * kAmu;
  return k;
}

struct Reference {
  const char *name;
  double energy, m2, mpw;
};
// thm_distortion_reference.py
const Reference kReference[] = {
    {"c12_qf", 0.8, 4.397970624489e-10, 8.185198991197e+00},  // l <= 22
    {"c12_qf", 1.6, 1.205400874483e-10, 8.368412822221e+00},  // l <= 21
    {"c12_qf", 2.6, 2.010216793842e-12, 8.480831554576e+00},  // l <= 19
    {"c12_cm90", 1.2, 1.642320447168e-08, 1.095809008805e+01},  // l <= 19
    {"c12_rmin", 1.2, 1.884942422829e-10, 1.381260751734e+00},  // l <= 21
    {"optical", 1.0, 2.705453160181e-04, 8.235304497384e+00},  // l <= 20
    {"optical", 2.5, 1.150227976121e-05, 8.479947738456e+00},  // l <= 18
    {"he3_qf", 0.6, 3.425574843087e+02, 6.583432872958e+01},  // l <= 80
    {"he3_cm60", 0.6, 2.072506035193e+01, 4.301580372331e+00},  // l <= 78
};

ThmExperiment Settings(const std::string &name) {
  ThmExperiment x;
  x.name = name;
  x.distortion = ThmExperiment::DIST_COULOMB;
  if (name == "c12_cm90") {
    x.angleKind = 2;
    x.angle = 90.0;
    x.boundYukawa = true;
  } else if (name == "c12_rmin") {
    x.boundRmin = 3.0;
  } else if (name == "he3_cm60") {
    x.angleKind = 2;
    x.angle = 60.0;
  } else if (name == "optical") {
    x.distortion = ThmExperiment::DIST_OPTICAL;
    const double aa[10] = {50.0, 5.5, 0.6, 20.0, 5.6, 0.6, 0.0, 0.0, 0.0, 5.5};
    const double sf[10] = {90.0, 3.4, 0.75, 0.0, 0.0, 0.0, 10.0, 3.8, 0.65, 3.8};
    x.opticalAA.kind = x.opticalSF.kind = 2;
    for (int t = 0; t < 10; t++) x.opticalAA.p[t] = aa[t], x.opticalSF.p[t] = sf[t];
  }
  return x;
}

// Builds on a small grid around E (E_ref = E).
std::string Build(ThmDistortion &d, const ThmExperiment &x, const ThmDistortion::Kinematics &k, double e) {
  return d.Build(x, k, e - 0.02, e + 0.02, e);
}

}  // namespace

int main(int argc, char **argv) {
  verbose = argc > 1 && std::strcmp(argv[1], "-v") == 0;

  std::printf("(c) |M|^2 and M_PW against the independent reference\n");
  for (const Reference &r : kReference) {
    std::string name = r.name;
    ThmDistortion d;
    std::string why = Build(d, Settings(name), name.compare(0, 3, "he3") == 0 ? He3() : C12(), r.energy);
    if (!why.empty()) {
      std::printf("  FAIL  %s: %s\n", name.c_str(), why.c_str());
      failures++;
      continue;
    }
    ThmDistortion::Point p = d.Evaluate(r.energy);
    char label[96];
    std::snprintf(label, sizeof label, "%s E=%.2f |M|^2", r.name, r.energy);
    Check(label, std::norm(p.m), r.m2, 2.0e-5);
    std::snprintf(label, sizeof label, "%s E=%.2f M_PW", r.name, r.energy);
    Check(label, p.mpw, r.mpw, 2.0e-6);
  }

  std::printf("(b) plane waves in both channels: the Yukawa Fourier transform, R = 1\n");
  for (double angle : {0.0, 70.0, 180.0}) {
    ThmExperiment x;
    x.distortion = ThmExperiment::DIST_OPTICAL;
    x.opticalAA.kind = x.opticalSF.kind = 0;
    x.boundYukawa = true;
    x.angleKind = 2;
    x.angle = angle;
    for (int sys = 0; sys < 2; sys++) {
      ThmDistortion d;
      double e = sys ? 0.6 : 1.5;
      std::string why = Build(d, x, sys ? He3() : C12(), e);
      if (!why.empty()) {
        std::printf("  FAIL  plane: %s\n", why.c_str());
        failures++;
        continue;
      }
      for (double de : {-0.013, 0.0, 0.017}) {
        ThmDistortion::Point p = d.Evaluate(e + de);
        double exact = 4.0 * M_PI / (d.kappa * d.kappa + p.q * p.q);
        char label[96];
        std::snprintf(label, sizeof label, "%s %g deg E=%.3f M", sys ? "he3" : "c12", angle, e + de);
        Check(label, std::sqrt(std::norm(p.m)), exact, 2.0e-6);
        std::snprintf(label, sizeof label, "%s %g deg E=%.3f M_PW", sys ? "he3" : "c12", angle, e + de);
        Check(label, p.mpw, exact, 2.0e-6);
        std::snprintf(label, sizeof label, "%s %g deg E=%.3f R", sys ? "he3" : "c12", angle, e + de);
        Check(label, d.R(p), 1.0, 2.0e-6);
      }
    }
  }

  std::printf("(e) optical potential with the nuclear part off = point Coulomb\n");
  {
    ThmExperiment coulomb = Settings("c12_qf"), off = coulomb, sphere = coulomb;
    off.distortion = sphere.distortion = ThmExperiment::DIST_OPTICAL;
    off.opticalAA.kind = off.opticalSF.kind = 2;  // all ten numbers 0
    sphere.opticalAA.kind = sphere.opticalSF.kind = 2;
    sphere.opticalAA.p[9] = sphere.opticalSF.p[9] = 0.02;  // uniform sphere, 0.02 fm
    for (double e : {0.9, 2.4}) {
      ThmDistortion dc, doff, ds;
      std::string why = Build(dc, coulomb, C12(), e) + Build(doff, off, C12(), e) + Build(ds, sphere, C12(), e);
      if (!why.empty()) {
        std::printf("  FAIL  %s\n", why.c_str());
        failures++;
        continue;
      }
      ThmDistortion::Point pc = dc.Evaluate(e), po = doff.Evaluate(e), ps = ds.Evaluate(e);
      char label[96];
      std::snprintf(label, sizeof label, "E=%.1f |M|^2 nuclear off", e);
      Check(label, std::norm(po.m), std::norm(pc.m), 1.0e-12);
      std::snprintf(label, sizeof label, "E=%.1f |M|^2 RC = 0.02 fm", e);
      Check(label, std::norm(ps.m), std::norm(pc.m), 1.0e-3);
    }
  }

  std::printf("the weight: interpolation between the nodes, R(E_ref) = 1\n");
  {
    ThmDistortion d;
    ThmExperiment x = Settings("c12_qf");
    std::string why = d.Build(x, C12(), 0.7, 2.7, 1.7);
    if (!why.empty()) {
      std::printf("  FAIL  %s\n", why.c_str());
      failures++;
    } else {
      Check("R(E_ref)", d.R(d.Evaluate(1.7)), 1.0, 1.0e-14);
      Check("weight at E_ref", d.Weight(1.7), 1.0, 1.0e-6);
      for (double e : {0.7137, 1.2551, 2.1049, 2.6953}) {
        char label[96];
        std::snprintf(label, sizeof label, "weight(%.4f) vs R", e);
        Check(label, d.Weight(e), d.R(d.Evaluate(e)), 1.0e-5);
      }
    }
  }

  if (failures) {
    std::printf("%d check(s) failed\n", failures);
    return 1;
  }
  std::printf("all checks passed\n");
  return 0;
}
