/*!
 * Reference check of the global optical potentials of the library
 * (physics/transfer/GlobalOptical.h) in the THM distorted waves
 * (ThmDistortion: opticalAA=<name>, opticalSF=<name>).
 *
 *  The ten numbers of every set against RIPL-3 / FRONT21 are the library's test
 *  (physics/tests/transfer_global_optical_test.cpp); here their use in the engine:
 *  (b) Elastic scattering d + 40Ca at 56 MeV with daehnick80 (no spin-orbit):
 *      sigma/sigma_R from the S matrix of the waves the distortion factor
 *      solves (the library's SolveStoredWaves) against scipy's DOP853 + mpmath, and against the
 *      measured ratio (Hatanaka et al. 1980, EXFOR E0682-022): within 0.2 dex
 *      at every angle (15-78 deg), 0.08 dex rms.
 *  (c) In the distortion factor: a global s + F potential is re-evaluated at
 *      E_sF of every energy (the amplitude equals that of the ten numbers
 *      at E_sF(E) written out); refusals (wrong projectile, outside the
 *      validity range) and :extrapolate (warned); the parser.
 *
 * Run:  tests/reference/thm_optical_test      (ctest: thm_optical)
 *       tests/reference/thm_optical_test -v   (print every value)
 */
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include <gsl/gsl_sf_gamma.h>

#include "Config.h"
#include "Constants.h"
#include "ThmDistortion.h"
#include "ThmExperiment.h"
#include "azure2/physics/transfer/GlobalOptical.h"

namespace {

bool verbose = false;
int failures = 0;

void Check(const std::string &what, double got, double want, double tol, double absTol = 0.0) {
  double diff = std::fabs(got - want);
  double rel = diff / std::max(std::fabs(want), 1.0e-300);
  bool ok = rel <= tol || diff <= absTol;
  if (!ok) failures++;
  if (!ok || verbose)
    std::printf("  %s  %-44s got %.12e want %.12e rel %.2e (tol %.0e)\n", ok ? "ok  " : "FAIL", what.c_str(), got,
                want, rel, tol);
}
void Expect(const std::string &what, bool ok, const std::string &detail = "") {
  if (!ok) failures++;
  if (!ok || verbose) std::printf("  %s  %s%s%s\n", ok ? "ok  " : "FAIL", what.c_str(), detail.empty() ? "" : ": ",
                                  detail.c_str());
}

// thm_optical_reference.py: theta_cm, sigma/sigma_R (DOP853), measured sigma/sigma_R.
struct Elastic {
  double theta, model, data;
};
const Elastic kElastic[] = {
    {14.72, 2.53602027e-01, 2.45116e-01}, {18.92, 3.70529879e-01, 4.70149e-01},
    {21.01, 8.25919761e-01, 9.27574e-01}, {23.11, 8.81823026e-01, 9.28617e-01},
    {31.48, 6.60691920e-01, 7.70816e-01}, {34.10, 1.56693212e+00, 1.64812e+00},
    {36.70, 1.92378399e+00, 1.88892e+00}, {44.51, 1.29716392e+00, 1.45748e+00},
    {49.69, 2.42726378e+00, 2.47953e+00}, {57.43, 1.71132250e+00, 1.79389e+00},
    {67.69, 1.30000623e+00, 1.50358e+00}, {77.87, 5.17349669e-01, 7.48311e-01},
};

namespace transfer = azure2::physics::transfer;

// The ten numbers V,R,a,W,RW,aW,WD,RD,aD,RC (MeV, fm) of set `name` for (zp, ap) on (zt, at) at elab (MeV).
void Evaluate(const char *name, int zp, int ap, int zt, int at, double elab, double p[10]) {
  const transfer::OpticalPotential o =
      transfer::GlobalOpticalSets()[transfer::GlobalOpticalIndex(name)].at(zp, ap, zt, at, elab);
  const transfer::WoodsSaxonTerm *terms[3] = {&o.real_volume, &o.imaginary_volume, &o.imaginary_surface};
  for (int t = 0; t < 3; t++) {
    p[3 * t] = terms[t]->depth_MeV;
    p[3 * t + 1] = transfer::RadiusOf(o, terms[t]->r0_fm);
    p[3 * t + 2] = terms[t]->diffuseness_fm;
  }
  p[9] = transfer::RadiusOf(o, o.coulomb_r0_fm);
}

const double kAmu = 931.49410242;  // WorkingCopy::BuildThmGroups
const double kD = 2.0135532134, kN = 1.0086649159, kP = 1.0072764675, kF19 = 18.9934689019;

// 19F + d at 55 MeV, neutron spectator (examples/f19_pag_thm): a = d, A = 19F, s = n, x = p.
ThmDistortion::Kinematics F19() {
  ThmDistortion::Kinematics k;
  k.Za = 1, k.ZA = 9, k.Zs = 0, k.Zx = 1;
  k.Aa = 2, k.AA = 19, k.As = 1, k.Ax = 1;
  k.ma = kD, k.mA = kF19, k.ms = kN, k.mx = kP;
  k.horseIsBeam = false;
  k.mBeam = kF19, k.mTarget = kD, k.beamEnergy = 55.0;
  k.bind = (k.mx + k.ms - k.ma) * kAmu;
  return k;
}

ThmExperiment Optical(const std::string &aa, const std::string &sf, std::string *why) {
  std::vector<ThmExperiment> xs;
  *why = ParseThmExperimentLine("experiment[X] beam=19F target=d spectator=n Ebeam=55 "
                                "distortion=optical opticalAA=" + aa + " opticalSF=" + sf,
                                xs);
  return xs.empty() ? ThmExperiment() : xs.front();
}

std::string Build(ThmDistortion &d, const ThmExperiment &x, double e, double lo, double hi) {
  d.dataLo = lo;
  d.dataHi = hi;
  return d.Build(x, F19(), e - 0.02, e + 0.02, e);
}

std::string Ten(const double p[10]) {
  std::string s;
  char b[32];
  for (int t = 0; t < 10; t++) {
    std::snprintf(b, sizeof b, "%.17g", p[t]);
    s += (t ? "," : "") + std::string(b);
  }
  return s;
}

}  // namespace

int main(int argc, char **argv) {
  verbose = argc > 1 && std::strcmp(argv[1], "-v") == 0;

  std::printf("(b) elastic d + 40Ca at 56 MeV, daehnick80: S matrix and sigma/sigma_R\n");
  {
    const double mCa = 39.9516979, elab = 56.0;
    ThmDistortion d;
    ThmDistortion::Channel c;
    c.kind = ThmDistortion::Channel::WOODS_SAXON;
    Evaluate("daehnick80", 1, 2, 20, 40, elab, c.p);
    c.Z1 = 1, c.Z2 = 20;
    c.muU = kD * mCa / (kD + mCa);
    c.mu = c.muU * uconv;
    const double ecm = elab * mCa / (kD + mCa);
    c.k = std::sqrt(2.0 * c.mu * ecm) / hbarc;
    c.eta = 20.0 * fstruc * c.mu / (hbarc * c.k);
    const int lmax = 45;
    std::vector<complex> S(lmax + 1);
    bool ok = true;
    try {
      // The waves as the distortion factor solves them (matched beyond the turning point).
      const transfer::DistortedWaves w =
          transfer::SolveStoredWaves(d.LibraryChannel(c, ecm), 0.02, 2, lmax, d.numerics.matching);
      for (int l = 0; l <= lmax; l++) S[l] = w.waves[l].S;
    } catch (const std::exception &) {
      ok = false;
    }
    Expect("waves normalized", ok);
    gsl_sf_result lnr, arg;
    gsl_sf_lngamma_complex_e(1.0, c.eta, &lnr, &arg);
    double sum2 = 0.0;
    int n = 0;
    for (const Elastic &e : kElastic) {
      const double t = e.theta * M_PI / 180.0, s2 = std::pow(std::sin(t / 2.0), 2), x = std::cos(t);
      const complex fc = -c.eta / (2.0 * c.k * s2) * std::polar(1.0, -c.eta * std::log(s2) + 2.0 * arg.val);
      complex fn(0.0, 0.0);
      double sigma = arg.val, pm = 1.0, pl = x;
      for (int l = 0; l <= lmax; l++) {
        if (l > 0) sigma += std::atan(c.eta / l);
        double P = l == 0 ? 1.0 : x;
        if (l >= 2) {
          double next = ((2.0 * l - 1.0) * x * pl - (l - 1.0) * pm) / l;
          pm = pl;
          pl = next;
          P = pl;
        }
        fn += (2.0 * l + 1.0) * std::polar(1.0, 2.0 * sigma) * (S[l] - 1.0) * P / (2.0 * complex(0.0, 1.0) * c.k);
      }
      const double ratio = std::norm(fc + fn) / std::norm(fc);
      char label[96];
      std::snprintf(label, sizeof label, "sigma/sigma_R %.2f deg vs DOP853", e.theta);
      Check(label, ratio, e.model, 2.0e-4);
      const double dex = std::log10(ratio / e.data);
      std::snprintf(label, sizeof label, "sigma/sigma_R %.2f deg vs data (|dex| <= 0.2)", e.theta);
      Expect(label, std::fabs(dex) <= 0.2, verbose ? std::to_string(dex) : "");
      sum2 += dex * dex;
      n++;
    }
    const double rms = std::sqrt(sum2 / n);
    Expect("rms log10(model/data) <= 0.08", rms <= 0.08, std::to_string(rms));
  }

  std::printf("(c) global potentials in the distortion factor (19F + d, 55 MeV)\n");
  {
    std::string why;
    // The s + F potential follows E_sF: the global amplitude equals the one of
    // the ten numbers written out at E_sF(E) (and a + A at E_aA).
    ThmExperiment g = Optical("ancai06", "kd03:extrapolate", &why);
    Expect("parse ancai06 / kd03:extrapolate", why.empty(), why);
    for (double e : {0.05, 0.85}) {
      ThmDistortion dg;
      why = Build(dg, g, e, 0.05, 0.85);
      Expect("build global", why.empty(), why);
      if (!why.empty()) continue;
      Expect("one warning (kd03 for 20Ne)", dg.warnings.size() == 1,
             dg.warnings.empty() ? "" : dg.warnings.front());
      double pa[10], ps[10];
      Evaluate("ancai06", 1, 2, 9, 19, dg.aa.LabEnergy(dg.eAA), pa);
      const double esf = dg.EsF(e), mF = kP + kF19;
      Evaluate("kd03", 0, 1, 10, 20, esf * (kN + mF) / mF, ps);
      ThmExperiment w = Optical(Ten(pa), Ten(ps), &why);
      ThmDistortion dw;
      why += Build(dw, w, e, 0.05, 0.85);
      Expect("build written-out", why.empty(), why);
      if (!why.empty()) continue;
      ThmDistortion::Point pg = dg.Evaluate(e), pw = dw.Evaluate(e);
      char label[96];
      std::snprintf(label, sizeof label, "E=%.2f |M|^2 global = written out", e);
      Check(label, std::norm(pg.m), std::norm(pw.m), 1.0e-6);
      // ... and not the numbers of the other end: the potential does move.
      Evaluate("kd03", 0, 1, 10, 20, dg.sf.LabEnergy(dg.EsF(0.9 - e)), ps);
      Expect("the s + F depths move across the data", std::fabs(ps[6] - dg.SfAt(e).p[6]) > 1.0e-3);
    }
    // Refusals and the range.
    ThmDistortion d1;
    why = Build(d1, Optical("kd03", "plane", &why), 0.5, 0.05, 0.85);
    Expect("kd03 for d + 19F refused", why.find("does not describe") != std::string::npos, why);
    ThmDistortion d2;
    why = Build(d2, Optical("daehnick80", "plane", &why), 0.5, 0.05, 0.85);
    Expect("daehnick80 for d + 19F at 5.8 MeV refused",
           why.find("outside its validity range") != std::string::npos &&
               why.find("daehnick80:extrapolate") != std::string::npos,
           why);
    ThmDistortion d3;
    why = Build(d3, Optical("ancai06", "kd03", &why), 0.5, 0.05, 0.85);
    Expect("kd03 for n + 20Ne refused (A = 20 < 24)", why.find("target A = 20") != std::string::npos, why);
    ThmDistortion d4;
    why = Build(d4, Optical("ancai06", "plane", &why), 0.5, 0.05, 0.85);
    Expect("ancai06 for d + 19F at 5.8 MeV accepted, no warning", why.empty() && d4.warnings.empty(), why);
    // Parser.
    std::vector<ThmExperiment> xs;
    why = ParseThmExperimentLine("experiment[X] distortion=optical opticalSF=kd04", xs);
    Expect("unknown name refused", why.find("kd03") != std::string::npos, why);
    why = ParseThmExperimentLine("experiment[X] distortion=optical opticalSF=kd03:extra", xs);
    Expect("unknown option refused", !why.empty(), why);
    why = ParseThmExperimentLine("experiment[X] distortion=optical opticalSF=1,2,3", xs);
    Expect("ten-number message kept", why.find("expected plane, coulomb or ten numbers") != std::string::npos, why);
  }

  if (failures) {
    std::printf("thm_optical: %d check(s) FAILED\n", failures);
    return 1;
  }
  std::printf("thm_optical: all checks passed\n");
  return 0;
}
