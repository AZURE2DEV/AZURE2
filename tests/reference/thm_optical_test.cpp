/*!
 * Reference check of the global optical potentials (ThmOptical,
 * src/ThmOptical.cpp) and of their use in the THM distorted waves
 * (ThmDistortion: opticalAA=<name>, opticalSF=<name>).
 *
 *  (a) The ten numbers V,R,a,W,RW,aW,WD,RD,aD,RC of every model at 23
 *      (projectile, target, lab energy) points against
 *      thm_optical_reference.py, which evaluates the RIPL-3 library entries
 *      (om-parameter-u.dat: An & Cai 6200, Koning & Delaroche 2405/5405,
 *      Becchetti & Greenlees 7100/8100, McFadden & Satchler 9100, Avrigeanu
 *      9600) by the library's own coefficient formulas, and the FRONT21
 *      (TWOFNR) code for Daehnick and Liang.  Includes the published d + 19F
 *      An & Cai numbers used before for examples/f19_pag_thm.
 *  (b) Elastic scattering d + 40Ca at 56 MeV with daehnick80 (no spin-orbit):
 *      sigma/sigma_R from the engine's S matrix (ThmDistortion::Wave,
 *      Numerov + COUL) against scipy's DOP853 + mpmath, and against the
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
#include "ThmOptical.h"

Config *g_config = nullptr;

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

struct Params {
  const char *model;
  int zp, ap, zt, at;
  double elab;
  double p[10];
};
// thm_optical_reference.py (RIPL-3 / FRONT21).
const Params kParams[] = {
    {"ancai06", 1, 2, 9, 19, 5.83, {92.5676139, 3.066238699, 0.7526218608, 1.466626, 3.581964152, 0.5881598351, 10.651602, 3.711647799, 0.6964409022, 3.476927348}},  // RIPL 6200
    {"ancai06", 1, 2, 12, 24, 2, {94.02329163, 3.31518301, 0.7553446892, 1.2284, 3.863971379, 0.5665500859, 10.7688, 3.999921854, 0.7098389467, 3.75850238}},  // RIPL 6200
    {"ancai06", 1, 2, 20, 40, 56, {82.02421478, 3.932024581, 0.7620913939, 4.5872, 4.562737221, 0.5130048107, 9.1164, 4.714215826, 0.7430370174, 4.456197317}},  // RIPL 6200
    {"ancai06", 1, 2, 82, 208, 100, {76.99507508, 6.817830942, 0.7936549009, 7.324, 7.831814739, 0.2625007863, 7.77, 8.055939511, 0.8983495125, 7.720264754}},  // RIPL 6200
    {"daehnick80", 1, 2, 20, 40, 56, {79.08627122, 4.001343715, 0.8042, 3.676040966, 4.531436259, 0.729396628, 9.979959034, 4.531436259, 0.729396628, 4.445937461}},  // FRONT21
    {"daehnick80", 1, 2, 28, 58, 20, {89.66548314, 4.52892567, 0.743, 0.498758334, 5.128911549, 0.7862461872, 12.22124167, 5.128911549, 0.7862461872, 5.032139633}},  // FRONT21
    {"daehnick80", 1, 2, 82, 208, 80, {79.87891912, 6.9322408, 0.845, 6.750264185, 7.850614581, 0.9047494496, 7.529735815, 7.850614581, 0.9047494496, 7.702489778}},  // FRONT21
    {"daehnick80", 1, 2, 9, 19, 5.8, {89.96006892, 3.122029929, 0.71886, 0.04147828561, 3.535632185, 0.7020729378, 12.30932171, 3.535632185, 0.7020729378, 3.468922143}},  // FRONT21
    {"kd03", 0, 1, 10, 20, 2.6, {53.37047263, 3.13392913, 0.674826, 0.3835113907, 3.13392913, 0.674826, 6.86586386, 3.52705041, 0.541288, 0}},  // RIPL 2405
    {"kd03", 0, 1, 26, 56, 14, {47.45208266, 4.583141938, 0.6694728, 1.118259206, 4.583141938, 0.6694728, 7.168952993, 4.903837657, 0.5353264, 0}},  // RIPL 2405
    {"kd03", 0, 1, 82, 208, 0.5, {47.7303101, 7.320197247, 0.6468704, 0.07572086866, 7.320197247, 0.6468704, 2.586696471, 7.397286765, 0.5101552, 0}},  // RIPL 2405
    {"kd03", 0, 1, 82, 208, 150, {13.88561425, 7.320197247, 0.6468704, 11.74180123, 7.320197247, 0.6468704, 0.7592637623, 7.397286765, 0.5101552, 0}},  // RIPL 2405
    {"kd03", 1, 1, 20, 40, 30, {46.5397075, 4.053875274, 0.671852, 2.946652563, 4.053875274, 0.671852, 6.417163507, 4.405560697, 0.53952, 4.395893121}},  // RIPL 5405
    {"kd03", 1, 1, 82, 208, 65, {37.50515578, 7.320197247, 0.6468704, 6.33616321, 7.320197247, 0.6468704, 5.302106915, 7.397286765, 0.626964, 7.226321535}},  // RIPL 5405
    {"kd03", 1, 1, 12, 24, 3, {56.31526791, 3.355698429, 0.6742312, 0.3142593949, 3.355698429, 0.6742312, 6.048115684, 3.740274332, 0.531192, 3.88496505}},  // RIPL 5405
    {"bg71", 1, 3, 40, 90, 15, {161.7388889, 5.377685696, 0.72, 28.82777778, 6.273966645, 0.84, 0, 0, 0, 5.825826171}},  // RIPL 7100
    {"bg71", 2, 3, 28, 58, 30, {148.5241379, 4.645051969, 0.72, 33.31724138, 5.419227297, 0.88, 0, 0, 0, 5.032139633}},  // RIPL 8100
    {"liang09", 2, 3, 8, 16, 10, {118.0948739, 2.977479936, 0.740955184, 0, 0, 0, 21.53873315, 3.032617646, 0.794753626, 3.248076467}},  // FRONT21
    {"liang09", 2, 3, 50, 120, 100, {108.3687635, 5.78982683, 0.8145389365, 17.4533, 6.958045141, 0.8812200889, 8.919704402, 5.87994698, 0.9532602666, 6.357894728}},  // FRONT21
    {"mcfadden66", 2, 4, 8, 16, 20, {185, 3.52777894, 0.52, 25, 3.52777894, 0.52, 0, 0, 0, 3.52777894}},  // RIPL 9100
    {"mcfadden66", 2, 4, 82, 208, 24.7, {185, 8.294988992, 0.52, 25, 8.294988992, 0.52, 0, 0, 0, 8.294988992}},  // RIPL 9100
    {"avrigeanu94", 2, 4, 20, 40, 20, {131.5323107, 4.257840107, 0.7879304089, 10.80556207, 5.369324473, 0.6236009621, 0, 0, 0, 4.257840107}},  // RIPL 9600
    {"avrigeanu94", 2, 4, 82, 208, 50, {172.4577483, 7.37661521, 0.7666375668, 12.53196341, 9.302237655, 0.5735001573, 0, 0, 0, 7.37661521}},  // RIPL 9600
};

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

const double kAmu = 931.49410242;  // EData::BuildThmGroups
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
  *why = ParseThmExperimentLine("experiment[X] segments=1 beam=19F target=d spectator=n Ebeam=55 "
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

  std::printf("(a) parameters against RIPL-3 / FRONT21\n");
  for (const Params &r : kParams) {
    int m = ThmGlobalOpticalIndex(r.model);
    Expect(std::string(r.model) + " known", m >= 0);
    if (m < 0) continue;
    Expect(std::string(r.model) + " is for its projectile", ThmGlobalOpticalFor(m, r.zp, r.ap));
    double p[10];
    ThmGlobalOpticalEvaluate(m, r.zp, r.ap, r.zt, r.at, r.elab, p);
    static const char *names[10] = {"V", "R", "a", "W", "RW", "aW", "WD", "RD", "aD", "RC"};
    for (int t = 0; t < 10; t++) {
      char label[96];
      std::snprintf(label, sizeof label, "%s (%d,%d)+(%d,%d) %g MeV %s", r.model, r.zp, r.ap, r.zt, r.at, r.elab,
                    names[t]);
      Check(label, p[t], r.p[t], 1.0e-9, 1.0e-12);
    }
  }
  Expect("kd03 is not a deuteron potential", !ThmGlobalOpticalFor(ThmGlobalOpticalIndex("kd03"), 1, 2));
  Expect("ancai06 is not an alpha potential", !ThmGlobalOpticalFor(ThmGlobalOpticalIndex("ancai06"), 2, 4));

  std::printf("(b) elastic d + 40Ca at 56 MeV, daehnick80: S matrix and sigma/sigma_R\n");
  {
    const double mCa = 39.9516979, elab = 56.0;
    ThmDistortion d;
    ThmDistortion::Channel c;
    c.kind = ThmDistortion::Channel::WOODS_SAXON;
    ThmGlobalOpticalEvaluate(ThmGlobalOpticalIndex("daehnick80"), 1, 2, 20, 40, elab, c.p);
    c.Z1 = 1, c.Z2 = 20;
    c.mu = kD * mCa / (kD + mCa) * uconv;
    const double ecm = elab * mCa / (kD + mCa);
    c.k = std::sqrt(2.0 * c.mu * ecm) / hbarc;
    c.eta = 20.0 * fstruc * c.mu / (hbarc * c.k);
    const int lmax = 45;
    std::vector<complex> S(lmax + 1);
    std::vector<complex> u;
    bool ok = true;
    for (int l = 0; l <= lmax && ok; l++) {
      complex T;
      ok = d.Wave(c, l, 0.02, 2, u, &T);
      S[l] = 1.0 + 2.0 * complex(0.0, 1.0) * T;
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
      ThmGlobalOpticalEvaluate(ThmGlobalOpticalIndex("ancai06"), 1, 2, 9, 19, dg.aa.LabEnergy(dg.eAA), pa);
      const double esf = dg.EsF(e), mF = kP + kF19;
      ThmGlobalOpticalEvaluate(ThmGlobalOpticalIndex("kd03"), 0, 1, 10, 20, esf * (kN + mF) / mF, ps);
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
      ThmGlobalOpticalEvaluate(ThmGlobalOpticalIndex("kd03"), 0, 1, 10, 20, dg.sf.LabEnergy(dg.EsF(0.9 - e)), ps);
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
    why = ParseThmExperimentLine("experiment[X] segments=1 distortion=optical opticalSF=kd04", xs);
    Expect("unknown name refused", why.find("kd03") != std::string::npos, why);
    why = ParseThmExperimentLine("experiment[X] segments=1 distortion=optical opticalSF=kd03:extra", xs);
    Expect("unknown option refused", !why.empty(), why);
    why = ParseThmExperimentLine("experiment[X] segments=1 distortion=optical opticalSF=1,2,3", xs);
    Expect("ten-number message kept", why.find("expected plane, coulomb or ten numbers") != std::string::npos, why);
  }

  if (failures) {
    std::printf("thm_optical: %d check(s) FAILED\n", failures);
    return 1;
  }
  std::printf("thm_optical: all checks passed\n");
  return 0;
}
