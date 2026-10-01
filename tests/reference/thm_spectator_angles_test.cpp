/*!
 * The spectator-direction window of a THM experiment
 * (ThmDistortion::AngleNodes, `spectatorAngles=`;
 * docs/source/theory/thm_implementation.rst, "Experimental acceptance"):
 * the accepted c.m. directions at one energy, as Gauss-Legendre nodes in
 * cos(theta_cm) with the acceptance as weight.
 *
 *  (a) The measure: the weights of a uniform window sum to the solid angle
 *      (over 2 pi) of the accepted c.m. directions, against a brute-force
 *      midpoint sum of the indicator theta_lab(theta_cm) in the window over
 *      4e6 bins of cos(theta_cm) -- c.m. windows, lab windows with one
 *      branch (gamma = V_cm/v_s < 1) and two (gamma > 1), the whole sphere
 *      (2), with the q cut of a ps window.
 *  (b) Every node lies inside the window (its lab angle, from
 *      tan theta_lab = sin theta/(cos theta + gamma), and its q).
 *  (c) A lab window of zero width on two branches: the branch weights are
 *      the limit of a shrinking window (1e-6 deg around it; [0, 0.01] at
 *      0 deg, where they are (gamma - 1)^2 : (gamma + 1)^2).
 *  (d) An acceptance table: a linear ramp in the c.m. or the lab angle
 *      integrates to the brute-force sum of the ramp.
 *
 * Kinematics: 12C(14N,d) at 30 MeV (Trojan horse = beam; the deuteron is
 * slower in the c.m. than the c.m. above E = 2.2 MeV) and the made-up
 * 18O(3He,d) at 115 MeV (Trojan horse = target).
 *
 * Run:  tests/reference/thm_spectator_angles_test      (ctest: thm_spectator_angles)
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

Config *g_config = nullptr;

namespace {

bool verbose = false;
int failures = 0;

void Check(const std::string &what, double got, double want, double tol, bool absolute = false) {
  double err = std::fabs(got - want) / (absolute ? 1.0 : std::max(std::fabs(want), 1.0e-300));
  bool ok = err <= tol;
  if (!ok) failures++;
  if (!ok || verbose)
    std::printf("  %s  %-58s got %.12e want %.12e %s %.2e (tol %.0e)\n", ok ? "ok  " : "FAIL", what.c_str(), got,
                want, absolute ? "abs" : "rel", err, tol);
}
void Truth(const std::string &what, bool ok) {
  if (!ok) failures++;
  if (!ok || verbose) std::printf("  %s  %s\n", ok ? "ok  " : "FAIL", what.c_str());
}

const double kAmu = 931.49410242;
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

struct Window {
  bool cm = false;
  double lo = 0.0, hi = 180.0;
  bool cut = false;
  double qlo = 0.0, qhi = 0.0;
  std::vector<double> t, w;  // table
};

// A ThmDistortion set up for a window (Setup; the s + F channel, vcm, beta,
// k_aA as the engine has them), angSlots for energies up to eHi.
bool Make(ThmDistortion &d, const ThmDistortion::Kinematics &k, const Window &win, int nodes, double eLo,
          double eHi) {
  ThmExperiment x;
  x.name = "test";
  x.distortion = ThmExperiment::DIST_COULOMB;
  x.angleWindow = win.t.empty() ? 1 : 2;
  x.angleCm = win.cm;
  x.angleMin = win.t.empty() ? win.lo : win.t.front();
  x.angleMax = win.t.empty() ? win.hi : win.t.back();
  x.angleTableT = win.t;
  x.angleTableW = win.w;
  x.angleNodes = nodes;
  if (win.cut) {
    x.psKind = ThmExperiment::PS_HULTHEN;
    x.psMin = win.qlo;
    x.psMax = win.qhi;
  }
  d.sf.mu = k.ms * (k.mx + k.mA) / (k.ms + k.mx + k.mA) * uconv;
  std::string why = d.Setup(x, k, eLo);
  if (!why.empty()) {
    std::printf("  FAIL  Setup: %s\n", why.c_str());
    failures++;
    return false;
  }
  d.SetAngleSlots(eHi);
  return true;
}

double Ks(const ThmDistortion &d, double e) { return std::sqrt(2.0 * d.sf.mu * d.EsF(e)) / hbarc; }
double Gamma(const ThmDistortion &d, double ks) { return d.vcm / (hbarc * ks / (d.kin.ms * uconv)); }

// Brute force: Int d cos(theta_cm) A(theta) [theta in the window, q in the cut].
double Brute(const ThmDistortion &d, const Window &win, double ks) {
  const int n = 4000000;
  const double g = Gamma(d, ks), kb = d.beta * d.aa.k;
  double sum = 0.0;
  for (int i = 0; i < n; i++) {
    const double c = -1.0 + (i + 0.5) * 2.0 / n;
    const double th = std::acos(c);
    const double tl = std::atan2(std::sin(th), c + g) * 180.0 / M_PI, tc = th * 180.0 / M_PI;
    const double t = win.cm ? tc : tl;
    double a = 0.0;
    if (win.t.empty()) {
      a = (t >= win.lo && t <= win.hi) ? 1.0 : 0.0;
    } else if (t >= win.t.front() && t <= win.t.back()) {
      size_t j = 1;
      while (j + 1 < win.t.size() && win.t[j] < t) j++;
      a = win.w[j - 1] + (t - win.t[j - 1]) / (win.t[j] - win.t[j - 1]) * (win.w[j] - win.w[j - 1]);
    }
    if (a == 0.0) continue;
    if (win.cut) {
      const double x = (d.kin.horseIsBeam ? 1.0 : -1.0) * c;
      const double q = std::sqrt(std::max(0.0, ks * ks + kb * kb - 2.0 * ks * kb * x)) * hbarc;
      if (q < win.qlo || q > win.qhi) continue;
    }
    sum += a;
  }
  return sum * 2.0 / n;
}

void Measure(const std::string &tag, const ThmDistortion::Kinematics &k, const Window &win, double e,
             double tol = 5.0e-6) {
  ThmDistortion d;
  if (!Make(d, k, win, 24, e - 0.1, e + 0.1)) return;
  const double ks = Ks(d, e), g = Gamma(d, ks), kb = d.beta * d.aa.k;
  std::vector<ThmDistortion::AngleNode> nodes;
  bool any = d.AngleNodes(ks, nodes);
  double sum = 0.0;
  bool inside = true;
  for (const ThmDistortion::AngleNode &n : nodes) {
    if (!(n.w > 0.0)) continue;
    sum += n.w;
    const double th = n.theta * M_PI / 180.0;
    const double tl = std::atan2(std::sin(th), std::cos(th) + g) * 180.0 / M_PI;
    const double t = win.cm ? n.theta : tl;
    const double lo = win.t.empty() ? win.lo : win.t.front(), hi = win.t.empty() ? win.hi : win.t.back();
    inside = inside && t >= lo - 1.0e-9 && t <= hi + 1.0e-9;
    const double x = (d.kin.horseIsBeam ? 1.0 : -1.0) * std::cos(th);
    const double q = std::sqrt(std::max(0.0, ks * ks + kb * kb - 2.0 * ks * kb * x));
    inside = inside && std::fabs(q - n.q) < 1.0e-12 && std::fabs(x - n.x) < 1.0e-12;
    if (win.cut) inside = inside && q * hbarc >= win.qlo - 1.0e-9 && q * hbarc <= win.qhi + 1.0e-9;
  }
  char head[160];
  std::snprintf(head, sizeof head, "%s, E = %.2f MeV, gamma = %.3f, %d slot(s)", tag.c_str(), e, g, d.angSlots);
  Truth(std::string(head) + ": accepted", any);
  Check(std::string(head) + ": measure", sum, Brute(d, win, ks), tol);
  Truth(std::string(head) + ": nodes inside (angle, q)", inside);
}

}  // namespace

int main(int argc, char **argv) {
  for (int i = 1; i < argc; i++)
    if (!std::strcmp(argv[i], "-v")) verbose = true;

  std::printf("(a, b) the measure of the accepted directions\n");
  Window w;
  w.cm = true, w.lo = 20.0, w.hi = 70.0;
  Measure("c12 cm:20-70", C12(), w, 1.0);
  w.cm = false, w.lo = 0.0, w.hi = 180.0;
  Measure("c12 lab 0-180 (one branch)", C12(), w, 0.8);
  Measure("c12 lab 0-180 (two branches)", C12(), w, 2.4);
  w.lo = 5.0, w.hi = 15.0;
  Measure("c12 lab 5-15 (one branch)", C12(), w, 0.8);
  Measure("c12 lab 5-15 (two branches)", C12(), w, 2.4);
  w.lo = 20.0, w.hi = 60.0;
  Measure("c12 lab 20-60 (beyond the largest angle part)", C12(), w, 2.4);
  w.lo = 0.0, w.hi = 30.0, w.cut = true, w.qlo = 0.0, w.qhi = 40.0;
  Measure("c12 lab 0-30 with q <= 40 MeV/c", C12(), w, 2.4);
  w = Window();
  w.cm = true, w.lo = 120.0, w.hi = 180.0, w.cut = true, w.qlo = 10.0, w.qhi = 60.0;
  Measure("he3 cm:120-180 with 10 <= q <= 60 MeV/c", He3(), w, 0.6, 2.0e-5);  // the bins at the cut edges
  w = Window();
  w.cm = false, w.lo = 0.0, w.hi = 180.0;
  Measure("he3 lab 0-180", He3(), w, 0.6);

  std::printf("(c) a lab window of zero width on two branches\n");
  for (double th : {0.0, 10.0}) {
    ThmDistortion d, dt;
    Window z, t;
    z.lo = z.hi = th;
    // (At 0 deg the backward measure of [0, delta] is ((gamma - 1) delta)^2/2:
    // delta = 0.01 deg keeps it well above the round-off of 1 - cos.)
    t.lo = th == 0.0 ? 0.0 : th - 1.0e-6, t.hi = th + (th == 0.0 ? 1.0e-2 : 1.0e-6);
    if (!Make(d, C12(), z, 4, 2.3, 2.5) || !Make(dt, C12(), t, 4, 2.3, 2.5)) continue;
    const double ks = Ks(d, 2.4);
    std::vector<ThmDistortion::AngleNode> nz, nt;
    d.AngleNodes(ks, nz);
    dt.AngleNodes(ks, nt);
    double z0 = 0.0, z1 = 0.0, t0 = 0.0, t1 = 0.0;
    for (int i = 0; i < 4; i++) z0 += nz[i].w, z1 += nz[4 + i].w, t0 += nt[i].w, t1 += nt[4 + i].w;
    char what[120];
    std::snprintf(what, sizeof what, "lab %g deg: backward/forward weight == shrinking window", th);
    Check(what, z1 / z0, t1 / t0, 1.0e-5);
    if (th == 0.0) {
      const double g = Gamma(d, ks);
      Check("lab 0 deg: (gamma - 1)^2/(gamma + 1)^2", z1 / z0, (g - 1.0) * (g - 1.0) / ((g + 1.0) * (g + 1.0)), 1.0e-12);
    }
  }

  std::printf("(d) an acceptance table\n");
  w = Window();
  // Smooth in cos(theta_cm) inside each interval (a kink of the table, or
  // theta = 0 or 180 deg inside it, would be integrated to O(h^2) only).
  w.cm = true, w.t = {100.0, 170.0}, w.w = {0.2, 1.0};
  Measure("he3 cm table (ramp)", He3(), w, 0.6, 2.0e-5);
  w.cm = false, w.t = {2.0, 25.0}, w.w = {1.0, 0.1};
  Measure("c12 lab table (ramp, two branches)", C12(), w, 2.4, 2.0e-5);

  if (failures) {
    std::printf("FAILED: %d check(s)\n", failures);
    return 1;
  }
  std::printf("all spectator-direction window checks passed\n");
  return 0;
}
