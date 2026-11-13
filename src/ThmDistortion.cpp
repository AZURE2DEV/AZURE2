#include "ThmDistortion.h"
#include "Config.h"
#include "ThmExperiment.h"
#include "ThmNumerics.h"
#include "azure2/physics/transfer/GlobalOptical.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <gsl/gsl_integration.h>

/*
 * Numerics (docs/source/theory/thm_implementation.rst, "Distortion factor R(E)";
 * docs/upgrade/transfer.md section 24).
 *
 * Radial grid.  r = r_sx on nodes r_i = i h, h = 0.02 fm (smaller if a local
 * wave number times the step exceeds 0.1, or 0.2/kappa); the a + A wave is
 * needed at beta r_i, so the library solves it on its own grid of step beta h
 * and reads it at the same index.  A lower cutoff rmin of the bound state is a
 * node.  The radial integrals run from rmin to r_end = rmin + 50/kappa by
 * Simpson's rule (an even number of steps): the integrand is phi_sx(r)
 * u^sF(r) u^aA(beta r), and phi decays as e^{-kappa r} while the regular waves
 * stay bounded outside their turning points and are small inside, so the
 * integrals converge absolutely; the l = 0 integrand at r_end, over kappa, is
 * reported relative to |M|.
 *
 * Waves, bound state and the sum over l: the library (ZeroRangeWeight.h).
 */

namespace {

namespace transfer = azure2::physics::transfer;

const double kGridStep = 0.01;  // MeV, the ln R grid
const int kMaxRadial = 400000;  // radial nodes at most

// The ten numbers V,R,a,W,RW,aW,WD,RD,aD,RC (MeV, fm) of a potential of the global table.
void TenNumbers(const transfer::OpticalPotential &o, double p[10]) {
  const transfer::WoodsSaxonTerm *terms[3] = {&o.real_volume, &o.imaginary_volume, &o.imaginary_surface};
  for (int t = 0; t < 3; t++) {
    p[3 * t] = terms[t]->depth_MeV;
    p[3 * t + 1] = transfer::RadiusOf(o, terms[t]->r0_fm);
    p[3 * t + 2] = terms[t]->diffuseness_fm;
  }
  p[9] = transfer::RadiusOf(o, o.coulomb_r0_fm);
}

}  // namespace

// A nucleus by name if the built-in table has it, else (Z,A).
static std::string NucleusName(int Z, int A) {
  if (const ThmNuclide *n = ThmNuclide::Find(Z, A)) return n->name;
  std::ostringstream s;
  s << "(Z,A)=(" << Z << "," << A << ")";
  return s.str();
}

void ThmDistortion::Channel::SetEnergy(double ecm) {
  if (global < 0) return;
  TenNumbers(transfer::GlobalOpticalSets()[global].at(Zp, Ap, Zt, At, LabEnergy(ecm)), p);
}

std::string ThmDistortion::Channel::Describe() const {
  std::ostringstream s;
  s.precision(6);
  if (kind == PLANE) return "plane wave";
  if (kind == POINT_COULOMB) {
    s << "point Coulomb (Z1 Z2 = " << Z1 * Z2 << ")";
    return s.str();
  }
  if (global >= 0)
    s << transfer::GlobalOpticalSets()[global].name << (extrapolate ? ":extrapolate" : "") << " (" << NucleusName(Zp, Ap)
      << " on " << NucleusName(Zt, At) << ") ";
  s << "Woods-Saxon V=" << p[0] << " R=" << p[1] << " a=" << p[2] << ", W=" << p[3] << " RW=" << p[4]
    << " aW=" << p[5] << ", WD=" << p[6] << " RD=" << p[7] << " aD=" << p[8] << ", Coulomb "
    << (p[9] > 0.0 ? "RC=" + ThmNumberText(p[9]) : std::string("point")) << " (Z1 Z2 = " << Z1 * Z2 << ")";
  return s.str();
}

transfer::WaveChannel ThmDistortion::LibraryChannel(const Channel &c, double ecm) const {
  transfer::WaveChannel w;
  w.plane_wave = c.kind == Channel::PLANE;
  w.pair.z1 = c.Z1;
  w.pair.z2 = c.Z2;
  w.pair.mu_u = c.muU;
  w.energy_MeV = ecm;
  if (c.kind == Channel::WOODS_SAXON) {
    w.potential.real_volume = {c.p[0], c.p[1], c.p[2]};
    w.potential.imaginary_volume = {c.p[3], c.p[4], c.p[5]};
    w.potential.imaginary_surface = {c.p[6], c.p[7], c.p[8]};
    w.potential.coulomb_r0_fm = c.p[9];
  }
  return w;
}

transfer::BoundStateRequest ThmDistortion::BoundRequest() const {
  transfer::BoundStateRequest q;
  q.pair.z1 = kin.Zs;
  q.pair.z2 = kin.Zx;
  q.pair.mu_u = kin.ms * kin.mx / (kin.ms + kin.mx);
  q.binding_MeV = kin.bind;
  q.step_fm = h;
  q.r_max_fm = (n - 1) * h;
  return q;
}

std::string ThmDistortion::Amplitudes(double energy, const std::vector<double> &x, std::vector<complex> &m,
                                      Point &p) const {
  transfer::ZeroRangeWeight w;
  try {
    w = transfer::ZeroRangeWeightSum(overlap, entrance, LibraryChannel(SfAt(energy), EsF(energy)), x, numerics);
  } catch (const std::exception &e) {
    return "the s + F channel at E = " + ThmNumberText(energy) + " MeV: " + e.what();
  }
  const double kb = beta * aa.k;
  m.assign(x.size(), complex(0.0, 0.0));
  double smallest = HUGE_VAL;
  for (size_t k = 0; k < x.size(); k++) {
    m[k] = 4.0 * M_PI / (p.ksf * kb) * w.sum[k];
    smallest = std::min(smallest, std::abs(w.sum[k]));
  }
  p.lmax = w.l_max;
  p.tail = smallest > 0.0 ? w.end_integrand / kappa / smallest : HUGE_VAL;
  return "";
}

std::string ThmDistortion::CheckEnergy(double energy) const {
  std::ostringstream why;
  double esf = EsF(energy);
  if (!(esf > 0.0)) {
    why << "at E = " << energy << " MeV the spectator has no energy left (E_sF = E_aA - B - E = " << eAA << " - "
        << kin.bind << " - " << energy << " MeV <= 0); check Ebeam";
    return why.str();
  }
  double ksf = std::sqrt(2.0 * sf.mu * esf) / hbarc;
  double eta = sf.kind == Channel::PLANE ? 0.0 : kin.Zs * (kin.Zx + kin.ZA) * fstruc * sf.mu / (hbarc * ksf);
  if (eta > 50.0) {
    why << "at E = " << energy << " MeV the spectator is too far below the s + F barrier (eta_sF = " << eta
        << " > 50)";
    return why.str();
  }
  if (angleKind == LAB) {
    double vs = hbarc * ksf / (kin.ms * uconv);
    double s = vcm / vs * std::sin(angle * M_PI / 180.0);
    if (s > 1.0) {
      why << "at E = " << energy << " MeV the lab angle " << angle
          << " deg is beyond the reach of the spectator (max " << std::asin(vs / vcm) * 180.0 / M_PI << " deg)";
      return why.str();
    }
  }
  return "";
}

double ThmDistortion::SpectatorCos(double ksf, double *thetaCm, bool *clamped) const {
  // Spectator direction in the c.m. (angle to the beam) and x = cos(k_sF, k_aA).
  double theta = 0.0;
  if (angleKind == QF) {
    theta = kin.horseIsBeam ? 0.0 : M_PI;
  } else if (angleKind == CM) {
    theta = angle * M_PI / 180.0;
  } else {
    double vs = hbarc * ksf / (kin.ms * uconv);
    double tl = angle * M_PI / 180.0;
    double s = vcm / vs * std::sin(tl);
    if (s > 1.0) {
      s = 1.0;
      if (clamped) *clamped = true;
    }
    theta = tl + std::asin(s);
  }
  if (thetaCm) *thetaCm = theta * 180.0 / M_PI;
  return (kin.horseIsBeam ? 1.0 : -1.0) * std::cos(theta);
}

double ThmDistortion::PlaneWaveSource(double q) const { return transfer::ZeroRangeMomentumAmplitude(overlap, q); }

ThmDistortion::Point ThmDistortion::Evaluate(double energy) const {
  if (angWindow) return EvaluateWindow(energy);
  Point p;
  p.energy = energy;
  p.esf = EsF(energy);
  if (!(p.esf > 0.0)) {
    p.why = CheckEnergy(energy);
    return p;
  }
  p.ksf = std::sqrt(2.0 * sf.mu * p.esf) / hbarc;
  p.etasf = sf.kind == Channel::PLANE ? 0.0 : kin.Zs * (kin.Zx + kin.ZA) * fstruc * sf.mu / (hbarc * p.ksf);
  p.x = SpectatorCos(p.ksf, &p.thetaCm);
  const double kb = beta * aa.k;
  p.q = std::sqrt(std::max(0.0, p.ksf * p.ksf + kb * kb - 2.0 * p.ksf * kb * p.x));
  p.mpw = PlaneWaveSource(p.q);
  std::vector<complex> m;
  p.why = Amplitudes(energy, {p.x}, m, p);
  if (!p.why.empty()) return p;
  p.m = m[0];
  p.ok = std::isfinite(std::abs(p.m)) && std::abs(p.m) > 0.0 && p.mpw != 0.0;
  if (!p.ok) {
    std::ostringstream why;
    why << "the amplitude at E = " << energy << " MeV is zero or not finite";
    p.why = why.str();
  }
  return p;
}

void ThmDistortion::SetAngleSlots(double eHi) {
  angSlots = 1;
  if (!angWindow || angCm) return;
  const double esf = EsF(eHi);
  if (!(esf > 0.0)) return;
  const double vs = hbarc * (std::sqrt(2.0 * sf.mu * esf) / hbarc) / (kin.ms * uconv);
  if (vcm > vs) angSlots = 2;
}

bool ThmDistortion::AngleNodes(double ksf, std::vector<AngleNode> &out) const {
  const int N = angNodes;
  out.assign((size_t)angSlots * N, AngleNode());
  const double d2r = M_PI / 180.0;
  const double vs = hbarc * ksf / (kin.ms * uconv);
  const double g = vs > 0.0 ? vcm / vs : HUGE_VAL;
  // The accepted c.m. intervals [t0, t1] (rad) before the q cut, and for a
  // zero-width lab window the weight of each branch, |d cos(theta_cm)/d theta_lab|.
  double iv[2][2], dw[2] = {1.0, 1.0};
  int nIv = 0;
  if (angCm) {
    iv[0][0] = angLo * d2r;
    iv[0][1] = angHi * d2r;
    nIv = 1;
  } else {
    const double tm = g > 1.0 ? std::asin(1.0 / g) : (g == 1.0 ? 0.5 * M_PI : M_PI);
    const double lo = angLo * d2r, hi = std::min(angHi * d2r, tm);
    auto as = [&](double t) { return std::asin(std::min(1.0, g * std::sin(t))); };
    if (lo <= hi) {
      iv[0][0] = lo + as(lo);
      iv[0][1] = hi + as(hi);
      nIv = 1;
      if (g > 1.0) {
        iv[1][0] = hi + M_PI - as(hi);
        iv[1][1] = lo + M_PI - as(lo);
        nIv = 2;
      }
      // d theta_cm/d theta_lab = 1 +- g cos(t)/sqrt(1 - g^2 sin^2 t) at t = lo
      // (at t = 0 both measures vanish like t dt, in the ratio (g + 1)^2 : (g - 1)^2).
      if (nIv == 2 && lo < 1.0e-9) {
        dw[0] = (g + 1.0) * (g + 1.0);
        dw[1] = (g - 1.0) * (g - 1.0);
      } else if (nIv == 2) {
        const double r = std::sqrt(std::max(1.0e-300, 1.0 - std::pow(g * std::sin(lo), 2)));
        dw[0] = std::fabs(std::sin(iv[0][0]) * (1.0 + g * std::cos(lo) / r));
        dw[1] = std::fabs(std::sin(iv[1][1]) * (1.0 - g * std::cos(lo) / r));
      }
    }
  }
  // The q cut as a range of cos(theta_cm).
  double cLo = -1.0, cHi = 1.0;
  const double kb = beta * aa.k;
  if (qCut && ksf > 0.0 && kb > 0.0) {
    const double ql = qCutLo / hbarc, qh = qCutHi / hbarc;
    const double xl = (ksf * ksf + kb * kb - qh * qh) / (2.0 * ksf * kb);
    const double xh = (ksf * ksf + kb * kb - ql * ql) / (2.0 * ksf * kb);
    cLo = kin.horseIsBeam ? xl : -xh;
    cHi = kin.horseIsBeam ? xh : -xl;
    cLo = std::max(cLo, -1.0);
    cHi = std::min(cHi, 1.0);
  }
  double c[2][2];
  bool open[2] = {false, false}, positive = false;
  for (int i = 0; i < nIv; i++) {
    c[i][0] = std::max(std::cos(std::min(M_PI, iv[i][1])), cLo);
    c[i][1] = std::min(std::cos(std::max(0.0, iv[i][0])), cHi);
    open[i] = c[i][0] <= c[i][1];
    positive = positive || (open[i] && c[i][1] > c[i][0]);
  }
  const double sgn = kin.horseIsBeam ? 1.0 : -1.0;
  bool any = false;
  for (int i = 0; i < nIv && i < angSlots; i++) {
    if (!open[i]) continue;
    const double half = 0.5 * (c[i][1] - c[i][0]), mid = 0.5 * (c[i][1] + c[i][0]);
    if (positive && !(half > 0.0)) continue;  // measure zero next to a finite interval
    for (int k = 0; k < N; k++) {
      AngleNode &n = out[(size_t)i * N + k];
      const double ct = half > 0.0 ? mid + half * angGx[k] : mid;
      const double theta = std::acos(std::max(-1.0, std::min(1.0, ct)));
      n.theta = theta * 180.0 / M_PI;
      n.x = sgn * ct;
      n.q = std::sqrt(std::max(0.0, ksf * ksf + kb * kb - 2.0 * ksf * kb * n.x));
      double a = 1.0;
      if (!angT.empty()) {
        // Acceptance at the table's angle: c.m., or the lab angle of this direction.
        double t = angCm ? n.theta : std::atan2(std::sin(theta), std::cos(theta) + g) * 180.0 / M_PI;
        t = std::max(angT.front(), std::min(angT.back(), t));
        size_t j = std::upper_bound(angT.begin(), angT.end(), t) - angT.begin();
        j = std::max<size_t>(1, std::min(j, angT.size() - 1));
        const double f = (t - angT[j - 1]) / (angT[j] - angT[j - 1]);
        a = angW[j - 1] + f * (angW[j] - angW[j - 1]);
      }
      n.w = (half > 0.0 ? half * angGw[k] : dw[i] / N) * a;
      any = any || n.w > 0.0;
    }
  }
  // Empty intervals copy the directions of the first open one (weight 0), so
  // that tables interpolated across energies stay finite.
  for (int i = 0; i < angSlots; i++) {
    bool empty = true;
    for (int k = 0; k < N; k++) empty = empty && out[(size_t)i * N + k].w == 0.0 && out[(size_t)i * N + k].q == 0.0;
    if (!empty) continue;
    for (int j = 0; j < angSlots; j++) {
      if (j == i) continue;
      bool full = false;
      for (int k = 0; k < N; k++) full = full || out[(size_t)j * N + k].q != 0.0;
      if (!full) continue;
      for (int k = 0; k < N; k++) {
        out[(size_t)i * N + k] = out[(size_t)j * N + k];
        out[(size_t)i * N + k].w = 0.0;
      }
      break;
    }
  }
  return any;
}

std::string ThmDistortion::CheckWindow(double energy) const {
  if (!angWindow) return "";
  const double esf = EsF(energy);
  if (!(esf > 0.0)) return CheckEnergy(energy);
  std::vector<AngleNode> nodes;
  if (AngleNodes(std::sqrt(2.0 * sf.mu * esf) / hbarc, nodes)) return "";
  std::ostringstream why;
  if (angAll) {
    const double ks = std::sqrt(2.0 * sf.mu * esf) / hbarc, kb = beta * aa.k;
    why << "at E = " << energy << " MeV the spectator momenta the kinematics reach, " << std::fabs(ks - kb) * hbarc
        << " to " << (ks + kb) * hbarc << " MeV/c, do not overlap the ps window [" << qCutLo << ", " << qCutHi
        << "] MeV/c";
    return why.str();
  }
  why << "at E = " << energy << " MeV no spectator direction of spectatorAngles=" << AngleText()
      << " is accepted";
  const double vs = hbarc * std::sqrt(2.0 * sf.mu * esf) / hbarc / (kin.ms * uconv);
  if (!angCm && vcm > vs)
    why << " (the largest lab angle the spectator reaches is " << std::asin(vs / vcm) * 180.0 / M_PI << " deg)";
  if (qCut) why << " with |p_s| in the ps window [" << qCutLo << ", " << qCutHi << "] MeV/c";
  return why.str();
}

std::string ThmDistortion::AngleText() const {
  if (!angWindow) return "";
  std::ostringstream t;
  t.precision(6);
  if (angAll) {
    t << "every direction with |p_s| in [" << qCutLo << ", " << qCutHi << "] MeV/c";
    return t.str();
  }
  t << (angCm ? "cm:" : "") << angLo << "-" << angHi;
  if (!angT.empty()) t << " (acceptance table, " << angT.size() << " rows)";
  return t.str();
}

ThmDistortion::Point ThmDistortion::EvaluateWindow(double energy) const {
  Point p;
  p.energy = energy;
  p.esf = EsF(energy);
  if (!(p.esf > 0.0)) {
    p.why = CheckEnergy(energy);
    return p;
  }
  p.ksf = std::sqrt(2.0 * sf.mu * p.esf) / hbarc;
  p.etasf = sf.kind == Channel::PLANE ? 0.0 : kin.Zs * (kin.Zx + kin.ZA) * fstruc * sf.mu / (hbarc * p.ksf);
  std::vector<AngleNode> nodes;
  if (!AngleNodes(p.ksf, nodes)) {
    p.why = CheckWindow(energy);
    return p;
  }
  // The active nodes.
  std::vector<AngleNode> act;
  for (const AngleNode &n : nodes)
    if (n.w > 0.0) act.push_back(n);
  const int nk = (int)act.size();
  // Plane-wave limits.
  std::vector<double> mpw(nk, 0.0);
  for (int k = 0; k < nk; k++) mpw[k] = PlaneWaveSource(act[k].q);
  // The radial integrals once, P_l(x_k) per node (the library's sum at every direction).
  std::vector<double> xk(nk);
  for (int k = 0; k < nk; k++) xk[k] = act[k].x;
  std::vector<complex> mk;
  p.why = Amplitudes(energy, xk, mk, p);
  if (!p.why.empty()) return p;
  double wsum = 0.0, m2 = 0.0, mpw2 = 0.0, th = 0.0, xs = 0.0, qs = 0.0;
  bool finite = true;
  for (int k = 0; k < nk; k++) {
    const complex m = mk[k];
    finite = finite && std::isfinite(std::abs(m)) && mpw[k] != 0.0;
    wsum += act[k].w;
    m2 += act[k].w * std::norm(m);
    mpw2 += act[k].w * mpw[k] * mpw[k];
    th += act[k].w * act[k].theta;
    xs += act[k].w * act[k].x;
    qs += act[k].w * act[k].q;
    if (k == 0) {
      p.m = m;
      p.mpw = mpw[k];
    }
  }
  p.nodes = nk;
  p.m2 = m2 / wsum;
  p.mpw2 = mpw2 / wsum;
  p.thetaCm = th / wsum;
  p.x = xs / wsum;
  p.q = qs / wsum;
  p.ok = finite && p.m2 > 0.0 && p.mpw2 > 0.0 && std::isfinite(p.m2);
  if (!p.ok) {
    std::ostringstream why;
    why << "the amplitude at E = " << energy << " MeV is zero or not finite";
    p.why = why.str();
  }
  return p;
}

std::string ThmDistortion::Setup(const ThmExperiment &x, const Kinematics &k, double eLo) {
  experiment = x.name;
  kin = k;
  kind = x.distortion == ThmExperiment::DIST_OPTICAL ? OPTICAL : COULOMB;
  angleKind = x.angleKind == 1 ? LAB : x.angleKind == 2 ? CM : QF;
  angle = x.angle;
  // The acceptance: spectatorAngles= and/or the |p_s| cut of a ps window;
  // a ps window alone accepts every direction whose q lies in it (angAll),
  // with psNodes nodes (one for a window of zero width).
  angWindow = x.angleWindow != 0 || x.psKind != ThmExperiment::PS_DELTA;
  angAll = angWindow && x.angleWindow == 0;
  angCm = angAll ? true : x.angleCm;
  angLo = angAll ? 0.0 : x.angleMin;
  angHi = angAll ? 180.0 : x.angleMax;
  angT = x.angleWindow == 2 ? x.angleTableT : std::vector<double>();
  angW = x.angleWindow == 2 ? x.angleTableW : std::vector<double>();
  angNodes = !angAll ? x.angleNodes : x.psMin == x.psMax ? 1 : x.psNodes;
  qCut = x.psKind != ThmExperiment::PS_DELTA;
  qCutLo = x.psMin;
  qCutHi = x.psMax;
  angGx.assign(angNodes, 0.0);
  angGw.assign(angNodes, 0.0);
  if (angWindow) {
    gsl_integration_glfixed_table *t = gsl_integration_glfixed_table_alloc(angNodes);
    for (int i = 0; i < angNodes; i++) gsl_integration_glfixed_point(-1.0, 1.0, i, &angGx[i], &angGw[i], t);
    gsl_integration_glfixed_table_free(t);
    if (angNodes == 1) angGw[0] = 2.0;
  }
  ratioPW = x.distortionRatioPW;
  yukawa = x.boundYukawa;
  if (!(k.bind > 0.0)) return "the Trojan horse is not bound: B(x+s) = " + ThmNumberText(k.bind) + " MeV";
  eAA = k.beamEnergy * k.mTarget / (k.mBeam + k.mTarget);
  vcm = std::sqrt(2.0 * k.mBeam * uconv * k.beamEnergy) / ((k.mBeam + k.mTarget) * uconv);
  beta = k.ms / k.ma;
  muSx = k.ms * k.mx / (k.ms + k.mx) * uconv;
  kappa = std::sqrt(2.0 * muSx * k.bind) / hbarc;
  etaB = k.Zs * k.Zx * fstruc * muSx / (hbarc * kappa);
  warnings.clear();
  const char *keys[2] = {"opticalAA", "opticalSF"}, *labels[2] = {"a + A", "s + F"};
  auto channel = [&](int which, const ThmExperiment::Optical &o, int z1, int a1, double m1, int z2, int a2,
                     double m2, Channel &c) -> std::string {
    c = Channel();
    c.kind = kind == COULOMB ? Channel::POINT_COULOMB
             : o.kind == 0   ? Channel::PLANE
             : o.kind == 1   ? Channel::POINT_COULOMB
                             : Channel::WOODS_SAXON;
    for (int t = 0; t < 10; t++) c.p[t] = o.p[t];
    c.Z1 = z1;
    c.Z2 = z2;
    c.muU = m1 * m2 / (m1 + m2);
    c.mu = c.muU * uconv;
    if (kind != OPTICAL || o.kind != 3) return "";
    // A global potential: its projectile is the partner it is made for.
    const transfer::GlobalOpticalSet &g = transfer::GlobalOpticalSets()[o.global];
    const std::string what = std::string(keys[which]) + "=" + g.name + " (" + g.reference + ", for " +
                             g.projectiles + "): the " + labels[which] + " channel is " + NucleusName(z1, a1) +
                             " + " + NucleusName(z2, a2);
    if (a1 <= 0 || a2 <= 0) return std::string(keys[which]) + "=" + g.name + ": the mass numbers are not known";
    const bool first = transfer::Describes(g, z1, a1);
    if (!first && !transfer::Describes(g, z2, a2))
      return what + ", which it does not describe (heavy-ion and other channels take the ten numbers "
                    "V,R,a,W,RW,aW,WD,RD,aD,RC)";
    c.global = o.global;
    c.extrapolate = o.extrapolate;
    c.Zp = first ? z1 : z2;
    c.Ap = first ? a1 : a2;
    c.mp = first ? m1 : m2;
    c.Zt = first ? z2 : z1;
    c.At = first ? a2 : a1;
    c.mt = first ? m2 : m1;
    return "";
  };
  std::string bad = channel(0, x.opticalAA, k.Za, k.Aa, k.ma, k.ZA, k.AA, k.mA, aa);
  if (bad.empty()) bad = channel(1, x.opticalSF, k.Zs, k.As, k.ms, k.Zx + k.ZA, k.Ax + k.AA, k.mx + k.mA, sf);
  if (!bad.empty()) return bad;
  aa.k = std::sqrt(2.0 * aa.mu * eAA) / hbarc;
  aa.eta = aa.kind == Channel::PLANE ? 0.0 : aa.Z1 * aa.Z2 * fstruc * aa.mu / (hbarc * aa.k);
  // Global potentials: a + A at E_aA; s + F follows E_sF, so the radial step
  // takes the largest depth over (0, E_sF(eLo)]; both checked against the
  // validity range over the data.
  const double lo = dataLo <= dataHi ? dataLo : eLo, hi = dataLo <= dataHi ? dataHi : eLo;
  for (int which = 0; which < 2; which++) {
    Channel &c = which == 0 ? aa : sf;
    if (c.global < 0) continue;
    const transfer::GlobalOpticalSet &g = transfer::GlobalOpticalSets()[c.global];
    double e0 = which == 0 ? eAA : EsF(hi), e1 = which == 0 ? eAA : EsF(lo);
    double l0 = c.LabEnergy(e0), l1 = c.LabEnergy(e1);
    if (which == 0) {
      c.SetEnergy(eAA);
      c.depthMax = std::fabs(c.p[0]) + std::fabs(c.p[3]) + std::fabs(c.p[6]);
    } else {
      const double top = std::max(EsF(eLo), e1);
      for (int i = 1; i <= 40; i++) {
        c.SetEnergy(top * i / 40.0);
        c.depthMax = std::max(c.depthMax, std::fabs(c.p[0]) + std::fabs(c.p[3]) + std::fabs(c.p[6]));
      }
      c.SetEnergy(0.5 * (e0 + e1));
    }
    const bool massOk = c.At >= g.a_min && c.At <= g.a_max;
    const bool energyOk = l0 >= g.e_min_MeV && l1 <= g.e_max_MeV;
    if (massOk && energyOk) continue;
    std::ostringstream why;
    why.precision(4);
    why << keys[which] << "=" << g.name << " (" << g.reference << ") is outside its validity range for "
        << NucleusName(c.Zp, c.Ap) << " + " << NucleusName(c.Zt, c.At) << ":";
    if (!massOk) why << " target A = " << c.At << " (valid " << g.a_min << "-" << g.a_max << ")";
    if (!energyOk) {
      why << (massOk ? "" : ",") << " lab energy of " << NucleusName(c.Zp, c.Ap) << " " << l0;
      if (which == 1 && l1 != l0) why << "-" << l1;
      why << " MeV over the data (valid " << g.e_min_MeV << "-" << g.e_max_MeV << " MeV)";
    }
    if (c.extrapolate) {
      warnings.push_back(why.str() + "; extrapolated as asked (:extrapolate).");
      continue;
    }
    return why.str() + ". Write " + g.name + ":extrapolate to use it there anyway (warned), or give the ten numbers";
  }

  // Radial step: h = 0.02 fm, or k_local h <= 0.1 in both channels.
  double kSF = std::sqrt(2.0 * sf.mu * (std::max(EsF(eLo), 0.0) + Depth(sf))) / hbarc;
  double kAA = std::sqrt(2.0 * aa.mu * (eAA + Depth(aa))) / hbarc;
  h = std::min(maxStep, std::min(phasePerStep / std::max(kSF, 1.0e-6), phasePerStep / std::max(beta * kAA, 1.0e-6)));
  h = std::min(h, stepKappa / kappa);
  rmin = x.boundRmin;
  i0 = 0;
  if (rmin > 0.0) {
    i0 = (int)std::ceil(rmin / h - 1.0e-9);
    h = rmin / i0;
  }
  return "";
}

double ThmDistortion::Depth(const Channel &c) {
  if (c.kind == Channel::WOODS_SAXON && c.global >= 0) return c.depthMax;
  return c.kind == Channel::WOODS_SAXON ? std::fabs(c.p[0]) + std::fabs(c.p[3]) + std::fabs(c.p[6]) : 0.0;
}

std::string ThmDistortion::ChannelText(int which) const {
  const Channel &c = which == 0 ? aa : sf;
  if (c.global < 0) return c.Describe();
  std::ostringstream t;
  t.precision(6);
  if (which == 0) {
    t << c.Describe() << " at E_lab = " << c.LabEnergy(eAA) << " MeV";
    return t.str();
  }
  // The s + F potential follows E_sF: at the two ends of the data.
  const double lo = dataLo <= dataHi ? dataLo : 0.0, hi = dataLo <= dataHi ? dataHi : 0.0;
  Channel c0 = SfAt(lo), c1 = SfAt(hi);
  t << c1.Describe() << " at E_lab = " << c.LabEnergy(EsF(hi)) << " MeV (E = " << hi << " MeV) to V=" << c0.p[0]
    << " R=" << c0.p[1] << " a=" << c0.p[2] << ", W=" << c0.p[3] << " RW=" << c0.p[4] << " aW=" << c0.p[5]
    << ", WD=" << c0.p[6] << " RD=" << c0.p[7] << " aD=" << c0.p[8] << " at E_lab = " << c.LabEnergy(EsF(lo))
    << " MeV (E = " << lo << " MeV)";
  return t.str();
}

ThmDistortion::Channel ThmDistortion::SfAt(double energy) const {
  Channel c = sf;
  const double esf = EsF(energy);
  c.k = std::sqrt(2.0 * c.mu * std::max(esf, 0.0)) / hbarc;
  c.eta = c.kind == Channel::PLANE ? 0.0 : kin.Zs * (kin.Zx + kin.ZA) * fstruc * c.mu / (hbarc * c.k);
  c.SetEnergy(esf);
  return c;
}

bool ThmDistortion::GlobalEnds(int which, double lo, double hi, double elab[2], double p[2][10]) const {
  const Channel &c = which == 0 ? aa : sf;
  if (c.global < 0) return false;
  for (int end = 0; end < 2; end++) {
    const double ecm = which == 0 ? eAA : EsF(end == 0 ? lo : hi);
    elab[end] = c.LabEnergy(ecm);
    TenNumbers(transfer::GlobalOpticalSets()[c.global].at(c.Zp, c.Ap, c.Zt, c.At, elab[end]), p[end]);
  }
  return true;
}

double ThmDistortion::Phi(double r) const {
  if (!(r > 0.0)) return 0.0;
  try {
    return transfer::AnalyticTailAt(Shape(), BoundRequest(), 1.0, rmin, r) / r;
  } catch (const std::exception &) {
    return 0.0;  // the Whittaker function failed (far beyond the double range)
  }
}

std::string ThmDistortion::Build(const ThmExperiment &x, const Kinematics &k, double eLo, double eHi,
                                 double eRefDefault) {
  if (!(dataLo <= dataHi)) {
    dataLo = eLo;
    dataHi = eHi;
  }
  std::string setup = Setup(x, k, eLo);
  if (!setup.empty()) return setup;
  const double rEnd = rmin + tailKappa / kappa;
  int intervals = (int)std::ceil((rEnd - rmin) / h);
  if (intervals % 2) intervals++;
  n = i0 + intervals + 1;
  if (n > kMaxRadial) return "the radial grid would need " + ThmNumberText(n) + " points (kappa too small)";
  // The bound state u_sx = r phi_sx (unit ANC) and the a + A side of the sum, once.
  numerics.c_i = beta;
  numerics.c_f = 1.0;
  try {
    overlap = transfer::RadialOverlapOf(transfer::AnalyticTail(Shape(), BoundRequest(), 1.0, rmin));
    entrance = transfer::PrepareZeroRangeEntrance(overlap, LibraryChannel(aa, eAA), numerics);
  } catch (const std::exception &e) {
    return std::string("the a + A channel: ") + e.what();
  }

  SetAngleSlots(eHi);
  std::ostringstream d;
  d.precision(6);
  d << (kind == COULOMB ? "coulomb" : "optical") << "; a + A: " << ChannelText(0) << "; s + F: " << ChannelText(1)
    << (angWindow ? "; spectator directions " + AngleText() + " (" + ThmNumberText(angNodes) + " nodes in cos theta_cm" +
                        (angSlots > 1 ? " per branch" : "") + (qCut ? ", |p_s| cut by the ps window" : "") +
                        "), acceptance-averaged |M|^2 and |M_PW|^2"
                  : "; spectator angle " +
                        (angleKind == QF ? std::string("qf (k_sF along k_aA)")
                                         : (angleKind == LAB ? "lab " : "cm ") + ThmNumberText(angle) + " deg"))
    << "; bound state " << (yukawa ? "yukawa" : "whittaker") << (rmin > 0.0 ? ", r >= " + ThmNumberText(rmin) + " fm" : "")
    << "; R = " << (ratioPW ? "(|M|^2/|M_PW|^2)(E) / (same)(E_ref)" : "|M(E)|^2/|M(E_ref)|^2");
  description = d.str();

  // Reference.
  eRef = x.hasDistortionRef ? x.distortionRef : eRefDefault;
  std::string why = CheckEnergy(eRef);
  if (!why.empty()) return "distortionRef: " + why;
  ref = Evaluate(eRef);
  if (!ref.ok) return "distortionRef: " + ref.why;
  refRatio = Ratio(ref);
  if (!(refRatio > 0.0) || !std::isfinite(refRatio)) return "the ratio at E_ref is zero or not finite";

  // ln R on the grid.
  int nodes = std::max(4, (int)std::ceil((eHi - eLo) / kGridStep - 1.0e-9) + 1);
  gridLo = eLo;
  gridStep = kGridStep;
  lnR.assign(nodes, 0.0);
  std::vector<double> tails(nodes, 0.0), pw(nodes, 0.0);
  std::vector<std::string> errors(nodes);
  // A window can accept no direction at grid energies beyond the data
  // (WorkingCopy refuses data points without one): those take the nearest value.
  std::vector<char> empty(nodes, 0);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < nodes; i++) {
    if (angWindow && !CheckWindow(gridLo + i * gridStep).empty() && EsF(gridLo + i * gridStep) > 0.0) {
      empty[i] = 1;
      continue;
    }
    Point p = Evaluate(gridLo + i * gridStep);
    if (!p.ok) {
      errors[i] = p.why;
      continue;
    }
    lnR[i] = std::log(Ratio(p) / refRatio);
    tails[i] = p.tail;
    pw[i] = p.mpw;
  }
  for (int i = 0; i < nodes; i++)
    if (!errors[i].empty()) return errors[i];
  if (angWindow) {
    int best = -1;
    for (int i = 0; i < nodes; i++)
      if (!empty[i]) best = i;
    if (best < 0) return "no spectator direction of the window is accepted on the grid";
    for (int i = 0; i < nodes; i++) {
      if (!empty[i]) continue;
      int near = -1;
      for (int j = 0; j < nodes; j++)
        if (!empty[j] && (near < 0 || std::abs(j - i) < std::abs(near - i))) near = j;
      lnR[i] = lnR[near];
      pw[i] = pw[near];
    }
  }
  tailWorst = std::max(ref.tail, *std::max_element(tails.begin(), tails.end()));
  // One direction: rho = |M|^2/M_PW^2 is singular where M_PW changes sign.  A
  // window of directions divides by <M_PW^2>, which a sign change of one
  // node's M_PW does not make vanish.
  pwSignChange = false;
  if (!angWindow)
    for (int i = 0; i < nodes; i++) pwSignChange = pwSignChange || (pw[i] > 0.0) != (ref.mpw > 0.0);
  return "";
}

double ThmDistortion::Weight(double energy, bool *outside) const {
  if (kind == TABLE) return (*table)(energy, outside);
  const int nodes = (int)lnR.size();
  double t = (energy - gridLo) / gridStep;
  bool out = t < 0.0 || t > nodes - 1.0;
  if (outside) *outside = out;
  if (t <= 0.0) return std::exp(lnR.front());
  if (t >= nodes - 1.0) return std::exp(lnR.back());
  // Cubic Lagrange through the four nearest nodes.
  int j = std::min(std::max((int)std::floor(t), 1), nodes - 3);
  double c[4];
  ThmCubicLagrange(t - j, c);  // nodes j-1, j, j+1, j+2 at s = -1, 0, 1, 2
  double v = c[0] * lnR[j - 1] + c[1] * lnR[j] + c[2] * lnR[j + 1] + c[3] * lnR[j + 2];
  return std::exp(v);
}
