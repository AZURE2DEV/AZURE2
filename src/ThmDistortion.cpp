#include "ThmDistortion.h"
#include "Config.h"
#include "GSLException.h"
#include "ThmExperiment.h"
#include "ThmOptical.h"
#include "cwfcomp_accurate.H"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>

/*
 * Numerics (docs/source/theory/thm_implementation.rst, "Distortion factor R(E)").
 *
 * Radial grid.  r = r_sx on nodes r_i = i h, h = 0.02 fm (smaller if a local
 * wave number times the step exceeds 0.1); the a + A wave is needed at
 * beta r_i, so it is integrated on its own grid of step beta h and read at the
 * same index.  A lower cutoff rmin of the bound state is a node.  The radial
 * integrals run from rmin to r_end = rmin + 50/kappa by Simpson's rule: the
 * integrand is phi_sx(r) u^sF(r) u^aA(beta r), and phi decays as e^{-kappa r}
 * while the regular waves stay bounded outside their turning points and
 * are small inside, so the integrals converge absolutely (no complex
 * rotation or damping is needed, unlike the Coulomb term C_l of the vertex);
 * the l = 0 integrand at r_end, over kappa, is reported relative to |M|.
 *
 * Waves.  u'' = [l(l+1)/r^2 + 2 mu (V_C + U)/hbar^2 - k^2] u, integrated
 * outward from the origin by Numerov's method (complex for an optical
 * potential), started from the power series r^{l+1}(1 + a1 r + a2 r^2 +
 * a3 r^3) at the first two nodes and rescaled on the way (the regular
 * solution is the growing one under a barrier, so the outward integration is
 * stable).  Beyond the potential and the turning point
 * (eta + sqrt(eta^2 + l(l+1)))/k it is matched at two nodes a quarter
 * wavelength apart to c1 F_l + c2 G_l, with F_l, G_l from COUL
 * (Coulomb_wave_functions_accurate, as CoulFunc), and divided by
 * N = c1 - i c2, which normalizes it to F_l + T_l H_l^+ (for point Coulomb
 * c2 = 0 up to the Numerov error, and u = F_l).
 *
 * Partial waves.  The a + A waves are computed once, for l = 0, 1, ... until
 * A_l = (2l+1) Int |phi u_l^aA| has fallen below 1e-17 of its largest value;
 * per energy the sum stops once sum_{l' > l} A_l' times max|u_l^sF| is below
 * 1e-13 |M|.
 */

namespace {

const double kGridStep = 0.01;  // MeV, the ln R grid

double WoodsSaxon(double r, double R, double a) {
  double x = (r - R) / a;
  return x > 700.0 ? 0.0 : 1.0 / (1.0 + std::exp(x));
}
// 4 e^x/(1 + e^x)^2 = -4 a d/dr of the Woods-Saxon form factor.
double SurfaceShape(double r, double R, double a) {
  double x = (r - R) / a;
  if (std::fabs(x) > 700.0) return 0.0;
  double e = std::exp(-std::fabs(x));
  return 4.0 * e / ((1.0 + e) * (1.0 + e));
}

// Coulomb phases sigma_l = arg Gamma(l + 1 + i eta), l = 0..lmax, as e^{i sigma_l}.
std::vector<complex> CoulombPhases(double eta, int lmax) {
  std::vector<complex> out(lmax + 1, complex(1.0, 0.0));
  if (eta == 0.0) return out;
  gsl_sf_result lnr, arg;
  gsl_sf_lngamma_complex_e(1.0, eta, &lnr, &arg);
  double sigma = arg.val;
  for (int l = 0; l <= lmax; l++) {
    if (l > 0) sigma += std::atan(eta / l);
    out[l] = std::polar(1.0, sigma);
  }
  return out;
}

// F_l, G_l at rho (eta real).
void CoulombFG(int l, double eta, double rho, double &F, double &G) {
  Coulomb_wave_functions_accurate coul(true, complex((double)l, 0.0), complex(eta, 0.0));
  complex f, df, g, dg;
  coul.compute(complex(rho, 0.0), f, df, g, dg);
  F = f.real();
  G = g.real();
}

std::string Number(double x) {
  std::ostringstream s;
  s.precision(6);
  s << x;
  return s.str();
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
  ThmGlobalOpticalEvaluate(global, Zp, Ap, Zt, At, LabEnergy(ecm), p);
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
    s << ThmGlobalOpticals()[global].name << (extrapolate ? ":extrapolate" : "") << " (" << NucleusName(Zp, Ap)
      << " on " << NucleusName(Zt, At) << ") ";
  s << "Woods-Saxon V=" << p[0] << " R=" << p[1] << " a=" << p[2] << ", W=" << p[3] << " RW=" << p[4]
    << " aW=" << p[5] << ", WD=" << p[6] << " RD=" << p[7] << " aD=" << p[8] << ", Coulomb "
    << (p[9] > 0.0 ? "RC=" + Number(p[9]) : std::string("point")) << " (Z1 Z2 = " << Z1 * Z2 << ")";
  return s.str();
}

bool ThmDistortion::Wave(const Channel &c, int l, double step, int nStore, std::vector<complex> &u,
                         complex *T) const {
  const double k = c.k;
  const bool plane = c.kind == Channel::PLANE;
  const bool ws = c.kind == Channel::WOODS_SAXON;
  const double eta = plane ? 0.0 : c.eta;
  const double fac = 2.0 * c.mu / (hbarc * hbarc);  // fm^-2 MeV^-1
  const double zz = plane ? 0.0 : c.Z1 * c.Z2 * fstruc * hbarc;  // MeV fm
  const double rc = ws ? c.p[9] : 0.0;
  const bool sphere = rc > 0.0;
  // Potential range (fm): where the Woods-Saxon terms are below e^-12.
  double rPot = rc;
  if (ws)
    for (int t = 0; t < 3; t++)
      if (c.p[3 * t] != 0.0) rPot = std::max(rPot, c.p[3 * t + 1] + 12.0 * c.p[3 * t + 2]);
  // w(r) = 2 mu (V_C + U)/hbar^2 without the point-Coulomb 1/r part.
  auto w = [&](double r) {
    complex v(0.0, 0.0);
    if (sphere) v += r >= rc ? zz / r : zz * (3.0 - r * r / (rc * rc)) / (2.0 * rc);
    else if (r > 0.0) v += zz / r;
    if (ws) {  // a term with depth 0 is absent (its R, a may be 0)
      if (c.p[0] != 0.0) v -= c.p[0] * WoodsSaxon(r, c.p[1], c.p[2]);
      if (c.p[3] != 0.0) v -= complex(0.0, c.p[3] * WoodsSaxon(r, c.p[4], c.p[5]));
      if (c.p[6] != 0.0) v -= complex(0.0, c.p[6] * SurfaceShape(r, c.p[7], c.p[8]));
    }
    return fac * v;
  };
  const double rtp = (eta + std::sqrt(eta * eta + l * (l + 1.0))) / k;
  const double quarter = 0.5 * M_PI / k;
  double rMatch = std::max(std::max((nStore - 1) * step, rPot), rtp) + std::max(2.0 * quarter, 1.0);
  int m1 = (int)std::ceil(rMatch / step);
  int d = std::max(2, (int)std::lround(quarter / step));
  if (m1 - d < nStore) m1 = nStore + d;
  const int nInt = m1 + 1;
  std::vector<complex> v(nInt), f(nInt);
  for (int i = 1; i < nInt; i++) {
    double r = i * step;
    complex g = k * k - l * (l + 1.0) / (r * r) - w(r);
    f[i] = 1.0 + step * step * g / 12.0;
  }
  // Series start: u = r^{l+1}(1 + a1 r + a2 r^2 + a3 r^3), with the point
  // Coulomb 2 c1/r (c1 = eta k) and the finite rest w0 at the origin.
  const double c1 = (sphere || plane) ? 0.0 : eta * k;
  complex w0 = fac * complex(sphere ? zz * 1.5 / rc : 0.0, 0.0);
  if (ws) {
    if (c.p[0] != 0.0) w0 -= fac * c.p[0] * WoodsSaxon(0.0, c.p[1], c.p[2]);
    if (c.p[3] != 0.0) w0 -= fac * complex(0.0, c.p[3] * WoodsSaxon(0.0, c.p[4], c.p[5]));
    if (c.p[6] != 0.0) w0 -= fac * complex(0.0, c.p[6] * SurfaceShape(0.0, c.p[7], c.p[8]));
  }
  complex a1 = c1 / (l + 1.0);
  complex a2 = (2.0 * c1 * a1 + w0 - k * k) / (4.0 * l + 6.0);
  complex a3 = (2.0 * c1 * a2 + (w0 - k * k) * a1) / (6.0 * l + 12.0);
  auto series = [&](double r) { return 1.0 + r * (a1 + r * (a2 + r * a3)); };
  v[0] = 0.0;
  v[1] = series(step);
  v[2] = std::pow(2.0, l + 1.0) * series(2.0 * step);
  for (int i = 2; i + 1 < nInt; i++) {
    v[i + 1] = ((12.0 - 10.0 * f[i]) * v[i] - f[i - 1] * v[i - 1]) / f[i + 1];
    if (std::abs(v[i + 1]) > 1.0e150)
      for (int j = 0; j <= i + 1; j++) v[j] *= 1.0e-150;
  }
  // Match c1 F + c2 G at m1 and m2 = m1 - d.
  const int m2 = m1 - d;
  double F1, G1, F2, G2;
  CoulombFG(l, eta, k * m1 * step, F1, G1);
  CoulombFG(l, eta, k * m2 * step, F2, G2);
  double det = F1 * G2 - F2 * G1;
  if (!(std::fabs(det) > 0.0) || !std::isfinite(det)) return false;
  complex cf = (v[m1] * G2 - v[m2] * G1) / det;
  complex cg = (F1 * v[m2] - F2 * v[m1]) / det;
  complex norm = cf - complex(0.0, 1.0) * cg;
  if (!(std::abs(norm) > 0.0) || !std::isfinite(std::abs(norm))) return false;
  u.assign(nStore, complex(0.0, 0.0));
  for (int i = 0; i < nStore; i++) u[i] = v[i] / norm;
  if (T) *T = cg / norm;
  return true;
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
  p.x = SpectatorCos(p.ksf, &p.thetaCm, &p.angleClamped);
  const double kb = beta * aa.k;
  p.q = std::sqrt(std::max(0.0, p.ksf * p.ksf + kb * kb - 2.0 * p.ksf * kb * p.x));
  // Plane-wave limit: the Fourier transform of phi at q.
  double mpw = 0.0;
  for (int i = i0; i < n; i++) {
    double r = i * h;
    double qr = p.q * r;
    double j0 = qr < 1.0e-6 ? 1.0 - qr * qr / 6.0 : std::sin(qr) / qr;
    mpw += simpson[i] * r * r * j0 * phi[i];
  }
  p.mpw = 4.0 * M_PI * mpw;
  // Partial waves.
  Channel c = SfAt(energy);
  const int lCap = (int)uAA.size() - 1;
  std::vector<complex> sigma = CoulombPhases(c.eta, lCap);
  std::vector<complex> u;
  complex sum(0.0, 0.0);
  double pPrev = 1.0, pl = 1.0;  // Legendre recurrence
  for (int l = 0; l <= lCap; l++) {
    if (l == 1) {
      pPrev = 1.0;
      pl = p.x;
    } else if (l > 1) {
      double next = ((2.0 * l - 1.0) * p.x * pl - (l - 1.0) * pPrev) / l;
      pPrev = pl;
      pl = next;
    }
    if (!Wave(c, l, h, n, u)) {
      std::ostringstream why;
      why << "the s + F wave l = " << l << " at E = " << energy << " MeV could not be normalized";
      p.why = why.str();
      return p;
    }
    complex integral(0.0, 0.0);
    double umax = 0.0;
    for (int i = i0; i < n; i++) {
      integral += simpson[i] * phi[i] * u[i] * uAA[l][i];
      umax = std::max(umax, std::abs(u[i]));
    }
    sum += (2.0 * l + 1.0) * sigma[l] * sigmaAA[l] * pl * integral;
    if (l == 0) p.tail = std::abs(phi[n - 1] * u[n - 1] * uAA[0][n - 1]) / kappa;
    p.lmax = l;
    if (l >= 1 && l < lCap && tailAA[l + 1] * std::max(2.0, umax) < 1.0e-13 * std::abs(sum)) break;
  }
  p.m = 4.0 * M_PI / (p.ksf * kb) * sum;
  p.tail = std::abs(sum) > 0.0 ? p.tail / std::abs(sum) : HUGE_VAL;
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
  const double kb = beta * aa.k;
  // Plane-wave limits.
  std::vector<double> mpw(nk, 0.0);
  for (int k = 0; k < nk; k++) {
    double s = 0.0;
    for (int i = i0; i < n; i++) {
      double r = i * h;
      double qr = act[k].q * r;
      double j0 = qr < 1.0e-6 ? 1.0 - qr * qr / 6.0 : std::sin(qr) / qr;
      s += simpson[i] * r * r * j0 * phi[i];
    }
    mpw[k] = 4.0 * M_PI * s;
  }
  // Partial waves: the radial integrals once, P_l(x_k) per node.
  Channel c = SfAt(energy);
  const int lCap = (int)uAA.size() - 1;
  std::vector<complex> sigma = CoulombPhases(c.eta, lCap);
  std::vector<complex> u, sum(nk, complex(0.0, 0.0));
  std::vector<double> pPrev(nk, 1.0), pl(nk, 1.0);
  for (int l = 0; l <= lCap; l++) {
    for (int k = 0; k < nk; k++) {
      if (l == 1) {
        pPrev[k] = 1.0;
        pl[k] = act[k].x;
      } else if (l > 1) {
        double next = ((2.0 * l - 1.0) * act[k].x * pl[k] - (l - 1.0) * pPrev[k]) / l;
        pPrev[k] = pl[k];
        pl[k] = next;
      }
    }
    if (!Wave(c, l, h, n, u)) {
      std::ostringstream why;
      why << "the s + F wave l = " << l << " at E = " << energy << " MeV could not be normalized";
      p.why = why.str();
      return p;
    }
    complex integral(0.0, 0.0);
    double umax = 0.0;
    for (int i = i0; i < n; i++) {
      integral += simpson[i] * phi[i] * u[i] * uAA[l][i];
      umax = std::max(umax, std::abs(u[i]));
    }
    double smallest = HUGE_VAL;
    for (int k = 0; k < nk; k++) {
      sum[k] += (2.0 * l + 1.0) * sigma[l] * sigmaAA[l] * pl[k] * integral;
      smallest = std::min(smallest, std::abs(sum[k]));
    }
    if (l == 0) p.tail = std::abs(phi[n - 1] * u[n - 1] * uAA[0][n - 1]) / kappa;
    p.lmax = l;
    if (l >= 1 && l < lCap && tailAA[l + 1] * std::max(2.0, umax) < 1.0e-13 * smallest) break;
  }
  double wsum = 0.0, m2 = 0.0, mpw2 = 0.0, th = 0.0, xs = 0.0, qs = 0.0, smallest = HUGE_VAL;
  bool finite = true;
  for (int k = 0; k < nk; k++) {
    const complex m = 4.0 * M_PI / (p.ksf * kb) * sum[k];
    finite = finite && std::isfinite(std::abs(m)) && mpw[k] != 0.0;
    smallest = std::min(smallest, std::abs(sum[k]));
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
  p.tail = smallest > 0.0 ? p.tail / smallest : HUGE_VAL;
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
  if (!(k.bind > 0.0)) return "the Trojan horse is not bound: B(x+s) = " + Number(k.bind) + " MeV";
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
    c.mu = m1 * m2 / (m1 + m2) * uconv;
    if (kind != OPTICAL || o.kind != 3) return "";
    // A global potential: its projectile is the partner it is made for.
    const ThmGlobalOptical &g = ThmGlobalOpticals()[o.global];
    const std::string what = std::string(keys[which]) + "=" + g.name + " (" + g.reference + ", for " +
                             g.projectiles + "): the " + labels[which] + " channel is " + NucleusName(z1, a1) +
                             " + " + NucleusName(z2, a2);
    if (a1 <= 0 || a2 <= 0) return std::string(keys[which]) + "=" + g.name + ": the mass numbers are not known";
    const bool first = ThmGlobalOpticalFor(o.global, z1, a1);
    if (!first && !ThmGlobalOpticalFor(o.global, z2, a2))
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
    const ThmGlobalOptical &g = ThmGlobalOpticals()[c.global];
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
    const bool massOk = c.At >= g.aMin && c.At <= g.aMax;
    const bool energyOk = l0 >= g.eMin && l1 <= g.eMax;
    if (massOk && energyOk) continue;
    std::ostringstream why;
    why.precision(4);
    why << keys[which] << "=" << g.name << " (" << g.reference << ") is outside its validity range for "
        << NucleusName(c.Zp, c.Ap) << " + " << NucleusName(c.Zt, c.At) << ":";
    if (!massOk) why << " target A = " << c.At << " (valid " << g.aMin << "-" << g.aMax << ")";
    if (!energyOk) {
      why << (massOk ? "" : ",") << " lab energy of " << NucleusName(c.Zp, c.Ap) << " " << l0;
      if (which == 1 && l1 != l0) why << "-" << l1;
      why << " MeV over the data (valid " << g.eMin << "-" << g.eMax << " MeV)";
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
  h = std::min(0.02, std::min(0.1 / std::max(kSF, 1.0e-6), 0.1 / std::max(beta * kAA, 1.0e-6)));
  h = std::min(h, 0.2 / kappa);
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
    ThmGlobalOpticalEvaluate(c.global, c.Zp, c.Ap, c.Zt, c.At, elab[end], p[end]);
  }
  return true;
}

double ThmDistortion::Phi(double r) const {
  if (!(r > 0.0) || r < rmin - 1.0e-12) return 0.0;
  if (yukawa || etaB == 0.0) return std::exp(-kappa * r) / r;
  gsl_sf_result U;
  int status;
  {
    GslQuiet quiet;  // status checked here; thread-safe (GSLException.h)
    status = gsl_sf_hyperg_U_e(1.0 + etaB, 2.0, 2.0 * kappa * r, &U);
  }
  return status == GSL_SUCCESS || status == GSL_EUNDRFLW ? 2.0 * kappa * std::exp(-kappa * r) * U.val : 0.0;
}

std::string ThmDistortion::Build(const ThmExperiment &x, const Kinematics &k, double eLo, double eHi,
                                 double eRefDefault) {
  if (!(dataLo <= dataHi)) {
    dataLo = eLo;
    dataHi = eHi;
  }
  std::string setup = Setup(x, k, eLo);
  if (!setup.empty()) return setup;
  const double rEnd = rmin + 50.0 / kappa;
  int intervals = (int)std::ceil((rEnd - rmin) / h);
  if (intervals % 2) intervals++;
  n = i0 + intervals + 1;
  if (n > 400000) return "the radial grid would need " + Number(n) + " points (kappa too small)";
  phi.assign(n, 0.0);
  simpson.assign(n, 0.0);
  {
    GslQuiet quiet;  // statuses checked below; thread-safe (GSLException.h)
    for (int i = std::max(i0, 1); i < n; i++) {
      double r = i * h;
      if (yukawa || etaB == 0.0) {
        phi[i] = std::exp(-kappa * r) / r;
      } else {
        gsl_sf_result U;
        int status = gsl_sf_hyperg_U_e(1.0 + etaB, 2.0, 2.0 * kappa * r, &U);
        phi[i] = status == GSL_SUCCESS || status == GSL_EUNDRFLW ? 2.0 * kappa * std::exp(-kappa * r) * U.val : 0.0;
      }
    }
  }
  for (int i = i0; i < n; i++) {
    int j = i - i0;
    simpson[i] = h / 3.0 * ((i == i0 || i == n - 1) ? 1.0 : (j % 2 ? 4.0 : 2.0));
  }

  // The a + A waves, once.
  uAA.clear();
  std::vector<double> amp;
  double ampMax = 0.0;
  const int kLCap = 400;
  // The waves are kept, (l + 1) n complex numbers: bound them before they
  // grow (n <= 400000 and l <= 400 alone would allow 2.6 GB).
  const double kMaxWaveMB = 512.0;
  for (int l = 0;; l++) {
    if (l > kLCap) return "the partial-wave sum does not converge by l = " + Number(kLCap);
    if ((l + 1.0) * n * sizeof(complex) / 1048576.0 > kMaxWaveMB)
      return "the a + A partial waves would need more than " + Number(kMaxWaveMB) + " MB (" + Number(n) +
             " radial points, l = " + Number(l) + ")";
    std::vector<complex> u;
    if (!Wave(aa, l, beta * h, n, u)) return "the a + A wave l = " + Number(l) + " could not be normalized";
    double a = 0.0;
    for (int i = i0; i < n; i++) a += simpson[i] * std::fabs(phi[i]) * std::abs(u[i]);
    a *= 2.0 * l + 1.0;
    uAA.push_back(u);
    amp.push_back(a);
    ampMax = std::max(ampMax, a);
    if (l >= 2 && a < 1.0e-17 * ampMax && a < amp[l - 1]) break;
  }
  tailAA.assign(amp.size() + 1, 0.0);
  for (size_t l = amp.size(); l-- > 0;) tailAA[l] = tailAA[l + 1] + amp[l];
  sigmaAA = CoulombPhases(aa.eta, (int)uAA.size() - 1);

  SetAngleSlots(eHi);
  std::ostringstream d;
  d.precision(6);
  d << (kind == COULOMB ? "coulomb" : "optical") << "; a + A: " << ChannelText(0) << "; s + F: " << ChannelText(1)
    << (angWindow ? "; spectator directions " + AngleText() + " (" + Number(angNodes) + " nodes in cos theta_cm" +
                        (angSlots > 1 ? " per branch" : "") + (qCut ? ", |p_s| cut by the ps window" : "") +
                        "), acceptance-averaged |M|^2 and |M_PW|^2"
                  : "; spectator angle " +
                        (angleKind == QF ? std::string("qf (k_sF along k_aA)")
                                         : (angleKind == LAB ? "lab " : "cm ") + Number(angle) + " deg"))
    << "; bound state " << (yukawa ? "yukawa" : "whittaker") << (rmin > 0.0 ? ", r >= " + Number(rmin) + " fm" : "")
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
  // (EData refuses data points without one): those take the nearest value.
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
  pwSignChange = false;
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
  double s = t - j;  // nodes j-1, j, j+1, j+2 at s = -1, 0, 1, 2
  double v = -s * (s - 1.0) * (s - 2.0) / 6.0 * lnR[j - 1] + (s + 1.0) * (s - 1.0) * (s - 2.0) / 2.0 * lnR[j] -
             (s + 1.0) * s * (s - 2.0) / 2.0 * lnR[j + 1] + (s + 1.0) * s * (s - 1.0) / 6.0 * lnR[j + 2];
  return std::exp(v);
}
