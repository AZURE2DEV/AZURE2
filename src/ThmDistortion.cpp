#include "ThmDistortion.h"
#include "Config.h"
#include "ThmExperiment.h"
#include "cwfcomp_accurate.H"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <gsl/gsl_errno.h>
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

std::string ThmDistortion::Channel::Describe() const {
  std::ostringstream s;
  s.precision(6);
  if (kind == PLANE) return "plane wave";
  if (kind == POINT_COULOMB) {
    s << "point Coulomb (Z1 Z2 = " << Z1 * Z2 << ")";
    return s.str();
  }
  s << "Woods-Saxon V=" << p[0] << " R=" << p[1] << " a=" << p[2] << ", W=" << p[3] << " RW=" << p[4]
    << " aW=" << p[5] << ", WD=" << p[6] << " RD=" << p[7] << " aD=" << p[8] << ", Coulomb "
    << (p[9] > 0.0 ? "RC=" + Number(p[9]) : std::string("point")) << " (Z1 Z2 = " << Z1 * Z2 << ")";
  return s.str();
}

bool ThmDistortion::Wave(const Channel &c, int l, double step, int nStore, std::vector<complex> &u) const {
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

ThmDistortion::Point ThmDistortion::Evaluate(double energy) const {
  Point p;
  p.energy = energy;
  p.esf = EsF(energy);
  if (!(p.esf > 0.0)) {
    p.why = CheckEnergy(energy);
    return p;
  }
  p.ksf = std::sqrt(2.0 * sf.mu * p.esf) / hbarc;
  p.etasf = sf.kind == Channel::PLANE ? 0.0 : kin.Zs * (kin.Zx + kin.ZA) * fstruc * sf.mu / (hbarc * p.ksf);
  // Spectator direction in the c.m. (angle to the beam) and x = cos(k_sF, k_aA).
  double theta = 0.0;
  if (angleKind == QF) {
    theta = kin.horseIsBeam ? 0.0 : M_PI;
  } else if (angleKind == CM) {
    theta = angle * M_PI / 180.0;
  } else {
    double vs = hbarc * p.ksf / (kin.ms * uconv);
    double tl = angle * M_PI / 180.0;
    double s = vcm / vs * std::sin(tl);
    if (s > 1.0) {
      s = 1.0;
      p.angleClamped = true;
    }
    theta = tl + std::asin(s);
  }
  p.thetaCm = theta * 180.0 / M_PI;
  p.x = (kin.horseIsBeam ? 1.0 : -1.0) * std::cos(theta);
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
  Channel c = sf;
  c.k = p.ksf;
  c.eta = p.etasf;
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

std::string ThmDistortion::Build(const ThmExperiment &x, const Kinematics &k, double eLo, double eHi,
                                 double eRefDefault) {
  experiment = x.name;
  kin = k;
  kind = x.distortion == ThmExperiment::DIST_OPTICAL ? OPTICAL : COULOMB;
  angleKind = x.angleKind == 1 ? LAB : x.angleKind == 2 ? CM : QF;
  angle = x.angle;
  ratioPW = x.distortionRatioPW;
  yukawa = x.boundYukawa;
  if (!(k.bind > 0.0)) return "the Trojan horse is not bound: B(x+s) = " + Number(k.bind) + " MeV";
  eAA = k.beamEnergy * k.mTarget / (k.mBeam + k.mTarget);
  vcm = std::sqrt(2.0 * k.mBeam * uconv * k.beamEnergy) / ((k.mBeam + k.mTarget) * uconv);
  beta = k.ms / k.ma;
  muSx = k.ms * k.mx / (k.ms + k.mx) * uconv;
  kappa = std::sqrt(2.0 * muSx * k.bind) / hbarc;
  etaB = k.Zs * k.Zx * fstruc * muSx / (hbarc * kappa);
  auto channel = [&](const ThmExperiment::Optical &o, int z1, int z2, double m1, double m2, Channel &c) {
    c = Channel();
    c.kind = kind == COULOMB ? Channel::POINT_COULOMB
             : o.kind == 0   ? Channel::PLANE
             : o.kind == 1   ? Channel::POINT_COULOMB
                             : Channel::WOODS_SAXON;
    for (int t = 0; t < 10; t++) c.p[t] = o.p[t];
    c.Z1 = z1;
    c.Z2 = z2;
    c.mu = m1 * m2 / (m1 + m2) * uconv;
  };
  channel(x.opticalAA, k.Za, k.ZA, k.ma, k.mA, aa);
  channel(x.opticalSF, k.Zs, k.Zx + k.ZA, k.ms, k.mx + k.mA, sf);
  aa.k = std::sqrt(2.0 * aa.mu * eAA) / hbarc;
  aa.eta = aa.kind == Channel::PLANE ? 0.0 : aa.Z1 * aa.Z2 * fstruc * aa.mu / (hbarc * aa.k);

  // Radial step: h = 0.02 fm, or k_local h <= 0.1 in both channels.
  auto depth = [](const Channel &c) {
    return c.kind == Channel::WOODS_SAXON ? std::fabs(c.p[0]) + std::fabs(c.p[3]) + std::fabs(c.p[6]) : 0.0;
  };
  double kSF = std::sqrt(2.0 * sf.mu * (std::max(EsF(eLo), 0.0) + depth(sf))) / hbarc;
  double kAA = std::sqrt(2.0 * aa.mu * (eAA + depth(aa))) / hbarc;
  h = std::min(0.02, std::min(0.1 / std::max(kSF, 1.0e-6), 0.1 / std::max(beta * kAA, 1.0e-6)));
  h = std::min(h, 0.2 / kappa);
  rmin = x.boundRmin;
  i0 = 0;
  if (rmin > 0.0) {
    i0 = (int)std::ceil(rmin / h - 1.0e-9);
    h = rmin / i0;
  }
  const double rEnd = rmin + 50.0 / kappa;
  int intervals = (int)std::ceil((rEnd - rmin) / h);
  if (intervals % 2) intervals++;
  n = i0 + intervals + 1;
  if (n > 400000) return "the radial grid would need " + Number(n) + " points (kappa too small)";
  phi.assign(n, 0.0);
  simpson.assign(n, 0.0);
  gsl_error_handler_t *old = gsl_set_error_handler_off();
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
  gsl_set_error_handler(old);
  for (int i = i0; i < n; i++) {
    int j = i - i0;
    simpson[i] = h / 3.0 * ((i == i0 || i == n - 1) ? 1.0 : (j % 2 ? 4.0 : 2.0));
  }

  // The a + A waves, once.
  uAA.clear();
  std::vector<double> amp;
  double ampMax = 0.0;
  const int kLCap = 400;
  for (int l = 0;; l++) {
    if (l > kLCap) return "the partial-wave sum does not converge by l = " + Number(kLCap);
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

  std::ostringstream d;
  d.precision(6);
  d << (kind == COULOMB ? "coulomb" : "optical") << "; a + A: " << aa.Describe() << "; s + F: " << sf.Describe()
    << "; spectator angle "
    << (angleKind == QF ? std::string("qf (k_sF along k_aA)")
                        : (angleKind == LAB ? "lab " : "cm ") + Number(angle) + " deg")
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
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < nodes; i++) {
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
