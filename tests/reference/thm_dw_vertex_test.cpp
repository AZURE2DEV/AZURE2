/*!
 * Reference check of the distorted-wave THM entrance vertex (ThmDwVertex,
 * src/ThmDwVertex.cpp; docs/source/theory/thm_implementation.rst,
 * "Distorted-wave entrance vertex"): the surface term of the prior-form DWBA,
 *
 *   V_lm(B) = (B - 1) S_lm(a) - a S_lm'(a),
 *   S(r) = Int d^3u phi(u) chi(-)*_{k_sF}(alpha r + u) chi(+)_{k_aA}(r + beta u),
 *
 * divided by 4 pi phi~(q), summed over m: the Gram matrix G of
 * (S_lm(a), a S_lm'(a)) with |M_l|^2 = c^+ G c, c = (B - 1, -1).
 *
 *  (b) Plane waves in both channels: S = phi~(q) e^{i p.r}, so G must be
 *      [[j^2, j rho j'], [rho j' j, (rho j')^2]] at rho = p a, p = |k_aA -
 *      alpha k_sF| -- the plane-wave vertex M_l = (B - 1) j_l - rho j_l' --
 *      and the ratios |M_l|^2/|M_0|^2 for any B equal the plane-wave ones
 *      (1e-8).  19F(d,n) (Trojan horse as target, qf), 12C(14N,d) (horse as
 *      beam, 40 deg in the c.m.: every m_l enters), and a ps window: every
 *      direction whose q lies in the cut, Gauss-Legendre in cos(theta_cm)
 *      with the fixed-E weight |phi~(q)|^2 d cos(theta_cm) (against an
 *      independent quadrature of phi~), one node of averaged G.
 *  (c) Point-Coulomb a + A wave for 19F(d,n) at 55 MeV: G for l = 0 and 1
 *      against thm_dw_vertex_reference.py, a direct three-dimensional
 *      quadrature of S(r) near r = a with mpmath Coulomb functions, projected
 *      on Y_l0 -- nothing shared with the engine's reduced amplitudes,
 *      Numerov/COUL waves or Lagrange tables.
 *  Also: ThmCG against AngCoeff::ClebGord; the Cholesky factor reproduces
 *  c^+ G c; interpolation between grid nodes; G is the same to the last bit
 *  on 1 and 4 threads.
 *
 * Run:  tests/reference/thm_dw_vertex_test      (ctest: thm_dw_vertex)
 *       tests/reference/thm_dw_vertex_test -v   (print every value)
 */
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <string>
#include <vector>

#include "AngCoeff.h"
#include "Config.h"
#include "Constants.h"
#include "ThmDwVertex.h"
#include "ThmExperiment.h"
#include "ThmLineshape.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// Each consumer of the engine defines this itself (see thm_coulomb_term_test).
Config *g_config = nullptr;

namespace {

bool verbose = false;
int failures = 0;

void Check(const std::string &what, double got, double want, double tol, double scale = 0.0) {
  double rel = std::fabs(got - want) / std::max(std::max(std::fabs(want), scale), 1.0e-300);
  bool ok = rel <= tol;
  if (!ok) failures++;
  if (!ok || verbose)
    std::printf("  %s  %-52s got %.12e want %.12e rel %.2e (tol %.0e)\n", ok ? "ok  " : "FAIL", what.c_str(), got,
                want, rel, tol);
}

const double kAmu = 931.49410242;  // EData::BuildThmGroups
const double kN14 = 13.9992355671, kC12 = 11.9967096429, kD = 2.0135532134, kP = 1.0072764675,
             kN = 1.0086649159, kF19 = 18.9934652;

ThmDistortion::Kinematics F19() {  // 19F + d at 55 MeV, d = p + n, spectator n
  ThmDistortion::Kinematics k;
  k.Za = 1, k.ZA = 9, k.Zs = 0, k.Zx = 1;
  k.ma = kD, k.mA = kF19, k.ms = kN, k.mx = kP;
  k.horseIsBeam = false;
  k.mBeam = kF19, k.mTarget = kD, k.beamEnergy = 55.0;
  k.bind = (k.mx + k.ms - k.ma) * kAmu;
  return k;
}
ThmDistortion::Kinematics C12() {  // 14N + 12C at 30 MeV, 14N = 12C + d, spectator d
  ThmDistortion::Kinematics k;
  k.Za = 7, k.ZA = 6, k.Zs = 1, k.Zx = 6;
  k.ma = kN14, k.mA = kC12, k.ms = kD, k.mx = kC12;
  k.horseIsBeam = true;
  k.mBeam = kN14, k.mTarget = kC12, k.beamEnergy = 30.0;
  k.bind = (k.mx + k.ms - k.ma) * kAmu;
  return k;
}

ThmExperiment Plane() {
  ThmExperiment x;
  x.name = "plane";
  x.distortion = ThmExperiment::DIST_OPTICAL;
  x.opticalAA.kind = x.opticalSF.kind = 0;
  return x;
}

double JlPrime(int l, double x) {  // j_l'(x) = l j_l/x - j_{l+1}
  return l * gsl_sf_bessel_jl(l, x) / x - gsl_sf_bessel_jl(l + 1, x);
}

// The plane-wave Gram entries at rho.
void PlaneGram(int l, double rho, double g[4]) {
  double j = gsl_sf_bessel_jl(l, rho), d = rho * JlPrime(l, rho);
  g[0] = j * j;
  g[1] = d * d;
  g[2] = j * d;
  g[3] = 0.0;
}

// Checks one node's Gram entries against plane waves at momentum p (fm^-1).
void CheckPlaneNode(const std::string &tag, const ThmDwVertex &v, const double *g, double p) {
  const int nl = (int)v.lvals.size();
  double ref0[4];
  PlaneGram(v.lvals[0], p * v.radius, ref0);
  for (int li = 0; li < nl; li++) {
    int l = v.lvals[li];
    double ref[4];
    PlaneGram(l, p * v.radius, ref);
    const double *gl = g + li * 4;
    const double scale = std::max(ref[0], ref[1]);
    char label[128];
    std::snprintf(label, sizeof label, "%s l=%d G11", tag.c_str(), l);
    Check(label, gl[0], ref[0], 1.0e-9, scale);
    std::snprintf(label, sizeof label, "%s l=%d G22", tag.c_str(), l);
    Check(label, gl[1], ref[1], 1.0e-9, scale);
    std::snprintf(label, sizeof label, "%s l=%d Re G12", tag.c_str(), l);
    Check(label, gl[2], ref[2], 1.0e-9, scale);
    std::snprintf(label, sizeof label, "%s l=%d Im G12", tag.c_str(), l);
    Check(label, gl[3], 0.0, 1.0e-9, scale);
    // Vertex ratios for real and complex boundaries.
    for (complex B : {complex(-1.72, 0.0), complex(0.35, 0.0), complex(-1.2, 0.4)}) {
      double dw = ThmDwVertex::Vertex2(gl, B) / ThmDwVertex::Vertex2(g, B);
      double pw = ThmDwVertex::Vertex2(ref, B) / ThmDwVertex::Vertex2(ref0, B);
      std::snprintf(label, sizeof label, "%s |M_%d/M_%d|^2 B=(%g,%g)", tag.c_str(), l, v.lvals[0], B.real(),
                    B.imag());
      Check(label, dw, pw, 1.0e-8);
    }
  }
}

// k_sF at E (fm^-1).
double Ksf(const ThmDwVertex &v, double e) { return std::sqrt(2.0 * v.dist.sf.mu * v.dist.EsF(e)) / hbarc; }

}  // namespace

int main(int argc, char **argv) {
  verbose = argc > 1 && std::strcmp(argv[1], "-v") == 0;
  const bool probe = argc > 1 && std::strcmp(argv[1], "-p") == 0;

  std::printf("ThmCG against AngCoeff::ClebGord\n");
  for (int j1 = 0; j1 <= 9; j1 += 3)
    for (int j2 = 0; j2 <= 8; j2 += 2)
      for (int J = std::abs(j1 - j2); J <= j1 + j2; J += 3)
        for (int m1 = -j1; m1 <= j1; m1 += 2) {
          int m2 = std::max(-j2, std::min(j2, 1 - m1));
          if (std::abs(m1 + m2) > J) continue;
          char label[96];
          std::snprintf(label, sizeof label, "<%d %d %d %d|%d %d>", j1, m1, j2, m2, J, m1 + m2);
          Check(label, ThmCG(j1, m1, j2, m2, J, m1 + m2), AngCoeff::ClebGord(j1, j2, J, m1, m2, m1 + m2), 1.0e-12,
                1.0);
        }
  {
    // Large L, small J (the vertex needs J = l <= ~10): orthonormality
    // sum_M <L1 M L2 -M|J 0><L1 M L2 -M|J' 0> = delta_JJ'.
    const int L1 = 60, L2 = 63;
    for (int J = 3; J <= 11; J += 2)
      for (int Jp = J; Jp <= 11; Jp += 2) {
        double sum = 0.0;
        for (int M = -L1; M <= L1; M++) sum += ThmCG(L1, M, L2, -M, J, 0) * ThmCG(L1, M, L2, -M, Jp, 0);
        char label[96];
        std::snprintf(label, sizeof label, "sum_M <60 M 63 -M|%d 0><..|%d 0>", J, Jp);
        Check(label, sum, J == Jp ? 1.0 : 0.0, 1.0e-12, 1.0);
      }
  }

  if (probe) {
    // The values of the Python prototype (scratch); not a test.
    ThmExperiment x;
    x.name = "probe";
    x.distortion = ThmExperiment::DIST_COULOMB;
    ThmDwVertex v;
    std::string why = v.Build(x, F19(), 5.136, {0, 1, 3}, 0.15, 0.9, {});
    std::printf("built: '%s' %.2f s, L <= %d, %d; %d x %d nodes\n", why.c_str(), v.buildSeconds, v.laMax, v.lsMax,
                v.uNodes, v.cNodes);
    for (double e : {0.2127, 0.3239, 0.8282}) {
      std::vector<double> w, q, g, gd;
      double qd, pd;
      v.Interpolate(e, w, q, g, gd, qd, pd);
      for (int li = 0; li < 3; li++)
        std::printf("E %.4f l %d q %.6f p %.6f G11 %.10e G22 %.10e G12 %.10e %+.10e\n", e, v.lvals[li], qd, pd,
                    gd[li * 4], gd[li * 4 + 1], gd[li * 4 + 2], gd[li * 4 + 3]);
    }
    return 0;
  }

  std::printf("(b) plane waves in both channels: the plane-wave vertex at p = |k_aA - alpha k_sF|\n");
  {
    // 19F(d,n): the Trojan horse is the target, qf; l = 0..4 at grid nodes.
    ThmDwVertex v;
    std::string why = v.Build(Plane(), F19(), 5.136, {0, 1, 2, 3, 4}, -0.3, 0.9, {});
    if (!why.empty()) {
      std::printf("  FAIL  19F plane: %s\n", why.c_str());
      failures++;
    } else {
      for (int i : {0, 17, 31, 60}) {
        double e = v.gridLo + i * v.gridStep;
        const int nl = (int)v.lvals.size();
        const double *g = &v.Gd[(size_t)i * nl * 4];
        double ks = Ksf(v, e), ka = v.dist.aa.k;
        double p = std::fabs(ka - v.alpha * ks);  // qf: k_sF along k_aA
        char tag[64];
        std::snprintf(tag, sizeof tag, "19F qf E=%.3f", e);
        Check(std::string(tag) + " p", v.pd[i], p, 1.0e-13);
        CheckPlaneNode(tag, v, g, p);
      }
    }
  }
  {
    // 12C(14N,d): the horse is the beam, 40 deg in the c.m.: m_l != 0 enters.
    ThmExperiment x = Plane();
    x.angleKind = 2;
    x.angle = 40.0;
    ThmDwVertex v;
    std::string why = v.Build(x, C12(), 6.0, {0, 2, 4, 6, 8}, 1.0, 2.0, {});
    if (!why.empty()) {
      std::printf("  FAIL  12C plane: %s\n", why.c_str());
      failures++;
    } else {
      for (int i : {3, 30}) {
        double e = v.gridLo + i * v.gridStep;
        const int nl = (int)v.lvals.size();
        double ks = Ksf(v, e), ka = v.dist.aa.k;
        double xsa = std::cos(40.0 * M_PI / 180.0);  // horse = beam: x = cos theta_cm
        double p = std::sqrt(ka * ka + v.alpha * v.alpha * ks * ks - 2.0 * ka * v.alpha * ks * xsa);
        char tag[64];
        std::snprintf(tag, sizeof tag, "12C 40deg E=%.3f", e);
        Check(std::string(tag) + " p", v.pd[i], p, 1.0e-13);
        CheckPlaneNode(tag, v, &v.Gd[(size_t)i * nl * 4], p);
      }
    }
  }
  {
    // A ps window: every direction whose q = |k_sF - beta k_aA| lies in the
    // cut, Gauss-Legendre in cos(theta_cm) with the fixed-E event weight
    // d cos(theta_cm) |phi~(q)|^2; the averaged G is one node.
    ThmExperiment x = Plane();
    x.psKind = ThmExperiment::PS_HULTHEN;
    x.psMin = 5.0;
    x.psMax = 40.0;
    x.psNodes = 8;
    ThmDwVertex v;
    std::string why = v.Build(x, F19(), 5.136, {0, 1, 3}, 0.0, 0.8, {0.1, 0.5});
    if (!why.empty() || !v.angles || v.nAng != x.psNodes) {
      std::printf("  FAIL  19F window: %s (angles %d, %d directions)\n", why.c_str(), (int)v.angles, v.nAng);
      failures++;
    } else {
      const int nl = (int)v.lvals.size();
      // phi~(q) = 4 pi Int r^2 j_0(q r) phi(r) dr (Simpson, independent of the engine's u rule).
      auto phit = [&](double q) {
        const int n = 40000;
        const double rEnd = v.dist.rmin + 60.0 / v.dist.kappa, h = rEnd / n;
        double s = 0.0;
        for (int i = 0; i <= n; i++) {
          double r = i * h, qr = q * r, j0 = qr < 1.0e-6 ? 1.0 - qr * qr / 6.0 : std::sin(qr) / qr;
          s += (i == 0 || i == n ? 1.0 : i % 2 ? 4.0 : 2.0) * r * r * j0 * v.dist.Phi(r);
        }
        return 4.0 * M_PI * s * h / 3.0;
      };
      gsl_integration_glfixed_table *t = gsl_integration_glfixed_table_alloc(x.psNodes);
      for (int i : {2, 25}) {
        double e = v.gridLo + i * v.gridStep;
        double ks = Ksf(v, e), ka = v.dist.aa.k, kb = v.beta * ka;
        // The cut as a range of x = k^_sF . k^_aA.
        double ql = x.psMin / hbarc, qh = x.psMax / hbarc;
        double xlo = std::max(-1.0, (ks * ks + kb * kb - qh * qh) / (2.0 * ks * kb));
        double xhi = std::min(1.0, (ks * ks + kb * kb - ql * ql) / (2.0 * ks * kb));
        std::vector<std::pair<double, double>> want;  // (q MeV/c, weight)
        std::vector<double> pk;
        double total = 0.0;
        for (int k = 0; k < x.psNodes; k++) {
          double xk, wk;
          gsl_integration_glfixed_point(xlo, xhi, k, &xk, &wk, t);
          double q = std::sqrt(ks * ks + kb * kb - 2.0 * ks * kb * xk), f = phit(q);
          want.push_back({q * hbarc, wk * f * f});
          total += wk * f * f;
          pk.push_back(std::sqrt(ka * ka + v.alpha * v.alpha * ks * ks - 2.0 * ka * v.alpha * ks * xk));
        }
        std::vector<std::pair<double, double>> got;
        for (int k = 0; k < v.nAng; k++) got.push_back({v.aq[(size_t)i * v.nAng + k], v.aw[(size_t)i * v.nAng + k]});
        std::sort(want.begin(), want.end());
        std::sort(got.begin(), got.end());
        for (int k = 0; k < x.psNodes; k++) {
          char tag[64];
          std::snprintf(tag, sizeof tag, "19F window E=%.2f node %d", e, k);
          Check(std::string(tag) + " q", got[k].first, want[k].first, 1.0e-12);
          Check(std::string(tag) + " weight", got[k].second, want[k].second / total, 1.0e-9);
        }
        // The averaged Gram matrix = the weight-averaged plane-wave ones at
        // p = |k_aA - alpha k_sF| of each direction (the engine's weights).
        std::vector<double> avg((size_t)nl * 4, 0.0);
        for (int k = 0; k < v.nAng; k++) {
          const double q = v.aq[(size_t)i * v.nAng + k] / hbarc;
          const double xk = (ks * ks + kb * kb - q * q) / (2.0 * ks * kb);
          const double p = std::sqrt(ka * ka + v.alpha * v.alpha * ks * ks - 2.0 * ka * v.alpha * ks * xk);
          for (int li = 0; li < nl; li++) {
            double g[4];
            PlaneGram(v.lvals[li], p * v.radius, g);
            for (int c = 0; c < 4; c++) avg[li * 4 + c] += v.aw[(size_t)i * v.nAng + k] * g[c];
          }
        }
        for (int li = 0; li < nl; li++) {
          const double scale = std::max(avg[li * 4], avg[li * 4 + 1]);
          for (int c = 0; c < 4; c++) {
            char label[96];
            std::snprintf(label, sizeof label, "19F window E=%.2f l=%d <G>[%d]", e, v.lvals[li], c);
            Check(label, v.G[(size_t)i * nl * 4 + c + li * 4], avg[li * 4 + c], 1.0e-9, scale);
          }
        }
        // The spectator relation p^2/2mu_xA = E + B + q^2/2mu_sx holds up to
        // the mass defect of the Trojan horse, O(B/m_a c^2) = 1.2e-3 (the
        // identity needs m_a = m_s + m_x; nuclear masses are used).
        const double muxA = kP * kF19 / (kP + kF19) * uconv, musx = kN * kP / (kN + kP) * uconv;
        for (int k : {0, x.psNodes - 1}) {
          const double q = want[k].first / hbarc;
          const double xk = (ks * ks + kb * kb - q * q) / (2.0 * ks * kb);
          const double p = std::sqrt(ka * ka + v.alpha * v.alpha * ks * ks - 2.0 * ka * v.alpha * ks * xk);
          char tag[64];
          std::snprintf(tag, sizeof tag, "19F window E=%.2f node %d", e, k);
          Check(std::string(tag) + " p^2/2mu_xA - q^2/2mu_sx",
                p * p * hbarc * hbarc / (2.0 * muxA) - q * q * hbarc * hbarc / (2.0 * musx), e + v.dist.kin.bind,
                4.0e-3);
        }
      }
      gsl_integration_glfixed_table_free(t);
    }
  }

  std::printf("(c) point Coulomb in d + 19F: against the direct quadrature of thm_dw_vertex_reference.py\n");
  {
    // l, E, G11, G22, Re G12, Im G12 (resolution changes 3e-8, 1e-8).
    struct Ref {
      int l;
      double e, g11, g22, g12r, g12i;
    };
    const Ref refs[] = {
        {0, 0.3239, 3.151001686349e-01, 2.409660963146e-02, 7.359989426075e-03, 8.682555104375e-02},
        {1, 0.2127, 3.952661733153e-02, 3.246030552278e-02, 3.556264493470e-02, -4.283031649614e-03},
    };
    ThmExperiment x;
    x.name = "coulomb";
    x.distortion = ThmExperiment::DIST_COULOMB;
    for (const Ref &r : refs) {
      ThmDwVertex v;
      // The energy is the third grid node: no interpolation.
      std::string why = v.Build(x, F19(), 5.136, {r.l}, r.e - 0.04, r.e + 0.06, {});
      if (!why.empty()) {
        std::printf("  FAIL  19F coulomb: %s\n", why.c_str());
        failures++;
        continue;
      }
      const double *g = &v.Gd[2 * 4];
      const double scale = std::max(r.g11, r.g22);
      char label[96];
      std::snprintf(label, sizeof label, "l=%d E=%.4f G11", r.l, r.e);
      Check(label, g[0], r.g11, 1.0e-6, scale);
      std::snprintf(label, sizeof label, "l=%d E=%.4f G22", r.l, r.e);
      Check(label, g[1], r.g22, 1.0e-6, scale);
      std::snprintf(label, sizeof label, "l=%d E=%.4f Re G12", r.l, r.e);
      Check(label, g[2], r.g12r, 1.0e-6, scale);
      std::snprintf(label, sizeof label, "l=%d E=%.4f Im G12", r.l, r.e);
      Check(label, g[3], r.g12i, 1.0e-6, scale);
    }
  }

  std::printf("the Cholesky factor and the interpolation\n");
  {
    ThmExperiment x;
    x.name = "coulomb";
    x.distortion = ThmExperiment::DIST_COULOMB;
    ThmDwVertex v;
    std::string why = v.Build(x, F19(), 5.136, {0, 1}, 0.1, 0.5, {});
    if (!why.empty()) {
      std::printf("  FAIL  19F coulomb: %s\n", why.c_str());
      failures++;
    } else {
      ThmDwVertex::At at;
      for (double e : {0.2127, 0.3239}) {
        v.Evaluate(e, at);
        std::vector<double> w, q, g, gd;
        double qd, pd;
        v.Interpolate(e, w, q, g, gd, qd, pd);
        for (int li = 0; li < 2; li++)
          for (complex B : {complex(-1.7, 0.0), complex(0.8, -0.3)}) {
            complex m0 = at.a[li * 2] * (B - 1.0) - at.d[li * 2], m1 = at.a[li * 2 + 1] * (B - 1.0) - at.d[li * 2 + 1];
            char label[96];
            std::snprintf(label, sizeof label, "E=%.4f l=%d |M0|^2+|M1|^2 = c+Gc", e, v.lvals[li]);
            // At the qf angle only m = 0 enters and G has rank one; the cubic
            // interpolation of its entries breaks that at the 1e-9 level, and
            // the factor clamps the (tiny negative) second pivot to 0.
            Check(label, std::norm(m0) + std::norm(m1), ThmDwVertex::Vertex2(&g[li * 4], B), 1.0e-8);
          }
      }
      // Between nodes the cubic interpolation of G against a direct build on
      // a grid shifted by half a step.
      ThmDwVertex h;
      h.Build(x, F19(), 5.136, {0, 1}, 0.11, 0.5, {});
      for (int i : {4, 9}) {
        double e = h.gridLo + i * h.gridStep;
        std::vector<double> w, q, g, gd;
        double qd, pd;
        v.Interpolate(e, w, q, g, gd, qd, pd);
        for (int li = 0; li < 2; li++) {
          char label[96];
          std::snprintf(label, sizeof label, "E=%.3f l=%d G11 interpolated", e, v.lvals[li]);
          Check(label, gd[li * 4], h.Gd[(size_t)i * 2 * 4 + li * 4], 1.0e-6);
          std::snprintf(label, sizeof label, "E=%.3f l=%d G22 interpolated", e, v.lvals[li]);
          Check(label, gd[li * 4 + 1], h.Gd[(size_t)i * 2 * 4 + li * 4 + 1], 1.0e-6);
        }
      }
    }
  }

  std::printf("the same vertex, bit for bit, on 1 and 4 threads\n");
  {
    // The sum over the (u, cos theta) points is taken in a fixed order (fixed
    // blocks, added in block order), so the thread count and the schedule do
    // not move the last bits of G (they did with the arrival-order sum).
    ThmExperiment x;
    x.name = "coulomb";
    x.distortion = ThmExperiment::DIST_COULOMB;
    std::vector<double> g1, gd1;
    bool same = true;
    std::string why;
    for (int threads : {1, 4, 1, 4}) {
#ifdef _OPENMP
      omp_set_num_threads(threads);
#endif
      ThmDwVertex v;
      why = v.Build(x, F19(), 5.136, {0, 1, 2}, 0.1, 0.5, {});
      if (!why.empty()) break;
      if (g1.empty()) {
        g1 = v.G;
        gd1 = v.Gd;
      } else {
        same = same && v.G.size() == g1.size() && v.Gd.size() == gd1.size() &&
               std::memcmp(v.G.data(), g1.data(), g1.size() * sizeof(double)) == 0 &&
               std::memcmp(v.Gd.data(), gd1.data(), gd1.size() * sizeof(double)) == 0;
      }
    }
    const bool ok = why.empty() && same && !g1.empty();
    if (!ok) failures++;
    std::printf("  %s  G and Gd identical on 1, 4, 1, 4 threads (%zu values)%s%s\n", ok ? "ok  " : "FAIL", g1.size() + gd1.size(),
                why.empty() ? "" : ": ", why.c_str());
  }

  if (failures) {
    std::printf("%d check(s) failed\n", failures);
    return 1;
  }
  std::printf("all checks passed\n");
  return 0;
}
