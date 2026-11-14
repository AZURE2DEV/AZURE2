// The Legendre polynomials of a point: P_l(cos theta) at its c.m. angle, or, for a target effect with
// finite-geometry attenuation (Q) coefficients, their lab-frame average over the detector (Gauss-Legendre rule,
// lab boost of the detected particle), with their getters and setters.
#include <cmath>
#include <vector>

#include "CNuc.h"
#include "WorkingCopy.h"
#include "Unit.h"
#include "Kinematics.h"
#include "PPair.h"
#include "TargetEffect.h"

/// The maximum order of the Legendre polynomials stored in the point object.

int Site::GetMaxLOrder() const {
  return legendreP_ ? (int)legendreP_->size() - 1 : -1;
}

/*!
 * Returns the Legendre polynomial specified by an order.
 */

double Site::GetLegendreP(int lOrder) const {
  return (*legendreP_)[lOrder];
}

/*!
 * Adds a Legendre polynomial to the vector.  Polynomials are added to the vector in the order as L=0,1,2,...
 */

void Site::AddLegendreP(double polynomial) {
  MutableLegendre().push_back(polynomial);
}

/*!
 * Clears all Legendre polynomials stored in the point.
 * This is needed when recalculating polynomials after energy shifts.
 */

void Site::ClearLegendrePolynomials() {
  // An empty list needs no copy of a shared one.
  if (legendreP_ && !legendreP_.SoleOwner())
    legendreP_.reset();
  else if (legendreP_)
    MutableLegendre().clear();
}

/*!
 * Gauss-Legendre nodes and weights on [-1, 1] (Newton iteration on P_n), built once.
 */
namespace {
struct GaussLegendreRule {
  std::vector<double> x, w;
  explicit GaussLegendreRule(int n) : x(n), w(n) {
    for (int i = 0; i < (n + 1) / 2; i++) {
      double z = cos(pi * (i + 0.75) / (n + 0.5)), dp = 1.0;
      for (int it = 0; it < 100; it++) {
        double p0 = 1.0, p1 = z;
        for (int k = 2; k <= n; k++) {
          const double p2 = ((2.0 * k - 1.0) * z * p1 - (k - 1.0) * p0) / k;
          p0 = p1;
          p1 = p2;
        }
        dp = n * (z * p1 - p0) / (z * z - 1.0);
        const double dz = p1 / dp;
        z -= dz;
        if (fabs(dz) < 1e-15) break;
      }
      x[i] = -z;
      x[n - 1 - i] = z;
      w[i] = w[n - 1 - i] = 2.0 / ((1.0 - z * z) * dp * dp);
    }
  }
};

const GaussLegendreRule &QTransformRule() {
  static const GaussLegendreRule rule(96);
  return rule;
}

/// P_0(x) .. P_n(x) into p (size n + 1).
// Kept out of line: inlined into CalcLegendreP (GCC 11, -O3), the T2 calculation at 1 thread ran 3.6 % more
// instructions and 3.8 % more cycles than with this call; callgrind put all of the difference in CalcLegendreP.
__attribute__((noinline)) void LegendreSeries(int n, double x, double *p) {
  p[0] = 1.0;
  if (n >= 1) p[1] = x;
  for (int l = 2; l <= n; l++) p[l] = (2.0 * l - 1.0) / l * x * p[l - 1] - (l - 1.0) / l * p[l - 2];  // as before 2026-10-25
}
}  // namespace

/*!
 * Lab-frame motion of the detected particle (particle 1 of the exit pair) for the finite-geometry attenuation
 * coefficients: with cos theta* the c.m. angle, cos theta_lab = gb (x + g) / sqrt(gb^2 (x + g)^2 + 1 - x^2) and
 * dOmega_lab/dOmega_cm = gb (1 + g x) / (gb^2 (g + x)^2 + 1 - x^2)^(3/2).  Non-relativistic: gb = 1 and
 * g^2 = m1 m3 E / (m2 m4 (E + Q)) (E the c.m. energy; the same g as ConvertLabAngle's a13/a24); relativistic=1:
 * the c.m. boost gamma and beta / beta_3* of RelativisticLabToCM.  Returns false where the transform is not used.
 */
bool Site::QLabBoostAt(double ti, const PPair *entrancePair, const PPair *exitPair, double &g, double &gb) const {
  const double tf = ti + entrancePair->GetSepE() + entrancePair->GetExE() - exitPair->GetSepE() - exitPair->GetExE();
  if (ti <= 0.0 || tf <= 0.0) return false;
  if (entrancePair->FunctionSettings().relativistic) {
    const double m1 = entrancePair->GetM(1) * uconv, m2 = entrancePair->GetM(2) * uconv;
    const double m3 = exitPair->GetM(1) * uconv, m4 = exitPair->GetM(2) * uconv;
    const double e1 = Kinematics::LabFromCM(m1, m2, ti);
    const double sqrtS = m1 + m2 + ti, eTot = e1 + m1 + m2;
    const double beta = sqrt(e1 * (e1 + 2.0 * m1)) / eTot;
    const double pf = Kinematics::K(m3, m4, tf) * hbarc, e3 = sqrt(m3 * m3 + pf * pf);
    g = beta / (pf / e3);
    gb = eTot / sqrtS;
  } else {
    g = sqrt(entrancePair->GetM(1) * exitPair->GetM(1) * ti / (entrancePair->GetM(2) * exitPair->GetM(2) * tf));
    gb = 1.0;
  }
  return std::isfinite(g) && g < 1.0;
}

/*!
 * Calculates Legendre polynomials up to a maximum order.  The polynomials are added, in order, to a vector.
 *
 * With attenuation coefficients Q_K (finite detector geometry, M. E. Rose, Phys. Rev. 91, 610 (1953)) the stored
 * factors become the coefficients of the R-matrix's c.m. Legendre coefficients c_L in the measured quantity:
 * - photon exits, UPOS (the detected radiation is the secondary gamma ray), analyzing powers, non-particle entrances
 *   and charged-particle elastic scattering: P_L(x) Q_L, Q as entered (for elastic scattering the Coulomb amplitude is
 *   not a finite Legendre series, so the lab transform below does not apply; Q on such data stays an approximation);
 * - a particle exit otherwise: the detector averages the LAB angular distribution, so Q_K multiplies its LAB Legendre
 *   coefficients b_K = sum_L M_KL c_L, M_KL = (2K+1)/2 int dx P_K(x_lab(x)) P_L(x) (x the c.m. cosine; the lab and c.m.
 *   solid-angle elements cancel).  Back at the point's c.m. angle, with J = dOmega_lab/dOmega_cm,
 *     factor_L = P_L(x) + J sum_{K < nQ} (Q_K - 1) M_KL P_K(x_lab),
 *   exact for Q_K = 1 beyond the coefficients given; Q_K = 1 for every K gives P_L(x).
 * Before 2026-10-25 the particle-exit case multiplied P_L by B = U Q, U an O(gamma^2) lab -> c.m. coefficient matrix
 * evaluated at gamma^2 instead of gamma: B is not the attenuation of the lab distribution (it changed the result even
 * with every Q_K = 1), a non-particle entrance gave every factor 0, and more than maxL coefficients read past the end.
 */

void Site::CalcLegendreP(WorkingCopy *data, int maxL, const CNuc *theCNuc, const TargetEffect *targetEffect) {
  if (theCNuc) {
    const int entrance = theCNuc->GetPairNumFromKey(this->GetEntranceKey());
    if (entrance > 0)
      entrance_symmetry_ = LegendreSymmetryFactor(*theCNuc->GetPair(entrance), GetEntranceKey() == GetExitKey());
  }
  // A slot holding Legendre coefficients has no angle: its readers' rows apply them (ValueFor).
  if (is_sub_point_ && AngleFree()) return;
  if (maxL >= 0) {
    std::vector<double> p(maxL + 1);
    this->LegendreRowAt(this->GetCMEnergy(), maxL, theCNuc, targetEffect, p.data());
    for (int l = 0; l <= maxL; l++) this->AddLegendreP(p[l]);
  }
}

bool Site::RowDependsOnEnergy(const CNuc *theCNuc, const TargetEffect *targetEffect) const {
  if (!targetEffect || targetEffect->NumQCoefficients() <= 0 || !theCNuc) return false;
  const PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetEntranceKey()));
  const PPair *exitPair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetExitKey()));
  const bool elastic = this->GetEntranceKey() == this->GetExitKey();
  const bool charged = entrancePair->GetZ(1) != 0 && entrancePair->GetZ(2) != 0;
  return entrancePair->GetPType() == 0 && exitPair->GetPType() == 0 && !this->IsUPOS() && !this->IsAnalyzingPower() &&
         !(elastic && charged);
}

void Site::LegendreRowAt(double cmEnergy, int maxL, const CNuc *theCNuc, const TargetEffect *targetEffect, double *p) const {
  this->LegendreRowAt(cmEnergy, this->GetCMAngle(), maxL, theCNuc, targetEffect, p);
}

void Site::LegendreRowAt(double cmEnergy, double cmAngle, int maxL, const CNuc *theCNuc, const TargetEffect *targetEffect,
                         double *p) const {
  const double x = cos(cmAngle * pi / 180.0);
  LegendreSeries(maxL, x, p);
  const int numQ = targetEffect ? targetEffect->NumQCoefficients() : 0;
  if (numQ > 0) {
    const PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetEntranceKey()));
    const PPair *exitPair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetExitKey()));
    double g = 0.0, gb = 1.0;
    const bool labTransform =
        this->RowDependsOnEnergy(theCNuc, targetEffect) && this->QLabBoostAt(cmEnergy, entrancePair, exitPair, g, gb);
    if (!labTransform) {
      for (int l = 0; l <= maxL && l < numQ; l++) p[l] *= targetEffect->GetQCoefficient(l);
    } else {
      const GaussLegendreRule &rule = QTransformRule();
      const int n = static_cast<int>(rule.x.size());
      auto toLab = [g, gb](double c) { return gb * (c + g) / sqrt(gb * gb * (c + g) * (c + g) + 1.0 - c * c); };
      const double jac = gb * (1.0 + g * x) / pow(gb * gb * (g + x) * (g + x) + 1.0 - x * x, 1.5);
      std::vector<double> pc((maxL + 1) * n), pl(numQ * n), pK(numQ);
      for (int i = 0; i < n; i++) {
        LegendreSeries(maxL, rule.x[i], &pc[i * (maxL + 1)]);
        LegendreSeries(numQ - 1, toLab(rule.x[i]), &pl[i * numQ]);
      }
      LegendreSeries(numQ - 1, toLab(x), pK.data());
      std::vector<double> add(maxL + 1, 0.0);
      for (int K = 0; K < numQ; K++) {
        const double dq = targetEffect->GetQCoefficient(K) - 1.0;
        if (dq == 0.0) continue;
        const double f = jac * dq * pK[K] * (2.0 * K + 1.0) / 2.0;
        for (int i = 0; i < n; i++) {
          const double wk = f * rule.w[i] * pl[i * numQ + K];
          for (int l = 0; l <= maxL; l++) add[l] += wk * pc[i * (maxL + 1) + l];
        }
      }
      for (int l = 0; l <= maxL; l++) p[l] += add[l];
    }
  }
}
