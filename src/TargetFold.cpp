// The apparatus folds of a target-effect point: the integration of its sub-point values over the target, the beam
// spread, the UDR and beam-profile kernels and the free-gas Doppler broadening (one quadrature, GridQuadrature),
// the sigma(E_bar) fallback outside the windows of a thin target, and the linear weights of the fold for the
// analytic gradient (TargetEffectWeights).
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <utility>
#include <vector>

#include "Config.h"
#include "LayerExit.h"
#include "WorkingCopy.h"
#include "Unit.h"
#include "TELinearForm.h"
#include "TargetEffect.h"

namespace {
// 2-point Gauss-Legendre over the sub-point grid of `point` (energies descending from sub-point 1; units of `data`),
// restricted to the c.m. window [lo, hi]: the cross section (through `rd`) and, for a yield, the stopping power are
// interpolated linearly between the two sub-points that bracket each node.
// Accumulates num += integral of wNum(E) sigma(E) [/ (stopping(E) 1e24)] and den += integral of
// wDen(E), so a kernel average is num / den and a yield is num alone: the one quadrature of every effect.
template <class V, class Read, class WNum, class WDen>
void GridQuadrature(const Unit *point, const WorkingCopy &data, Read rd, double lo, double hi, bool divideByStopping,
                    WNum wNum, WDen wDen, V &num, double &den) {
  const double x1 = -1.0 / sqrt(3.0), x2 = 1.0 / sqrt(3.0);
  const int n = point->NumSubPoints();
  for (int j = 1; j < n; j++) {
    const double Ea = point->GetSubPoint(data, j)->GetCMEnergy();
    const double Eb = point->GetSubPoint(data, j + 1)->GetCMEnergy();
    if (Eb >= hi || Ea <= lo) continue;
    const double a = std::min(Ea, hi), b = std::max(Eb, lo);
    const double mid = 0.5 * (a + b), half = 0.5 * (a - b);
    if (half < 1.0e-12) continue;
    for (int gp = 0; gp < 2; gp++) {
      const double e = mid + half * (gp == 0 ? x1 : x2);
      const double frac = (Ea - e) / (Ea - Eb);   // 0 at sub-point j, 1 at sub-point j + 1
      V sigma = (1.0 - frac) * rd(j) + frac * rd(j + 1);
      double f = 1.0;
      if (divideByStopping) {
        double stopping = (1.0 - frac) * point->GetSubPoint(data, j)->GetStoppingPower() + frac * point->GetSubPoint(data, j + 1)->GetStoppingPower();
        f = 1.0 / stopping / 1.0e24;
      }
      num += (half * wNum(e) * f) * sigma;
      den += half * wDen(e);
    }
  }
}
}  // namespace

/*!
 * Integrates the sub-points of this point into its fitted value for the experimental effect it carries.
 *
 * Gaussian beam spread (a constant sigma or an expression, both through TargetEffect::SigmaCM):
 *   Y(E0) = int K(E'; E0) sigma(E') dE' / int K(E'; E0) dE',  K a Gaussian of width sigma_b(E0) about E0,
 * normalised by its own integral over the sub-point window, so a constant cross section is reproduced
 * exactly however the window is truncated (before 2026-10 the constant-sigma form was not normalised and
 * came out 0.27 % low, and the expression form took its width at E' rather than at the point's energy).
 *
 * Target integration:  Y(E0) = (1/N) int_{E0-T}^{E0} <sigma>(E_d) / epsilon(E_d) dE_d  over the depth, where
 * <sigma>(E_d) is sigma itself, or its Gaussian average about E_d of width sqrt(sigma_b(E_d)^2 +
 * sigma_straggling(E0 - E_d)^2) when there is a beam spread and/or straggling, normalised as above (the quadrature
 * averages sigma/epsilon over that Gaussian, which is <sigma>/epsilon(E_d) to the curvature of 1/epsilon across it).
 * Energies and epsilon are c.m. (lab epsilon times the point's c.m./lab factor): dE_cm/eps_cm = dE_lab/eps_lab.
 * An infinitely thick target (density 0) and a yield on a layer (TargetEffect::IsLayerYield) are not divided by N:
 *   Y(E0) = 1e-24 int_{E_out}^{E0} <sigma>(E_d) / epsilon(E_d) dE_d   reactions per incident particle (sigma in b,
 * epsilon in MeV cm^2 per active nucleus), with E_out = 1 keV for the thick target and, for the layer, the energy at
 * which the integral of 1/epsilon on the same nodes reaches the working copy's N (LayerExit.h).
 *
 * UDR (time-of-flight) and beam-profile kernels follow; both are normalised averages as well.
 * 2-point Gauss-Legendre on every sub-point interval throughout (GridQuadrature), exact for cubics and
 * indifferent to the adaptive grid's spacing.
 */

template <class V, class Read>
V Unit::IntegrateTargetEffectImpl(const WorkingCopy &data, const Config &configure, Read rd) {
  const TargetEffect *targetEffect = data.Data().GetTargetEffect(this->GetTargetEffectNum());
  if (!targetEffect->IsDoppler() || doppler_applied_) return this->IntegrateKernelsImpl<V>(data, configure, rd);
  // Free-gas Doppler broadening first (SAMMY's order), then the line's kernel on the broadened values.
  const int n = this->NumSubPoints();
  if (n < 2) return V(0.0);
  std::vector<std::pair<int, double>> w;
  auto broaden = [&](double e) {
    this->DopplerWeights(data, e, w);
    V s = 0.0;
    for (const auto &t : w) s += t.second * rd(t.first);
    return s;
  };
  if (!targetEffect->HasResolutionStage()) return broaden(this->GetCMEnergy());
  std::vector<V> b(n);
  for (int j = 1; j <= n; j++) b[j - 1] = broaden(this->GetSubPoint(data, j)->GetCMEnergy());
  return this->IntegrateKernelsImpl<V>(data, configure, [&b](int j) { return b[j - 1]; });
}

void Unit::DopplerWeights(const WorkingCopy &data, double eCM, std::vector<std::pair<int, double>> &w,
                            std::vector<std::pair<int, double>> *dwdT) const {
  w.clear();
  if (dwdT) dwdT->clear();
  const TargetEffect *te = data.Data().GetTargetEffect(this->GetTargetEffectNum());
  const int n = this->NumSubPoints();
  const double c = cm_factor_;
  // the present temperature (K) of the working copy, a fit parameter or the line's
  const double T = data.EffectValues(this->GetTargetEffectNum()).dopplerTemperature;
  const double u = te->DopplerUAt(c, T);
  if (n < 2 || !(u > 0.0) || !(eCM > 0.0)) return;
  // The window (layout temperature) and the pieces (nominal temperature) do not move with T: only the kernel
  // values below do, so sigma_D is smooth in a fitted T.  For a fixed line both are the line's T.
  const double u0 = te->DopplerUAt(c, te->GetDopplerNominal());
  double lo, hi;
  te->DopplerWindow(eCM, c, lo, hi);
  const double v = std::sqrt(eCM / c);
  const double norm = 1.0 / (u * std::sqrt(M_PI) * v * v);
  // 4-point Gauss-Legendre on [-1, 1]
  static const double gx[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526};
  static const double gw[4] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538};
  auto energy = [this, &data](int j) { return this->GetSubPoint(data, j)->GetCMEnergy(); };
  // sub-points run from high to low energy: the first interval (j, j+1) whose lower end lies below hi
  int first = 1, last = n - 1;
  {
    int a = 1, b = n;  // energy(a) >= ... >= energy(b)
    while (b - a > 1) {
      const int m = (a + b) / 2;
      (energy(m) > hi ? a : b) = m;
    }
    first = a;
  }
  std::map<int, double> acc, dacc;
  // d(norm k)/dT = norm [-k/u + dk/du] du/dT, du/dT = u/(2T)
  const double dudT = (T > 0.0) ? 0.5 * u / T : 0.0;
  const double u2 = u * u, u3 = u2 * u;
  for (int j = first; j <= last; j++) {
    const double Ea = energy(j), Eb = energy(j + 1);
    if (Ea <= lo) break;
    if (Eb >= hi) continue;
    const double a = std::min(Ea, hi), b = std::max(Eb, lo);
    if (!(a > b) || !(Ea > Eb)) continue;
    const double wa = std::sqrt(a / c), wb = std::sqrt(b / c);
    const int pieces = std::max(1, static_cast<int>(std::ceil((wa - wb) / (0.25 * u0))));
    const double hw = (wa - wb) / pieces;
    double toA = 0.0, toB = 0.0, dToA = 0.0, dToB = 0.0;
    for (int p = 0; p < pieces; p++) {
      const double mid = wb + (p + 0.5) * hw, half = 0.5 * hw;
      for (int g = 0; g < 4; g++) {
        const double x = mid + half * gx[g];
        const double k = (std::exp(-(v - x) * (v - x) / (u * u)) - std::exp(-(v + x) * (v + x) / (u * u))) * x * x;
        const double frac = (Ea - c * x * x) / (Ea - Eb);  // 0 at sub-point j, 1 at sub-point j + 1
        toA += half * gw[g] * k * (1.0 - frac);
        toB += half * gw[g] * k * frac;
        if (dwdT) {
          const double dm2 = (v - x) * (v - x), dp2 = (v + x) * (v + x);
          const double dk = -k / u + x * x * 2.0 * (dm2 * std::exp(-dm2 / u2) - dp2 * std::exp(-dp2 / u2)) / u3;
          dToA += half * gw[g] * dk * (1.0 - frac);
          dToB += half * gw[g] * dk * frac;
        }
      }
    }
    acc[j] += norm * toA;
    acc[j + 1] += norm * toB;
    if (dwdT) {
      dacc[j] += norm * dToA * dudT;
      dacc[j + 1] += norm * dToB * dudT;
    }
  }
  w.assign(acc.begin(), acc.end());
  if (dwdT) dwdT->assign(dacc.begin(), dacc.end());
}

double Unit::DopplerTemperatureDerivative(const WorkingCopy &data, const Config &configure,
                                            const std::vector<double> &values) {
  const TargetEffect *targetEffect = data.Data().GetTargetEffect(this->GetTargetEffectNum());
  if (!targetEffect->IsDoppler() || doppler_applied_) return 0.0;
  const int n = this->NumSubPoints();
  if (n < 2 || (int)values.size() < n) return 0.0;
  std::vector<std::pair<int, double>> w, dw;
  auto dBroaden = [&](double e) {
    this->DopplerWeights(data, e, w, &dw);
    double s = 0.0;
    for (const auto &t : dw) s += t.second * values[t.first - 1];
    return s;
  };
  if (!targetEffect->HasResolutionStage()) return dBroaden(this->GetCMEnergy());
  // The kernel after the Doppler step is linear in its sub-point values: apply it to d(sigma_D)/dT.
  std::vector<double> b(n);
  for (int j = 1; j <= n; j++) b[j - 1] = dBroaden(this->GetSubPoint(data, j)->GetCMEnergy());
  return this->IntegrateKernelsImpl<double>(data, configure, [&b](int j) { return b[j - 1]; });
}

template <class V, class Read>
V Unit::TargetFallbackImpl(const WorkingCopy &data, Read rd) {
  const int n = this->NumSubPoints();
  if (n <= 0) return V(0.0);
  if (n == 1) return rd(1);
  const double e = this->GetCMEnergy() - target_fallback_depth_;
  if (target_fallback_only_ && n == 3) {
    // three fallback nodes spanning the E_bar a fitted density can reach: the quadratic through them
    double x[3];
    for (int k = 0; k < 3; k++) x[k] = this->GetSubPoint(data, k + 1)->GetCMEnergy();
    const double lo = std::min(x[0], x[2]), hi = std::max(x[0], x[2]);
    const double ec = std::min(std::max(e, lo), hi);
    V v = 0.0;
    for (int k = 0; k < 3; k++) {
      double l = 1.0;
      for (int m = 0; m < 3; m++)
        if (m != k) l *= (ec - x[m]) / (x[k] - x[m]);
      v += l * rd(k + 1);
    }
    return v;
  }
  // a full grid (descending energies): linear between the two nodes around E_bar, the end node outside it
  if (e >= this->GetSubPoint(data, 1)->GetCMEnergy()) return rd(1);
  if (e <= this->GetSubPoint(data, n)->GetCMEnergy()) return rd(n);
  for (int j = 1; j < n; j++) {
    const double Ea = this->GetSubPoint(data, j)->GetCMEnergy(), Eb = this->GetSubPoint(data, j + 1)->GetCMEnergy();
    if (e <= Ea && e >= Eb) {
      const double frac = (Ea > Eb) ? (Ea - e) / (Ea - Eb) : 0.0;
      return (1.0 - frac) * rd(j) + frac * rd(j + 1);
    }
  }
  return rd(n);
}

template <class V, class Read>
V Unit::IntegrateKernelsImpl(const WorkingCopy &data, const Config &configure, Read rd) {
  V yield = 0.0;
  const TargetEffect *targetEffect = data.Data().GetTargetEffect(this->GetTargetEffectNum());
  const int numPoints = this->NumSubPoints();
  // A thin-target point outside the windows of an energy-restricted integration, or one the automatic
  // tolerance let off: sigma at the depth-averaged energy E_bar (SetTargetFallbackOnly, TargetFallbackScope).
  const bool thinTarget = targetEffect->IsTargetIntegration() && !targetEffect->IsThickTarget();
  if (thinTarget && (target_fallback_only_ || target_fallback_now_)) return this->TargetFallbackImpl<V>(data, rd);
  if (numPoints < 2) return V(0.0);
  const double eHigh = this->GetSubPoint(data, 1)->GetCMEnergy();
  const double eLow = this->GetSubPoint(data, numPoints)->GetCMEnergy();
  const double centroid = this->GetCMEnergy();
  const double range = targetEffect->convolutionRange;
  auto one = [](double) { return 1.0; };
  // sigma(E) / (epsilon(E) 1e24) interpolated at one energy, for a depth node with no spread to average over
  auto yieldAt = [&](double e) {
    for (int j = 1; j < numPoints; j++) {
      double Ea = this->GetSubPoint(data, j)->GetCMEnergy(), Eb = this->GetSubPoint(data, j + 1)->GetCMEnergy();
      if (e <= Ea && e >= Eb) {
        double frac = (Ea - e) / (Ea - Eb);
        V sigma = (1.0 - frac) * rd(j) + frac * rd(j + 1);
        double stopping = (1.0 - frac) * this->GetSubPoint(data, j)->GetStoppingPower() + frac * this->GetSubPoint(data, j + 1)->GetStoppingPower();
        return (1.0 / stopping / 1.0e24) * sigma;
      }
    }
    return V(0.0);
  };
  // Areal density (atoms/cm^2) of target traversed from the grid's top node down to c.m. energy e, int dE / epsilon
  // over the sub-points (trapezoid in 1/epsilon; epsilon in c.m. units as the sub-points carry it, and dE_cm /
  // eps_cm = dE_lab / eps_lab).  For Bohr straggling.  A node with a non-positive stopping power contributes nothing.
  std::vector<double> arealCum;
  auto invEps = [](double eps) { return (eps > 0.0) ? 1.0 / eps : 0.0; };
  if (targetEffect->IsBohrStraggling()) {
    arealCum.assign(numPoints + 1, 0.0);
    for (int j = 1; j < numPoints; j++) {
      double Ea = this->GetSubPoint(data, j)->GetCMEnergy(), Eb = this->GetSubPoint(data, j + 1)->GetCMEnergy();
      double ia = invEps(this->GetSubPoint(data, j)->GetStoppingPower()), ib = invEps(this->GetSubPoint(data, j + 1)->GetStoppingPower());
      arealCum[j + 1] = arealCum[j] + 0.5 * (ia + ib) * (Ea - Eb);
    }
  }
  auto arealTo = [&](double e) {
    if (arealCum.empty() || e >= eHigh) return 0.0;
    for (int j = 1; j < numPoints; j++) {
      double Ea = this->GetSubPoint(data, j)->GetCMEnergy(), Eb = this->GetSubPoint(data, j + 1)->GetCMEnergy();
      if (e <= Ea && e >= Eb) {
        double frac = (Ea - e) / (Ea - Eb);
        double ia = invEps(this->GetSubPoint(data, j)->GetStoppingPower()), ib = invEps(this->GetSubPoint(data, j + 1)->GetStoppingPower());
        double ie = (1.0 - frac) * ia + frac * ib;
        return arealCum[j] + 0.5 * (ia + ie) * (Ea - e);
      }
    }
    return arealCum[numPoints];
  };

  if (targetEffect->IsTargetIntegration()) {
    const double surface = centroid;
    double back = centroid - this->GetTargetThickness();
    V outer = 0.0;
    const bool noSpread = !targetEffect->IsGaussian() && !targetEffect->IsStraggling();
    // A yield on a layer (reactions per incident particle, docs/upgrade/energy_grid.md section 5.2): the depth
    // integral runs from the top of the integration down to E_out, found on these nodes with this quadrature
    // (LayerExit.h), so that int_{E_out}^{top} dE/epsilon is the areal density the working copy holds.
    if (targetEffect->IsLayerYield()) {
      const double top = noSpread ? eHigh : surface;
      auto energy = [&](int j) { return this->GetSubPoint(data, j)->GetCMEnergy(); };
      auto stopping = [&](int j) { return this->GetSubPoint(data, j)->GetStoppingPower(); };
      if (!LayerExitEnergy(numPoints, energy, stopping, top, data.EffectValues(this->GetTargetEffectNum()).density, back))
        return V(std::numeric_limits<double>::quiet_NaN());  // refused at set-up (WorkingCopy::SetUpTargetEffects)
    }
    if (noSpread) {
      // The grid was built to span the target (WorkingCopy::GridExtent), so it is integrated whole, as it
      // always was: clipping it to [E0 - T, E0] here would drop a point whose segment energy shift
      // moved E0 away from the grid's own centre (13N, segment 1, lowest point).
      // A fitted density (TargetEffect::HasDensityFit) lays the grid out for the largest density the fit may
      // reach, so the depth integral is cut at the back of the target, E0 - T(N): smooth in N, grid unchanged.
      double unused = 0.0;
      const double lo = (targetEffect->HasDensityFit() || targetEffect->IsLayerYield()) ? std::max(eLow, back) : eLow;
      GridQuadrature<V>(this, data, rd, lo, eHigh, true, one, one, outer, unused);
    } else {
      const double x1 = -1.0 / sqrt(3.0), x2 = 1.0 / sqrt(3.0);
      for (int i = 1; i < numPoints; i++) {
        double Ea = this->GetSubPoint(data, i)->GetCMEnergy(), Eb = this->GetSubPoint(data, i + 1)->GetCMEnergy();
        if (Eb >= surface || Ea <= back) continue;
        double a = std::min(Ea, surface), b = std::max(Eb, back);
        double mid = 0.5 * (a + b), half = 0.5 * (a - b);
        if (half < 1.0e-10) continue;
        for (int gp = 0; gp < 2; gp++) {
          double eDepth = mid + half * (gp == 0 ? x1 : x2);
          double beam = targetEffect->IsGaussian() ? targetEffect->SigmaCM(eDepth, configure, cm_factor_) : 0.0;
          if (!(beam > 0.)) beam = 0.0;
          double straggle = targetEffect->IsBohrStraggling()
                                ? targetEffect->BohrStragglingSigma(arealTo(eDepth) - arealTo(surface), cm_factor_, bohr_z1_, bohr_z2_)
                                : targetEffect->StragglingSigma(surface - eDepth);
          double sigmaEff = std::sqrt(beam * beam + straggle * straggle);
          V inner = 0.0;
          if (sigmaEff > 1.0e-10) {
            double den = 0.0;
            auto k = [&](double e) { return TargetEffect::GaussianKernel(e, eDepth, sigmaEff); };
            GridQuadrature<V>(this, data, rd, std::max(eLow, eDepth - range * sigmaEff), std::min(eHigh, eDepth + range * sigmaEff), true, k, k, inner, den);
            inner = (den > 0.0) ? inner / den : V(0.0);
          } else {
            inner = yieldAt(eDepth);
          }
          outer += half * inner;
        }
      }
    }
    // a thick target (density 0) and a layer yield are the yield per incident particle, 1e-24 int sigma/epsilon dE;
    // a thin target's quantity is the target-averaged cross section, that divided by 1e-24 N
    yield = (targetEffect->IsThickTarget() || targetEffect->IsLayerYield())
                ? outer
                : outer / (data.EffectValues(this->GetTargetEffectNum()).density * 1.0e-24);
    // The smooth edge of a window: blend with sigma(E_bar), not with sigma at the beam energy.  Done here, in
    // the linear kernel, so the gradient weights (TargetEffectWeights) carry it too.
    if (thinTarget && target_blend_weight_ < 1.0)
      yield = target_blend_weight_ * yield + (1.0 - target_blend_weight_) * this->TargetFallbackImpl<V>(data, rd);
  } else if (targetEffect->IsGaussian()) {
    const double s = targetEffect->SigmaCM(centroid, configure, cm_factor_);
    if (!(s > 0.)) {
      // no spread: the bare value at the point, interpolated on the grid
      for (int j = 1; j < numPoints; j++) {
        double Ea = this->GetSubPoint(data, j)->GetCMEnergy(), Eb = this->GetSubPoint(data, j + 1)->GetCMEnergy();
        if (centroid <= Ea && centroid >= Eb) {
          double frac = (Ea - centroid) / (Ea - Eb);
          return (1.0 - frac) * rd(j) + frac * rd(j + 1);
        }
      }
      return V(0.0);
    }
    auto k = [&](double e) { return TargetEffect::GaussianKernel(e, centroid, s); };
    V num = 0.0;
    double den = 0.0;
    GridQuadrature<V>(this, data, rd, eLow, eHigh, false, k, k, num, den);
    yield = (den > 0.0) ? num / den : V(0.0);
  } else if (targetEffect->IsUdr()) {
    // User-defined numerical resolution function (SAMMY UDR).  The kernel
    // R(tau) is the density of time-of-flight delays at this point's energy;
    // a sub-point at true lab energy E' is seen at the delay
    //   tau(E') = t(E0) - t(E'),
    // so the fitted value is  sum_i w_i sigma_i / sum_i w_i  with
    //   w_i = int R(tau) phi_i(tau) dtau,
    // phi_i the hat function of sub-point i on the delay axis: the exact
    // integral of the piecewise-linear kernel against the linearly
    // interpolated cross section, as in SAMMY's Udr_Resb.  Each product is a
    // quadratic on every piece between breakpoints of either grid, so
    // Simpson's rule evaluates it exactly.
    const DataSet::PointSetUp *setUp = data.Data().SetUpOf(*this);
    const UdrKernel *kernel = setUp ? setUp->udr.get() : nullptr;
    int numPoints = this->NumSubPoints();
    if (!kernel || kernel->tau.size() < 2 || numPoints < 2) {
      return V(0.0);
    }
    const double mass = setUp->udrMass;
    const double path = targetEffect->GetUdrFlightPath();
    const double cmPerLab = (lab_energy_ != 0.0) ? cm_energy_ / lab_energy_ : 1.0;
    const double t0 = UdrKinematics::TimeOfFlight(lab_energy_, mass, path);
    // Sub-points run from high to low energy, i.e. from large to small delay;
    // collect them in ascending delay.
    std::vector<double> tauSub(numPoints);
    std::vector<V> sigmaSub(numPoints);
    for (int i = 1; i <= numPoints; i++) {
      const Unit *sp = this->GetSubPoint(data, i);
      tauSub[numPoints - i] = t0 - UdrKinematics::TimeOfFlight(sp->GetCMEnergy() / cmPerLab, mass, path);
      sigmaSub[numPoints - i] = rd(i);
    }
    std::vector<double> w(numPoints, 0.0);
    const std::vector<double> &T = kernel->tau;
    const std::vector<double> &R = kernel->r;
    size_t j = 0;  // sub-point interval [tauSub[j], tauSub[j+1]]
    for (size_t k = 0; k + 1 < T.size(); k++) {
      double ta = T[k], tb = T[k + 1];
      if (tb <= tauSub.front() || ta >= tauSub.back()) continue;
      while (j + 1 < static_cast<size_t>(numPoints) && tauSub[j + 1] <= ta) j++;
      double lo = std::max(ta, tauSub.front());
      double hiAll = std::min(tb, tauSub.back());
      size_t jj = j;
      while (lo < hiAll && jj + 1 < static_cast<size_t>(numPoints)) {
        double hi = std::min(hiAll, tauSub[jj + 1]);
        if (hi > lo) {
          double sa = tauSub[jj], sb = tauSub[jj + 1], ds = sb - sa;
          auto rOf = [&](double x) { return R[k] + (R[k + 1] - R[k]) * (x - ta) / (tb - ta); };
          double mid = 0.5 * (lo + hi);
          double fl = rOf(lo), fm = rOf(mid), fh = rOf(hi);
          double h = (hi - lo) / 6.0;
          // hat of jj:  (sb - x)/ds ;  hat of jj+1:  (x - sa)/ds
          w[jj] += h * (fl * (sb - lo) + 4.0 * fm * (sb - mid) + fh * (sb - hi)) / ds;
          w[jj + 1] += h * (fl * (lo - sa) + 4.0 * fm * (mid - sa) + fh * (hi - sa)) / ds;
        }
        lo = hi;
        if (lo >= tauSub[jj + 1]) jj++;
      }
    }
    V numerator = 0.0;
    double denominator = 0.0;
    for (int i = 0; i < numPoints; i++) {
      numerator += w[i] * sigmaSub[i];
      denominator += w[i];
    }
    yield = (denominator > 0.0) ? numerator / denominator : V(0.0);
  } else if (targetEffect->IsBeamProfile()) {
    // Beam-profile kernel (photodissociation in a broad, skewed gamma beam
    // with an event-by-event reconstructed energy, cf. Haverson 2026 App. A):
    //
    //   Y = sum_i G(E_i) W(E_i) D(E_i) sigma(E_i) / sum_i G(E_i) W(E_i),
    //
    // G  the absolute beam energy profile (a sum of skewed Gaussians, not
    //    centred on the data point),
    // W  the fraction of the Gaussian detector resolution (sigma s) that
    //    lands inside the point's energy window [a, b], shifted with the
    //    segment's energy shift; 1 when the point has no window,
    // D  the detailed-balance factor sigma(gamma,a)/sigma(a,gamma) relative
    //    to its value at the point's own (unshifted) energy; 1 unless
    //    requested.  The measurement spreads its luminosity over G*W and
    //    averages sigma(gamma,a) = f_db * sigma(a,gamma); the datum is that
    //    average divided by f_db at the point's energy.  So D weights the
    //    numerator only: in the denominator too it would replace f_db(E0) by
    //    its beam average and cancel its own normalisation (a 1-5 % error
    //    when the kernel's centroid is tens of keV off the point energy).
    const double s = cm_factor_ * targetEffect->GetBeamTpcSigma();  // the line's lab sigma, this pair's frame
    const bool window = this->HasBinWindow();
    // The energy shift moved the point (and the sub-point grid); the window
    // is stored unshifted, so shift it by the same c.m. amount.
    double cmPerLab = (lab_energy_ != 0.0) ? cm_energy_ / lab_energy_ : 1.0;
    double originalCM = this->GetOriginalEnergy() * cmPerLab;
    double delta = cm_energy_ - originalCM;
    double a = bin_low_cm_ + delta;
    double b = bin_high_cm_ + delta;
    const bool photo = targetEffect->IsBeamPhotodissociation() && photo_mcn_ > 0.0;
    // g(E) = E / E_gamma(E)^2 is the energy dependence of the detailed-balance
    // factor mu c^2 E / E_gamma^2, with the photon energy including the
    // recoil of the compound nucleus: E_gamma^2/(2 M) - E_gamma + (E + Q) = 0.
    auto dbWeight = [&](double e) {
      double sum = e + photo_qgamma_;
      double arg = 1.0 - 2.0 * sum / photo_mcn_;
      double egamma = (arg > 0.0) ? photo_mcn_ * (1.0 - std::sqrt(arg)) : sum;
      return (egamma > 0.0) ? e / (egamma * egamma) : 0.0;
    };
    double dbNorm = photo ? dbWeight(originalCM) : 1.0;
    auto kernel = [&](double e) {
      double k = targetEffect->BeamProfileWeight(e / cm_factor_);  // the profile is lab; e is c.m.
      if (window) {
        if (s > 0.0)
          k *= 0.5 * (std::erf((b - e) / (s * std::sqrt(2.0))) - std::erf((a - e) / (s * std::sqrt(2.0))));
        else
          k *= (e >= a && e <= b) ? 1.0 : 0.0;
      }
      return k;
    };
    auto weighted = [&](double e) {
      return kernel(e) * ((photo && dbNorm > 0.0) ? dbWeight(e) / dbNorm : 1.0);
    };
    V numerator = 0.0;
    double denominator = 0.0;
    GridQuadrature<V>(this, data, rd, eLow, eHigh, false, weighted, kernel, numerator, denominator);
    yield = (denominator > 0.0) ? numerator / denominator : V(0.0);
  }
  return yield;
}
void Unit::IntegrateTargetEffect(const WorkingCopy &data, const Config &configure, azure2::Results &results) {
  this->SetFitCrossSection(results, this->IntegrateTargetEffectImpl<double>(data, configure, [this, &data, &results](int j) {
    return this->GetSubPoint(data, j)->GetFitCrossSection(results);
  }));
}
double Unit::IntegrateSubPointValues(const WorkingCopy &data, const Config &configure, const std::vector<double> &values) {
  return this->IntegrateTargetEffectImpl<double>(data, configure, [&values](int j) { return values[j - 1]; });
}

bool Unit::TargetEffectWeights(const WorkingCopy &data, const Config &configure, std::vector<double> &c) {
  const int n = this->NumSubPoints();
  c.assign(n > 0 ? n : 0, 0.0);
  TELinearForm y = this->IntegrateTargetEffectImpl<TELinearForm>(data, configure, [](int j) { return TELinearForm::Unit(j); });
  if (!y.Valid()) return false;
  for (const auto &t : y.Terms()) {
    if (t.first < 1 || t.first > n) return false;
    c[t.first - 1] += t.second;
  }
  return true;
}
