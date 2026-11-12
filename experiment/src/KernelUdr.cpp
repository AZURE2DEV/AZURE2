// The user-defined numerical resolution kernel (SAMMY's UDR, Kernels.h): a time-of-flight delay distribution per
// point, built at set-up from the detector's table (KernelPoint::udr).
#include <algorithm>
#include <vector>

#include "azure2/experiment/Kernels.h"
#include "azure2/experiment/LinearForm.h"
#include "azure2/experiment/UdrKernel.h"

namespace azure2 {
namespace experiment {
namespace {

/*!
 * The kernel R(tau) is the density of time-of-flight delays at this point's energy; a node at true lab energy E' is
 * seen at the delay tau(E') = t(E0) - t(E'), so the folded value is sum_i w_i sigma_i / sum_i w_i with
 * w_i = int R(tau) phi_i(tau) dtau, phi_i the hat function of node i on the delay axis: the exact integral of the
 * piecewise-linear kernel against the linearly interpolated cross section, as in SAMMY's Udr_Resb.  Each product is a
 * quadratic on every piece between breakpoints of either grid, so Simpson's rule evaluates it exactly.
 */
template <class V>
V FoldUdr(const KernelParameters &p, const KernelPoint &point, const KernelNodes &nodes, const V *values) {
  const int numPoints = nodes.n;
  if (numPoints < 2) return V(0.0);
  const UdrKernel *kernel = point.udr;
  if (!kernel || kernel->tau.size() < 2) return V(0.0);
  const double mass = point.udrMass;
  const double path = p.flightPath;
  const double cmPerLab = (point.labEnergy != 0.0) ? point.energy / point.labEnergy : 1.0;
  const double t0 = UdrKinematics::TimeOfFlight(point.labEnergy, mass, path);
  // Nodes run from high to low energy, i.e. from large to small delay; collect them in ascending delay.
  std::vector<double> tauSub(numPoints);
  std::vector<V> sigmaSub(numPoints);
  for (int i = 1; i <= numPoints; i++) {
    tauSub[numPoints - i] = t0 - UdrKinematics::TimeOfFlight(nodes.Energy(i) / cmPerLab, mass, path);
    sigmaSub[numPoints - i] = values[i - 1];
  }
  std::vector<double> w(numPoints, 0.0);
  const std::vector<double> &T = kernel->tau;
  const std::vector<double> &R = kernel->r;
  size_t j = 0;  // node interval [tauSub[j], tauSub[j+1]]
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
  return (denominator > 0.0) ? numerator / denominator : V(0.0);
}

/*!
 * A delay tau maps the point's nominal flight time t0 to a true flight time t0 - tau: the kernel's support in delay
 * becomes a window in true energy, extending above the point for positive delays (moderation) and below it for
 * negative ones.
 */
void UdrSupport(const KernelParameters &p, const SupportPoint &s, double &high, double &low) {
  const UdrKernel &kernel = *s.udr;
  const double mass = s.udrMass;
  double t0 = UdrKinematics::TimeOfFlight(s.labEnergy, mass, p.flightPath);
  double eHighLab = UdrKinematics::EnergyFromTime(t0 - kernel.tauMax, mass, p.flightPath);
  double eLowLab = UdrKinematics::EnergyFromTime(t0 - kernel.tauMin, mass, p.flightPath);
  double cmPerLab = s.energy / s.labEnergy;
  if (eHighLab <= s.labEnergy) eHighLab = s.labEnergy * 1.0001;
  if (eLowLab <= 0.0 || eLowLab >= s.labEnergy) eLowLab = s.labEnergy * 0.9999;
  high = eHighLab * cmPerLab;
  low = std::max(eLowLab * cmPerLab, 0.001);
}

}  // namespace

Kernel UdrKernelRow() {
  Kernel k;
  k.name = "udr";
  k.keyword = "detector.udr";
  k.fold = FoldUdr<double>;
  k.weights = FoldUdr<LinearForm>;
  k.thinLimit = AverageThinLimit;
  k.support = UdrSupport;
  k.doc = "SAMMY's user-defined time-of-flight resolution (UDR): a tabulated delay distribution convolved with the "
          "burst and the channel width";
  return k;
}

}  // namespace experiment
}  // namespace azure2
