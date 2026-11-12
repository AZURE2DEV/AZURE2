// The beam-profile kernel (Kernels.h): photodissociation in a broad, skewed gamma beam with an event-by-event
// reconstructed energy (cf. Haverson 2026 App. A).
#include <algorithm>
#include <cmath>

#include "azure2/experiment/KernelQuadrature.h"
#include "azure2/experiment/Kernels.h"
#include "azure2/experiment/LinearForm.h"

namespace azure2 {
namespace experiment {
namespace {
constexpr double kPi = 3.14159265358979323846;
}

/*!
 * The weighted sum of the skewed-Gaussian components, G(E|xi,omega,alpha) = exp(-((E-xi)/omega)^2/2) /(omega
 * sqrt(2 pi)) * (1 + erf(alpha (E-xi)/(omega sqrt 2))), each optionally zeroed outside its mean +- nCut standard
 * deviations.
 */
double BeamProfileWeight(const KernelParameters &p, double energy) {
  const double pi = kPi;
  double weight = 0.0;
  for (const BeamProfileComponent &c : *p.profile) {
    if (!(c.omega > 0.0)) continue;
    double z = (energy - c.xi) / c.omega;
    if (p.truncation > 0.0) {
      double delta = c.alpha / std::sqrt(1.0 + c.alpha * c.alpha);
      double mean = c.xi + c.omega * delta * std::sqrt(2.0 / pi);
      double sd = c.omega * std::sqrt(1.0 - 2.0 * delta * delta / pi);
      if (std::fabs(energy - mean) > p.truncation * sd) continue;
    }
    weight += c.weight * std::exp(-0.5 * z * z) / (c.omega * std::sqrt(2.0 * pi)) * (1.0 + std::erf(c.alpha * z / std::sqrt(2.0)));
  }
  return weight;
}

/// The union of xi +- 5 omega over the components (tighter if the profile is truncated), lab.
void BeamProfileSupport(const KernelParameters &p, double &low, double &high) {
  const double pi = kPi;
  low = 1.0e300;
  high = -1.0e300;
  for (const BeamProfileComponent &c : *p.profile) {
    double lo = c.xi - 5.0 * c.omega;
    double hi = c.xi + 5.0 * c.omega;
    if (p.truncation > 0.0) {
      double delta = c.alpha / std::sqrt(1.0 + c.alpha * c.alpha);
      double mean = c.xi + c.omega * delta * std::sqrt(2.0 / pi);
      double sd = c.omega * std::sqrt(1.0 - 2.0 * delta * delta / pi);
      lo = std::max(lo, mean - p.truncation * sd);
      hi = std::min(hi, mean + p.truncation * sd);
    }
    low = std::min(low, lo);
    high = std::max(high, hi);
  }
}

namespace {

/*!
 *   Y = sum_i G(E_i) W(E_i) D(E_i) sigma(E_i) / sum_i G(E_i) W(E_i),
 *
 * G  the absolute beam energy profile (a sum of skewed Gaussians, not centred on the data point),
 * W  the fraction of the Gaussian detector resolution (sigma s) that lands inside the point's energy window [a, b],
 *    shifted with the segment's energy shift; 1 when the point has no window,
 * D  the detailed-balance factor sigma(gamma,a)/sigma(a,gamma) relative to its value at the point's own (unshifted)
 *    energy; 1 unless requested.  The measurement spreads its luminosity over G*W and averages sigma(gamma,a) =
 *    f_db * sigma(a,gamma); the datum is that average divided by f_db at the point's energy.  So D weights the
 *    numerator only: in the denominator too it would replace f_db(E0) by its beam average and cancel its own
 *    normalisation (a 1-5 % error when the kernel's centroid is tens of keV off the point energy).
 */
template <class V>
V FoldBeamProfile(const KernelParameters &p, const KernelPoint &point, const KernelNodes &nodes, const V *values) {
  const int numPoints = nodes.n;
  if (numPoints < 2) return V(0.0);
  const double eHigh = nodes.High();
  const double eLow = nodes.Low();
  const double cmFactor = point.cmFactor;
  const double s = cmFactor * p.tpcSigma;  // the record's lab sigma, this pair's frame
  const bool window = point.binLow == point.binLow && point.binHigh == point.binHigh;
  // The energy shift moved the point (and the nodes); the window is stored unshifted, so shift it by the same c.m.
  // amount.
  double cmPerLab = (point.labEnergy != 0.0) ? point.energy / point.labEnergy : 1.0;
  double originalCM = point.nominalLabEnergy * cmPerLab;
  double delta = point.energy - originalCM;
  double a = point.binLow + delta;
  double b = point.binHigh + delta;
  const bool photo = p.detailedBalance && point.photoMass > 0.0;
  // g(E) = E / E_gamma(E)^2 is the energy dependence of the detailed-balance factor mu c^2 E / E_gamma^2, with the
  // photon energy including the recoil of the compound nucleus: E_gamma^2/(2 M) - E_gamma + (E + Q) = 0.
  auto dbWeight = [&](double e) {
    double sum = e + point.photoQ;
    double arg = 1.0 - 2.0 * sum / point.photoMass;
    double egamma = (arg > 0.0) ? point.photoMass * (1.0 - std::sqrt(arg)) : sum;
    return (egamma > 0.0) ? e / (egamma * egamma) : 0.0;
  };
  double dbNorm = photo ? dbWeight(originalCM) : 1.0;
  auto kernel = [&](double e) {
    double k = BeamProfileWeight(p, e / cmFactor);  // the profile is lab; e is c.m.
    if (window) {
      if (s > 0.0)
        k *= 0.5 * (std::erf((b - e) / (s * std::sqrt(2.0))) - std::erf((a - e) / (s * std::sqrt(2.0))));
      else
        k *= (e >= a && e <= b) ? 1.0 : 0.0;
    }
    return k;
  };
  auto weighted = [&](double e) { return kernel(e) * ((photo && dbNorm > 0.0) ? dbWeight(e) / dbNorm : 1.0); };
  V numerator = 0.0;
  double denominator = 0.0;
  GridQuadrature<V>(nodes, values, eLow, eHigh, false, weighted, kernel, numerator, denominator);
  return (denominator > 0.0) ? numerator / denominator : V(0.0);
}

/*!
 * Where the beam is (the profile is absolute, not centred on the point), narrowed to the point's energy window plus
 * the resolution tails when it has one: 4 s, or for a THM (HOES) point the 5 s of its Gaussian fold (a narrow peak
 * 4-5 s from a point dominates the tail there).  The floor is the lowest node energy (a THM point may sit below the
 * entrance threshold).
 */
void BeamProfileGridSupport(const KernelParameters &p, const SupportPoint &sp, double &high, double &low) {
  BeamProfileSupport(p, low, high);  // lab
  low *= sp.cmFactor;
  high *= sp.cmFactor;
  double s = sp.cmFactor * p.tpcSigma;
  double tails = sp.thm ? kThmConvolutionRange : 4.0;
  if (sp.binLow == sp.binLow && sp.binHigh == sp.binHigh) {
    low = std::max(low, sp.binLow - tails * s);
    high = std::min(high, sp.binHigh + tails * s);
  }
  if (low < sp.minimumEnergy) low = sp.minimumEnergy;
  if (high <= low) high = low + 0.001;
}

}  // namespace

Kernel BeamProfileKernel() {
  Kernel k;
  k.name = "beam_profile";
  k.keyword = "beam.profile";
  k.fold = FoldBeamProfile<double>;
  k.weights = FoldBeamProfile<LinearForm>;
  k.support = BeamProfileGridSupport;
  // no thinLimit: the support is at least 1 keV wide (BeamProfileGridSupport), never a single node
  k.doc = "an absolute beam energy profile (skewed Gaussians), the detector's Gaussian resolution inside the point's "
          "energy window, and optionally the detailed-balance weight of a photodissociation measurement";
  return k;
}

}  // namespace experiment
}  // namespace azure2
