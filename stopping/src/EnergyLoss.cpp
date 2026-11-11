#include "azure2/stopping/EnergyLoss.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "azure2/stopping/Stopping.h"

namespace azure2 {
namespace stopping {

namespace {

// Gauss-Kronrod 7-15 nodes (positive half, the last is 0) and weights; Gauss weights at the odd Kronrod nodes.
constexpr std::array<double, 8> kXgk = {0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
                                        0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
                                        0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
                                        0.207784955007898467600689403773245, 0.0};
constexpr std::array<double, 8> kWgk = {0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
                                        0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
                                        0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
                                        0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
constexpr std::array<double, 4> kWg = {0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
                                       0.381830050505118944950369775488975, 0.417959183673469387755102040816327};
constexpr int kMaxDepth = 60;

/// One Gauss-Kronrod 7-15 panel on [a, b]: the Kronrod value and |Kronrod - Gauss|.
std::pair<double, double> Panel(const std::function<double(double)> &g, double a, double b) {
  const double c = 0.5 * (a + b);
  const double h = 0.5 * (b - a);
  const double fc = g(c);
  double kronrod = kWgk[7] * fc;
  double gauss = kWg[3] * fc;
  for (int j = 0; j < 7; ++j) {
    const double f1 = g(c - h * kXgk[j]);
    const double f2 = g(c + h * kXgk[j]);
    kronrod += kWgk[j] * (f1 + f2);
    if (j % 2 == 1) gauss += kWg[j / 2] * (f1 + f2);
  }
  return {kronrod * h, std::fabs((kronrod - gauss) * h)};
}

double Adaptive(const std::function<double(double)> &g, double a, double b, double abs_tol, int depth) {
  const auto [value, error] = Panel(g, a, b);
  if (error <= abs_tol || std::fabs(value) <= 1e-300) return value;
  if (depth >= kMaxDepth) {
    std::ostringstream o;
    o << "stopping: quadrature did not converge on [" << a << ", " << b << "] (sqrt(MeV)) after " << kMaxDepth << " bisections";
    throw std::runtime_error(o.str());
  }
  const double m = 0.5 * (a + b);
  return Adaptive(g, a, m, 0.5 * abs_tol, depth + 1) + Adaptive(g, m, b, 0.5 * abs_tol, depth + 1);
}

/// 1/epsilon in cm^2 per (MeV atom) x atoms: dN/dE = 1 / (epsilon [eV/(1e15 atoms/cm^2)] x 1e-21 MeV cm^2).
double InverseStopping(const StoppingFunction &epsilon, double energy_MeV) {
  const double e = epsilon(energy_MeV);
  if (!(e > 0.0)) {
    std::ostringstream o;
    o << "stopping: epsilon(" << energy_MeV << " MeV) = " << e << " is not positive";
    throw std::domain_error(o.str());
  }
  return 1.0 / (e * kStoppingToMeVcm2);
}

}  // namespace

double IntegrateInSqrtE(const std::function<double(double)> &f, double a_MeV, double b_MeV, double rel_tol,
                        const std::vector<double> &breaks_MeV) {
  if (!(a_MeV >= 0.0) || !(b_MeV >= a_MeV)) throw std::invalid_argument("stopping: integral needs 0 <= a <= b");
  if (a_MeV == b_MeV) return 0.0;
  const auto g = [&f](double u) { return u > 0.0 ? 2.0 * u * f(u * u) : 0.0; };
  // the pieces in u, cut at the break points strictly inside (a, b)
  std::vector<double> nodes{std::sqrt(a_MeV)};
  for (double e : breaks_MeV) {
    if (e > a_MeV && e < b_MeV) nodes.push_back(std::sqrt(e));
  }
  nodes.push_back(std::sqrt(b_MeV));
  std::sort(nodes.begin(), nodes.end());
  double first = 0.0;
  for (std::size_t k = 0; k + 1 < nodes.size(); ++k) first += std::fabs(Panel(g, nodes[k], nodes[k + 1]).first);
  const double piece_tol = std::max(rel_tol * first / static_cast<double>(nodes.size() - 1), 1e-300);
  double sum = 0.0;
  for (std::size_t k = 0; k + 1 < nodes.size(); ++k) sum += Adaptive(g, nodes[k], nodes[k + 1], piece_tol, 0);
  return sum;
}

double ArealDensityBetween(const StoppingFunction &epsilon, double energy_in_MeV, double energy_out_MeV,
                           const std::vector<double> &breaks_MeV) {
  if (!(energy_out_MeV >= 0.0) || !(energy_in_MeV >= energy_out_MeV)) {
    throw std::invalid_argument("stopping: areal density between E_in and E_out needs E_in >= E_out >= 0");
  }
  return IntegrateInSqrtE([&epsilon](double e) { return InverseStopping(epsilon, e); }, energy_out_MeV, energy_in_MeV,
                          kQuadratureRelTol, breaks_MeV);
}

double Range(const StoppingFunction &epsilon, double energy_MeV, const std::vector<double> &breaks_MeV) {
  return ArealDensityBetween(epsilon, energy_MeV, 0.0, breaks_MeV);
}

double EnergyAfter(const StoppingFunction &epsilon, double energy_in_MeV, double areal_density,
                   const std::vector<double> &breaks_MeV) {
  if (!(areal_density >= 0.0)) throw std::invalid_argument("stopping: areal density must be >= 0");
  if (!(energy_in_MeV > 0.0)) throw std::invalid_argument("stopping: energy must be > 0");
  if (areal_density == 0.0) return energy_in_MeV;
  const double range = Range(epsilon, energy_in_MeV, breaks_MeV);
  if (areal_density >= range) {
    std::ostringstream o;
    o << "stopping: the projectile stops: areal density " << areal_density << " atoms/cm^2 >= range " << range;
    throw std::domain_error(o.str());
  }
  // g(E_out) = N(E_in, E_out) - N: g(0) > 0, g(E_in) < 0, g' = -1/epsilon(E_out)
  double lo = 0.0;
  double hi = energy_in_MeV;
  double e = std::clamp(energy_in_MeV - areal_density * epsilon(energy_in_MeV) * kStoppingToMeVcm2, 0.0, energy_in_MeV);
  for (int iteration = 0; iteration < 200; ++iteration) {
    const double g = ArealDensityBetween(epsilon, energy_in_MeV, e, breaks_MeV) - areal_density;
    if (g > 0.0) lo = e; else hi = e;
    double next = e + g * epsilon(e) * kStoppingToMeVcm2;
    if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
    if (std::fabs(next - e) <= 1e-13 * energy_in_MeV || hi - lo <= 1e-13 * energy_in_MeV) return next;
    e = next;
  }
  throw std::runtime_error("stopping: energy after the layer did not converge");
}

StoppingFunction MaterialStopping(const Material &material, const Projectile &projectile) {
  return [material, projectile](double e) { return Stopping(material, projectile, e).total; };
}

double ArealDensityBetween(const Material &material, const Projectile &projectile, double energy_in_MeV,
                           double energy_out_MeV) {
  return ArealDensityBetween(MaterialStopping(material, projectile), energy_in_MeV, energy_out_MeV,
                             StoppingBreakpoints(material, projectile));
}

double EnergyAfter(const Material &material, const Projectile &projectile, double energy_in_MeV, double areal_density) {
  return EnergyAfter(MaterialStopping(material, projectile), energy_in_MeV, areal_density,
                     StoppingBreakpoints(material, projectile));
}

double Range(const Material &material, const Projectile &projectile, double energy_MeV) {
  return Range(MaterialStopping(material, projectile), energy_MeV, StoppingBreakpoints(material, projectile));
}

double ArealDensityScaleDerivative(const StoppingFunction &epsilon, double energy_in_MeV, double energy_out_MeV,
                                   const std::vector<double> &breaks_MeV) {
  return -ArealDensityBetween(epsilon, energy_in_MeV, energy_out_MeV, breaks_MeV);
}

double EnergyAfterScaleDerivative(const StoppingFunction &epsilon, double energy_in_MeV, double areal_density,
                                  const std::vector<double> &breaks_MeV) {
  const double e_out = EnergyAfter(epsilon, energy_in_MeV, areal_density, breaks_MeV);
  return -areal_density * epsilon(e_out) * kStoppingToMeVcm2;
}

double RangeScaleDerivative(const StoppingFunction &epsilon, double energy_MeV, const std::vector<double> &breaks_MeV) {
  return -Range(epsilon, energy_MeV, breaks_MeV);
}

}  // namespace stopping
}  // namespace azure2
