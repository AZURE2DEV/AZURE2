// The free-gas Doppler broadening (SAMMY's FGM, Kernels.h): a broadening of the node values applied before the
// effect's resolution kernel.
#include <algorithm>
#include <cmath>
#include <map>
#include <utility>
#include <vector>

#include "azure2/experiment/Kernels.h"
#include "azure2/experiment/LinearForm.h"

namespace azure2 {
namespace experiment {
namespace {
constexpr double kPi = 3.14159265358979323846;
}

double DopplerU(double cmFactor, double temperature) {
  if (!(cmFactor > 0.0) || !(cmFactor < 1.0)) return 0.0;
  const double kT = 8.617333262e-11 * temperature;     // MeV
  return std::sqrt(kT * (1.0 - cmFactor) / cmFactor);  // m/M = (1 - c)/c for c = M/(m+M)
}

void DopplerWindow(const KernelParameters &p, double eCM, double cmFactor, double &lo, double &hi) {
  lo = hi = eCM;
  const double u = DopplerU(cmFactor, p.dopplerLayout);
  if (!(u > 0.0) || !(eCM > 0.0)) return;
  // the kernel is a Gaussian in velocity, exp(-(v-w)^2/u^2): standard deviation u/sqrt(2)
  const double v = std::sqrt(eCM / cmFactor), dv = kConvolutionRange * u / std::sqrt(2.0);
  hi = cmFactor * (v + dv) * (v + dv);
  lo = (v > dv) ? cmFactor * (v - dv) * (v - dv) : 0.0;
}

namespace {

/*!
 *   sigma_D(E) = 1/(u sqrt(pi) v^2) int_0^inf [exp(-(v-w)^2/u^2) - exp(-(v+w)^2/u^2)] w^2 sigma(w^2) dw,
 * v = sqrt(E), w = sqrt(E'), E the lab-equivalent energy (E_cm over the point's c.m./lab factor), u^2 = kT m/M, sigma
 * linear in E between nodes: the weights w with sigma_D(eCM) = sum w.second * value(node w.first), 4-point
 * Gauss-Legendre on pieces no wider than u0/4.  The window (layout temperature) and the pieces (nominal temperature)
 * do not move with T: only the kernel values do, so sigma_D is smooth in a fitted T; dwdT, if given, receives dw/dT
 * (per K) on the same nodes.
 */
void DopplerWeights(const KernelParameters &p, const KernelPoint &point, const KernelNodes &nodes, double eCM,
                    std::vector<std::pair<int, double>> &w, std::vector<std::pair<int, double>> *dwdT) {
  w.clear();
  if (dwdT) dwdT->clear();
  const int n = nodes.n;
  const double c = point.cmFactor;
  const double T = point.temperature;
  const double u = DopplerU(c, T);
  if (n < 2 || !(u > 0.0) || !(eCM > 0.0)) return;
  const double u0 = DopplerU(c, p.dopplerNominal);
  double lo, hi;
  DopplerWindow(p, eCM, c, lo, hi);
  const double v = std::sqrt(eCM / c);
  const double norm = 1.0 / (u * std::sqrt(kPi) * v * v);
  // 4-point Gauss-Legendre on [-1, 1]
  static const double gx[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526};
  static const double gw[4] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538};
  auto energy = [&nodes](int j) { return nodes.Energy(j); };
  // nodes run from high to low energy: the first interval (j, j+1) whose lower end lies below hi
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
    for (int q = 0; q < pieces; q++) {
      const double mid = wb + (q + 0.5) * hw, half = 0.5 * hw;
      for (int g = 0; g < 4; g++) {
        const double x = mid + half * gx[g];
        const double k = (std::exp(-(v - x) * (v - x) / (u * u)) - std::exp(-(v + x) * (v + x) / (u * u))) * x * x;
        const double frac = (Ea - c * x * x) / (Ea - Eb);  // 0 at node j, 1 at node j + 1
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

/// The broadened values at the energies e[0..m): out[r] = sum of the weights at e[r] times the node values.
template <class V>
void BroadenDoppler(const KernelParameters &p, const KernelPoint &point, const KernelNodes &nodes, const V *values,
                    const double *e, int m, V *out) {
  std::vector<std::pair<int, double>> w;
  for (int r = 0; r < m; r++) {
    DopplerWeights(p, point, nodes, e[r], w, nullptr);
    V s = 0.0;
    for (const auto &t : w) s += t.second * values[t.first - 1];
    out[r] = s;
  }
}

/// d/dT of BroadenDoppler<double>.
void BroadenDopplerDT(const KernelParameters &p, const KernelPoint &point, const KernelNodes &nodes,
                      const double *values, const double *e, int m, double *out) {
  std::vector<std::pair<int, double>> w, dw;
  for (int r = 0; r < m; r++) {
    DopplerWeights(p, point, nodes, e[r], w, &dw);
    double s = 0.0;
    for (const auto &t : dw) s += t.second * values[t.first - 1];
    out[r] = s;
  }
}

/// The resolution kernel's interval widened by the Doppler width at each end: every node the kernel reaches needs
/// its own Doppler window inside the grid.
void DopplerSupport(const KernelParameters &p, const SupportPoint &s, double &high, double &low) {
  double lo, hi, unused;
  DopplerWindow(p, high, s.cmFactor, unused, hi);
  DopplerWindow(p, low, s.cmFactor, lo, unused);
  high = std::max(high, hi);
  low = std::max(std::min(low, lo), 0.001);
}

}  // namespace

Kernel DopplerKernel() {
  Kernel k;
  k.name = "doppler";
  k.stage = Kernel::Stage::kBroadening;
  k.keyword = "target.doppler free_gas";
  k.parameters = {{"temperature", "temperature", "K", true, "the effective temperature of the target atoms"}};
  k.broaden = BroadenDoppler<double>;
  k.broadenForm = BroadenDoppler<LinearForm>;
  k.broadenDT = BroadenDopplerDT;
  k.broadeningWeights = DopplerWeights;
  k.thinLimit = AverageThinLimit;
  k.support = DopplerSupport;
  k.doc = "SAMMY's free-gas Doppler broadening of the cross section, applied before the resolution kernel";
  return k;
}

}  // namespace experiment
}  // namespace azure2
