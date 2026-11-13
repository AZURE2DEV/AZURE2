#include "azure2/physics/transfer/ZeroRangeWeight.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "azure2/physics/transfer/RadialNumerics.h"

namespace azure2 {
namespace physics {
namespace transfer {

namespace {

using cplx = std::complex<double>;
const double kPi = std::acos(-1.0);

void CheckOverlap(const RadialOverlap &o, const ZeroRangeWeightRequest &q) {
  if (o.l != 0) throw std::invalid_argument("zero-range weight: the overlap must have l = 0");
  if (!(o.step_fm > 0.0) || o.u.size() < 8) throw std::invalid_argument("zero-range weight: an empty overlap");
  if (static_cast<int>(o.u.size()) - FirstNode(o) < 4)
    throw std::invalid_argument("zero-range weight: fewer than 4 nodes beyond r_min");
  if (!(q.c_i > 0.0) || !(q.c_f > 0.0)) throw std::invalid_argument("zero-range weight: scale factors must be positive");
}

/// u(r_p)/r_p on the overlap's grid (0 at the origin).
std::vector<double> PhiOf(const RadialOverlap &o) {
  std::vector<double> phi(o.u.size(), 0.0);
  for (std::size_t p = 1; p < o.u.size(); ++p) phi[p] = o.u[p] / (p * o.step_fm);
  return phi;
}

}  // namespace

ZeroRangeEntrance PrepareZeroRangeEntrance(const RadialOverlap &overlap, const WaveChannel &entrance,
                                           const ZeroRangeWeightRequest &q) {
  CheckOverlap(overlap, q);
  const int n = static_cast<int>(overlap.u.size()), first = FirstNode(overlap);
  const double h = overlap.step_fm;
  const std::vector<double> phi = PhiOf(overlap);
  const RadialChannel channel = ChannelFor(entrance, q.c_i * h, n, q.l_cap, q.matching);
  ZeroRangeEntrance e;
  e.k_fm = channel.k_fm;
  std::vector<double> amplitude, f(n - first);
  double largest = 0.0;
  for (int l = 0;; ++l) {
    if (l > q.l_cap) throw std::runtime_error("the partial-wave sum does not converge by l = " + std::to_string(q.l_cap));
    if ((l + 1.0) * n * sizeof(cplx) / 1048576.0 > q.max_megabytes)
      throw std::runtime_error("the entrance partial waves would need more than " +
                               std::to_string(static_cast<int>(q.max_megabytes)) + " MB (" + std::to_string(n) +
                               " radial points, l = " + std::to_string(l) + ")");
    PartialWave pw;
    try {
      pw = channel.plane_wave ? SolvePartialWave(channel, l, 2 * l, WaveMatching(), n)
                              : SolvePartialWave(channel, l, 2 * l, MatchingFor(channel, l, n, q.matching), n);
    } catch (const std::runtime_error &) {
      throw std::runtime_error("the entrance wave l = " + std::to_string(l) + " could not be normalized");
    }
    for (int p = first; p < n; ++p) f[p - first] = std::fabs(phi[p]) * std::abs(pw.u[p]);
    const double a = (2.0 * l + 1.0) * Simpson(f, h);
    e.u.push_back(std::move(pw.u));
    amplitude.push_back(a);
    largest = std::max(largest, a);
    if (l >= 2 && a < q.entrance_cut * largest && a < amplitude[l - 1]) break;
  }
  e.remaining.assign(amplitude.size() + 1, 0.0);
  for (std::size_t l = amplitude.size(); l-- > 0;) e.remaining[l] = e.remaining[l + 1] + amplitude[l];
  e.sigma = CoulombPhases(channel, static_cast<int>(e.u.size()) - 1);
  return e;
}

ZeroRangeWeight ZeroRangeWeightSum(const RadialOverlap &overlap, const ZeroRangeEntrance &entrance,
                                   const WaveChannel &exit, const std::vector<double> &cos_theta,
                                   const ZeroRangeWeightRequest &q) {
  CheckOverlap(overlap, q);
  const int n = static_cast<int>(overlap.u.size()), first = FirstNode(overlap);
  const double h = overlap.step_fm;
  if (entrance.u.empty() || static_cast<int>(entrance.u.front().size()) != n)
    throw std::invalid_argument("zero-range weight: the entrance was prepared on another grid");
  const int l_top = static_cast<int>(entrance.u.size()) - 1;
  const std::vector<double> phi = PhiOf(overlap);
  const RadialChannel channel = ChannelFor(exit, q.c_f * h, n, l_top, q.matching);
  const std::vector<double> sigma = CoulombPhases(channel, l_top);
  // Plane waves come from one pass of the Riccati-Bessel recurrence for every l.
  const DistortedWaves plane =
      channel.plane_wave ? SolveStoredWaves(exit, q.c_f * h, n, l_top, q.matching) : DistortedWaves();
  const int nk = static_cast<int>(cos_theta.size());
  ZeroRangeWeight out;
  out.sum.assign(nk, 0.0);
  std::vector<double> p_prev(nk, 1.0), p_l(nk, 1.0);  // Legendre recurrence per direction
  std::vector<cplx> f(n - first);
  int stop = -1;
  for (int l = 0; l <= l_top; ++l) {
    for (int k = 0; k < nk; ++k) {
      if (l == 1) {
        p_prev[k] = 1.0;
        p_l[k] = cos_theta[k];
      } else if (l > 1) {
        const double next = ((2.0 * l - 1.0) * cos_theta[k] * p_l[k] - (l - 1.0) * p_prev[k]) / l;
        p_prev[k] = p_l[k];
        p_l[k] = next;
      }
    }
    PartialWave pw;
    try {
      pw = channel.plane_wave ? plane.waves[l]
                              : SolvePartialWave(channel, l, 2 * l, MatchingFor(channel, l, n, q.matching), n);
    } catch (const std::runtime_error &) {
      throw std::runtime_error("the exit wave l = " + std::to_string(l) + " could not be normalized");
    }
    double u_max = 0.0;
    for (int p = first; p < n; ++p) {
      f[p - first] = phi[p] * pw.u[p] * entrance.u[l][p];
      u_max = std::max(u_max, std::abs(pw.u[p]));
    }
    const cplx radial = Simpson(f, h);
    const cplx phase = std::polar(1.0, sigma[l]) * std::polar(1.0, entrance.sigma[l]);
    double smallest = HUGE_VAL;
    for (int k = 0; k < nk; ++k) {
      out.sum[k] += (2.0 * l + 1.0) * phase * p_l[k] * radial;
      smallest = std::min(smallest, std::abs(out.sum[k]));
    }
    if (l == 0) out.end_integrand = std::abs(f.back());
    out.l_max = l;
    if (stop < 0 && l >= 1 && l < l_top &&
        entrance.remaining[l + 1] * std::max(2.0, u_max) < q.sum_tolerance * smallest)
      stop = l;
    if (stop >= 0 && l >= stop + q.extra_partial_waves) break;
  }
  const double scale = std::sqrt(4.0 * kPi) / (entrance.k_fm * channel.k_fm * q.c_i * q.c_f);
  for (const cplx &s : out.sum) out.O.push_back(scale * s);
  return out;
}

double ZeroRangeMomentumAmplitude(const RadialOverlap &overlap, double q) {
  const int n = static_cast<int>(overlap.u.size()), first = FirstNode(overlap);
  std::vector<double> f(n - first, 0.0);
  for (int p = first; p < n; ++p) {
    const double r = p * overlap.step_fm, qr = q * r;
    const double j0 = qr < 1.0e-6 ? 1.0 - qr * qr / 6.0 : std::sin(qr) / qr;
    f[p - first] = r * overlap.u[p] * j0;
  }
  return 4.0 * kPi * Simpson(f, overlap.step_fm);
}

}  // namespace transfer
}  // namespace physics
}  // namespace azure2
