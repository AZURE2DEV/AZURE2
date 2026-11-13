#include "azure2/physics/transfer/Transfer.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "azure2/physics/AngularCoupling.h"
#include "azure2/physics/Constants.h"
#include "azure2/physics/transfer/RadialNumerics.h"
#include "azure2/physics/transfer/TransferAmplitudes.h"
#include "azure2/physics/transfer/TransferKernel.h"

namespace azure2 {
namespace physics {
namespace transfer {

namespace {

using cplx = std::complex<double>;
const double kPi = std::acos(-1.0);

void CheckReaction(const TransferReaction &r, const DistortedWaves &in, const DistortedWaves &out) {
  if (!(r.mass_a_u > 0.0) || !(r.mass_A_u > 0.0) || !(r.mass_s_u > 0.0) || !(r.mass_F_u > 0.0))
    throw std::invalid_argument("transfer: masses must be positive");
  if (!(r.mass_s_u < r.mass_a_u) || !(r.mass_A_u < r.mass_F_u))
    throw std::invalid_argument("transfer: m_s < m_a and m_A < m_F");
  const Overlap *ov[2] = {&r.projectile, &r.residual};
  for (const Overlap *o : ov) {
    if (!triangle(2 * o->radial.l, r.two_s_x, o->radial.two_j) ||
        !triangle(o->radial.two_j, o->two_J_core, o->two_J_composite))
      throw std::invalid_argument("transfer: an overlap's angular momenta do not couple");
    if (o->radial.u.empty() || !(o->radial.step_fm > 0.0)) throw std::invalid_argument("transfer: an empty overlap");
  }
  (void)in;
  (void)out;
}

/// R_{wf wi} = int dR u_wf(c_f R) u_F(R) u_wi(c_i R) / R (fm^-1/2 times the unit of u_F) on the residual overlap's
/// grid, from its first node (r_min of a cut tail) as far as both waves reach, for every pair of waves (Eq. (6.7)).
std::vector<cplx> ZeroRangeRadialIntegrals(const RadialOverlap &state, const DistortedWaves &in,
                                           const DistortedWaves &out, double c_i, double c_f) {
  const double hb = state.step_fm;
  const double r_in = (Wave(in, 0).size() - 1) * in.step_fm, r_out = (Wave(out, 0).size() - 1) * out.step_fm;
  const int n = std::min<int>(static_cast<int>(state.u.size()),
                              static_cast<int>(std::floor(std::min(r_in / c_i, r_out / c_f) / hb)) + 1);
  const std::size_t nwi = in.waves.size(), nwf = out.waves.size();
  std::vector<std::vector<cplx>> ui(nwi), uf(nwf);
  for (std::size_t k = 0; k < nwi; ++k) {
    ui[k].resize(n);
    for (int p = 0; p < n; ++p) ui[k][p] = InterpolateOnGrid(Wave(in, k), in.step_fm, c_i * p * hb);
  }
  for (std::size_t k = 0; k < nwf; ++k) {
    uf[k].resize(n);
    for (int p = 0; p < n; ++p) uf[k][p] = InterpolateOnGrid(Wave(out, k), out.step_fm, c_f * p * hb);
  }
  const int first = FirstNode(state);  // 0, or the node of a cut tail's r_min
  std::vector<cplx> R(nwf * nwi, 0.0), f(n);
  for (std::size_t kf = 0; kf < nwf; ++kf)
    for (std::size_t ki = 0; ki < nwi; ++ki) {
      const int lf = out.waves[kf].l, li = in.waves[ki].l;
      if ((li + lf + state.l) % 2 != 0 || li < std::abs(lf - state.l) || li > lf + state.l) continue;
      f[0] = 0.0;
      for (int p = 1; p < n; ++p) f[p] = uf[kf][p] * state.u[p] * ui[ki][p] / (p * hb);
      R[kf * nwi + ki] = Simpson(std::vector<cplx>(f.begin() + first, f.end()), hb) / (c_i * c_f);
    }
  return R;
}

/// sqrt((2l_f+1)(2l_i+1)/(4 pi (2l+1))) <l_f 0 l_i 0 | l 0>: the angular integral of Eq. (6.7).
double Gaunt0(int lf, int li, int l) {
  return std::sqrt((2 * lf + 1.0) * (2 * li + 1.0) / (4.0 * kPi * (2 * l + 1.0))) *
         clebsch_gordan(2 * lf, 0, 2 * li, 0, 2 * l, 0);
}

/// The observables at each angle from the radial integrals and the prefactor of Eq. (19.5).
TransferObservables ObservablesOf(const TransferReaction &r, const DistortedWaves &in, const DistortedWaves &out,
                                  const KernelIntegrals &K, double prefactor, const std::vector<double> &theta_deg) {
  const AngularCoefficients c = AngularCoefficientsOf(r, in, out, K, prefactor);
  TransferObservables obs;
  for (double th : theta_deg) {
    const AmplitudeTable t = AmplitudesAt(r, c, th);
    obs.cross_section.push_back(CrossSectionOf(r, in, out, t));
    obs.analysing_power_iT11.push_back(VectorAnalysingPowerOf(t));
  }
  return obs;
}

}  // namespace

std::vector<cplx> ZeroRangeOrbitalIntegral(const RadialOverlap &state, const DistortedWaves &in, const DistortedWaves &out,
                                           double c_i, double c_f, double theta_deg) {
  if (!(c_i > 0.0) || !(c_f > 0.0)) throw std::invalid_argument("zero range: scale factors must be positive");
  if (in.spin_orbit || out.spin_orbit)
    throw std::invalid_argument("zero range orbital integral: waves without a spin-orbit term (one per l)");
  const int lF = state.l;
  const int li_max = MaxL(in), lf_max = MaxL(out);
  const std::vector<cplx> radial = ZeroRangeRadialIntegrals(state, in, out, c_i, c_f);
  const std::vector<double> y = SphericalHarmonicsAtPhiZero(lf_max, lF, std::cos(theta_deg * kPi / 180.0));
  std::vector<cplx> O(2 * lF + 1, 0.0);
  for (int lf = 0; lf <= lf_max; ++lf)
    for (int li = std::abs(lf - lF); li <= std::min(li_max, lf + lF); ++li) {
      if ((li + lf + lF) % 2 != 0) continue;
      const cplx phase = IPower(li - lf) * std::exp(cplx(0.0, in.sigma[li] + out.sigma[lf])) *
                         std::sqrt((2 * li + 1.0) / (4.0 * kPi));
      for (int m = -lF; m <= lF; ++m) {
        if (std::abs(m) > lf) continue;
        O[m + lF] += phase * Ylm(y, lF + 1, lf, m) * Gaunt0(lf, li, lF) *
                     clebsch_gordan(2 * lf, 2 * m, 2 * li, 0, 2 * lF, 2 * m) * radial[lf * (li_max + 1) + li];
      }
    }
  const double pref = 16.0 * kPi * kPi / (in.k_fm * out.k_fm);
  for (cplx &o : O) o *= pref;
  return O;
}

TransferObservables ZeroRangeObservables(const TransferReaction &r, double d0, const DistortedWaves &in,
                                         const DistortedWaves &out, const std::vector<double> &theta_deg) {
  CheckReaction(r, in, out);
  if (r.projectile.radial.l != 0) throw std::invalid_argument("zero range: the projectile overlap must have l = 0");
  const double alpha = r.mass_A_u / r.mass_F_u;
  const int lF = r.residual.radial.l;
  // The zero range as the kernel table of Eq. (6.4) with l_a = 0, L = l_F (section 19): I = D0 g R.
  const std::vector<cplx> radial = ZeroRangeRadialIntegrals(r.residual.radial, in, out, 1.0, alpha);
  KernelIntegrals K;
  K.L_lo = lF;
  K.nL = 1;
  K.n_exit = static_cast<int>(out.waves.size());
  K.n_entrance = static_cast<int>(in.waves.size());
  K.I.assign(radial.size(), 0.0);
  for (int kf = 0; kf < K.n_exit; ++kf)
    for (int ki = 0; ki < K.n_entrance; ++ki) {
      const cplx v = radial[static_cast<std::size_t>(kf) * K.n_entrance + ki];
      if (v != 0.0) K.I[K.Index(0, kf, ki)] = d0 * Gaunt0(out.waves[kf].l, in.waves[ki].l, lF) * v;
    }
  return ObservablesOf(r, in, out, K, 16.0 * kPi * kPi / (in.k_fm * out.k_fm), theta_deg);
}

std::vector<double> ZeroRangeCrossSection(const TransferReaction &r, double d0, const DistortedWaves &in,
                                          const DistortedWaves &out, const std::vector<double> &theta_deg) {
  return ZeroRangeObservables(r, d0, in, out, theta_deg).cross_section;
}

TransferObservables FiniteRangeObservables(const TransferReaction &r, const DistortedWaves &in,
                                           const DistortedWaves &out, const TransferGrid &grid,
                                           const std::vector<double> &theta_deg) {
  CheckReaction(r, in, out);
  const double alpha = r.mass_A_u / r.mass_F_u, beta = r.mass_s_u / r.mass_a_u, jac = 1.0 - alpha * beta;
  const KernelIntegrals K = FiniteRangeIntegrals(r, in, out, grid);
  return ObservablesOf(r, in, out, K, 16.0 * kPi * kPi / (in.k_fm * out.k_fm * jac * jac * jac), theta_deg);
}

std::vector<double> FiniteRangeCrossSection(const TransferReaction &r, const DistortedWaves &in,
                                            const DistortedWaves &out, const TransferGrid &grid,
                                            const std::vector<double> &theta_deg) {
  return FiniteRangeObservables(r, in, out, grid, theta_deg).cross_section;
}

}  // namespace transfer
}  // namespace physics
}  // namespace azure2
