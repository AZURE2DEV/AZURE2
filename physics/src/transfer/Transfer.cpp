#include "azure2/physics/transfer/Transfer.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "azure2/physics/AngularCoupling.h"
#include "azure2/physics/Constants.h"
#include "azure2/physics/transfer/RadialNumerics.h"
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
  if (in.spin_orbit || out.spin_orbit)
    throw std::invalid_argument("transfer: distorted waves with a spin-orbit term are not implemented");
  const Overlap *ov[2] = {&r.projectile, &r.residual};
  for (const Overlap *o : ov) {
    if (!triangle(2 * o->radial.l, r.two_s_x, o->radial.two_j) ||
        !triangle(o->radial.two_j, o->two_J_core, o->two_J_composite))
      throw std::invalid_argument("transfer: an overlap's angular momenta do not couple");
    if (o->radial.u.empty() || !(o->radial.step_fm > 0.0)) throw std::invalid_argument("transfer: an empty overlap");
  }
}

/*!
 * sum over M_a, M_A, M_s, M_F of |T|^2 (Eq. (6.5) of the note) from the orbital amplitudes O[(m_F + l_F)(2 l_a + 1) +
 * m_a + l_a]: the overlaps ((l, s_x) j, J_core; J_composite) contracted over the projection of s_x.
 */
double SpinSum(const TransferReaction &r, const std::vector<cplx> &O) {
  const Overlap &pa = r.projectile, &pf = r.residual;
  const int la = pa.radial.l, lF = pf.radial.l, ja = pa.radial.two_j, jF = pf.radial.two_j;
  const int sx = r.two_s_x;
  double sum = 0.0;
  for (int Ma = -pa.two_J_composite; Ma <= pa.two_J_composite; Ma += 2)
    for (int Ms = -pa.two_J_core; Ms <= pa.two_J_core; Ms += 2) {
      const int mua = Ma - Ms;
      if (std::abs(mua) > ja) continue;
      const double ca = clebsch_gordan(ja, mua, pa.two_J_core, Ms, pa.two_J_composite, Ma);
      if (ca == 0.0) continue;
      for (int MF = -pf.two_J_composite; MF <= pf.two_J_composite; MF += 2)
        for (int MA = -pf.two_J_core; MA <= pf.two_J_core; MA += 2) {
          const int muF = MF - MA;
          if (std::abs(muF) > jF) continue;
          const double cF = clebsch_gordan(jF, muF, pf.two_J_core, MA, pf.two_J_composite, MF);
          if (cF == 0.0) continue;
          cplx t = 0.0;
          for (int sig = -sx; sig <= sx; sig += 2) {
            const int ma2 = mua - sig, mF2 = muF - sig;  // 2 m_a, 2 m_F
            if (std::abs(ma2) > 2 * la || std::abs(mF2) > 2 * lF) continue;
            const double c = clebsch_gordan(2 * la, ma2, sx, sig, ja, mua) * clebsch_gordan(2 * lF, mF2, sx, sig, jF, muF);
            t += c * O[(mF2 / 2 + lF) * (2 * la + 1) + ma2 / 2 + la];
          }
          sum += std::norm(ca * cF * t);
        }
    }
  return sum * std::pow(pa.amplitude * pf.amplitude, 2);
}

/// Eq. (6.6): mb/sr from sum |T|^2 (MeV^2 fm^6).
double CrossSection(const TransferReaction &r, const DistortedWaves &in, const DistortedWaves &out, double spin_sum) {
  const double mu_i = r.mass_a_u * r.mass_A_u / (r.mass_a_u + r.mass_A_u) * kAtomicMassUnit_MeV;
  const double mu_f = r.mass_s_u * r.mass_F_u / (r.mass_s_u + r.mass_F_u) * kAtomicMassUnit_MeV;
  const double hc2 = kHbarC_MeV_fm * kHbarC_MeV_fm;
  const double factor = mu_i * mu_f / std::pow(2.0 * kPi * hc2, 2) * (out.k_fm / in.k_fm) /
                        ((r.projectile.two_J_composite + 1.0) * (r.residual.two_J_core + 1.0));
  return 10.0 * factor * spin_sum;
}

}  // namespace

std::vector<cplx> ZeroRangeOrbitalIntegral(const RadialOverlap &state, const DistortedWaves &in, const DistortedWaves &out,
                                           double c_i, double c_f, double theta_deg) {
  if (!(c_i > 0.0) || !(c_f > 0.0)) throw std::invalid_argument("zero range: scale factors must be positive");
  const int lF = state.l;
  const int li_max = MaxL(in), lf_max = MaxL(out);
  const double hb = state.step_fm;
  // Radial integrals R_{l_f l_i} = int dR u_lf(c_f R) u_F(R) u_li(c_i R) / R on the bound state's grid, as far as
  // both waves reach.
  const double r_in = (Wave(in, 0).size() - 1) * in.step_fm, r_out = (Wave(out, 0).size() - 1) * out.step_fm;
  const int n = std::min<int>(static_cast<int>(state.u.size()),
                              static_cast<int>(std::floor(std::min(r_in / c_i, r_out / c_f) / hb)) + 1);
  const std::vector<double> y = SphericalHarmonicsAtPhiZero(lf_max, lF, std::cos(theta_deg * kPi / 180.0));
  std::vector<cplx> O(2 * lF + 1, 0.0);
  std::vector<cplx> f(n);
  std::vector<std::vector<cplx>> ui(li_max + 1), uf(lf_max + 1);
  for (int l = 0; l <= li_max; ++l) {
    ui[l].resize(n);
    for (int k = 0; k < n; ++k) ui[l][k] = InterpolateOnGrid(Wave(in, l), in.step_fm, c_i * k * hb);
  }
  for (int l = 0; l <= lf_max; ++l) {
    uf[l].resize(n);
    for (int k = 0; k < n; ++k) uf[l][k] = InterpolateOnGrid(Wave(out, l), out.step_fm, c_f * k * hb);
  }
  for (int lf = 0; lf <= lf_max; ++lf)
    for (int li = std::abs(lf - lF); li <= std::min(li_max, lf + lF); ++li) {
      if ((li + lf + lF) % 2 != 0) continue;
      f[0] = 0.0;
      for (int k = 1; k < n; ++k) f[k] = uf[lf][k] * state.u[k] * ui[li][k] / (k * hb);
      const cplx radial = Simpson(f, hb) / (c_i * c_f);
      const double gaunt0 = std::sqrt((2 * lf + 1.0) * (2 * li + 1.0) / (4.0 * kPi * (2 * lF + 1.0))) *
                            clebsch_gordan(2 * lf, 0, 2 * li, 0, 2 * lF, 0);
      const cplx phase = IPower(li - lf) * std::exp(cplx(0.0, in.sigma[li] + out.sigma[lf])) *
                         std::sqrt((2 * li + 1.0) / (4.0 * kPi));
      for (int m = -lF; m <= lF; ++m) {
        if (std::abs(m) > lf) continue;
        O[m + lF] += phase * Ylm(y, lF + 1, lf, m) * gaunt0 * clebsch_gordan(2 * lf, 2 * m, 2 * li, 0, 2 * lF, 2 * m) * radial;
      }
    }
  const double pref = 16.0 * kPi * kPi / (in.k_fm * out.k_fm);
  for (cplx &o : O) o *= pref;
  return O;
}

std::vector<double> ZeroRangeCrossSection(const TransferReaction &r, double d0, const DistortedWaves &in,
                                          const DistortedWaves &out, const std::vector<double> &theta_deg) {
  CheckReaction(r, in, out);
  if (r.projectile.radial.l != 0) throw std::invalid_argument("zero range: the projectile overlap must have l = 0");
  const double alpha = r.mass_A_u / r.mass_F_u;
  std::vector<double> xs;
  for (double th : theta_deg) {
    std::vector<cplx> O = ZeroRangeOrbitalIntegral(r.residual.radial, in, out, 1.0, alpha, th);
    for (cplx &o : O) o *= d0;
    xs.push_back(CrossSection(r, in, out, SpinSum(r, O)));
  }
  return xs;
}

std::vector<double> FiniteRangeCrossSection(const TransferReaction &r, const DistortedWaves &in,
                                            const DistortedWaves &out, const TransferGrid &grid,
                                            const std::vector<double> &theta_deg) {
  CheckReaction(r, in, out);
  const double alpha = r.mass_A_u / r.mass_F_u, beta = r.mass_s_u / r.mass_a_u, jac = 1.0 - alpha * beta;
  const int la = r.projectile.radial.l, lF = r.residual.radial.l;
  const KernelIntegrals K = FiniteRangeIntegrals(r, in, out, grid);
  const int Mmax = la + lF;
  // O_{m_F m_a}(theta) (Eq. (6.3)) and the cross section.
  const double pref = 16.0 * kPi * kPi / (in.k_fm * out.k_fm * jac * jac * jac);
  std::vector<double> xs;
  for (double th : theta_deg) {
    const std::vector<double> yk = SphericalHarmonicsAtPhiZero(K.lf_max, Mmax, std::cos(th * kPi / 180.0));
    std::vector<cplx> O((2 * lF + 1) * (2 * la + 1), 0.0);
    for (int mF = -lF; mF <= lF; ++mF)
      for (int ma = -la; ma <= la; ++ma) {
        const int M = ma - mF;
        cplx o = 0.0;
        for (int iL = 0; iL < K.nL; ++iL) {
          const int L = K.L_lo + iL;
          if (std::abs(M) > L) continue;
          const double cv = clebsch_gordan(2 * lF, -2 * mF, 2 * la, 2 * ma, 2 * L, 2 * M);
          if (cv == 0.0) continue;
          for (int lf = std::abs(M); lf <= K.lf_max; ++lf)
            for (int li = std::abs(lf - L); li <= std::min(K.li_max, lf + L); ++li) {
              if ((li + lf + la + lF) % 2 != 0) continue;
              o += cv * IPower(li - lf) * std::exp(cplx(0.0, in.sigma[li] + out.sigma[lf])) *
                   std::sqrt((2 * li + 1.0) / (4.0 * kPi)) * Ylm(yk, Mmax + 1, lf, M) *
                   clebsch_gordan(2 * lf, 2 * M, 2 * li, 0, 2 * L, 2 * M) * K.I[K.Index(iL, lf, li)];
            }
        }
        O[(mF + lF) * (2 * la + 1) + ma + la] = ((mF % 2 != 0) ? -1.0 : 1.0) * pref * o;
      }
    xs.push_back(CrossSection(r, in, out, SpinSum(r, O)));
  }
  return xs;
}

}  // namespace transfer
}  // namespace physics
}  // namespace azure2
