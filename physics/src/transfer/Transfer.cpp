#include "azure2/physics/transfer/Transfer.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "azure2/physics/AngularCoupling.h"
#include "azure2/physics/Constants.h"
#include "azure2/physics/transfer/RadialNumerics.h"

namespace azure2 {
namespace physics {
namespace transfer {

namespace {

using cplx = std::complex<double>;
const double kPi = std::acos(-1.0);

/// i^n for any integer n.
cplx IPower(int n) {
  static const cplx table[4] = {{1.0, 0.0}, {0.0, 1.0}, {-1.0, 0.0}, {0.0, -1.0}};
  return table[((n % 4) + 4) % 4];
}

/// Y_lm(theta, 0) from a table of SphericalHarmonicsAtPhiZero (stride m_max + 1), any sign of m.
double Ylm(const std::vector<double> &table, int stride, int l, int m) {
  const double y = table[l * stride + std::abs(m)];
  return (m < 0 && (m % 2 != 0)) ? -y : y;
}

/// The highest l solved in a set of waves.
int MaxL(const DistortedWaves &w) {
  int l_max = 0;
  for (const PartialWave &pw : w.waves) l_max = std::max(l_max, pw.l);
  return l_max;
}

/// u_l of a spinless set of waves (one per l), checked.
const std::vector<cplx> &Wave(const DistortedWaves &w, int l) {
  const PartialWave &pw = WaveOf(w, l, 2 * l);
  if (pw.u.empty()) throw std::invalid_argument("transfer: the distorted waves were solved without keeping u");
  return pw.u;
}

void CheckReaction(const TransferReaction &r, const DistortedWaves &in, const DistortedWaves &out) {
  if (!(r.mass_a_u > 0.0) || !(r.mass_A_u > 0.0) || !(r.mass_s_u > 0.0) || !(r.mass_F_u > 0.0))
    throw std::invalid_argument("transfer: masses must be positive");
  if (!(r.mass_s_u < r.mass_a_u) || !(r.mass_A_u < r.mass_F_u))
    throw std::invalid_argument("transfer: m_s < m_a and m_A < m_F");
  if (in.spin_orbit || out.spin_orbit)
    throw std::invalid_argument("transfer: distorted waves with a spin-orbit term are not implemented");
  const Overlap *ov[2] = {&r.projectile, &r.residual};
  for (const Overlap *o : ov) {
    if (!triangle(2 * o->state.l, r.two_s_x, o->state.two_j) ||
        !triangle(o->state.two_j, o->two_J_core, o->two_J_composite))
      throw std::invalid_argument("transfer: an overlap's angular momenta do not couple");
  }
}

/*!
 * sum over M_a, M_A, M_s, M_F of |T|^2 (Eq. (6.5) of the note) from the orbital amplitudes O[(m_F + l_F)(2 l_a + 1) +
 * m_a + l_a]: the overlaps ((l, s_x) j, J_core; J_composite) contracted over the projection of s_x.
 */
double SpinSum(const TransferReaction &r, const std::vector<cplx> &O) {
  const Overlap &pa = r.projectile, &pf = r.residual;
  const int la = pa.state.l, lF = pf.state.l, ja = pa.state.two_j, jF = pf.state.two_j;
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

std::vector<cplx> ZeroRangeOrbitalIntegral(const BoundState &state, const DistortedWaves &in, const DistortedWaves &out,
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
  if (r.projectile.state.l != 0) throw std::invalid_argument("zero range: the projectile overlap must have l = 0");
  const double alpha = r.mass_A_u / r.mass_F_u;
  std::vector<double> xs;
  for (double th : theta_deg) {
    std::vector<cplx> O = ZeroRangeOrbitalIntegral(r.residual.state, in, out, 1.0, alpha, th);
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
  const BoundState &sa = r.projectile.state, &sF = r.residual.state;
  const int la = sa.l, lF = sF.l;
  const int L_lo = std::abs(lF - la), L_hi = lF + la, nL = L_hi - L_lo + 1, Mmax = L_hi, nM = 2 * Mmax + 1;
  const int li_max = MaxL(in), lf_max = MaxL(out);
  const double h = grid.step_fm;
  auto stride_of = [&](const DistortedWaves &w) {
    const double s = h / w.step_fm;
    const int k = static_cast<int>(std::lround(s));
    if (k < 1 || std::fabs(s - k) > 1e-9) throw std::invalid_argument("transfer: the grid step must be a multiple of the waves' steps");
    return k;
  };
  const int si = stride_of(in), sf = stride_of(out);
  const int n = static_cast<int>(std::lround(grid.r_max_fm / h)) + 1;
  if ((n - 1) * si >= static_cast<int>(Wave(in, 0).size()) || (n - 1) * sf >= static_cast<int>(Wave(out, 0).size()))
    throw std::invalid_argument("transfer: r_max beyond the matching radius of the waves");

  // V_sx u_a on the projectile overlap's grid.
  std::vector<double> vu(sa.u.size());
  const bool nuclear = r.interaction == PostInteraction::NuclearPart;
  for (std::size_t k = 0; k < vu.size(); ++k)
    vu[k] = (sa.potential_MeV[k] - (nuclear ? sa.coulomb_MeV[k] : 0.0)) * sa.u[k];

  // <l_F m1 l_a m2 | L M>, by L, m1, m2.
  auto cg_index = [&](int iL, int m1, int m2) { return (iL * (2 * lF + 1) + m1 + lF) * (2 * la + 1) + m2 + la; };
  std::vector<double> cg_vertex(nL * (2 * lF + 1) * (2 * la + 1), 0.0);
  for (int iL = 0; iL < nL; ++iL)
    for (int m1 = -lF; m1 <= lF; ++m1)
      for (int m2 = -la; m2 <= la; ++m2)
        cg_vertex[cg_index(iL, m1, m2)] = clebsch_gordan(2 * lF, 2 * m1, 2 * la, 2 * m2, 2 * (L_lo + iL), 2 * (m1 + m2));

  // acc[L][l_f][M][l_i] = sum over the (r_i, r_f) grid of h^2 r_i r_f u_lf(r_f) u_li(r_i) h^L_{lf M}(r_f, r_i),
  // h^L_{lf M} = int dx Y_lf,M(x) G_LM(x) (Eq. (6.4) of the note).
  auto acc_index = [&](int iL, int lf, int M, int li) {
    return ((static_cast<std::size_t>(iL) * (lf_max + 1) + lf) * nM + M + Mmax) * (li_max + 1) + li;
  };
  std::vector<cplx> acc(static_cast<std::size_t>(nL) * (lf_max + 1) * nM * (li_max + 1), 0.0);
  std::vector<double> hker(static_cast<std::size_t>(nL) * nM * (lf_max + 1));
  std::vector<double> gx, gw;
  GaussLegendre(grid.angle_nodes, gx, gw);
  std::vector<double> G(nL * nM);
  const double reach = jac * grid.vertex_range_fm;
  std::vector<const std::vector<cplx> *> wi(li_max + 1), wf(lf_max + 1);
  for (int l = 0; l <= li_max; ++l) wi[l] = &Wave(in, l);
  for (int l = 0; l <= lf_max; ++l) wf[l] = &Wave(out, l);
  std::vector<cplx> ufv(lf_max + 1), uiv(li_max + 1);

  for (int i = 1; i < n; ++i) {
    const double ri = i * h;
    for (int l = 0; l <= li_max; ++l) uiv[l] = (*wi[l])[i * si];
    const int f_lo = std::max(1, static_cast<int>(std::ceil((alpha * ri - reach) / h)));
    const int f_hi = std::min(n - 1, static_cast<int>(std::floor((alpha * ri + reach) / h)));
    for (int f = f_lo; f <= f_hi; ++f) {
      const double rf = f * h;
      const double x_min = (rf * rf + alpha * alpha * ri * ri - reach * reach) / (2.0 * alpha * ri * rf);
      if (x_min >= 1.0) continue;
      const double x_lo = std::max(-1.0, x_min);
      std::fill(hker.begin(), hker.end(), 0.0);
      for (int g = 0; g < grid.angle_nodes; ++g) {
        const double x = x_lo + (1.0 - x_lo) * 0.5 * (gx[g] + 1.0), wx = gw[g] * 0.5 * (1.0 - x_lo);
        const double sn = std::sqrt(std::max(0.0, 1.0 - x * x));
        // r_i along z, r_f = rf (sin, 0, x); R = (r_i - beta r_f)/jac, rho = (r_f - alpha r_i)/jac.
        const double Rx = -beta * rf * sn / jac, Rz = (ri - beta * rf * x) / jac;
        const double px = rf * sn / jac, pz = (rf * x - alpha * ri) / jac;
        const double R = std::max(std::hypot(Rx, Rz), 1e-12), rho = std::max(std::hypot(px, pz), 1e-12);
        if (rho > grid.vertex_range_fm) continue;
        const double radial = InterpolateOnGrid(sF.u, sF.step_fm, R) / R * InterpolateOnGrid(vu, sa.step_fm, rho) / rho;
        if (radial == 0.0) continue;
        const std::vector<double> yR = SphericalHarmonicsAtPhiZero(lF, lF, Rz / R);  // phi_R = pi: (-1)^m
        const std::vector<double> yr = SphericalHarmonicsAtPhiZero(la, la, pz / rho);  // phi_rho = 0
        std::fill(G.begin(), G.end(), 0.0);
        for (int iL = 0; iL < nL; ++iL)
          for (int m1 = -lF; m1 <= lF; ++m1)
            for (int m2 = -la; m2 <= la; ++m2) {
              const double c = cg_vertex[cg_index(iL, m1, m2)];
              if (c == 0.0) continue;
              const double sign = (m1 % 2 != 0) ? -1.0 : 1.0;
              G[iL * nM + m1 + m2 + Mmax] += c * sign * Ylm(yR, lF + 1, lF, m1) * Ylm(yr, la + 1, la, m2) * radial;
            }
        const std::vector<double> yf = SphericalHarmonicsAtPhiZero(lf_max, Mmax, x);
        for (int iL = 0; iL < nL; ++iL) {
          const int L = L_lo + iL;
          for (int M = -L; M <= L; ++M) {
            const double gv = G[iL * nM + M + Mmax] * wx;
            if (gv == 0.0) continue;
            double *hk = &hker[(static_cast<std::size_t>(iL) * nM + M + Mmax) * (lf_max + 1)];
            for (int lf = std::abs(M); lf <= lf_max; ++lf) hk[lf] += gv * Ylm(yf, Mmax + 1, lf, M);
          }
        }
      }
      const double weight = h * h * ri * rf;
      for (int l = 0; l <= lf_max; ++l) ufv[l] = (*wf[l])[f * sf] * weight;
      for (int iL = 0; iL < nL; ++iL) {
        const int L = L_lo + iL;
        for (int M = -L; M <= L; ++M) {
          const double *hk = &hker[(static_cast<std::size_t>(iL) * nM + M + Mmax) * (lf_max + 1)];
          for (int lf = std::abs(M); lf <= lf_max; ++lf) {
            if (hk[lf] == 0.0) continue;
            const cplx hv = hk[lf] * ufv[lf];
            for (int li = std::abs(lf - L); li <= std::min(li_max, lf + L); ++li)
              if ((li + lf + la + lF) % 2 == 0) acc[acc_index(iL, lf, M, li)] += hv * uiv[li];
          }
        }
      }
    }
  }

  // I^L_{lf li} (Eq. (6.4)), then O_{m_F m_a}(theta) (Eq. (6.3)) and the cross section.
  std::vector<cplx> I(static_cast<std::size_t>(nL) * (lf_max + 1) * (li_max + 1), 0.0);
  auto I_index = [&](int iL, int lf, int li) { return (static_cast<std::size_t>(iL) * (lf_max + 1) + lf) * (li_max + 1) + li; };
  for (int iL = 0; iL < nL; ++iL) {
    const int L = L_lo + iL;
    for (int lf = 0; lf <= lf_max; ++lf)
      for (int li = std::abs(lf - L); li <= std::min(li_max, lf + L); ++li) {
        if ((li + lf + la + lF) % 2 != 0) continue;
        cplx s = 0.0;
        for (int M = -std::min(L, lf); M <= std::min(L, lf); ++M)
          s += clebsch_gordan(2 * lf, 2 * M, 2 * li, 0, 2 * L, 2 * M) * acc[acc_index(iL, lf, M, li)];
        I[I_index(iL, lf, li)] = 8.0 * kPi * kPi / (2 * L + 1.0) * std::sqrt((2 * li + 1.0) / (4.0 * kPi)) * s;
      }
  }
  const double pref = 16.0 * kPi * kPi / (in.k_fm * out.k_fm * jac * jac * jac);
  std::vector<double> xs;
  for (double th : theta_deg) {
    const std::vector<double> yk = SphericalHarmonicsAtPhiZero(lf_max, Mmax, std::cos(th * kPi / 180.0));
    std::vector<cplx> O((2 * lF + 1) * (2 * la + 1), 0.0);
    for (int mF = -lF; mF <= lF; ++mF)
      for (int ma = -la; ma <= la; ++ma) {
        const int M = ma - mF;
        cplx o = 0.0;
        for (int iL = 0; iL < nL; ++iL) {
          const int L = L_lo + iL;
          if (std::abs(M) > L) continue;
          const double cv = clebsch_gordan(2 * lF, -2 * mF, 2 * la, 2 * ma, 2 * L, 2 * M);
          if (cv == 0.0) continue;
          for (int lf = std::abs(M); lf <= lf_max; ++lf)
            for (int li = std::abs(lf - L); li <= std::min(li_max, lf + L); ++li) {
              if ((li + lf + la + lF) % 2 != 0) continue;
              o += cv * IPower(li - lf) * std::exp(cplx(0.0, in.sigma[li] + out.sigma[lf])) *
                   std::sqrt((2 * li + 1.0) / (4.0 * kPi)) * Ylm(yk, Mmax + 1, lf, M) *
                   clebsch_gordan(2 * lf, 2 * M, 2 * li, 0, 2 * L, 2 * M) * I[I_index(iL, lf, li)];
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
