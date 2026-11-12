#include "azure2/physics/transfer/TransferKernel.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "azure2/physics/AngularCoupling.h"
#include "azure2/physics/transfer/RadialNumerics.h"

namespace azure2 {
namespace physics {
namespace transfer {

namespace {
using cplx = std::complex<double>;
const double kPi = std::acos(-1.0);
}  // namespace

cplx IPower(int n) {
  static const cplx table[4] = {{1.0, 0.0}, {0.0, 1.0}, {-1.0, 0.0}, {0.0, -1.0}};
  return table[((n % 4) + 4) % 4];
}

int MaxL(const DistortedWaves &w) {
  int l_max = 0;
  for (const PartialWave &pw : w.waves) l_max = std::max(l_max, pw.l);
  return l_max;
}

const std::vector<cplx> &Wave(const DistortedWaves &w, std::size_t k) {
  const PartialWave &pw = w.waves.at(k);
  if (pw.u.empty()) throw std::invalid_argument("transfer: the distorted waves were solved without keeping u");
  return pw.u;
}

double Ylm(const std::vector<double> &table, int stride, int l, int m) {
  const double y = table[l * stride + std::abs(m)];
  return (m < 0 && (m % 2 != 0)) ? -y : y;
}

KernelIntegrals FiniteRangeIntegrals(const TransferReaction &r, const DistortedWaves &in, const DistortedWaves &out,
                                     const TransferGrid &grid) {
  const double alpha = r.mass_A_u / r.mass_F_u, beta = r.mass_s_u / r.mass_a_u, jac = 1.0 - alpha * beta;
  const RadialOverlap &sa = r.projectile.radial, &sF = r.residual.radial;
  const int la = sa.l, lF = sF.l;
  KernelIntegrals K;
  K.L_lo = std::abs(lF - la);
  K.nL = lF + la - K.L_lo + 1;
  K.n_exit = static_cast<int>(out.waves.size());
  K.n_entrance = static_cast<int>(in.waves.size());
  const int L_lo = K.L_lo, nL = K.nL, Mmax = lF + la, nM = 2 * Mmax + 1, lf_max = MaxL(out), li_max = MaxL(in);
  const int nwf = K.n_exit, nwi = K.n_entrance;
  // The waves of each l: positions first[l] ... first[l + 1] - 1 (one per l without a spin-orbit term, one per j with).
  auto first_of = [](const DistortedWaves &w, int l_max) {
    std::vector<int> first(l_max + 2, 0);
    for (const PartialWave &pw : w.waves) ++first[pw.l + 1];
    for (int l = 0; l <= l_max; ++l) first[l + 1] += first[l];
    for (std::size_t k = 0; k < w.waves.size(); ++k)
      if (k > 0 && w.waves[k].l < w.waves[k - 1].l) throw std::invalid_argument("transfer: waves not ordered by l");
    return first;
  };
  const std::vector<int> first_f = first_of(out, lf_max), first_i = first_of(in, li_max);
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

  // The interaction (Eq. (17.1)): one binding term, V u tabulated on its overlap's grid; the other terms are
  // evaluated at each node and multiply u_F(r_xA) u_a(r_sx).
  std::vector<double> vu_a, vu_F;
  std::vector<const InteractionTerm *> others;
  for (const InteractionTerm &t : r.interaction) {
    if (t.kind == InteractionKind::ProjectileBinding || t.kind == InteractionKind::ProjectileBindingNuclear) {
      if (!vu_a.empty() || !vu_F.empty()) throw std::invalid_argument("transfer: more than one binding term");
      vu_a = BindingSource(t, sa);
    } else if (t.kind == InteractionKind::ResidualBinding) {
      if (!vu_a.empty() || !vu_F.empty()) throw std::invalid_argument("transfer: more than one binding term");
      vu_F = BindingSource(t, sF);
    } else {
      others.push_back(&t);
    }
  }
  if (vu_a.empty() && vu_F.empty()) throw std::invalid_argument("transfer: the interaction has no binding term");

  // <l_F m1 l_a m2 | L M>, by L, m1, m2.
  auto cg_index = [&](int iL, int m1, int m2) { return (iL * (2 * lF + 1) + m1 + lF) * (2 * la + 1) + m2 + la; };
  std::vector<double> cg_vertex(nL * (2 * lF + 1) * (2 * la + 1), 0.0);
  for (int iL = 0; iL < nL; ++iL)
    for (int m1 = -lF; m1 <= lF; ++m1)
      for (int m2 = -la; m2 <= la; ++m2)
        cg_vertex[cg_index(iL, m1, m2)] = clebsch_gordan(2 * lF, 2 * m1, 2 * la, 2 * m2, 2 * (L_lo + iL), 2 * (m1 + m2));

  // acc[L][wf][M][wi] = sum over the (r_i, r_f) grid of h^2 r_i r_f u_wf(r_f) u_wi(r_i) h^L_{lf M}(r_f, r_i),
  // h^L_{lf M} = int dx Y_lf,M(x) G_LM(x) (Eq. (6.4) of the note); wf, wi the waves of l_f, l_i.
  auto acc_index = [&](int iL, int wf, int M, int wi) {
    return ((static_cast<std::size_t>(iL) * nwf + wf) * nM + M + Mmax) * nwi + wi;
  };
  std::vector<cplx> acc(static_cast<std::size_t>(nL) * nwf * nM * nwi, 0.0);
  std::vector<cplx> hker(static_cast<std::size_t>(nL) * nM * (lf_max + 1));
  std::vector<double> gx, gw;
  GaussLegendre(grid.angle_nodes, gx, gw);
  std::vector<cplx> G(nL * nM);
  const double reach = jac * grid.vertex_range_fm;
  std::vector<const std::vector<cplx> *> wi(nwi), wf(nwf);
  for (int k = 0; k < nwi; ++k) wi[k] = &Wave(in, k);
  for (int k = 0; k < nwf; ++k) wf[k] = &Wave(out, k);
  std::vector<cplx> ufv(nwf), uiv(nwi);

  for (int i = 1; i < n; ++i) {
    const double ri = i * h;
    for (int k = 0; k < nwi; ++k) uiv[k] = (*wi[k])[i * si];
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
        cplx radial = 0.0;
        if (!vu_a.empty()) radial = InterpolateOnGrid(sF.u, sF.step_fm, R) * InterpolateOnGrid(vu_a, sa.step_fm, rho);
        else radial = InterpolateOnGrid(vu_F, sF.step_fm, R) * InterpolateOnGrid(sa.u, sa.step_fm, rho);
        if (!others.empty()) {
          const NodeDistances d{rho, R, std::hypot(Rx + px, Rz + pz), rf, ri};
          cplx v = 0.0;
          for (const InteractionTerm *t : others) v += TermValue(*t, d);
          radial += InterpolateOnGrid(sF.u, sF.step_fm, R) * InterpolateOnGrid(sa.u, sa.step_fm, rho) * v;
        }
        if (radial == 0.0) continue;
        radial /= R * rho;
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
            const cplx gv = G[iL * nM + M + Mmax] * wx;
            if (gv == 0.0) continue;
            cplx *hk = &hker[(static_cast<std::size_t>(iL) * nM + M + Mmax) * (lf_max + 1)];
            for (int lf = std::abs(M); lf <= lf_max; ++lf) hk[lf] += gv * Ylm(yf, Mmax + 1, lf, M);
          }
        }
      }
      const double weight = h * h * ri * rf;
      for (int k = 0; k < nwf; ++k) ufv[k] = (*wf[k])[f * sf] * weight;
      for (int iL = 0; iL < nL; ++iL) {
        const int L = L_lo + iL;
        for (int M = -L; M <= L; ++M) {
          const cplx *hk = &hker[(static_cast<std::size_t>(iL) * nM + M + Mmax) * (lf_max + 1)];
          for (int lf = std::abs(M); lf <= lf_max; ++lf) {
            if (hk[lf] == 0.0) continue;
            for (int kf = first_f[lf]; kf < first_f[lf + 1]; ++kf) {
              const cplx hv = hk[lf] * ufv[kf];
              for (int li = std::abs(lf - L); li <= std::min(li_max, lf + L); ++li)
                if ((li + lf + la + lF) % 2 == 0)
                  for (int ki = first_i[li]; ki < first_i[li + 1]; ++ki) acc[acc_index(iL, kf, M, ki)] += hv * uiv[ki];
            }
          }
        }
      }
    }
  }

  // I^L of each pair of waves (Eq. (6.4)).
  K.I.assign(static_cast<std::size_t>(nL) * nwf * nwi, 0.0);
  for (int iL = 0; iL < nL; ++iL) {
    const int L = L_lo + iL;
    for (int kf = 0; kf < nwf; ++kf)
      for (int ki = 0; ki < nwi; ++ki) {
        const int lf = out.waves[kf].l, li = in.waves[ki].l;
        if ((li + lf + la + lF) % 2 != 0 || li < std::abs(lf - L) || li > lf + L) continue;
        cplx s = 0.0;
        for (int M = -std::min(L, lf); M <= std::min(L, lf); ++M)
          s += clebsch_gordan(2 * lf, 2 * M, 2 * li, 0, 2 * L, 2 * M) * acc[acc_index(iL, kf, M, ki)];
        K.I[K.Index(iL, kf, ki)] = 8.0 * kPi * kPi / (2 * L + 1.0) * std::sqrt((2 * li + 1.0) / (4.0 * kPi)) * s;
      }
  }
  return K;
}

}  // namespace transfer
}  // namespace physics
}  // namespace azure2
