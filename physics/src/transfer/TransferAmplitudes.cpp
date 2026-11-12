#include "azure2/physics/transfer/TransferAmplitudes.h"

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

/// 2j of the distorted waves the wave at position k stands for: its own j with a spin-orbit term; every j of
/// (l, s) without one (one radial wave for all j: the sums over j of Eq. (19.5) are then the identity (19.7)).
std::vector<int> TwoJsOf(const DistortedWaves &w, std::size_t k, int two_s) {
  const PartialWave &pw = w.waves[k];
  if (w.spin_orbit) return {pw.two_j};
  std::vector<int> js;
  for (int two_j = std::abs(2 * pw.l - two_s); two_j <= 2 * pw.l + two_s; two_j += 2) js.push_back(two_j);
  return js;
}

/// Position of the projection set (m_F, m_a, m', m, mu', mu) (projections of spins as 2m).
std::size_t SetIndex(const AngularCoefficients &c, int mF, int ma, int two_mp, int two_m, int two_mup, int two_mu) {
  std::size_t k = static_cast<std::size_t>(mF + c.lF);
  k = k * (2 * c.la + 1) + (ma + c.la);
  k = k * (c.two_Ja + 1) + (two_mp + c.two_Ja) / 2;
  k = k * (c.two_Ja + 1) + (two_m + c.two_Ja) / 2;
  k = k * (c.two_Js + 1) + (two_mup + c.two_Js) / 2;
  return k * (c.two_Js + 1) + (two_mu + c.two_Js) / 2;
}

std::size_t SetCount(const AngularCoefficients &c) {
  return static_cast<std::size_t>(2 * c.lF + 1) * (2 * c.la + 1) * (c.two_Ja + 1) * (c.two_Ja + 1) * (c.two_Js + 1) *
         (c.two_Js + 1);
}

/// <m | S_y | m~> of spin s (2s = two_s), times i: (S_+ - S_-)/2.
double ISy(int two_s, int two_m, int two_mt) {
  const double s = 0.5 * two_s, mt = 0.5 * two_mt;
  if (two_m == two_mt + 2) return 0.5 * std::sqrt(s * (s + 1.0) - mt * (mt + 1.0));
  if (two_m == two_mt - 2) return -0.5 * std::sqrt(s * (s + 1.0) - mt * (mt - 1.0));
  return 0.0;
}

}  // namespace

AngularCoefficients AngularCoefficientsOf(const TransferReaction &r, const DistortedWaves &in, const DistortedWaves &out,
                                          const KernelIntegrals &K, double prefactor) {
  AngularCoefficients c;
  c.la = r.projectile.radial.l;
  c.lF = r.residual.radial.l;
  c.two_Ja = r.projectile.two_J_composite;
  c.two_Js = r.projectile.two_J_core;
  if (in.spin_orbit && in.two_s != c.two_Ja)
    throw std::invalid_argument("transfer: the entrance waves' spin is not that of the projectile overlap (J_a)");
  if (out.spin_orbit && out.two_s != c.two_Js)
    throw std::invalid_argument("transfer: the exit waves' spin is not that of the projectile's core (J_s)");
  c.lf_max = MaxL(out);
  c.n_max = c.lf_max;
  const int nn = 2 * c.n_max + 1;
  c.C.assign(SetCount(c) * (c.lf_max + 1) * nn, 0.0);
  const int la = c.la, lF = c.lF, Ja = c.two_Ja, Js = c.two_Js;
  for (int iL = 0; iL < K.nL; ++iL) {
    const int L = K.L_lo + iL;
    for (int kf = 0; kf < K.n_exit; ++kf)
      for (int ki = 0; ki < K.n_entrance; ++ki) {
        const cplx I = K.I[K.Index(iL, kf, ki)];
        if (I == 0.0) continue;
        const int lf = out.waves[kf].l, li = in.waves[ki].l;
        // Eqs. (19.4)-(19.5): the phase of the partial waves and R^L = (-1)^l_i sqrt((2L+1)/(2l_f+1)) I^L.
        const cplx base = prefactor * I * IPower(li - lf) * std::exp(cplx(0.0, in.sigma[li] + out.sigma[lf])) *
                          std::sqrt((2 * li + 1.0) / (4.0 * kPi)) * ((li % 2 != 0) ? -1.0 : 1.0) *
                          std::sqrt((2 * L + 1.0) / (2 * lf + 1.0));
        for (int ji : TwoJsOf(in, ki, Ja))
          for (int jf : TwoJsOf(out, kf, Js))
            for (int m = -Ja; m <= Ja; m += 2) {
              const double e0 = clebsch_gordan(2 * li, 0, Ja, m, ji, m);
              if (e0 == 0.0) continue;
              for (int mp = -Ja; mp <= Ja; mp += 2) {
                const int mi2 = m - mp;  // 2 m_i: the orbital projection of the entrance wave
                if (std::abs(mi2) > 2 * li) continue;
                const double ei = e0 * clebsch_gordan(2 * li, mi2, Ja, mp, ji, m);
                if (ei == 0.0) continue;
                for (int mF = -lF; mF <= lF; ++mF)
                  for (int ma = -la; ma <= la; ++ma) {
                    const int M = ma - mF;
                    if (std::abs(M) > L) continue;
                    const double cv = ((mF % 2 != 0) ? -1.0 : 1.0) * clebsch_gordan(2 * lF, -2 * mF, 2 * la, 2 * ma, 2 * L, 2 * M);
                    if (cv == 0.0) continue;
                    const int np2 = mi2 + 2 * M;  // 2 n': the orbital projection of the exit wave at the vertex
                    if (std::abs(np2) > 2 * lf) continue;
                    const double cr = clebsch_gordan(2 * li, mi2, 2 * L, 2 * M, 2 * lf, np2);
                    if (cr == 0.0) continue;
                    for (int mup = -Js; mup <= Js; mup += 2) {
                      const int N2 = np2 + mup;
                      if (std::abs(N2) > jf) continue;
                      const double f1 = clebsch_gordan(2 * lf, np2, Js, mup, jf, N2);
                      if (f1 == 0.0) continue;
                      for (int mu = -Js; mu <= Js; mu += 2) {
                        const int n2 = N2 - mu;
                        if (std::abs(n2) > 2 * lf) continue;
                        const double f0 = clebsch_gordan(2 * lf, n2, Js, mu, jf, N2);
                        if (f0 == 0.0) continue;
                        const std::size_t k = (SetIndex(c, mF, ma, mp, m, mup, mu) * (c.lf_max + 1) + lf) * nn +
                                              n2 / 2 + c.n_max;
                        c.C[k] += base * (ei * cv * cr * f1 * f0);
                      }
                    }
                  }
              }
            }
      }
  }
  return c;
}

AmplitudeTable AmplitudesAt(const TransferReaction &r, const AngularCoefficients &c, double theta_deg) {
  const int nn = 2 * c.n_max + 1;
  const std::vector<double> y = SphericalHarmonicsAtPhiZero(c.lf_max, c.n_max, std::cos(theta_deg * kPi / 180.0));
  // Q of Eq. (19.5) for every projection set.
  std::vector<cplx> Q(SetCount(c), 0.0);
  for (std::size_t s = 0; s < Q.size(); ++s) {
    cplx q = 0.0;
    for (int lf = 0; lf <= c.lf_max; ++lf)
      for (int n = -lf; n <= lf; ++n) {
        const cplx v = c.C[(s * (c.lf_max + 1) + lf) * nn + n + c.n_max];
        if (v != 0.0) q += v * Ylm(y, c.n_max + 1, lf, n);
      }
    Q[s] = q;
  }
  // Eq. (19.6): the overlaps ((l, s_x) j, J_core; J_composite) at the vertex, with the projections m' of a and mu'
  // of s there, contracted over the projection sigma of s_x.
  const Overlap &pa = r.projectile, &pf = r.residual;
  const int la = c.la, lF = c.lF, ja = pa.radial.two_j, jF = pf.radial.two_j, sx = r.two_s_x;
  AmplitudeTable t;
  t.two_Ja = pa.two_J_composite, t.two_Js = pa.two_J_core, t.two_JA = pf.two_J_core, t.two_JF = pf.two_J_composite;
  t.T.assign(static_cast<std::size_t>(t.two_Js + 1) * (t.two_JF + 1) * (t.two_Ja + 1) * (t.two_JA + 1), 0.0);
  const double amplitude = pa.amplitude * pf.amplitude;
  for (int MF = -t.two_JF; MF <= t.two_JF; MF += 2)
    for (int MA = -t.two_JA; MA <= t.two_JA; MA += 2) {
      const int muF = MF - MA;
      if (std::abs(muF) > jF) continue;
      const double cF = clebsch_gordan(jF, muF, t.two_JA, MA, t.two_JF, MF);
      if (cF == 0.0) continue;
      for (int mp = -t.two_Ja; mp <= t.two_Ja; mp += 2)
        for (int mup = -t.two_Js; mup <= t.two_Js; mup += 2) {
          const int mua = mp - mup;
          if (std::abs(mua) > ja) continue;
          const double ca = clebsch_gordan(ja, mua, t.two_Js, mup, t.two_Ja, mp);
          if (ca == 0.0) continue;
          for (int sig = -sx; sig <= sx; sig += 2) {
            const int ma2 = mua - sig, mF2 = muF - sig;  // 2 m_a, 2 m_F
            if (std::abs(ma2) > 2 * la || std::abs(mF2) > 2 * lF) continue;
            const double v = amplitude * ca * cF * clebsch_gordan(2 * la, ma2, sx, sig, ja, mua) *
                             clebsch_gordan(2 * lF, mF2, sx, sig, jF, muF);
            if (v == 0.0) continue;
            for (int m = -t.two_Ja; m <= t.two_Ja; m += 2)
              for (int mu = -t.two_Js; mu <= t.two_Js; mu += 2)
                t.T[t.Index(mu, MF, m, MA)] += v * Q[SetIndex(c, mF2 / 2, ma2 / 2, mp, m, mup, mu)];
          }
        }
    }
  return t;
}

double CrossSectionOf(const TransferReaction &r, const DistortedWaves &in, const DistortedWaves &out,
                      const AmplitudeTable &t) {
  double sum = 0.0;
  for (const cplx &v : t.T) sum += std::norm(v);
  const double mu_i = r.mass_a_u * r.mass_A_u / (r.mass_a_u + r.mass_A_u) * kAtomicMassUnit_MeV;
  const double mu_f = r.mass_s_u * r.mass_F_u / (r.mass_s_u + r.mass_F_u) * kAtomicMassUnit_MeV;
  const double hc2 = kHbarC_MeV_fm * kHbarC_MeV_fm;
  const double factor = mu_i * mu_f / std::pow(2.0 * kPi * hc2, 2) * (out.k_fm / in.k_fm) /
                        ((t.two_Ja + 1.0) * (t.two_JA + 1.0));
  return 10.0 * factor * sum;  // Eq. (6.6), fm^2 -> mb
}

double VectorAnalysingPowerOf(const AmplitudeTable &t) {
  if (t.two_Ja == 0) return 0.0;
  const double s = 0.5 * t.two_Ja;
  double num = 0.0, den = 0.0;
  for (int mu = -t.two_Js; mu <= t.two_Js; mu += 2)
    for (int MF = -t.two_JF; MF <= t.two_JF; MF += 2)
      for (int MA = -t.two_JA; MA <= t.two_JA; MA += 2)
        for (int m = -t.two_Ja; m <= t.two_Ja; m += 2) {
          const cplx a = t.T[t.Index(mu, MF, m, MA)];
          den += std::norm(a);
          for (int mt = -t.two_Ja; mt <= t.two_Ja; mt += 2) {
            const double isy = ISy(t.two_Ja, m, mt);
            if (isy != 0.0) num += (a * std::conj(t.T[t.Index(mu, MF, mt, MA)]) * cplx(0.0, -isy)).real();
          }
        }
  // Eq. (19.8): iT11 = sqrt(3/(2 s (s+1))) Tr(T S_y T^+) / Tr(T T^+).
  return den > 0.0 ? std::sqrt(3.0 / (2.0 * s * (s + 1.0))) * num / den : 0.0;
}

}  // namespace transfer
}  // namespace physics
}  // namespace azure2
