#include "azure2/physics/transfer/SourceProjection.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "azure2/physics/AngularCoupling.h"
#include "azure2/physics/transfer/RadialNumerics.h"

namespace azure2 {
namespace physics {
namespace transfer {

namespace {

using cplx = std::complex<double>;

constexpr int kTableMargin = 8;  // nodes of a wave table beyond R_max (the six-point stencil stays inside)
constexpr cplx kIPower[4] = {cplx(1, 0), cplx(0, 1), cplx(-1, 0), cplx(0, -1)};  // i^n

// The waves L = 0 ... l_max on n nodes of step h with their phases e^{i sigma_L}: the library's stored waves (matched
// beyond the turning point; Riccati-Bessel functions for a plane wave).
void WaveTable(const WaveChannel &channel, double step, int n, int l_max, const MatchingBeyondTurningPoint &rule,
               std::vector<std::vector<cplx>> &u, std::vector<cplx> &phase) {
  DistortedWaves w = SolveStoredWaves(channel, step, n, l_max, rule);
  u.assign(l_max + 1, std::vector<cplx>());
  phase.assign(l_max + 1, cplx(1.0, 0.0));
  for (int L = 0; L <= l_max; L++) {
    u[L] = std::move(w.waves[L].u);
    phase[L] = std::polar(1.0, w.sigma[L]);
  }
}

// Gauss-Legendre nodes and weights on [lo, hi] appended to x, w.
void AppendGaussLegendre(int n, double lo, double hi, std::vector<double> &x, std::vector<double> &w) {
  std::vector<double> gx, gw;
  GaussLegendre(n, gx, gw);
  const double half = 0.5 * (hi - lo), mid = 0.5 * (hi + lo);
  for (int i = 0; i < n; i++) {
    x.push_back(mid + half * gx[i]);
    w.push_back(half * gw[i]);
  }
}

// Value and derivative of the wave table u at the stencil, divided by R: f = u/R, f' = u'/R - u/R^2.
inline void ReadWave(const SixPointStencil &st, const std::vector<cplx> &u, double R, cplx &f, cplx &df) {
  cplx v(0.0, 0.0), dv(0.0, 0.0);
  for (int i = 0; i < 6; i++) {
    v += st.value[i] * u[st.first + i];
    dv += st.derivative[i] * u[st.first + i];
  }
  f = v / R;
  df = dv / R - v / (R * R);
}

}  // namespace

SourceProjectionPlan PrepareSourceProjection(const SourceSetup &setup, const SourceProjectionRequest &request) {
  if (setup.l_values.empty() || !std::is_sorted(setup.l_values.begin(), setup.l_values.end()) ||
      std::adjacent_find(setup.l_values.begin(), setup.l_values.end()) != setup.l_values.end() ||
      setup.l_values.front() < 0)
    throw std::invalid_argument("source projection: the l values must be ascending, distinct and >= 0");
  if (!(setup.radius_fm > 0.0) || !(setup.alpha > 0.0) || !(setup.beta > 0.0) || !(setup.overlap.kappa_fm > 0.0) ||
      !(setup.entrance_step_fm > 0.0) || !(setup.exit_step_fm > 0.0) || !(setup.exit_k_max_fm >= 0.0) ||
      !setup.overlap.phi)
    throw std::invalid_argument("source projection: a, alpha, beta, kappa and the steps must be positive, phi given");
  SourceProjectionPlan p;
  p.setup = setup;
  p.request = request;
  const double pi = std::acos(-1.0);
  const double a = setup.radius_fm, alpha = setup.alpha, beta = setup.beta;
  const double rmin = setup.overlap.r_min_fm;
  p.u_max_fm = rmin + request.tail_per_kappa / setup.overlap.kappa_fm;
  p.entrance_k_fm = rho(setup.entrance.pair, 1.0, setup.entrance.energy_MeV);
  const double ka = p.entrance_k_fm, ks_max = setup.exit_k_max_fm;

  // Partial waves and the (u, cos theta_u) quadrature.
  const double rA = a + beta * p.u_max_fm, rS = alpha * a + p.u_max_fm;
  p.l_entrance_max = (int)std::ceil(ka * rA) + request.partial_wave_margin;
  p.l_exit_max = (int)std::ceil(ks_max * rS) + request.partial_wave_margin;
  if (p.l_entrance_max > request.l_cap || p.l_exit_max > request.l_cap)
    throw std::length_error("the partial-wave sums would need L up to " +
                            std::to_string(std::max(p.l_entrance_max, p.l_exit_max)) + " (kappa too small)");
  {
    const double width = std::min(request.panel_max_fm, request.panel_phase_rad / (ks_max + beta * ka));
    const int panels = std::max(2, (int)std::ceil((p.u_max_fm - rmin) / width));
    for (int k = 0; k < panels; k++) {
      const double lo = rmin + (p.u_max_fm - rmin) * k / panels, hi = rmin + (p.u_max_fm - rmin) * (k + 1) / panels;
      AppendGaussLegendre(request.nodes_per_panel, lo, hi, p.u_nodes, p.u_weights);
    }
  }
  const int cos_nodes =
      (int)std::ceil(request.angle_nodes_per_l * std::max(p.l_entrance_max, p.l_exit_max)) + request.angle_nodes_extra;
  AppendGaussLegendre(cos_nodes, -1.0, 1.0, p.cos_nodes, p.cos_weights);
  p.phi_at_u.resize(p.u_nodes.size());
  for (std::size_t i = 0; i < p.u_nodes.size(); i++) p.phi_at_u[i] = setup.overlap.phi(p.u_nodes[i]);

  // Triples (l; Ls, La): |Ls - La| <= l <= Ls + La, Ls + La + l even.
  for (int li = 0; li < (int)setup.l_values.size(); li++) {
    const int l = setup.l_values[li];
    for (int Ls = 0; Ls <= p.l_exit_max; Ls++)
      for (int La = std::abs(Ls - l); La <= std::min(p.l_entrance_max, Ls + l); La++) {
        if ((Ls + La + l) % 2) continue;
        SourceTriple t;
        t.l_index = li;
        t.l = l;
        t.l_exit = Ls;
        t.l_entrance = La;
        const double pre = std::sqrt((2.0 * Ls + 1.0) / (2.0 * l + 1.0));
        for (int nu = 0; nu <= std::min(l, La); nu++)
          t.kernel.push_back(pre * (nu % 2 ? -1.0 : 1.0) * (nu ? 2.0 : 1.0) *
                             clebsch_gordan_integer(Ls, 0, La, nu, l, nu));
        for (int m = 0; m <= std::min(l, Ls); m++)
          t.direction.push_back(clebsch_gordan_integer(Ls, m, La, 0, l, m) * std::sqrt((2.0 * La + 1.0) / (4.0 * pi)));
        p.triples.push_back(std::move(t));
      }
  }

  // The entrance waves, once.
  p.entrance_nodes = (int)std::ceil(rA / setup.entrance_step_fm) + kTableMargin;
  p.exit_nodes = (int)std::ceil(rS / setup.exit_step_fm) + kTableMargin;
  WaveTable(setup.entrance, setup.entrance_step_fm, p.entrance_nodes, p.l_entrance_max, request.matching,
            p.entrance_waves, p.entrance_phases);
  return p;
}

SourceExitWaves SolveSourceExit(const SourceProjectionPlan &plan, const WaveChannel &exit) {
  SourceExitWaves x;
  x.k_fm = rho(exit.pair, 1.0, exit.energy_MeV);
  WaveTable(exit, plan.setup.exit_step_fm, plan.exit_nodes, plan.l_exit_max, plan.request.matching, x.waves, x.phases);
  return x;
}

std::size_t SourceBlockSize(const SourceProjectionPlan &plan, int energies) {
  return plan.triples.size() * (std::size_t)energies * 2;
}

int SourceSumBlocks(const SourceProjectionPlan &plan, int energies) {
  const std::size_t points = plan.u_nodes.size() * plan.cos_nodes.size();
  const std::size_t bytes = (std::size_t)(plan.request.sum_megabytes * 1048576.0);
  const std::size_t per_block = std::max<std::size_t>(1, SourceBlockSize(plan, energies) * sizeof(cplx));
  return (int)std::max<std::size_t>(1, std::min<std::size_t>({points, (std::size_t)plan.request.sum_blocks,
                                                               bytes / per_block}));
}

void SumSourceBlock(const SourceProjectionPlan &plan, const std::vector<SourceExitWaves> &exits, int block, int blocks,
                    cplx *acc) {
  const double pi = std::acos(-1.0);
  const double a = plan.setup.radius_fm, alpha = plan.setup.alpha, beta = plan.setup.beta;
  const int nb = (int)exits.size();
  const int la_max = plan.l_entrance_max, ls_max = plan.l_exit_max;
  const int l_top = plan.setup.l_values.back();
  const int stride = l_top + 1;
  const int c_nodes = (int)plan.cos_nodes.size();
  const long long points = (long long)plan.u_nodes.size() * c_nodes;
  const int n_triples = (int)plan.triples.size();
  std::vector<double> yl, dyl, ya, dya;
  std::vector<cplx> fa(la_max + 1), dfa(la_max + 1);
  std::vector<cplx> fs((std::size_t)nb * (ls_max + 1)), dfs((std::size_t)nb * (ls_max + 1));
  const int end = (int)(points * (block + 1) / blocks);
  for (int ip = (int)(points * block / blocks); ip < end; ip++) {
    const int iu = ip / c_nodes, ic = ip % c_nodes;
    const double u = plan.u_nodes[iu], ct = plan.cos_nodes[ic], st = std::sqrt(std::max(0.0, 1.0 - ct * ct));
    const double wgt = 2.0 * pi * u * u * plan.phi_at_u[iu] * plan.u_weights[iu] * plan.cos_weights[ic];
    if (wgt == 0.0) continue;
    // R_s = alpha r z^ + u and R_a = r z^ + beta u at r = a, u in the xz plane.
    const double xs = u * st, zs = alpha * a + u * ct, xa = beta * u * st, za = a + beta * u * ct;
    const double Rs = std::hypot(xs, zs), Ra = std::hypot(xa, za);
    if (!(Rs > 0.0) || !(Ra > 0.0)) continue;
    const double ths = std::atan2(xs, zs), tha = std::atan2(xa, za);
    const double cs = zs / Rs, ss = xs / Rs, ca = za / Ra, sa = xa / Ra;
    SphericalHarmonicsWithDerivative(l_top, l_top, ths, yl, dyl);
    SphericalHarmonicsWithDerivative(la_max, l_top, tha - ths, ya, dya);
    const SixPointStencil sta = SixPointLagrange(Ra, plan.setup.entrance_step_fm, plan.entrance_nodes);
    for (int L = 0; L <= la_max; L++) ReadWave(sta, plan.entrance_waves[L], Ra, fa[L], dfa[L]);
    const SixPointStencil sts = SixPointLagrange(Rs, plan.setup.exit_step_fm, plan.exit_nodes);
    for (int b = 0; b < nb; b++)
      for (int L = 0; L <= ls_max; L++)
        ReadWave(sts, exits[b].waves[L], Rs, fs[b * (ls_max + 1) + L], dfs[b * (ls_max + 1) + L]);
    const double dthsdr = -alpha * ss / Rs, dthadr = -sa / Ra;
    for (int t = 0; t < n_triples; t++) {
      const SourceTriple &tr = plan.triples[t];
      double K = 0.0, Ks = 0.0, Ka = 0.0;  // K and dK/dtheta_s, dK/dtheta_a
      for (int nu = 0; nu < (int)tr.kernel.size(); nu++) {
        const double yln = yl[tr.l * stride + nu], dyln = dyl[tr.l * stride + nu];
        const double yan = ya[tr.l_entrance * stride + nu], dyan = dya[tr.l_entrance * stride + nu];
        K += tr.kernel[nu] * yln * yan;
        Ka += tr.kernel[nu] * yln * dyan;
        Ks += tr.kernel[nu] * (dyln * yan - yln * dyan);
      }
      const cplx fA = fa[tr.l_entrance], dfA = dfa[tr.l_entrance];
      const cplx A1 = wgt * fA * K;
      const cplx A2 = wgt * alpha * cs * fA * K;
      const cplx A3 = wgt * (fA * (Ks * dthsdr + Ka * dthadr) + ca * dfA * K);
      cplx *out = &acc[(std::size_t)t * nb * 2];
      for (int b = 0; b < nb; b++) {
        const cplx f = fs[b * (ls_max + 1) + tr.l_exit], df = dfs[b * (ls_max + 1) + tr.l_exit];
        out[2 * b] += A1 * f;
        out[2 * b + 1] += A2 * df + A3 * f;
      }
    }
  }
}

std::vector<SourceAmplitudes> CombineSourceBlocks(const SourceProjectionPlan &plan,
                                                  const std::vector<SourceExitWaves> &exits,
                                                  const std::vector<cplx> &block_sums, int blocks) {
  const double pi = std::acos(-1.0);
  const int nb = (int)exits.size(), n_triples = (int)plan.triples.size();
  const std::size_t size = SourceBlockSize(plan, nb);
  if (block_sums.size() != size * (std::size_t)blocks)
    throw std::invalid_argument("source projection: the block sums have the wrong size");
  std::vector<cplx> h(size, cplx(0.0, 0.0));
  for (int blk = 0; blk < blocks; blk++)
    for (std::size_t i = 0; i < size; i++) h[i] += block_sums[(std::size_t)blk * size + i];
  std::vector<SourceAmplitudes> out(nb);
  for (int b = 0; b < nb; b++) {
    out[b].value.resize(n_triples);
    out[b].derivative.resize(n_triples);
    for (int t = 0; t < n_triples; t++) {
      const SourceTriple &tr = plan.triples[t];
      // (4pi)^2/(k_f k_i) (-i)^Ls i^La e^{i(sigma^f + sigma^i)} sqrt(4pi/(2l+1)).
      cplx ph = exits[b].phases[tr.l_exit] * plan.entrance_phases[tr.l_entrance];
      const int n = ((tr.l_entrance - tr.l_exit) % 4 + 4) % 4;
      ph *= kIPower[n] * (16.0 * pi * pi / (exits[b].k_fm * plan.entrance_k_fm)) *
            std::sqrt(4.0 * pi / (2.0 * tr.l + 1.0));
      out[b].value[t] = ph * h[((std::size_t)t * nb + b) * 2];
      out[b].derivative[t] = ph * h[((std::size_t)t * nb + b) * 2 + 1];
    }
  }
  return out;
}

std::vector<SourceAmplitudes> ProjectSource(const SourceProjectionPlan &plan, const std::vector<WaveChannel> &exits) {
  std::vector<SourceAmplitudes> out;
  const int pass = std::max(1, plan.request.energies_per_pass);
  for (std::size_t e0 = 0; e0 < exits.size(); e0 += pass) {
    std::vector<SourceExitWaves> waves;
    for (std::size_t e = e0; e < std::min(exits.size(), e0 + pass); e++)
      waves.push_back(SolveSourceExit(plan, exits[e]));
    const int nb = (int)waves.size(), blocks = SourceSumBlocks(plan, nb);
    const std::size_t size = SourceBlockSize(plan, nb);
    std::vector<cplx> sums(size * blocks, cplx(0.0, 0.0));
    for (int blk = 0; blk < blocks; blk++) SumSourceBlock(plan, waves, blk, blocks, &sums[(std::size_t)blk * size]);
    for (SourceAmplitudes &x : CombineSourceBlocks(plan, waves, sums, blocks)) out.push_back(std::move(x));
  }
  return out;
}

SourceMultipoles SourceMultipolesAt(const SourceProjectionPlan &plan, const SourceAmplitudes &amp, double cos_fi) {
  const int nl = (int)plan.setup.l_values.size(), l_top = plan.setup.l_values.back();
  const double th = std::acos(std::max(-1.0, std::min(1.0, cos_fi)));
  std::vector<double> ys, dys;
  SphericalHarmonicsWithDerivative(plan.l_exit_max, l_top, th, ys, dys);
  SourceMultipoles s;
  s.value.resize(nl);
  s.derivative.resize(nl);
  for (int li = 0; li < nl; li++) {
    const int l = plan.setup.l_values[li];
    for (int m = 0; m <= l; m++) {
      cplx sm(0.0, 0.0), dm(0.0, 0.0);
      for (std::size_t t = 0; t < plan.triples.size(); t++) {
        const SourceTriple &tr = plan.triples[t];
        if (tr.l_index != li || m > tr.l_exit) continue;
        const double c = tr.direction[m] * ys[tr.l_exit * (l_top + 1) + m];
        if (c == 0.0) continue;
        sm += c * amp.value[t];
        dm += c * amp.derivative[t];
      }
      s.value[li].push_back(sm);
      s.derivative[li].push_back(dm);
    }
  }
  return s;
}

double SourceMomentumAmplitude(const SourceProjectionPlan &plan, double q) {
  const double pi = std::acos(-1.0);
  double phit = 0.0;
  for (std::size_t i = 0; i < plan.u_nodes.size(); i++) {
    const double z = q * plan.u_nodes[i];
    const double j0 = z < 1.0e-4 ? 1.0 - z * z / 6.0 : std::sin(z) / z;
    phit += plan.u_weights[i] * plan.u_nodes[i] * plan.u_nodes[i] * plan.phi_at_u[i] * j0;
  }
  phit *= 4.0 * pi;
  return phit;
}

}  // namespace transfer
}  // namespace physics
}  // namespace azure2
