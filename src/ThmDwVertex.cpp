#include "ThmDwVertex.h"
#include "ThmExperiment.h"
#include "ThmLineshape.h"
#include "ThmNumerics.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "azure2/physics/transfer/SourceProjection.h"
#include "azure2/support/FloatingPointPrecision.h"

/*
 * Numerics (docs/source/theory/thm_implementation.rst, "Distorted-wave
 * entrance vertex"; docs/upgrade/transfer.md section 25).  The projection of
 * the source at r = a -- the reduced amplitudes per (l, Ls, La) by the 2D
 * quadrature in (u, cos theta_u), and S_lm(a), S_lm'(a) at a direction -- is
 * the library's (physics/transfer SourceProjection.h).  Here: the energy
 * grid and its passes (the exit waves of a pass and the blocks of the sum over
 * the points run on the threads; the library fixes the order of the sum), the
 * directions, the normalisation by 4 pi phi~(q), the Gram matrix and its
 * interpolation.
 */

namespace {

const double kGridStep = 0.02;  // MeV

}  // namespace

int ThmDwVertex::LIndex(int l) const {
  for (size_t i = 0; i < lvals.size(); i++)
    if (lvals[i] == l) return (int)i;
  return -1;
}

void ThmDwVertex::Factor(const double g[4], complex a[2], complex d[2]) {
  const double g11 = std::max(g[0], 0.0), g22 = std::max(g[1], 0.0);
  const complex g12(g[2], g[3]);
  a[0] = a[1] = d[0] = d[1] = complex(0.0, 0.0);
  if (g11 >= g22) {
    if (!(g11 > 0.0)) return;
    double r = std::sqrt(g11);
    a[0] = r;
    d[0] = g12 / r;
    d[1] = std::sqrt(std::max(0.0, g22 - std::norm(g12) / g11));
  } else {
    // Pivoted: G' = (g22, conj g12; g12, g11) for c' = (-1, B - 1).
    double r = std::sqrt(g22);
    d[0] = r;
    a[0] = std::conj(g12) / r;
    a[1] = std::sqrt(std::max(0.0, g11 - std::norm(g12) / g22));
  }
}

double ThmDwVertex::Vertex2(const double g[4], complex B) {
  const complex c1 = B - 1.0;
  const complex g12(g[2], g[3]);
  // c^+ G c with c = (c1, -1): |c1|^2 g11 + g22 - 2 Re(conj(c1) g12).
  return std::norm(c1) * g[0] + g[1] - 2.0 * std::real(std::conj(c1) * g12);
}

double ThmDwVertex::Vertex2Factored(const double g[4], complex B) {
  complex a[2], d[2];
  Factor(g, a, d);
  return std::norm(a[0] * (B - 1.0) - d[0]) + std::norm(a[1] * (B - 1.0) - d[1]);
}

std::string ThmDwVertex::Build(const ThmExperiment &x, const ThmDistortion::Kinematics &k, double a,
                               const std::vector<int> &ls, double eLo, double eHi,
                               const std::vector<double> &points) {
  const auto t0 = std::chrono::steady_clock::now();
  experiment = x.name;
  radius = a;
  lvals = ls;
  std::sort(lvals.begin(), lvals.end());
  lvals.erase(std::unique(lvals.begin(), lvals.end()), lvals.end());
  if (lvals.empty()) return "the entrance pair has no channel";
  if (!(a > 0.0)) return "the channel radius of the entrance pair is not positive";
  const int nl = (int)lvals.size();

  gridLo = eLo;
  gridStep = kGridStep;
  nE = std::max(4, (int)std::ceil((eHi - eLo) / kGridStep - 1.0e-9) + 1);
  if (!points.empty()) {
    dist.dataLo = *std::min_element(points.begin(), points.end());
    dist.dataHi = *std::max_element(points.begin(), points.end());
  }
  std::string why = dist.Setup(x, k, eLo);
  if (!why.empty()) return why;
  // Every data point must be reachable, as for R(E) (WorkingCopy): a lab angle the
  // spectator cannot reach would otherwise be clamped (SpectatorCos) to a
  // different direction.  Grid energies beyond the data may be out of reach.
  for (double e : points) {
    std::string w = dist.CheckEnergy(e);
    if (!w.empty()) return w;
  }
  for (int e = 0; e < nE; e++) {
    std::string w = dist.CheckEnergy(eLo + e * gridStep);
    if (!w.empty() && dist.angleKind != ThmDistortion::LAB)
      return "the vertex grid, " + ThmNumberText(eLo) + " to " + ThmNumberText(eLo + (nE - 1) * gridStep) +
             " MeV (the data and up to 0.3 MeV beyond, for the folding sub-points), cannot be built: " + w;
  }
  alpha = k.mA / (k.mx + k.mA);
  beta = k.ms / k.ma;
  const double ka = dist.aa.k;

  // The acceptance: the accepted spectator directions (spectatorAngles=
  // and/or the |p_s| cut of a ps window; ThmDistortion::AngleNodes).
  angles = dist.angWindow;
  nNodes = 1;
  nAng = 0;
  if (angles) {
    // The model is linear in G (sum over the Cholesky components of
    // |X a_k + Y d_k|^2 = |X|^2 G11 + |Y|^2 G22 + 2 Re(X* Y G12)), so the
    // average over the directions is one vertex with the averaged G: one node.
    dist.SetAngleSlots(eHi);
    nAng = dist.angSlots * dist.angNodes;
  }

  // The source's projection at r = a: the plan once (quadrature, triples, a + A waves), then per pass of
  // exit energies the s + F waves and the blocks of the sum over the points, each on the threads.
  namespace transfer = azure2::physics::transfer;
  std::vector<double> ksE(nE), esfE(nE);
  double ksMax = 0.0;
  for (int e = 0; e < nE; e++) {
    esfE[e] = dist.EsF(eLo + e * gridStep);
    ksE[e] = std::sqrt(2.0 * dist.sf.mu * esfE[e]) / hbarc;
    ksMax = std::max(ksMax, ksE[e]);
  }
  // The steps of the wave tables: 0.02 fm or 0.1 over the local wave number with the potential's depth.
  auto tableStep = [](const ThmDistortion::Channel &c, double e) {
    double kl = std::sqrt(2.0 * c.mu * (std::max(e, 0.0) + ThmDistortion::Depth(c))) / hbarc;
    return std::min(0.02, 0.1 / std::max(kl, 1.0e-6));
  };
  transfer::SourceSetup setup;
  setup.overlap.phi = [this](double r) { return dist.Phi(r); };
  setup.overlap.r_min_fm = dist.rmin;
  setup.overlap.kappa_fm = dist.kappa;
  setup.entrance = dist.LibraryChannel(dist.aa, dist.eAA);
  setup.entrance_step_fm = tableStep(dist.aa, dist.eAA);
  setup.exit_step_fm = tableStep(dist.sf, *std::max_element(esfE.begin(), esfE.end()));
  setup.exit_k_max_fm = ksMax;
  setup.alpha = alpha;
  setup.beta = beta;
  setup.radius_fm = a;
  setup.l_values = lvals;
  transfer::SourceProjectionRequest request;
  request.matching = dist.numerics.matching;
  transfer::SourceProjectionPlan plan;
  try {
    plan = transfer::PrepareSourceProjection(setup, request);
  } catch (const std::length_error &e) {
    return e.what();
  } catch (const std::exception &e) {
    return std::string("the a + A channel: ") + e.what();
  }
  laMax = plan.l_entrance_max;
  lsMax = plan.l_exit_max;
  uMax = plan.u_max_fm;
  uNodes = (int)plan.u_nodes.size();
  cNodes = (int)plan.cos_nodes.size();
  std::vector<transfer::SourceAmplitudes> amplitudes(nE);
  for (int e0 = 0; e0 < nE; e0 += request.energies_per_pass) {
    const int nb = std::min(request.energies_per_pass, nE - e0);
    std::vector<transfer::SourceExitWaves> exits(nb);
    std::vector<std::string> errors(nb);
#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < nb; b++) {
      const azure2::support::ExtendedPrecisionScope fpPrecision;
      const double e = eLo + (e0 + b) * gridStep;
      try {
        exits[b] = transfer::SolveSourceExit(plan, dist.LibraryChannel(dist.SfAt(e), dist.EsF(e)));
      } catch (const std::exception &x) {
        errors[b] = "the s + F channel at E = " + ThmNumberText(e) + " MeV: " + x.what();
      }
    }
    for (const std::string &err : errors)
      if (!err.empty()) return err;
    const int nBlocks = transfer::SourceSumBlocks(plan, nb);
    const size_t blockSize = transfer::SourceBlockSize(plan, nb);
    std::vector<complex> blockSums((size_t)nBlocks * blockSize, complex(0.0, 0.0));
#pragma omp parallel for schedule(dynamic, 1)
    for (int blk = 0; blk < nBlocks; blk++) {
      const azure2::support::ExtendedPrecisionScope fpPrecision;
      transfer::SumSourceBlock(plan, exits, blk, nBlocks, &blockSums[(size_t)blk * blockSize]);
    }
    std::vector<transfer::SourceAmplitudes> pass = transfer::CombineSourceBlocks(plan, exits, blockSums, nBlocks);
    for (int b = 0; b < nb; b++) amplitudes[e0 + b] = std::move(pass[b]);
  }

  // Per energy and node: s_m, d_m and the Gram matrix.
  w.assign((size_t)nE * nNodes, 0.0);
  q.assign((size_t)nE * nNodes, 0.0);
  G.assign((size_t)nE * nNodes * nl * 4, 0.0);
  Gd.assign((size_t)nE * nl * 4, 0.0);
  aw.assign((size_t)nE * nAng, 0.0);
  aq.assign((size_t)nE * nAng, 0.0);
  ath.assign((size_t)nE * nAng, 0.0);
  qd.assign(nE, 0.0);
  pd.assign(nE, 0.0);
  thd.assign(nE, 0.0);
  valid.assign(nE, 1);
  // The sign of phi~(q) in the spectatorAngle direction, per grid energy.
  std::vector<signed char> phiSign(nE, 0);
  for (int e = 0; e < nE; e++) {
    const double ks = ksE[e];
    // The Gram entries at an angle theta_sa between k_sF and k_aA, for q (fm^-1).
    auto gram = [&](double xsa, double qq, double *g) {
      const double phit = transfer::SourceMomentumAmplitude(plan, qq);
      const transfer::SourceMultipoles mp = transfer::SourceMultipolesAt(plan, amplitudes[e], xsa);
      for (int li = 0; li < nl; li++) {
        const int l = lvals[li];
        double g11 = 0.0, g22 = 0.0;
        complex g12(0.0, 0.0);
        for (int m = 0; m <= l; m++) {
          complex sm = mp.value[li][m], dm = mp.derivative[li][m];
          sm /= 4.0 * M_PI * phit;
          dm *= a / (4.0 * M_PI * phit);
          double mult = m ? 2.0 : 1.0;
          g11 += mult * std::norm(sm);
          g22 += mult * std::norm(dm);
          g12 += mult * std::conj(sm) * dm;
        }
        const double f = 4.0 * M_PI / (2.0 * l + 1.0);
        g[li * 4 + 0] = f * g11;
        g[li * 4 + 1] = f * g22;
        g[li * 4 + 2] = f * g12.real();
        g[li * 4 + 3] = f * g12.imag();
      }
      return phit;
    };
    // The spectatorAngle direction.
    {
      double th = 0.0;
      double xsa = dist.SpectatorCos(ks, &th);
      const double kb = beta * ka;
      qd[e] = std::sqrt(std::max(0.0, ks * ks + kb * kb - 2.0 * ks * kb * xsa));
      pd[e] = std::sqrt(std::max(0.0, ka * ka + alpha * alpha * ks * ks - 2.0 * ka * alpha * ks * xsa));
      thd[e] = std::acos(std::max(-1.0, std::min(1.0, xsa))) * 180.0 / M_PI;
      double phit = gram(xsa, qd[e], &Gd[(size_t)e * nl * 4]);
      if (!(std::fabs(phit) > 0.0) || !std::isfinite(phit))
        return "the plane-wave source phi~(q) vanishes at E = " + ThmNumberText(eLo + e * gridStep) + " MeV";
      phiSign[e] = phit > 0.0 ? 1 : -1;
    }
    if (angles) {
      // The accepted directions: event weight d cos(theta_cm) x acceptance x
      // |phi~(q)|^2 (the data are divided by |phi|^2 averaged over the same
      // events, and G is normalized by phi~ of its own node).
      std::vector<ThmDistortion::AngleNode> an;
      if (!dist.AngleNodes(ks, an)) {
        valid[e] = 0;
        continue;
      }
      double total = 0.0;
      std::vector<double> gk((size_t)nAng * nl * 4, 0.0);
      for (int i = 0; i < nAng; i++) {
        const ThmDistortion::AngleNode &n = an[i];
        double *ai = &aw[(size_t)e * nAng + i];
        aq[(size_t)e * nAng + i] = n.q * hbarc;
        ath[(size_t)e * nAng + i] = n.theta;
        if (!(n.w > 0.0)) continue;  // an empty branch
        double phit = gram(n.x, n.q, &gk[(size_t)i * nl * 4]);
        if (!(std::fabs(phit) > 0.0) || !std::isfinite(phit))
          return "the plane-wave source phi~(q) vanishes at E = " + ThmNumberText(eLo + e * gridStep) + " MeV, q = " +
                 ThmNumberText(n.q * hbarc) + " MeV/c";
        *ai = n.w * phit * phit;
        total += *ai;
      }
      if (!(total > 0.0) || !std::isfinite(total)) {
        valid[e] = 0;
        continue;
      }
      double qm = 0.0;
      for (int i = 0; i < nAng; i++) {
        const double wi = aw[(size_t)e * nAng + i] /= total;
        qm += wi * aq[(size_t)e * nAng + i];
        for (int c = 0; c < nl * 4; c++) G[(size_t)e * nl * 4 + c] += wi * gk[(size_t)i * nl * 4 + c];
      }
      w[e] = 1.0;
      q[e] = qm;
      continue;
    }
    w[e] = 1.0;
    q[e] = qd[e] * hbarc;
    std::copy(&Gd[(size_t)e * nl * 4], &Gd[(size_t)e * nl * 4] + nl * 4, &G[(size_t)e * nl * 4]);
  }

  // One direction: G is normalized by phi~(q)^2 of that direction, so a node
  // of phi~ between two grid energies is a pole of G that the cubic
  // interpolation runs across.  (With a window of directions the average is
  // weighted by phi~^2 and stays finite.)  Said where it reaches the data:
  // within two grid steps, the span of the interpolation.
  if (!angles && !points.empty()) {
    const double lo = dist.dataLo - 2.0 * gridStep, hi = dist.dataHi + 2.0 * gridStep;
    for (int e = 0; e + 1 < nE; e++) {
      const double e0 = gridLo + e * gridStep, e1 = e0 + gridStep;
      if (e1 < lo || e0 > hi || phiSign[e] == 0 || phiSign[e + 1] == 0 || phiSign[e] == phiSign[e + 1]) continue;
      dist.warnings.push_back("the plane-wave source phi~(q) of the vertex changes sign between E = " + ThmNumberText(e0) +
                              " and " + ThmNumberText(e1) + " MeV (q = " + ThmNumberText(qd[e] * hbarc) + " to " +
                              ThmNumberText(qd[e + 1] * hbarc) +
                              " MeV/c in the spectatorAngle direction): a node of the momentum distribution, "
                              "where the vertex, normalized by phi~(q), is singular and interpolated across the "
                              "pole; another spectatorAngle or boundState= moves it");
      break;
    }
  }

  // Every data point needs a reachable window; grid energies without one take
  // the nearest valid entry (they are beyond the data).
  for (double e : points) {
    double t = (e - gridLo) / gridStep;
    int i0 = std::max(0, std::min(nE - 1, (int)std::floor(t))), i1 = std::min(nE - 1, i0 + 1);
    if (!valid[i0] || !valid[i1]) {
      std::string why = dist.CheckWindow(e);
      return why.empty() ? "at E = " + ThmNumberText(e) + " MeV the spectator-direction window is out of reach next to it"
                         : why;
    }
  }
  for (int e = 0; e < nE; e++) {
    if (valid[e]) continue;
    int best = -1;
    for (int f = 0; f < nE; f++)
      if (valid[f] && (best < 0 || std::abs(f - e) < std::abs(best - e))) best = f;
    if (best < 0) return "the spectator-direction window is out of reach at every energy";
    for (int i = 0; i < nAng; i++) {
      aw[(size_t)e * nAng + i] = aw[(size_t)best * nAng + i];
      aq[(size_t)e * nAng + i] = aq[(size_t)best * nAng + i];
      ath[(size_t)e * nAng + i] = ath[(size_t)best * nAng + i];
    }
    for (int i = 0; i < nNodes; i++) {
      w[(size_t)e * nNodes + i] = w[(size_t)best * nNodes + i];
      q[(size_t)e * nNodes + i] = q[(size_t)best * nNodes + i];
      for (int c = 0; c < nl * 4; c++)
        G[((size_t)e * nNodes + i) * nl * 4 + c] = G[((size_t)best * nNodes + i) * nl * 4 + c];
    }
  }

  std::ostringstream d;
  d.precision(6);
  d << "vertexModel=dw: surface term of the prior-form DWBA; a + A: " << dist.ChannelText(0)
    << "; s + F: " << dist.ChannelText(1) << "; bound state " << (dist.yukawa ? "yukawa" : "whittaker")
    << (dist.rmin > 0.0 ? ", r >= " + ThmNumberText(dist.rmin) + " fm" : "") << "; "
    << (angles ? "spectator directions " + dist.AngleText() + " (" + ThmNumberText(dist.angNodes) +
                     " nodes in cos theta_cm" + (dist.angSlots > 1 ? " per branch" : "") +
                     (dist.qCut && !dist.angAll ? ", |p_s| cut by the ps window" : "") +
                     "; weight d cos theta_cm x acceptance x |phi~(q)|^2)"
               : std::string("spectator angle ") +
                     (dist.angleKind == ThmDistortion::QF
                          ? "qf (k_sF along k_aA)"
                          : (dist.angleKind == ThmDistortion::LAB ? "lab " : "cm ") + ThmNumberText(dist.angle) + " deg"));
  description = d.str();
  buildSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return "";
}

bool ThmDwVertex::Reached(double energy) const {
  if (nE == 0) return false;
  const double t = (energy - gridLo) / gridStep;
  if (t < 0.0 || t > nE - 1.0) return false;
  const int i0 = std::max(0, std::min(nE - 1, (int)std::floor(t))), i1 = std::min(nE - 1, i0 + 1);
  return valid.empty() || (valid[i0] && valid[i1]);
}

void ThmDwVertex::AngleNodesAt(double energy, std::vector<double> &weight, std::vector<double> &qk,
                               std::vector<double> &theta) const {
  weight.clear();
  qk.clear();
  theta.clear();
  if (!angles || nE == 0) return;
  const int e = std::max(0, std::min(nE - 1, (int)std::lround((energy - gridLo) / gridStep)));
  for (int i = 0; i < nAng; i++) {
    weight.push_back(aw[(size_t)e * nAng + i]);
    qk.push_back(aq[(size_t)e * nAng + i]);
    theta.push_back(ath[(size_t)e * nAng + i]);
  }
}

void ThmDwVertex::Interpolate(double energy, std::vector<double> &weight, std::vector<double> &qk,
                              std::vector<double> &g, std::vector<double> &gd, double &qDelta, double &pDelta,
                              bool *outside) const {
  const int nl = (int)lvals.size();
  double t = (energy - gridLo) / gridStep;
  bool out = t < 0.0 || t > nE - 1.0;
  if (outside) *outside = out;
  t = std::max(0.0, std::min(nE - 1.0, t));
  int j = std::min(std::max((int)std::floor(t), 1), nE - 3);
  double c[4];
  ThmCubicLagrange(t - j, c);
  auto lag = [&](const std::vector<double> &v, size_t stride, size_t off) {
    double r = 0.0;
    for (int i = 0; i < 4; i++) r += c[i] * v[(size_t)(j - 1 + i) * stride + off];
    return r;
  };
  weight.assign(nNodes, 0.0);
  qk.assign(nNodes, 0.0);
  g.assign((size_t)nNodes * nl * 4, 0.0);
  gd.assign((size_t)nl * 4, 0.0);
  double total = 0.0;
  for (int k = 0; k < nNodes; k++) {
    weight[k] = std::max(0.0, lag(w, nNodes, k));
    total += weight[k];
    qk[k] = lag(q, nNodes, k);
    for (int m = 0; m < nl * 4; m++) g[(size_t)k * nl * 4 + m] = lag(G, (size_t)nNodes * nl * 4, (size_t)k * nl * 4 + m);
  }
  if (total > 0.0)
    for (double &v : weight) v /= total;
  for (int m = 0; m < nl * 4; m++) gd[m] = lag(Gd, (size_t)nl * 4, m);
  qDelta = lag(qd, 1, 0);
  pDelta = lag(pd, 1, 0);
}

void ThmDwVertex::Evaluate(double energy, At &out) const {
  std::vector<double> qk, g, gd;
  double qq, pp;
  Interpolate(energy, out.weight, qk, g, gd, qq, pp, &out.outside);
  const int nl = (int)lvals.size();
  out.nodes = nNodes;
  out.nl = nl;
  out.a.assign((size_t)nNodes * nl * 2, complex(0.0, 0.0));
  out.d.assign((size_t)nNodes * nl * 2, complex(0.0, 0.0));
  for (int k = 0; k < nNodes; k++)
    for (int li = 0; li < nl; li++) Factor(&g[((size_t)k * nl + li) * 4], &out.a[((size_t)k * nl + li) * 2], &out.d[((size_t)k * nl + li) * 2]);
}
