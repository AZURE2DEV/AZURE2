// level_parameters (LevelParameters.h) on a synthetic System with analytic channel functions, using physics/ alone.
// The engine comparison on T1-T3 is tests/upgrade/level_parameters_from_cnuc_test.cpp.
//
//  1. Brune -> formal against the defining identity: the formal level matrix with constant B,
//       (A^-1)_{lm} = (E_l - E) delta - sum_c g_lc g_mc (S_c(E) - B_c + i P_c(E)),
//     and Brune's, (A^-1)_{lm} = (E~_l - E) delta - sum_c g~ g~ (S_c(E) + i P_c(E)) + sum_c g~ g~ s_c(l, m; E) with
//     s_c(l, l) = S_c(E~_l) and s_c(l, m) = [S_c(E~_l)(E - E~_m) - S_c(E~_m)(E - E~_l)]/(E~_l - E~_m), give the same
//     G(E) = g^T A g at every energy E (Brune, PRC 66, 044611 (2002)); the shift function is non-linear in E.  An
//     error in N or M, in the rank assignment or in B_c fails it at > 1e-3.  Also for levels listed out of energy
//     order and for Park input.
//  2. Round trips: observed -> Brune -> observed; Brune -> Park -> Brune; formal -> Brune -> formal (iteration
//     tolerance 1e-13 MeV).
//  3. Convention invariance: internal amplitudes -> each convention's columns -> back; the convention matrix is
//     orthogonal; observed widths written in jj columns give, through level_parameters, the same formal parameters
//     as the internal columns.
//  4. Gaps: a jj column dropped, and a block that lacks a channel spin, are reported and refused.
//  5. The formal basis (HANDOFF decision 10): B_c from the lowest level, so the formal parameters do not depend on the
//     listing order and the lowest level's formal energy is its Brune energy; each formal level oriented like its
//     Brune level (formal -> Brune returns the input signs); widths too large for the radii refused.
//  7. observed_levels against hand values: Gamma, C with an external part inside |.|^2 (one turning the sign), the
//     amplitude flag, N (Brune) and J (Park), a static M1 moment.
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include "azure2/physics/LevelParameters.h"

using namespace azure2::physics;
using cplx = std::complex<double>;

namespace {

int failures = 0;
long checks = 0;

void check(const std::string &what, bool ok, const std::string &detail = "") {
  checks++;
  if (!ok) {
    std::printf("  FAIL  %s %s\n", what.c_str(), detail.c_str());
    failures++;
  }
}

void report(const std::string &what, double value, double tol) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.3g (tolerance %.0g)", value, tol);
  std::printf("  %-62s %s\n", what.c_str(), buf);
  check(what, value <= tol, buf);
}

// ---- the synthetic system: 8Be-like, p + 7Li (1/2+, 3/2-), alpha + alpha, gamma to a 2+ state at 3.03 MeV.
System MakeSystem(const std::vector<int> &channels, const std::vector<double> &E, const std::vector<double> &gamma) {
  std::vector<Nucleus> nuclei = {{4, 8}};
  Partition p;
  p.Z = {{1, 3}};
  p.mass_u = {{1.00782503223, 7.016003434}};
  p.two_j = {{1, 3}};
  p.parity = {{1, -1}};
  p.separation_MeV = 17.255;
  Partition a;
  a.Z = {{2, 2}};
  a.mass_u = {{4.00260325413, 4.00260325413}};
  a.separation_MeV = -0.0918;
  Partition g;
  g.kind = PartitionKind::Photon;
  g.Z = {{0, 4}};
  g.mass_u = {{0.0, 8.00530510}};
  g.two_j = {{2, 4}};
  g.excitation_MeV = {{0.0, 3.03}};
  std::vector<Channel> all = {{0, ChannelKind::Particle, 1, 2, 4.5},   // c0: p, l = 1, s = 1
                              {0, ChannelKind::Particle, 1, 4, 4.5},   // c1: p, l = 1, s = 2
                              {1, ChannelKind::Particle, 2, 0, 5.0},   // c2: alpha, l = 2
                              {2, ChannelKind::Magnetic, 1, 4, 0.0},   // c3: M1 to the 2+ state
                              {0, ChannelKind::Particle, 3, 2, 4.5},   // c4: p, l = 3, s = 1
                              {0, ChannelKind::Particle, 3, 4, 4.5},   // c5: p, l = 3, s = 2
                              {1, ChannelKind::Particle, 0, 0, 5.0}};  // c6: alpha, l = 0
  std::vector<double> B(channels.size(), 0.0);
  JpiBlock blk(4, 1, 0, channels, B, E, gamma);
  return System(LevelBasis::Brune, nuclei, {p, a, g}, all, {blk});
}

// Analytic channel functions: S non-linear in E_c, P ~ sqrt(E_c) above threshold, constant ANC factor below.
ChannelFunctions Synthetic(const System &sys) {
  return [&sys](std::size_t b, std::size_t k, double E) {
    const Channel &ch = sys.channel(sys.block(b).channel(k));
    const Partition &p = sys.partition(ch.partition);
    const double Ec = E - (p.excitation_MeV[0] + p.excitation_MeV[1]) - p.separation_MeV;
    ChannelValues v;
    if (ch.kind == ChannelKind::Particle) {
      v.S = -0.8 + 0.06 * Ec + 0.004 * ch.l * Ec * Ec;
      v.dS_dE = 0.06 + 0.008 * ch.l * Ec;
      v.P = (Ec > 0.0) ? 0.05 * std::sqrt(Ec) * (1 + ch.l) : 0.3;
    } else {
      v.P = std::pow(std::fabs(Ec) / 197.3269804, 2 * ch.l + 1);
    }
    return v;
  };
}

double OpenP(const System &sys, const ChannelFunctions &cf, std::size_t k, double E) {
  const Channel &ch = sys.channel(sys.block(0).channel(k));
  if (ch.kind != ChannelKind::Particle) return 0.0;
  const Partition &p = sys.partition(ch.partition);
  return (E - p.separation_MeV > 0.0) ? cf(0, k, E).P : 0.0;
}

std::vector<cplx> Inverse(std::vector<cplx> A, std::size_t n) {
  std::vector<cplx> I(n * n, 0.0);
  for (std::size_t i = 0; i < n; i++) I[i * n + i] = 1.0;
  for (std::size_t c = 0; c < n; c++) {
    std::size_t piv = c;
    for (std::size_t r = c + 1; r < n; r++)
      if (std::abs(A[r * n + c]) > std::abs(A[piv * n + c])) piv = r;
    for (std::size_t j = 0; j < n; j++) {
      std::swap(A[c * n + j], A[piv * n + j]);
      std::swap(I[c * n + j], I[piv * n + j]);
    }
    const cplx d = A[c * n + c];
    for (std::size_t j = 0; j < n; j++) {
      A[c * n + j] /= d;
      I[c * n + j] /= d;
    }
    for (std::size_t r = 0; r < n; r++)
      if (r != c) {
        const cplx f = A[r * n + c];
        for (std::size_t j = 0; j < n; j++) {
          A[r * n + j] -= f * A[c * n + j];
          I[r * n + j] -= f * I[c * n + j];
        }
      }
  }
  return I;
}

// G(E) = gamma^T A(E) gamma, formal (constant B) or Brune form.
std::vector<cplx> GMatrix(const System &sys, const ChannelFunctions &cf, const std::vector<double> &E,
                          const std::vector<double> &g, const std::vector<double> *B, double energy) {
  const std::size_t n = E.size(), nc = sys.block(0).num_channels();
  std::vector<cplx> Ainv(n * n, 0.0);
  for (std::size_t l = 0; l < n; l++)
    for (std::size_t m = 0; m < n; m++) {
      cplx a = (l == m) ? cplx(E[l] - energy) : cplx(0.0);
      for (std::size_t k = 0; k < nc; k++) {
        if (sys.channel(sys.block(0).channel(k)).kind != ChannelKind::Particle) continue;
        const double gg = g[l * nc + k] * g[m * nc + k];
        const ChannelValues f = cf(0, k, energy);
        a -= gg * cplx(f.S - (B ? (*B)[k] : 0.0), OpenP(sys, cf, k, energy));
        if (!B) {
          const double Sl = cf(0, k, E[l]).S, Sm = cf(0, k, E[m]).S;
          a += gg * ((l == m) ? Sl : (Sl * (energy - E[m]) - Sm * (energy - E[l])) / (E[l] - E[m]));
        }
      }
      Ainv[l * n + m] = a;
    }
  const std::vector<cplx> A = Inverse(Ainv, n);
  std::vector<cplx> G(nc * nc, 0.0);
  for (std::size_t c = 0; c < nc; c++)
    for (std::size_t d = 0; d < nc; d++)
      for (std::size_t l = 0; l < n; l++)
        for (std::size_t m = 0; m < n; m++) G[c * nc + d] += g[l * nc + c] * A[l * n + m] * g[m * nc + d];
  return G;
}

double GDiff(const System &sys, const ChannelFunctions &cf, const std::vector<double> &Eb, const std::vector<double> &gb,
             const std::vector<double> &Ef, const std::vector<double> &gf, const std::vector<double> &B) {
  double worst = 0.0;
  for (double energy : {16.0, 17.6, 18.3, 19.1, 21.5}) {
    const std::vector<cplx> a = GMatrix(sys, cf, Ef, gf, &B, energy), b = GMatrix(sys, cf, Eb, gb, nullptr, energy);
    double scale = 0.0, d = 0.0;
    for (std::size_t i = 0; i < a.size(); i++) {
      scale = std::fmax(scale, std::abs(b[i]));
      d = std::fmax(d, std::abs(a[i] - b[i]));
    }
    worst = std::fmax(worst, d / scale);
  }
  return worst;
}

double MaxRel(const std::vector<double> &a, const std::vector<double> &b) {
  double scale = 0.0, d = 0.0;
  for (std::size_t i = 0; i < a.size(); i++) {
    scale = std::fmax(scale, std::fabs(b[i]));
    d = std::fmax(d, std::fabs(a[i] - b[i]));
  }
  return scale > 0.0 ? d / scale : d;
}

// Largest element-wise relative difference (for inputs spanning orders of magnitude).
double MaxElemRel(const std::vector<double> &a, const std::vector<double> &b) {
  double d = 0.0;
  for (std::size_t i = 0; i < a.size(); i++) d = std::fmax(d, std::fabs(a[i] - b[i]) / std::fmax(std::fabs(b[i]), 1e-300));
  return d;
}

}  // namespace

int main() {
  // Block 2+: c0 p (l=1, s=1), c1 p (l=1, s=2), c2 alpha (l=2), c3 M1.  Levels: 16.63 (bound for p), 18.15, 19.07,
  // listed out of energy order to exercise the rank assignment.
  const std::vector<int> chans = {0, 1, 2, 3};
  const std::vector<double> E = {18.15, 16.63, 19.07};
  // Observed: p width (eV) or ANC (fm^-1/2), alpha width (eV), Gamma_gamma (eV).
  const std::vector<double> observed = {-40000.0, 9000.0, 120000.0, 0.8,   // 18.15
                                        2.1, -0.9, 1200.0, 0.05,          // 16.63: p channels bound, ANCs
                                        150000.0, 60000.0, -400000.0, 2.0};  // 19.07
  const std::size_t nc = chans.size(), nl = E.size();
  const System sys = MakeSystem(chans, E, std::vector<double>(nl * nc, 0.0));
  const ChannelFunctions cf = Synthetic(sys);
  const std::vector<ProjectChannel> own = internal_channels(sys, 0);
  ProjectLevels obs;
  obs.E = E;
  obs.value = observed;

  std::printf("1. Brune -> formal against the level-matrix identity\n");
  const BasisAmplitudes brune = amplitudes_from_observed(sys, 0, own, obs, LevelBasis::Brune, cf);
  check("no warnings", brune.warnings.empty());
  const std::vector<double> B = boundary_values(sys, 0, E, cf);
  report("B_c = S_c(lowest level, 16.63, listed second)", std::fabs(B[0] - cf(0, 0, 16.63).S) + std::fabs(B[3] - B[0]), 0.0);
  const LevelSet formal = formal_from_brune(sys, 0, E, brune.gamma, B, cf);
  report("G(E) formal vs Brune, 5 energies (relative)", GDiff(sys, cf, E, brune.gamma, formal.E, formal.gamma, B), 1e-12);
  {
    // B_c = S_c(E~_low) makes E~_low an eigenvalue of the formal problem (Brune 2002): the level keeps its energy.
    report("lowest level: formal E = Brune E (relative)", std::fabs(formal.E[1] - E[1]) / E[1], 1e-14);
    bool oriented = true;
    for (std::size_t l = 0; l < nl; l++) {
      double dot = 0.0;
      for (std::size_t k = 0; k < nc; k++) dot += formal.gamma[l * nc + k] * brune.gamma[l * nc + k];
      oriented = oriented && dot > 0.0;
    }
    check("each formal level points the way of its Brune level (gamma . gamma~ > 0)", oriented);
    // The same block listed in energy order (16.63, 18.15, 19.07): the same formal levels, permuted.  With B_c at the
    // first listed level the formal energies differ.
    const std::vector<std::size_t> order = {1, 0, 2};
    std::vector<double> Es(nl), gs(nl * nc);
    for (std::size_t l = 0; l < nl; l++) {
      Es[l] = E[order[l]];
      for (std::size_t k = 0; k < nc; k++) gs[l * nc + k] = brune.gamma[order[l] * nc + k];
    }
    const std::vector<double> Bs = boundary_values(sys, 0, Es, cf);
    const LevelSet fs = formal_from_brune(sys, 0, Es, gs, Bs, cf);
    std::vector<double> back(nl * nc), backE(nl);
    for (std::size_t l = 0; l < nl; l++) {
      backE[order[l]] = fs.E[l];
      for (std::size_t k = 0; k < nc; k++) back[order[l] * nc + k] = fs.gamma[l * nc + k];
    }
    report("listing order: B_c unchanged", MaxRel(Bs, B), 0.0);
    report("listing order: formal energies (relative)", MaxRel(backE, formal.E), 1e-13);
    report("listing order: formal amplitudes, signs included (relative to max)", MaxRel(back, formal.gamma), 1e-13);
    std::vector<double> wide = brune.gamma;
    for (double &g : wide) g *= 30.0;
    bool threw = false;
    try {
      formal_from_brune(sys, 0, E, wide, B, cf);
    } catch (const std::domain_error &) {
      threw = true;
    }
    check("widths too large for the radii (M indefinite): refused (domain_error)", threw);
  }
  {
    // Strong mixing: two alpha channels (l = 2, 0), levels at 17.0, 18.0 and 18.05 MeV, the upper two with large
    // widths.  Here the eigenvector of the 18.05 MeV level has its largest component on another level, so taking the
    // eigensolver's sign (largest component positive) reverses that formal level; checked by applying that mutation.
    const std::vector<double> Em = {17.0, 18.0, 18.05}, gm = {0.3, 0.1, 0.8, 0.2, -0.3, 0.5};
    const System sm = MakeSystem({2, 6}, Em, std::vector<double>(6, 0.0));
    const ChannelFunctions cm = Synthetic(sm);
    const std::vector<double> Bm = boundary_values(sm, 0, Em, cm);
    const LevelSet fm = formal_from_brune(sm, 0, Em, gm, Bm, cm);
    bool oriented = true;
    for (std::size_t l = 0; l < 3; l++) oriented = oriented && fm.gamma[2 * l] * gm[2 * l] + fm.gamma[2 * l + 1] * gm[2 * l + 1] > 0.0;
    check("strong mixing: each formal level points the way of its Brune level", oriented);
    report("strong mixing: G(E) formal vs Brune (relative)", GDiff(sm, cm, Em, gm, fm.E, fm.gamma, Bm), 1e-12);
    FormalToBruneOptions tightm;
    tightm.energy_tolerance_MeV = 1e-13;
    const LevelSet bm = brune_from_formal(sm, 0, fm.E, fm.gamma, Bm, cm, tightm);
    report("strong mixing: Brune -> formal -> Brune, signs included", MaxRel(bm.gamma, gm) + MaxRel(bm.E, Em), 1e-11);
  }
  {
    const BasisAmplitudes park = amplitudes_from_observed(sys, 0, own, obs, LevelBasis::Park, cf);
    const std::vector<double> back = brune_from_park(sys, 0, E, park.gamma, cf);
    report("Park observed -> Brune = Brune observed", MaxRel(back, brune.gamma), 1e-15);
    std::vector<double> J;
    park_from_brune(sys, 0, E, brune.gamma, cf, &J);
    report("Park norm J from inputs = 1/(1 + sum g~^2 dS/dE)", MaxRel(park.park_norm, J), 1e-15);
  }

  std::printf("2. round trips\n");
  {
    const ProjectLevels again = observed_from_amplitudes(sys, 0, own, E, brune.gamma, {}, LevelBasis::Brune, cf);
    report("observed -> Brune -> observed (element-wise relative)", MaxElemRel(again.value, observed), 1e-14);
    const BasisAmplitudes park = amplitudes_from_observed(sys, 0, own, obs, LevelBasis::Park, cf);
    const ProjectLevels againP = observed_from_amplitudes(sys, 0, own, E, park.gamma, {}, LevelBasis::Park, cf);
    report("observed -> Park -> observed (element-wise relative)", MaxElemRel(againP.value, observed), 1e-14);
    const std::vector<double> p = park_from_brune(sys, 0, E, brune.gamma, cf);
    report("Brune -> Park -> Brune", MaxRel(brune_from_park(sys, 0, E, p, cf), brune.gamma), 1e-15);
    FormalToBruneOptions tight;
    tight.energy_tolerance_MeV = 1e-13;
    const LevelSet b2 = brune_from_formal(sys, 0, formal.E, formal.gamma, B, cf, tight);
    report("formal -> Brune: energies", MaxRel(b2.E, E), 1e-13);
    report("formal -> Brune: amplitudes, signs included", MaxRel(b2.gamma, brune.gamma), 1e-11);
    const LevelSet f2 = formal_from_brune(sys, 0, b2.E, b2.gamma, B, cf);
    report("formal -> Brune -> formal: energies", MaxRel(f2.E, formal.E), 1e-12);
    report("formal -> Brune -> formal: amplitudes, signs included", MaxRel(f2.gamma, formal.gamma), 1e-11);
    const LevelSet b3 = brune_from_formal(sys, 0, formal.E, formal.gamma, B, cf);  // the engine's tolerance 1e-6 MeV
    report("formal -> Brune at the engine's tolerance (1e-6 MeV): energies", MaxRel(b3.E, E), 1e-6);
  }

  std::printf("3. convention invariance\n");
  std::vector<std::pair<std::string, Convention> > conventions;
  {
    Convention c;
    c.spin_order = SpinOrder::HeavyFirst;
    conventions.push_back({"LS, s = j2 + j1", c});
    c = Convention();
    c.orbital_order = OrbitalOrder::OrbitalFirst;
    c.orbital_phase = OrbitalPhase::IPowL;
    c.magnetic_sign = -1;
    conventions.push_back({"LS, J = l + s, i^l, M sign -", c});
    for (SpinOrder so : {SpinOrder::LightFirst, SpinOrder::HeavyFirst})
      for (OrbitalOrder oo : {OrbitalOrder::SpinFirst, OrbitalOrder::OrbitalFirst}) {
        c = Convention();
        c.coupling = CouplingScheme::JJ;
        c.spin_order = so;
        c.orbital_order = oo;
        conventions.push_back({std::string("jj, ") + (so == SpinOrder::LightFirst ? "light" : "heavy") +
                                   (oo == OrbitalOrder::SpinFirst ? ", (j_x l)" : ", (l j_x)"),
                               c});
      }
  }
  for (const auto &nc_ : conventions) {
    const Convention &conv = nc_.second;
    const std::vector<ProjectChannel> cols = project_channels(sys, 0, conv);
    check(nc_.first + ": no gaps", recoupling_gaps(sys, 0, conv, cols).empty());
    const ConventionTransform T = convention_transform(sys, 0, conv, cols);
    double orth = 0.0;
    for (std::size_t a = 0; a < T.num_internal(); a++)
      for (std::size_t c = 0; c < T.num_internal(); c++) {
        double s = 0.0;
        for (std::size_t m = 0; m < T.num_project(); m++) s += T(a, m) * T(c, m);
        orth = std::fmax(orth, std::fabs(s - (a == c ? 1.0 : 0.0)));
      }
    report(nc_.first + ": T T^T = 1", orth, 1e-15);
    const std::vector<double> gp = T.to_project(brune.gamma, nl);
    report(nc_.first + ": internal -> project -> internal", MaxRel(T.to_internal(gp, nl), brune.gamma), 1e-15);
    // Widths written in the project's columns, then level_parameters in that convention.
    const ProjectLevels pobs = observed_from_amplitudes(sys, 0, cols, E, gp, {}, LevelBasis::Brune, cf);
    Convention withBasis = conv;
    withBasis.basis = LevelBasis::Brune;
    const InternalParameters ip = level_parameters(sys, 0, withBasis, cols, pobs, cf);
    report(nc_.first + ": Brune amplitudes from the project's widths", MaxRel(ip.basis_gamma, brune.gamma), 1e-14);
    report(nc_.first + ": formal energies", MaxRel(ip.E, formal.E), 1e-13);
    report(nc_.first + ": formal amplitudes, signs included", MaxRel(ip.gamma, formal.gamma), 1e-12);
    if (conv.coupling == CouplingScheme::JJ) {
      double mixed = 0.0;
      for (std::size_t m = 0; m < T.num_project(); m++)
        for (std::size_t a = 0; a < T.num_internal(); a++)
          if (std::fabs(T(a, m)) > 1e-12 && std::fabs(std::fabs(T(a, m)) - 1.0) > 1e-12) mixed = 1.0;
      check(nc_.first + ": the jj matrix mixes the two channel spins (not a relabeling)", mixed == 1.0);
    }
  }
  {
    Convention internal;
    const InternalParameters ip = level_parameters(sys, 0, internal, own, obs, cf);
    report("internal convention: level_parameters = the steps by hand (bit for bit)",
           MaxRel(ip.gamma, formal.gamma) + MaxRel(ip.E, formal.E) + MaxRel(ip.basis_gamma, brune.gamma), 0.0);
  }

  std::printf("4. gaps are reported\n");
  {
    Convention jj;
    jj.coupling = CouplingScheme::JJ;
    std::vector<ProjectChannel> cols = project_channels(sys, 0, jj);
    cols.erase(cols.begin());  // drop j = 1/2 of the p, l = 1 group
    const std::vector<std::string> gaps = recoupling_gaps(sys, 0, jj, cols);
    check("dropped jj column reported", !gaps.empty());
    for (const std::string &g : gaps) std::printf("    reported: %s\n", g.c_str());
    bool threw = false;
    try {
      convention_transform(sys, 0, jj, cols);
    } catch (const std::invalid_argument &) {
      threw = true;
    }
    check("dropped jj column refused", threw);
    // A block with only s = 2 of the p, l = 3 pair (J = 2 allows s = 1 and 2): jj needs the missing s = 1.
    const System lack = MakeSystem({5, 2}, {18.0}, {0.1, 0.2});
    const std::vector<ProjectChannel> full = project_channels(lack, 0, jj);
    const std::vector<std::string> g2 = recoupling_gaps(lack, 0, jj, full);
    check("block lacking a channel spin reported", !g2.empty());
    for (const std::string &g : g2) std::printf("    reported: %s\n", g.c_str());
    // LS relabelings never need a missing channel.
    Convention heavy;
    heavy.spin_order = SpinOrder::HeavyFirst;
    check("LS heavy-first on the incomplete block: no gap", recoupling_gaps(lack, 0, heavy, project_channels(lack, 0, heavy)).empty());
  }

  std::printf("6. diagnostics for callers\n");
  {
    // No amplitude-input column, so N = 1 and Park's J = D/2 for every level (Eq. 2.3).
    double worst = 0.0;
    bool layout = brune.half_denominator.size() == nl && brune.denominator_terms.size() == nl * nc;
    for (std::size_t l = 0; l < nl && layout; l++) {
      worst = std::fmax(worst, std::fabs(brune.half_denominator[l] - brune.park_norm[l]));
      double sum = 0.0;
      for (std::size_t k = 0; k < nc; k++) {
        const DenominatorTerm &t = brune.denominator_terms[l * nc + k];
        layout = layout && t.present == (k < 3) && brune.notes[l * nc + k] == ColumnNote::None;
        if (t.present) {
          sum += t.term;
          layout = layout && t.P > 0.0;
        }
      }
      worst = std::fmax(worst, std::fabs(brune.half_denominator[l] - (1.0 - sum / 2.0)));
    }
    check("denominator terms: particle width and ANC columns only, no notes", layout);
    report("D/2 = 1 - sum of the terms / 2 = Park's J (no amplitude input)", worst, 1e-15);
    ProjectLevels flagged = obs;
    flagged.is_amplitude.assign(observed.size(), 0);
    flagged.is_amplitude[3] = 1;  // the M1 column of level 18.15
    const BasisAmplitudes f = amplitudes_from_observed(sys, 0, own, flagged, LevelBasis::Brune, cf);
    check("amplitude flag on a photon column: noted and read as a width",
          f.notes[3] == ColumnNote::AmplitudeFlagIgnored && f.gamma == brune.gamma && f.warnings.size() == 1);

    FormalToBruneOptions one;
    one.energy_tolerance_MeV = 1e-13;
    one.max_iterations = 1;
    const LevelSet stuck = brune_from_formal(sys, 0, formal.E, formal.gamma, B, cf, one);
    // Level 1 is the lowest one, where B_c = S_c(E): its formal and Brune energies coincide, the first iteration
    // moves it by rounding only and it converges at once.  The other two stop after one iteration.
    bool kept = stuck.status[1] == LevelSet::Status::Converged && stuck.residual_MeV[1] == 0.0, shifts = true;
    for (std::size_t l : {std::size_t(0), std::size_t(2)})
      kept = kept && stuck.status[l] == LevelSet::Status::NotConverged && stuck.residual_MeV[l] > 1e-13 &&
          stuck.E[l] == formal.E[l];
    for (std::size_t l = 0; l < nl; l++)
      for (std::size_t k = 0; k < nc; k++)
        shifts = shifts && ((l != 1 && k < 3) == (stuck.boundary_shift[l * nc + k] != 0.0));
    check("not converged: status, step > tolerance, formal values kept (levels 0, 2; level 1 converged)", kept);
    check("boundary shifts of the last iteration: particle channels of the moving levels only", shifts);
    const LevelSet b3 = brune_from_formal(sys, 0, formal.E, formal.gamma, B, cf);
    bool converged = true;
    for (std::size_t l = 0; l < nl; l++) converged = converged && b3.status[l] == LevelSet::Status::Converged;
    check("converged levels report Converged", converged);
  }

  std::printf("7. observed quantities from total amplitudes (observed_levels, cnuc_split.md section 8)\n");
  {
    // Hand values: N = 1 + sum_particle g^2 S', W = s 2 P |g + g_ext|^2 / N, s = sign Re(g + g_ext); value by kind.
    const std::vector<double> g = {0.31, -0.12, 0.45, 0.004,    // 18.15
                                   0.20, -0.09, 0.15, 0.002,    // 16.63: p channels bound
                                   -0.28, 0.17, 0.52, -0.003};  // 19.07
    std::vector<ExternalWidth> ext(g.size());
    ext[3] = {true, cplx(-0.0015, 0.0021)};
    ext[11] = {true, cplx(0.0040, -0.0007)};  // turns the sign of the total: Re(-0.003 + 0.004) > 0
    std::vector<unsigned char> rwa(g.size(), 0);
    rwa[2] = 1;  // the alpha channel of 18.15 given as an amplitude
    double worst = 0.0;
    bool kinds = true;
    for (LevelBasis basis : {LevelBasis::Brune, LevelBasis::Park}) {
      const ObservedLevels o = observed_levels(sys, 0, own, E, g, ext, rwa, basis, cf);
      for (std::size_t l = 0; l < nl; l++) {
        double sum = 0.0;
        for (std::size_t k = 0; k < 3; k++) sum += g[l * nc + k] * g[l * nc + k] * cf(0, k, E[l]).dS_dE;
        const double N = (basis == LevelBasis::Park) ? 1.0 : 1.0 + sum, J = 1.0 - sum;
        worst = std::max(worst, std::fabs(o.norm[l] - ((basis == LevelBasis::Park) ? J : N)));
        for (std::size_t k = 0; k < nc; k++) {
          const std::size_t i = l * nc + k;
          const cplx tot = g[i] + (ext[i].present ? ext[i].value : cplx(0.0, 0.0));
          const double s = tot.real() < 0.0 ? -1.0 : 1.0;
          const double W = s * 2.0 * cf(0, k, E[l]).P * std::norm(tot) / N;
          const bool bound = (k < 2 && E[l] < 17.255);
          double v = bound ? s * std::sqrt(std::fabs(W)) : W * 1e6;
          if (rwa[i]) v = (basis == LevelBasis::Park) ? g[i] / std::sqrt(J) : g[i];
          kinds = kinds && o.kind[i] == (bound ? ObservedKind::ANC : ObservedKind::Width);
          worst = std::max(worst, std::fabs(o.W[i] - W) / std::fabs(W));
          worst = std::max(worst, std::fabs(o.value[i] - v) / std::fabs(v));
        }
      }
    }
    report("Gamma, C, gamma_ext inside |.|^2, N and J vs hand values (relative)", worst, 1e-15);
    check("kinds: ANC below the p threshold, Width above and for the photon", kinds);
    // A static M1 moment (photon channel to the level's own state): mu = s (2 P g^2 / N)^1/2.
    const ChannelFunctions moment = [&](std::size_t b, std::size_t k, double En) {
      ChannelValues v = cf(b, k, En);
      if (k == 3) {
        v.static_moment = true;
        v.P = 1.7;
      }
      return v;
    };
    const ObservedLevels m = observed_levels(sys, 0, own, E, g, {}, {}, LevelBasis::Brune, moment);
    double sum = 0.0;
    for (std::size_t k = 0; k < 3; k++) sum += g[k] * g[k] * cf(0, k, E[0]).dS_dE;
    const double mu = std::sqrt(2.0 * 1.7 * g[3] * g[3] / (1.0 + sum));
    check("static M1: kind MagneticMoment", m.kind[3] == ObservedKind::MagneticMoment);
    report("static M1: mu vs hand value (relative)", std::fabs(m.value[3] - mu) / mu, 1e-15);
  }

  std::printf("%s: %ld checks, %d failure(s)\n", failures ? "FAIL" : "PASS", checks, failures);
  return failures ? 1 : 0;
}
