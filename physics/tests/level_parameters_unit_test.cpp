// level_parameters (LevelParameters.h) on a synthetic System with analytic channel functions, using physics/ alone.
// The engine comparison on T1-T3 is tests/upgrade/level_parameters_from_cnuc_test.cpp.
//
//  1. Brune -> formal against the defining identity: the formal level matrix with constant B,
//       (A^-1)_{lm} = (E_l - E) delta - sum_c g_lc g_mc (S_c(E) - B_c + i P_c(E)),
//     and Brune's, (A^-1)_{lm} = (E~_l - E) delta - sum_c g~ g~ (S_c(E) + i P_c(E)) + sum_c g~ g~ s_c(l, m; E) with
//     s_c(l, l) = S_c(E~_l) and s_c(l, m) = [S_c(E~_l)(E - E~_m) - S_c(E~_m)(E - E~_l)]/(E~_l - E~_m), give the same
//     G(E) = g^T A g at every energy E (Brune, PRC 66, 044611 (2002)); the shift function is non-linear in E.  An
//     error in N or M, in the rank assignment or in B_c fails it at > 1e-3.  Also for levels listed out of energy
//     order, for B_c from the lowest level, and for Park input.
//  2. Round trips: observed -> Brune -> observed; Brune -> Park -> Brune; formal -> Brune -> formal (iteration
//     tolerance 1e-13 MeV).
//  3. Convention invariance: internal amplitudes -> each convention's columns -> back; the convention matrix is
//     orthogonal; observed widths written in jj columns give, through level_parameters, the same formal parameters
//     as the internal columns.
//  4. Gaps: a jj column dropped, and a block that lacks a channel spin, are reported and refused.
//  5. The proposals (B_c from the lowest level; largest |gamma| positive) do what they say and are off by default.
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
                              {0, ChannelKind::Particle, 3, 4, 4.5}};  // c5: p, l = 3, s = 2
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

// Largest difference up to a common sign of each level (row).
double MaxRelUpToLevelSign(const std::vector<double> &a, const std::vector<double> &b, std::size_t nc) {
  std::vector<double> c = a;
  for (std::size_t l = 0; l < a.size() / nc; l++) {
    double dot = 0.0;
    for (std::size_t k = 0; k < nc; k++) dot += a[l * nc + k] * b[l * nc + k];
    if (dot < 0.0)
      for (std::size_t k = 0; k < nc; k++) c[l * nc + k] = -c[l * nc + k];
  }
  return MaxRel(c, b);
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
  const std::vector<double> B = boundary_values(sys, 0, E, BoundaryRule::FirstListedLevel, cf);
  report("B_c = S_c(first listed level)", std::fabs(B[0] - cf(0, 0, E[0]).S) + std::fabs(B[3] - B[0]), 0.0);
  const LevelSet formal = formal_from_brune(sys, 0, E, brune.gamma, B, cf);
  report("G(E) formal vs Brune, 5 energies (relative)", GDiff(sys, cf, E, brune.gamma, formal.E, formal.gamma, B), 1e-12);
  const std::vector<double> Blow = boundary_values(sys, 0, E, BoundaryRule::LowestLevel, cf);
  report("B_c lowest level = S_c(16.63)", std::fabs(Blow[0] - cf(0, 0, 16.63).S), 0.0);
  const LevelSet formalLow = formal_from_brune(sys, 0, E, brune.gamma, Blow, cf, FormalSignRule::LargestPositive);
  report("G(E) with B_c from the lowest level", GDiff(sys, cf, E, brune.gamma, formalLow.E, formalLow.gamma, Blow), 1e-12);
  check("formal energies depend on B (the proposal changes parameters)", std::fabs(formalLow.E[0] - formal.E[0]) > 1e-6);
  {
    bool positive = true;
    for (std::size_t l = 0; l < nl; l++) {
      std::size_t big = 0;
      for (std::size_t k = 1; k < nc; k++)
        if (std::fabs(formalLow.gamma[l * nc + k]) > std::fabs(formalLow.gamma[l * nc + big])) big = k;
      positive = positive && formalLow.gamma[l * nc + big] > 0.0;
    }
    check("LargestPositive: largest |gamma| of each formal level is positive", positive);
    Convention def;
    check("proposals off by default", def.boundary == BoundaryRule::FirstListedLevel && def.formal_sign == FormalSignRule::AsComputed);
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
    report("formal -> Brune: amplitudes (up to the sign of a level)", MaxRelUpToLevelSign(b2.gamma, brune.gamma, nc), 1e-11);
    const LevelSet f2 = formal_from_brune(sys, 0, b2.E, b2.gamma, B, cf);
    report("formal -> Brune -> formal: energies", MaxRel(f2.E, formal.E), 1e-12);
    report("formal -> Brune -> formal: amplitudes (up to level sign)", MaxRelUpToLevelSign(f2.gamma, formal.gamma, nc), 1e-11);
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
    report(nc_.first + ": formal amplitudes (up to level sign)", MaxRelUpToLevelSign(ip.gamma, formal.gamma, nc), 1e-12);
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

  std::printf("%s: %ld checks, %d failure(s)\n", failures ? "FAIL" : "PASS", checks, failures);
  return failures ? 1 : 0;
}
