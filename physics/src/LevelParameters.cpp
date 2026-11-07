// Basis conversions of level_parameters (LevelParameters.h; docs/upgrade/level_parameters.md sections 2-4).  Each
// function reproduces the engine's arithmetic, in the engine's order of operations, with the file and line it follows
// named in a comment; the eigenproblems are SymmetricEigen.h's instead of GSL's.
#include "azure2/physics/LevelParameters.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include "azure2/physics/SymmetricEigen.h"

namespace azure2 {
namespace physics {

namespace {

constexpr double kEVPerMeV = 1e6;  // the <levels> widths are in eV, the amplitudes in MeV^1/2

bool IsPhoton(ChannelKind k) { return k == ChannelKind::Electric || k == ChannelKind::Magnetic; }
bool IsBeta(ChannelKind k) { return k == ChannelKind::Fermi || k == ChannelKind::GamowTeller; }

// A particle channel of a particle partition: the channels whose shift functions enter the level matrix.
bool IsParticle(const System &s, int c) {
  const Channel &ch = s.channel(c);
  return ch.kind == ChannelKind::Particle && s.partition(ch.partition).kind == PartitionKind::Particle;
}

// E_c = E_x - (e_0 + e_1) - S_alpha (MeV, c.m.), in the engine's order: CNuc.cpp TransformIn,
// theLevel->GetE() - pair->GetExE() - pair->GetSepE(), with GetExE() = e_light + e_heavy.
double ChannelEnergy(const System &s, int partition, double E) {
  const Partition &p = s.partition(partition);
  return E - (p.excitation_MeV[0] + p.excitation_MeV[1]) - p.separation_MeV;
}

std::string Label(const System &s, std::size_t b, std::size_t lambda) {
  const JpiBlock &blk = s.block(b);
  std::ostringstream o;
  o << "J^pi = " << blk.two_J() << "/2" << (blk.parity() > 0 ? "+" : "-") << " (block " << b << ") level " << lambda;
  if (lambda < blk.num_levels()) o << " (E = " << blk.E(lambda) << " MeV)";
  return o.str();
}

// For column m, the block channel the channel functions are read from: same partition, kind and l (and, for
// photon and beta channels, the same label).  P_c and S_c depend on (alpha, l) only.
std::vector<std::size_t> FunctionChannels(const System &s, std::size_t b, const std::vector<ProjectChannel> &columns) {
  const JpiBlock &blk = s.block(b);
  std::vector<std::size_t> fk(columns.size());
  for (std::size_t m = 0; m < columns.size(); m++) {
    const ProjectChannel &col = columns[m];
    bool found = false;
    for (std::size_t k = 0; k < blk.num_channels() && !found; k++) {
      const Channel &ch = s.channel(blk.channel(k));
      if (ch.partition != col.partition || ch.kind != col.kind || ch.l != col.l) continue;
      if (col.kind != ChannelKind::Particle && ch.two_s != col.two_x) continue;
      fk[m] = k;
      found = true;
    }
    if (!found) {
      std::ostringstream o;
      o << "level_parameters: block " << b << " has no channel of partition " << col.partition << " with l = " << col.l
        << " for project column " << m;
      throw std::invalid_argument(o.str());
    }
  }
  return fk;
}

void CheckLayout(const System &s, std::size_t b, std::size_t nColumns, const std::vector<double> &E,
                 const std::vector<double> &gamma, const char *who) {
  if (b >= s.num_blocks()) throw std::out_of_range(std::string(who) + ": block index out of range");
  if (gamma.size() != E.size() * nColumns)
    throw std::invalid_argument(std::string(who) + ": " + std::to_string(E.size()) + " levels x " +
                                std::to_string(nColumns) + " columns need as many amplitudes, got " +
                                std::to_string(gamma.size()));
}

// Rank of each energy in ascending order, ties in listed order (CNuc.cpp AscendingRanks).
std::vector<int> AscendingRanks(const std::vector<double> &energies) {
  std::vector<int> order(energies.size());
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](int a, int c) { return energies[a] < energies[c]; });
  std::vector<int> rank(energies.size());
  for (int k = 0; k < (int)order.size(); k++) rank[order[k]] = k;
  return rank;
}

}  // namespace

std::vector<ProjectChannel> internal_channels(const System &system, std::size_t b) {
  const JpiBlock &blk = system.block(b);
  std::vector<ProjectChannel> out;
  for (std::size_t k = 0; k < blk.num_channels(); k++) {
    const Channel &ch = system.channel(blk.channel(k));
    out.push_back(ProjectChannel{ch.partition, ch.kind, ch.l, ch.two_s});
  }
  return out;
}

// ------------------------------------------------------------------------------------------ observed <-> Brune

BasisAmplitudes amplitudes_from_observed(const System &system, std::size_t b, const std::vector<ProjectChannel> &columns,
                                         const ProjectLevels &levels, LevelBasis basis, const ChannelFunctions &cf) {
  if (basis == LevelBasis::Formal)
    throw std::invalid_argument("amplitudes_from_observed: observed widths convert to Brune or Park amplitudes; "
                                "formal ones follow from formal_from_brune");
  const std::size_t n = columns.size();
  CheckLayout(system, b, n, levels.E, levels.value, "amplitudes_from_observed");
  if (!levels.is_amplitude.empty() && levels.is_amplitude.size() != levels.value.size())
    throw std::invalid_argument("amplitudes_from_observed: is_amplitude has the wrong size");
  const std::vector<std::size_t> fk = FunctionChannels(system, b, columns);
  const JpiBlock &blk = system.block(b);
  const bool park = (basis == LevelBasis::Park);
  const bool halfIntegerJ = (blk.two_J() % 2 != 0);

  BasisAmplitudes out;
  out.gamma.assign(levels.value.size(), 0.0);
  out.park_norm.assign(levels.E.size(), 1.0);
  // CNuc.cpp TransformIn, first loop (lines 1200-1460 at e54e6fcc).
  for (std::size_t la = 0; la < levels.E.size(); la++) {
    const double E = levels.E[la];
    std::vector<double> tg(n, 0.0), pene(n, 1.0);
    std::vector<bool> negative(n, false), passThrough(n, false);
    double denom = 2.0;  // 2 (1 - sum over width channels of Gamma dS/dE / (2P))
    double numer = 1.0;  // 1 + sum over amplitude channels of gamma^2 dS/dE
    for (std::size_t m = 0; m < n; m++) {
      const ProjectChannel &col = columns[m];
      const double v = levels.value[la * n + m];
      const bool rwa = !levels.is_amplitude.empty() && levels.is_amplitude[la * n + m];
      negative[m] = (v < 0.0);
      if (col.kind == ChannelKind::Particle) {
        const double Ec = ChannelEnergy(system, col.partition, E);
        if (rwa) {
          const ChannelValues f = cf(b, fk[m], E);
          tg[m] = std::fabs(v);
          pene[m] = 1.0;
          passThrough[m] = true;
          numer += (v * v) * f.dS_dE;
        } else if (Ec > 0.0) {
          const ChannelValues f = cf(b, fk[m], E);
          tg[m] = std::fabs(v) / kEVPerMeV;
          double P = f.P;
          if (!(P > 0.0)) {
            // Just above a Coulomb threshold P = 0 in double precision: gamma = 0 (exact for an input of 0).
            if (tg[m] != 0.0)
              out.warnings.push_back(Label(system, b, la) + ", column " + std::to_string(m) +
                                     ": P = 0 at the level energy; the observed width cannot be converted, gamma = 0");
            tg[m] = 0.0;
            P = 1.0;
          } else {
            denom -= tg[m] / P * f.dS_dE;
          }
          pene[m] = P;
        } else {
          tg[m] = v * v;  // ANC^2 (fm^-1)
          if (Ec == 0.0) {
            if (tg[m] != 0.0)
              out.warnings.push_back(Label(system, b, la) + ", column " + std::to_string(m) +
                                     ": the level is at the channel threshold, where an ANC defines no width; gamma = 0");
            tg[m] = 0.0;
            pene[m] = 1.0;
            continue;
          }
          const ChannelValues f = cf(b, fk[m], E);
          pene[m] = f.P;
          denom -= tg[m] / pene[m] * f.dS_dE;
        }
      } else if (IsPhoton(col.kind)) {
        if (rwa)
          out.warnings.push_back(Label(system, b, la) + ", column " + std::to_string(m) +
                                 ": an amplitude input is supported for particle channels only; read as a width");
        const ChannelValues f = cf(b, fk[m], E);
        if (f.static_moment) {
          // A static moment: the input is the moment itself; its sign is reversed for half-integer J.
          int sign = (v < 0.0) ? -1 : 1;
          if (halfIntegerJ) sign = -sign;
          const double scaled = v * f.input_scale;
          tg[m] = scaled * scaled;
          pene[m] = f.P;
          negative[m] = (sign < 0);
        } else {
          tg[m] = std::fabs(v) / kEVPerMeV;
          double P = f.P;
          if (P < 1e-16) P = 1e-16;
          pene[m] = P;
        }
      } else {  // beta: the input is the feeding amplitude
        tg[m] = std::fabs(v);
        pene[m] = 1.0;
      }
    }
    // Park's J = 1 - sum gamma_Park^2 dS/dE = denom / (2 numer) (CNuc.cpp, comment before parkNorm).
    const double parkNorm = denom / (2.0 * numer);
    out.park_norm[la] = parkNorm;
    if (park && parkNorm <= 0.0) {
      std::ostringstream o;
      o << Label(system, b, la) << ": Park norm J = " << parkNorm << " is not positive";
      out.warnings.push_back(o.str());
    } else if (!park && denom < 0.0) {
      std::ostringstream o;
      o << Label(system, b, la) << ": 1 - sum_c Gamma_c (dS_c/dE)/(2 P_c) = " << denom / 2.0
        << " < 0: no reduced width reproduces these observed widths";
      out.warnings.push_back(o.str());
    }
    if (park) {
      denom = 2.0;
      numer = 1.0;
    }
    for (std::size_t m = 0; m < n; m++) {
      double g = tg[m];
      if (!IsBeta(columns[m].kind)) {
        if (!passThrough[m])
          g = std::sqrt(std::fabs(g / pene[m] * numer / denom));
        else if (park)
          g *= std::sqrt(std::fabs(parkNorm));
      } else if (park) {
        g *= std::sqrt(std::fabs(parkNorm));
      }
      if (negative[m]) g = -g;
      out.gamma[la * n + m] = g;
    }
  }
  return out;
}

ProjectLevels observed_from_amplitudes(const System &system, std::size_t b, const std::vector<ProjectChannel> &columns,
                                       const std::vector<double> &E, const std::vector<double> &gamma,
                                       const std::vector<unsigned char> &is_amplitude, LevelBasis basis,
                                       const ChannelFunctions &cf) {
  if (basis == LevelBasis::Formal) throw std::invalid_argument("observed_from_amplitudes: give Brune or Park amplitudes");
  const std::size_t n = columns.size();
  CheckLayout(system, b, n, E, gamma, "observed_from_amplitudes");
  if (!is_amplitude.empty() && is_amplitude.size() != gamma.size())
    throw std::invalid_argument("observed_from_amplitudes: is_amplitude has the wrong size");
  const std::vector<std::size_t> fk = FunctionChannels(system, b, columns);
  const bool park = (basis == LevelBasis::Park);
  const bool halfIntegerJ = (system.block(b).two_J() % 2 != 0);
  ProjectLevels out;
  out.form = ProjectLevels::Form::Observed;
  out.E = E;
  out.value.assign(gamma.size(), 0.0);
  out.is_amplitude = is_amplitude;
  for (std::size_t la = 0; la < E.size(); la++) {
    // Brune: N = 1 + sum_c gamma_c^2 dS_c/dE, Gamma_c = 2 P_c gamma_c^2 / N.  Park: J = 1 - sum gamma^2 dS/dE,
    // Gamma_c = 2 P_c gamma_c^2, and an amplitude-form or beta column is reported as Brune's, gamma / sqrt(J).
    double sum = 0.0;
    std::vector<ChannelValues> f(n);
    for (std::size_t m = 0; m < n; m++) {
      f[m] = cf(b, fk[m], E[la]);
      if (columns[m].kind == ChannelKind::Particle) sum += gamma[la * n + m] * gamma[la * n + m] * f[m].dS_dE;
    }
    const double N = park ? 1.0 : 1.0 + sum;
    const double J = 1.0 - sum;
    if (park && !(J > 0.0)) throw std::domain_error(Label(system, b, la) + ": Park norm J is not positive");
    for (std::size_t m = 0; m < n; m++) {
      const double g = gamma[la * n + m];
      const double sign = (g < 0.0) ? -1.0 : 1.0;
      const bool rwa = !is_amplitude.empty() && is_amplitude[la * n + m];
      double v = 0.0;
      if (columns[m].kind == ChannelKind::Particle) {
        const double Ec = ChannelEnergy(system, columns[m].partition, E[la]);
        if (rwa)
          v = park ? g / std::sqrt(J) : g;
        else if (Ec > 0.0)
          v = (f[m].P > 0.0) ? sign * 2.0 * f[m].P * g * g / N * kEVPerMeV : 0.0;
        else if (Ec < 0.0)
          v = sign * std::sqrt(2.0 * f[m].P * g * g / N);
      } else if (IsPhoton(columns[m].kind)) {
        if (f[m].static_moment) {
          const double s = halfIntegerJ ? -sign : sign;
          v = s * std::sqrt(2.0 * f[m].P * g * g / N) / f[m].input_scale;
        } else {
          const double P = (f[m].P < 1e-16) ? 1e-16 : f[m].P;
          v = sign * 2.0 * P * g * g / N * kEVPerMeV;
        }
      } else {
        v = park ? g / std::sqrt(J) : g;
      }
      out.value[la * n + m] = v;
    }
  }
  return out;
}

// --------------------------------------------------------------------------------------------- boundary values

std::vector<double> boundary_values(const System &system, std::size_t b, const std::vector<double> &E,
                                    BoundaryRule rule, const ChannelFunctions &cf) {
  const JpiBlock &blk = system.block(b);
  const std::size_t nc = blk.num_channels();
  std::vector<double> B(nc, 0.0);
  if (E.empty()) return B;  // a hard-sphere block: B drops out (CNuc.cpp CalcBoundaryConditions)
  std::size_t first = 0;
  if (rule == BoundaryRule::LowestLevel)
    for (std::size_t la = 1; la < E.size(); la++)
      if (E[la] < E[first]) first = la;
  // CNuc.cpp CalcBoundaryConditions (lines 1975-2003 at e54e6fcc): S_c at the chosen level for particle channels;
  // a non-particle channel copies channel 0's value at that point (it is never used).
  for (std::size_t k = 0; k < nc; k++) {
    if (system.partition(system.channel(blk.channel(k)).partition).kind == PartitionKind::Particle)
      B[k] = cf(b, k, E[first]).S;
    else
      B[k] = B[0];
  }
  return B;
}

// ---------------------------------------------------------------------------------------------- external capture

void subtract_external_capture(const System &system, std::size_t b, std::vector<double> &gamma,
                               const std::vector<ExternalWidth> &ext, bool ignore_zero_widths,
                               std::vector<std::string> &warnings) {
  const JpiBlock &blk = system.block(b);
  const std::size_t nc = blk.num_channels();
  if (gamma.size() % (nc ? nc : 1) != 0 || ext.size() != gamma.size())
    throw std::invalid_argument("subtract_external_capture: ext must have the layout of gamma");
  const std::size_t nl = nc ? gamma.size() / nc : 0;
  // CNuc.cpp TransformIn, second loop (lines 1488-1515 at e54e6fcc).
  for (std::size_t la = 0; la < nl; la++)
    for (std::size_t k = 0; k < nc; k++) {
      if (!IsPhoton(system.channel(blk.channel(k)).kind) || !ext[la * nc + k].present) continue;
      double &g = gamma[la * nc + k];
      if (ignore_zero_widths && std::fabs(g) < 1.0e-8) continue;
      const double re = ext[la * nc + k].value.real(), im = ext[la * nc + k].value.imag();
      if (g * g >= im * im) {
        g = (g < 0.0) ? -std::sqrt(g * g - im * im) - re : std::sqrt(g * g - im * im) - re;
      } else {
        warnings.push_back(Label(system, b, la) + ", channel " + std::to_string(k) +
                           ": imaginary part of the external amplitude exceeds the total; internal amplitude = -Re");
        g = -re;
      }
    }
}

// ------------------------------------------------------------------------------------------- Brune <-> formal

LevelSet formal_from_brune(const System &system, std::size_t b, const std::vector<double> &E,
                           const std::vector<double> &gamma, const std::vector<double> &B, const ChannelFunctions &cf,
                           FormalSignRule sign) {
  const JpiBlock &blk = system.block(b);
  const std::size_t nc = blk.num_channels();
  CheckLayout(system, b, nc, E, gamma, "formal_from_brune");
  if (B.size() != nc) throw std::invalid_argument("formal_from_brune: one B_c per block channel");
  const std::size_t n = E.size();
  LevelSet out;
  out.iterations.assign(n, 0);
  if (n == 0) return out;
  for (std::size_t a = 0; a < n; a++)
    for (std::size_t c = a + 1; c < n; c++)
      if (E[a] == E[c])
        throw std::invalid_argument("formal_from_brune: " + Label(system, b, a) + " and level " + std::to_string(c) +
                                    " have the same energy; the Brune transformation needs distinct energies");
  std::vector<bool> particle(nc);
  for (std::size_t k = 0; k < nc; k++) particle[k] = system.channel(blk.channel(k)).kind == ChannelKind::Particle;
  std::vector<double> shifts(n * nc, 0.0);  // S_c(E_lambda)
  for (std::size_t la = 0; la < n; la++)
    for (std::size_t k = 0; k < nc; k++)
      if (particle[k]) shifts[la * nc + k] = cf(b, k, E[la]).S;
  // CNuc.cpp TransformIn, second loop (lines 1532-1560 at e54e6fcc); the matrices of Brune (2002):
  //   M_{mu la} = delta - (1 - delta) sum_c (S_mu - S_la)/(E_mu - E_la) gamma_la gamma_mu
  //   N_{mu la} = delta (E_la + sum_c (S_la - B) gamma_la^2)
  //             + (1 - delta) sum_c ((E_mu S_la - E_la S_mu)/(E_mu - E_la) - B) gamma_la gamma_mu
  std::vector<double> Nm(n * n, 0.0), Mm(n * n, 0.0);
  for (std::size_t mu = 0; mu < n; mu++)
    for (std::size_t la = 0; la < n; la++) {
      if (la == mu) {
        Mm[mu * n + la] = 1.0;
        double sum = E[la];
        for (std::size_t k = 0; k < nc; k++)
          if (particle[k]) sum += (shifts[la * nc + k] - B[k]) * (gamma[la * nc + k] * gamma[la * nc + k]);
        Nm[mu * n + la] = sum;
      } else {
        double mSum = 0.0, nSum = 0.0;
        for (std::size_t k = 0; k < nc; k++) {
          if (!particle[k]) continue;
          const double gl = gamma[la * nc + k], gm = gamma[mu * nc + k];
          const double Sl = shifts[la * nc + k], Sm = shifts[mu * nc + k];
          mSum += (Sm - Sl) / (E[mu] - E[la]) * gl * gm;
          nSum += ((E[mu] * Sl - E[la] * Sm) / (E[mu] - E[la]) - B[k]) * gl * gm;
        }
        Mm[mu * n + la] = -mSum;
        Nm[mu * n + la] = nSum;
      }
    }
  const EigenResult eig = generalized_symmetric_eigen(Nm, Mm, n);
  if (!eig.ok) {
    std::ostringstream o;
    o << "No standard R-matrix (constant boundary conditions) corresponds to the levels of J^pi = " << blk.two_J()
      << "/2" << (blk.parity() > 0 ? "+" : "-") << " (block " << b
      << "): the matrix 1 - sum gamma gamma' dS/dE of the transformation is not positive definite";
    throw std::domain_error(o.str());
  }
  // Each level keeps its position and receives the eigenvalue of the rank of its Brune energy (CNuc.cpp 1565-1600).
  const std::vector<int> rank = AscendingRanks(E);
  out.E.assign(n, 0.0);
  out.gamma.assign(n * nc, 0.0);
  for (std::size_t la = 0; la < n; la++) {
    const std::size_t k = (std::size_t)rank[la];
    out.E[la] = eig.values[k];
    for (std::size_t c = 0; c < nc; c++) {
      double sum = 0.0;
      for (std::size_t mu = 0; mu < n; mu++) sum += eig.vectors[mu * n + k] * gamma[mu * nc + c];
      out.gamma[la * nc + c] = sum;
    }
  }
  if (sign == FormalSignRule::LargestPositive)
    for (std::size_t la = 0; la < n; la++) {
      std::size_t big = 0;
      for (std::size_t c = 1; c < nc; c++)
        if (std::fabs(out.gamma[la * nc + c]) > std::fabs(out.gamma[la * nc + big])) big = c;
      if (nc > 0 && out.gamma[la * nc + big] < 0.0)
        for (std::size_t c = 0; c < nc; c++) out.gamma[la * nc + c] = -out.gamma[la * nc + c];
    }
  return out;
}

LevelSet brune_from_formal(const System &system, std::size_t b, const std::vector<double> &E,
                           const std::vector<double> &gamma, const std::vector<double> &B, const ChannelFunctions &cf,
                           const FormalToBruneOptions &options) {
  const JpiBlock &blk = system.block(b);
  const std::size_t nc = blk.num_channels();
  CheckLayout(system, b, nc, E, gamma, "brune_from_formal");
  if (B.size() != nc) throw std::invalid_argument("brune_from_formal: one B_c per block channel");
  const std::size_t n = E.size();
  std::vector<bool> particle(nc);
  for (std::size_t k = 0; k < nc; k++) particle[k] = system.channel(blk.channel(k)).kind == ChannelKind::Particle;
  LevelSet out;
  out.E.assign(n, 0.0);
  out.gamma.assign(n * nc, 0.0);
  out.iterations.assign(n, 0);
  const std::vector<int> rank = AscendingRanks(E);
  // CNuc.cpp TransformOut (lines 2563-2745 at e54e6fcc), one level at a time.
  for (std::size_t lev = 0; lev < n; lev++) {
    std::vector<double> tempE(n);
    std::vector<double> tempG(n * nc);
    for (std::size_t i = 0; i < n; i++) {  // ascending formal energy
      tempE[rank[i]] = E[i];
      for (std::size_t k = 0; k < nc; k++) tempG[rank[i] * nc + k] = gamma[i * nc + k];
    }
    const std::size_t me = (std::size_t)rank[lev];
    std::vector<double> tempB = B;
    int iteration = 1;
    bool done = false;
    while (iteration <= options.max_iterations && !done) {
      std::vector<double> diff(nc, 0.0);
      for (std::size_t k = 0; k < nc; k++)
        if (particle[k]) {
          const double newB = cf(b, k, tempE[me]).S;
          diff[k] = newB - tempB[k];
          tempB[k] = newB;
        }
      std::vector<double> C(n * n);
      for (std::size_t mu = 0; mu < n; mu++)
        for (std::size_t mup = 0; mup < n; mup++) {
          double chanSum = 0.0;
          for (std::size_t k = 0; k < nc; k++)
            if (particle[k]) chanSum += diff[k] * tempG[mu * nc + k] * tempG[mup * nc + k];
          C[mu * n + mup] = (mu == mup) ? tempE[mu] - chanSum : -chanSum;
        }
      const EigenResult eig = symmetric_eigen(C, n);
      if (!eig.ok) {
        out.warnings.push_back(Label(system, b, lev) + ": could not transform the level (eigenvalue problem failed); "
                                                       "the formal values are kept");
        tempE[me] = E[lev];
        for (std::size_t k = 0; k < nc; k++) tempG[me * nc + k] = gamma[lev * nc + k];
        break;
      }
      if (std::fabs(eig.values[me] - tempE[me]) <= options.energy_tolerance_MeV) done = true;
      std::vector<double> newG(n * nc);
      for (std::size_t mu = 0; mu < n; mu++)
        for (std::size_t k = 0; k < nc; k++) {
          double s = 0.0;
          for (std::size_t mup = 0; mup < n; mup++) s += eig.vectors[mup * n + mu] * tempG[mup * nc + k];
          newG[mu * nc + k] = s;
        }
      tempE = eig.values;
      tempG.swap(newG);
      if (!done) {
        if (iteration == options.max_iterations) {
          std::ostringstream o;
          o << Label(system, b, lev) << ": not transformed after " << options.max_iterations
            << " iterations; the formal values are kept";
          out.warnings.push_back(o.str());
          tempE[me] = E[lev];
          for (std::size_t k = 0; k < nc; k++) tempG[me * nc + k] = gamma[lev * nc + k];
        }
        iteration++;
      }
    }
    out.E[lev] = tempE[me];
    for (std::size_t k = 0; k < nc; k++) out.gamma[lev * nc + k] = tempG[me * nc + k];
    out.iterations[lev] = iteration;
  }
  return out;
}

// --------------------------------------------------------------------------------------------- Brune <-> Park

namespace {
std::vector<double> ScaleLevels(const System &system, std::size_t b, const std::vector<double> &E,
                                const std::vector<double> &gamma, const ChannelFunctions &cf, bool toPark,
                                std::vector<double> *park_norm) {
  const JpiBlock &blk = system.block(b);
  const std::size_t nc = blk.num_channels();
  CheckLayout(system, b, nc, E, gamma, toPark ? "park_from_brune" : "brune_from_park");
  std::vector<double> out = gamma;
  if (park_norm) park_norm->assign(E.size(), 1.0);
  // CNuc.cpp ConvertAmplitudeBasis (lines 3422-3474 at e54e6fcc).
  for (std::size_t la = 0; la < E.size(); la++) {
    double sum = 0.0;  // sum_c gamma_c^2 dS_c/dE at E_lambda
    for (std::size_t k = 0; k < nc; k++) {
      if (!IsParticle(system, blk.channel(k))) continue;
      const double g = gamma[la * nc + k];
      if (g == 0.0) continue;
      sum += g * g * cf(b, k, E[la]).dS_dE;
    }
    double factor, J;
    if (toPark) {
      J = 1.0 / (1.0 + sum);
      factor = std::sqrt(1.0 / (1.0 + sum));
    } else {
      J = 1.0 - sum;
      if (J <= 0.0) {
        std::ostringstream o;
        o << Label(system, b, la) << ": Park norm J = " << J << " is not positive; no Brune parameters correspond";
        throw std::domain_error(o.str());
      }
      factor = 1.0 / std::sqrt(J);
    }
    if (park_norm) (*park_norm)[la] = J;
    for (std::size_t k = 0; k < nc; k++) out[la * nc + k] = gamma[la * nc + k] * factor;
  }
  return out;
}
}  // namespace

std::vector<double> park_from_brune(const System &system, std::size_t b, const std::vector<double> &E,
                                    const std::vector<double> &gamma, const ChannelFunctions &cf,
                                    std::vector<double> *park_norm) {
  return ScaleLevels(system, b, E, gamma, cf, true, park_norm);
}

std::vector<double> brune_from_park(const System &system, std::size_t b, const std::vector<double> &E,
                                    const std::vector<double> &gamma, const ChannelFunctions &cf,
                                    std::vector<double> *park_norm) {
  return ScaleLevels(system, b, E, gamma, cf, false, park_norm);
}

// ------------------------------------------------------------------------------------------------ composite

InternalParameters level_parameters(const System &system, std::size_t b, const Convention &convention,
                                    const std::vector<ProjectChannel> &columns, const ProjectLevels &levels,
                                    const ChannelFunctions &cf, const LevelParameterOptions &options) {
  const JpiBlock &blk = system.block(b);
  const std::size_t nc = blk.num_channels();
  const std::size_t nl = levels.E.size();
  InternalParameters out;

  // 1. Amplitudes of the project's basis, in the project's columns.
  std::vector<double> projectGamma;
  if (levels.form == ProjectLevels::Form::Observed) {
    out.basis = (convention.basis == LevelBasis::Park) ? LevelBasis::Park : LevelBasis::Brune;
    BasisAmplitudes a = amplitudes_from_observed(system, b, columns, levels, out.basis, cf);
    projectGamma.swap(a.gamma);
    out.park_norm.swap(a.park_norm);
    out.warnings.insert(out.warnings.end(), a.warnings.begin(), a.warnings.end());
  } else {
    out.basis = convention.basis;
    CheckLayout(system, b, columns.size(), levels.E, levels.value, "level_parameters");
    projectGamma = levels.value;
  }

  // 2. Recoupling to the block's channels.  For the internal coupling with the block's own columns T is the identity
  //    and is skipped, so the numbers pass unchanged (bit for bit).
  bool identity = is_internal_coupling(convention) && columns.size() == nc;
  const std::vector<ProjectChannel> own = internal_channels(system, b);
  for (std::size_t k = 0; identity && k < nc; k++)
    identity = own[k].partition == columns[k].partition && own[k].kind == columns[k].kind &&
        own[k].l == columns[k].l && own[k].two_x == columns[k].two_x;
  std::vector<double> g = identity ? projectGamma
                                   : convention_transform(system, b, convention, columns).to_internal(projectGamma, nl);

  // 3. External capture, 4. B_c.
  if (options.external) subtract_external_capture(system, b, g, *options.external, options.ignore_zero_widths, out.warnings);
  out.basis_E = levels.E;
  out.basis_gamma = g;
  out.B = boundary_values(system, b, levels.E, convention.boundary, cf);

  // 5. Formal parameters.
  if (!options.formal) return out;
  if (out.basis == LevelBasis::Formal) {
    out.E = levels.E;
    out.gamma = g;
    return out;
  }
  const std::vector<double> brune = (out.basis == LevelBasis::Park) ? brune_from_park(system, b, levels.E, g, cf) : g;
  LevelSet f = formal_from_brune(system, b, levels.E, brune, out.B, cf, convention.formal_sign);
  out.E.swap(f.E);
  out.gamma.swap(f.gamma);
  return out;
}

}  // namespace physics
}  // namespace azure2
