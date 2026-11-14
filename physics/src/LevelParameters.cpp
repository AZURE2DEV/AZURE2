// Basis conversions of level_parameters (LevelParameters.h; docs/upgrade/level_parameters.md sections 2-4).  Each
// function reproduces the arithmetic of the engine code it replaced, with the file and line it follows (at e54e6fcc)
// named in a comment; the eigenproblems are SymmetricEigen.h's.  The formal basis takes B_c at the lowest level and
// orients each formal level like its Brune level (HANDOFF decision 10; cnuc_split.md section 7).
#include "azure2/physics/LevelParameters.h"

#include <algorithm>
#include <cmath>
#include <limits>
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

// E_c = (E_x - S_alpha) - (e_0 + e_1) (MeV, c.m.): the energy above the partition's ground-state threshold less the
// excitation of the pair, in the order the engine forms it at a data point (Unit.cpp, inEnergy - GetSepE() -
// GetExE(), GetExE() = e_light + e_heavy) and in the adapter's channel functions (cnuc_split.md section 7).
double ChannelEnergy(const System &s, int partition, double E) {
  const Partition &p = s.partition(partition);
  return E - p.separation_MeV - (p.excitation_MeV[0] + p.excitation_MeV[1]);
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
  out.half_denominator.assign(levels.E.size(), 1.0);
  out.denominator_terms.assign(levels.value.size(), DenominatorTerm());
  out.notes.assign(levels.value.size(), ColumnNote::None);
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
            if (tg[m] != 0.0) {
              out.notes[la * n + m] = ColumnNote::ZeroPenetrability;
              out.warnings.push_back(Label(system, b, la) + ", column " + std::to_string(m) +
                                     ": P = 0 at the level energy; the observed width cannot be converted, gamma = 0");
            }
            tg[m] = 0.0;
            P = 1.0;
          } else {
            denom -= tg[m] / P * f.dS_dE;
            out.denominator_terms[la * n + m] = DenominatorTerm{true, tg[m] / P * f.dS_dE, P, f.dS_dE};
          }
          pene[m] = P;
        } else {
          tg[m] = v * v;  // ANC^2 (fm^-1)
          if (Ec == 0.0) {
            if (tg[m] != 0.0) {
              out.notes[la * n + m] = ColumnNote::AtThreshold;
              out.warnings.push_back(Label(system, b, la) + ", column " + std::to_string(m) +
                                     ": the level is at the channel threshold, where an ANC defines no width; gamma = 0");
            }
            tg[m] = 0.0;
            pene[m] = 1.0;
            continue;
          }
          const ChannelValues f = cf(b, fk[m], E);
          pene[m] = f.P;
          denom -= tg[m] / pene[m] * f.dS_dE;
          out.denominator_terms[la * n + m] = DenominatorTerm{true, tg[m] / pene[m] * f.dS_dE, pene[m], f.dS_dE};
        }
      } else if (IsPhoton(col.kind)) {
        if (rwa) {
          out.notes[la * n + m] = ColumnNote::AmplitudeFlagIgnored;
          out.warnings.push_back(Label(system, b, la) + ", column " + std::to_string(m) +
                                 ": an amplitude input is supported for particle channels only; read as a width");
        }
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
        if (rwa) {
          out.notes[la * n + m] = ColumnNote::AmplitudeFlagIgnored;
          out.warnings.push_back(Label(system, b, la) + ", column " + std::to_string(m) +
                                 ": an amplitude input is supported for particle channels only; ignored");
        }
        tg[m] = std::fabs(v);
        pene[m] = 1.0;
      }
    }
    // Park's J = 1 - sum gamma_Park^2 dS/dE = denom / (2 numer) (CNuc.cpp, comment before parkNorm).
    const double parkNorm = denom / (2.0 * numer);
    out.park_norm[la] = parkNorm;
    out.half_denominator[la] = denom / 2.0;
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

ObservedLevels observed_levels(const System &system, std::size_t b, const std::vector<ProjectChannel> &columns,
                               const std::vector<double> &E, const std::vector<double> &gamma,
                               const std::vector<ExternalWidth> &ext, const std::vector<unsigned char> &is_amplitude,
                               LevelBasis basis, const ChannelFunctions &cf) {
  if (basis == LevelBasis::Formal) throw std::invalid_argument("observed_levels: give Brune or Park amplitudes");
  const std::size_t n = columns.size();
  CheckLayout(system, b, n, E, gamma, "observed_levels");
  if (!is_amplitude.empty() && is_amplitude.size() != gamma.size())
    throw std::invalid_argument("observed_levels: is_amplitude has the wrong size");
  if (!ext.empty() && ext.size() != gamma.size()) throw std::invalid_argument("observed_levels: ext has the wrong size");
  const std::vector<std::size_t> fk = FunctionChannels(system, b, columns);
  const bool park = (basis == LevelBasis::Park);
  const bool halfIntegerJ = (system.block(b).two_J() % 2 != 0);
  ObservedLevels out;
  out.kind.assign(gamma.size(), ObservedKind::Width);
  out.W.assign(gamma.size(), 0.0);
  out.value.assign(gamma.size(), 0.0);
  out.norm.assign(E.size(), 1.0);
  // cnuc_split.md section 8, Eq. (8.1); the arithmetic of CNuc::TransformOut's second loop (9d50667f, 2423-2516).
  for (std::size_t la = 0; la < E.size(); la++) {
    double sum = 0.0;  // sum over particle columns of gamma_int^2 dS/dE
    double J = 1.0;    // Park's J, accumulated as CompoundSystem::CalcShiftFunctions does
    std::vector<ChannelValues> f(n);
    for (std::size_t m = 0; m < n; m++) {
      f[m] = cf(b, fk[m], E[la]);
      if (columns[m].kind == ChannelKind::Particle) {
        const double g = gamma[la * n + m];
        sum += f[m].dS_dE * (g * g);
        J -= (g * g) * f[m].dS_dE;
      }
    }
    const double N = park ? 1.0 : 1.0 + sum;
    out.norm[la] = park ? J : N;
    for (std::size_t m = 0; m < n; m++) {
      const std::size_t i = la * n + m;
      const ProjectChannel &col = columns[m];
      const double g = gamma[i];
      const std::complex<double> tot = (!ext.empty() && ext[i].present) ? g + ext[i].value : std::complex<double>(g, 0.0);
      const double s = (tot.real() < 0.0) ? -1.0 : 1.0;
      const double mod2 = tot.real() * tot.real() + tot.imag() * tot.imag();
      const bool rwa = !is_amplitude.empty() && is_amplitude[i];
      if (IsBeta(col.kind)) {
        out.kind[i] = ObservedKind::Beta;
        out.W[i] = park ? tot.real() / std::sqrt(std::fabs(J)) : tot.real();
        out.value[i] = out.W[i];
        continue;
      }
      double P = f[m].P, sign = s;
      if (col.kind == ChannelKind::Particle) {
        const double Ec = ChannelEnergy(system, col.partition, E[la]);
        out.kind[i] = (Ec < 0.0) ? ObservedKind::ANC : ObservedKind::Width;
      } else if (f[m].static_moment) {
        if (halfIntegerJ) sign = -sign;  // the engine's static-moment factor changes sign for half-integer J
        if (col.kind == ChannelKind::Magnetic && col.l == 1) out.kind[i] = ObservedKind::MagneticMoment;
        else if (col.kind == ChannelKind::Electric && col.l == 2) out.kind[i] = ObservedKind::QuadrupoleMoment;
      }
      const double W = sign * 2.0 * mod2 * P / N;
      out.W[i] = W;
      const double r = (W < 0.0) ? -std::sqrt(std::fabs(W)) : std::sqrt(std::fabs(W));
      switch (out.kind[i]) {
        case ObservedKind::Width: out.value[i] = W * kEVPerMeV; break;
        case ObservedKind::ANC:
        case ObservedKind::MagneticMoment: out.value[i] = r; break;
        case ObservedKind::QuadrupoleMoment: out.value[i] = r / f[m].input_scale; break;
        case ObservedKind::Beta: break;
      }
      // A column given as an amplitude is reported as one, in Brune's normalization.
      if (rwa && col.kind == ChannelKind::Particle) out.value[i] = (park && J > 0.0) ? g / std::sqrt(J) : g;
    }
  }
  return out;
}

ProjectLevels observed_from_amplitudes(const System &system, std::size_t b, const std::vector<ProjectChannel> &columns,
                                       const std::vector<double> &E, const std::vector<double> &gamma,
                                       const std::vector<unsigned char> &is_amplitude, LevelBasis basis,
                                       const ChannelFunctions &cf) {
  ProjectLevels out;
  out.form = ProjectLevels::Form::Observed;
  out.E = E;
  out.value = observed_levels(system, b, columns, E, gamma, {}, is_amplitude, basis, cf).value;
  out.is_amplitude = is_amplitude;
  return out;
}

// --------------------------------------------------------------------------------------------- boundary values

std::vector<double> boundary_values(const System &system, std::size_t b, const std::vector<double> &E,
                                    const ChannelFunctions &cf) {
  const JpiBlock &blk = system.block(b);
  const std::size_t nc = blk.num_channels();
  std::vector<double> B(nc, 0.0);
  if (E.empty()) return B;  // a hard-sphere block: B drops out (CompoundSystemBuild.cpp CalcBoundaryConditions)
  std::size_t lowest = 0;
  for (std::size_t la = 1; la < E.size(); la++)
    if (E[la] < E[lowest]) lowest = la;
  // CNuc.cpp CalcBoundaryConditions (lines 1975-2003 at e54e6fcc), at the lowest level instead of the first listed:
  // S_c there for particle channels; a non-particle channel copies channel 0's value (it is never used).
  for (std::size_t k = 0; k < nc; k++) {
    if (system.partition(system.channel(blk.channel(k)).partition).kind == PartitionKind::Particle)
      B[k] = cf(b, k, E[lowest]).S;
    else
      B[k] = B[0];
  }
  return B;
}

// ---------------------------------------------------------------------------------------------- external capture

void subtract_external_capture(const System &system, std::size_t b, std::vector<double> &gamma,
                               const std::vector<ExternalWidth> &ext, std::vector<std::string> &warnings, std::vector<std::size_t> *exceeded) {
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
      const double re = ext[la * nc + k].value.real(), im = ext[la * nc + k].value.imag();
      if (g * g >= im * im) {
        g = (g < 0.0) ? -std::sqrt(g * g - im * im) - re : std::sqrt(g * g - im * im) - re;
      } else {
        warnings.push_back(Label(system, b, la) + ", channel " + std::to_string(k) +
                           ": imaginary part of the external amplitude exceeds the total; internal amplitude = -Re");
        if (exceeded) exceeded->push_back(la * nc + k);
        g = -re;
      }
    }
}

// ------------------------------------------------------------------------------------------- Brune <-> formal

LevelSet formal_from_brune(const System &system, std::size_t b, const std::vector<double> &E,
                           const std::vector<double> &gamma, const std::vector<double> &B, const ChannelFunctions &cf) {
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
  // Each level keeps its position and receives the eigenvalue of the rank of its Brune energy (CNuc.cpp 1565-1600),
  // and the orientation of its Brune level: the eigenproblem leaves the sign of an eigenvector free, and the one taken
  // makes the coefficient a_{lambda k} of the level's own Brune state positive.  The formal amplitudes then carry the
  // signs of the input wherever the levels mix weakly, and a level that other levels' photon channels end on (an
  // external-capture final state) keeps the sign those amplitudes were given relative to it.
  const std::vector<int> rank = AscendingRanks(E);
  out.E.assign(n, 0.0);
  out.gamma.assign(n * nc, 0.0);
  for (std::size_t la = 0; la < n; la++) {
    const std::size_t k = (std::size_t)rank[la];
    const double orientation = (eig.vectors[la * n + k] < 0.0) ? -1.0 : 1.0;
    out.E[la] = eig.values[k];
    for (std::size_t c = 0; c < nc; c++) {
      double sum = 0.0;
      for (std::size_t mu = 0; mu < n; mu++) sum += eig.vectors[mu * n + k] * gamma[mu * nc + c];
      out.gamma[la * nc + c] = orientation * sum;
    }
  }
  return out;
}

namespace {
/// brune_from_formal in an evaluation, a level that cannot be transformed: NaN energy, amplitudes and row of T, so that
/// every observable that uses it is not finite (level_parameters.md section 9).
void UndefinedLevel(std::vector<double> &tempE, std::vector<double> &tempG, std::vector<double> &W, std::size_t n,
                    std::size_t nc, std::size_t me) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  tempE[me] = nan;
  for (std::size_t k = 0; k < nc; k++) tempG[me * nc + k] = nan;
  if (!W.empty())
    for (std::size_t q = 0; q < n; q++) W[me * n + q] = nan;
}

// A level whose formal values are kept is its own formal level: its row of the transform is a unit vector.
void KeepFormalRow(std::vector<double> &W, std::size_t n, std::size_t me) {
  for (std::size_t q = 0; q < n; q++) W[me * n + q] = (q == me) ? 1.0 : 0.0;
}

/// g(x) = x - e_r(x), e_r the eigenvalue of rank r of A(x) = diag(E) - sum_c (S_c(x) - B_c) gamma_c gamma_c^T, the
/// formal levels in rank order (level_parameters.md section 10).  NaN when S or the eigenproblem fails.
struct RootFunction {
  const System &system;
  std::size_t b;
  const std::vector<double> &E, &gamma, &B;
  const ChannelFunctions &cf;
  const std::vector<bool> &particle;
  std::size_t r;
  EigenResult eig;
  std::vector<double> shift;
  int calls = 0;

  double operator()(double x) {
    calls++;
    const std::size_t n = E.size(), nc = B.size();
    shift.assign(nc, 0.0);
    for (std::size_t k = 0; k < nc; k++)
      if (particle[k]) shift[k] = cf(b, k, x).S - B[k];
    std::vector<double> C(n * n);
    for (std::size_t mu = 0; mu < n; mu++)
      for (std::size_t mup = 0; mup < n; mup++) {
        double chanSum = 0.0;
        for (std::size_t k = 0; k < nc; k++)
          if (particle[k]) chanSum += shift[k] * gamma[mu * nc + k] * gamma[mup * nc + k];
        C[mu * n + mup] = (mu == mup) ? E[mu] - chanSum : -chanSum;
      }
    for (double c : C)
      if (!std::isfinite(c)) return std::numeric_limits<double>::quiet_NaN();
    eig = symmetric_eigen(C, n);
    return eig.ok ? x - eig.values[r] : std::numeric_limits<double>::quiet_NaN();
  }
};

/// The root of g by bisection on a bracket: the two iterates if g changes sign between them, else expanded from x0
/// (the formal energy).  On success f.eig and f.shift belong to the root.
bool BracketedRoot(RootFunction &f, double xa, double xb, double x0, double &root) {
  double lo = 0.0, hi = 0.0, glo = 0.0, ghi = 0.0;
  bool found = false;
  if (std::isfinite(xa) && std::isfinite(xb)) {
    const double ga = f(xa), gb = f(xb);
    if (std::isfinite(ga) && std::isfinite(gb) && (ga <= 0.0) != (gb <= 0.0)) {
      found = true;
      if (ga <= 0.0) lo = xa, glo = ga, hi = xb, ghi = gb;
      else lo = xb, glo = gb, hi = xa, ghi = ga;
    }
  }
  if (!found) {
    double x = x0, g = f(x0), width = 1e-3;
    if (!std::isfinite(g)) return false;
    const double dir = (g < 0.0) ? 1.0 : -1.0;
    for (int i = 0; i < 64 && !found; i++, width *= 2.0) {
      const double x2 = x + dir * width, g2 = f(x2);
      if (!std::isfinite(g2)) return false;
      if ((g2 <= 0.0) != (g <= 0.0)) {
        found = true;
        if (dir > 0.0) lo = x, glo = g, hi = x2, ghi = g2;
        else lo = x2, glo = g2, hi = x, ghi = g;
      }
      x = x2, g = g2;
    }
    if (!found) return false;
  }
  // lo: g <= 0, hi: g > 0.  Bisection to full double precision.
  for (int i = 0; i < 2100; i++) {
    const double mid = lo + 0.5 * (hi - lo);
    if (mid == lo || mid == hi) break;
    const double gm = f(mid);
    if (!std::isfinite(gm)) return false;
    if (gm <= 0.0) lo = mid, glo = gm;
    else hi = mid, ghi = gm;
  }
  root = (std::fabs(glo) <= std::fabs(ghi)) ? lo : hi;
  return std::isfinite(f(root));
}
}  // namespace

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
  out.status.assign(n, LevelSet::Status::Converged);
  out.residual_MeV.assign(n, 0.0);
  out.boundary_shift.assign(n * nc, 0.0);
  if (options.want_transform) out.transform.assign(n * n, 0.0);
  const std::vector<int> rank = AscendingRanks(E);
  // CNuc.cpp TransformOut (lines 2563-2745 at e54e6fcc), one level at a time.
  for (std::size_t lev = 0; lev < n; lev++) {
    std::vector<double> tempE(n);
    std::vector<double> tempG(n * nc);
    for (std::size_t i = 0; i < n; i++) {  // ascending formal energy
      tempE[rank[i]] = E[i];
      for (std::size_t k = 0; k < nc; k++) tempG[rank[i] * nc + k] = gamma[i * nc + k];
    }
    const std::vector<double> rankedE = tempE, rankedG = tempG;  // the formal levels in rank order
    const std::size_t me = (std::size_t)rank[lev];
    // W: the rows of tempG as combinations of the formal levels in rank order (tempG = W gamma_ranked), carried only
    // when the transform is asked for (external_capture.md section 3).
    std::vector<double> W;
    if (options.want_transform) {
      W.assign(n * n, 0.0);
      for (std::size_t i = 0; i < n; i++) W[i * n + i] = 1.0;
    }
    std::vector<double> tempB = B;
    int iteration = 1, reportedIterations = -1;
    // Section 10: the root of E - e(E) on a bracket (xa, xb, else expanded from the formal energy); the level takes the
    // eigenvector of A(root) in the formal levels, oriented on its own level.  false: no bracket, nothing changed.
    auto bracketed = [&](double xa, double xb, std::ostringstream &o) {
      RootFunction g{system, b, rankedE, rankedG, B, cf, particle, me, {}, {}};
      double root = 0.0;
      if (!options.bracket || !BracketedRoot(g, xa, xb, rankedE[me], root)) return false;
      out.status[lev] = LevelSet::Status::Bracketed;
      out.residual_MeV[lev] = std::fabs(root - g.eig.values[me]);
      for (std::size_t k = 0; k < nc; k++) out.boundary_shift[lev * nc + k] = g.shift[k];
      const double orientation = (g.eig.vectors[me * n + me] < 0.0) ? -1.0 : 1.0;
      tempE = g.eig.values;
      for (std::size_t k = 0; k < nc; k++) {
        double sum = 0.0;
        for (std::size_t mup = 0; mup < n; mup++) sum += g.eig.vectors[mup * n + me] * rankedG[mup * nc + k];
        tempG[me * nc + k] = orientation * sum;
      }
      if (!W.empty())
        for (std::size_t q = 0; q < n; q++) W[me * n + q] = orientation * g.eig.vectors[q * n + me];
      reportedIterations = iteration + g.calls;
      o << "; the fixed point repels here: the level is the root of E - e(E), bracketed and bisected (" << g.calls
        << " evaluations, residual " << out.residual_MeV[lev] << " MeV)";
      return true;
    };
    bool done = false;
    while (iteration <= options.max_iterations && !done) {
      std::vector<double> diff(nc, 0.0);
      for (std::size_t k = 0; k < nc; k++)
        if (particle[k]) {
          const double newB = cf(b, k, tempE[me]).S;
          diff[k] = newB - tempB[k];
          tempB[k] = newB;
        }
      for (std::size_t k = 0; k < nc; k++) out.boundary_shift[lev * nc + k] = diff[k];
      std::vector<double> C(n * n);
      for (std::size_t mu = 0; mu < n; mu++)
        for (std::size_t mup = 0; mup < n; mup++) {
          double chanSum = 0.0;
          for (std::size_t k = 0; k < nc; k++)
            if (particle[k]) chanSum += diff[k] * tempG[mu * nc + k] * tempG[mup * nc + k];
          C[mu * n + mup] = (mu == mup) ? tempE[mu] - chanSum : -chanSum;
        }
      const EigenResult eig = symmetric_eigen(C, n);
      bool finite = true;
      for (double c : C) finite = finite && std::isfinite(c);
      if (!eig.ok && !finite) {
        // Section 10, rule 5: the iterates overflowed (a repelling fixed point), the bracket comes from the formal energy.
        std::ostringstream o;
        o << Label(system, b, lev) << ": the iteration diverged after " << iteration << " iterations";
        const double nan = std::numeric_limits<double>::quiet_NaN();
        if (bracketed(nan, nan, o)) {
          out.warnings.push_back(o.str());
          break;
        }
      }
      if (!eig.ok) {
        out.status[lev] = LevelSet::Status::EigenFailed;
        if (options.evaluation) {
          out.warnings.push_back(Label(system, b, lev) + ": could not transform the level (eigenvalue problem failed); "
                                                         "the model is undefined here (NaN)");
          UndefinedLevel(tempE, tempG, W, n, nc, me);
        } else {
          out.warnings.push_back(Label(system, b, lev) + ": could not transform the level (eigenvalue problem failed); "
                                                         "the formal values are kept");
          tempE[me] = E[lev];
          for (std::size_t k = 0; k < nc; k++) tempG[me * nc + k] = gamma[lev * nc + k];
          if (!W.empty()) KeepFormalRow(W, n, me);
        }
        break;
      }
      const double step = std::fabs(eig.values[me] - tempE[me]);  // MeV
      if (step <= std::max(options.energy_tolerance_MeV, options.relative_tolerance * std::fabs(tempE[me]))) done = true;
      // Each eigenvector oriented like formal_from_brune's: its component on its own level positive.
      std::vector<double> newG(n * nc);
      for (std::size_t mu = 0; mu < n; mu++) {
        const double orientation = (eig.vectors[mu * n + mu] < 0.0) ? -1.0 : 1.0;
        for (std::size_t k = 0; k < nc; k++) {
          double s = 0.0;
          for (std::size_t mup = 0; mup < n; mup++) s += eig.vectors[mup * n + mu] * tempG[mup * nc + k];
          newG[mu * nc + k] = orientation * s;
        }
      }
      if (!W.empty()) {
        std::vector<double> newW(n * n);
        for (std::size_t mu = 0; mu < n; mu++) {
          const double orientation = (eig.vectors[mu * n + mu] < 0.0) ? -1.0 : 1.0;
          for (std::size_t q = 0; q < n; q++) {
            double s = 0.0;
            for (std::size_t mup = 0; mup < n; mup++) s += eig.vectors[mup * n + mu] * W[mup * n + q];
            newW[mu * n + q] = orientation * s;
          }
        }
        W.swap(newW);
      }
      const double previous = tempE[me];
      tempE = eig.values;
      tempG.swap(newG);
      if (!done) {
        if (iteration == options.max_iterations) {
          out.residual_MeV[lev] = step;
          std::ostringstream o;
          o << Label(system, b, lev) << ": not transformed after " << options.max_iterations << " iterations (last step "
            << step << " MeV)";
          const bool atFloor = options.evaluation && step <= options.rounding_bound_MeV;
          if (!atFloor && bracketed(previous, tempE[me], o)) {
          } else if (atFloor) {
            // At the rounding floor: the last iterate is accurate to the step and continuous in the parameters.
            out.status[lev] = LevelSet::Status::RoundingLimited;
            o << "; the last iterate is used";
          } else if (options.evaluation) {
            out.status[lev] = LevelSet::Status::NotConverged;
            o << "; the model is undefined here (NaN)";
            UndefinedLevel(tempE, tempG, W, n, nc, me);
          } else {
            out.status[lev] = LevelSet::Status::NotConverged;
            o << "; the formal values are kept";
            tempE[me] = E[lev];
            for (std::size_t k = 0; k < nc; k++) tempG[me * nc + k] = gamma[lev * nc + k];
            if (!W.empty()) KeepFormalRow(W, n, me);
          }
          out.warnings.push_back(o.str());
        }
        iteration++;
      }
    }
    out.E[lev] = tempE[me];
    for (std::size_t k = 0; k < nc; k++) out.gamma[lev * nc + k] = tempG[me * nc + k];
    out.iterations[lev] = (reportedIterations >= 0) ? reportedIterations : iteration;
    if (!W.empty())
      for (std::size_t i = 0; i < n; i++) out.transform[lev * n + i] = W[me * n + (std::size_t)rank[i]];
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

  // 3. External capture, 4. B_c (lowest level).
  if (options.external) subtract_external_capture(system, b, g, *options.external, out.warnings);
  out.basis_E = levels.E;
  out.basis_gamma = g;
  out.B = boundary_values(system, b, levels.E, cf);

  // 5. Formal parameters.
  if (!options.formal) return out;
  if (out.basis == LevelBasis::Formal) {
    out.E = levels.E;
    out.gamma = g;
    return out;
  }
  const std::vector<double> brune = (out.basis == LevelBasis::Park) ? brune_from_park(system, b, levels.E, g, cf) : g;
  LevelSet f = formal_from_brune(system, b, levels.E, brune, out.B, cf);
  out.E.swap(f.E);
  out.gamma.swap(f.gamma);
  return out;
}

}  // namespace physics
}  // namespace azure2
