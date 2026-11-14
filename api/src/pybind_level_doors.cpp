// The two doors of the level parameters for pyazr.endf (docs/upgrade/endf.md section 8.2): one J^pi block of
// particle channels, given by plain pair records and numbers, converted by the engine's own functions of
// physics/LevelParameters.h with the engine's channel functions (ChannelFunc, WhitFunc on PPair objects built from
// the records, as CompoundSystem::Fill builds them).
//
//   formal_levels_of_block    the block's levels in the project's basis (Brune or Park amplitudes, or the observed
//                             widths and ANCs) -> the formal levels with the engine's B_c = S_c(E_1)
//                             (level_parameters: what a session with use_brune=False holds);
//   observed_levels_of_block  formal levels with the caller's B_c (any boundary condition, e.g. an ENDF file's) ->
//                             the Brune levels and their observed widths and ANCs (brune_from_formal, then
//                             observed_levels).
// use_coulomb_table and use_gsl_coulomb are the session options of the same names (the engine's defaults: true,
// false), so that the levels are those of a session with the same options.
//
// Only particle channels: photon and beta channels do not enter the level matrix (LevelParameters.h,
// boundary_values) and are not asked for.  No physics in this file: it builds the System and the PPair objects and
// calls physics/.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "CoulombTable.h"
#include "EngineChannelFunctions.h"
#include "PPair.h"
#include "azure2/io/LevelRows.h"
#include "azure2/physics/LevelParameters.h"
#include "azure2/physics/System.h"
#include "azure2/support/Environment.h"

namespace py = pybind11;
namespace phys = azure2::physics;

namespace {

/// A particle pair as pyazr's pair records give it (Session.pairs, pyazr.kinds()["pair_fields"]).
struct PairFields {
  int key = 0;
  int z1 = 0, z2 = 0;
  double m1 = 0., m2 = 0.;          // atomic masses (u)
  double j1 = 0., j2 = 0.;
  int pi1 = 1, pi2 = 1;
  double light_excitation = 0.;     // MeV
  double heavy_excitation = 0.;     // MeV
  double separation = 0.;           // MeV
  double radius = 0.;               // fm
};

double Number(const py::dict &d, const char *name) {
  if (!d.contains(name)) throw std::invalid_argument(std::string("a pair record without the field '") + name + "'");
  return py::cast<double>(d[name]);
}

PairFields PairOf(const py::dict &d) {
  const std::string kind = d.contains("pair_kind") ? py::cast<std::string>(d["pair_kind"]) : "particle";
  PairFields p;
  p.key = static_cast<int>(Number(d, "key"));
  if (kind != "particle")
    throw std::invalid_argument("pair " + std::to_string(p.key) + " is a " + kind +
                                " pair; the level doors take particle pairs only");
  p.z1 = static_cast<int>(Number(d, "Z1"));
  p.z2 = static_cast<int>(Number(d, "Z2"));
  p.m1 = Number(d, "M1");
  p.m2 = Number(d, "M2");
  p.j1 = Number(d, "J1");
  p.j2 = Number(d, "J2");
  p.pi1 = static_cast<int>(Number(d, "parity1"));
  p.pi2 = static_cast<int>(Number(d, "parity2"));
  p.light_excitation = d.contains("light_excitation") ? Number(d, "light_excitation") : 0.0;
  p.heavy_excitation = Number(d, "excitation") - p.light_excitation;
  p.separation = Number(d, "sep_energy");
  p.radius = Number(d, "channel_radius");
  return p;
}

/// The pair as the engine builds it from a <levels> row (PPair(LevelRow, settings)), with the run's default
/// settings of its functions: the Coulomb table where it covers the pair (as a session), the environment's switches.
std::unique_ptr<PPair> EnginePair(const PairFields &p, bool use_coulomb_table) {
  std::ostringstream row;
  row.precision(17);
  // J pi E fix aa ir s l id active chfix gamma j1 pi1 j2 pi2 e2 m1 m2 z1 z2 Sent S j3 pi3 e3 ptype a g1 g2 ecmask e1
  row << 0.5 << ' ' << 1 << ' ' << 0.0 << " 1 1 " << p.key << ' ' << 0.5 << " 0 1 1 1 0.0 " << p.j1 << ' ' << p.pi1
      << ' ' << p.j2 << ' ' << p.pi2 << ' ' << p.heavy_excitation << ' ' << p.m1 << ' ' << p.m2 << ' ' << p.z1 << ' '
      << p.z2 << ' ' << p.separation << ' ' << p.separation << " 0 0 0.0 0 " << p.radius << " 0.0 0.0 0 "
      << p.light_excitation;
  std::istringstream in(row.str());
  azure2::io::LevelRow line(in);
  if (in.fail()) throw std::invalid_argument("pair " + std::to_string(p.key) + ": its fields do not make a row");
  const azure2::support::Environment &env = azure2::support::ProcessEnvironment();
  PairFunctionSettings settings;
  settings.useCoulombTable = use_coulomb_table;
  settings.coulombTableDefault = CoulombTable::DefaultPathIn(env);
  settings.boundShiftNumerical = env.Set(azure2::support::Switch::kBoundShiftNumerical);
  settings.noNeutralClosedForm = env.Set(azure2::support::Switch::kNoNeutralClosedForm);
  settings.shiftDerivative1eV = env.Set(azure2::support::Switch::kShiftDerivative1eV);
  return std::make_unique<PPair>(line, settings);
}

/// One block of particle channels: the System, the engine's pairs of its partitions and the block's J.
struct Block {
  std::vector<std::unique_ptr<PPair>> pairs;  // by partition
  std::unique_ptr<phys::System> system;
  double J = 0.;
  bool gsl = false;  // GSL's Coulomb functions (the session option use_gsl_coul)
};

/// channels: (pair key, l, 2s) per block channel; E (MeV, E_x), gamma row-major (levels x channels), B per channel.
Block MakeBlock(const py::list &pair_records, int two_J, int parity, const std::vector<std::vector<int>> &channels,
                const std::vector<double> &E, const std::vector<double> &gamma, const std::vector<double> &B,
                phys::LevelBasis basis, bool use_coulomb_table, bool use_gsl_coulomb) {
  std::map<int, PairFields> byKey;
  for (const py::handle &h : pair_records) {
    const py::dict d = py::cast<py::dict>(h);
    if (d.contains("pair_kind") && py::cast<std::string>(d["pair_kind"]) != "particle") continue;
    const PairFields p = PairOf(d);
    byKey[p.key] = p;
  }
  Block out;
  out.J = two_J / 2.0;
  out.gsl = use_gsl_coulomb;
  std::map<int, int> partitionOfKey;
  std::vector<phys::Partition> partitions;
  std::vector<phys::Channel> chans;
  std::vector<int> blockChannels;
  int nucleusZ = -1, nucleusA = -1;
  for (const std::vector<int> &c : channels) {
    if (c.size() != 3) throw std::invalid_argument("a channel is (pair key, l, 2s)");
    auto it = byKey.find(c[0]);
    if (it == byKey.end())
      throw std::invalid_argument("a channel names pair " + std::to_string(c[0]) + ", which is not a particle pair "
                                  "of the records given");
    const PairFields &p = it->second;
    if (!partitionOfKey.count(p.key)) {
      const int z = p.z1 + p.z2, a = static_cast<int>(std::floor(p.m1 + p.m2 + 0.5));
      if (nucleusZ < 0) nucleusZ = z, nucleusA = a;
      if (z != nucleusZ || a != nucleusA)
        throw std::invalid_argument("the channels of one block name pairs of two compound nuclei");
      phys::Partition part;
      part.kind = phys::PartitionKind::Particle;
      part.nucleus = 0;
      part.Z = {{p.z1, p.z2}};
      part.mass_u = {{p.m1, p.m2}};
      part.two_j = {{static_cast<int>(std::lround(2 * p.j1)), static_cast<int>(std::lround(2 * p.j2))}};
      part.parity = {{p.pi1, p.pi2}};
      part.excitation_MeV = {{p.light_excitation, p.heavy_excitation}};
      part.separation_MeV = p.separation;
      partitionOfKey[p.key] = static_cast<int>(partitions.size());
      partitions.push_back(part);
      out.pairs.push_back(EnginePair(p, use_coulomb_table));
    }
    phys::Channel ch;
    ch.partition = partitionOfKey[p.key];
    ch.kind = phys::ChannelKind::Particle;
    ch.l = c[1];
    ch.two_s = c[2];
    ch.radius_fm = p.radius;
    blockChannels.push_back(static_cast<int>(chans.size()));
    chans.push_back(ch);
  }
  const std::size_t nc = channels.size(), nl = E.size();
  if (gamma.size() != nl * nc)
    throw std::invalid_argument("gamma has " + std::to_string(gamma.size()) + " values; levels x channels is " +
                                std::to_string(nl * nc));
  if (B.size() != nc) throw std::invalid_argument("B has one value per channel");
  phys::JpiBlock block(two_J, parity, 0, blockChannels, B, E, gamma);
  phys::Convention convention;
  convention.basis = basis;
  out.system = std::make_unique<phys::System>(convention, std::vector<phys::Nucleus>{{nucleusZ, nucleusA}}, partitions,
                                              chans, std::vector<phys::JpiBlock>{block});
  return out;
}

/// The engine's channel functions of a particle channel at a level energy E (MeV, E_x), as the engine's adapter
/// EngineChannelFunctions (src/EngineChannelFunctions.cpp) forms them: E_c = (E - S) - e; S, dS/dE; P above threshold,
/// mu a / (hbar^2 W^2) below it, 0 at threshold.
phys::ChannelFunctions ParticleChannelFunctions(const Block &blk) {
  const Block *b = &blk;
  return [b](std::size_t, std::size_t k, double E) {
    const phys::Channel &ch = b->system->channel(b->system->block(0).channel(k));
    const PPair *pair = b->pairs.at(ch.partition).get();
    const double Ec = E - pair->GetSepE() - pair->GetExE();
    const int l = ch.l;
    const int twoJ = pair->SpinOrbitTwoJ(b->J);
    return ParticleChannelValues(*pair, l, twoJ, Ec, b->gsl);
  };
}

phys::LevelBasis BasisOf(const std::string &name) {
  if (name == "brune") return phys::LevelBasis::Brune;
  if (name == "park") return phys::LevelBasis::Park;
  throw std::invalid_argument("basis is 'brune' or 'park', not '" + name + "'");
}

const char *StatusName(phys::LevelSet::Status s) {
  switch (s) {
    case phys::LevelSet::Status::Converged: return "converged";
    case phys::LevelSet::Status::EigenFailed: return "eigen_failed";
    case phys::LevelSet::Status::NotConverged: return "not_converged";
    case phys::LevelSet::Status::RoundingLimited: return "rounding_limited";
    case phys::LevelSet::Status::Bracketed: return "bracketed";
  }
  return "unknown";
}

const char *KindName(phys::ObservedKind k) {
  switch (k) {
    case phys::ObservedKind::Width: return "width";
    case phys::ObservedKind::ANC: return "anc";
    case phys::ObservedKind::MagneticMoment: return "magnetic_moment";
    case phys::ObservedKind::QuadrupoleMoment: return "quadrupole_moment";
    case phys::ObservedKind::Beta: return "beta";
  }
  return "unknown";
}

py::dict FormalLevels(const py::list &pairs, int two_J, int parity, const std::vector<std::vector<int>> &channels,
                      const std::vector<double> &E, const std::vector<double> &values,
                      const std::vector<int> &is_amplitude, const std::string &form, const std::string &basis,
                      bool use_coulomb_table, bool use_gsl_coulomb) {
  const phys::LevelBasis b = BasisOf(basis);
  const std::vector<double> zeroB(channels.size(), 0.0);
  const std::vector<double> zeroG(E.size() * channels.size(), 0.0);
  Block blk = MakeBlock(pairs, two_J, parity, channels, E, zeroG, zeroB, b, use_coulomb_table, use_gsl_coulomb);
  phys::ProjectLevels levels;
  if (form == "observed") {
    levels.form = phys::ProjectLevels::Form::Observed;
    levels.is_amplitude.assign(is_amplitude.begin(), is_amplitude.end());
    if (!levels.is_amplitude.empty() && levels.is_amplitude.size() != values.size())
      throw std::invalid_argument("is_amplitude has one flag per value or none");
  } else if (form == "amplitudes") {
    levels.form = phys::ProjectLevels::Form::Amplitudes;
  } else {
    throw std::invalid_argument("form is 'observed' or 'amplitudes', not '" + form + "'");
  }
  levels.E = E;
  levels.value = values;
  const phys::InternalParameters p =
      phys::level_parameters(*blk.system, 0, blk.system->convention(), phys::internal_channels(*blk.system, 0), levels,
                             ParticleChannelFunctions(blk));
  py::dict out;
  out["E"] = p.E;
  out["gamma"] = p.gamma;
  out["B"] = p.B;
  out["basis_E"] = p.basis_E;
  out["basis_gamma"] = p.basis_gamma;
  out["warnings"] = p.warnings;
  return out;
}

py::dict ObservedLevels(const py::list &pairs, int two_J, int parity, const std::vector<std::vector<int>> &channels,
                        const std::vector<double> &E, const std::vector<double> &gamma, const std::vector<double> &B,
                        double energy_tolerance_MeV, bool use_coulomb_table, bool use_gsl_coulomb) {
  Block blk = MakeBlock(pairs, two_J, parity, channels, E, gamma, B, phys::LevelBasis::Brune, use_coulomb_table,
                        use_gsl_coulomb);
  const phys::ChannelFunctions cf = ParticleChannelFunctions(blk);
  phys::FormalToBruneOptions options;
  options.energy_tolerance_MeV = energy_tolerance_MeV;
  const phys::LevelSet brune = phys::brune_from_formal(*blk.system, 0, E, gamma, B, cf, options);
  const std::vector<phys::ProjectChannel> columns = phys::internal_channels(*blk.system, 0);
  const phys::ObservedLevels obs = phys::observed_levels(*blk.system, 0, columns, brune.E, brune.gamma, {}, {},
                                                         phys::LevelBasis::Brune, cf);
  py::dict out;
  out["E"] = brune.E;
  out["gamma"] = brune.gamma;
  out["value"] = obs.value;
  std::vector<std::string> kinds, status;
  for (phys::ObservedKind k : obs.kind) kinds.push_back(KindName(k));
  for (phys::LevelSet::Status s : brune.status) status.push_back(StatusName(s));
  out["kind"] = kinds;
  out["status"] = status;
  out["iterations"] = brune.iterations;
  out["norm"] = obs.norm;
  out["warnings"] = brune.warnings;
  return out;
}

}  // namespace

void BindLevelDoors(py::module_ &m) {
  m.def("formal_levels_of_block", &FormalLevels, py::arg("pairs"), py::arg("two_J"), py::arg("parity"),
        py::arg("channels"), py::arg("E"), py::arg("values"), py::arg("is_amplitude"), py::arg("form"),
        py::arg("basis"), py::arg("use_coulomb_table") = true, py::arg("use_gsl_coulomb") = false,
        "One J^pi block of particle channels in the project's basis -> its formal levels (physics/ level_parameters, "
        "the engine's channel functions).  pairs: pair records (dicts with key, pair_kind, Z1, Z2, M1, M2, J1, J2, "
        "parity1, parity2, excitation, light_excitation, sep_energy, channel_radius); channels: (pair key, l, 2s); E: "
        "level energies (MeV, E_x); values: row-major levels x channels, form 'observed' (widths in eV, ANCs in "
        "fm^-1/2, or amplitudes where is_amplitude) or 'amplitudes' (MeV^1/2 of `basis`, 'brune' or 'park').  Returns "
        "{E, gamma (formal, MeV^1/2), B (S_c at the lowest level), basis_E, basis_gamma, warnings}; ValueError for "
        "records that make no block, std::domain_error (ValueError) when no formal parameters correspond.");
  m.def("observed_levels_of_block", &ObservedLevels, py::arg("pairs"), py::arg("two_J"), py::arg("parity"),
        py::arg("channels"), py::arg("E"), py::arg("gamma"), py::arg("B"),
        py::arg("energy_tolerance_MeV") = phys::FormalToBruneOptions().energy_tolerance_MeV,
        py::arg("use_coulomb_table") = true, py::arg("use_gsl_coulomb") = false,
        "One J^pi block of particle channels given by formal levels with the caller's boundary values B_c -> the "
        "Brune levels (brune_from_formal, its iteration stopped when a level moves by <= energy_tolerance_MeV, the "
        "engine's 1e-6 MeV by default) and their observed values (observed_levels): {E (MeV, E_x), gamma (Brune, "
        "MeV^1/2), value (row-major: width in eV for an open channel, ANC in fm^-1/2 for a closed one, signed), kind "
        "('width' or 'anc' per value), status per level (converged, bracketed, not_converged, ...), iterations, norm "
        "(Brune's N per level), warnings}.");
}
