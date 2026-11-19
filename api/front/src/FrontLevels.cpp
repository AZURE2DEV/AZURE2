// The facade's part for the Levels and Particle Pairs tabs (docs/upgrade/qt_gui.md section 11): the columns of a
// <levels> row by name, what a row's width column is and its channel's quantities at the level energy, a channel's
// radius from <channelRadii>, and what names a pair.  Every number comes from the engine's own functions: the pair
// the engine builds from a row (PPair(io::LevelRow)), ParticleChannelValues (src/EngineChannelFunctions.cpp, the one
// implementation of a particle channel's functions), physics::wigner_limit and single_particle_width
// (physics/Widths.h), experiment::ChannelRadiusFor and the engine's build of the records.
#include "azure2/front/Engine.h"

#include <cmath>
#include <limits>
#include <sstream>

#include "EngineChannelFunctions.h"
#include "PPair.h"
#include "azure2/experiment/Build.h"
#include "azure2/experiment/ChannelRadii.h"
#include "azure2/experiment/ProjectReader.h"
#include "azure2/io/LevelRowText.h"
#include "azure2/io/LevelRows.h"
#include "azure2/physics/Widths.h"

namespace azure2 {
namespace front {

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();

/// The row as the engine's reader takes it: the tokens written, joined by blanks (an optional column the row does not
/// write is "", Project::Levels).
bool ReadRow(const LevelRow &row, std::unique_ptr<io::LevelRow> &out, std::string &why) {
  std::string text;
  for (const std::string &t : row.tokens) {
    if (t.empty()) continue;
    if (!text.empty()) text += ' ';
    text += t;
  }
  std::istringstream in(text);
  out = std::make_unique<io::LevelRow>(in);
  if (in.fail()) {
    why = "the row does not read as a <levels> row";
    return false;
  }
  return true;
}

void SetKind(ChannelQuantities &q, const char *kind, const char *symbol, const char *unit, const char *meaning) {
  q.kind = kind;
  q.symbol = symbol;
  q.unit = unit;
  q.meaning = meaning;
}

}  // namespace

const std::vector<std::string> &LevelColumnNames() { return io::LevelColumnNames(); }

int LevelColumnIndex(const std::string &name) { return io::LevelColumnIndex(name); }

ChannelQuantities ChannelAtLevel(const LevelRow &row, double radius_fm, bool relativistic) {
  ChannelQuantities q;
  q.channelEnergy = q.P = q.S = q.dSdE = q.wignerLimit = q.singleParticleWidth = q.reducedWidth = q.wignerRatio = kNaN;
  std::unique_ptr<io::LevelRow> line;
  if (!ReadRow(row, line, q.error)) return q;
  const io::LevelRow &r = *line;
  q.value = r.gamma();
  q.particle = r.pType() == 0;
  q.radius = radius_fm > 0 ? radius_fm : r.chRad();
  // E_c = (E_x - S) - e, as the engine forms it at a level (EngineChannelFunctions) and at a data point
  q.channelEnergy = r.levelE() - r.sepE() - (r.e2() + r.e1());
  if (r.pType() == 20) {
    SetKind(q, "beta", "B", "", "beta-decay feeding amplitude (Brune normalization)");
    return q;
  }
  if (r.pType() == 10) {
    // A photon channel to the level's own state is a static moment (src/EngineChannelFunctions.cpp: |E - E_f| < 1 keV,
    // the same J and parity); its multipole is M1 or E2 by the parities (the engine's E/M rule).
    const int L = r.l();  // io::LevelRow halves the l2 column
    const bool electric = (L % 2 != 0) == (r.pi2() * r.levelPi() == -1);
    const bool own = std::fabs(r.levelE() - r.e2()) < 1.e-3 && r.levelJ() == r.j2() && r.levelPi() == r.pi2();
    if (own && !electric && L == 1)
      SetKind(q, "magnetic_moment", "μ", "μ_N", "static magnetic dipole moment of the level");
    else if (own && electric && L == 2)
      SetKind(q, "quadrupole_moment", "Q", "b", "static electric quadrupole moment of the level");
    else
      SetKind(q, "photon", "Γ", "eV", "partial gamma width");
    return q;
  }
  if (r.pType() != 0) {
    q.error = "pair type " + std::to_string(r.pType()) + " is none the engine knows";
    return q;
  }
  if (r.gammaIsRWA() == 1)
    SetKind(q, "amplitude", "γ", "MeV^1/2", "reduced-width amplitude (gammaIsRWA)");
  else if (q.channelEnergy > 0)
    SetKind(q, "open", "Γ", "eV", "partial width of the open channel");
  else if (q.channelEnergy < 0)
    SetKind(q, "closed", "C", "fm^-1/2", "asymptotic normalization coefficient (ANC) of the closed channel");
  else
    SetKind(q, "threshold", "Γ", "eV", "the level is at the channel's threshold: no width or ANC is defined");

  PairFunctionSettings settings;
  settings.relativistic = relativistic;
  const PPair pair(r, settings);
  const int twoJ = pair.SpinOrbitTwoJ(r.levelJ());
  const int l = r.l();  // io::LevelRow halves the l2 column
  const azure2::physics::ChannelValues v = ParticleChannelValues(pair, l, twoJ, q.radius, q.channelEnergy);
  q.P = v.P;
  q.S = v.S;
  q.dSdE = v.dS_dE;
  q.wignerLimit = azure2::physics::wigner_limit(pair.GetRedMass(), q.radius);
  if (q.channelEnergy > 0) {
    azure2::physics::ChannelKinematics kin;
    kin.z1 = r.z1();
    kin.z2 = r.z2();
    kin.mu_u = pair.GetRedMass();
    q.singleParticleWidth = azure2::physics::single_particle_width(kin, l, q.radius, q.channelEnergy);
  }
  // the observed reduced width: Gamma = 2 P gamma^2 (Gamma in MeV), C^2 = 2 P gamma^2 below threshold (the factor P of
  // ChannelValues for each, physics/LevelParameters.h); an amplitude is gamma itself
  if (q.kind == "amplitude") q.reducedWidth = q.value * q.value;
  else if (q.kind == "open" && q.P > 0) q.reducedWidth = std::fabs(q.value) * 1.e-6 / (2. * q.P);
  else if (q.kind == "closed" && q.P > 0) q.reducedWidth = q.value * q.value / (2. * q.P);
  if (std::isfinite(q.reducedWidth) && q.wignerLimit > 0) q.wignerRatio = q.reducedWidth / q.wignerLimit;
  return q;
}

ChannelRadiusInfo ChannelRadius(const Project &project, int pairKey, int l, int twoS, int twoJ, double pairRadius) {
  ChannelRadiusInfo out;
  out.radius = pairRadius;
  if (project.Records("channelRadii").empty()) return out;
  experiment::Format3ReadOptions read;
  read.fileName = project.Path();
  read.projectDir = project.Directory();
  read.readDataFiles = false;
  read.dataFiles = experiment::DataFiles::NotUsed;
  read.setUpChecks = false;
  const experiment::Format3Project built = experiment::ReadFormat3Text(project.Text(), read);
  const experiment::ChannelRadiusEntry *e =
      experiment::ChannelRadiusFor(built.built.channelRadii, pairKey, l, twoS, twoJ, &out.ambiguous);
  if (e) {
    out.radius = e->radius_fm;
    out.own = true;
    out.record = e->record;
  }
  return out;
}

std::vector<std::string> PairUses(const Project &project, int pairKey) {
  std::vector<std::string> out;
  experiment::Format3ReadOptions read;
  read.fileName = project.Path();
  read.projectDir = project.Directory();
  read.readDataFiles = false;
  read.dataFiles = experiment::DataFiles::NotUsed;
  read.setUpChecks = false;
  read.allRecordsActive = true;
  const experiment::Format3Project built = experiment::ReadFormat3Text(project.Text(), read);
  auto names = [&project](const std::string &section, const std::string &kind) {
    std::vector<std::string> n;
    for (const Record &r : project.Records(section))
      if (r.kind == kind) n.push_back(r.name);
    return n;
  };
  const std::vector<std::string> data = names("data", "dataset"), extrapolations = names("extrapolations", "extrapolation");
  for (const auto &i : built.built.data.datasets.items) {
    const bool names = i.spec.reaction.entranceKey == pairKey || i.spec.reaction.exitKey == pairKey;
    if (names && i.spec.key >= 1 && i.spec.key <= (int)data.size()) out.push_back("data: " + data[i.spec.key - 1]);
  }
  for (const auto &i : built.built.extrapolations.extrapolations.items) {
    const bool names = i.spec.reaction.entranceKey == pairKey || i.spec.reaction.exitKey == pairKey;
    if (names && i.spec.key >= 1 && i.spec.key <= (int)extrapolations.size())
      out.push_back("extrapolations: " + extrapolations[i.spec.key - 1]);
  }
  return out;
}

}  // namespace front
}  // namespace azure2
