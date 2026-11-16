// A level's radiative decay as Gamma_gamma, b_f, delta_f between the parameter vector and the engine's amplitudes
// (GammaDecayMap.h; docs/upgrade/gamma_parametrisation.md sections 2.2-2.6).
#include "GammaDecayMap.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <sstream>

#include "AChannel.h"
#include "ALevel.h"
#include "AZUREGrad.h"
#include "AZUREParams.h"
#include "CompoundSystem.h"
#include "ExternalCaptureAmplitudes.h"
#include "EngineChannelFunctions.h"
#include "ChannelFunc.h"
#include "Config.h"
#include "JGroup.h"
#include "PPair.h"
#include "azure2/experiment/LevelTables.h"
#include "azure2/experiment/ProjectReader.h"
#include "azure2/experiment/Schema.h"
#include "azure2/io/LevelRows.h"
#include "azure2/io/RecordReader.h"
#include "azure2/physics/ExternalCapture.h"
#include "azure2/physics/LevelParameters.h"

namespace ex = azure2::experiment;
namespace phys = azure2::physics;

namespace {

const char *const kSection = "gammaDecays";
/// Width of the wall on the sum of the explicit branchings above 1 (section 2.6).
constexpr double kWallWidth = 1e-3;
/// Step of the central difference of the photon amplitudes in the level energy (MeV), Eq. 2.5.
constexpr double kEnergyStep = 1e-4;
constexpr double kEVPerMeV = 1e6;
/// A branching below this is reported at its bound 0.
constexpr double kAtBound = 1e-6;

void Print(const Config &configure, const ex::Diagnostic &d) {
  configure.outStream << (d.severity == ex::Severity::Warning ? "WARNING: " : "ERROR: ") << "<gammaDecays> line " << d.line
                      << " [" << d.code << "]: " << d.message << std::endl;
}

bool IsPhotonChannel(const CompoundSystem &compound, const JGroup *group, int ch) {
  return compound.GetPair(group->GetChannel(ch)->GetPairNum())->GetPType() == 10;
}

bool IsParticleChannel(const CompoundSystem &compound, const JGroup *group, int ch) {
  return compound.GetPair(group->GetChannel(ch)->GetPairNum())->GetPType() == 0;
}

/// The full-vector index of level la's energy in group j; its amplitudes follow it, in channel order
/// (LevelParameterLayout.h).
int LevelEnergyIndex(const CompoundSystem &compound, int j, int la) {
  return compound.ParameterLayout().Of(j, la).energy;
}

/// The numbers of record level L as the file gives them, slot order (Gamma_gamma in eV).
std::vector<double> SpecNumbers(const GammaDecayLevel &L) {
  std::vector<double> q;
  for (const GammaDecayLevel::Slot &slot : L.slots) {
    switch (slot.part) {
      case GammaDecayLevel::Part::TotalWidth: q.push_back(L.spec.totalWidth.value * kEVPerMeV); break;
      case GammaDecayLevel::Part::Branching: q.push_back(L.spec.finals[slot.final].branching.value); break;
      case GammaDecayLevel::Part::Mixing: q.push_back(L.spec.finals[slot.final].mixing.value); break;
    }
  }
  return q;
}

std::string Fmt(double x) {
  char buf[32];
  std::snprintf(buf, sizeof buf, "%.6g", x);
  return buf;
}

}  // namespace

std::shared_ptr<const GammaDecaySet> GammaDecaySet::Read(const CompoundSystem &compound, const Config &configure,
                                                         const std::vector<azure2::io::LevelRow> &rows) {
  auto set = std::make_shared<GammaDecaySet>();
  std::unique_ptr<std::istream> section = configure.SectionStream(kSection);
  if (!section) return set;
  std::vector<std::string> lines;
  std::string line;
  std::getline(*section, line);  // <gammaDecays>
  while (std::getline(*section, line) && line.find("</gammaDecays>") == std::string::npos) lines.push_back(line);
  const azure2::io::RecordReadResult read = azure2::io::ReadRecords(lines, 1, ex::Schema::KindsOfSection(kSection));
  std::vector<ex::Diagnostic> diagnostics;
  for (const azure2::io::RecordDiagnostic &d : read.diagnostics) diagnostics.push_back(ex::FromRecordDiagnostic(d, kSection));
  const ex::ValidationResult validated = ex::Schema::Validate({{kSection, read.set.records}});
  diagnostics.insert(diagnostics.end(), validated.diagnostics.begin(), validated.diagnostics.end());
  bool ok = true;
  for (const ex::Diagnostic &d : diagnostics) ok = ok && d.severity == ex::Severity::Warning;
  ex::GammaDecayBuild built;
  if (ok) {
    built = ex::BuildGammaDecays(validated, ex::LevelTablesOfRows(rows));
    diagnostics.insert(diagnostics.end(), built.diagnostics.begin(), built.diagnostics.end());
    ok = built.Ok();
  }
  for (const ex::Diagnostic &d : diagnostics) Print(configure, d);
  if (!ok) return nullptr;
  if (built.specs.empty()) return set;

  auto refuse = [&](const ex::GammaDecaySpec &s, const std::string &why) {
    configure.outStream << "ERROR: <gammaDecays> record '" << s.name << "' (level " << s.levelId << "): " << why << std::endl;
    return nullptr;
  };
  set->park_ = (configure.paramMask & Config::USE_PARK_FORMALISM);
  set->gslCoulomb_ = (configure.paramMask & Config::USE_GSL_COULOMB_FUNC);
  set->rmc_ = (configure.paramMask & Config::USE_RMC_FORMALISM);
  const bool formal = !(configure.paramMask & Config::USE_BRUNE_FORMALISM);
  const bool transform = (configure.paramMask & Config::TRANSFORM_PARAMETERS);
  const bool external = (configure.paramMask & Config::USE_EXTERNAL_CAPTURE);
  const bool ignoreZero = (configure.paramMask & Config::IGNORE_ZERO_WIDTHS);
  set->configure_ = &configure;
  try {
    set->structure_ = compound.SystemStructureHandle();
  } catch (const std::exception &e) {
    configure.outStream << "ERROR: <gammaDecays>: " << e.what() << std::endl;
    return nullptr;
  }
  const SystemOfCompoundMap &map = set->structure_->map;
  for (const ex::GammaDecaySpec &s : built.specs) {
    if (formal)
      return refuse(s, "the formal basis (--no-brune) mixes the photon amplitudes of all levels of a J^pi group, so a "
                       "level's total radiative width and branchings cannot be its own parameters; use the Brune or "
                       "Park basis (gamma_parametrisation.md section 2.4)");
    if (!transform)
      return refuse(s, "--no-transform reads the level parameters as amplitudes; the record gives observed widths");
    GammaDecayLevel L;
    L.spec = s;
    const JGroup *group = compound.GetJGroup(s.group);
    const ALevel *level = group->GetLevel(s.level);
    if (!group->IsInRMatrix() || !level->IsInRMatrix()) return refuse(s, "the level is not in the R-matrix");
    for (std::size_t b = 0; b < map.groupOfBlock.size(); b++)
      if (map.groupOfBlock[b] == s.group) L.block = b;
    L.channelOfBlockChannel = map.channelOfBlockChannel.at(L.block);
    if ((int)L.channelOfBlockChannel.size() != group->NumChannels())
      return refuse(s, "the level's J^pi group has channels outside the R-matrix");
    for (int ch = 1; ch <= group->NumChannels(); ch++)
      if (IsPhotonChannel(compound, group, ch)) L.photonChannels.push_back(ch);
    L.energyIndex = LevelEnergyIndex(compound, s.group, s.level);
    L.levelNumber = s.level;
    for (int j = 1; j < s.group; j++) L.levelNumber += compound.GetJGroup(j)->NumLevels();
    for (const int ch : L.photonChannels) {
      const AChannel *c = group->GetChannel(ch);
      const PPair *pair = compound.GetPair(c->GetPairNum());
      if (std::fabs(level->FileE() - pair->GetHeavyExE()) < 1.e-3 && group->GetJ() == pair->GetJ(2) && group->GetPi() == pair->GetPi(2))
        return refuse(s, "a photon channel to the level's own state is a static moment, not a decay");
    }
    // The record's numbers, in slot order: Gamma_gamma, the branchings but the remainder, the mixing ratios.
    L.slots.push_back({GammaDecayLevel::Part::TotalWidth, -1});
    for (std::size_t f = 0; f < s.finals.size(); f++)
      if ((int)f != s.remainder) L.slots.push_back({GammaDecayLevel::Part::Branching, (int)f});
    for (std::size_t f = 0; f < s.finals.size(); f++)
      if (s.finals[f].L.size() == 2) L.slots.push_back({GammaDecayLevel::Part::Mixing, (int)f});
    if (L.slots.size() != L.photonChannels.size())
      return refuse(s, "internal error: " + std::to_string(L.slots.size()) + " numbers for " +
                           std::to_string(L.photonChannels.size()) + " photon channels (report it)");
    // The partial widths of the record, which CompoundSystem::Fill makes the level's input widths (eV, WriteInputWidths).
    const std::vector<double> widths = PartialWidths(L, SpecNumbers(L));
    // The external part of each slot (section 2.7): the final state whose capture covers the channel's multipolarity,
    // none under --ignore-externals for a width given as 0 (TransformIn's rule on the input).
    L.external.assign(L.photonChannels.size(), GammaDecayLevel::External());
    for (std::size_t k = 0; k < L.photonChannels.size() && external; k++) {
      if (ignoreZero && widths[k] == 0.0) continue;
      const std::pair<int, int> f = ExternalCaptureFinalState(compound, group->GetChannel(L.photonChannels[k]));
      if (f.first)
        L.external[k] = {f.first, f.second, LevelEnergyIndex(compound, f.first, f.second),
                         compound.GetJGroup(f.first)->NumChannels()};
    }
    set->levels_.push_back(L);
  }
  return set;
}

int GammaDecaySet::RecordOf(int j, int la) const {
  for (std::size_t r = 0; r < levels_.size(); r++)
    if (levels_[r].spec.group == j && levels_[r].spec.level == la) return (int)r;
  return -1;
}

int GammaDecaySet::SlotOf(int record, int ch) const {
  if (record < 0) return -1;
  const std::vector<int> &pc = levels_[record].photonChannels;
  for (std::size_t k = 0; k < pc.size(); k++)
    if (pc[k] == ch) return (int)k;
  return -1;
}

bool GammaDecaySet::OwnsSlot(int j, int la, int ch) const { return SlotOf(RecordOf(j, la), ch) >= 0; }

bool GammaDecaySet::AddParameter(AZUREParams &params, const std::vector<std::vector<double>> &inputs, int j, int la,
                                 int ch, int levelNumber) const {
  const int r = RecordOf(j, la);
  if (SlotOf(r, ch) < 0) return false;
  const GammaDecayLevel &L = levels_[r];
  const GammaDecayLevel::Slot &slot = L.slots[SlotOf(r, ch)];
  const ex::GammaDecaySpec &s = L.spec;
  ROOT::Minuit2::MnUserParameters &p = params.GetMinuitParams();
  std::string name = s.name, part, final;
  ex::GammaNumber number;
  switch (slot.part) {
    case GammaDecayLevel::Part::TotalWidth:
      part = "total_width";
      number = s.totalWidth;
      break;
    case GammaDecayLevel::Part::Branching:
      part = "branching";
      final = s.finals[slot.final].label;
      number = s.finals[slot.final].branching;
      break;
    case GammaDecayLevel::Part::Mixing:
      part = "mixing";
      final = s.finals[slot.final].label;
      number = s.finals[slot.final].mixing;
      break;
  }
  name += "." + part + (final.empty() ? "" : "." + final);
  const double value = inputs.at(r).at(SlotOf(r, ch));
  const double step = slot.part == GammaDecayLevel::Part::TotalWidth ? (value != 0.0 ? 0.1 * value : 1e-3)
                      : slot.part == GammaDecayLevel::Part::Branching ? 0.05
                                                                      : 0.05 * std::fabs(value) + 0.01;
  params.Add(name, value, step, ParameterRecord::Radiative(j, la, ch, levelNumber, s.name, part, final));
  if (slot.part == GammaDecayLevel::Part::TotalWidth) p.SetLowerLimit(name, 0.0);
  if (slot.part == GammaDecayLevel::Part::Branching) p.SetLimits(name, 0.0, 1.0);
  const bool mixingLine = slot.part != GammaDecayLevel::Part::Mixing || s.finals[slot.final].hasMixing;
  const bool listed = slot.part != GammaDecayLevel::Part::Branching ||
                      s.finals[slot.final].role == ex::GammaFinalState::Role::Listed;
  if (!(number.fit && mixingLine && listed)) p.Fix(name);
  return true;
}

std::vector<std::vector<double>> GammaDecaySet::FileNumbers() const {
  std::vector<std::vector<double>> numbers;
  for (const GammaDecayLevel &L : levels_) numbers.push_back(SpecNumbers(L));
  return numbers;
}

std::vector<std::vector<double>> GammaDecaySet::NumbersAt(const vector_r &p) const {
  std::vector<std::vector<double>> numbers;
  for (const GammaDecayLevel &L : levels_) {
    numbers.emplace_back();
    for (const int ch : L.photonChannels) numbers.back().push_back(p[L.energyIndex + ch]);
  }
  return numbers;
}

void GammaDecaySet::WriteInputWidths(LevelInputs &inputs) const {
  const std::vector<std::vector<double>> &numbers = inputs.GammaDecayNumbers();
  for (std::size_t r = 0; r < levels_.size(); r++) {
    const GammaDecayLevel &L = levels_[r];
    const std::vector<double> widths = PartialWidths(L, numbers[r]);
    for (std::size_t k = 0; k < L.photonChannels.size(); k++)
      inputs.SetGamma(L.spec.group, L.spec.level, L.photonChannels[k], widths[k]);
  }
}

vector_r GammaDecaySet::WithPartialWidths(const vector_r &observed) const {
  vector_r p = observed;
  const std::vector<std::vector<double>> numbers = NumbersAt(observed);
  for (std::size_t r = 0; r < levels_.size(); r++) {
    const GammaDecayLevel &L = levels_[r];
    const std::vector<double> widths = PartialWidths(L, numbers[r]);
    for (std::size_t k = 0; k < L.photonChannels.size(); k++) p[L.energyIndex + L.photonChannels[k]] = widths[k];
  }
  return p;
}

double GammaDecaySet::ValueOr(const CompoundSystem &compound, const LevelState &state, int j, int la, int ch, double otherwise) const {
  const int r = RecordOf(j, la), slot = SlotOf(r, ch);
  if (slot < 0 || !state.HasGammaDecayFill(r)) return otherwise;
  return state.GammaDecayFillOf(r).q.at(slot);
}

std::vector<double> GammaDecaySet::PartialWidths(const GammaDecayLevel &L, const std::vector<double> &q,
                                                 std::vector<double> *dGdq, double *excess) {
  const ex::GammaDecaySpec &s = L.spec;
  const std::size_t m = L.slots.size(), nf = s.finals.size();
  // Gamma_gamma (eV), b_f and delta_f with db_f/dq (the remainder and the scaling of section 2.6).
  double total = 0.0;
  std::vector<double> b(nf, 0.0), delta(nf, 0.0);
  std::vector<int> deltaSlot(nf, -1), bSlot(nf, -1);
  int totalSlot = 0;
  for (std::size_t i = 0; i < m; i++) {
    const GammaDecayLevel::Slot &slot = L.slots[i];
    if (slot.part == GammaDecayLevel::Part::TotalWidth) {
      total = q[i];
      totalSlot = (int)i;
    } else if (slot.part == GammaDecayLevel::Part::Branching) {
      b[slot.final] = q[i];
      bSlot[slot.final] = (int)i;
    } else {
      delta[slot.final] = q[i];
      deltaSlot[slot.final] = (int)i;
    }
  }
  double sum = 0.0;
  for (std::size_t f = 0; f < nf; f++)
    if ((int)f != s.remainder) sum += b[f];
  if (excess) *excess = sum - 1.0;
  std::vector<double> db(nf * m, 0.0);  // db_f / dq_i
  const bool scaled = sum > 1.0;
  for (std::size_t f = 0; f < nf; f++) {
    if ((int)f == s.remainder) continue;
    if (scaled) {
      for (std::size_t g = 0; g < nf; g++)
        if ((int)g != s.remainder) db[f * m + bSlot[g]] = ((f == g ? 1.0 : 0.0) - b[f] / sum) / sum;
    } else {
      db[f * m + bSlot[f]] = 1.0;
    }
  }
  std::vector<double> bEff = b;
  if (scaled)
    for (std::size_t f = 0; f < nf; f++) bEff[f] = (int)f == s.remainder ? 0.0 : b[f] / sum;
  if (s.remainder >= 0) {
    bEff[s.remainder] = scaled ? 0.0 : 1.0 - sum;
    if (!scaled)
      for (std::size_t f = 0; f < nf; f++)
        if ((int)f != s.remainder) db[s.remainder * m + bSlot[f]] = -1.0;
  }
  // Eq. (2.3): Gamma_L = sigma Gamma_gamma b / (1 + delta^2), Gamma_{L+1} = -sigma Gamma_gamma b delta |delta| / (1 +
  // delta^2) (delta of Krane and Steffen: widths of equal sign are delta < 0, section 2.3).
  std::vector<double> width(m, 0.0);
  if (dGdq) dGdq->assign(m * m, 0.0);
  for (std::size_t k = 0; k < m; k++) {
    const int ch = L.photonChannels[k];
    for (std::size_t f = 0; f < nf; f++)
      for (std::size_t c = 0; c < s.finals[f].channels.size(); c++) {
        if (s.finals[f].channels[c] != ch) continue;
        const double d = delta[f], w = 1.0 + d * d, sigma = s.finals[f].sign;
        const double a = (c == 0 ? sigma / w : -sigma * d * std::fabs(d) / w) + 0.0;  // Gamma_k / (Gamma_gamma b); +0, not -0
        width[k] = total * bEff[f] * a;
        if (!dGdq) continue;
        (*dGdq)[k * m + totalSlot] += bEff[f] * a;
        for (std::size_t i = 0; i < m; i++) (*dGdq)[k * m + i] += total * a * db[f * m + i];
        if (deltaSlot[f] >= 0) {
          const double da = c == 0 ? -sigma * 2.0 * d / (w * w) : -sigma * 2.0 * std::fabs(d) / (w * w);
          (*dGdq)[k * m + deltaSlot[f]] += total * bEff[f] * da;
        }
      }
  }
  return width;
}

namespace {

/// The photon amplitudes (MeV^1/2, slot order) of record level L at energy E (MeV) with the record's values q and
/// the full vector p's particle amplitudes: physics/ amplitudes_from_observed with the particle columns as amplitudes
/// (Brune: n_lambda from them; Park: 0, the photon amplitude is (Gamma / 2P)^1/2), the photon columns as widths.
std::vector<double> PhotonAmplitudes(const CompoundSystem &compound, const SystemOfCompound &structure, bool park, bool gsl, bool rmc,
                                     const GammaDecayLevel &L, double E, const std::vector<double> &q,
                                     const std::vector<double> &particle) {
  const std::vector<double> widths = GammaDecaySet::PartialWidths(L, q);
  const std::size_t n = L.channelOfBlockChannel.size();
  phys::ProjectLevels in;
  in.form = phys::ProjectLevels::Form::Observed;
  in.E = {E};
  in.value.assign(n, 0.0);
  in.is_amplitude.assign(n, 1);
  std::vector<int> slotOfColumn(n, -1);
  for (std::size_t k = 0; k < n; k++) {
    const int ch = L.channelOfBlockChannel[k];
    for (std::size_t s = 0; s < L.photonChannels.size(); s++)
      if (L.photonChannels[s] == ch) slotOfColumn[k] = (int)s;
    if (slotOfColumn[k] >= 0) {
      in.value[k] = widths[slotOfColumn[k]];
      in.is_amplitude[k] = 0;
    } else {
      in.value[k] = park ? 0.0 : particle[ch - 1];
    }
  }
  const phys::ChannelFunctions cf = EngineChannelFunctions(compound, structure.map, gsl, rmc);
  const phys::BasisAmplitudes out =
      phys::amplitudes_from_observed(structure.system, L.block, phys::internal_channels(structure.system, L.block), in,
                                     park ? phys::LevelBasis::Park : phys::LevelBasis::Brune, cf);
  std::vector<double> amplitude(L.photonChannels.size(), 0.0);
  for (std::size_t k = 0; k < n; k++)
    if (slotOfColumn[k] >= 0) amplitude[slotOfColumn[k]] = out.gamma[k];
  return amplitude;
}

/// N_f^(-1/2) of final state (j, la) at energy E (MeV) and amplitudes gamma (per channel of its group), with
/// dN_f/dgamma_d per channel: n0 = 1 (Brune) or Park's J_f = 1 - sum_c gamma_c^2 S'_c(E) over its particle channels,
/// as CompoundSystem::CalcShiftFunctions forms it (section 2.7).
double FinalSqrtNF(const CompoundSystem &compound, bool park, bool gsl, int j, int la, double E, const std::vector<double> &gamma,
                   std::vector<double> *dN) {
  const JGroup *group = compound.GetJGroup(j);
  const ALevel *level = group->GetLevel(la);
  double n0 = 1.0;
  if (dN) dN->assign(group->NumChannels(), 0.0);
  for (int ch = 1; ch <= group->NumChannels(); ch++) {
    const AChannel *channel = group->GetChannel(ch);
    const PPair *pair = compound.GetPair(channel->GetPairNum());
    if (dN && ch <= level->NumNFIntegrals())
      (*dN)[ch - 1] += 2.0 * gamma[ch - 1] * phys::bound_state_coefficient(channel->GetRadius(), pair->GetRedMass()) *
                       level->GetNFIntegral(ch);
    if (!park || channel->GetRadType() != 'P' || pair->GetPType() != 0) continue;
    ChannelFunc cf(pair, gsl, channel->GetRadius());
    const double dS = cf.ShiftDerivative(channel->GetL(), pair->SpinOrbitTwoJ(group->GetJ()),
                                         E - pair->GetSepE() - pair->GetExE());
    n0 -= std::pow(gamma[ch - 1], 2.0) * dS;
    if (dN) (*dN)[ch - 1] -= 2.0 * gamma[ch - 1] * dS;
  }
  return SqrtNFFactor(compound, group, level, n0, gamma);
}

}  // namespace

/// The map of one record at a point (GammaDecaySet::Map): the fill's parts and the engine's amplitudes.
struct GammaDecaySet::Point {
  GammaDecayFill fill;              ///< without q
  std::vector<double> internal;     ///< the engine's photon amplitudes (MeV^1/2), slot order
};

GammaDecaySet::Point GammaDecaySet::Map(const CompoundSystem &compound, const GammaDecayLevel &L, double E,
                                        const std::vector<double> &q, const std::vector<double> &particle,
                                        const std::vector<std::pair<double, std::vector<double>>> &finals,
                                        bool slopes) const {
  Point point;
  GammaDecayFill &fill = point.fill;
  fill.total = PhotonAmplitudes(compound, *structure_, park_, gslCoulomb_, rmc_, L, E, q, particle);
  point.internal = fill.total;
  const std::size_t m = L.photonChannels.size();
  fill.external.assign(m, std::complex<double>(0.0, 0.0));
  fill.floor.assign(m, 0);
  fill.dLevel.assign(m, {});
  fill.dFinal.assign(m, {});
  if (!L.HasExternal()) return point;
  // gamma^ext at the current parameters, subtracted as TransformIn does (physics/ subtract_external_capture) in the
  // level's row of its block.
  const JGroup *group = compound.GetJGroup(L.spec.group);
  const std::size_t n = L.channelOfBlockChannel.size();
  std::vector<double> row(n, 0.0);
  std::vector<phys::ExternalWidth> ext(n);
  std::vector<int> slotOfColumn(n, -1);
  for (std::size_t k = 0; k < n; k++)
    for (std::size_t s = 0; s < m; s++)
      if (L.photonChannels[s] == L.channelOfBlockChannel[k]) slotOfColumn[k] = (int)s;
  for (std::size_t k = 0; k < n; k++) {
    const int s = slotOfColumn[k];
    if (s < 0) continue;
    row[k] = fill.total[s];
    const GammaDecayLevel::External &e = L.external[s];
    if (!e.group) continue;
    ExternalCaptureInputs in;
    in.levelEnergy = E;
    in.levelGamma = particle;
    in.finalEnergy = finals[s].first;
    in.finalGamma = finals[s].second;
    std::vector<double> dN;
    in.finalSqrtNF = FinalSqrtNF(compound, park_, gslCoulomb_, e.group, e.level, in.finalEnergy, in.finalGamma,
                                 slopes ? &dN : nullptr);
    ExternalCaptureSlopes d;
    const AChannel *channel = group->GetChannel(L.photonChannels[s]);
    fill.external[s] = ExternalReducedWidthAt(compound, group, channel, {e.group, e.level}, in, *configure_,
                                              slopes ? &d : nullptr);
    ext[k].present = true;
    ext[k].value = fill.external[s];
    if (!slopes) continue;
    // d gamma^ext / d gamma~_fd = sum_c d gamma^ext / d s_fc ds_fc / d gamma~_fd, s_fc = gamma~_fc N_f^-1/2 K_c:
    // ds_fc / d gamma~_fd = delta_cd N_f^-1/2 K_c - s_fc N_f^-1 dN_f/d gamma~_fd / 2 (Eq. 2.8).
    const ALevel *finalLevel = compound.GetJGroup(e.group)->GetLevel(e.level);
    const double r = in.finalSqrtNF;
    fill.dLevel[s] = d.level;
    fill.dFinal[s].assign(in.finalGamma.size(), std::complex<double>(0.0, 0.0));
    for (std::size_t c = 0; c < in.finalGamma.size(); c++) {
      if (d.finalScale[c] == std::complex<double>(0.0, 0.0)) continue;
      const double K = finalLevel->GetECConversionFactor((int)c + 1);
      const double sc = in.finalGamma[c] * r * K;
      for (std::size_t dd = 0; dd < in.finalGamma.size(); dd++) {
        const double ds = (c == dd ? r * K : 0.0) - 0.5 * sc * r * r * dN[dd];
        fill.dFinal[s][dd] += d.finalScale[c] * ds;
      }
    }
  }
  std::vector<std::string> warnings;
  std::vector<std::size_t> exceeded;
  phys::subtract_external_capture(structure_->system, L.block, row, ext, warnings, &exceeded);
  for (std::size_t i : exceeded)
    if (slotOfColumn[i] >= 0) fill.floor[slotOfColumn[i]] = 1;
  for (std::size_t k = 0; k < n; k++)
    if (slotOfColumn[k] >= 0) point.internal[slotOfColumn[k]] = row[k];
  return point;
}

vector_r GammaDecaySet::ToEngine(const CompoundSystem &compound, LevelState &state, const vector_r &values) const {
  vector_r p = values;
  const std::vector<std::vector<double>> numbers = NumbersAt(values);
  for (std::size_t r = 0; r < levels_.size(); r++) {
    const GammaDecayLevel &L = levels_[r];
    const int nc = compound.GetJGroup(L.spec.group)->NumChannels();
    std::vector<double> particle(nc);
    const std::vector<double> &q = numbers[r];
    for (int ch = 1; ch <= nc; ch++) particle[ch - 1] = p[L.energyIndex + ch];
    // The final states of the external parts, as the vector has them (their particle amplitudes are no record slot).
    std::vector<std::pair<double, std::vector<double>>> finals(L.photonChannels.size());
    for (std::size_t s = 0; s < L.photonChannels.size(); s++) {
      const GammaDecayLevel::External &e = L.external[s];
      if (!e.group) continue;
      finals[s].first = values[e.energyIndex];
      for (int ch = 1; ch <= compound.GetJGroup(e.group)->NumChannels(); ch++)
        finals[s].second.push_back(values[e.energyIndex + ch]);
    }
    Point point = Map(compound, L, p[L.energyIndex], q, particle, finals, L.HasExternal());
    for (std::size_t s = 0; s < L.photonChannels.size(); s++) p[L.energyIndex + L.photonChannels[s]] = point.internal[s];
    bool floor = false;
    for (unsigned char f : point.fill.floor) floor = floor || f;
    if (floor && Config::FirstTimeIn(configure_->printedOnce, "gamma_decay: below the external-capture floor"))
      std::cerr << "WARNING: <gammaDecays> record '" << L.spec.name << "': a partial width is below 2P (Im gamma^ext)^2 / n,"
                << " which no internal amplitude reaches; the observed width is that floor (gamma_int = -Re gamma^ext, as"
                << " TransformIn; gamma_parametrisation.md section 2.7; further such warnings are not shown)" << std::endl;
    point.fill.q = q;
    state.SetGammaDecayFill((int)r, std::move(point.fill));
  }
  return p;
}

void GammaDecaySet::PullBack(const CompoundSystem &compound, const LevelState &state, const GradAccum &accum, const ParamIndexMap &map,
                             vector_r &gradFull) const {
  const int nFull = (int)gradFull.size();
  auto add = [&](int index, double v) {
    if (index >= 0 && index < nFull) gradFull[index] += v;
  };
  for (std::size_t r = 0; r < levels_.size(); r++) {
    const GammaDecayLevel &L = levels_[r];
    const int j = L.spec.group, la = L.spec.level;
    const JGroup *group = compound.GetJGroup(j);
    const ALevel *level = group->GetLevel(la);
    if (!state.HasGammaDecayFill((int)r)) continue;
    const GammaDecayFill &fill = state.GammaDecayFillOf((int)r);
    const std::vector<double> &q = fill.q;
    const std::vector<double> &amp = fill.total;  // g_k, the amplitude of the record's observed width
    const std::size_t m = L.photonChannels.size();
    std::vector<double> gbar(m), particle(group->NumChannels());
    // d gamma_int / d g = R and d gamma_int / d Im gamma^ext = B (Eq. 2.8): 1 and 0 without external part, 0 and 0
    // on the floor (gamma_int = -Re gamma^ext there).
    std::vector<double> R(m, 1.0), B(m, 0.0);
    for (std::size_t k = 0; k < m; k++) {
      gbar[k] = accum.gamma[j - 1][la - 1][L.photonChannels[k] - 1];
      if (!L.external[k].group) continue;
      const double g = amp[k], b = fill.external[k].imag(), root = std::sqrt(std::max(g * g - b * b, 0.0));
      R[k] = (fill.floor[k] || !(root > 0.0)) ? 0.0 : std::fabs(g) / root;
      B[k] = (fill.floor[k] || !(root > 0.0)) ? 0.0 : (g < 0.0 ? 1.0 : -1.0) * b / root;
    }
    bool any = false;
    for (double g : gbar) any = any || g != 0.0;
    if (!any) continue;
    for (int ch = 1; ch <= group->NumChannels(); ch++) particle[ch - 1] = state.Gamma(j, la, ch);
    // Gamma_gamma and the branchings: dg_k/dq_i = g_k / (2 Gamma_k) dGamma_k/dq_i (0 where Gamma_k = 0: g ~ b^1/2
    // has no finite derivative at b = 0), times R_k.
    std::vector<double> dGdq;
    const std::vector<double> widths = PartialWidths(L, q, &dGdq);
    auto slotOf = [&](int ch) {
      for (std::size_t k = 0; k < m; k++)
        if (L.photonChannels[k] == ch) return k;
      return m;
    };
    for (std::size_t i = 0; i < m; i++) {
      double v = 0.0;
      if (L.slots[i].part == GammaDecayLevel::Part::Mixing) {
        // A mixing ratio, in the amplitudes' own form, regular through delta = 0 (the width form is 0/0 there):
        // g_L ~ (1 + d^2)^-1/2, so dg_L/dd = -g_L d / (1 + d^2); g_L+1 = K d (1 + d^2)^-1/2, so dg_L+1/dd =
        // K (1 + d^2)^-3/2, with K = 2^1/2 g_L+1 at d = 1 (the map at d = 1, the other numbers as they are).
        const ex::GammaFinalState &f = L.spec.finals[L.slots[i].final];
        const std::size_t kL = slotOf(f.channels[0]), kL1 = slotOf(f.channels[1]);
        const double d = q[i], w = 1.0 + d * d;
        std::vector<double> unit = q;
        unit[i] = 1.0;
        const double K = std::sqrt(2.0) * PhotonAmplitudes(compound, *structure_, park_, gslCoulomb_, rmc_, L,
                                                           state.E(j, la), unit, particle)[kL1];
        v = -gbar[kL] * R[kL] * amp[kL] * d / w + gbar[kL1] * R[kL1] * K / (w * std::sqrt(w));
      } else {
        for (std::size_t k = 0; k < m; k++)
          if (widths[k] != 0.0) v += gbar[k] * R[k] * amp[k] / (2.0 * widths[k]) * dGdq[k * m + i];
      }
      add(map.GammaIndexFlat(j, la, L.photonChannels[i]), v);
    }
    // The particle amplitudes (Brune): g_k^2 = |Gamma_k| n / (2 P_k), n = 1 + sum_c gamma_c^2 S'_c.
    const double E = state.E(j, la);
    if (!park_) {
      const phys::ChannelFunctions cf = EngineChannelFunctions(compound, structure_->map, gslCoulomb_, rmc_);
      double n = 1.0;
      std::vector<std::pair<int, double>> terms;  // (channel, gamma_c S'_c)
      for (std::size_t k = 0; k < L.channelOfBlockChannel.size(); k++) {
        const int ch = L.channelOfBlockChannel[k];
        if (!IsParticleChannel(compound, group, ch)) continue;
        const double dS = cf(L.block, k, E).dS_dE;
        n += particle[ch - 1] * particle[ch - 1] * dS;
        terms.push_back({ch, particle[ch - 1] * dS});
      }
      double sumGP = 0.0;
      for (std::size_t k = 0; k < m; k++) sumGP += gbar[k] * R[k] * amp[k];
      for (const auto &t : terms) add(map.GammaIndexFlat(j, la, t.first), sumGP * t.second / n);
    }
    // The external parts (Eq. 2.8): -Re d gamma^ext + B Im d gamma^ext, in the level's particle amplitudes and in
    // the final states' amplitudes.
    for (std::size_t k = 0; k < m; k++) {
      const GammaDecayLevel::External &e = L.external[k];
      if (!e.group || gbar[k] == 0.0) continue;
      for (std::size_t c = 0; c < fill.dLevel[k].size(); c++) {
        const std::complex<double> d = fill.dLevel[k][c];
        if (d != std::complex<double>(0.0, 0.0)) add(map.GammaIndexFlat(j, la, (int)c + 1), gbar[k] * (-d.real() + B[k] * d.imag()));
      }
      for (std::size_t c = 0; c < fill.dFinal[k].size(); c++) {
        const std::complex<double> d = fill.dFinal[k][c];
        if (d != std::complex<double>(0.0, 0.0))
          add(map.GammaIndexFlat(e.group, e.level, (int)c + 1), gbar[k] * (-d.real() + B[k] * d.imag()));
      }
    }
    // The energies, by central differences of the whole map (Eq. 2.5, gamma^ext included), each when it is free: the
    // level's, and that of each final state of an external part (it is free when the state lies above a particle
    // threshold of its nucleus, CompoundSystem::FillMnParams; section 2.7).
    auto isFree = [&](int index) { return index >= 0 && index < map.NumFull() && map.FullToPacked(index) >= 0; };
    std::vector<std::pair<double, std::vector<double>>> finals(m);
    for (std::size_t s = 0; s < m; s++) {
      const GammaDecayLevel::External &e = L.external[s];
      if (!e.group) continue;
      finals[s].first = state.E(e.group, e.level);
      for (int ch = 1; ch <= compound.GetJGroup(e.group)->NumChannels(); ch++)
        finals[s].second.push_back(state.Gamma(e.group, e.level, ch));
    }
    auto difference = [&](double upE, double downE, const std::vector<std::pair<double, std::vector<double>>> &upF,
                          const std::vector<std::pair<double, std::vector<double>>> &downF) {
      const std::vector<double> up = Map(compound, L, upE, q, particle, upF, false).internal;
      const std::vector<double> down = Map(compound, L, downE, q, particle, downF, false).internal;
      double d = 0.0;
      for (std::size_t k = 0; k < m; k++) d += gbar[k] * (up[k] - down[k]) / (2.0 * kEnergyStep);
      return d;
    };
    const int ei = map.EIndexFlat(j, la);
    if (isFree(ei)) add(ei, difference(E + kEnergyStep, E - kEnergyStep, finals, finals));
    std::vector<int> done;
    for (std::size_t s = 0; s < m; s++) {
      const GammaDecayLevel::External &e = L.external[s];
      if (!e.group) continue;
      const int fi = map.EIndexFlat(e.group, e.level);
      if (!isFree(fi) || std::find(done.begin(), done.end(), fi) != done.end()) continue;
      done.push_back(fi);
      std::vector<std::pair<double, std::vector<double>>> upF = finals, downF = finals;
      for (std::size_t t = 0; t < m; t++)
        if (L.external[t].group == e.group && L.external[t].level == e.level) {
          upF[t].first += kEnergyStep;
          downF[t].first -= kEnergyStep;
        }
      add(fi, difference(E, E, upF, downF));
    }
  }
}

std::vector<int> GammaDecaySet::InputIndices(const ParamIndexMap &map, int r) const {
  const GammaDecayLevel &L = levels_.at(r);
  std::vector<int> indices{map.EIndexFlat(L.spec.group, L.spec.level)};
  const int nc = (int)L.channelOfBlockChannel.size();
  for (int ch = 1; ch <= nc; ch++) indices.push_back(map.GammaIndexFlat(L.spec.group, L.spec.level, ch));
  for (const GammaDecayLevel::External &e : L.external) {
    if (!e.group) continue;
    // its energy (free when the state lies above a particle threshold of its nucleus) and its channels: the
    // amplitudes N_f and s_fc read
    indices.push_back(map.EIndexFlat(e.group, e.level));
    for (int ch = 1; ch <= e.channels; ch++) indices.push_back(map.GammaIndexFlat(e.group, e.level, ch));
  }
  return indices;
}

std::vector<GammaDecaySet::PenaltyRow> GammaDecaySet::PenaltyRows(const CompoundSystem &compound, const LevelState &state) const {
  std::vector<PenaltyRow> rows;
  for (std::size_t r = 0; r < levels_.size(); r++) {
    const GammaDecayLevel &L = levels_[r];
    const ex::GammaDecaySpec &s = L.spec;
    if (!state.HasGammaDecayFill((int)r)) continue;
    const std::vector<double> &q = state.GammaDecayFillOf((int)r).q;
    auto index = [&](std::size_t i) { return L.energyIndex + L.photonChannels[i]; };
    double sum = 0.0;
    std::vector<std::pair<int, double>> sumTerms;
    for (std::size_t i = 0; i < L.slots.size(); i++) {
      const GammaDecayLevel::Slot &slot = L.slots[i];
      const ex::GammaNumber &number = slot.part == GammaDecayLevel::Part::TotalWidth ? s.totalWidth
                                      : slot.part == GammaDecayLevel::Part::Branching ? s.finals[slot.final].branching
                                                                                      : s.finals[slot.final].mixing;
      const double scale = slot.part == GammaDecayLevel::Part::TotalWidth ? kEVPerMeV : 1.0;
      if (slot.part == GammaDecayLevel::Part::Branching) {
        sum += q[i];
        sumTerms.push_back({index(i), 1.0});
      }
      if (!number.hasPrior) continue;
      const double w = number.priorWidth * scale;
      rows.push_back({(q[i] - number.priorCentre * scale) / w, {{index(i), 1.0 / w}}});
    }
    if (s.hasRemainderPrior && s.remainder >= 0) {
      PenaltyRow row;
      row.r = (1.0 - sum - s.remainderPriorCentre) / s.remainderPriorWidth;
      for (const auto &t : sumTerms) row.d.push_back({t.first, -1.0 / s.remainderPriorWidth});
      rows.push_back(row);
    }
    // The wall on a negative remainder: always one row, 0 inside, so that the number of rows does not move.
    PenaltyRow wall;
    if (s.remainder >= 0 && sum > 1.0) {
      wall.r = (sum - 1.0) / kWallWidth;
      for (const auto &t : sumTerms) wall.d.push_back({t.first, 1.0 / kWallWidth});
    }
    rows.push_back(wall);
  }
  return rows;
}

double GammaDecaySet::Penalty(const CompoundSystem &compound, const LevelState &state) const {
  double chi2 = 0.0;
  for (const PenaltyRow &row : PenaltyRows(compound, state)) chi2 += row.r * row.r;
  return chi2;
}

void GammaDecaySet::AddPenaltyGradient(const CompoundSystem &compound, const LevelState &state, std::vector<std::pair<int, double>> &grad) const {
  for (const PenaltyRow &row : PenaltyRows(compound, state))
    for (const auto &d : row.d) grad.push_back({d.first, 2.0 * row.r * d.second});
}

std::string GammaDecaySet::Report(const CompoundSystem &compound, const LevelState &state, const ObservedLevelTable &table) const {
  std::ostringstream out;
  if (levels_.empty()) return "";
  out << "RADIATIVE DECAYS (<gammaDecays>: Gamma_gamma, branchings, mixing ratios of Krane and Steffen; Gamma_gamma is "
         "not part of the resonance's total width)" << std::endl
      << std::endl;
  for (std::size_t r = 0; r < levels_.size(); r++) {
    const GammaDecayLevel &L = levels_[r];
    const ex::GammaDecaySpec &s = L.spec;
    if (!state.HasGammaDecayFill((int)r)) continue;
    const std::vector<double> &q = state.GammaDecayFillOf((int)r).q;
    double excess = 0.0;
    const std::vector<double> widths = PartialWidths(L, q, nullptr, &excess);
    out << s.name << ": level " << s.levelId << ", E_x = " << Fmt(s.levelEnergy) << " MeV" << std::endl;
    out << "  Gamma_gamma = " << Fmt(q[0]) << " eV" << std::endl;
    for (std::size_t f = 0; f < s.finals.size(); f++) {
      const ex::GammaFinalState &final = s.finals[f];
      // The branching the widths have (the remainder 1 - sum of the others; with a sum above 1, as scaled).
      double b = 0.0, delta = 0.0;
      for (std::size_t i = 0; i < L.slots.size(); i++)
        if (L.slots[i].final == (int)f && L.slots[i].part == GammaDecayLevel::Part::Mixing) delta = q[i];
      for (std::size_t k = 0; k < L.photonChannels.size(); k++)
        for (const int ch : final.channels)
          if (L.photonChannels[k] == ch && q[0] > 0.0) b += std::fabs(widths[k]) / q[0];
      out << "  " << final.label << ": b = " << Fmt(b) << ((int)f == s.remainder ? " (remainder)" : "");
      if (final.L.size() == 2) out << ", delta = " << Fmt(delta);
      out << "; Gamma =";
      for (std::size_t c = 0; c < final.channels.size(); c++)
        for (std::size_t k = 0; k < L.photonChannels.size(); k++)
          if (L.photonChannels[k] == final.channels[c]) out << " " << Fmt(widths[k]) << " eV (L = " << final.L[c] << ")";
      out << std::endl;
      // A channel on the external-capture floor: the engine's observed width (TransformOut) is not the record's.
      const GammaDecayFill &fill = state.GammaDecayFillOf((int)r);
      for (std::size_t k = 0; k < L.photonChannels.size(); k++)
        for (const int ch : final.channels)
          if (L.photonChannels[k] == ch && k < fill.floor.size() && fill.floor[k])
            out << "  note: the width " << Fmt(widths[k]) << " eV to " << final.label << " is below the floor of external "
                << "capture, 2P (Im gamma^ext)^2 / n; the observed width is "
                << Fmt(table.Level(s.group, s.level).value[ch - 1]) << " eV (section 2.7)"
                << std::endl;
      // MINUIT reaches a limit through its sine transformation, to about 1e-7 here (gamma_parametrisation.md 2.6).
      if (b < kAtBound && (final.role == ex::GammaFinalState::Role::Listed || (int)f == s.remainder))
        out << "  note: the branching to " << final.label << " is at its bound 0 (b < " << kAtBound << ")" << std::endl;
    }
    if (s.remainder >= 0 && excess > 0.0)
      out << "  note: the explicit branchings sum to " << Fmt(1.0 + excess)
          << " > 1: the remainder is 0 and they are scaled to sum to 1 (the wall of the objective pulls them back)" << std::endl;
    out << std::endl;
  }
  return out.str();
}
