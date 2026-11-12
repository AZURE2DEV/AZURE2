// A level's radiative decay as Gamma_gamma, b_f, delta_f between the parameter vector and the engine's amplitudes
// (GammaDecayMap.h; docs/upgrade/gamma_parametrisation.md sections 2.2-2.6).
#include "GammaDecayMap.h"

#include <cmath>
#include <cstdio>
#include <sstream>

#include "AChannel.h"
#include "ALevel.h"
#include "AZUREGrad.h"
#include "AZUREParams.h"
#include "CNuc.h"
#include "CNucLevelParameters.h"
#include "Config.h"
#include "JGroup.h"
#include "PPair.h"
#include "azure2/experiment/LevelTables.h"
#include "azure2/experiment/ProjectReader.h"
#include "azure2/experiment/Schema.h"
#include "azure2/io/LevelRows.h"
#include "azure2/io/RecordReader.h"
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

void Print(const Config &configure, const ex::Diagnostic &d) {
  configure.outStream << (d.severity == ex::Severity::Warning ? "WARNING: " : "ERROR: ") << "<gammaDecays> line " << d.line
                      << " [" << d.code << "]: " << d.message << std::endl;
}

bool IsPhotonChannel(CNuc &compound, JGroup *group, int ch) {
  return compound.GetPair(group->GetChannel(ch)->GetPairNum())->GetPType() == 10;
}

bool IsParticleChannel(CNuc &compound, JGroup *group, int ch) {
  return compound.GetPair(group->GetChannel(ch)->GetPairNum())->GetPType() == 0;
}

std::string Fmt(double x) {
  char buf[32];
  std::snprintf(buf, sizeof buf, "%.6g", x);
  return buf;
}

}  // namespace

std::shared_ptr<const GammaDecaySet> GammaDecaySet::Read(CNuc &compound, const Config &configure,
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
  try {
    set->structure_ = std::make_shared<const CNucSystem>(BuildSystemFromCNuc(compound, configure, CNucSystemContent::Structure));
  } catch (const std::exception &e) {
    configure.outStream << "ERROR: <gammaDecays>: " << e.what() << std::endl;
    return nullptr;
  }
  const CNucSystemMap &map = set->structure_->map;
  for (const ex::GammaDecaySpec &s : built.specs) {
    if (formal)
      return refuse(s, "the formal basis (--no-brune) mixes the photon amplitudes of all levels of a J^pi group, so a "
                       "level's total radiative width and branchings cannot be its own parameters; use the Brune or "
                       "Park basis (gamma_parametrisation.md section 2.4)");
    if (!transform)
      return refuse(s, "--no-transform reads the level parameters as amplitudes; the record gives observed widths");
    GammaDecayLevel L;
    L.spec = s;
    JGroup *group = compound.GetJGroup(s.group);
    ALevel *level = group->GetLevel(s.level);
    if (!group->IsInRMatrix() || !level->IsInRMatrix()) return refuse(s, "the level is not in the R-matrix");
    for (std::size_t b = 0; b < map.groupOfBlock.size(); b++)
      if (map.groupOfBlock[b] == s.group) L.block = b;
    L.channelOfBlockChannel = map.channelOfBlockChannel.at(L.block);
    if ((int)L.channelOfBlockChannel.size() != group->NumChannels())
      return refuse(s, "the level's J^pi group has channels outside the R-matrix");
    for (int ch = 1; ch <= group->NumChannels(); ch++)
      if (IsPhotonChannel(compound, group, ch)) L.photonChannels.push_back(ch);
    // The vector's layout (CNuc::FillMnParams): per level its energy, then one entry per channel of its group.
    int index = 0, number = 1;
    for (int j = 1; j < s.group; j++) {
      index += compound.GetJGroup(j)->NumLevels() * (1 + compound.GetJGroup(j)->NumChannels());
      number += compound.GetJGroup(j)->NumLevels();
    }
    L.energyIndex = index + (s.level - 1) * (1 + group->NumChannels());
    L.levelNumber = number + s.level - 1;
    for (const int ch : L.photonChannels) {
      AChannel *c = group->GetChannel(ch);
      PPair *pair = compound.GetPair(c->GetPairNum());
      if (std::fabs(level->GetE() - pair->GetHeavyExE()) < 1.e-3 && group->GetJ() == pair->GetJ(2) && group->GetPi() == pair->GetPi(2))
        return refuse(s, "a photon channel to the level's own state is a static moment, not a decay");
      // External capture into the final state by this multipolarity (CNucExternalCapture.cpp, ExternalReducedWidth):
      // its observed width would need gamma^ext at every evaluation (section 2.4).
      if (!external) continue;
      const char type = c->GetRadType();
      const int Lm = c->GetL();
      for (int j = 1; j <= compound.NumJGroups(); j++)
        for (int la = 1; la <= compound.GetJGroup(j)->NumLevels(); la++) {
          const ALevel *f = compound.GetJGroup(j)->GetLevel(la);
          if (!f->IsECLevel() || f->GetECPairNum() != c->GetPairNum()) continue;
          const unsigned char mask = f->GetECMultMask();
          if ((type == 'E' && Lm == 1 && (mask & isE1)) || (type == 'M' && Lm == 1 && (mask & isM1)) ||
              (type == 'E' && Lm == 2 && (mask & isE2)))
            return refuse(s, "a final state (pair " + std::to_string(pair->GetPairKey()) + ") has external capture by this multipolarity in "
                             "this run; write the level's partial widths in <levels> instead (gamma_parametrisation.md "
                             "section 2.4)");
        }
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
    // The partial widths of the record as the level's input widths (eV): TransformIn converts them.
    std::vector<double> q;
    for (const GammaDecayLevel::Slot &slot : L.slots) {
      switch (slot.part) {
        case GammaDecayLevel::Part::TotalWidth: q.push_back(s.totalWidth.value * kEVPerMeV); break;
        case GammaDecayLevel::Part::Branching: q.push_back(s.finals[slot.final].branching.value); break;
        case GammaDecayLevel::Part::Mixing: q.push_back(s.finals[slot.final].mixing.value); break;
      }
    }
    const std::vector<double> widths = PartialWidths(L, q);
    for (std::size_t k = 0; k < L.photonChannels.size(); k++) level->SetGamma(L.photonChannels[k], widths[k]);
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

bool GammaDecaySet::AddParameter(AZUREParams &params, int j, int la, int ch, int levelNumber) const {
  const int r = RecordOf(j, la);
  if (SlotOf(r, ch) < 0) return false;
  const GammaDecayLevel &L = levels_[r];
  const GammaDecayLevel::Slot &slot = L.slots[SlotOf(r, ch)];
  const ex::GammaDecaySpec &s = L.spec;
  ROOT::Minuit2::MnUserParameters &p = params.GetMinuitParams();
  std::string name = s.name, part, final;
  ex::GammaNumber number;
  double scale = 1.0;
  switch (slot.part) {
    case GammaDecayLevel::Part::TotalWidth:
      part = "total_width";
      number = s.totalWidth;
      scale = kEVPerMeV;
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
  const double value = number.value * scale;
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

double GammaDecaySet::ValueOr(const CNuc &compound, int j, int la, int ch, double otherwise) const {
  const int r = RecordOf(j, la), slot = SlotOf(r, ch);
  if (slot < 0 || !compound.HasGammaDecayValues(r)) return otherwise;
  return compound.GammaDecayValues(r).at(slot);
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
std::vector<double> PhotonAmplitudes(CNuc &compound, const CNucSystem &structure, bool park, bool gsl, bool rmc,
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
  const phys::ChannelFunctions cf = ChannelFunctionsFromCNuc(compound, structure.map, gsl, rmc);
  const phys::BasisAmplitudes out =
      phys::amplitudes_from_observed(structure.system, L.block, phys::internal_channels(structure.system, L.block), in,
                                     park ? phys::LevelBasis::Park : phys::LevelBasis::Brune, cf);
  std::vector<double> amplitude(L.photonChannels.size(), 0.0);
  for (std::size_t k = 0; k < n; k++)
    if (slotOfColumn[k] >= 0) amplitude[slotOfColumn[k]] = out.gamma[k];
  return amplitude;
}

}  // namespace

vector_r GammaDecaySet::ToEngine(CNuc &compound, const vector_r &values) const {
  vector_r p = values;
  for (std::size_t r = 0; r < levels_.size(); r++) {
    const GammaDecayLevel &L = levels_[r];
    const int nc = compound.GetJGroup(L.spec.group)->NumChannels();
    std::vector<double> q, particle(nc);
    for (const int ch : L.photonChannels) q.push_back(p[L.energyIndex + ch]);
    for (int ch = 1; ch <= nc; ch++) particle[ch - 1] = p[L.energyIndex + ch];
    const std::vector<double> amplitude =
        PhotonAmplitudes(compound, *structure_, park_, gslCoulomb_, rmc_, L, p[L.energyIndex], q, particle);
    for (std::size_t s = 0; s < L.photonChannels.size(); s++) p[L.energyIndex + L.photonChannels[s]] = amplitude[s];
    compound.SetGammaDecayValues((int)r, q);
  }
  return p;
}

void GammaDecaySet::PullBack(CNuc &compound, const GradAccum &accum, const ParamIndexMap &map, vector_r &gradFull) const {
  const int nFull = (int)gradFull.size();
  auto add = [&](int index, double v) {
    if (index >= 0 && index < nFull) gradFull[index] += v;
  };
  for (std::size_t r = 0; r < levels_.size(); r++) {
    const GammaDecayLevel &L = levels_[r];
    const int j = L.spec.group, la = L.spec.level;
    JGroup *group = compound.GetJGroup(j);
    ALevel *level = group->GetLevel(la);
    if (!compound.HasGammaDecayValues((int)r)) continue;
    const std::vector<double> &q = compound.GammaDecayValues((int)r);
    const std::size_t m = L.photonChannels.size();
    std::vector<double> gbar(m), amp(m), particle(group->NumChannels());
    for (std::size_t k = 0; k < m; k++) {
      gbar[k] = accum.gamma[j - 1][la - 1][L.photonChannels[k] - 1];
      amp[k] = level->GetFitGamma(L.photonChannels[k]);
    }
    bool any = false;
    for (double g : gbar) any = any || g != 0.0;
    if (!any) continue;
    for (int ch = 1; ch <= group->NumChannels(); ch++) particle[ch - 1] = level->GetFitGamma(ch);
    // The record's numbers: dp_k/dq_i = p_k / (2 Gamma_k) dGamma_k/dq_i (0 where Gamma_k = 0).
    std::vector<double> dGdq;
    const std::vector<double> widths = PartialWidths(L, q, &dGdq);
    for (std::size_t i = 0; i < m; i++) {
      double v = 0.0;
      for (std::size_t k = 0; k < m; k++)
        if (widths[k] != 0.0) v += gbar[k] * amp[k] / (2.0 * widths[k]) * dGdq[k * m + i];
      add(map.GammaIndexFlat(j, la, L.photonChannels[i]), v);
    }
    // The particle amplitudes (Brune): p_k^2 = |Gamma_k| n / (2 P_k), n = 1 + sum_c gamma_c^2 S'_c.
    const double E = level->GetFitE();
    if (!park_) {
      const phys::ChannelFunctions cf = ChannelFunctionsFromCNuc(compound, structure_->map, gslCoulomb_, rmc_);
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
      for (std::size_t k = 0; k < m; k++) sumGP += gbar[k] * amp[k];
      for (const auto &t : terms) add(map.GammaIndexFlat(j, la, t.first), sumGP * t.second / n);
    }
    // The level energy, by a central difference of the map (Eq. 2.5).
    const std::vector<double> up = PhotonAmplitudes(compound, *structure_, park_, gslCoulomb_, rmc_, L, E + kEnergyStep, q, particle);
    const std::vector<double> down = PhotonAmplitudes(compound, *structure_, park_, gslCoulomb_, rmc_, L, E - kEnergyStep, q, particle);
    double dE = 0.0;
    for (std::size_t k = 0; k < m; k++) dE += gbar[k] * (up[k] - down[k]) / (2.0 * kEnergyStep);
    add(map.EIndexFlat(j, la), dE);
  }
}

std::vector<GammaDecaySet::PenaltyRow> GammaDecaySet::PenaltyRows(const CNuc &compound) const {
  std::vector<PenaltyRow> rows;
  for (std::size_t r = 0; r < levels_.size(); r++) {
    const GammaDecayLevel &L = levels_[r];
    const ex::GammaDecaySpec &s = L.spec;
    if (!compound.HasGammaDecayValues((int)r)) continue;
    const std::vector<double> &q = compound.GammaDecayValues((int)r);
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

double GammaDecaySet::Penalty(const CNuc &compound) const {
  double chi2 = 0.0;
  for (const PenaltyRow &row : PenaltyRows(compound)) chi2 += row.r * row.r;
  return chi2;
}

void GammaDecaySet::AddPenaltyGradient(const CNuc &compound, std::vector<std::pair<int, double>> &grad) const {
  for (const PenaltyRow &row : PenaltyRows(compound))
    for (const auto &d : row.d) grad.push_back({d.first, 2.0 * row.r * d.second});
}

std::string GammaDecaySet::Report(const CNuc &compound) const {
  std::ostringstream out;
  if (levels_.empty()) return "";
  out << "RADIATIVE DECAYS (<gammaDecays>: Gamma_gamma, branchings, mixing ratios of Krane and Steffen; Gamma_gamma is "
         "not part of the resonance's total width)" << std::endl
      << std::endl;
  for (std::size_t r = 0; r < levels_.size(); r++) {
    const GammaDecayLevel &L = levels_[r];
    const ex::GammaDecaySpec &s = L.spec;
    if (!compound.HasGammaDecayValues((int)r)) continue;
    const std::vector<double> &q = compound.GammaDecayValues((int)r);
    double excess = 0.0;
    const std::vector<double> widths = PartialWidths(L, q, nullptr, &excess);
    out << s.name << ": level " << s.levelId << ", E_x = " << Fmt(s.levelEnergy) << " MeV" << std::endl;
    out << "  Gamma_gamma = " << Fmt(q[0]) << " eV" << std::endl;
    for (std::size_t f = 0; f < s.finals.size(); f++) {
      const ex::GammaFinalState &final = s.finals[f];
      double b = 0.0, delta = 0.0;
      for (std::size_t i = 0; i < L.slots.size(); i++) {
        if (L.slots[i].final != (int)f) continue;
        if (L.slots[i].part == GammaDecayLevel::Part::Branching) b = q[i];
        if (L.slots[i].part == GammaDecayLevel::Part::Mixing) delta = q[i];
      }
      if ((int)f == s.remainder) {
        b = 1.0;
        for (std::size_t i = 0; i < L.slots.size(); i++)
          if (L.slots[i].part == GammaDecayLevel::Part::Branching) b -= q[i];
      }
      out << "  " << final.label << ": b = " << Fmt(b) << ((int)f == s.remainder ? " (remainder)" : "");
      if (final.L.size() == 2) out << ", delta = " << Fmt(delta);
      out << "; Gamma =";
      for (std::size_t c = 0; c < final.channels.size(); c++)
        for (std::size_t k = 0; k < L.photonChannels.size(); k++)
          if (L.photonChannels[k] == final.channels[c]) out << " " << Fmt(widths[k]) << " eV (L = " << final.L[c] << ")";
      out << std::endl;
      if (b <= 1e-9 && (final.role == ex::GammaFinalState::Role::Listed || (int)f == s.remainder))
        out << "  note: the branching to " << final.label << " is at its bound 0" << std::endl;
    }
    if (s.remainder >= 0 && excess > 0.0)
      out << "  note: the explicit branchings sum to " << Fmt(1.0 + excess)
          << " > 1: the remainder is 0 and they are scaled to sum to 1 (the wall of the objective pulls them back)" << std::endl;
    out << std::endl;
  }
  return out.str();
}
