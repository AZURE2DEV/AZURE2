// The levels' input values (include/LevelInputs.h) and what an evaluation writes (include/LevelState.h).
//  (a) After CompoundSystem::Fill the input values are the <levels> rows' (ALevel::FileE, FileGamma, FileZeroWidth), but for a
//      gamma_decay record's photon channels, which hold the partial widths of the record's file numbers; the record
//      numbers are the file's; nothing is transformed yet.
//  (b) After CompoundSystem::Initialize the transformation has run (no error) and written the input values only: every level's
//      <levels> values are those of a compound only filled.  In the formal basis a level of several in a group moves
//      its input energy.
//  (c) A physical vector (FillCompoundFromParamsPhysical) is written into the input values at the layout's slots.
//  (d) What a fill computes beside the levels is the state's: a copy of the LevelState keeps the gamma_decay records'
//      fills of its vector when the compound is filled at another; the state has a capture-mask projector for every
//      slot of the set-up.
//  (e) The THM line shape's test for a ground-state moment reads the level's file energy, not one a vector moved.
// Arguments: .azr files, each run in the Brune and the formal basis.  Reads them in place and writes nothing.
#include <cmath>
#include <cstdio>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <Minuit2/MnUserParameters.h>

#include "AZUREParams.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "ECAmplitudeCache.h"
#include "GammaDecayMap.h"
#include "JGroup.h"
#include "ThmLineshape.h"

namespace {

int failures = 0;

int groundStateMoments = 0;  // the cases of (e)
void Check(bool ok, const std::string &what) {
  if (!ok) {
    std::printf("  FAIL  %s\n", what.c_str());
    failures++;
  }
}

struct Loaded {
  std::ostringstream sink;
  Config cfg{sink};
  CompoundSystem compound;
  LevelState state;  // the levels, as the run holds them
  LevelInputs inputs;
};

std::unique_ptr<Loaded> Load(const std::string &file, bool brune, bool initialize) {
  auto m = std::make_unique<Loaded>();
  m->cfg.configfile = file;
  if (m->cfg.ReadConfigFile() < 0) return nullptr;
  m->cfg.fileCheckMask = m->cfg.screenCheckMask = 0;
  if (!brune) m->cfg.paramMask &= ~(Config::USE_PARK_FORMALISM | Config::USE_BRUNE_FORMALISM);
  m->cfg.ResolveMatrixRoute();
  InitializeECAmplitudeCache(m->cfg);
  if (m->compound.Fill(m->cfg, m->state, m->inputs) == -1) return nullptr;
  if (initialize) m->compound.Initialize(m->cfg, m->state, m->inputs);
  return m;
}

void OnProject(const std::string &file, bool brune) {
  const std::string tag = file + (brune ? " (Brune)" : " (formal)");
  std::unique_ptr<Loaded> filled = Load(file, brune, false);
  std::unique_ptr<Loaded> initialized = Load(file, brune, true);
  if (!filled || !initialized) {
    // A gamma_decay record is refused in the formal basis: nothing to check there.
    Check(!brune, tag + ": cannot load");
    return;
  }
  const CompoundSystem &cn = filled->compound;
  const LevelInputs &in = filled->inputs;
  const GammaDecaySet &records = cn.GetGammaDecays();
  // (a)
  Check(!in.TransformedIn() && in.TransformError().empty(), tag + ": filled, not transformed");
  Check(in.GammaDecayNumbers() == records.FileNumbers(), tag + ": the records' numbers are the file's");
  for (int j = 1; j <= cn.NumJGroups(); j++)
    for (int la = 1; la <= cn.GetJGroup(j)->NumLevels(); la++) {
      const ALevel *level = cn.GetJGroup(j)->GetLevel(la);
      const std::string where = tag + ": level " + std::to_string(j) + "/" + std::to_string(la);
      Check(in.E(j, la) == level->FileE(), where + " energy");
      const int r = records.RecordOf(j, la);
      std::vector<double> widths;
      if (r >= 0) widths = GammaDecaySet::PartialWidths(records.Levels()[r], records.FileNumbers()[r]);
      for (int ch = 1; ch <= cn.GetJGroup(j)->NumChannels(); ch++) {
        const int slot = records.SlotOf(r, ch);
        const double expected = slot >= 0 ? widths[slot] : level->FileGamma(ch);
        Check(in.Gamma(j, la, ch) == expected, where + " channel " + std::to_string(ch) + " width");
        Check(filled->state.ZeroInputWidth(j, la, ch) == level->FileZeroWidth(ch), where + " channel " + std::to_string(ch) +
                                                                            " zero flag");
      }
    }
  // (b)
  const CompoundSystem &ci = initialized->compound;
  Check(initialized->inputs.TransformedIn() && initialized->inputs.TransformError().empty(), tag + ": transformed, no error");
  bool moved = false;
  for (int j = 1; j <= ci.NumJGroups(); j++)
    for (int la = 1; la <= ci.GetJGroup(j)->NumLevels(); la++) {
      const ALevel *a = ci.GetJGroup(j)->GetLevel(la), *b = cn.GetJGroup(j)->GetLevel(la);
      Check(a->FileE() == b->FileE(), tag + ": the level's file energy is not written by the transformation");
      for (int ch = 1; ch <= ci.GetJGroup(j)->NumChannels(); ch++)
        Check(a->FileGamma(ch) == b->FileGamma(ch), tag + ": the level's file width is not written");
      moved = moved || initialized->inputs.E(j, la) != a->FileE();
    }
  bool several = false;
  for (int j = 1; j <= ci.NumJGroups(); j++) several = several || ci.GetJGroup(j)->NumLevels() > 1;
  if (!brune && several) Check(moved, tag + ": a formal input energy differs from the file's");
  // (c)
  vector_r v(ci.ParameterLayout().Size());
  for (std::size_t i = 0; i < v.size(); i++) v[i] = 1.0 + 0.01 * (double)i;
  initialized->compound.FillCompoundFromParamsPhysical(v, initialized->state, initialized->inputs);
  const vector_r p = records.Empty() ? v : records.WithPartialWidths(v);
  for (const LevelParameterLayout::Level &slots : ci.ParameterLayout().Levels()) {
    Check(initialized->inputs.E(slots.group, slots.level) == p[slots.energy], tag + ": physical energy");
    for (int ch = 1; ch <= (int)slots.amplitude.size(); ch++)
      Check(initialized->inputs.Gamma(slots.group, slots.level, ch) == p[slots.amplitude[ch - 1]], tag + ": physical width");
  }
  if (!records.Empty()) Check(initialized->inputs.GammaDecayNumbers() == records.NumbersAt(v), tag + ": physical numbers");
  // (d) on a compound initialized from the file
  std::unique_ptr<Loaded> fresh = Load(file, brune, true);
  CompoundSystem &cf = fresh->compound;
  AZUREParams params;
  cf.FillMnParams(params, fresh->state, fresh->inputs, &fresh->cfg);
  const vector_r start = params.GetMinuitParams().Params();
  vector_r shifted = start;
  for (const GammaDecayLevel &L : records.Levels())
    for (const int ch : L.photonChannels) shifted[L.energyIndex + ch] *= 1.1;
  cf.FillCompoundFromParams(start, fresh->state);
  const LevelState atStart = fresh->state;
  cf.FillCompoundFromParams(shifted, fresh->state);
  for (std::size_t r = 0; r < records.Levels().size(); r++) {
    Check(atStart.HasGammaDecayFill((int)r) && fresh->state.HasGammaDecayFill((int)r), tag + ": record filled");
    if (!atStart.HasGammaDecayFill((int)r) || !fresh->state.HasGammaDecayFill((int)r)) continue;
    Check(atStart.GammaDecayFillOf((int)r).q == records.NumbersAt(start)[r], tag + ": the copy keeps its fill");
    Check(fresh->state.GammaDecayFillOf((int)r).q == records.NumbersAt(shifted)[r], tag + ": the state has the new fill");
  }
  for (int j = 1; j <= cf.NumJGroups(); j++)
    for (int ch = 1; ch <= cf.GetJGroup(j)->NumChannels(); ch++) {
      const int slot = cf.FormalCaptureMaskSlot(j, ch);
      if (slot < 0) continue;
      const std::size_t n = cf.FormalCaptureMaskOf(slot).levels.size();
      const std::size_t size = fresh->state.CaptureProjector(slot).size();
      Check(size == 0 || size == n * n, tag + ": a projector over the slot's levels");
    }
  // (e) The THM line shape tells a ground-state moment by the level's file energy: after a physical vector that moves
  // every level by 10 keV, a photon channel of a level that is its pair's final state still adds no width.
  vector_r shiftedE(ci.ParameterLayout().Size());
  for (std::size_t i = 0; i < shiftedE.size(); i++) shiftedE[i] = 1.0;
  for (const LevelParameterLayout::Level &slots : ci.ParameterLayout().Levels())
    shiftedE[slots.energy] = ci.GetJGroup(slots.group)->GetLevel(slots.level)->FileE() + 0.01;
  initialized->compound.FillCompoundFromParamsPhysical(shiftedE, initialized->state, initialized->inputs);
  for (int j = 1; j <= ci.NumJGroups(); j++)
    for (int la = 1; la <= ci.GetJGroup(j)->NumLevels(); la++)
      for (int ch = 1; ch <= ci.GetJGroup(j)->NumChannels(); ch++) {
        const JGroup *g = ci.GetJGroup(j);
        const AChannel *c = g->GetChannel(ch);
        const PPair *pair = ci.GetPair(c->GetPairNum());
        if ((c->GetRadType() != 'E' && c->GetRadType() != 'M') || g->GetJ() != pair->GetJ(2) ||
            g->GetPi() != pair->GetPi(2) || std::fabs(g->GetLevel(la)->FileE() - pair->GetExE()) >= 1.0e-3)
          continue;
        groundStateMoments++;
        std::vector<double> gammas(g->NumChannels(), 0.0);
        gammas[ch - 1] = 1.0;
        Check(ThmLevelWidth(&ci, j, la, g->GetLevel(la)->FileE() + 0.01, gammas,
                            initialized->cfg) == 0.0,
              tag + ": a ground-state moment adds no THM width when the level moves");
      }
}

}  // namespace

int main(int argc, char **argv) {
  for (int a = 1; a < argc; a++) {
    OnProject(argv[a], true);
    OnProject(argv[a], false);
  }
  std::printf("ground-state moments checked: %d\n", groundStateMoments);
  std::printf("%s: %d failure(s)\n", failures ? "FAIL" : "ok", failures);
  return failures ? 1 : 0;
}
