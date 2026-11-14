// The levels' input values (include/LevelInputs.h) and what an evaluation writes (include/LevelState.h).
//  (a) After CNuc::Fill the input values are the <levels> rows' (ALevel::FileE, FileGamma, FileZeroWidth), but for a
//      gamma_decay record's photon channels, which hold the partial widths of the record's file numbers; the record
//      numbers are the file's; nothing is transformed yet.
//  (b) After CNuc::Initialize the transformation has run (no error) and written the input values only: every level's
//      <levels> values are those of a compound only filled.  In the formal basis a level of several in a group moves
//      its input energy.
//  (c) A physical vector (FillCompoundFromParamsPhysical) is written into the input values at the layout's slots.
// Arguments: .azr files, each run in the Brune and the formal basis.  Reads them in place and writes nothing.
#include <cstdio>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "CNuc.h"
#include "Config.h"
#include "ECAmplitudeCache.h"
#include "GammaDecayMap.h"
#include "JGroup.h"

namespace {

int failures = 0;

void Check(bool ok, const std::string &what) {
  if (!ok) {
    std::printf("  FAIL  %s\n", what.c_str());
    failures++;
  }
}

struct Loaded {
  std::ostringstream sink;
  Config cfg{sink};
  CNuc compound;
};

std::unique_ptr<Loaded> Load(const std::string &file, bool brune, bool initialize) {
  auto m = std::make_unique<Loaded>();
  m->cfg.configfile = file;
  if (m->cfg.ReadConfigFile() < 0) return nullptr;
  m->cfg.fileCheckMask = m->cfg.screenCheckMask = 0;
  if (!brune) m->cfg.paramMask &= ~(Config::USE_PARK_FORMALISM | Config::USE_BRUNE_FORMALISM);
  m->cfg.ResolveMatrixRoute();
  InitializeECAmplitudeCache(m->cfg);
  if (m->compound.Fill(m->cfg) == -1) return nullptr;
  if (initialize) m->compound.Initialize(m->cfg);
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
  const CNuc &cn = filled->compound;
  const LevelInputs &in = cn.Inputs();
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
        Check(in.ZeroInputWidth(j, la, ch) == level->FileZeroWidth(ch), where + " channel " + std::to_string(ch) +
                                                                            " zero flag");
      }
    }
  // (b)
  const CNuc &ci = initialized->compound;
  Check(ci.Inputs().TransformedIn() && ci.Inputs().TransformError().empty(), tag + ": transformed, no error");
  bool moved = false;
  for (int j = 1; j <= ci.NumJGroups(); j++)
    for (int la = 1; la <= ci.GetJGroup(j)->NumLevels(); la++) {
      const ALevel *a = ci.GetJGroup(j)->GetLevel(la), *b = cn.GetJGroup(j)->GetLevel(la);
      Check(a->FileE() == b->FileE(), tag + ": the level's file energy is not written by the transformation");
      for (int ch = 1; ch <= ci.GetJGroup(j)->NumChannels(); ch++)
        Check(a->FileGamma(ch) == b->FileGamma(ch), tag + ": the level's file width is not written");
      moved = moved || ci.Inputs().E(j, la) != a->FileE();
    }
  bool several = false;
  for (int j = 1; j <= ci.NumJGroups(); j++) several = several || ci.GetJGroup(j)->NumLevels() > 1;
  if (!brune && several) Check(moved, tag + ": a formal input energy differs from the file's");
  // (c)
  vector_r v(ci.ParameterLayout().Size());
  for (std::size_t i = 0; i < v.size(); i++) v[i] = 1.0 + 0.01 * (double)i;
  initialized->compound.FillCompoundFromParamsPhysical(v);
  const vector_r p = records.Empty() ? v : records.WithPartialWidths(v);
  for (const LevelParameterLayout::Level &slots : ci.ParameterLayout().Levels()) {
    Check(ci.Inputs().E(slots.group, slots.level) == p[slots.energy], tag + ": physical energy");
    for (int ch = 1; ch <= (int)slots.amplitude.size(); ch++)
      Check(ci.Inputs().Gamma(slots.group, slots.level, ch) == p[slots.amplitude[ch - 1]], tag + ": physical width");
  }
  if (!records.Empty()) Check(ci.Inputs().GammaDecayNumbers() == records.NumbersAt(v), tag + ": physical numbers");
}

}  // namespace

int main(int argc, char **argv) {
  for (int a = 1; a < argc; a++) {
    OnProject(argv[a], true);
    OnProject(argv[a], false);
  }
  std::printf("%s: %d failure(s)\n", failures ? "FAIL" : "ok", failures);
  return failures ? 1 : 0;
}
