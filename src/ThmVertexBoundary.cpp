#include "ThmVertexBoundary.h"
#include "RunScratch.h"
#include "AChannel.h"
#include "ALevel.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "Segment.h"
#include "JGroup.h"
#include "PPair.h"

ThmVertexBoundary::Rule ThmVertexBoundary::RuleOf(const Config &configure) {
  switch (configure.thm.vertex) {
    case Config::ThmOptions::ON_SHELL: return ON_SHELL;
    case Config::ThmOptions::CONSTANT: return LOWEST_LEVEL_SHIFT;
    case Config::ThmOptions::PER_LEVEL:
      return (configure.paramMask & Config::USE_BRUNE_FORMALISM) ? LEVEL_SHIFT : CHANNEL_CONSTANT;
  }
  return CHANNEL_CONSTANT;
}

ThmVertexBoundary::ThmVertexBoundary(const Config &configure, const PPair *pair, const JGroup *jgroup,
                                     const LevelState &state, int j, int channel)
    : rule_(RuleOf(configure)),
      state_(&state),
      j_(j),
      channel_(channel),
      l_(jgroup->GetChannel(channel)->GetL()),
      two_j_(pair->SpinOrbitTwoJ(jgroup->GetJ())),
      lowest_(0),
      fixed_(jgroup->GetChannel(channel)->GetBoundaryCondition()),
      channelFunc_(pair) {
  double eMin = 0.0;
  for (int la = 1; la <= jgroup->NumLevels(); la++) {
    const ALevel *level = jgroup->GetLevel(la);
    if (!level->IsInRMatrix()) continue;
    if (!lowest_ || state.E(j, la) < eMin) {
      eMin = state.E(j, la);
      lowest_ = la;
    }
  }
  // (A J group without levels in the R matrix has no level term; its
  // level-free term, if any, takes S_c at the excitation energy 0.)
  if (rule_ == LOWEST_LEVEL_SHIFT)
    fixed_ = ShiftAtLevelEnergy(pair, l_, two_j_, eMin,
                                &configure.scratch->ForThisThread().thmShift);
}

double ThmVertexBoundary::Level(int la) const {
  return rule_ == LEVEL_SHIFT ? state_->Shift(j_, la, channel_) : fixed_;
}

complex ThmVertexBoundary::OnShellAt(double e) {
  return complex(channelFunc_.Shift(l_, two_j_, e), channelFunc_.Penetrability(l_, two_j_, e));
}

double ThmVertexBoundary::ShiftAtLevelEnergy(const PPair *pair, int l, int two_j, double levelEnergy,
                                             ThmShiftMemo *memo) {
  const double resonanceEnergy = levelEnergy - (pair->GetSepE() + pair->GetExE());
  // Continuous through threshold, and at it (ChannelFunc).
  if (!memo) return ChannelFunc(pair).Shift(l, two_j, resonanceEnergy);
  const int z1 = pair->GetZ(1), z2 = pair->GetZ(2);
  const double redmass = pair->GetRedMass(), radius = pair->GetChRad();
  const std::optional<azure2::physics::HybridPotential> &potential = pair->GetHybridPotential();
  for (int i = 0; i < memo->filled; i++) {
    const ThmShiftMemo::Entry &m = memo->entries[i];
    if (m.z1 == z1 && m.z2 == z2 && m.redmass == redmass && m.radius == radius && m.potential == potential &&
        m.l == l && m.two_j == two_j && m.e == resonanceEnergy)
      return m.s;
  }
  const double s = ChannelFunc(pair).Shift(l, two_j, resonanceEnergy);
  memo->entries[memo->next] = ThmShiftMemo::Entry{z1, z2, redmass, radius, potential, l, two_j, resonanceEnergy, s};
  memo->next = (memo->next + 1) % ThmShiftMemo::kSize;
  if (memo->filled < ThmShiftMemo::kSize) memo->filled++;
  return s;
}

bool ThmVertexAllowsPair(const CompoundSystem *compound, const Segment &segment, int key, bool isTest, const Config &configure) {
  if (!segment.IsTHM() || !compound->IsPairKey(segment.GetEntranceKey())) return true;
  const PPair *pair = compound->GetPair(compound->GetPairNumFromKey(segment.GetEntranceKey()));
  if (!pair->GetHybridPotential()) return true;
  configure.outStream << "ERROR (E_THM_HYBRID_POTENTIAL): " << (isTest ? "test segment #" : "segment #") << key
                      << " (THM): entrance pair " << segment.GetEntranceKey()
                      << " has a hybrid nuclear potential (<potential>); the THM vertex uses point-Coulomb waves and "
                         "is not built with a potential.  Switch the potential off for this pair."
                      << std::endl;
  return false;
}
