#include "ThmVertexBoundary.h"
#include "CoulombTable.h"
#include "AChannel.h"
#include "ALevel.h"
#include "CNuc.h"
#include "Config.h"
#include "ESegment.h"
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

ThmVertexBoundary::ThmVertexBoundary(const Config &configure, PPair *pair, JGroup *jgroup, int channel)
    : rule_(RuleOf(configure)),
      channel_(channel),
      l_(jgroup->GetChannel(channel)->GetL()),
      two_j_(pair->SpinOrbitTwoJ(jgroup->GetJ())),
      lowest_(0),
      fixed_(jgroup->GetChannel(channel)->GetBoundaryCondition()),
      channelFunc_(pair, !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC)) {
  double eMin = 0.0;
  for (int la = 1; la <= jgroup->NumLevels(); la++) {
    ALevel *level = jgroup->GetLevel(la);
    if (!level->IsInRMatrix()) continue;
    if (!lowest_ || level->GetFitE() < eMin) {
      eMin = level->GetFitE();
      lowest_ = la;
    }
  }
  // (A J group without levels in the R matrix has no level term; its
  // level-free term, if any, takes S_c at the excitation energy 0.)
  if (rule_ == LOWEST_LEVEL_SHIFT)
    fixed_ = ShiftAtLevelEnergy(pair, l_, two_j_, eMin, !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC));
}

double ThmVertexBoundary::Level(const ALevel *level) const {
  return rule_ == LEVEL_SHIFT ? level->GetShiftFunction(channel_) : fixed_;
}

complex ThmVertexBoundary::OnShellAt(double e) const {
  return complex(channelFunc_.Shift(l_, two_j_, e), channelFunc_.Penetrability(l_, two_j_, e));
}

/*
 * The memo is keyed by what S_c depends on -- charges, reduced mass, channel
 * radius, the hybrid-potential state of the pair, l, the resonance energy and
 * the Coulomb routine -- not by the PPair address: the memo outlives a CNuc
 * (thread_local), a later CNuc in the same process (pyazr sessions one after
 * the other, a radius scan) can get the same address back from the allocator,
 * and a pair's radius can change in place.  Keyed by the address, a session at
 * 6.1 fm read S_c of an earlier 4.1 fm session (18O(p,a) narrow levels: chi2
 * 35.99 instead of 208.6 for the same parameters).  Within one evaluation
 * every point asks for the same few (pair, l, energy) triples.
 */
double ThmVertexBoundary::ShiftAtLevelEnergy(PPair *pair, int l, int two_j, double levelEnergy, bool useGSL) {
  struct Entry {
    int z1, z2; double redmass, radius; std::optional<azure2::physics::HybridPotential> potential; int l, two_j; double e; bool gsl; const void *table; double s;
  };
  static const int kMemo = 32;
  thread_local Entry memo[kMemo];
  thread_local int filled = 0, next = 0;
  const int z1 = pair->GetZ(1), z2 = pair->GetZ(2);
  const double redmass = pair->GetRedMass(), radius = pair->GetChRad();
  const std::optional<azure2::physics::HybridPotential> &potential = pair->GetHybridPotential();
  // The Coulomb table in use (CoulFunc and ShftFunc read S from it where it
  // covers the channel), so a session with and one without it do not share values.
  const void *table = (g_config && g_config->useCoulombTable) ? CoulombTable::Shared(g_config->coulombTablePath) : nullptr;
  double resonanceEnergy = levelEnergy - (pair->GetSepE() + pair->GetExE());
  for (int i = 0; i < filled; i++) {
    const Entry &m = memo[i];
    if (m.z1 == z1 && m.z2 == z2 && m.redmass == redmass && m.radius == radius &&
        m.potential == potential && m.l == l && m.two_j == two_j && m.e == resonanceEnergy && m.gsl == useGSL && m.table == table)
      return m.s;
  }
  // Continuous through threshold, and at it (ChannelFunc).
  double s = ChannelFunc(pair, useGSL).Shift(l, two_j, resonanceEnergy);
  memo[next] = Entry{z1, z2, redmass, radius, potential, l, two_j, resonanceEnergy, useGSL, table, s};
  next = (next + 1) % kMemo;
  if (filled < kMemo) filled++;
  return s;
}

bool ThmVertexAllowsPair(CNuc *compound, const ESegment &segment, int key, bool isTest, const Config &configure) {
  if (!segment.IsTHM() || !compound->IsPairKey(segment.GetEntranceKey())) return true;
  PPair *pair = compound->GetPair(compound->GetPairNumFromKey(segment.GetEntranceKey()));
  if (!pair->GetHybridPotential()) return true;
  configure.outStream << "ERROR (E_THM_HYBRID_POTENTIAL): " << (isTest ? "test segment #" : "segment #") << key
                      << " (THM): entrance pair " << segment.GetEntranceKey()
                      << " has a hybrid nuclear potential (<potential>); the THM vertex uses point-Coulomb waves and "
                         "is not built with a potential.  Switch the potential off for this pair."
                      << std::endl;
  return false;
}
