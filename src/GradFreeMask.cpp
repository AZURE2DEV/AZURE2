// BuildFreeMask (GradFreeMask.h), out of AZUREGrad.cpp (the size ratchet).
#include "GradFreeMask.h"

#include <vector>

#include "CompoundSystem.h"

// The free level parameters of the packed vector, per J-group (FormalBasis).  An isospin-derived slot is fixed in the
// packed vector but its derivative reaches free parameters through the chain rule of GradAccum::Scatter, so it counts
// as free when any of its sources does.
AMatrixFunc::FreeMask BuildFreeMask(const CompoundSystem *compound, const ParamIndexMap &pmap) {
  AMatrixFunc::FreeMask mask;
  const int nJ = compound->NumJGroups();
  mask.e.assign(nJ, {});
  mask.g.assign(nJ, {});
  std::vector<char> needed(pmap.NumFull(), 0);
  for (int f = 0; f < pmap.NumFull(); f++) needed[f] = pmap.FullToPacked(f) >= 0;
  for (const IsospinTie &t : pmap.IsospinTies()) {
    if (t.slot < 0 || t.slot >= pmap.NumFull() || needed[t.slot]) continue;
    for (const auto &src : t.sources)
      if (src.first >= 0 && src.first < pmap.NumFull() && needed[src.first]) needed[t.slot] = 1;
  }
  auto isFree = [&needed](int full) { return full >= 0 && full < (int)needed.size() && needed[full]; };
  for (int j = 1; j <= nJ; j++) {
    const JGroup *jg = compound->GetJGroup(j);
    const int nl = jg->NumLevels(), nch = jg->NumChannels();
    mask.e[j - 1].assign(nl, 0);
    mask.g[j - 1].assign(nl, std::vector<char>(nch, 0));
    for (int la = 1; la <= nl; la++) {
      mask.e[j - 1][la - 1] = isFree(pmap.EIndexFlat(j, la));
      for (int ch = 1; ch <= nch; ch++) mask.g[j - 1][la - 1][ch - 1] = isFree(pmap.GammaIndexFlat(j, la, ch));
      // A gamma_decay level: its photon amplitudes move with any free number of the level (the record's, its energy,
      // its particle amplitudes), so their cotangents are needed for the pull-back (GammaDecayMap.h).
      // With external capture also with the final states' amplitudes (GammaDecaySet::InputIndices).
      const GammaDecaySet &decays = compound->GetGammaDecays();
      const int record = decays.RecordOf(j, la);
      if (record < 0) continue;
      bool anyFree = false;
      for (const int index : decays.InputIndices(pmap, record)) anyFree = anyFree || isFree(index);
      if (anyFree)
        for (const int ch : decays.Levels()[record].photonChannels) mask.g[j - 1][la - 1][ch - 1] = 1;
    }
  }
  return mask;
}
