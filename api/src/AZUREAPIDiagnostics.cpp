// AZUREAPI, the diagnostics of the external region: the external-capture radial integrals on a requested energy
// grid (GetECIntegrals) and the counters of the input-file cache (GetCacheStats).  One unit because neither is
// part of a calculation: they read the model as set up and return numbers for inspection (pyazr's
// ec_integrals(), cache_stats()).
#include <vector>

#include "AZUREAPI.h"
#include "ALevel.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "ECIntegral.h"
#include "InputFiles.h"
#include "JGroup.h"
#include "PPair.h"

// ---------------------------------------------------------------------------
//  Diagnostics: the external region, and the caches that make it affordable
// ---------------------------------------------------------------------------

/*!
 * External-capture radial integrals on a requested energy grid.
 *
 * Walks exactly the pathway structure Unit::CalculateECAmplitudes walks --
 * every EC level whose final pair matches, every spin group of the exit to it in
 * the pathway table, every external pathway within it -- and evaluates the integral at each
 * requested energy instead of at a data point's energy.  The quantum numbers of
 * each pathway come back with it, so no bookkeeping is needed on the caller's
 * side.
 */
vector_r AZUREAPI::GetECIntegrals(const vector_r &request) const {
  vector_r out;
  if (request.size() < 2) return out;

  int pairKey = static_cast<int>(std::lround(request[0]));
  int nE = static_cast<int>(std::lround(request[1]));
  if (nE <= 0 || request.size() < static_cast<size_t>(2 + nE)) return out;
  std::vector<double> energies(request.begin() + 2, request.begin() + 2 + nE);

  if (!compound_ || !compound_->IsPairKey(pairKey)) return out;
  int aa = compound_->GetPairNumFromKey(pairKey);
  const PPair *entrancePair = compound_->GetPair(aa);
  if (!entrancePair || !entrancePair->IsEntrance()) return out;

  // One block per pathway: 6 descriptors then 2*nE numbers.
  std::vector<vector_r> blocks;

  for (int j = 1; j <= compound_->NumJGroups(); j++) {
    for (int la = 1; la <= compound_->GetJGroup(j)->NumLevels(); la++) {
      const ALevel *ecLevel = compound_->GetJGroup(j)->GetLevel(la);
      if (!ecLevel->IsECLevel()) continue;
      int ir = ecLevel->GetECPairNum();
      if (ir < 1 || ir > compound_->NumPairs()) continue;
      // The exit to the final state's photon pair; none when no pathway leads there from this entrance.
      const azure2::EntrancePathways &P = compound_->Pathways().Entrance(aa);
      const int exit = P.ExitOf(ir);
      if (!exit) continue;

      for (int k = 1; k <= P.exits[exit - 1].num_groups; k++) {
        const azure2::SpinGroupRow &theKGroup = P.Group(P.exits[exit - 1], k);
        for (int ecm = 1; ecm <= theKGroup.num_external; ecm++) {
          const azure2::ExternalRow &g = P.External(theKGroup, ecm);

          const AChannel *finalChannel =
              compound_->GetJGroup(j)->GetChannel(g.final_channel);
          const PPair *finalPair = compound_->GetPair(finalChannel->GetPairNum());

          int liValue;
          double siValue;
          if (g.channel_capture) {
            const azure2::InternalRow &cc = P.internal[g.cc_row];
            liValue = compound_->GetJGroup(cc.group)->GetChannel(cc.chp)->GetL();
            siValue = compound_->GetJGroup(cc.group)->GetChannel(cc.chp)->GetS();
          } else {
            liValue = g.l;
            siValue = theKGroup.s;
          }

          vector_r block;
          block.push_back(static_cast<double>(liValue));
          block.push_back(static_cast<double>(finalChannel->GetL()));
          block.push_back(2.0 * siValue);
          block.push_back(2.0 * finalChannel->GetS());
          block.push_back(static_cast<double>(g.mult));
          block.push_back(g.rad_type == 'E' ? 1.0 : 0.0);

          ECIntegral theECIntegral(finalPair, configure());
          // The final state's capture energy, as the engine's integrals use it (ALevel::GetECEnergy, set by
          // CompoundSystem::SetUpECFinalBlocks): in the formal basis GetE() is the formal energy, not that one.
          double levelEnergy = ecLevel->GetECEnergy();
          for (int i = 0; i < nE; ++i) {
            double inEnergy = energies[i] + entrancePair->GetSepE() + entrancePair->GetExE();
            complex v(0., 0.);
            if (energies[i] > 0.0) {
              try {
                v = theECIntegral(liValue, finalChannel->GetL(),
                                  siValue, finalChannel->GetS(),
                                  g.J, compound_->GetJGroup(j)->GetJ(),
                                  g.mult, g.rad_type,
                                  inEnergy, levelEnergy,
                                  g.channel_capture);
              } catch (...) {
                v = complex(0., 0.);
              }
            }
            block.push_back(real(v));
            block.push_back(imag(v));
          }
          blocks.push_back(std::move(block));
        }
      }
    }
  }

  out.push_back(static_cast<double>(blocks.size()));
  out.push_back(static_cast<double>(nE));
  for (const vector_r &b : blocks) out.insert(out.end(), b.begin(), b.end());
  return out;
}

/*!
 * Six zeros (the counters of the deleted Coulomb memo), then the input-file counters of the process
 * (InputFiles: files read from disk, opens served from memory, files held).
 */
vector_r AZUREAPI::GetCacheStats() const {
  vector_r out(9, 0.0);
  const InputFiles::Stats in = InputFiles::GetStats();
  out[6] = static_cast<double>(in.reads);
  out[7] = static_cast<double>(in.served);
  out[8] = static_cast<double>(in.files);
  // out[0..5] were the counters of the process-wide Coulomb memo (CoulFuncCache), deleted with
  // docs/upgrade/coulomb.md: they stay 0, so that the layout and pyazr's cache_stats() keys do not change.
  return out;
}
