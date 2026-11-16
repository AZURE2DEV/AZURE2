// The external-capture amplitudes of a point: one amplitude per capture pathway (spin group k and external pathway m
// of the exit in the pathway table, PathwayTable.h) at the point's
// energy and at each of its sub-points, their storage, and their value at an energy an energy shift has moved the
// point to (from the energy tables, or interpolated between the segment's points).
#include "RunScratch.h"
#include <cmath>
#include <complex>
#include <iostream>

#include "AChannel.h"
#include "ALevel.h"
#include "CompoundSystem.h"
#include "ChannelFunc.h"
#include "ChannelTable.h"
#include "Config.h"
#include "CoulFunc.h"
#include "ECAmplitudeCache.h"
#include "ECIntegral.h"
#include "ECTable.h"
#include "ExternalCaptureExit.h"
#include "WorkingCopy.h"
#include "Unit.h"
#include "JGroup.h"
#include "Kinematics.h"
#include "PPair.h"
#include "PathwayTable.h"

static complex ECPathwayAmplitude(const CompoundSystem *, const Config &, int entrancePair, const azure2::EntrancePathways &,
                                  const azure2::SpinGroupRow &, const azure2::ExternalRow &, double levelEnergy,
                                  double cmEnergy, const ChannelTable *table = nullptr);

/// Returns the external capture amplitude of external pathway m (1-based) of spin group k (1-based) of the point's
/// exit, as stored at set-up.

complex Site::GetECAmplitude(int kGroupNum, int ecMGroupNum) const {
  return ec_amplitudes_[kGroupNum - 1][ecMGroupNum - 1];
}

/*!
 * Returns the external capture amplitude at the point's current energy, which
 * an energy shift may have moved away from the energy the amplitude was
 * computed at.
 *
 * With the energy tables (Config::useECTable, the default) a moved point is
 * evaluated at its new energy, from the same tables as every other amplitude,
 * and the result is kept until the energy moves again.  The final level's
 * energy is the one the point's own amplitudes were computed with: the energy
 * of the set-up (ALevel::GetECEnergy, written once by CompoundSystem::Initialize), which a
 * fit or an MCMC that moves the level does not change, so shifted and unshifted
 * segments stay consistent and the memo holds no history of the copy
 * (tests/upgrade/history_test.cpp, project ec).
 *
 * With --no-ec-table the amplitude is interpolated linearly between the
 * energies of the segment's own points (ECAmplitudeCache) and held flat beyond
 * its ends, as before the tables.
 */
complex Site::GetECAmplitudeWithShift(int kGroupNum, int ecMGroupNum, const CompoundSystem *compound, const Config &configure) {
  // Get current energy including shifts
  double currentEnergy = this->GetCMEnergy();

  if (configure.useECTable) {
    const size_t k = kGroupNum - 1, m = ecMGroupNum - 1;
    if (k < ec_energies_.size() && m < ec_energies_[k].size() && ec_energies_[k][m] == currentEnergy)
      return ec_amplitudes_[k][m];
    // Sized from the requested indices: a point without stored amplitudes (a sub-point or a fit's
    // working copy) has an empty ec_amplitudes_, and sizing from it indexed past the end.
    if (values_.ec_shifted.size() <= k) values_.ec_shifted.resize(k + 1);
    if (values_.ec_shifted[k].size() <= m) values_.ec_shifted[k].resize(m + 1);
    ShiftedECAmplitude &memo = values_.ec_shifted[k][m];
    if (memo.valid && memo.energy == currentEnergy) return memo.amplitude;

    const int aa = compound->GetPairNumFromKey(this->GetEntranceKey());
    const azure2::EntrancePathways &P = compound->Pathways().Entrance(aa);
    const azure2::SpinGroupRow &group =
        P.Group(ExternalCaptureExit(*compound, aa, compound->GetPairNumFromKey(this->GetExitKey()), segment_key_), kGroupNum);
    const azure2::ExternalRow &pathway = P.External(group, ecMGroupNum);
    const double levelEnergy = compound->GetJGroup(pathway.final_group)->GetLevel(pathway.final_level)->GetECEnergy();
    memo.amplitude = ECPathwayAmplitude(compound, configure, aa, P, group, pathway, levelEnergy, currentEnergy);
    memo.energy = currentEnergy;
    memo.valid = true;
    return memo.amplitude;
  }

  // The run's ECAmplitudeCache for interpolation at any energy
  if (ECAmplitudeCache *cache = RunECAmplitudes(configure)) {
    ECAmplitudeCache::AmplitudeKey key;
    key.kGroupNum = kGroupNum;
    key.ecMGroupNum = ecMGroupNum;
    key.entranceKey = this->GetEntranceKey();
    key.exitKey = this->GetExitKey();
    // Use segment key directly - same as used in WorkingCopy::CalculateECAmplitudes
    key.segmentKey = segment_key_;

    // Simple interpolation from cache
    complex interpolatedAmplitude = cache->GetInterpolatedAmplitude(key, currentEnergy);
    return interpolatedAmplitude;
  }

  // The point's own amplitude when the run has no cache
  return GetECAmplitude(kGroupNum, ecMGroupNum);
}

/*!
 * The external-capture amplitude of one external pathway (a row of spin group
 * `group` of the pathways from pair `entrance`) at c.m. energy cmEnergy, for its
 * final level at energy levelEnergy (MeV): the entrance Coulomb and hard-sphere phases,
 * the photon-energy factor and the radial integral (from ECTable, or directly
 * with --no-ec-table).  CalculateECAmplitudes calls it at the point's energy;
 * GetECAmplitudeWithShift at a shifted one.
 */
static complex ECPathwayAmplitude(const CompoundSystem *compound, const Config &configure, int entrance,
                                  const azure2::EntrancePathways &P, const azure2::SpinGroupRow &group,
                                  const azure2::ExternalRow &pathway, double levelEnergy, double cmEnergy,
                                  const ChannelTable *table) {
  const PPair *entrancePair = compound->GetPair(entrance);
  const int j = pathway.final_group;
  // A closed entrance channel -- at or below threshold, or so close above it
  // that P = 0 in double precision (ChannelFunc) -- captures nothing: the
  // amplitude carries the entrance Coulomb wave, whose routines return NaN or
  // garbage there (e.g. a Gaussian-folding sub-point at E_cm <= 0).
  if (ChannelFunc(entrancePair).IsClosed(cmEnergy))
    return complex(0.0, 0.0);
  double inEnergy = cmEnergy + entrancePair->GetSepE() + entrancePair->GetExE();
  // entrance Phase Calculations;
  CoulFunc theCoulombFunction(entrancePair);
  // Every pathway of this point with the same l asks for the same exp(-i phi) at a_c (CoulFunc::SetRadialMemo).
  if (!(configure.paramMask & Config::CALCULATE_REACTION_RATE))
    theCoulombFunction.SetRadialMemo(&configure.scratch->ForThisThread().ec.radial);

  double eta = Kinematics::Channel(entrancePair).Eta(cmEnergy);
  complex expCP(1.0, 0.0);
  for (int ll = 1; ll <= pathway.l; ll++)
    expCP *= complex((double)ll / sqrt(pow((double)ll, 2.0) + pow(eta, 2.0)),
                     eta / sqrt(pow((double)ll, 2.0) + pow(eta, 2.0)));
  // exp(-i phi) at a_c from the point's channel table where it holds the open channel, at E_c = (E - S) - e
  // (the table's expression; within one ulp of cmEnergy); otherwise computed.
  const int entranceTwoJ = entrancePair->SpinOrbitTwoJ(pathway.J);
  const ChannelTable::Row *row = table ? table->Find({entrance, pathway.l, entranceTwoJ,
                                                      (inEnergy - entrancePair->GetSepE()) - entrancePair->GetExE()})
                                       : nullptr;
  complex expHSP = (row && row->P > 0.0) ? row->expHSP
                                         : theCoulombFunction.HardSpherePhase(pathway.l, entranceTwoJ, entrancePair->GetChRad(), cmEnergy);

  double sqrtGammaPene = pow((inEnergy - levelEnergy) / hbarc, pathway.mult + 0.5);

  const AChannel *theFinalChannel = compound->GetJGroup(j)->GetChannel(pathway.final_channel);
  const PPair *theFinalPair = compound->GetPair(theFinalChannel->GetPairNum());
  int theInitialLValue;
  double theInitialSValue;
  if (pathway.channel_capture) {
    // the exit channel of the internal pathway the capture passes through
    const azure2::InternalRow &through = P.internal[pathway.cc_row];
    theInitialLValue = compound->GetJGroup(through.group)->GetChannel(through.chp)->GetL();
    theInitialSValue = compound->GetJGroup(through.group)->GetChannel(through.chp)->GetS();
  } else {
    theInitialLValue = pathway.l;
    theInitialSValue = group.s;
  }

  complex integrals;
  if (configure.useECTable) {
    integrals = ECTable::Integral(theFinalPair, configure, theInitialLValue, theFinalChannel->GetL(),
                                  theInitialSValue, theFinalChannel->GetS(),
                                  pathway.J, compound->GetJGroup(j)->GetJ(), pathway.mult, pathway.rad_type,
                                  inEnergy, levelEnergy, pathway.channel_capture, table,
                                  theFinalChannel->GetPairNum());
  } else {
    ECIntegral theECIntegral(theFinalPair, configure);
    integrals = theECIntegral(theInitialLValue, theFinalChannel->GetL(),
                              theInitialSValue, theFinalChannel->GetS(),
                              pathway.J, compound->GetJGroup(j)->GetJ(), pathway.mult, pathway.rad_type,
                              inEnergy, levelEnergy, pathway.channel_capture);
  }

  // calculate the total radial integral
  return expCP * expHSP * sqrtGammaPene * integrals;
}

/*!
 * Calculates the external capture amplitudes for the data point.
 * The amplitudes are calculated for all reaction pathways with corresponding
 * entrance and exit pairs.
 */

void Site::CalculateECAmplitudes(WorkingCopy *data, const CompoundSystem *compound, const Config &configure, const ChannelTable *table) {
  int aa = compound->GetPairNumFromKey(this->GetEntranceKey());
  if (compound->GetPair(aa)->GetPType() == 20) return;
  if (compound->GetPair(aa)->IsEntrance()) {
    const azure2::EntrancePathways &P = compound->Pathways().Entrance(aa);
    for (int j = 1; j <= compound->NumJGroups(); j++) {
      for (int la = 1; la <= compound->GetJGroup(j)->NumLevels(); la++) {
        if (compound->GetJGroup(j)->GetLevel(la)->IsECLevel()) {
          const ALevel *ecLevel = compound->GetJGroup(j)->GetLevel(la);
          int ir = compound->GetPairNumFromKey(this->GetExitKey());
          if (ecLevel->GetECPairNum() == ir) {
            const azure2::ExitRow &exit = ExternalCaptureExit(*compound, aa, ir, segment_key_);
            for (int k = 1; k <= exit.num_groups; k++) {
              const azure2::SpinGroupRow &group = P.Group(exit, k);
              for (int ecm = 1; ecm <= group.num_external; ecm++) {
                complex ecAmplitude = ECPathwayAmplitude(compound, configure, aa, P, group, P.External(group, ecm),
                                                         ecLevel->GetECEnergy(), this->GetCMEnergy(), table);
                this->AddECAmplitude(k, ecm, ecAmplitude, configure);
              }
            }
          }
        }
      }
    }
  }
}

/*!
 * Adds an external capture amplitude with reference to a specified reaction pathway.
 */

void Site::AddECAmplitude(int kGroupNum, int ecMGroupNum, complex ecAmplitude, double energy, const Config &configure) {
  vector_c d_amp;
  vector_r d_energy;
  while (kGroupNum > ec_amplitudes_.size()) {
    ec_amplitudes_.push_back(d_amp);
    ec_energies_.push_back(d_energy);
  }
  ec_amplitudes_[kGroupNum - 1].push_back(ecAmplitude);
  ec_energies_[kGroupNum - 1].push_back(energy);
  assert(ecMGroupNum == ec_amplitudes_[kGroupNum - 1].size());
  assert(ec_amplitudes_[kGroupNum - 1].size() == ec_energies_[kGroupNum - 1].size());

  // Add to shared cache for interpolation across different units.  The energy
  // tables evaluate shifted points directly (GetECAmplitudeWithShift), so the
  // cache is only filled without them.
  if (ECAmplitudeCache *cache = configure.useECTable ? nullptr : RunECAmplitudes(configure)) {
    ECAmplitudeCache::AmplitudeKey key;
    key.kGroupNum = kGroupNum;
    key.ecMGroupNum = ecMGroupNum;
    key.entranceKey = this->GetEntranceKey();
    key.exitKey = this->GetExitKey();
    key.segmentKey = segment_key_;
    cache->AddAmplitude(key, energy, ecAmplitude);

    // A non-finite amplitude (never cached, see ECAmplitudeCache::AddAmplitude) is replaced by the value
    // interpolated from the neighbouring cached energies.
    if (!std::isfinite(std::real(ecAmplitude)) || !std::isfinite(std::imag(ecAmplitude))) {
      std::cout << "WARNING: Replacing EC Amplitude with interpolated value from cache for kGroup "
                << kGroupNum << ", ecMGroup " << ecMGroupNum
                << " at energy " << energy << " MeV." << std::endl;
      complex interpAmplitude = cache->GetInterpolatedAmplitude(key, energy);
      std::cout << "         Original Amplitude: " << ecAmplitude
                << ", Interpolated Amplitude: " << interpAmplitude << std::endl;
      ec_amplitudes_[kGroupNum - 1].back() = interpAmplitude;
    }
  }
}

void Site::AddECAmplitude(int kGroupNum, int ecMGroupNum, complex ecAmplitude, const Config &configure) {
  // Default behavior - use current CM energy
  double currentEnergy = this->GetCMEnergy();
  AddECAmplitude(kGroupNum, ecMGroupNum, ecAmplitude, currentEnergy, configure);
}

/*!
 * Clears all EC amplitudes and energies. Used when cloning to avoid duplicate accumulation.
 */
void Site::ClearECAmplitudes() {
  ec_amplitudes_.clear();
  ec_energies_.clear();
  values_.ec_shifted.clear();
}
