// The external-capture amplitudes of a point: one amplitude per capture pathway (KGroup, ECMGroup) at the point's
// energy and at each of its sub-points, their storage, and their value at an energy an energy shift has moved the
// point to (from the energy tables, or interpolated between the segment's points).
#include <cmath>
#include <complex>
#include <iostream>

#include "AChannel.h"
#include "ALevel.h"
#include "CNuc.h"
#include "ChannelFunc.h"
#include "ChannelTable.h"
#include "Config.h"
#include "CoulFunc.h"
#include "Decay.h"
#include "ECAmplitudeCache.h"
#include "ECIntegral.h"
#include "ECMGroup.h"
#include "ECTable.h"
#include "ExternalCaptureExit.h"
#include "WorkingCopy.h"
#include "Unit.h"
#include "JGroup.h"
#include "KGroup.h"
#include "Kinematics.h"
#include "PPair.h"

static complex ECPathwayAmplitude(CNuc *, const Config &, PPair *, int, double, KGroup *, ECMGroup *, double,
                                  const ChannelTable *table = nullptr);

/// Returns the external capture amplitude for a given external reaction pathway
/// specified by positions in the KGroup and subsequent ECMGroup vectors.

complex Unit::GetECAmplitude(int kGroupNum, int ecMGroupNum) const {
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
 * of the set-up (ALevel::GetECEnergy, written once by CNuc::Initialize), which a
 * fit or an MCMC that moves the level does not change, so shifted and unshifted
 * segments stay consistent and the memo holds no history of the copy
 * (tests/upgrade/history_test.cpp, project ec).
 *
 * With --no-ec-table the amplitude is interpolated linearly between the
 * energies of the segment's own points (ECAmplitudeCache) and held flat beyond
 * its ends, as before the tables.
 */
complex Unit::GetECAmplitudeWithShift(int kGroupNum, int ecMGroupNum, CNuc *theCNuc, const Config &configure) const {
  // Get current energy including shifts
  double currentEnergy = this->GetCMEnergy();

  if (configure.useECTable) {
    const size_t k = kGroupNum - 1, m = ecMGroupNum - 1;
    if (k < ec_energies_.size() && m < ec_energies_[k].size() && ec_energies_[k][m] == currentEnergy)
      return ec_amplitudes_[k][m];
    // Sized from the requested indices: a point without stored amplitudes (a sub-point or a fit's
    // working copy) has an empty ec_amplitudes_, and sizing from it indexed past the end.
    if (ec_shifted_.size() <= k) ec_shifted_.resize(k + 1);
    if (ec_shifted_[k].size() <= m) ec_shifted_[k].resize(m + 1);
    ShiftedECAmplitude &memo = ec_shifted_[k][m];
    if (memo.valid && memo.energy == currentEnergy) return memo.amplitude;

    const int aa = theCNuc->GetPairNumFromKey(this->GetEntranceKey());
    PPair *entrancePair = theCNuc->GetPair(aa);
    int ir = theCNuc->GetPairNumFromKey(this->GetExitKey());
    KGroup *theKGroup =
        entrancePair->GetDecay(ExternalCaptureExit(*theCNuc, aa, ir, segment_key_))->GetKGroup(kGroupNum);
    ECMGroup *theECMGroup = theKGroup->GetECMGroup(ecMGroupNum);
    int j = theECMGroup->GetJGroupNum();
    const double levelEnergy = theCNuc->GetJGroup(j)->GetLevel(theECMGroup->GetLevelNum())->GetECEnergy();
    memo.amplitude = ECPathwayAmplitude(theCNuc, configure, entrancePair, j, levelEnergy, theKGroup, theECMGroup,
                                        currentEnergy);
    memo.energy = currentEnergy;
    memo.valid = true;
    return memo.amplitude;
  }

  // Use ECAmplitudeCache for interpolation at any energy
  if (g_ecAmplitudeCache) {
    ECAmplitudeCache::AmplitudeKey key;
    key.kGroupNum = kGroupNum;
    key.ecMGroupNum = ecMGroupNum;
    key.entranceKey = this->GetEntranceKey();
    key.exitKey = this->GetExitKey();
    // Use segment key directly - same as used in WorkingCopy::CalculateECAmplitudes
    key.segmentKey = segment_key_;

    // Simple interpolation from cache
    complex interpolatedAmplitude = g_ecAmplitudeCache->GetInterpolatedAmplitude(key, currentEnergy);
    return interpolatedAmplitude;
  }

  // Fallback to original cached amplitude if no global cache available
  return GetECAmplitude(kGroupNum, ecMGroupNum);
}

/*!
 * The external-capture amplitude of one pathway (KGroup, ECMGroup of the decay
 * to the point's exit pair) at c.m. energy cmEnergy, for the final level of
 * J-group j at energy levelEnergy: the entrance Coulomb and hard-sphere phases,
 * the photon-energy factor and the radial integral (from ECTable, or directly
 * with --no-ec-table).  CalculateECAmplitudes calls it at the point's energy;
 * GetECAmplitudeWithShift at a shifted one.
 */
static complex ECPathwayAmplitude(CNuc *theCNuc, const Config &configure, PPair *entrancePair, int j,
                                  double levelEnergy, KGroup *theKGroup, ECMGroup *theECMGroup,
                                  double cmEnergy,
                                  const ChannelTable *table) {
  // A closed entrance channel -- at or below threshold, or so close above it
  // that P = 0 in double precision (ChannelFunc) -- captures nothing: the
  // amplitude carries the entrance Coulomb wave, whose routines return NaN or
  // garbage there (e.g. a Gaussian-folding sub-point at E_cm <= 0).
  if (ChannelFunc(entrancePair, !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC)).IsClosed(cmEnergy))
    return complex(0.0, 0.0);
  double inEnergy = cmEnergy + entrancePair->GetSepE() + entrancePair->GetExE();
  // entrance Phase Calculations;
  CoulFunc theCoulombFunction(entrancePair, !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC));
  // Every pathway of this point with the same l asks for the same exp(-i phi) at a_c (CoulFunc::SetRadialScratch).
  if (!(configure.paramMask & Config::CALCULATE_REACTION_RATE)) theCoulombFunction.SetRadialScratch(true);

  double eta = Kinematics::Channel(entrancePair).Eta(cmEnergy);
  complex expCP(1.0, 0.0);
  for (int ll = 1; ll <= theECMGroup->GetL(); ll++)
    expCP *= complex((double)ll / sqrt(pow((double)ll, 2.0) + pow(eta, 2.0)),
                     eta / sqrt(pow((double)ll, 2.0) + pow(eta, 2.0)));
  // exp(-i phi) at a_c from the point's channel table where it holds the open channel, at E_c = (E - S) - e
  // (the table's expression; within one ulp of cmEnergy); otherwise computed.
  int entranceNum = 0;
  for (int p = 1; p <= theCNuc->NumPairs() && !entranceNum; p++)
    if (theCNuc->GetPair(p) == entrancePair) entranceNum = p;
  const int entranceTwoJ = entrancePair->SpinOrbitTwoJ(theECMGroup->GetJ());
  const ChannelTable::Row *row = table ? table->Find({entranceNum, theECMGroup->GetL(), entranceTwoJ,
                                                      (inEnergy - entrancePair->GetSepE()) - entrancePair->GetExE()})
                                       : nullptr;
  complex expHSP = (row && row->P > 0.0) ? row->expHSP
                                         : theCoulombFunction.HardSpherePhase(theECMGroup->GetL(), entranceTwoJ, entrancePair->GetChRad(), cmEnergy);

  double sqrtGammaPene = pow((inEnergy - levelEnergy) / hbarc, theECMGroup->GetMult() + 0.5);

  AChannel *theFinalChannel = theCNuc->GetJGroup(j)->GetChannel(theECMGroup->GetFinalChannel());
  PPair *theFinalPair = theCNuc->GetPair(theFinalChannel->GetPairNum());
  int theInitialLValue;
  double theInitialSValue;
  if (theECMGroup->IsChannelCapture()) {
    MGroup *theChanCapMGroup = entrancePair->GetDecay(theECMGroup->GetChanCapDecay())->GetKGroup(theECMGroup->GetChanCapKGroup())->GetMGroup(theECMGroup->GetChanCapMGroup());
    theInitialLValue = theCNuc->GetJGroup(theChanCapMGroup->GetJNum())->GetChannel(theChanCapMGroup->GetChpNum())->GetL();
    theInitialSValue = theCNuc->GetJGroup(theChanCapMGroup->GetJNum())->GetChannel(theChanCapMGroup->GetChpNum())->GetS();
  } else {
    theInitialLValue = theECMGroup->GetL();
    theInitialSValue = theKGroup->GetS();
  }

  complex integrals;
  if (configure.useECTable) {
    integrals = ECTable::Integral(theFinalPair, configure, theInitialLValue, theFinalChannel->GetL(),
                                  theInitialSValue, theFinalChannel->GetS(),
                                  theECMGroup->GetJ(), theCNuc->GetJGroup(j)->GetJ(),
                                  theECMGroup->GetMult(), theECMGroup->GetRadType(),
                                  inEnergy, levelEnergy, theECMGroup->IsChannelCapture(), table,
                                  theFinalChannel->GetPairNum());
  } else {
    ECIntegral theECIntegral(theFinalPair, configure);
    integrals = theECIntegral(theInitialLValue, theFinalChannel->GetL(),
                              theInitialSValue, theFinalChannel->GetS(),
                              theECMGroup->GetJ(), theCNuc->GetJGroup(j)->GetJ(),
                              theECMGroup->GetMult(), theECMGroup->GetRadType(),
                              inEnergy, levelEnergy,
                              theECMGroup->IsChannelCapture());
  }

  // calculate the total radial integral
  return expCP * expHSP * sqrtGammaPene * integrals;
}

/*!
 * Calculates the external capture amplitudes for the data point.
 * The amplitudes are calculated for all reaction pathways with corresponding
 * entrance and exit pairs.
 */

void Unit::CalculateECAmplitudes(WorkingCopy *data, CNuc *theCNuc, const Config &configure, const ChannelTable *table) {
  int aa = theCNuc->GetPairNumFromKey(this->GetEntranceKey());
  if (theCNuc->GetPair(aa)->GetPType() == 20) return;
  if (theCNuc->GetPair(aa)->IsEntrance()) {
    PPair *entrancePair = theCNuc->GetPair(aa);
    for (int j = 1; j <= theCNuc->NumJGroups(); j++) {
      for (int la = 1; la <= theCNuc->GetJGroup(j)->NumLevels(); la++) {
        if (theCNuc->GetJGroup(j)->GetLevel(la)->IsECLevel()) {
          ALevel *ecLevel = theCNuc->GetJGroup(j)->GetLevel(la);
          int ir = theCNuc->GetPairNumFromKey(this->GetExitKey());
          if (ecLevel->GetECPairNum() == ir) {
            Decay *decay = entrancePair->GetDecay(ExternalCaptureExit(*theCNuc, aa, ir, segment_key_));
            for (int k = 1; k <= decay->NumKGroups(); k++) {
              KGroup *theKGroup = decay->GetKGroup(k);
              for (int ecm = 1; ecm <= theKGroup->NumECMGroups(); ecm++) {
                ECMGroup *theECMGroup = theKGroup->GetECMGroup(ecm);
                complex ecAmplitude = ECPathwayAmplitude(theCNuc, configure, entrancePair, j, ecLevel->GetECEnergy(),
                                                         theKGroup, theECMGroup, this->GetCMEnergy(), table);
                this->AddECAmplitude(k, ecm, ecAmplitude);
              }
            }
          }
        }
      }
    }
  }
  for (int i = 1; i <= this->NumSubPoints(); i++) {
    this->GetSubPoint(*data, i)->CalculateECAmplitudes(data, theCNuc, configure, table);
  }
}

/*!
 * Adds an external capture amplitude with reference to a specified reaction pathway.
 */

void Unit::AddECAmplitude(int kGroupNum, int ecMGroupNum, complex ecAmplitude, double energy) {
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
  if (g_ecAmplitudeCache && !(g_config && g_config->useECTable)) {
    ECAmplitudeCache::AmplitudeKey key;
    key.kGroupNum = kGroupNum;
    key.ecMGroupNum = ecMGroupNum;
    key.entranceKey = this->GetEntranceKey();
    key.exitKey = this->GetExitKey();
    key.segmentKey = segment_key_;
    g_ecAmplitudeCache->AddAmplitude(key, energy, ecAmplitude);

    // A non-finite amplitude (never cached, see ECAmplitudeCache::AddAmplitude) is replaced by the value
    // interpolated from the neighbouring cached energies.
    if (!std::isfinite(std::real(ecAmplitude)) || !std::isfinite(std::imag(ecAmplitude))) {
      std::cout << "WARNING: Replacing EC Amplitude with interpolated value from cache for kGroup "
                << kGroupNum << ", ecMGroup " << ecMGroupNum
                << " at energy " << energy << " MeV." << std::endl;
      complex interpAmplitude = g_ecAmplitudeCache->GetInterpolatedAmplitude(key, energy);
      std::cout << "         Original Amplitude: " << ecAmplitude
                << ", Interpolated Amplitude: " << interpAmplitude << std::endl;
      ec_amplitudes_[kGroupNum - 1].back() = interpAmplitude;
    }
  }
}

void Unit::AddECAmplitude(int kGroupNum, int ecMGroupNum, complex ecAmplitude) {
  // Default behavior - use current CM energy
  double currentEnergy = this->GetCMEnergy();
  AddECAmplitude(kGroupNum, ecMGroupNum, ecAmplitude, currentEnergy);
}

/*!
 * Clears all EC amplitudes and energies. Used when cloning to avoid duplicate accumulation.
 */
void Unit::ClearECAmplitudes() {
  ec_amplitudes_.clear();
  ec_energies_.clear();
  ec_shifted_.clear();
}
