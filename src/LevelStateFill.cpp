#include <cmath>
#include <exception>
#include <mutex>
#include "azure2/support/FloatingPointPrecision.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "CompoundSystem.h"
#include "ChannelFunc.h"
#include "ExternalCaptureAmplitudes.h"
#include "Config.h"
#include "CoulFunc.h"
#include "RunScratch.h"
#include "ShftFunc.h"

// The members of CompoundSystem that write its LevelState at a parameter vector: the fit values, the shift functions and what
// follows from them.

/*!
 * Fills the level state from a parameter vector (FillMnParams' order), with everything derived from it but the
 * shift functions (CalcShiftFunctions).
 */

void CompoundSystem::FillCompoundFromParams(const vector_r &params, LevelState &state) const {
  // IsospinTies.h: the derived members of the isospin levels follow the free parameters of the vector.
  // A vector that holds the isospin parameters sets them; a shorter one (the level block alone) keeps the state's.
  vector_r tied;
  if (isospin_->Active()) {
    if ((int)params.size() >= this->NumParameters()) state.SetIsospinValues(isospin_->ValuesIn(params));
    tied = params;
    isospin_->Apply(tied, state.IsospinValues());
  }
  const vector_r &untransformed = tied.empty() ? params : tied;
  // GammaDecayMap.h: the records' numbers become amplitudes; what the map computed is the state's.
  const vector_r &p = gammaDecays_->Empty() ? untransformed : gammaDecays_->ToEngine(*this, state, untransformed);
  for (const LevelParameterLayout::Level &slots : this->ParameterLayout().Levels()) {
    const int j = slots.group, la = slots.level;
    state.SetE(j, la, p[slots.energy]);
    std::vector<double> g;
    for (int ch = 1; ch <= (int)slots.amplitude.size(); ch++) {
      const double gamma = p[slots.amplitude[ch - 1]];
      state.SetGamma(j, la, ch, gamma);
      state.SetECGamma(j, la, ch, gamma);
      g.push_back(gamma);
    }
    state.SetSqrtNF(j, la, SqrtNFFactor(*this, this->GetJGroup(j), this->GetJGroup(j)->GetLevel(la), 1.0, g));
  }
  if (ecFinalBlocks_) this->ECFinalStatesFromFormal(state);
  if (captureMaskSetup_) this->FormalCaptureMasksFromParams(state);
}

/*!
 * This function is called for each iteration to calculate the shift
 * functions at new level energies when the Brune parametrization is used.
 *
 * In the Park formalism, and for a level coincident with another of its
 * J-group in every basis, the first and second energy derivatives at the level
 * energies are numerical (ChannelFunc::ShiftDerivative, CoulFunc::PEShift_d2E,
 * ShftFunc), four
 * and five more shift functions per channel and level, most of a serial
 * evaluation at new parameters.  They are taken after the shift functions, for
 * all levels at once and in parallel; the Coulomb functions are computed at
 * every call (physics/coulomb, no memo), so the values do not depend on the
 * number of threads.  Channels of a level that differ only in
 * the channel spin share one derivative, and the working copies of a fit share
 * them through the run's exact-argument memo (RunScratch::shiftDerivatives).
 */

void CompoundSystem::CalcShiftFunctions(const Config &configure, LevelState &state, ShiftOrders orders) const {
  const bool park = (configure.paramMask & Config::USE_PARK_FORMALISM);
  struct Derivative {
    int j, la;
    const PPair *pair;
    int ch, l, two_j;  ///< two_j: the channel's 2J key (PPair::SpinOrbitTwoJ)
    double radius;     ///< the channel's radius a_c (fm)
    double energy;
    bool first, second;
    double d1, d2;
    int source;  ///< an earlier entry of the same level, pair, l and radius (another channel spin), or -1
  };
  std::vector<Derivative> derivatives;
  struct Moved {
    const JGroup *group;
    int j, la;
    double energy;
    bool recompute, first, second;
  };
  std::vector<Moved> levels;
  for (int j = 1; j <= this->NumJGroups(); j++) {
    if (this->GetJGroup(j)->IsInRMatrix()) {
      const JGroup *theJGroup = this->GetJGroup(j);
      for (int la = 1; la <= theJGroup->NumLevels(); la++) {
        const ALevel *theLevel = theJGroup->GetLevel(la);
        if (theLevel->IsInRMatrix()) {
          const double levelEnergy = state.E(j, la);
          // The shift functions depend on the level energy alone, and a
          // minimizer moves one parameter at a time, so most levels have not
          // moved since the last call.  Park's J below depends on the widths
          // too and is always refreshed.
          const bool recompute = !state.ShiftsValidAt(j, la, levelEnergy);
          // A level on top of another of its J-group (azure2::physics::levels_coincide) takes the limits of the
          // divided differences between them (LevelPairShift), in S' and, for the gradient, S'': in every basis.
          bool coincident = false;
          for (int other = 1; other <= theJGroup->NumLevels() && !coincident; other++)
            coincident = other != la && theJGroup->GetLevel(other)->IsInRMatrix() &&
                         azure2::physics::levels_coincide(levelEnergy, state.E(j, other));
          const bool first = (park && recompute) || (coincident && !state.FirstDerivativesValidAt(j, la, levelEnergy));
          // d^2S/dE^2 (Park) enters only the gradient: the numerical second
          // derivative costs five shift functions per channel, so a chi-squared
          // alone leaves it out, and a later call that asks for it fills it in.
          // Coincident levels (rare) take it at once.
          const bool second = ((park && orders == ShiftOrders::kWithSecondDerivatives) || coincident) &&
                              !state.SecondDerivativesValidAt(j, la, levelEnergy);
          const int firstOfLevel = (int)derivatives.size();
          for (int ch = 1; ch <= theJGroup->NumChannels(); ch++) {
            const AChannel *theChannel = theJGroup->GetChannel(ch);
            const PPair *thePair = this->GetPair(theChannel->GetPairNum());
            if (thePair->GetPType() == 0) {
              if (!recompute && !first && !second) continue;
              int lValue = theChannel->GetL();
              const bool derivativeChannel = (park || coincident) && theChannel->GetRadType() == 'P';
              double resonanceEnergy = levelEnergy - (thePair->GetSepE() + thePair->GetExE());
              if (recompute) {
                // S as everywhere else (ChannelFunc: continuous through threshold).
                ChannelFunc channelFunc(thePair, theChannel->GetRadius());
                state.SetShift(j, la, ch, channelFunc.Shift(lValue, thePair->SpinOrbitTwoJ(theJGroup->GetJ()), resonanceEnergy));
              }
              if (derivativeChannel && (first || second)) {
                // Channels that differ only in the channel spin share the derivative (at one radius).
                int source = -1;
                for (int k = firstOfLevel; k < (int)derivatives.size() && source < 0; k++)
                  if (derivatives[k].source < 0 && derivatives[k].pair == thePair && derivatives[k].l == lValue && derivatives[k].two_j == thePair->SpinOrbitTwoJ(theJGroup->GetJ()) &&
                      derivatives[k].radius == theChannel->GetRadius())
                    source = k;
                derivatives.push_back({j, la, thePair, ch, lValue, thePair->SpinOrbitTwoJ(theJGroup->GetJ()), theChannel->GetRadius(),
                                       levelEnergy, first, second, 0.0, 0.0, source});
              }
            } else {
              state.SetShift(j, la, ch, state.Shift(j, 1, 1));
            }
          }
          levels.push_back({theJGroup, j, la, levelEnergy, recompute, first, second});
        }
      }
    }
  }

  // Park: dS/dE and d^2S/dE^2, in parallel (see above); a GSL error is
  // carried out of the parallel region and rethrown.
  const int numDerivatives = (int)derivatives.size();
  azure2::ShiftDerivativeMemo &memo = configure.scratch->shiftDerivatives;
  std::exception_ptr error;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) if (numDerivatives > 1 && !omp_in_parallel())
#endif
  for (int k = 0; k < numDerivatives; k++) {
    const azure2::support::ExtendedPrecisionScope fpPrecision;
    try {
      Derivative &d = derivatives[k];
      if (d.source >= 0) continue;
      const double threshold = d.pair->GetSepE() + d.pair->GetExE();
      const double resonanceEnergy = d.energy - threshold;
      CoulFunc theCoulombFunction(d.pair);
      const azure2::ShiftDerivativeMemo::Key key{d.pair->GetZ(1), d.pair->GetZ(2), d.l, d.two_j, azure2::physics::PotentialKey(theCoulombFunction.HybridPotential()),
                              theCoulombFunction.UsesTable(), d.pair->FunctionSettings().relativistic, d.pair->GetRedMass(),
                              d.pair->GetM(1) + d.pair->GetM(2), d.radius, threshold, d.energy};
      bool need1 = d.first, need2 = d.second;
      {
        std::lock_guard<std::mutex> lock(memo.mutex);
        auto it = memo.values.find(key);
        if (it != memo.values.end()) {
          if (need1 && it->second.has1) {
            d.d1 = it->second.d1;
            need1 = false;
          }
          if (need2 && it->second.has2) {
            d.d2 = it->second.d2;
            need2 = false;
          }
        }
      }
      if (!need1 && !need2) continue;
      // dS/dE as everywhere else (ChannelFunc, also BuildShiftDerivTable's
      // under Brune); d^2S/dE^2 from the branch the level is on.
      if (need1) {
        ChannelFunc channelFunc(d.pair, d.radius);
        d.d1 = channelFunc.ShiftDerivative(d.l, d.two_j, resonanceEnergy);
      }
      if (need2) {
        if (resonanceEnergy < 0.0) {
          ShftFunc theShiftFunction(d.pair, d.radius);
          d.d2 = theShiftFunction.EnergySecondDerivative(d.l, d.two_j, d.energy);
        } else {
          d.d2 = theCoulombFunction.PEShift_d2E(d.l, d.two_j, d.radius, resonanceEnergy);
        }
      }
      std::lock_guard<std::mutex> lock(memo.mutex);
      if (memo.values.size() >= azure2::ShiftDerivativeMemo::kMax) memo.values.clear();
      azure2::ShiftDerivativeMemo::Value &v = memo.values[key];
      if (need1) {
        v.has1 = true;
        v.d1 = d.d1;
      }
      if (need2) {
        v.has2 = true;
        v.d2 = d.d2;
      }
    } catch (...) {
#ifdef _OPENMP
#pragma omp critical(calc_shift_functions_error)
#endif
      if (!error) error = std::current_exception();
    }
  }
  if (error) std::rethrow_exception(error);
  for (Derivative &d : derivatives) {
    if (d.source >= 0) {
      d.d1 = derivatives[d.source].d1;
      d.d2 = derivatives[d.source].d2;
    }
    if (d.first) state.SetShiftDerivative(d.j, d.la, d.ch, d.d1);
    if (d.second) state.SetShiftSecondDerivative(d.j, d.la, d.ch, d.d2);
  }

  for (const Moved &m : levels) {
    const ALevel *theLevel = m.group->GetLevel(m.la);
    const JGroup *theJGroup = m.group;
    if (m.recompute) state.SetShiftsEnergy(m.j, m.la, m.energy);
    if (m.first) state.SetFirstDerivativesEnergy(m.j, m.la, m.energy);
    if (m.second) state.SetSecondDerivativesEnergy(m.j, m.la, m.energy);
    if (park) {
      // Park (2021) Eq. (28): J = 1 - sum_c gamma_c^2 (dS_c/dE) at E_lambda.
      double parkNorm = 1.0;
      for (int ch = 1; ch <= theJGroup->NumChannels(); ch++) {
        const AChannel *theChannel = theJGroup->GetChannel(ch);
        if (theChannel->GetRadType() == 'P' && this->GetPair(theChannel->GetPairNum())->GetPType() == 0)
          parkNorm -= pow(state.Gamma(m.j, m.la, ch), 2.0) * state.ShiftDerivative(m.j, m.la, ch);
      }
      state.SetParkNorm(m.j, m.la, parkNorm);
      // N_f^{1/2} multiplies the level's reduced width amplitude wherever
      // it serves as an external-capture final state.  It is defined with
      // Brune's amplitudes, gamma_Park / sqrt(J); folding the 1/sqrt(J)
      // in here keeps gamma * N_f^{1/2} the same number in both
      // parametrizations.
      std::vector<double> g;
      for (int ch = 1; ch <= theJGroup->NumChannels(); ch++) g.push_back(state.Gamma(m.j, m.la, ch));
      state.SetSqrtNF(m.j, m.la, SqrtNFFactor(*this, theJGroup, theLevel, parkNorm, g));
    }
  }
}

/*!
 * The fit energies and amplitudes of `from`, nothing derived from them.
 */

void CompoundSystem::CopyFitEnergiesAndAmplitudes(const LevelState &from, LevelState &to) const {
  for (int j = 1; j <= this->NumJGroups(); j++)
    for (int la = 1; la <= this->GetJGroup(j)->NumLevels(); la++) {
      to.SetE(j, la, from.E(j, la));
      for (int ch = 1; ch <= this->GetJGroup(j)->NumChannels(); ch++) to.SetGamma(j, la, ch, from.Gamma(j, la, ch));
    }
}
