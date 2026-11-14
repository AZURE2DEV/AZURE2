#include <cmath>
#include <exception>
#include <mutex>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "CNuc.h"
#include "ChannelFunc.h"
#include "CNucExternalCapture.h"
#include "Config.h"
#include "CoulFunc.h"
#include "RunScratch.h"
#include "ShftFunc.h"

// The members of CNuc that write its LevelState at a parameter vector: the fit values, the shift functions and what
// follows from them.

/*!
 * Fills the level state from a parameter vector (FillMnParams' order), with everything derived from it but the
 * shift functions (CalcShiftFunctions).
 */

void CNuc::FillCompoundFromParams(const vector_r &params) {
  const vector_r &p = gammaDecays_->Empty() ? params : gammaDecays_->ToEngine(*this, params);  // GammaDecayMap.h
  int i = 0;
  for (int j = 1; j <= this->NumJGroups(); j++) {
    for (int la = 1; la <= this->Group(j)->NumLevels(); la++) {
      ALevel *level = this->Group(j)->GetLevel(la);
      levelState_.SetE(j, la, p[i]);
      i++;
      std::vector<double> g;
      for (int ch = 1; ch <= this->Group(j)->NumChannels(); ch++) {
        levelState_.SetGamma(j, la, ch, p[i]);
        level->SetECGamma(ch, p[i]);
        g.push_back(p[i]);
        i++;
      }
      level->SetSqrtNFFactor(SqrtNFFactor(*this, this->Group(j), level, 1.0, g));
    }
  }
  if (ecFinalBlocks_) this->ECFinalStatesFromFormal();
  if (captureMaskSetup_) this->FormalCaptureMasksFromParams();
}

/*!
 * This function is called for each iteration to calculate the shift
 * functions at new level energies when the Brune parametrization is used.
 *
 * In the Park formalism the first and second energy derivatives at the level
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

void CNuc::CalcShiftFunctions(const Config &configure, bool secondDerivatives) {
  const bool park = (configure.paramMask & Config::USE_PARK_FORMALISM);
  struct Derivative {
    int j, la;
    PPair *pair;
    int ch, l, two_j;  ///< two_j: the channel's 2J key (PPair::SpinOrbitTwoJ)
    double energy;
    bool first, second;
    double d1, d2;
    int source;  ///< an earlier entry of the same level, pair and l (another channel spin), or -1
  };
  std::vector<Derivative> derivatives;
  struct Moved {
    JGroup *group;
    int j, la;
    double energy;
    bool recompute, second;
  };
  std::vector<Moved> levels;
  for (int j = 1; j <= this->NumJGroups(); j++) {
    if (this->Group(j)->IsInRMatrix()) {
      JGroup *theJGroup = this->Group(j);
      for (int la = 1; la <= theJGroup->NumLevels(); la++) {
        ALevel *theLevel = theJGroup->GetLevel(la);
        if (theLevel->IsInRMatrix()) {
          const double levelEnergy = levelState_.E(j, la);
          // The shift functions depend on the level energy alone, and a
          // minimizer moves one parameter at a time, so most levels have not
          // moved since the last call.  Park's J below depends on the widths
          // too and is always refreshed.
          const bool recompute = !levelState_.ShiftsValidAt(j, la, levelEnergy);
          // d^2S/dE^2 (Park) enters only the gradient: the numerical second
          // derivative costs five shift functions per channel, so a chi-squared
          // alone leaves it out, and a later call that asks for it fills it in.
          const bool second = park && secondDerivatives && !levelState_.SecondDerivativesValidAt(j, la, levelEnergy);
          const int firstOfLevel = (int)derivatives.size();
          for (int ch = 1; ch <= theJGroup->NumChannels(); ch++) {
            AChannel *theChannel = theJGroup->GetChannel(ch);
            PPair *thePair = this->Pair(theChannel->GetPairNum());
            if (thePair->GetPType() == 0) {
              if (!recompute && !second) continue;
              int lValue = theChannel->GetL();
              const bool parkP = park && theChannel->GetRadType() == 'P';
              double resonanceEnergy = levelEnergy - (thePair->GetSepE() + thePair->GetExE());
              if (recompute) {
                // S as everywhere else (ChannelFunc: continuous through threshold).
                ChannelFunc channelFunc(thePair, !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC));
                levelState_.SetShift(j, la, ch, channelFunc.Shift(lValue, thePair->SpinOrbitTwoJ(theJGroup->GetJ()), resonanceEnergy));
              }
              if (parkP && (recompute || second)) {
                // Channels that differ only in the channel spin share the derivative.
                int source = -1;
                for (int k = firstOfLevel; k < (int)derivatives.size() && source < 0; k++)
                  if (derivatives[k].source < 0 && derivatives[k].pair == thePair && derivatives[k].l == lValue && derivatives[k].two_j == thePair->SpinOrbitTwoJ(theJGroup->GetJ())) source = k;
                derivatives.push_back({j, la, thePair, ch, lValue, thePair->SpinOrbitTwoJ(theJGroup->GetJ()), levelEnergy, recompute, second, 0.0, 0.0, source});
              }
            } else {
              levelState_.SetShift(j, la, ch, levelState_.Shift(j, 1, 1));
            }
          }
          levels.push_back({theJGroup, j, la, levelEnergy, recompute, second});
        }
      }
    }
  }

  // Park: dS/dE and d^2S/dE^2, in parallel (see above); a GSL error is
  // carried out of the parallel region and rethrown.
  const bool gslCoulomb = !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC);
  const int numDerivatives = (int)derivatives.size();
  azure2::ShiftDerivativeMemo &memo = configure.scratch->shiftDerivatives;
  std::exception_ptr error;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) if (numDerivatives > 1 && !omp_in_parallel())
#endif
  for (int k = 0; k < numDerivatives; k++) {
    try {
      Derivative &d = derivatives[k];
      if (d.source >= 0) continue;
      const double threshold = d.pair->GetSepE() + d.pair->GetExE();
      const double resonanceEnergy = d.energy - threshold;
      CoulFunc theCoulombFunction(d.pair, gslCoulomb);
      const azure2::ShiftDerivativeMemo::Key key{d.pair->GetZ(1), d.pair->GetZ(2), d.l, d.two_j, gslCoulomb, azure2::physics::PotentialKey(theCoulombFunction.HybridPotential()),
                              theCoulombFunction.UsesTable(), d.pair->FunctionSettings().relativistic, d.pair->GetRedMass(),
                              d.pair->GetM(1) + d.pair->GetM(2), d.pair->GetChRad(), threshold, d.energy};
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
        ChannelFunc channelFunc(d.pair, gslCoulomb);
        d.d1 = channelFunc.ShiftDerivative(d.l, d.two_j, resonanceEnergy);
      }
      if (need2) {
        if (resonanceEnergy < 0.0) {
          ShftFunc theShiftFunction(d.pair);
          d.d2 = theShiftFunction.EnergySecondDerivative(d.l, d.two_j, d.energy);
        } else {
          d.d2 = theCoulombFunction.PEShift_d2E(d.l, d.two_j, d.pair->GetChRad(), resonanceEnergy);
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
    if (d.first) levelState_.SetShiftDerivative(d.j, d.la, d.ch, d.d1);
    if (d.second) levelState_.SetShiftSecondDerivative(d.j, d.la, d.ch, d.d2);
  }

  for (const Moved &m : levels) {
    ALevel *theLevel = m.group->GetLevel(m.la);
    JGroup *theJGroup = m.group;
    if (m.recompute) levelState_.SetShiftsEnergy(m.j, m.la, m.energy);
    if (m.second) levelState_.SetSecondDerivativesEnergy(m.j, m.la, m.energy);
    if (park) {
      // Park (2021) Eq. (28): J = 1 - sum_c gamma_c^2 (dS_c/dE) at E_lambda.
      double parkNorm = 1.0;
      for (int ch = 1; ch <= theJGroup->NumChannels(); ch++) {
        AChannel *theChannel = theJGroup->GetChannel(ch);
        if (theChannel->GetRadType() == 'P' && this->Pair(theChannel->GetPairNum())->GetPType() == 0)
          parkNorm -= pow(levelState_.Gamma(m.j, m.la, ch), 2.0) * levelState_.ShiftDerivative(m.j, m.la, ch);
      }
      levelState_.SetParkNorm(m.j, m.la, parkNorm);
      // N_f^{1/2} multiplies the level's reduced width amplitude wherever
      // it serves as an external-capture final state.  It is defined with
      // Brune's amplitudes, gamma_Park / sqrt(J); folding the 1/sqrt(J)
      // in here keeps gamma * N_f^{1/2} the same number in both
      // parametrizations.
      std::vector<double> g;
      for (int ch = 1; ch <= theJGroup->NumChannels(); ch++) g.push_back(levelState_.Gamma(m.j, m.la, ch));
      theLevel->SetSqrtNFFactor(SqrtNFFactor(*this, theJGroup, theLevel, parkNorm, g));
    }
  }
}

/*!
 * The fit energies and amplitudes of `from`, nothing derived from them.
 */

void CNuc::CopyFitEnergiesAndAmplitudes(const LevelState &from) {
  for (int j = 1; j <= this->NumJGroups(); j++)
    for (int la = 1; la <= this->Group(j)->NumLevels(); la++) {
      levelState_.SetE(j, la, from.E(j, la));
      for (int ch = 1; ch <= this->Group(j)->NumChannels(); ch++) levelState_.SetGamma(j, la, ch, from.Gamma(j, la, ch));
    }
}
