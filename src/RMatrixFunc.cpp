#include "EPoint.h"
#include "LevelCoincidence.h"
#include "CNuc.h"
#include "MatrixInv.h"
#include "RMatrixFunc.h"
#include <assert.h>
#include <cmath>
#include <limits>

namespace {

/// A level closer than this to the energy (MeV) puts the evaluation at a pole of the R-matrix.
const double kRPoleWindow = 1.0e-6;

/*!
 * \f$ [1-RL]^{-1}R \f$ of a group whose R-matrix \f$ R = \gamma^T M^{-1} \gamma \f$ has a pole at this energy,
 * in the level form \f$ \gamma^T (M - \gamma L \gamma^T)^{-1} \gamma \f$ (push-through identity), which stays finite
 * there; also \f$ [1-RL]^{-1} = 1 + [1-RL]^{-1} R L \f$.  M is the matrix whose inverse gives R (diagonal
 * \f$ E_\lambda - E \f$ in the standard formalism, Brune's real level matrix otherwise), g[level][channel] the
 * reduced width amplitudes of the levels kept in M, L the channels' \f$ L \f$.
 */
void PoleLevelForm(const matrix_c &M, const std::vector<std::vector<double>> &g, const vector_c &L,
                   matrix_c &rlInvR, matrix_c &rlInv) {
  const size_t nl = M.size(), nc = L.size();
  matrix_c aInv(M);
  for (size_t la = 0; la < nl; la++)
    for (size_t lap = 0; lap < nl; lap++)
      for (size_t c = 0; c < nc; c++) aInv[la][lap] -= g[la][c] * L[c] * g[lap][c];
  MatrixInv a(aInv);
  rlInvR.assign(nc, vector_c(nc, complex(0.0, 0.0)));
  rlInv.assign(nc, vector_c(nc, complex(0.0, 0.0)));
  for (size_t c = 0; c < nc; c++)
    for (size_t cp = 0; cp < nc; cp++) {
      complex x(0.0, 0.0);
      for (size_t la = 0; la < nl; la++)
        for (size_t lap = 0; lap < nl; lap++) x += g[la][c] * a.inverse()[la][lap] * g[lap][cp];
      rlInvR[c][cp] = x;
    }
  for (size_t c = 0; c < nc; c++)
    for (size_t cp = 0; cp < nc; cp++) rlInv[c][cp] = (c == cp ? 1.0 : 0.0) + rlInvR[c][cp] * L[cp];
}

}  // namespace

/*!
 * The RMatrixFunc object is created with reference to a CNuc object.
 */

RMatrixFunc::RMatrixFunc(CNuc *compound, const Config &configure) :
  compound_(compound),
  configure_(configure) {}

/*!
 * Returns an R-Matrix element specified by positions in the JGroup and AChannel vectors.
 */

complex RMatrixFunc::GetRMatrixElement(int jGroupNum, int channelNum, int channelPrimeNum) const {
  return r_matrices_[jGroupNum - 1][channelNum - 1][channelPrimeNum - 1];
}

/*!
 * Returns an \f$ [1-RL] \f$ Matrix element specified by positions in the JGroup and AChannel vectors.
 */

complex RMatrixFunc::GetRLMatrixElement(int jGroupNum, int channelNum, int channelPrimeNum) const {
  return rl_matrices_[jGroupNum - 1][channelNum - 1][channelPrimeNum - 1];
}

/*!
 * Returns a \f$ [1-RL]^{-1} \f$ Matrix element specified by positions in the JGroup and AChannel vectors.
 */

complex RMatrixFunc::GetRLInvMatrixElement(int jGroupNum, int channelNum, int channelPrimeNum) const {
  return rl_inv_matrices_[jGroupNum - 1][channelNum - 1][channelPrimeNum - 1];
}

/*!
 * Returns a \f$ [1-RL]^{-1}R \f$ Matrix element specified by positions in the JGroup and AChannel vectors.
 */

complex RMatrixFunc::GetRLInvRMatrixElement(int jGroupNum, int channelNum, int channelPrimeNum) const {
  return rl_inv_r_matrices_[jGroupNum - 1][channelNum - 1][channelPrimeNum - 1];
}

/*!
 * Returns an entire \f$ [1-RL] \f$ Matrix specified by a position in the JGroup vector.
 */

matrix_c *RMatrixFunc::GetJSpecRLMatrix(int jGroupNum) {
  matrix_c *b = &rl_matrices_[jGroupNum - 1];
  return b;
}

/*!
 * Clears all matrices associated with the RMatrixFunc object.
 */

void RMatrixFunc::ClearMatrices() {
  r_matrices_.clear();
  rl_matrices_.clear();
  rl_inv_matrices_.clear();
  rl_inv_r_matrices_.clear();
  tmatrix_.clear();
  ec_tmatrix_.clear();
  pole_rl_inv_r_.clear();
  pole_rl_inv_.clear();
  InvalidateTMatrixProducts();
}

/*!
 * This function creates the \f$ [1-RL] \f$ and \f$ R \f$ Matrices from the CNuc object.
 */

void RMatrixFunc::FillMatrices(EPoint *point) {
  double inEnergy;
  if (compound()->GetPair(compound()->GetPairNumFromKey(point->GetEntranceKey()))->GetPType() == 20)
    inEnergy = point->GetCMEnergy() +
        compound()->GetPair(compound()->GetPairNumFromKey(point->GetExitKey()))->GetSepE() +
        compound()->GetPair(compound()->GetPairNumFromKey(point->GetExitKey()))->GetExE();
  else
    inEnergy = point->GetCMEnergy() +
        compound()->GetPair(compound()->GetPairNumFromKey(point->GetEntranceKey()))->GetSepE() +
        compound()->GetPair(compound()->GetPairNumFromKey(point->GetEntranceKey()))->GetExE();
  // Groups of another compound nucleus are not part of this point's matrices
  // (EPoint::CalcEDependentValues left their rows empty).
  activeNucleus_ = compound()->GetPair(compound()->GetPairNumFromKey(point->GetEntranceKey()))->GetNucleus();
  for (int j = 1; j <= compound()->NumJGroups(); j++) {
    if (compound()->GetJGroup(j)->IsInRMatrix() && compound()->GetJGroup(j)->GetNucleus() == activeNucleus_) {
      JGroup *theJGroup = compound()->GetJGroup(j);
      // The R-matrix has a pole at the energy of every level with a width (standard formalism: the formal energy;
      // Brune: the observed one, where the level's row of the real level matrix vanishes).  An energy on or next to
      // one (the adaptive sub-point grid of a target effect places a node exactly there) is solved in the level
      // form, which is finite there.  Before 2026-11-01 the standard formalism gave NaN there and Brune's an
      // ill-conditioned inverse (tests/13N capture with target integration: NaN above the 1/2+ resonance, 1.2 %
      // off at it).
      std::vector<int> coupled;  // levels with a width, in R-matrix order
      bool atPole = false;
      for (int la = 1; la <= theJGroup->NumLevels(); la++) {
        if (!theJGroup->GetLevel(la)->IsInRMatrix()) continue;
        bool hasWidth = false;
        for (int ch = 1; ch <= theJGroup->NumChannels() && !hasWidth; ch++)
          hasWidth = (theJGroup->GetLevel(la)->GetFitGamma(ch) != 0.0);
        if (!hasWidth) continue;
        coupled.push_back(la);
        if (std::fabs(theJGroup->GetLevel(la)->GetFitE() - inEnergy) < kRPoleWindow) atPole = true;
      }
      matrix_c poleM;  // the matrix whose inverse gives R, over the coupled levels (filled only at a pole)
      if (configure().paramMask & Config::USE_BRUNE_FORMALISM) {
        matrix_c qMatrixInverse;
        for (int la = 1; la <= compound()->GetJGroup(j)->NumLevels(); la++) {
          if (!compound()->GetJGroup(j)->GetLevel(la)->IsInRMatrix()) continue;
          ALevel *level = compound()->GetJGroup(j)->GetLevel(la);
          vector_c tempVector;
          for (int lap = 1; lap <= compound()->GetJGroup(j)->NumLevels(); lap++) {
            if (!compound()->GetJGroup(j)->GetLevel(lap)->IsInRMatrix()) continue;
            ALevel *levelp = compound()->GetJGroup(j)->GetLevel(lap);
            complex sum = (la == lap) ? complex(level->GetFitE() - inEnergy, 0.0) : complex(0.0, 0.0);
            for (int ch = 1; ch <= compound()->GetJGroup(j)->NumChannels(); ch++) {
              if (compound()->GetJGroup(j)->GetChannel(ch)->GetRadType() != 'P') continue;
              double gammaCh = level->GetFitGamma(ch);
              double gammaChp = levelp->GetFitGamma(ch);
              complex channelShift = point->GetLoElement(j, ch) + complex(compound()->GetJGroup(j)->GetChannel(ch)->GetBoundaryCondition(), -1. * pow(point->GetSqrtPenetrability(j, ch), 2.));
              sum -= gammaCh * gammaChp * channelShift;
              if (la == lap)
                sum += gammaCh * gammaChp * level->GetShiftFunction(ch);
              else if (LevelsCoincide(level->GetFitE(), levelp->GetFitE()))
                sum += gammaCh * gammaChp * (level->GetShiftFunction(ch) + level->GetShiftDerivative(ch) * (inEnergy - level->GetFitE()));
              else
                sum += gammaCh * gammaChp *
                    (level->GetShiftFunction(ch) * (inEnergy - levelp->GetFitE()) - levelp->GetShiftFunction(ch) * (inEnergy - level->GetFitE())) /
                    (level->GetFitE() - levelp->GetFitE());
            }
            // A level without widths (zero row and column) at E = E_lambda: as in
            // AMatrixFunc::FillLevelMatrix, a unit diagonal keeps the matrix regular;
            // the level contributes nothing either way.
            if (la == lap && sum == complex(0.0, 0.0)) {
              bool decoupled = true;
              for (int ch = 1; ch <= compound()->GetJGroup(j)->NumChannels() && decoupled; ch++)
                decoupled = (level->GetFitGamma(ch) == 0.0);
              if (decoupled) sum = 1.0;
            }
            tempVector.push_back(sum);
          }
          qMatrixInverse.push_back(tempVector);
        }
        // qMatrixInverse's rows and columns follow the in-R-matrix levels only (a bound state kept just as an
        // external-capture final state is not one of them): row of level la.
        std::vector<int> rowOf;
        for (int la = 1, r = 0; la <= theJGroup->NumLevels(); la++)
          rowOf.push_back(theJGroup->GetLevel(la)->IsInRMatrix() ? r++ : -1);
        if (atPole) {
          // Keep the coupled levels (a level without widths has no off-diagonal elements, so leaving it out is
          // exact).  R itself is not finite here and is not formed: Q^{-1} is singular at an exact pole.
          for (int la : coupled) {
            vector_c row;
            for (int lap : coupled) row.push_back(qMatrixInverse[rowOf[la - 1]][rowOf[lap - 1]]);
            poleM.push_back(row);
          }
        }
        // A hard-sphere group has no levels: its level matrix is empty and its R-matrix is zero.
        const matrix_c unit(1, vector_c(1, complex(1.0, 0.0)));
        MatrixInv matrixInv(qMatrixInverse.empty() || atPole ? unit : qMatrixInverse);
        const double notFormed = std::numeric_limits<double>::quiet_NaN();
        for (int ch = 1; ch <= compound()->GetJGroup(j)->NumChannels(); ch++) {
          for (int chp = 1; chp <= compound()->GetJGroup(j)->NumChannels(); chp++) {
            complex rTemp(0.0, 0.0);
            for (int la = 1; la <= compound()->GetJGroup(j)->NumLevels() && !atPole; la++) {
              if (!compound()->GetJGroup(j)->GetLevel(la)->IsInRMatrix()) continue;
              complex temp(0.0, 0.0);
              for (int lap = 1; lap <= compound()->GetJGroup(j)->NumLevels(); lap++) {
                if (!compound()->GetJGroup(j)->GetLevel(lap)->IsInRMatrix()) continue;
                // Rows by in-R-matrix position: before 2026-11-01 by level number, out of range (undefined) when a
                // capture final state came before an R-matrix level of the group.
                temp += matrixInv.inverse()[rowOf[la - 1]][rowOf[lap - 1]] *
                        compound()->GetJGroup(j)->GetLevel(lap)->GetFitGamma(chp);
              }
              rTemp += temp * compound()->GetJGroup(j)->GetLevel(la)->GetFitGamma(ch);
            }
            if (atPole) rTemp = complex(notFormed, notFormed);  // not finite; the level form is used instead
            this->AddRMatrixElement(j, ch, chp, rTemp);
            double tempPene = 0.;
            if (compound()->GetJGroup(j)->GetChannel(chp)->GetRadType() == 'P')
              tempPene = pow(point->GetSqrtPenetrability(j, chp), 2.);
            if (ch == chp)
              this->AddRLMatrixElement(j, ch, chp, 1. - complex(0., 1.) * rTemp * tempPene);
            else
              this->AddRLMatrixElement(j, ch, chp, -complex(0., 1.) * rTemp * tempPene);
          }
        }
      } else {
        if (atPole)
          for (size_t i = 0; i < coupled.size(); i++) {
            ALevel *level = theJGroup->GetLevel(coupled[i]);
            double gammaSum = 0.;
            if (configure().paramMask & Config::USE_RMC_FORMALISM)
              for (int chpp = 1; chpp <= theJGroup->NumChannels(); chpp++)
                if (theJGroup->GetChannel(chpp)->GetRadType() == 'M' || theJGroup->GetChannel(chpp)->GetRadType() == 'E')
                  gammaSum += pow(level->GetFitGamma(chpp), 2.0);
            vector_c row(coupled.size(), complex(0.0, 0.0));
            row[i] = level->GetFitE() - inEnergy - complex(0.0, 1.0) * gammaSum;
            poleM.push_back(row);
          }
        for (int ch = 1; ch <= compound()->GetJGroup(j)->NumChannels(); ch++) {
          for (int chp = 1; chp <= compound()->GetJGroup(j)->NumChannels(); chp++) {
            complex sum(0.0, 0.0);
            for (int la = 1; la <= compound()->GetJGroup(j)->NumLevels(); la++) {
              if (compound()->GetJGroup(j)->GetLevel(la)->IsInRMatrix()) {
                ALevel *level = compound()->GetJGroup(j)->GetLevel(la);
                double gammaCh = level->GetFitGamma(ch);
                double gammaChp = level->GetFitGamma(chp);
                // A zero width contributes zero, also at E = E_lambda (not 0/0).
                if (gammaCh == 0.0 || gammaChp == 0.0) continue;
                double resenergy = level->GetFitE();
                double gammaSum = 0.;
                if (configure().paramMask & Config::USE_RMC_FORMALISM)
                  for (int chpp = 1; chpp <= compound()->GetJGroup(j)->NumChannels(); chpp++)
                    if (compound()->GetJGroup(j)->GetChannel(chpp)->GetRadType() == 'M' ||
                        compound()->GetJGroup(j)->GetChannel(chpp)->GetRadType() == 'E')
                      gammaSum += pow(compound()->GetJGroup(j)->GetLevel(la)->GetFitGamma(chpp), 2.0);
                sum += gammaCh * gammaChp / (resenergy - inEnergy - complex(0.0, 1.0) * gammaSum);
              }
            }
            this->AddRMatrixElement(j, ch, chp, sum);
            complex loElement = point->GetLoElement(j, chp);
            if (ch == chp)
              this->AddRLMatrixElement(j, ch, chp, 1.0 - sum * loElement);
            else
              this->AddRLMatrixElement(j, ch, chp, -sum * loElement);
          }
        }
      }
      if (atPole) {
        // L as each branch puts it into [1-RL]: the channel's L^o in the standard formalism; i P of a particle
        // channel (0 otherwise) in Brune's, whose R already holds the shift.
        const bool brune = configure().paramMask & Config::USE_BRUNE_FORMALISM;
        vector_c L;
        for (int ch = 1; ch <= theJGroup->NumChannels(); ch++) {
          if (!brune)
            L.push_back(point->GetLoElement(j, ch));
          else if (theJGroup->GetChannel(ch)->GetRadType() == 'P')
            L.push_back(complex(0.0, pow(point->GetSqrtPenetrability(j, ch), 2.)));
          else
            L.push_back(complex(0.0, 0.0));
        }
        std::vector<std::vector<double>> g;
        for (int la : coupled) {
          std::vector<double> row;
          for (int ch = 1; ch <= theJGroup->NumChannels(); ch++) row.push_back(theJGroup->GetLevel(la)->GetFitGamma(ch));
          g.push_back(row);
        }
        if (pole_rl_inv_r_.size() < (size_t)j) {
          pole_rl_inv_r_.resize(j);
          pole_rl_inv_.resize(j);
        }
        PoleLevelForm(poleM, g, L, pole_rl_inv_r_[j - 1], pole_rl_inv_[j - 1]);
      }
    }
  }
}


/*!
 * This function inverts the \f$ [1-RL] \f$ matrices and creates the \f$ [1-RL]^{-1}R \f$ matrices.
 */

void RMatrixFunc::InvertMatrices() {
  for (int j = 1; j <= compound()->NumJGroups(); j++) {
    if (compound()->GetJGroup(j)->IsInRMatrix() && compound()->GetJGroup(j)->GetNucleus() != activeNucleus_) {
      // Another nucleus's group: keep the per-group vectors aligned with an empty entry.
      this->AddRLInvMatrix(matrix_c());
      continue;
    }
    if (compound()->GetJGroup(j)->IsInRMatrix() && (size_t)j <= pole_rl_inv_r_.size() && !pole_rl_inv_r_[j - 1].empty()) {
      // At a pole of R (FillMatrices): [1-RL] is not finite; its inverse and [1-RL]^{-1}R come from the level form.
      this->AddRLInvMatrix(pole_rl_inv_[j - 1]);
      for (int ch = 1; ch <= compound()->GetJGroup(j)->NumChannels(); ch++)
        for (int chp = 1; chp <= compound()->GetJGroup(j)->NumChannels(); chp++)
          this->AddRLInvRMatrixElement(j, ch, chp, pole_rl_inv_r_[j - 1][ch - 1][chp - 1]);
      continue;
    }
    if (compound()->GetJGroup(j)->IsInRMatrix()) {
      matrix_c *theRLMatrix = this->GetJSpecRLMatrix(j);
      MatrixInv matrixInv(*theRLMatrix);
      this->AddRLInvMatrix(matrixInv.inverse());
      for (int ch = 1; ch <= compound()->GetJGroup(j)->NumChannels(); ch++) {
        for (int chp = 1; chp <= compound()->GetJGroup(j)->NumChannels(); chp++) {
          complex rlinvrElement(0., 0.);
          for (int chpp = 1; chpp <= compound()->GetJGroup(j)->NumChannels(); chpp++) {
            rlinvrElement += this->GetRLInvMatrixElement(j, ch, chpp) *
                this->GetRMatrixElement(j, chpp, chp);
          }
          this->AddRLInvRMatrixElement(j, ch, chp, rlinvrElement);
        }
      }
    }
  }
}

/*!
 * This function calculates the T-Matrix for each reaction pathways based on the \f$ [1-RL]^{-1}R \f$ matrix.
 */

void RMatrixFunc::CalculateTMatrix(EPoint *point) {
  int aa = compound()->GetPairNumFromKey(point->GetEntranceKey());
  int irEnd;
  int irStart;
  bool isRMC = false;
  if ((configure().paramMask & Config::USE_RMC_FORMALISM) &&
      compound()->GetPair(compound()->GetPairNumFromKey(point->GetExitKey()))->GetPType() == 10) {
    irStart = 1;
    irEnd = compound()->GetPair(aa)->NumDecays();
    isRMC = true;
  } else {
    irStart = 0;
    while (irStart < compound()->GetPair(aa)->NumDecays()) {
      irStart++;
      if (compound()->GetPair(aa)->GetDecay(irStart)->GetPairNum() ==
          compound()->GetPairNumFromKey(point->GetExitKey())) break;
    }
    irEnd = irStart;
  }
  for (int ir = irStart; ir <= irEnd; ir++) {
    Decay *theDecay = compound()->GetPair(aa)->GetDecay(ir);
    for (int k = 1; k <= theDecay->NumKGroups(); k++) {
      for (int m = 1; m <= theDecay->GetKGroup(k)->NumMGroups(); m++) {
        MGroup *theMGroup = theDecay->GetKGroup(k)->GetMGroup(m);
        JGroup *theJGroup = compound()->GetJGroup(theMGroup->GetJNum());
        AChannel *entranceChannel = theJGroup->GetChannel(theMGroup->GetChNum());
        AChannel *exitChannel = theJGroup->GetChannel(theMGroup->GetChpNum());
        complex uphase = point->GetExpCoulombPhase(theMGroup->GetJNum(), theMGroup->GetChNum()) *
            point->GetExpHardSpherePhase(theMGroup->GetJNum(), theMGroup->GetChNum()) *
            point->GetExpCoulombPhase(theMGroup->GetJNum(), theMGroup->GetChpNum()) *
            point->GetExpHardSpherePhase(theMGroup->GetJNum(), theMGroup->GetChpNum());
        complex umatrix = 2.0 * std::complex<double>(0.0, 1.0) *
            point->GetSqrtPenetrability(theMGroup->GetJNum(), theMGroup->GetChNum()) *
            point->GetSqrtPenetrability(theMGroup->GetJNum(), theMGroup->GetChpNum()) *
            this->GetRLInvRMatrixElement(theMGroup->GetJNum(),
                                         theMGroup->GetChpNum(),
                                         theMGroup->GetChNum());
        complex tphase = point->GetExpCoulombPhase(theMGroup->GetJNum(), theMGroup->GetChNum()) *
            point->GetExpCoulombPhase(theMGroup->GetJNum(), theMGroup->GetChNum());
        complex tmatrix;
        if (isRMC)
          this->AddTMatrixElement(k, m, complex(0.0, -1.0) * umatrix, ir);
        else {
          if (theMGroup->GetChNum() == theMGroup->GetChpNum()) {
            tmatrix = tphase - uphase * (1.0 + umatrix);
          } else
            tmatrix = -uphase * umatrix;
          this->AddTMatrixElement(k, m, tmatrix);
        }
      }
      for (int m = 1; m <= theDecay->GetKGroup(k)->NumECMGroups(); m++) {
        ECMGroup *theECMGroup = theDecay->GetKGroup(k)->GetECMGroup(m);
        ALevel *finalLevel = compound()->GetJGroup(theECMGroup->GetJGroupNum())->GetLevel(theECMGroup->GetLevelNum());
        double ecNormParam = finalLevel->GetFitGamma(theECMGroup->GetFinalChannel()) *
            finalLevel->GetSqrtNFFactor() * finalLevel->GetECConversionFactor(theECMGroup->GetFinalChannel());
        // Use energy-shift aware EC amplitude calculation if energy shifts are active
        complex ecAmplitude;
        if (configure().paramMask & Config::USE_EXTERNAL_CAPTURE) {
          // Use the new method that accounts for energy shifts through interpolation
          ecAmplitude = point->GetECAmplitudeWithShift(k, m, compound(), configure());
        } else {
          // Use cached amplitude for backward compatibility
          ecAmplitude = point->GetECAmplitude(k, m);
        }
        complex tmatrix = ecNormParam * ecAmplitude;
        if (theECMGroup->IsChannelCapture()) {
          MGroup *chanMGroup = compound()->GetPair(aa)->GetDecay(theECMGroup->GetChanCapDecay())->GetKGroup(theECMGroup->GetChanCapKGroup())->GetMGroup(theECMGroup->GetChanCapMGroup());
          AChannel *chanEntranceChannel = compound()->GetJGroup(chanMGroup->GetJNum())->GetChannel(chanMGroup->GetChNum());
          AChannel *chanExitChannel = compound()->GetJGroup(chanMGroup->GetJNum())->GetChannel(chanMGroup->GetChpNum());
          complex umatrix = 2.0 * std::complex<double>(0.0, 1.0) *
              point->GetSqrtPenetrability(chanMGroup->GetJNum(), chanMGroup->GetChNum()) *
              this->GetRLInvRMatrixElement(chanMGroup->GetJNum(),
                                           chanMGroup->GetChpNum(),
                                           chanMGroup->GetChNum());
          tmatrix = tmatrix * umatrix;
        }
        this->AddECTMatrixElement(k, m, tmatrix);
      }
    }
  }
}

/*!
 * This function adds an R-Matrix element specified by positions in the JGroup and AChannel vectors.
 */

void RMatrixFunc::AddRMatrixElement(int jGroupNum, int channelNum, int channelPrimeNum, complex matrixElement) {
  matrix_c e;
  vector_c f;
  while (jGroupNum > r_matrices_.size()) r_matrices_.push_back(e);
  while (channelNum > r_matrices_[jGroupNum - 1].size()) r_matrices_[jGroupNum - 1].push_back(f);
  r_matrices_[jGroupNum - 1][channelNum - 1].push_back(matrixElement);
  assert(channelPrimeNum = r_matrices_[jGroupNum - 1][channelNum - 1].size());
}

/*!
 * This function adds a \f$ [1-RL] \f$ matrix element specified by positions in the JGroup and AChannel vectors.
 */

void RMatrixFunc::AddRLMatrixElement(int jGroupNum, int channelNum, int channelPrimeNum, complex matrixElement) {
  matrix_c e;
  vector_c f;
  while (jGroupNum > rl_matrices_.size()) rl_matrices_.push_back(e);
  while (channelNum > rl_matrices_[jGroupNum - 1].size()) rl_matrices_[jGroupNum - 1].push_back(f);
  rl_matrices_[jGroupNum - 1][channelNum - 1].push_back(matrixElement);
  assert(channelPrimeNum = rl_matrices_[jGroupNum - 1][channelNum - 1].size());
}

/*!
 * This function adds a \f$ [1-RL]^{-1}R \f$ matrix element specified by positions in the JGroup and AChannel vectors.
 */

void RMatrixFunc::AddRLInvRMatrixElement(int jGroupNum, int channelNum, int channelPrimeNum, complex matrixElement) {
  matrix_c e;
  vector_c f;
  while (jGroupNum > rl_inv_r_matrices_.size()) rl_inv_r_matrices_.push_back(e);
  while (channelNum > rl_inv_r_matrices_[jGroupNum - 1].size()) rl_inv_r_matrices_[jGroupNum - 1].push_back(f);
  rl_inv_r_matrices_[jGroupNum - 1][channelNum - 1].push_back(matrixElement);
  assert(channelPrimeNum = rl_inv_r_matrices_[jGroupNum - 1][channelNum - 1].size());
}

/*!
 * This function adds an entire \f$ [1-RL]^{-1} \f$ matrix to a vector.
 */

void RMatrixFunc::AddRLInvMatrix(matrix_c matrix) {
  rl_inv_matrices_.push_back(matrix);
}
