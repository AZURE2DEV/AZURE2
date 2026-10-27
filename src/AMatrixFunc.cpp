#include "AMatrixFunc.h"
#include "Profile.h"
#include "AngCoeff.h"
#include "PolarizationFunc.h"
#include "CNuc.h"
#include "Config.h"
#include "EPoint.h"
#include "AZUREGrad.h"
#include "Decay.h"
#include "KGroup.h"
#include "KLGroup.h"
#include "MGroup.h"
#include "Interference.h"
#include "JGroup.h"
#include "ALevel.h"
#include "AChannel.h"
#include "PPair.h"
#include <assert.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif

/*!
 * The AMatrixFunc object is created with reference to a CNuc object.
 */

AMatrixFunc::AMatrixFunc(CNuc *compound, const Config &configure) :
  configure_(&configure),
  compound_(compound),
  cached_max_levels_(0),
  cached_max_channels_(0) {}

/*!
 * Returns the A-Matrix of a J-Group, materializing it from the stored
 * factorization on first use.  Returns an empty matrix if the J-Group has no
 * usable factorization.
 */

const matrix_c &AMatrixFunc::GetAMatrix(int jGroupNum) const {
  static const matrix_c emptyMatrix;
  if (jGroupNum < 1 || jGroupNum > (int)solvers_.size()) return emptyMatrix;
  if (route_[jGroupNum - 1] == kRouteChannel) {
    // A = V A_f V^T, one column per unit vector.  The cache is part of the
    // object's per-energy state, like the factorization's own inverse.
    AMatrixFunc *self = const_cast<AMatrixFunc *>(this);
    ChannelState &st = self->chan_[jGroupNum - 1];
    if (!st.inverseValid) {
      const int n = (*formal_[jGroupNum - 1]).n;
      st.inverse.assign(n, vector_c(n, complex(0.0, 0.0)));
      vector_c unit(n), col(n);
      for (int c = 0; c < n; c++) {
        std::fill(unit.begin(), unit.end(), complex(0.0, 0.0));
        unit[c] = 1.0;
        self->ApplyA(jGroupNum, unit.data(), col.data());
        for (int r = 0; r < n; r++) st.inverse[r][c] = col[r];
      }
      st.inverseValid = true;
    }
    return st.inverse;
  }
  if (!solver_valid_[jGroupNum - 1]) return emptyMatrix;
  return solvers_[jGroupNum - 1].Inverse();
}

/*!
 * Returns an A-Matrix element specified by positions in the JGroup and ALevel vectors.
 * The level numbers are the original ones; they are mapped onto the compacted
 * indices of the levels actually included in the R-Matrix.
 */

complex AMatrixFunc::GetAMatrixElement(int jGroupNum, int lambdaNum, int muNum) const {
  if (jGroupNum < 1 || jGroupNum > (int)solvers_.size()) return complex(0.0, 0.0);
  const bool channel = (route_[jGroupNum - 1] == kRouteChannel);
  if (!channel && !solver_valid_[jGroupNum - 1]) return complex(0.0, 0.0);
  const std::vector<int> &map = level_active_index_[jGroupNum - 1];
  if (lambdaNum < 1 || lambdaNum >= (int)map.size()) return complex(0.0, 0.0);
  if (muNum < 1 || muNum >= (int)map.size()) return complex(0.0, 0.0);
  const int row = map[lambdaNum];
  const int col = map[muNum];
  if (row == 0 || col == 0) return complex(0.0, 0.0);
  if (channel) return this->GetAMatrix(jGroupNum)[row - 1][col - 1];
  return solvers_[jGroupNum - 1].Inverse()[row - 1][col - 1];
}

/*!
 * Returns a pointer to an entire A-Matrix specified by a position in the JGroup vector.
 */

matrix_c *AMatrixFunc::GetJSpecAInvMatrix(int jGroupNum) {
  matrix_c *b = &a_inv_matrices_[jGroupNum - 1];
  return b;
}

/*!
 * Clears all matrices associated with the AMatrixFunc object.
 */

void AMatrixFunc::ClearMatrices() {
  tbar_decay_ = nullptr;
  tmatrix_.clear();
  ec_tmatrix_.clear();
  InvalidateTMatrixProducts();
  chan_cap_cache_.clear();

  const int numJGroups = compound()->NumJGroups();

  // Grow the per-JGroup storage rather than rebuilding it, so an instance reused
  // across energy points reallocates nothing after its first point.
  if ((int)a_inv_matrices_.size() != numJGroups) {
    a_inv_matrices_.resize(numJGroups);
    solvers_.resize(numJGroups);
    solver_valid_.resize(numJGroups);
    level_active_index_.resize(numJGroups);
    gamma_vectors_.resize(numJGroups);
    gamma_valid_.resize(numJGroups);
    channel_solves_.resize(numJGroups);
    channel_solve_valid_.resize(numJGroups);
    raw_solves_.resize(numJGroups);
    raw_solve_valid_.resize(numJGroups);
    route_.resize(numJGroups);
    formal_.resize(numJGroups);
    formal_jgroup_.assign(numJGroups, nullptr);
    chan_.resize(numJGroups);
    active_levels_.resize(numJGroups);
  }

  all_levels_active_ = true;
  for (int j = 1; j <= numJGroups; ++j) {
    solver_valid_[j - 1] = 0;
    route_[j - 1] = kRouteNone;
    if (!compound()->GetJGroup(j)->IsInRMatrix()) {
      // keep empty placeholders to preserve external indexing
      level_active_index_[j - 1].clear();
      a_inv_matrices_[j - 1].clear();
      gamma_valid_[j - 1].clear();
      channel_solve_valid_[j - 1].clear();
      raw_solve_valid_[j - 1].clear();
      continue;
    }

    JGroup *jg = compound()->GetJGroup(j);
    const int numLevels = jg->NumLevels();
    const int numChannels = jg->NumChannels();

    // Build original->active map (1-based for readability in existing code)
    std::vector<int> &map = level_active_index_[j - 1];
    std::vector<int> &levels = active_levels_[j - 1];
    map.assign(numLevels + 1, 0);
    levels.clear();
    int act = 0;
    for (int la = 1; la <= numLevels; ++la) {
      if (jg->GetLevel(la)->IsInRMatrix()) {
        map[la] = ++act;
        levels.push_back(la);
      } else
        all_levels_active_ = false;
    }
    // The level matrix itself is sized and zeroed only for the level route
    // (FillLevelMatrix).

    gamma_vectors_[j - 1].resize(numChannels);
    gamma_valid_[j - 1].assign(numChannels, 0);
    channel_solves_[j - 1].resize(numChannels);
    channel_solve_valid_[j - 1].assign(numChannels, 0);
    raw_solves_[j - 1].resize(numChannels);
    raw_solve_valid_[j - 1].assign(numChannels, 0);
  }
}


/*!
 * This function creates the inverted A-Matrix from the parameters in the CNuc object.
 */

void AMatrixFunc::FillMatrices(EPoint *point) {
  AZUREProfile::Timer profileTimer(AZUREProfile::kFill);
  const double inEnergy = this->InEnergy(point);
  // Groups of another compound nucleus are not part of this point's matrices
  // (EPoint::CalcEDependentValues left their rows empty).
  activeNucleus_ = compound()->GetPair(compound()->GetPairNumFromKey(point->GetEntranceKey()))->GetNucleus();
  const bool channelSpace = !level_only_ && ChannelSpace::Enabled();
  for (int j = 1; j <= compound()->NumJGroups(); j++) {
    JGroup *jGroup = compound()->GetJGroup(j);
    if (!jGroup->IsInRMatrix() || jGroup->GetNucleus() != activeNucleus_) continue;
    route_[j - 1] = kRouteLevel;
    if (channelSpace) {
      if (!formal_[j - 1] || formal_jgroup_[j - 1] != jGroup) {
        formal_[j - 1] = ChannelSpace::Cached(jGroup);
        formal_jgroup_[j - 1] = jGroup;
      }
      ChannelSpace::FormalLevels &f = (*formal_[j - 1]);
      ChannelSpace::Update(jGroup, active_levels_[j - 1], configure(), f);
      if (f.usable) {
        if (this->SolveChannelSpace(j, point, inEnergy)) {
          route_[j - 1] = kRouteChannel;
          AZUREProfile::Count(AZUREProfile::kChannelRoute);
          static const bool verify = std::getenv("AZURE_CHANNEL_VERIFY") != nullptr;
          if (verify) this->VerifyChannelSpace(j, point, inEnergy);
          continue;
        }
      } else {
        AZUREProfile::Count(f.status == 1 ? AZUREProfile::kLevelRouteSize : AZUREProfile::kLevelRouteSetup);
      }
    } else {
      AZUREProfile::Count(AZUREProfile::kLevelRouteForced);
    }
    this->FillLevelMatrix(j, point, inEnergy);
  }
}

/*!
 * The level matrix of J-group j at one energy (the level route).
 */

void AMatrixFunc::FillLevelMatrix(int j, EPoint *point, double inEnergy) {
  const bool park = (configure().paramMask & Config::USE_PARK_FORMALISM);
  {
    const int act = (int)active_levels_[j - 1].size();
    matrix_c &M = a_inv_matrices_[j - 1];
    M.resize(act);
    for (int r = 0; r < act; ++r) M[r].assign(act, complex(0.0, 0.0));
  }
  // Cache JGroup pointer to avoid repeated calls
  JGroup *jGroup = compound()->GetJGroup(j);
  int numLevels = jGroup->NumLevels();
  int numChannels = jGroup->NumChannels();

  // Ensure pre-allocated buffers are large enough - resize only when needed
  if (numLevels > cached_max_levels_ || numChannels > cached_max_channels_) {
    cached_max_levels_ = std::max(numLevels, cached_max_levels_);
    cached_max_channels_ = std::max(numChannels, cached_max_channels_);

    levelGammas_.resize(cached_max_levels_ + 1);
    levelEnergies_.resize(cached_max_levels_ + 1);
    shiftFunctions_.resize(cached_max_levels_ + 1);

    for (int la = 0; la <= cached_max_levels_; la++) {
      levelGammas_[la].resize(cached_max_channels_ + 1);
      shiftFunctions_[la].resize(cached_max_channels_ + 1);
    }
  }

  // Pre-cache gamma values for all levels and channels using pre-allocated buffers
  for (int la = 1; la <= numLevels; la++) {
    if (jGroup->GetLevel(la)->IsInRMatrix()) {
      ALevel *level = jGroup->GetLevel(la);
      levelEnergies_[la] = level->GetFitE();

      for (int ch = 1; ch <= numChannels; ch++) {
        levelGammas_[la][ch] = level->GetFitGamma(ch);
        shiftFunctions_[la][ch] = level->GetShiftFunction(ch);
      }
    }
  }

  for (int la = 1; la <= numLevels; la++) {
    if (jGroup->GetLevel(la)->IsInRMatrix()) {
      ALevel *level = jGroup->GetLevel(la);
      for (int lap = 1; lap <= numLevels; lap++) {
        if (jGroup->GetLevel(lap)->IsInRMatrix()) {
          ALevel *levelp = jGroup->GetLevel(lap);
          complex sum(0.0, 0.0);
          double parkOverlap = 0.0;
          for (int ch = 1; ch <= numChannels; ch++) {
            double gammaCh = levelGammas_[la][ch];
            double gammaChp = levelGammas_[lap][ch];

            // Early termination for effectively zero gamma values
            if (fabs(gammaCh) < 1.0e-12 || fabs(gammaChp) < 1.0e-12) continue;

            complex loElement = point->GetLoElement(j, ch);
            sum += gammaCh * gammaChp * loElement;

            // Cache channel pointer and radiation type to avoid repeated calls
            AChannel *channel = jGroup->GetChannel(ch);
            char radType = channel->GetRadType();

            if ((radType == 'M' || radType == 'E') &&
                la == lap &&
                (configure().paramMask & Config::USE_RMC_FORMALISM))
              sum += complex(0.0, 1.0) * gammaCh * gammaChp;
            if (park && radType == 'P') {
              // Park (2021) Eq. (22) with B_{lambda c} = S_c(E_lambda), Eq. (26):
              // the kernel is S_c(E) - S_c(E_lambda) + i P_c, where loElement
              // is S_c(E) - B_c + i P_c.  Eq. (9) gives the off-diagonal overlap.
              sum += gammaCh * gammaChp * (channel->GetBoundaryCondition() - shiftFunctions_[la][ch]);
              if (la != lap)
                parkOverlap -= gammaCh * gammaChp * (shiftFunctions_[la][ch] - shiftFunctions_[lap][ch]) /
                    (levelEnergies_[la] - levelEnergies_[lap]);
            } else if ((configure().paramMask & Config::USE_BRUNE_FORMALISM) && radType == 'P') {
              sum += gammaCh * gammaChp * channel->GetBoundaryCondition();
              if (la == lap)
                sum -= gammaCh * gammaChp * shiftFunctions_[la][ch];
              else
                sum -= gammaCh * gammaChp *
                    (shiftFunctions_[la][ch] * (inEnergy - levelEnergies_[lap]) - shiftFunctions_[lap][ch] * (inEnergy - levelEnergies_[la])) /
                    (levelEnergies_[la] - levelEnergies_[lap]);
            }
          }
          if (park) {
            // Park (2021) Eqs. (21)-(22): A^-1 = (E_lambda - E) J - gamma (S - B_lambda + iP) gamma^T,
            // with the diagonal overlap J_{lambda lambda} of Eq. (28).
            if (la == lap) parkOverlap = level->GetParkNorm();
            this->AddAInvMatrixElement(j, la, lap, (levelEnergies_[la] - inEnergy) * parkOverlap - sum);
          } else if (la == lap) {
            double resenergy = levelEnergies_[la];
            this->AddAInvMatrixElement(j, la, lap, resenergy - inEnergy - sum);
          } else
            this->AddAInvMatrixElement(j, la, lap, -sum);
        }
      }
    }
  }
}

/*!
 * This function factors the inverse A-Matrix of each J-Group.  The A-Matrix
 * itself is not formed: CalculateTMatrix needs only bilinear forms of it, which
 * the factorization supplies with one triangular solve per channel, and the
 * callers that do need individual elements get them from GetAMatrixElement,
 * which materializes the inverse on demand.
 */

void AMatrixFunc::InvertMatrices() {
  AZUREProfile::Timer profileTimer(AZUREProfile::kFactor);
  for (int j = 1; j <= compound()->NumJGroups(); j++) {
    if (compound()->GetJGroup(j)->IsInRMatrix() && compound()->GetJGroup(j)->GetNucleus() == activeNucleus_) {
      // Anything derived from the previous solve is now stale.
      gamma_valid_[j - 1].assign(gamma_valid_[j - 1].size(), 0);
      channel_solve_valid_[j - 1].assign(channel_solve_valid_[j - 1].size(), 0);
      raw_solve_valid_[j - 1].assign(raw_solve_valid_[j - 1].size(), 0);
      if (route_[j - 1] != kRouteLevel) continue;  // solved in channel space by FillMatrices
      matrix_c *theAInvMatrix = this->GetJSpecAInvMatrix(j);
      // Add validation to catch corrupted matrices before processing
      if (theAInvMatrix->empty()) {
        continue;  // Skip empty matrices
      }
      solver_valid_[j - 1] = solvers_[j - 1].Decompose(*theAInvMatrix) ? 1 : 0;
    }
  }
  chan_cap_cache_.clear();
}

/*!
 * Returns the reduced widths of one channel over the active levels of a J-Group.
 * Widths below 1e-12 are zeroed so that the bilinear form reproduces, term for
 * term, the level pairs the explicit double loop used to skip.
 */

const vector_r &AMatrixFunc::GetGammaVector(int jGroupNum, int chNum) {
  vector_r &gammas = gamma_vectors_[jGroupNum - 1][chNum - 1];
  if (!gamma_valid_[jGroupNum - 1][chNum - 1]) {
    JGroup *jGroup = compound()->GetJGroup(jGroupNum);
    const std::vector<int> &map = level_active_index_[jGroupNum - 1];
    gammas.assign(active_levels_[jGroupNum - 1].size(), 0.0);
    for (int la = 1; la < (int)map.size(); la++) {
      const int row = map[la];
      if (row == 0) continue;
      double gamma = jGroup->GetLevel(la)->GetFitGamma(chNum);
      if (fabs(gamma) < 1.0e-12) continue;
      gammas[row - 1] = gamma;
    }
    gamma_valid_[jGroupNum - 1][chNum - 1] = 1;
  }
  return gammas;
}

/*!
 * Returns A times the reduced-width vector of one channel, cached for the
 * duration of the energy point.  Every pathway sharing an exit channel reuses
 * the same solve, which is what replaces the old per-pathway double loop.
 */

const vector_c &AMatrixFunc::GetChannelSolve(int jGroupNum, int chNum) {
  vector_c &solve = channel_solves_[jGroupNum - 1][chNum - 1];
  if (!channel_solve_valid_[jGroupNum - 1][chNum - 1]) {
    AZUREProfile::Timer profileTimer(AZUREProfile::kChannelSolve);
    const vector_r &gammas = this->GetGammaVector(jGroupNum, chNum);
    const int n = (int)gammas.size();
    solve.resize(n);
    for (int i = 0; i < n; i++) solve[i] = gammas[i];
    if (n > 0) solvers_[jGroupNum - 1].Solve(&solve[0]);
    channel_solve_valid_[jGroupNum - 1][chNum - 1] = 1;
  }
  return solve;
}

/*!
 * The level matrix is complex symmetric, so A is too and the bilinear form may
 * contract either index with the cached solve.
 */

complex AMatrixFunc::GetUBilinear(int jGroupNum, int chNum, int chpNum) {
  if (jGroupNum < 1 || jGroupNum > (int)solvers_.size()) return complex(0.0, 0.0);
  if (route_[jGroupNum - 1] == kRouteChannel) {
    // F = G^T A G is symmetric: one column per exit channel serves every entrance.
    const ChannelSpace::FormalLevels &f = (*formal_[jGroupNum - 1]);
    const int a = f.pos[chNum - 1], b = f.pos[chpNum - 1];
    if (a < 0 || b < 0) return complex(0.0, 0.0);
    return this->FColumn(jGroupNum, b)[a];
  }
  if (!solver_valid_[jGroupNum - 1]) return complex(0.0, 0.0);
  // Fill the solve cache first: it populates the gamma cache for chpNum, so
  // taking the chNum reference afterwards keeps it clear that nothing can move
  // out from under it.
  const vector_c &solve = this->GetChannelSolve(jGroupNum, chpNum);
  const vector_r &gammas = this->GetGammaVector(jGroupNum, chNum);
  complex sum(0.0, 0.0);
  for (size_t i = 0; i < solve.size(); i++) sum += gammas[i] * solve[i];
  return sum;
}

/*!
 * Channel capture excludes levels with a negligible width in the internal
 * channel from both indices, which is a mask on the level set rather than the
 * per-channel cutoff, so this path builds its vectors fresh.  Results are
 * memoized because the same pathway recurs across k-groups.
 */

complex AMatrixFunc::GetChannelCaptureBilinear(int jGroupNum, int chNum, int chpNum,
                                               int maskChannel) {
  if (jGroupNum < 1 || jGroupNum > (int)solvers_.size()) return complex(0.0, 0.0);
  for (size_t i = 0; i < chan_cap_cache_.size(); i++) {
    const ChanCapEntry &entry = chan_cap_cache_[i];
    if (entry.jGroupNum == jGroupNum && entry.chNum == chNum &&
        entry.chpNum == chpNum && entry.maskChannel == maskChannel) return entry.value;
  }
  complex value(0.0, 0.0);
  const bool channel = (route_[jGroupNum - 1] == kRouteChannel);
  if (channel || solver_valid_[jGroupNum - 1]) {
    JGroup *jGroup = compound()->GetJGroup(jGroupNum);
    const std::vector<int> &map = level_active_index_[jGroupNum - 1];
    const int n = (int)active_levels_[jGroupNum - 1].size();
    chan_cap_entrance_.assign(n, complex(0.0, 0.0));
    chan_cap_exit_.assign(n, complex(0.0, 0.0));
    for (int la = 1; la < (int)map.size(); la++) {
      const int row = map[la];
      if (row == 0) continue;
      ALevel *level = jGroup->GetLevel(la);
      if (maskChannel && fabs(level->GetFitGamma(maskChannel)) < 1.0e-8) continue;
      chan_cap_entrance_[row - 1] = level->GetFitGamma(chNum);
      chan_cap_exit_[row - 1] = level->GetFitGamma(chpNum);
    }
    value = channel ? this->ChannelBilinear(jGroupNum, chan_cap_entrance_.data(), chan_cap_exit_.data())
                    : solvers_[jGroupNum - 1].Bilinear(chan_cap_entrance_, chan_cap_exit_);
  }
  ChanCapEntry entry;
  entry.jGroupNum = jGroupNum;
  entry.chNum = chNum;
  entry.chpNum = chpNum;
  entry.maskChannel = maskChannel;
  entry.value = value;
  chan_cap_cache_.push_back(entry);
  return value;
}

/*!
 * This function calculates the T-Matrix for each reaction pathway based on the A-Matrix.
 */

void AMatrixFunc::CalculateTMatrix(EPoint *point) {
  AZUREProfile::Timer profileTimer(AZUREProfile::kTMatrix);
  tbar_decay_ = nullptr;
  // Cache frequently accessed values to avoid repeated function calls
  int entranceKey = point->GetEntranceKey();
  int exitKey = point->GetExitKey();
  int aa = compound()->GetPairNumFromKey(entranceKey);
  int exitPairNum = compound()->GetPairNumFromKey(exitKey);

  int irEnd;
  int irStart;
  bool isRMC = false;
  if ((configure().paramMask & Config::USE_RMC_FORMALISM) &&
      compound()->GetPair(exitPairNum)->GetPType() == 10) {
    irStart = 1;
    irEnd = compound()->GetPair(aa)->NumDecays();
    isRMC = true;
  } else {
    irStart = 0;
    int numDecays = compound()->GetPair(aa)->NumDecays();
    while (irStart < numDecays) {
      irStart++;
      if (compound()->GetPair(aa)->GetDecay(irStart)->GetPairNum() == exitPairNum) break;
    }
    irEnd = irStart;
  }
  for (int ir = irStart; ir <= irEnd; ir++) {
    Decay *theDecay = compound()->GetPair(aa)->GetDecay(ir);
    for (int k = 1; k <= theDecay->NumKGroups(); k++) {
      for (int m = 1; m <= theDecay->GetKGroup(k)->NumMGroups(); m++) {
        MGroup *theMGroup = theDecay->GetKGroup(k)->GetMGroup(m);
        // Cache frequently accessed values to avoid repeated function calls
        int jNum = theMGroup->GetJNum();
        int chNum = theMGroup->GetChNum();
        int chpNum = theMGroup->GetChpNum();

        JGroup *theJGroup = compound()->GetJGroup(jNum);
        AChannel *entranceChannel = theJGroup->GetChannel(chNum);
        AChannel *exitChannel = theJGroup->GetChannel(chpNum);

        // Cache phase calculations
        complex coulombPhaseEn = point->GetExpCoulombPhase(jNum, chNum);
        complex hardSpherePhaseEn = point->GetExpHardSpherePhase(jNum, chNum);
        complex coulombPhaseEx = point->GetExpCoulombPhase(jNum, chpNum);
        complex hardSpherePhaseEx = point->GetExpHardSpherePhase(jNum, chpNum);
        complex sqrtPenEn = point->GetSqrtPenetrability(jNum, chNum);
        complex sqrtPenEx = point->GetSqrtPenetrability(jNum, chpNum);

        complex uphase = coulombPhaseEn * hardSpherePhaseEn * coulombPhaseEx * hardSpherePhaseEx;
        complex umatrix = 2.0 * complex(0.0, 1.0) * sqrtPenEn * sqrtPenEx *
            this->GetUBilinear(jNum, chNum, chpNum);

        complex tphase = coulombPhaseEn * coulombPhaseEn;
        complex tmatrix;
        if (isRMC)
          this->AddTMatrixElement(k, m, complex(0.0, -1.0) * umatrix, ir);
        else {
          if (chNum == chpNum) {
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
          int internalChannel = theECMGroup->GetIntChannelNum();
          MGroup *chanMGroup = compound()->GetPair(aa)->GetDecay(theECMGroup->GetChanCapDecay())->GetKGroup(theECMGroup->GetChanCapKGroup())->GetMGroup(theECMGroup->GetChanCapMGroup());
          int maskChannel = (configure().paramMask & Config::IGNORE_ZERO_WIDTHS) ? internalChannel : 0;
          complex umatrix = 2.0 * complex(0.0, 1.0) *
              point->GetSqrtPenetrability(chanMGroup->GetJNum(), chanMGroup->GetChNum()) *
              this->GetChannelCaptureBilinear(chanMGroup->GetJNum(), chanMGroup->GetChNum(),
                                              chanMGroup->GetChpNum(), maskChannel);
          tmatrix = tmatrix * umatrix;
        }
        this->AddECTMatrixElement(k, m, tmatrix);
      }
    }
  }
}

/*!
 * This function adds an inverse A-Matrix element specified by positions in the JGroup and ALevel vectors.
 */

void AMatrixFunc::AddAInvMatrixElement(int jGroupNum, int lambdaNum, int muNum, complex aMatrixElement) {
  // Basic bounds
  if (jGroupNum < 1 || jGroupNum > (int)a_inv_matrices_.size()) return;

  // If this J-group wasn’t included, skip
  if (level_active_index_.empty() || level_active_index_[jGroupNum - 1].empty()) return;

  // Translate original level indices to compact active indices
  const auto &map = level_active_index_[jGroupNum - 1];
  if (lambdaNum < 1 || lambdaNum >= (int)map.size()) return;
  if (muNum < 1 || muNum >= (int)map.size()) return;

  const int row = map[lambdaNum];
  const int col = map[muNum];
  if (row == 0 || col == 0) return;  // either level is inactive; logic says: ignore

  // Validate value
  if (!std::isfinite(aMatrixElement.real()) || !std::isfinite(aMatrixElement.imag())) return;

  // Direct assignment into pre-sized dense matrix
  a_inv_matrices_[jGroupNum - 1][row - 1][col - 1] = aMatrixElement;
}

/*!
 * Reverse-mode adjoint of one energy point's cross section (angle-integrated or
 * differential), accumulating dlnL/dE and dlnL/dgamma into `accum`.
 *
 * Adjoint convention (consistent throughout, see PLAN.md "Gotchas"):
 *   For a complex intermediate z and the real scalar lnL we store the cotangent
 *     zbar := dlnL/dRe(z) + i dlnL/dIm(z).
 *   - A linear map z = c * a (c complex constant) back-propagates as
 *         abar += conj(c) * zbar.
 *   - The contribution of z to a *real* parameter p is
 *         dlnL/dp += Re( conj(zbar) * dz/dp ).
 *   - The holomorphic inverse A = M^{-1} back-propagates as
 *         Mbar = -A^H * Abar * A^H.
 *
 * Forward chain reproduced here (per point, standard A-matrix):
 *   T(k,m) = [tphase] - uphase ( [1] + umatrix )         ([..] only if chNum==chpNum)
 *   umatrix = sum_{l,l'} 2i sqrtP_en sqrtP_ex gamma_en,l gamma_ex,l' A_j[l][l']
 *   A_j = M_j^{-1},   M_j[r][c] = (E_r - inE) delta_rc - sum_ch gamma_r,ch gamma_c,ch L_ch
 * with the cross-section -> T step depending on whether the point is
 * angle-integrated (quadratic |T_temp|^2 sum) or differential (Legendre
 * interference sum + elastic Coulomb interference).  Returns false (touching
 * nothing) for configurations outside the supported standard path.
 */
bool AMatrixFunc::PointAdjoint(EPoint *point, double fitBar, GradAccum &accum,
                               const vector_matrix_r *shiftDeriv, int xsComponent,
                               bool ayCrossSection) {
  TMatrixCotangent bar;
  if (!this->PointTBar(point, fitBar, bar, xsComponent, ayCrossSection)) return false;
  this->TBarAdjoint(point, bar, accum, shiftDeriv);
  return true;
}

/*!
 * Cross section -> T: the T-matrix cotangents of one point's observable
 * (see PointAdjoint).  Reads the T-matrix held by this object and the
 * angle-dependent quantities of `point`.
 */
bool AMatrixFunc::PointTBar(EPoint *point, double fitBar, TMatrixCotangent &bar, int xsComponent,
                            bool ayCrossSection) {
  AZUREProfile::Timer profileTimer(AZUREProfile::kGradTBar);
  const complex I(0.0, 1.0);

  // --- Common gating. ---
  if (configure().paramMask & Config::USE_RMC_FORMALISM) return false;
  if (point->IsAngularDist() && !point->IsCoefficient()) return false;  // angular-dist test points
  // xsComponent < 0 asks for A_0 = sigma_int / 4 pi of an angular-distribution point: the angle-integrated
  // cotangent scaled by 1 / 4 pi (GradTargetEffectAdjoint, normalised coefficients under a target effect).
  const bool a0Mode = (xsComponent < 0);
  if (a0Mode) xsComponent = 0;
  // E1/E2 component selection only applies to the angle-integrated capture XS.
  if (xsComponent != 0 && (point->IsDifferential() || point->IsPhase())) return false;

  int aa = compound()->GetPairNumFromKey(point->GetEntranceKey());
  int exitPairNum = compound()->GetPairNumFromKey(point->GetExitKey());
  // Beta-delayed particle emission is supported below; its observable is a plain
  // incoherent sum over pathways, so only the T cotangents differ from the
  // cross-section case. E1/E2 component selection has no meaning for it.
  const bool isBetaDelayed = (compound()->GetPair(aa)->GetPType() == 20);
  if (isBetaDelayed && xsComponent != 0) return false;

  // Active-level compaction unsupported here: require every level of every
  // in-R-matrix JGroup to itself be in the R matrix (so original idx == active;
  // ClearMatrices records it).
  if (!all_levels_active_) return false;

  // Locate the decay matching the exit pair (mirror CalculateCrossSection).
  int ir = 0;
  while (ir < compound()->GetPair(aa)->NumDecays()) {
    ir++;
    if (compound()->GetPair(aa)->GetDecay(ir)->GetPairNum() == exitPairNum) break;
  }
  Decay *theDecay = compound()->GetPair(aa)->GetDecay(ir);

  const double geom = point->GetGeometricalFactor();
  const double i1i2 = compound()->GetPair(aa)->GetI1I2Factor();
  const bool isDiff = point->IsDifferential();
  const bool isPhase = point->IsPhase();
  if (isPhase && aa != exitPairNum) return false;  // phase output is for elastic only
  // The cross section carried by an analyzing-power point is the differential
  // one at the point's angle (GenMatrixFunc computes it before A_y and a
  // target-effect sub-point keeps it); anything else is not supported.
  if (ayCrossSection && (!point->IsAnalyzingPower() || !isDiff || isPhase || xsComponent != 0)) return false;

  // Electron screening (GenMatrixFunc::CalculateCrossSection): a screened observable is f times the bare
  // one, so the bare one's cotangent is f fitBar, and dY/dU_e = k Y per MeV, 1e-6 k Y per eV, with Y the
  // screened value this point holds from its forward calculation.  A phase shift and an A_y value are
  // not screened (an A_y sub-point's cross section is: ayCrossSection).
  bar.screen.clear();
  if (point->GetScreeningParam() > 0 && point->GetScreeningK() != 0.0 && !isPhase && !point->IsCoefficient() &&
      !(point->IsAnalyzingPower() && !ayCrossSection)) {
    double y = point->GetFitCrossSection();
    if (xsComponent == 1) y = point->GetFitE1CrossSection();
    else if (xsComponent == 2) y = point->GetFitE2CrossSection();
    bar.screen.push_back({point->GetScreeningParam(), fitBar * 1.0e-6 * point->GetScreeningK() * y});
    fitBar *= point->ScreeningFactor();
  }

  // === Cross-section -> T:  build the T-matrix cotangents. ===
  // tBar  : cotangents of the internal T(k,m).
  // ecBar : cotangents of the external-capture T(k,ecm).
  int nK = theDecay->NumKGroups();
  std::vector<vector_c> &tBar = bar.t, &ecBar = bar.ec;
  tBar.assign(nK, vector_c());
  ecBar.assign(nK, vector_c());
  for (int k = 1; k <= nK; k++) {
    tBar[k - 1].assign(theDecay->GetKGroup(k)->NumMGroups(), complex(0.0, 0.0));
    ecBar[k - 1].assign(theDecay->GetKGroup(k)->NumECMGroups(), complex(0.0, 0.0));
  }

  // Keep only the contributions feeding the requested cross-section component
  // (0 = full, 1 = E1, 2 = E2); E1/E2 select electric channels of that L.
  auto includeInternal = [&](MGroup *mg) -> bool {
    if (xsComponent == 0) return true;
    AChannel *ex = compound()->GetJGroup(mg->GetJNum())->GetChannel(mg->GetChpNum());
    return ex->GetRadType() == 'E' && ex->GetL() == xsComponent;
  };
  auto includeEC = [&](ECMGroup *ecg) -> bool {
    if (xsComponent == 0) return true;
    return ecg->GetRadType() == 'E' && ecg->GetMult() == xsComponent;
  };

  // The angle-integrated cotangent as a function of its own cotangent `bar`: the chain below calls it with
  // fitBar; a normalised Legendre coefficient a_k = N_k / sigma_int calls it with -fitBar a_k / sigma_int.
  auto integratedBar = [&](double bar) {
    // ---- Angle-integrated: model = (1/100) sum_k sum_temp w_temp |T_temp|^2,
    //      w_temp = geom (2J+1) I1I2,  T_temp = sum_m T(k,m) sharing (J,l,l'). ----
    for (int k = 1; k <= nK; k++) {
      // Mirror the angle-integrated UPOS over-counting guard.
      if (aa != exitPairNum && compound()->GetPair(exitPairNum)->GetPType() == 0) {
        if (theDecay->GetKGroup(k)->GetSp() != theDecay->GetKGroup(k)->GetSp2()) continue;
      }
      int numM = theDecay->GetKGroup(k)->NumMGroups();
      int numEC = theDecay->GetKGroup(k)->NumECMGroups();
      struct Group {
        double j;
        int l;
        int lp;
        complex tsum;
      };
      std::vector<Group> groups;
      std::vector<int> mGroupIdx(numM + 1, -1);
      std::vector<int> ecGroupIdx(numEC + 1, -1);
      auto findGroup = [&](double jv, int lv, int lpv) -> int {
        for (int g = 0; g < (int)groups.size(); g++)
          if (groups[g].j == jv && groups[g].l == lv && groups[g].lp == lpv) return g;
        groups.push_back(Group{jv, lv, lpv, complex(0.0, 0.0)});
        return (int)groups.size() - 1;
      };
      // Internal pathways T(k,m), grouped by (J,l,l').
      for (int m = 1; m <= numM; m++) {
        MGroup *mg = theDecay->GetKGroup(k)->GetMGroup(m);
        if (!includeInternal(mg)) continue;
        JGroup *jg = compound()->GetJGroup(mg->GetJNum());
        int gi = findGroup(jg->GetJ(),
                           jg->GetChannel(mg->GetChNum())->GetL(),
                           jg->GetChannel(mg->GetChpNum())->GetL());
        groups[gi].tsum += this->GetTMatrixElement(k, m);
        mGroupIdx[m] = gi;
      }
      // External-capture pathways T(k,ecm), grouped by (J, L, mult) and merged
      // with the matching internal group (forward: same temp T-matrix list).
      for (int m = 1; m <= numEC; m++) {
        ECMGroup *ecg = theDecay->GetKGroup(k)->GetECMGroup(m);
        if (!includeEC(ecg)) continue;
        int gi = findGroup(ecg->GetJ(), ecg->GetL(), ecg->GetMult());
        groups[gi].tsum += this->GetECTMatrixElement(k, m);
        ecGroupIdx[m] = gi;
      }
      std::vector<complex> tBarGroup(groups.size());
      for (int g = 0; g < (int)groups.size(); g++) {
        double w = geom * (2.0 * groups[g].j + 1.0) * i1i2;
        tBarGroup[g] = bar * (2.0 / 100.0) * w * groups[g].tsum;
      }
      // accumulate (not assign): a normalised Legendre coefficient adds this to its order-k cotangents
      for (int m = 1; m <= numM; m++)
        if (mGroupIdx[m] >= 0) tBar[k - 1][m - 1] += tBarGroup[mGroupIdx[m]];
      for (int m = 1; m <= numEC; m++)
        if (ecGroupIdx[m] >= 0) ecBar[k - 1][m - 1] += tBarGroup[ecGroupIdx[m]];
    }
  };

  if (point->IsTotal()) {
    // ---- Total cross section by the optical theorem (GenMatrixFunc::CalculateBareCrossSection):
    //        model = (geom i1i2 / 100) sum (2J+1) 2 Re T(k,m) over the diagonal pathways,
    //      so d(model)/dRe T = (geom i1i2 / 100) 2 (2J+1), d(model)/dIm T = 0, and the cotangent (dL/dRe + i dL/dIm)
    //      is real. ----
    if (xsComponent != 0 || isDiff) return false;
    for (int k = 1; k <= nK; k++) {
      if (theDecay->GetKGroup(k)->IsCrossSpin()) continue;
      for (int m = 1; m <= theDecay->GetKGroup(k)->NumMGroups(); m++) {
        MGroup *mg = theDecay->GetKGroup(k)->GetMGroup(m);
        if (mg->GetChNum() != mg->GetChpNum()) continue;
        const double jValue = compound()->GetJGroup(mg->GetJNum())->GetJ();
        tBar[k - 1][m - 1] += fitBar * geom * i1i2 / 100. * 2. * (2. * jValue + 1.);
      }
    }
  } else if (isBetaDelayed) {
    // ---- Beta-delayed particle emission. GenMatrixFunc forms
    //        sum = sum_{k,m} 25 |T(k,m)|^2 ,   model = Re(sum)/100
    //      with no geometrical or spin factor and no external-capture term --
    //      an incoherent sum over pathways, with no interference between them.
    //      Hence d(model)/dT*(k,m) = (25/100) T(k,m), and the factor of two is
    //      the cotangent convention used throughout this function. ----
    for (int k = 1; k <= nK; k++)
      for (int m = 1; m <= theDecay->GetKGroup(k)->NumMGroups(); m++)
        tBar[k - 1][m - 1] += fitBar * (2.0 * 25.0 / 100.0) * this->GetTMatrixElement(k, m);
  } else if (point->IsAnalyzingPower() && !ayCrossSection) {
    // ---- Vector analyzing power. A_y is built from the channel-spin amplitude
    //      matrix, which is linear in T(k,m) with the coefficients AddPathway
    //      applies, so the whole derivative is a contraction against those same
    //      coefficients. The Coulomb amplitude carries no parameter dependence
    //      and drops out. ----
    if (point->NumSubPoints() > 0) return false;  // see the note below
    if (compound()->GetPair(aa)->GetPType() != 0) return false;

    if (compound()->GetPair(exitPairNum)->GetPType() == 10) {
      // ---- Capture analyzing power (Seyler and Weller, Phys. Rev. C 20
      //      (1979) 453). The observable is a ratio of two forms that are each
      //      bilinear in the T-matrix elements,
      //        A_y = N/D,  N = sum bk Re(i T1 T2*) P_k^1,
      //                    D = sum ak Re(T1 T2*)   P_k,
      //      so with d Re(T1 T2*)/dT1* = T2/2 and d Re(i T1 T2*)/dT1* = -i T2/2
      //      (and their partners with 1 <-> 2, i -> -i),
      //        dA_y/dT1* = [ -i bk P_k^1 T2 - A_y ak P_k T2 ] / (2 D),
      //        dA_y/dT2* = [ +i bk P_k^1 T1 - A_y ak P_k T1 ] / (2 D),
      //      and the cotangent convention of this function supplies the 2. ----
      const int maxL = point->GetMaxLOrder();
      compound()->CalcCaptureAnalyzingPower(aa, ir, maxL);
      const int nTerms = theDecay->NumCaptureAyTerms();
      if (nTerms == 0) return false;

      const double xCos = cos(point->GetCMAngle() * pi / 180.0);
      std::vector<double> assocP(maxL + 1, 0.0);
      for (int k = 1; k <= maxL; k++) assocP[k] = AngCoeff::LegendreP1(k, xCos);

      auto element = [&](int kg, int p, bool isEC) -> complex {
        return isEC ? this->GetECTMatrixElement(kg, p)
                    : this->GetTMatrixElement(kg, p);
      };
      double num = 0.0, den = 0.0;
      for (int i = 1; i <= nTerms; i++) {
        const CaptureAyTerm *t = theDecay->GetCaptureAyTerm(i);
        if (t->kOrder > maxL) continue;
        complex prod = element(t->kGroup1, t->path1, t->isEC1) * conj(element(t->kGroup2, t->path2, t->isEC2));
        if (t->ak != 0.0) den += t->ak * real(prod) * point->GetLegendreP(t->kOrder);
        if (t->bk != 0.0) num += t->bk * (-imag(prod)) * assocP[t->kOrder];
      }
      if (fabs(den) < 1.e-300) return false;
      const double ay = num / den;

      const complex I2(0.0, 1.0);
      for (int i = 1; i <= nTerms; i++) {
        const CaptureAyTerm *t = theDecay->GetCaptureAyTerm(i);
        if (t->kOrder > maxL) continue;
        const complex T1 = element(t->kGroup1, t->path1, t->isEC1);
        const complex T2 = element(t->kGroup2, t->path2, t->isEC2);
        const double pk = point->GetLegendreP(t->kOrder);
        const double pk1 = assocP[t->kOrder];
        const complex d1 = (-I2 * t->bk * pk1 * T2 - ay * t->ak * pk * T2) / den;
        const complex d2 = (I2 * t->bk * pk1 * T1 - ay * t->ak * pk * T1) / den;
        if (t->isEC1)
          ecBar[t->kGroup1 - 1][t->path1 - 1] += fitBar * d1;
        else
          tBar[t->kGroup1 - 1][t->path1 - 1] += fitBar * d1;
        if (t->isEC2)
          ecBar[t->kGroup2 - 1][t->path2 - 1] += fitBar * d2;
        else
          tBar[t->kGroup2 - 1][t->path2 - 1] += fitBar * d2;
      }
    } else if (compound()->GetPair(exitPairNum)->GetPType() != 0) {
      return false;
    } else {
      Polarization::AmplitudeMatrix M(compound(), point, aa, exitPairNum);
      for (int k = 1; k <= nK; k++) {
        if (theDecay->GetKGroup(k)->IsCrossSpin()) continue;  // as the forward sum
        for (int m = 1; m <= theDecay->GetKGroup(k)->NumMGroups(); m++) {
          MGroup *mg = theDecay->GetKGroup(k)->GetMGroup(m);
          M.AddPathway(mg->GetJNum(), mg->GetChNum(), mg->GetChpNum(),
                       this->GetTMatrixElement(k, m));
        }
      }
      if (aa == exitPairNum) M.AddCoulomb();
      if (M.size() == 0) return false;

      const std::vector<complex> bar = M.AnalyzingPowerBar();
      for (int k = 1; k <= nK; k++) {
        if (theDecay->GetKGroup(k)->IsCrossSpin()) continue;
        for (int m = 1; m <= theDecay->GetKGroup(k)->NumMGroups(); m++) {
          MGroup *mg = theDecay->GetKGroup(k)->GetMGroup(m);
          if (!includeInternal(mg)) continue;
          tBar[k - 1][m - 1] += fitBar * M.PathwayAdjoint(mg->GetJNum(), mg->GetChNum(), mg->GetChpNum(), bar);
        }
      }
    }
  } else if (isPhase) {
    // ---- Phase shift: model = (90/pi) * arg(U) [+ const for identical pairs],
    //      U = sum_{matching k,m} (expCP^2 - T(k,m)) / expCP^2, matched by
    //      (J, l, elastic channel) = (segment J, segment L). ----
    double segJ = point->GetJ();
    int segL = point->GetL();
    struct PM {
      int k;
      int m;
      complex invExpCP2;
    };
    std::vector<PM> matches;
    complex U(0.0, 0.0);
    for (int k = 1; k <= nK; k++) {
      for (int m = 1; m <= theDecay->GetKGroup(k)->NumMGroups(); m++) {
        MGroup *mg = theDecay->GetKGroup(k)->GetMGroup(m);
        JGroup *jg = compound()->GetJGroup(mg->GetJNum());
        AChannel *en = jg->GetChannel(mg->GetChNum());
        AChannel *ex = jg->GetChannel(mg->GetChpNum());
        if (jg->GetJ() == segJ && en->GetL() == segL && en == ex) {
          complex expCP = point->GetExpCoulombPhase(mg->GetJNum(), mg->GetChNum());
          complex expCP2 = expCP * expCP;
          U += (expCP2 - this->GetTMatrixElement(k, m)) / expCP2;
          matches.push_back(PM{k, m, 1.0 / expCP2});
        }
      }
    }
    double absU2 = std::norm(U);
    if (absU2 > 0.0) {
      // phase = (90/pi) arg(U) -> Ubar = fitBar (90/pi)(-Im U + i Re U)/|U|^2.
      complex Ubar = fitBar * (90.0 / pi) *
          complex(-std::imag(U), std::real(U)) / absU2;
      // U = sum [1 - T/expCP^2]  ->  dU/dT = -1/expCP^2.
      for (const PM &pm : matches)
        tBar[pm.k - 1][pm.m - 1] += conj(-pm.invExpCP2) * Ubar;
    }
  } else if (!isDiff && (!point->IsCoefficient() || a0Mode)) {
    integratedBar(a0Mode ? fitBar / (4.0 * pi) : fitBar);
  } else {
    // ---- Differential: model = (Re(CT)+Re(RT)+Re(IT))/100. ----
    //   RT = (geom I1I2 rtFactor/pi) Re( sum_{kL,inter} weight T1 conj(T2) P_L )
    //   IT = (geom itFactor/sqrt(pi)) Re( i sum_{k,m: en==ex} statSpin coulombAmp conj(T) P_l )  (elastic)
    //   CT = const (no T dependence).
    // Supports the standard case (entrance==exit pair, or gamma exit) with
    // external-capture interferences (RR/ER/RE/EE), and the inelastic-particle
    // UPOS / normal-angular-distribution case.
    PPair *exitPairDiff = compound()->GetPair(exitPairNum);
    double rtFactor = 1.0, itFactor = 1.0;
    if (aa == exitPairNum && compound()->GetPair(aa)->IsIdentical()) {
      rtFactor = 4.0;
      itFactor = 2.0;
    }
    double cR = geom * i1i2 * rtFactor / (pi * 100.0);
    // An angular-coefficient point (ESegment isDiff 9): the forward value is
    //   a_k = (4 geom i1i2 / (100 sigma_int)) sum_{L=k} z T1 T2*   (normalised)   or   A_k = a_k sigma_int / 4 pi,
    // i.e. the differential model's Legendre coefficient of order k, coeff_k = cR sum_{L=k} z T1 T2*, times
    // 4 pi / (rtFactor sigma_int) (normalised) or 1 / rtFactor (absolute).  So the differential machinery below
    // with a Kronecker delta in place of P_l(cos theta) gives the N_k part; the normalised one then adds the
    // sigma_int cotangent -fitBar a_k / sigma_int (integratedBar).  No Coulomb terms: elastic pairs are refused.
    const bool coef = point->IsCoefficient();
    const int kOrd = coef ? point->GetCoefficientOrder() : -1;
    const double sigmaInt = coef ? 4.0 * pi * point->GetAngDistA0() : 0.0;
    if (coef && (point->IsUPOS() || !(sigmaInt > 0.0))) return false;
    const double dfac = coef ? (point->GetCoefficientAbsolute() ? 1.0 / rtFactor : 4.0 * pi / (rtFactor * sigmaInt)) : 0.0;
    auto pl = [&](int l) -> double { return coef ? ((l == kOrd) ? dfac : 0.0) : point->GetLegendreP(l); };

    // Fetch a T value (internal 'R' or external-capture 'E') and route a
    // cotangent contribution to the right accumulator.
    auto Tval = [&](int K, int idx, char rad) -> complex {
      return (rad == 'E') ? this->GetECTMatrixElement(K, idx)
                          : this->GetTMatrixElement(K, idx);
    };
    auto Tadd = [&](int K, int idx, char rad, complex v) {
      if (rad == 'E')
        ecBar[K - 1][idx - 1] += v;
      else
        tBar[K - 1][idx - 1] += v;
    };
    // model += W Re(T1 conj T2)  ->  T1bar += fitBar W T2 ;  T2bar += fitBar W T1.
    auto addQuad = [&](int K, int m1, char r1, int m2, char r2, double W) {
      complex T1 = Tval(K, m1, r1), T2 = Tval(K, m2, r2);
      Tadd(K, m1, r1, fitBar * W * T2);
      Tadd(K, m2, r2, fitBar * W * T1);
    };

    const bool inelasticParticle = (aa != exitPairNum && exitPairDiff->GetPType() == 0);
    // Away from UPOS the sum is linear in the Legendre polynomials with
    // angle-independent coefficients: the cotangent of each pathway is a sum
    // over orders of products taken once per T-matrix (TBarProducts).
    if (!(inelasticParticle && point->IsUPOS()) && this->TBarProducts(theDecay, inelasticParticle)) {
      const complex scale = fitBar * cR;
      const std::vector<InterferencePathway> &paths = theDecay->InterferencePathways();
      const int nl = theDecay->NumInterferenceOrders();
      tbar_legendre_.assign(nl, 0.0);
      const int nP = std::min(nl, point->GetMaxLOrder() + 1);
      for (int l = 0; l < nP; l++) tbar_legendre_[l] = pl(l);
      for (size_t p = 0; p < paths.size(); p++) {
        const complex *q = &tbar_q_[p * nl];
        complex s(0.0, 0.0);
        for (int l = 0; l < nl; l++) s += q[l] * tbar_legendre_[l];
        if (s != complex(0.0, 0.0)) Tadd(paths[p].k, paths[p].m, paths[p].rad, scale * s);
      }
    } else
    for (int kL = 1; kL <= theDecay->NumKLGroups(); kL++) {
      int K = theDecay->GetKLGroup(kL)->GetK();
      int lOrder = theDecay->GetKLGroup(kL)->GetLOrder();
      double legP = pl(lOrder);
      KGroup *kg = theDecay->GetKGroup(K);
      for (int inter = 1; inter <= theDecay->GetKLGroup(kL)->NumInterferences(); inter++) {
        Interference *itf = theDecay->GetKLGroup(kL)->GetInterference(inter);
        const std::string &type = itf->GetInterferenceType();  // "RR","ER","RE","EE"
        if (type.size() != 2) return false;
        char rad1 = type[0], rad2 = type[1];  // 'R' internal, 'E' external
        int m1 = itf->GetM1(), m2 = itf->GetM2();

        if (inelasticParticle) {
          // Inelastic particle exit: pick the normal or UPOS sub-case.
          double sp1 = kg->GetSp(), sp2 = kg->GetSp2();
          int lp1 = compound()->GetJGroup(kg->GetMGroup(m1)->GetJNum())->GetChannel(kg->GetMGroup(m1)->GetChpNum())->GetL();
          int lp2 = compound()->GetJGroup(kg->GetMGroup(m2)->GetJNum())->GetChannel(kg->GetMGroup(m2)->GetChpNum())->GetL();
          if (sp1 == sp2 && !point->IsUPOS()) {
            addQuad(K, m1, rad1, m2, rad2, cR * itf->GetZ1Z2() * legP);
          }
          if (lp1 == lp2 && point->IsUPOS() && type == "RR") {
            // R_L coefficient (constant w.r.t. parameters), even lOrder only.
            double finalL = (double)point->GetSecondaryDecayL();
            double Ic = point->GetIc();
            double j2f = compound()->GetPair(exitPairNum)->GetJ(2);  // the exit pair, not the decay index (GenMatrixFunc)
            double delta = point->GetDelta();
            double R_L = 0.0;
            if ((int)lOrder % 2 == 0) {
              R_L = this->GetRk(j2f, finalL, finalL, Ic, lOrder);
              if (delta != 0.0) {
                double finalLp = finalL + 1.0;
                double R_LLp = this->GetRk(j2f, finalL, finalLp, Ic, lOrder);
                double R_LpLp = this->GetRk(j2f, finalLp, finalLp, Ic, lOrder);
                if (R_LpLp != 0.0)
                  R_L = (R_L + 2.0 * delta * R_LLp + delta * delta * R_LpLp) / (1.0 + delta * delta);
              }
            }
            double W = cR * itf->GetZ1Z2_UPOS() * std::pow(2.0 * lOrder + 1.0, 0.5) / 4.0 * R_L * legP;
            addQuad(K, m1, 'R', m2, 'R', W);
          }
        } else {
          // Standard case (entrance==exit pair, or gamma exit).
          addQuad(K, m1, rad1, m2, rad2, cR * itf->GetZ1Z2() * legP);
        }
      }
    }

    if (coef && !point->GetCoefficientAbsolute()) {
      // a_k = N_k / sigma_int: the sigma_int cotangent through the angle-integrated form
      integratedBar(-fitBar * point->GetFitCrossSection() / sigmaInt);
    }

    if (!coef && aa == exitPairNum && theDecay->TablesBuilt()) {  // elastic Coulomb-nuclear interference
      const complex cI = I * geom * itFactor / (std::sqrt(pi) * 100.0);
      for (const DiagonalPathway &d : theDecay->DiagonalPathways()) {
        const complex B = cI * d.statSpin * point->GetCoulombAmplitude(d.s) * point->GetLegendreP(d.l);
        tBar[d.k - 1][d.m - 1] += fitBar * B;
      }
    } else if (!coef && aa == exitPairNum) {  // elastic Coulomb-nuclear interference
      complex cI = I * geom * itFactor / (std::sqrt(pi) * 100.0);
      for (int k = 1; k <= nK; k++) {
        for (int m = 1; m <= theDecay->GetKGroup(k)->NumMGroups(); m++) {
          MGroup *mg = theDecay->GetKGroup(k)->GetMGroup(m);
          AChannel *en = compound()->GetJGroup(mg->GetJNum())->GetChannel(mg->GetChNum());
          AChannel *ex = compound()->GetJGroup(mg->GetJNum())->GetChannel(mg->GetChpNum());
          if (en == ex) {
            int l = en->GetL();
            // Per channel spin, as in the forward sum: only an identical pair
            // with spin has a Coulomb amplitude that depends on it.
            complex B = cI * mg->GetStatSpinFactor() * point->GetCoulombAmplitude(en->GetS()) * point->GetLegendreP(l);
            // model += Re(B conj T)  ->  Tbar += fitBar B.
            tBar[k - 1][m - 1] += fitBar * B;
          }
        }
      }
    }
  }

  return true;
}

/*!
 * T -> level energies and reduced widths (see PointAdjoint), at the energy of
 * `energyPoint`, whose T-matrix this object holds.
 */
void AMatrixFunc::TBarAdjoint(EPoint *energyPoint, const TMatrixCotangent &bar, GradAccum &accum,
                              const vector_matrix_r *shiftDeriv) {
  AZUREProfile::Timer profileTimer(AZUREProfile::kGradBackprop);
  const complex I(0.0, 1.0);
  for (const auto &s : bar.screen) accum.AddScreening(s.first, s.second);
  EPoint *point = energyPoint;
  // Brune formalism: the boundary/shift terms change the M-construction kernel
  // (handled in the structural gamma step via cached shift values) and add
  // explicit E_lambda dependence through S(E_lambda) (handled in the structural
  // energy step using shiftDeriv = dS/dE).
  // Park: the level matrix is (E_r - E) J[r][c] - sum_d g_r g_c (L_d + B_d - S_r,d),
  // with the overlap J of Park's Eqs. (9) and (28) depending on the widths and,
  // through dS/dE at E_lambda, on the level energies.  Its derivatives come from
  // the per-level values CNuc::CalcShiftFunctions stores (S, S', S'').
  const bool park = (configure().paramMask & Config::USE_PARK_FORMALISM);
  int aa = compound()->GetPairNumFromKey(point->GetEntranceKey());
  int exitPairNum = compound()->GetPairNumFromKey(point->GetExitKey());
  int ir = 0;
  while (ir < compound()->GetPair(aa)->NumDecays()) {
    ir++;
    if (compound()->GetPair(aa)->GetDecay(ir)->GetPairNum() == exitPairNum) break;
  }
  Decay *theDecay = compound()->GetPair(aa)->GetDecay(ir);
  bool hasEC = false;
  for (int k = 1; k <= theDecay->NumKGroups(); k++)
    if (theDecay->GetKGroup(k)->NumECMGroups() > 0) {
      hasEC = true;
      break;
    }
  const int nK = theDecay->NumKGroups();
  const std::vector<vector_c> &tBar = bar.t;
  const std::vector<vector_c> &ecBar = bar.ec;

  // Without external capture the level matrix enters only through bilinear
  // forms, whose adjoint is formed from rank-one terms (BilinearAdjoint).
  // AZURE_GRAD_DENSE=1 keeps the dense A^H Abar A^H products below.
  static const bool dense = [] {
    const char *e = std::getenv("AZURE_GRAD_DENSE");
    return e && *e && std::string(e) != "0";
  }();
  if (!hasEC && !dense) {
    std::vector<std::vector<BilinearTerm>> terms;
    this->CollectBilinearTerms(point, bar, terms);
    for (int j = 1; j <= (int)terms.size(); j++)
      if (!terms[j - 1].empty()) this->BilinearAdjoint(point, j, terms[j - 1], accum, shiftDeriv);
    return;
  }
  const double inEnergy = this->InEnergy(point);
  (void)exitPairNum;
  (void)nK;
  // Per-JGroup A-matrix cotangent accumulators (active == original indexing).
  std::vector<matrix_c> aBar(compound()->NumJGroups());
  std::vector<bool> aBarUsed(compound()->NumJGroups(), false);
  const complex I2 = 2.0 * I;
  auto ensureABar = [&](int jNum, int numLevels) {
    if (!aBarUsed[jNum - 1]) {
      aBar[jNum - 1].assign(numLevels, vector_c(numLevels, complex(0.0, 0.0)));
      aBarUsed[jNum - 1] = true;
    }
  };

  // === External-capture T -> parameters. ===
  // Forward:  T_ec(k,ecm) = gamma_fc * sqrtNF * conv * ecAmplitude  [* umatrix_cc],
  //   gamma_fc = final level's reduced width for the final channel (explicit),
  //   sqrtNF   = final level's 1/sqrt(NF) factor (depends on its gammas),
  //   conv     = EC conversion factor (constant),
  //   ecAmplitude treated constant w.r.t. E/gamma (energy shifts handled by FD).
  // For channel capture an extra factor umatrix_cc = sum 2i sqrtP gamma gamma A
  // over an internal pathway multiplies T_ec; it back-propagates to that
  // pathway's A-matrix and gammas just like the standard umatrix.
  if (hasEC) {
    const bool useEC = (configure().paramMask & Config::USE_EXTERNAL_CAPTURE);
    const bool ignoreZero = (configure().paramMask & Config::IGNORE_ZERO_WIDTHS);
    for (int k = 1; k <= nK; k++) {
      int numEC = theDecay->GetKGroup(k)->NumECMGroups();
      for (int m = 1; m <= numEC; m++) {
        complex eb = ecBar[k - 1][m - 1];
        if (eb == complex(0.0, 0.0)) continue;
        ECMGroup *ecg = theDecay->GetKGroup(k)->GetECMGroup(m);
        int jf = ecg->GetJGroupNum();
        int lf = ecg->GetLevelNum();
        int fc = ecg->GetFinalChannel();
        ALevel *finalLevel = compound()->GetJGroup(jf)->GetLevel(lf);
        JGroup *jgf = compound()->GetJGroup(jf);
        double sqrtNF = finalLevel->GetSqrtNFFactor();
        double conv = finalLevel->GetECConversionFactor(fc);
        complex ecAmp = useEC ? point->GetECAmplitudeWithShift(k, m, compound(), configure())
                              : point->GetECAmplitude(k, m);

        // base = gamma_fc * sqrtNF * conv * ecAmp  (the part without umatrix_cc).
        double gamma_fc = finalLevel->GetFitGamma(fc);
        complex base = gamma_fc * sqrtNF * conv * ecAmp;

        // "effective amplitude" carrying the umatrix_cc factor for the
        // gamma_fc / sqrtNF derivatives (= ecAmp for direct capture).
        complex effAmp = ecAmp;

        if (ecg->IsChannelCapture()) {
          int internalChannel = ecg->GetIntChannelNum();
          MGroup *ccMg = compound()->GetPair(aa)->GetDecay(ecg->GetChanCapDecay())->GetKGroup(ecg->GetChanCapKGroup())->GetMGroup(ecg->GetChanCapMGroup());
          int jcc = ccMg->GetJNum();
          int ccChN = ccMg->GetChNum();
          int ccChpN = ccMg->GetChpNum();
          JGroup *jgcc = compound()->GetJGroup(jcc);
          int ccLevels = jgcc->NumLevels();
          complex sqrtPcc = point->GetSqrtPenetrability(jcc, ccChN);  // single sqrtP

          // Forward umatrix_cc (mirrors the IGNORE_ZERO_WIDTHS skips).
          complex umatrix_cc(0.0, 0.0);
          for (int la = 1; la <= ccLevels; la++) {
            if (!jgcc->GetLevel(la)->IsInRMatrix()) continue;
            if (internalChannel && ignoreZero &&
                std::fabs(jgcc->GetLevel(la)->GetFitGamma(internalChannel)) < 1.0e-8) continue;
            for (int lap = 1; lap <= ccLevels; lap++) {
              if (!jgcc->GetLevel(lap)->IsInRMatrix()) continue;
              if (internalChannel && ignoreZero &&
                  std::fabs(jgcc->GetLevel(lap)->GetFitGamma(internalChannel)) < 1.0e-8) continue;
              umatrix_cc += I2 * sqrtPcc *
                  jgcc->GetLevel(la)->GetFitGamma(ccChN) *
                  jgcc->GetLevel(lap)->GetFitGamma(ccChpN) *
                  this->GetAMatrixElement(jcc, la, lap);
            }
          }
          effAmp = ecAmp * umatrix_cc;

          // umatrix_cc cotangent:  T_ec = base * umatrix_cc -> ubar = conj(base)*eb.
          complex uBarcc = conj(base) * eb;
          ensureABar(jcc, ccLevels);
          for (int la = 1; la <= ccLevels; la++) {
            if (!jgcc->GetLevel(la)->IsInRMatrix()) continue;
            if (internalChannel && ignoreZero &&
                std::fabs(jgcc->GetLevel(la)->GetFitGamma(internalChannel)) < 1.0e-8) continue;
            double gEn = jgcc->GetLevel(la)->GetFitGamma(ccChN);
            for (int lap = 1; lap <= ccLevels; lap++) {
              if (!jgcc->GetLevel(lap)->IsInRMatrix()) continue;
              if (internalChannel && ignoreZero &&
                  std::fabs(jgcc->GetLevel(lap)->GetFitGamma(internalChannel)) < 1.0e-8) continue;
              double gEx = jgcc->GetLevel(lap)->GetFitGamma(ccChpN);
              complex Aval = this->GetAMatrixElement(jcc, la, lap);
              complex w = I2 * sqrtPcc * gEn * gEx;  // weight of A[la][lap]
              aBar[jcc - 1][la - 1][lap - 1] += conj(w) * uBarcc;
              accum.AddGamma(jcc, la, ccChN, std::real(conj(uBarcc) * (I2 * sqrtPcc * gEx * Aval)));
              accum.AddGamma(jcc, lap, ccChpN, std::real(conj(uBarcc) * (I2 * sqrtPcc * gEn * Aval)));
            }
          }
        }

        // Explicit gamma_fc:  d T_ec / d gamma_fc = sqrtNF * conv * effAmp.
        accum.AddGamma(jf, lf, fc, std::real(conj(eb) * (sqrtNF * conv * effAmp)));

        // sqrtNF couples to every NF channel d of the final level:
        //   nFSum = 1 + sum_d A_d gamma_d^2,  A_d = 2 chRad redMass uconv/hbarc^2 NF_d,
        //   sqrtNF = nFSum^{-1/2},  d sqrtNF / d gamma_d = -A_d gamma_d sqrtNF^3.
        // Park: sqrtNF = (J + sum_d A_d gamma_d^2)^{-1/2} with J = 1 - sum_d gamma_d^2 S'_d
        // over the particle channels (CNuc::CalcShiftFunctions), so every particle
        // channel contributes  d sqrtNF / d gamma_d = (S'_d - A_d) gamma_d sqrtNF^3.
        complex pre = gamma_fc * conv * effAmp;  // d T_ec / d sqrtNF
        int numNF = finalLevel->NumNFIntegrals();
        double sqrtNF3 = sqrtNF * sqrtNF * sqrtNF;
        const bool parkFinal = park && finalLevel->IsInRMatrix();
        int numD = parkFinal ? jgf->NumChannels() : numNF;
        for (int d = 1; d <= numD; d++) {
          PPair *pd = compound()->GetPair(jgf->GetChannel(d)->GetPairNum());
          double Ad = (d <= numNF) ? 2.0 * pd->GetChRad() * pd->GetRedMass() * uconv / (hbarc * hbarc) * finalLevel->GetNFIntegral(d)
                                   : 0.0;
          if (parkFinal && jgf->GetChannel(d)->GetRadType() == 'P' && pd->GetPType() == 0)
            Ad -= finalLevel->GetShiftDerivative(d);
          if (Ad == 0.0) continue;
          double gd = finalLevel->GetFitGamma(d);
          double dSqrtNF = -Ad * gd * sqrtNF3;
          accum.AddGamma(jf, lf, d, std::real(conj(eb) * (pre * dSqrtNF)));
        }
      }
    }
  }

  // === T -> U -> A (+ explicit gamma), back-propagating tBar over all pathways. ===
  for (int k = 1; k <= nK; k++) {
    int numM = theDecay->GetKGroup(k)->NumMGroups();
    for (int m = 1; m <= numM; m++) {
      complex tb = tBar[k - 1][m - 1];
      if (tb == complex(0.0, 0.0)) continue;
      MGroup *mg = theDecay->GetKGroup(k)->GetMGroup(m);
      int jNum = mg->GetJNum();
      int chNum = mg->GetChNum();
      int chpNum = mg->GetChpNum();
      JGroup *jg = compound()->GetJGroup(jNum);

      complex uphase = point->GetExpCoulombPhase(jNum, chNum) *
          point->GetExpHardSpherePhase(jNum, chNum) *
          point->GetExpCoulombPhase(jNum, chpNum) *
          point->GetExpHardSpherePhase(jNum, chpNum);
      complex sqrtPenEn = point->GetSqrtPenetrability(jNum, chNum);
      complex sqrtPenEx = point->GetSqrtPenetrability(jNum, chpNum);

      // d T / d umatrix = -uphase  ->  ubar = conj(-uphase) * tBar.
      complex uBar = conj(-uphase) * tb;

      int numLevels = jg->NumLevels();
      if (!aBarUsed[jNum - 1]) {
        aBar[jNum - 1].assign(numLevels, vector_c(numLevels, complex(0.0, 0.0)));
        aBarUsed[jNum - 1] = true;
      }
      complex pen2 = 2.0 * I * sqrtPenEn * sqrtPenEx;  // common factor of w_{l,l'}

      for (int la = 1; la <= numLevels; la++) {
        if (!jg->GetLevel(la)->IsInRMatrix()) continue;
        double gammaEn = jg->GetLevel(la)->GetFitGamma(chNum);
        for (int lap = 1; lap <= numLevels; lap++) {
          if (!jg->GetLevel(lap)->IsInRMatrix()) continue;
          double gammaEx = jg->GetLevel(lap)->GetFitGamma(chpNum);
          complex Aval = this->GetAMatrixElement(jNum, la, lap);
          complex w = pen2 * gammaEn * gammaEx;  // weight of A[la][lap]

          // A cotangent:  Abar += conj(w) * ubar.
          aBar[jNum - 1][la - 1][lap - 1] += conj(w) * uBar;

          // Explicit gamma paths (umatrix depends on gamma_en and gamma_ex too):
          //   d umatrix / d gamma_en,la += pen2 * gamma_ex * A
          //   d umatrix / d gamma_ex,lap += pen2 * gamma_en * A
          accum.AddGamma(jNum, la, chNum, std::real(conj(uBar) * (pen2 * gammaEx * Aval)));
          accum.AddGamma(jNum, lap, chpNum, std::real(conj(uBar) * (pen2 * gammaEn * Aval)));
        }
      }
    }
  }

  // === A -> M lemma + structural E / gamma, once per JGroup. ===
  for (int jNum = 1; jNum <= compound()->NumJGroups(); jNum++) {
    if (!aBarUsed[jNum - 1]) continue;
    JGroup *jg = compound()->GetJGroup(jNum);
    const matrix_c &A = this->GetAMatrix(jNum);
    const matrix_c &Ab = aBar[jNum - 1];
    int n = (int)A.size();
    if (n == 0) continue;

    // Conjugate transpose A^H.
    matrix_c Ah(n, vector_c(n, complex(0.0, 0.0)));
    for (int r = 0; r < n; r++)
      for (int c = 0; c < n; c++)
        Ah[r][c] = conj(A[c][r]);

    // Mbar = -A^H * Abar * A^H.
    matrix_c tmp(n, vector_c(n, complex(0.0, 0.0)));  // tmp = Abar * A^H
    for (int r = 0; r < n; r++)
      for (int c = 0; c < n; c++) {
        complex s(0.0, 0.0);
        for (int t = 0; t < n; t++) s += Ab[r][t] * Ah[t][c];
        tmp[r][c] = s;
      }
    matrix_c Mbar(n, vector_c(n, complex(0.0, 0.0)));
    for (int r = 0; r < n; r++)
      for (int c = 0; c < n; c++) {
        complex s(0.0, 0.0);
        for (int t = 0; t < n; t++) s += Ah[r][t] * tmp[t][c];
        Mbar[r][c] = -s;
      }
    this->LevelMatrixAdjoint(point, jNum, Mbar, inEnergy, accum, shiftDeriv);
  }
}

void AMatrixFunc::TBarAdjointMany(const std::vector<EPoint *> &points, const std::vector<TMatrixCotangent> &bars,
                                  GradAccum &accum, const vector_matrix_r *shiftDeriv) {
  if (points.empty()) return;
  static const bool dense = [] {
    const char *e = std::getenv("AZURE_GRAD_DENSE");
    return e && *e && std::string(e) != "0";
  }();
  std::vector<std::vector<BilinearTerm>> merged(compound()->NumJGroups()), terms;
  for (size_t q = 0; q < points.size() && q < bars.size(); q++) {
    if (dense || !this->CollectBilinearTerms(points[q], bars[q], terms)) {
      this->TBarAdjoint(points[q], bars[q], accum, shiftDeriv);
      continue;
    }
    for (const auto &s : bars[q].screen) accum.AddScreening(s.first, s.second);
    for (int j = 1; j <= (int)terms.size(); j++)
      for (const BilinearTerm &t : terms[j - 1]) {
        bool found = false;
        for (BilinearTerm &m : merged[j - 1])
          if (m.ch == t.ch && m.chp == t.chp) {
            m.beta += t.beta;
            found = true;
            break;
          }
        if (!found) merged[j - 1].push_back(t);
      }
  }
  AZUREProfile::Timer profileTimer(AZUREProfile::kGradBackprop);
  for (int j = 1; j <= (int)merged.size(); j++)
    if (!merged[j - 1].empty()) this->BilinearAdjoint(points[0], j, merged[j - 1], accum, shiftDeriv);
}

/*!
 * Beam energy in the compound system, as FillMatrices takes it.
 */
double AMatrixFunc::InEnergy(EPoint *point) {
  if (compound()->GetPair(compound()->GetPairNumFromKey(point->GetEntranceKey()))->GetPType() == 20) {
    int ep = compound()->GetPairNumFromKey(point->GetExitKey());
    return point->GetCMEnergy() + compound()->GetPair(ep)->GetSepE() + compound()->GetPair(ep)->GetExE();
  }
  int aa = compound()->GetPairNumFromKey(point->GetEntranceKey());
  return point->GetCMEnergy() + compound()->GetPair(aa)->GetSepE() + compound()->GetPair(aa)->GetExE();
}

/*!
 * The level-matrix part of the adjoint of J-group jNum: given the cotangent
 * Mbar of its level matrix, the structural energy and reduced-width gradients
 * of the standard, Brune and Park level matrices (exactly mirroring
 * FillMatrices).
 */
void AMatrixFunc::LevelMatrixAdjoint(EPoint *energyPoint, int jNum, const matrix_c &Mbar, double inEnergy,
                                     GradAccum &accum, const vector_matrix_r *shiftDeriv, const complex *lambda) {
  AZUREProfile::Timer profileTimer(AZUREProfile::kGradLevelMatrix);
  EPoint *point = energyPoint;
  const bool brune = (configure().paramMask & Config::USE_BRUNE_FORMALISM);
  const bool park = (configure().paramMask & Config::USE_PARK_FORMALISM);
  JGroup *jg = compound()->GetJGroup(jNum);
  const int n = (int)Mbar.size();
  if (n == 0) return;
  // Structural energy gradient.  The (E_r - inE) delta term gives
  //   d M[a][a]/d E_a = 1  ->  dlnL/dE_a += Re(Mbar[a][a]).
  for (int la = 1; la <= n; la++)
    accum.AddE(jNum, la, std::real(Mbar[la - 1][la - 1]));

  // Brune adds explicit E_lambda dependence through the boundary/shift kernel
  //   M[r][c] = ... - sum_{ch P} g_r,ch g_c,ch K_ch[r][c],
  //   K_ch[r][c] = L + B - D_ch[r][c],  D diagonal = S_r,  off-diag = Q_{rc}.
  // so dlnL/dE_nu += sum_{r,c} Re(conj(Mbar[r][c]) * sum_ch g_r g_c dD/dE_nu).
  // dD/dE only couples to E_r and E_c; using S'(E) = shiftDeriv.
  if (park) {
    // Park:  M[r][c] = (E_r - inE) J[r][c] - sum_d g_r,d g_c,d K_d[r][c],
    //   K_d[r][c] = L_d + B_d - S_r,d   (row-indexed shift; particle channels),
    //   J[r][r]   = 1 - sum_d g_r,d^2 S'_r,d,
    //   J[r][c]   = -sum_d g_r,d g_c,d (S_r,d - S_c,d)/(E_r - E_c).
    // The delta term above already gave d M[r][r]/dE_r = 1; here it is J_rr
    // instead, so the difference J_rr - 1 is added, plus the J and S_r,d
    // derivatives.  M is real-coefficient in E, so only Re(Mbar) enters.
    int numChannels = jg->NumChannels();
    auto isParticle = [&](int d) {
      return jg->GetChannel(d)->GetRadType() == 'P' &&
             compound()->GetPair(jg->GetChannel(d)->GetPairNum())->GetPType() == 0;
    };
    for (int r = 1; r <= n; r++) {
      ALevel *lr = jg->GetLevel(r);
      double Er = lr->GetFitE();
      for (int c = 1; c <= n; c++) {
        double reM = std::real(Mbar[r - 1][c - 1]);
        if (reM == 0.0) continue;
        if (r == c) {
          double dJ = 0.0, sumS1 = 0.0;
          for (int d = 1; d <= numChannels; d++) {
            if (!isParticle(d)) continue;
            double g2 = lr->GetFitGamma(d) * lr->GetFitGamma(d);
            if (g2 == 0.0) continue;
            dJ -= g2 * lr->GetShiftSecondDerivative(d);
            sumS1 += g2 * lr->GetShiftDerivative(d);
          }
          accum.AddE(jNum, r, reM * ((lr->GetParkNorm() - 1.0) + (Er - inEnergy) * dJ + sumS1));
        } else {
          ALevel *lc = jg->GetLevel(c);
          double Ec = lc->GetFitE();
          double den = Er - Ec;
          double Jrc = 0.0, dJdEr = 0.0, dJdEc = 0.0, sumS1 = 0.0;
          for (int d = 1; d <= numChannels; d++) {
            if (!isParticle(d)) continue;
            double g = lr->GetFitGamma(d) * lc->GetFitGamma(d);
            if (g == 0.0) continue;
            double Sr = lr->GetShiftFunction(d), Sc = lc->GetShiftFunction(d);
            double dSr = lr->GetShiftDerivative(d), dSc = lc->GetShiftDerivative(d);
            double q = (Sr - Sc) / den;
            Jrc -= g * q;
            dJdEr -= g * (dSr - q) / den;
            dJdEc -= g * (-dSc + q) / den;
            sumS1 += g * dSr;
          }
          accum.AddE(jNum, r, reM * (Jrc + (Er - inEnergy) * dJdEr + sumS1));
          accum.AddE(jNum, c, reM * (Er - inEnergy) * dJdEc);
        }
      }
    }

    // Structural gamma gradient, Park.  For M[a][c] and M[c][a] (both carry
    // Mbar) the kernel is row-indexed, K_d[a][c] = L_d + B_d - S_a,d and
    // K_d[c][a] = L_d + B_d - S_c,d, and J contributes
    //   d J[a][a]/d g_a,d = -2 g_a,d S'_a,d,
    //   d J[a][c]/d g_a,d = d J[c][a]/d g_a,d = -g_c,d (S_a,d - S_c,d)/(E_a - E_c).
    for (int d = 1; d <= numChannels; d++) {
      const bool pd = isParticle(d);
      // lambda[d-1] stands for L_d + B_d (particle) or L_d.
      complex Ld = lambda ? lambda[d - 1] : point->GetLoElement(jNum, d);
      double Bd = (pd && !lambda) ? jg->GetChannel(d)->GetBoundaryCondition() : 0.0;
      for (int a = 1; a <= n; a++) {
        ALevel *la_ = jg->GetLevel(a);
        double Ea = la_->GetFitE();
        double ga = la_->GetFitGamma(d);
        double Sa = pd ? la_->GetShiftFunction(d) : 0.0;
        double grad = 0.0;
        for (int c = 1; c <= n; c++) {
          ALevel *lc_ = jg->GetLevel(c);
          double gc = lc_->GetFitGamma(d);
          if (c == a) {
            complex Kaa = pd ? (Ld + Bd - Sa) : Ld;
            grad += -2.0 * ga * std::real(conj(Mbar[a - 1][a - 1]) * Kaa);
            if (pd) grad += std::real(Mbar[a - 1][a - 1]) * (Ea - inEnergy) * (-2.0 * ga * la_->GetShiftDerivative(d));
          } else {
            if (gc == 0.0) continue;
            double Sc = pd ? lc_->GetShiftFunction(d) : 0.0;
            complex Kac = pd ? (Ld + Bd - Sa) : Ld;
            complex Kca = pd ? (Ld + Bd - Sc) : Ld;
            grad += -gc * std::real(conj(Mbar[a - 1][c - 1]) * Kac + conj(Mbar[c - 1][a - 1]) * Kca);
            if (pd) {
              double Ec = lc_->GetFitE();
              double delta = (Sa - Sc) / (Ea - Ec);
              grad += -gc * delta * (std::real(Mbar[a - 1][c - 1]) * (Ea - inEnergy) +
                                     std::real(Mbar[c - 1][a - 1]) * (Ec - inEnergy));
            }
          }
        }
        accum.AddGamma(jNum, a, d, grad);
      }
    }
    return;
  }

  // Standard and Brune level matrices.  The level data are gathered once into
  // flat arrays: the loops below are n^2 per channel, and the accessor calls
  // and divisions they used to make per element dominated the gradient.
  const int nch = jg->NumChannels();
  // E, G, S, S G and 1/(E_r - E_c) off the diagonal, kept per parameter set.
  const ChannelSpace::LevelArrays &levelArrays = ChannelSpace::CurrentLevelArrays(jg);
  const std::vector<double> &E = levelArrays.E, &G = levelArrays.G, &S = levelArrays.S, &inv = levelArrays.inv;
  std::vector<double> &reSym = lma_resym_, &imSym = lma_imsym_;
  std::vector<char> &bruneP = lma_brunep_;
  bruneP.assign(nch, 0);
  bool anyBruneP = false;
  for (int ch = 0; ch < nch; ch++)
    if (brune && jg->GetChannel(ch + 1)->GetRadType() == 'P') bruneP[ch] = anyBruneP = true;

  // Brune adds explicit E_lambda dependence through the boundary/shift kernel
  //   M[r][c] = ... - sum_{ch P} g_r,ch g_c,ch K_ch[r][c],
  //   K_ch[r][c] = L + B - D_ch[r][c],  D diagonal = S_r,  off-diag = Q_{rc},
  //   Q_rc = (S_r (inE - E_c) - S_c (inE - E_r)) / (E_r - E_c),
  // so dlnL/dE_nu += sum_{r,c} Re(Mbar[r][c]) sum_ch g_r g_c dD/dE_nu; dD/dE
  // only couples to E_r and E_c, with S'(E) = shiftDeriv.
  if (brune && shiftDeriv && anyBruneP) {
    const std::vector<std::vector<double>> &sd = (*shiftDeriv)[jNum - 1];
    std::vector<double> &dE = lma_de_;
    dE.assign(n, 0.0);
    for (int ch = 0; ch < nch; ch++) {
      if (!bruneP[ch]) continue;
      for (int r = 0; r < n; r++) {
        const double gr = G[(size_t)r * nch + ch];
        if (gr == 0.0) continue;
        const double Sr = S[(size_t)r * nch + ch];
        const double dSr = sd[r][ch];
        const double ar = inEnergy - E[r];
        const complex *mrow = &Mbar[r][0];
        // D = S_r ;  dD/dE_r = S'_r.
        dE[r] += std::real(mrow[r]) * gr * gr * dSr;
        for (int c = 0; c < n; c++) {
          if (c == r) continue;
          const double g = gr * G[(size_t)c * nch + ch];
          if (g == 0.0) continue;
          const double Sc = S[(size_t)c * nch + ch];
          const double dSc = sd[c][ch];
          const double ac = inEnergy - E[c];
          const double iv = inv[(size_t)r * n + c];
          const double Q = (Sr * ac - Sc * ar) * iv;
          const double dQdEr = (dSr * ac + Sc - Q) * iv;
          const double dQdEc = (-Sr - dSc * ar + Q) * iv;
          const double w = std::real(mrow[c]) * g;
          dE[r] += w * dQdEr;
          dE[c] += w * dQdEc;
        }
      }
    }
    for (int r = 0; r < n; r++) accum.AddE(jNum, r + 1, dE[r]);
  }

  // Structural gamma gradient.  M[r][c] = (E_r-inE) delta - sum_ch g_r,ch g_c,ch K_ch[r][c],
  // so  dlnL/d gamma_{a,d} += -sum_c Re( conj(Mbar[a][c]+Mbar[c][a]) * K_d[a][c] ) * gamma_{c,d}.
  // The kernel K_d[r][c] is the channel-d contribution to the level matrix:
  //   non-Brune (or non-'P'):  K_d = L_d  (the Lo-matrix element, indep. of r,c)
  //   Brune 'P' channel:       K_d[r][c] = L_d + B_d - D_d[r][c],
  //     D_d[r][c] = (r==c) ? S_r,d
  //                        : (S_r,d (inE-E_c) - S_c,d (inE-E_r)) / (E_r - E_c)
  // exactly mirroring FillMatrices.  With sym = Mbar + Mbar^T and D real,
  //   -Re(conj(sym) K) = -Re(sym) Re(L + B) - Im(sym) Im(L) + Re(sym) D,
  // and the D term splits into
  //   S_a sum_{c!=a} Re(sym_ac) (inE-E_c)/(E_a-E_c) g_c
  //   - (inE-E_a) sum_{c!=a} Re(sym_ac)/(E_a-E_c) S_c g_c  + Re(sym_aa) S_a g_a,
  // so every channel is a few row-times-matrix products.
  reSym.resize((size_t)n * n);
  imSym.resize((size_t)n * n);
  for (int a = 0; a < n; a++)
    for (int c = 0; c < n; c++) {
      const complex sym = Mbar[a][c] + Mbar[c][a];
      reSym[(size_t)a * n + c] = std::real(sym);
      imSym[(size_t)a * n + c] = std::imag(sym);
    }
  const std::vector<double> &SG = levelArrays.SG;
  std::vector<double> &acc = lma_acc_;
  acc.resize((size_t)4 * nch);
  double *sRe = &acc[0], *sIm = &acc[nch], *t1 = &acc[2 * (size_t)nch], *t2 = &acc[3 * (size_t)nch];
  std::vector<complex> &Lconst = lma_l_;
  Lconst.resize(nch);
  for (int d = 0; d < nch; d++) {
    if (lambda) {
      Lconst[d] = lambda[d];
      continue;
    }
    complex Ld = point->GetLoElement(jNum, d + 1);
    if (bruneP[d]) Ld += jg->GetChannel(d + 1)->GetBoundaryCondition();
    Lconst[d] = Ld;
  }
  for (int a = 0; a < n; a++) {
    std::fill(acc.begin(), acc.end(), 0.0);
    const double *rs = &reSym[(size_t)a * n], *is = &imSym[(size_t)a * n];
    for (int c = 0; c < n; c++) {
      const double *gc = &G[(size_t)c * nch];
      const double wr = rs[c], wi = is[c];
      for (int d = 0; d < nch; d++) {
        sRe[d] += wr * gc[d];
        sIm[d] += wi * gc[d];
      }
      if (anyBruneP && c != a) {
        const double iv = inv[(size_t)a * n + c];
        const double w1 = wr * (inEnergy - E[c]) * iv, w2 = wr * iv;
        const double *sgc = &SG[(size_t)c * nch];
        for (int d = 0; d < nch; d++) {
          t1[d] += w1 * gc[d];
          t2[d] += w2 * sgc[d];
        }
      }
    }
    for (int d = 0; d < nch; d++) {
      double g = -(std::real(Lconst[d]) * sRe[d] + std::imag(Lconst[d]) * sIm[d]);
      if (bruneP[d]) {
        const double Sa = S[(size_t)a * nch + d];
        g += Sa * t1[d] - (inEnergy - E[a]) * t2[d] + rs[a] * Sa * G[(size_t)a * nch + d];
      }
      accum.AddGamma(jNum, a + 1, d + 1, g);
    }
  }
}

/*!
 * A times the reduced-width vector of one channel, the widths taken as they are
 * (the forward bilinear forms drop widths below 1e-12, GetGammaVector; the
 * adjoint, like the dense one, differentiates every width).  Cached for the
 * current energy.
 */
const vector_c &AMatrixFunc::GetRawChannelSolve(int jGroupNum, int chNum) {
  vector_c &solve = raw_solves_[jGroupNum - 1][chNum - 1];
  if (!raw_solve_valid_[jGroupNum - 1][chNum - 1]) {
    AZUREProfile::Timer profileTimer(AZUREProfile::kChannelSolve);
    JGroup *jGroup = compound()->GetJGroup(jGroupNum);
    const std::vector<int> &map = level_active_index_[jGroupNum - 1];
    const int n = (int)active_levels_[jGroupNum - 1].size();
    solve.assign(n, complex(0.0, 0.0));
    for (int la = 1; la < (int)map.size(); la++) {
      const int row = map[la];
      if (row == 0) continue;
      solve[row - 1] = jGroup->GetLevel(la)->GetFitGamma(chNum);
    }
    if (n > 0) {
      if (route_[jGroupNum - 1] == kRouteChannel) {
        // The widths as they are equal the ones the channel route was built
        // from unless one is below the 1e-12 cutoff: A g_c then comes from the
        // channel solves already at hand.
        const ChannelSpace::FormalLevels &f = (*formal_[jGroupNum - 1]);
        bool same = true;
        for (int k = 0; k < n && same; k++) same = (std::real(solve[k]) == f.gcut[(size_t)k * f.nch + chNum - 1]);
        if (same && !f.pencil) {
          const vector_c &uf = this->FormalSolve(jGroupNum, chNum - 1);
          if (f.identity)
            solve = uf;
          else
            for (int r = 0; r < n; r++) {
              const double *vr = &f.V[(size_t)r * n];
              complex sum(0.0, 0.0);
              for (int k = 0; k < n; k++) sum += vr[k] * uf[k];
              solve[r] = sum;
            }
        } else if (same) {
          // Pencil: A g_c = Z X (e_c + z),  (1 - Lambda R) z = Lambda R e_c.
          ChannelState &st = chan_[jGroupNum - 1];
          const int m = (int)f.act.size(), b = f.pos[chNum - 1];
          std::fill(solve.begin(), solve.end(), complex(0.0, 0.0));
          if (b >= 0) {
            vector_c &z = st.wz;
            z.resize(m);
            for (int a = 0; a < m; a++) z[a] = st.lamAct[a] * st.R[(size_t)a * m + b];
            this->SolveBordered(jGroupNum, nullptr, z.data());
            z[b] += 1.0;
            vector_c &y = st.w4;
            y.resize(n);
            for (int k = 0; k < n; k++) {
              complex sum(0.0, 0.0);
              for (int a = 0; a < m; a++) sum += st.X[(size_t)k * m + a] * z[a];
              y[k] = sum;
            }
            for (int r = 0; r < n; r++) {
              const double *zr = &f.Z[(size_t)r * n];
              complex sum(0.0, 0.0);
              for (int k = 0; k < n; k++) sum += zr[k] * y[k];
              solve[r] = sum;
            }
          }
        } else {
          vector_c &tmp = chan_[jGroupNum - 1].w3;
          tmp = solve;
          this->ApplyA(jGroupNum, tmp.data(), solve.data());
        }
      } else
        solvers_[jGroupNum - 1].Solve(&solve[0]);
    }
    raw_solve_valid_[jGroupNum - 1][chNum - 1] = 1;
  }
  return solve;
}

/*!
 * The channel kernel Lambda_d of J-group j at the energy of `point`, as
 * LevelMatrixAdjoint takes it: the Lo-matrix element, plus the boundary
 * condition for the particle channels of the Brune and Park level matrices.
 */
const std::vector<complex> &AMatrixFunc::KernelLambda(EPoint *point, int j) {
  JGroup *jg = compound()->GetJGroup(j);
  const bool park = (configure().paramMask & Config::USE_PARK_FORMALISM);
  const bool brune = (configure().paramMask & Config::USE_BRUNE_FORMALISM);
  const int nch = jg->NumChannels();
  kernel_scratch_.resize(nch);
  for (int d = 1; d <= nch; d++) {
    AChannel *ch = jg->GetChannel(d);
    complex l = point->GetLoElement(j, d);
    const bool particle = park ? (ch->GetRadType() == 'P' && compound()->GetPair(ch->GetPairNum())->GetPType() == 0)
                               : (brune && ch->GetRadType() == 'P');
    if (particle) l += ch->GetBoundaryCondition();
    kernel_scratch_[d - 1] = l;
  }
  return kernel_scratch_;
}

/*!
 * T -> bilinear forms.  Each internal pathway T(k,m) of J-group j is
 *   T = [tphase] - uphase (1 + 2i sqrtP_en sqrtP_ex b),   b = g_en^T A g_ex,
 * so its cotangent tbar gives the bilinear form the cotangent
 *   beta = conj(2i sqrtP_en sqrtP_ex) conj(-uphase) tbar.
 * Pathways of one J-group sharing a channel pair add their betas.
 */
bool AMatrixFunc::CollectBilinearTerms(EPoint *energyPoint, const TMatrixCotangent &bar,
                                       std::vector<std::vector<BilinearTerm>> &terms) {
  EPoint *point = energyPoint;
  const complex I(0.0, 1.0);
  int aa = compound()->GetPairNumFromKey(point->GetEntranceKey());
  int exitPairNum = compound()->GetPairNumFromKey(point->GetExitKey());
  int ir = 0;
  while (ir < compound()->GetPair(aa)->NumDecays()) {
    ir++;
    if (compound()->GetPair(aa)->GetDecay(ir)->GetPairNum() == exitPairNum) break;
  }
  Decay *theDecay = compound()->GetPair(aa)->GetDecay(ir);
  terms.assign(compound()->NumJGroups(), std::vector<BilinearTerm>());
  for (int k = 1; k <= theDecay->NumKGroups(); k++)
    if (theDecay->GetKGroup(k)->NumECMGroups() > 0) return false;
  for (int k = 1; k <= theDecay->NumKGroups() && k <= (int)bar.t.size(); k++) {
    const int numM = theDecay->GetKGroup(k)->NumMGroups();
    for (int m = 1; m <= numM && m <= (int)bar.t[k - 1].size(); m++) {
      const complex tb = bar.t[k - 1][m - 1];
      if (tb == complex(0.0, 0.0)) continue;
      MGroup *mg = theDecay->GetKGroup(k)->GetMGroup(m);
      const int jNum = mg->GetJNum();
      const int chNum = mg->GetChNum();
      const int chpNum = mg->GetChpNum();
      complex uphase = point->GetExpCoulombPhase(jNum, chNum) * point->GetExpHardSpherePhase(jNum, chNum) *
          point->GetExpCoulombPhase(jNum, chpNum) * point->GetExpHardSpherePhase(jNum, chpNum);
      complex uBar = conj(-uphase) * tb;
      complex pen2 = 2.0 * I * point->GetSqrtPenetrability(jNum, chNum) * point->GetSqrtPenetrability(jNum, chpNum);
      const complex beta = conj(pen2) * uBar;
      std::vector<BilinearTerm> &list = terms[jNum - 1];
      bool merged = false;
      for (BilinearTerm &t : list)
        if (t.ch == chNum && t.chp == chpNum) {
          t.beta += beta;
          merged = true;
          break;
        }
      if (!merged) list.push_back(BilinearTerm{chNum, chpNum, beta});
    }
  }
  return true;
}

/*!
 * Adjoint of sum_t beta_t g_ch^T A g_chp over the terms of J-group jNum.
 *
 * Explicit widths:  d b / d g_{la,ch} = (A g_chp)_la,  d b / d g_{la,chp} = (A^T g_ch)_la,
 * and A is complex symmetric (the level matrix is), so both come from the
 * solves u_c = A g_c.  Level matrix:  A = M^{-1}, so with Abar = sum_t beta_t
 * g_ch g_chp^T,
 *   Mbar = -A^H Abar A^H = -sum_t beta_t conj(u_ch) conj(u_chp)^T,
 * the same matrix the dense path forms with two n^3 products; its structural
 * part is LevelMatrixAdjoint.
 */
void AMatrixFunc::BilinearAdjoint(EPoint *energyPoint, int jNum, const std::vector<BilinearTerm> &terms,
                                  GradAccum &accum, const vector_matrix_r *shiftDeriv) {
  if (terms.empty()) return;
  if (jNum < 1 || jNum > (int)solvers_.size()) return;
  const bool channel = (route_[jNum - 1] == kRouteChannel);
  if (!channel && !solver_valid_[jNum - 1]) return;
  const int n = (int)active_levels_[jNum - 1].size();
  if (n == 0) return;
  // In channel space the adjoint is taken in the formal representation: for
  // the standard level matrix that is the level matrix itself; for Brune and
  // Park only a scalar gradient pass can use it (its formal cotangents are
  // resolved once at the end of the pass, ResolveFormalCotangents).
  if (channel) {
    const ChannelSpace::FormalLevels &f = (*formal_[jNum - 1]);
    if (f.identity || (accum.formal && f.adjointOK)) {
      this->FormalBilinearAdjoint(jNum, terms, accum);
      return;
    }
  }
  // All levels are in the R matrix on this path (PointTBar), so active index
  // = level index.
  if (accum.formal) {
    // A scalar pass: Mbar = -sum_t beta_t conj(u_ch) conj(u_chp)^T enters the
    // structural adjoint through K - E N, linearly in E, so only Re(Mbar) and
    // E Re(Mbar) are summed here (resolved once per pass), and the kernel
    // G Lambda(E) G^T is differentiated now, through the rank-one terms:
    //   gammabar_{l d} -= Re(conj(Lambda_d) ((Mbar + Mbar^T) g_d)_l).
    JGroup *jg = compound()->GetJGroup(jNum);
    const int nch = jg->NumChannels();
    const double inEnergy = this->InEnergy(energyPoint);
    std::vector<double> &K = accum.lmK[jNum - 1], &N = accum.lmN[jNum - 1];
    if (!accum.lmUsed[jNum - 1]) {
      K.assign((size_t)n * n, 0.0);
      N.assign((size_t)n * n, 0.0);
      accum.lmUsed[jNum - 1] = 1;
    }
    accum.formalUsed = true;
    const std::vector<complex> &lam = this->KernelLambda(energyPoint, jNum);
    const std::vector<double> &G = ChannelSpace::CurrentLevelArrays(jg).G;
    vector_c ga(nch), gb(nch);
    for (const BilinearTerm &t : terms) {
      const vector_c &uEn = this->GetRawChannelSolve(jNum, t.ch);
      const vector_c &uEx = this->GetRawChannelSolve(jNum, t.chp);
      const complex cb = conj(t.beta);
      for (int la = 1; la <= n; la++) {
        accum.AddGamma(jNum, la, t.ch, std::real(cb * uEx[la - 1]));
        accum.AddGamma(jNum, la, t.chp, std::real(cb * uEn[la - 1]));
      }
      // conj(u)^T g_d for both solves.
      std::fill(ga.begin(), ga.end(), complex(0.0, 0.0));
      std::fill(gb.begin(), gb.end(), complex(0.0, 0.0));
      for (int la = 0; la < n; la++) {
        const complex ca = conj(uEn[la]), cbx = conj(uEx[la]);
        const double *g = &G[(size_t)la * nch];
        for (int d = 0; d < nch; d++) {
          ga[d] += ca * g[d];
          gb[d] += cbx * g[d];
        }
      }
      for (int la = 0; la < n; la++) {
        const complex a = -t.beta * conj(uEn[la]), b = -t.beta * conj(uEx[la]);
        // (Mbar g_d)_l = a_l (conj(uEx) . g_d);  (Mbar^T g_d)_l = b_l (conj(uEn) . g_d).
        for (int d = 0; d < nch; d++)
          accum.gamma[jNum - 1][la][d] -= std::real(conj(lam[d]) * (a * gb[d] + b * ga[d]));
        double *kr = &K[(size_t)la * n], *nr = &N[(size_t)la * n];
        for (int c = 0; c < n; c++) {
          const double v = std::real(a * conj(uEx[c]));
          kr[c] += v;
          nr[c] += inEnergy * v;
        }
      }
    }
    return;
  }
  // A row of the residual Jacobian or of the model gradients: the rank-two
  // derivatives of the level matrix by free parameter (LevelAdjoint).
  if (accum.mask && this->LevelAdjoint(energyPoint, jNum, terms, *accum.mask, shiftDeriv,
                                       [&](int l, int ch, double v) {
                                         if (ch < 0)
                                           accum.e[jNum - 1][l] += v;
                                         else
                                           accum.gamma[jNum - 1][l][ch] += v;
                                       }))
    return;
  matrix_c &Mbar = mbar_scratch_;
  Mbar.resize(n);
  for (int r = 0; r < n; r++) Mbar[r].assign(n, complex(0.0, 0.0));
  for (const BilinearTerm &t : terms) {
    const vector_c &uEn = this->GetRawChannelSolve(jNum, t.ch);
    const vector_c &uEx = this->GetRawChannelSolve(jNum, t.chp);
    const complex cb = conj(t.beta);
    for (int la = 1; la <= n; la++) {
      accum.AddGamma(jNum, la, t.ch, std::real(cb * uEx[la - 1]));
      accum.AddGamma(jNum, la, t.chp, std::real(cb * uEn[la - 1]));
    }
    for (int r = 0; r < n; r++) {
      const complex a = -t.beta * conj(uEn[r]);
      vector_c &row = Mbar[r];
      for (int c = 0; c < n; c++) row[c] += a * conj(uEx[c]);
    }
  }
  this->LevelMatrixAdjoint(energyPoint, jNum, Mbar, this->InEnergy(energyPoint), accum, shiftDeriv);
}

/*!
 * One diagonal bilinear form g_ch^T A g_ch with cotangent mBar (the scattering length's M_cc).  The
 * general BilinearAdjoint takes every level of the group as active, so groups with a level left out
 * of the R-matrix are refused here, as is the RMC formalism.
 */
bool AMatrixFunc::BilinearAdjoint(EPoint *point, int jGroupNum, int chNum, complex mBar, GradAccum &accum,
                                  const vector_matrix_r *shiftDeriv) {
  if (configure().paramMask & Config::USE_RMC_FORMALISM) return false;
  if (jGroupNum < 1 || jGroupNum > (int)solvers_.size()) return false;
  // Solved in channel space (FillMatrices) or as a level matrix (InvertMatrices).
  if (route_[jGroupNum - 1] != kRouteChannel && !solver_valid_[jGroupNum - 1]) return false;
  JGroup *jg = compound()->GetJGroup(jGroupNum);
  if (!jg->IsInRMatrix()) return false;
  for (int la = 1; la <= jg->NumLevels(); la++)
    if (!jg->GetLevel(la)->IsInRMatrix()) return false;
  this->BilinearAdjoint(point, jGroupNum, std::vector<BilinearTerm>{BilinearTerm{chNum, chNum, mBar}}, accum,
                        shiftDeriv);
  return true;
}

/*!
 * The per-energy part of the channel route for J-group j: the kernel
 * \f$\Lambda\f$, \f$D^{-1}\f$, \f$R = \Gamma^T D^{-1}\Gamma\f$ over the coupled channels and
 * the factors of \f$1 - R\Lambda\f$.  Returns false (the level route is then
 * taken) near a formal pole, where \f$R\f$ is dominated by one level and the
 * Woodbury form loses digits the level matrix keeps, or when \f$1 - R\Lambda\f$
 * is ill-conditioned.
 */
bool AMatrixFunc::SolveChannelSpace(int j, EPoint *point, double inEnergy) {
  const ChannelSpace::FormalLevels &f = (*formal_[j - 1]);
  ChannelState &st = chan_[j - 1];
  const int n = f.n, nch = f.nch, m = (int)f.act.size();
  st.inverseValid = false;
  st.lam.resize(nch);
  for (int c = 0; c < nch; c++) {
    complex l = point->GetLoElement(j, c + 1);
    l += f.lamShift[c];
    st.lam[c] = l;
  }
  st.lamAct.resize(m);
  for (int a = 0; a < m; a++) st.lamAct[a] = st.lam[f.act[a]];
  st.defl.clear();
  if (f.pencil) {
    if (!this->SolvePencil(j, inEnergy)) return false;
  } else if (!this->SolveFormal(j, inEnergy)) return false;
  // The bordered system of the deflated levels p and the channels:
  //   S = [ D_p        -Gamma_p    ]
  //       [ -Lambda Gamma_p^T   1 - Lambda R ]
  // (just 1 - Lambda R without deflated levels).
  const int kd = (int)st.defl.size(), ns = kd + m;
  st.lu.assign((size_t)ns * ns, complex(0.0, 0.0));
  for (int p = 0; p < kd; p++) {
    const int lv = st.defl[p];
    st.lu[(size_t)p * ns + p] = (f.e[lv] - inEnergy) + f.elo[lv];
    for (int a = 0; a < m; a++) {
      const double g = f.gamma[(size_t)lv * m + a];
      st.lu[(size_t)p * ns + kd + a] = -g;
      st.lu[(size_t)(kd + a) * ns + p] = -st.lamAct[a] * g;
    }
  }
  for (int a = 0; a < m; a++)
    for (int b = 0; b < m; b++)
      st.lu[(size_t)(kd + a) * ns + kd + b] =
          (a == b ? complex(1.0, 0.0) : complex(0.0, 0.0)) - st.lamAct[a] * st.R[(size_t)a * m + b];
  st.piv.resize(ns);
  if (ns > 0 && !(ChannelSpace::LUFactor(st.lu.data(), st.piv.data(), ns) >= ChannelSpace::kMinPivotRatio)) {
    AZUREProfile::Count(AZUREProfile::kLevelRoutePivot);
    return false;
  }
  st.F.resize((size_t)m * m);
  st.fcol.assign(m, 0);
  st.usolve.resize(nch);
  st.uvalid.assign(nch, 0);
  return true;
}

/*!
 * R of the pencil route: \f$R = G^T (K - E N)^{-1} G\f$ from a real LU of the
 * pencil.  Refused near a formal pole, where (K - E N) is nearly singular: the
 * estimated relative error of R (round-off times the condition of the pencil)
 * times the size of \f$R\Lambda\f$ must stay below kMaxPoleError.
 */
bool AMatrixFunc::SolvePencil(int j, double inEnergy) {
  const ChannelSpace::FormalLevels &f = (*formal_[j - 1]);
  ChannelState &st = chan_[j - 1];
  const int n = f.n, m = (int)f.act.size();
  // H - E T, formed in extended precision: H and E T are much larger than
  // their difference near a formal pole.  Upper Hessenberg.
  st.P.assign((size_t)n * n, 0.0);
  const long double e = inEnergy;
  double pmax = 0.0;
  for (int i = 0; i < n; i++)
    for (int k = std::max(0, i - 1); k < n; k++) {
      const double v = (double)(f.H[(size_t)i * n + k] - e * f.T[(size_t)i * n + k]);
      st.P[(size_t)i * n + k] = v;
      pmax = std::max(pmax, std::fabs(v));
    }
  st.hswap.resize(n);
  const double ratio = ChannelSpace::HessenbergLU(st.P.data(), st.hswap.data(), n);
  if (!(ratio > 0.0)) {
    AZUREProfile::Count(AZUREProfile::kLevelRoutePole);
    return false;
  }
  // X = (H - E T)^-1 Qt G;  R = (Z^T G)^T X  (= G^T (K - E N)^-1 G).
  st.X.assign(f.gq.begin(), f.gq.end());
  ChannelSpace::HessenbergSolve(st.P.data(), st.hswap.data(), n, st.X.data(), m);
  st.R.assign((size_t)m * m, complex(0.0, 0.0));
  double rl = 0.0, x2max = 0.0;
  for (int a = 0; a < m; a++) {
    double x2 = 0.0;
    for (int k = 0; k < n; k++) x2 += st.X[(size_t)k * m + a] * st.X[(size_t)k * m + a];
    x2max = std::max(x2max, x2 * std::abs(st.lamAct[a]));
    for (int b = a; b < m; b++) {
      double sab = 0.0, sba = 0.0;
      for (int k = 0; k < n; k++) {
        sab += f.gz[(size_t)k * m + a] * st.X[(size_t)k * m + b];
        sba += f.gz[(size_t)k * m + b] * st.X[(size_t)k * m + a];
      }
      const double r = 0.5 * (sab + sba);
      st.R[(size_t)a * m + b] = st.R[(size_t)b * m + a] = r;
      rl = std::max(rl, std::fabs(r) * std::max(std::abs(st.lamAct[a]), std::abs(st.lamAct[b])));
    }
  }
  // Round-off of R Lambda: the pencil is known to a few ulps of its entries,
  // dP, and dR = -G^T P^-1 dP P^-1 G, so the column of channel a of R Lambda
  // is uncertain by about ulps eps |P| |P^-1 G_a|^2 |Lambda_a| -- the analogue
  // of eta on the formal route (P^-1 G = Z X).
  const double err = ChannelSpace::kEnergyUlps * std::numeric_limits<double>::epsilon() * pmax * x2max;
  if (!(rl <= ChannelSpace::kMaxPoleRatio) || !(err <= ChannelSpace::kMaxPoleError)) {
    AZUREProfile::Count(AZUREProfile::kLevelRoutePole);
    return false;
  }
  return true;
}

/*!
 * R of the formal route: \f$R = \Gamma^T D^{-1} \Gamma\f$.
 */
bool AMatrixFunc::SolveFormal(int j, double inEnergy) {
  const ChannelSpace::FormalLevels &f = (*formal_[j - 1]);
  ChannelState &st = chan_[j - 1];
  const int n = f.n, m = (int)f.act.size();
  st.dinv.resize(n);
  st.isDefl.assign(n, 0);
  // Levels near a pole are deflated: kept as unknowns of the bordered system
  // (SolveChannelSpace) rather than summed into R.  Level k contributes
  // w_k / (e_k - E) to R, w_k = Gamma_k Gamma_k^T, which the kernel turns into
  // terms of size s_k = |Gamma_k|^2_Lambda / |e_k - E|.  Near a formal pole
  // two things would go wrong that the level matrix does not suffer from: the
  // large terms cancel in (1 - Lambda R)^-1 (kappa_k = s_k), and the absolute
  // error of e_k - E (round-off of E and of the formal energy e_k, a few ulps
  // of the energies involved) becomes a relative error of s_k
  // (eta_k = s_k de / |e_k - E|).
  const double de = ChannelSpace::kEnergyUlps * std::numeric_limits<double>::epsilon() * (f.emax + std::fabs(inEnergy));
  for (int k = 0; k < n; k++) {
    const complex d = (f.e[k] - inEnergy) + f.elo[k];
    const double *g = &f.gamma[(size_t)k * m];
    double w = 0.0;
    for (int a = 0; a < m; a++) w += g[a] * g[a] * std::abs(st.lamAct[a]);
    st.dinv[k] = (d == complex(0.0, 0.0)) ? complex(0.0, 0.0) : 1.0 / d;
    const double inv = std::abs(st.dinv[k]);
    const double kappa = w * inv;
    if (w > 0.0 && (d == complex(0.0, 0.0) || !(kappa <= ChannelSpace::kDeflateRatio) ||
                    !(kappa * de * inv <= ChannelSpace::kMaxPoleError))) {
      // A deflated level is solved with the channels; the round-off of that
      // solve moves e_k - E by a few ulps of the (unit-scale) bordered system,
      // which matters within the level's width w_k: relative error
      // ~ eps w_k / (|e_k - E|^2 + w_k^2).  Narrower than that, the level
      // matrix is taken.
      const double ad = std::abs(d);
      if (!(ChannelSpace::kEnergyUlps * std::numeric_limits<double>::epsilon() * w / (ad * ad + w * w) <=
            ChannelSpace::kMaxPoleError)) {
        AZUREProfile::Count(AZUREProfile::kLevelRoutePole);
        return false;
      }
      st.isDefl[k] = 1;
      st.defl.push_back(k);
    }
  }
  if ((int)st.defl.size() > ChannelSpace::kMaxDeflated) {
    AZUREProfile::Count(AZUREProfile::kLevelRoutePole);
    return false;
  }
  if (!st.defl.empty()) AZUREProfile::Count(AZUREProfile::kDeflated);
  // R = Gamma^T D^-1 Gamma, symmetric; real unless the formal energies are complex.
  st.R.assign((size_t)m * m, complex(0.0, 0.0));
  if (!f.complexE) {
    std::vector<double> &rr = st.rreal;
    rr.assign((size_t)m * m, 0.0);
    for (int k = 0; k < n; k++) {
      if (st.isDefl[k]) continue;
      const double w = std::real(st.dinv[k]);
      const double *g = &f.gamma[(size_t)k * m];
      for (int a = 0; a < m; a++) {
        const double ga = g[a] * w;
        if (ga == 0.0) continue;
        double *row = &rr[(size_t)a * m];
        for (int b = a; b < m; b++) row[b] += ga * g[b];
      }
    }
    for (int a = 0; a < m; a++)
      for (int b = a; b < m; b++) st.R[(size_t)a * m + b] = st.R[(size_t)b * m + a] = rr[(size_t)a * m + b];
  } else {
    for (int k = 0; k < n; k++) {
      if (st.isDefl[k]) continue;
      const complex w = st.dinv[k];
      const double *g = &f.gamma[(size_t)k * m];
      for (int a = 0; a < m; a++) {
        if (g[a] == 0.0) continue;
        const complex ga = g[a] * w;
        complex *row = &st.R[(size_t)a * m];
        for (int b = a; b < m; b++) row[b] += ga * g[b];
      }
    }
    for (int a = 0; a < m; a++)
      for (int b = a + 1; b < m; b++) st.R[(size_t)b * m + a] = st.R[(size_t)a * m + b];
  }
  return true;
}

/*!
 * Solves the bordered system of SolveChannelSpace for the right-hand side
 * (wp over the deflated levels, wz over the channels), in place: on return wp
 * holds y_p and wz holds z = Lambda Gamma^T y.
 */
void AMatrixFunc::SolveBordered(int j, complex *wp, complex *wz) {
  ChannelState &st = chan_[j - 1];
  const int kd = (int)st.defl.size(), m = (int)(*formal_[j - 1]).act.size(), ns = kd + m;
  if (ns == 0) return;
  vector_c &b = st.w3b;
  b.resize(ns);
  for (int p = 0; p < kd; p++) b[p] = wp[p];
  for (int a = 0; a < m; a++) b[kd + a] = wz[a];
  ChannelSpace::LUSolve(st.lu.data(), st.piv.data(), ns, b.data());
  for (int p = 0; p < kd; p++) wp[p] = b[p];
  for (int a = 0; a < m; a++) wz[a] = b[kd + a];
}

/*!
 * Column b of \f$F = \Gamma^T A_f \Gamma = G^T A G\f$: with the deflated levels p,
 * \f$F_b = \Gamma_p^T y_p + R_b + R z\f$ from the bordered system with
 * right-hand side \f$(\Gamma_{p b}, \Lambda R_b)\f$ (without deflated levels
 * \f$(1 - R\Lambda)^{-1} R_b\f$).
 */
const complex *AMatrixFunc::FColumn(int j, int b) {
  ChannelState &st = chan_[j - 1];
  const ChannelSpace::FormalLevels &f = (*formal_[j - 1]);
  const int m = (int)f.act.size(), kd = (int)st.defl.size();
  complex *col = &st.F[(size_t)b * m];
  if (!st.fcol[b]) {
    vector_c &yp = st.wp, &z = st.wz;
    yp.resize(kd);
    z.resize(m);
    for (int p = 0; p < kd; p++) yp[p] = f.gamma[(size_t)st.defl[p] * m + b];
    for (int a = 0; a < m; a++) z[a] = st.lamAct[a] * st.R[(size_t)a * m + b];
    this->SolveBordered(j, yp.data(), z.data());
    for (int a = 0; a < m; a++) {
      complex s = st.R[(size_t)a * m + b];
      for (int c = 0; c < m; c++) s += st.R[(size_t)a * m + c] * z[c];
      for (int p = 0; p < kd; p++) s += f.gamma[(size_t)st.defl[p] * m + a] * yp[p];
      col[a] = s;
    }
    st.fcol[b] = 1;
  }
  return col;
}

/*!
 * \f$y = A_f w\f$ on the formal route (formal-level coordinates): with
 * \f$t = \Gamma_r^T D_r^{-1} w_r\f$ the bordered system with right-hand side
 * \f$(w_p, \Lambda t)\f$ gives \f$y_p\f$ and \f$z\f$, and
 * \f$y_r = D_r^{-1}(w_r + \Gamma_r z)\f$.
 */
void AMatrixFunc::FormalApply(int j, const complex *w, complex *y) {
  const ChannelSpace::FormalLevels &f = (*formal_[j - 1]);
  ChannelState &st = chan_[j - 1];
  const int n = f.n, m = (int)f.act.size(), kd = (int)st.defl.size();
  vector_c &yp = st.wp, &z = st.wz;
  yp.resize(kd);
  z.assign(m, complex(0.0, 0.0));
  for (int p = 0; p < kd; p++) yp[p] = w[st.defl[p]];
  for (int k = 0; k < n; k++) {
    if (st.isDefl[k]) continue;
    const complex t = st.dinv[k] * w[k];
    if (t == complex(0.0, 0.0)) continue;
    const double *g = &f.gamma[(size_t)k * m];
    for (int a = 0; a < m; a++) z[a] += g[a] * t;
  }
  for (int a = 0; a < m; a++) z[a] *= st.lamAct[a];
  this->SolveBordered(j, yp.data(), z.data());
  for (int k = 0; k < n; k++) {
    if (st.isDefl[k]) continue;
    const double *g = &f.gamma[(size_t)k * m];
    complex s = w[k];
    for (int a = 0; a < m; a++) s += g[a] * z[a];
    y[k] = st.dinv[k] * s;
  }
  for (int p = 0; p < kd; p++) y[st.defl[p]] = yp[p];
}

/*!
 * \f$A u\f$ in channel space: \f$V A_f V^T u\f$ on the formal route
 * (FormalApply), \f$P^{-1}u + P^{-1}G z\f$ with \f$(1 - \Lambda R) z = \Lambda G^T P^{-1} u\f$
 * on the pencil route (P = K - E N).
 */
void AMatrixFunc::ApplyA(int j, const complex *u, complex *out) {
  const ChannelSpace::FormalLevels &f = (*formal_[j - 1]);
  ChannelState &st = chan_[j - 1];
  const int n = f.n, m = (int)f.act.size();
  vector_c &w = st.w1, &q = st.w2;
  if (f.pencil) {
    // P = Qt^T (H - E T) Z^T:  P^-1 u = Z (H - E T)^-1 Qt u,  P^-1 G = Z X,
    // G^T P^-1 u = X^T Z^T u.
    std::vector<double> &re = st.rreal;
    re.assign((size_t)n * 4, 0.0);
    double *qu = re.data(), *zu = re.data() + 2 * (size_t)n;  // n x 2 each: Qt u and Z^T u
    for (int i = 0; i < n; i++) {
      const double *qt = &f.Qt[(size_t)i * n];
      double qr = 0.0, qi = 0.0, zr = 0.0, zi = 0.0;
      for (int k = 0; k < n; k++) {
        qr += qt[k] * std::real(u[k]);
        qi += qt[k] * std::imag(u[k]);
        const double zz = f.Z[(size_t)k * n + i];
        zr += zz * std::real(u[k]);
        zi += zz * std::imag(u[k]);
      }
      qu[2 * i] = qr;
      qu[2 * i + 1] = qi;
      zu[2 * i] = zr;
      zu[2 * i + 1] = zi;
    }
    ChannelSpace::HessenbergSolve(st.P.data(), st.hswap.data(), n, qu, 2);
    q.resize(m);
    for (int a = 0; a < m; a++) {
      complex s(0.0, 0.0);
      for (int k = 0; k < n; k++) s += st.X[(size_t)k * m + a] * complex(zu[2 * k], zu[2 * k + 1]);
      q[a] = st.lamAct[a] * s;
    }
    this->SolveBordered(j, nullptr, q.data());
    vector_c &y = st.w1;
    y.resize(n);
    for (int k = 0; k < n; k++) {
      complex s(qu[2 * k], qu[2 * k + 1]);
      for (int a = 0; a < m; a++) s += st.X[(size_t)k * m + a] * q[a];
      y[k] = s;
    }
    for (int r = 0; r < n; r++) {
      const double *zr = &f.Z[(size_t)r * n];
      complex s(0.0, 0.0);
      for (int k = 0; k < n; k++) s += zr[k] * y[k];
      out[r] = s;
    }
    return;
  }
  w.assign(n, complex(0.0, 0.0));
  if (f.identity)
    for (int k = 0; k < n; k++) w[k] = u[k];
  else
    for (int r = 0; r < n; r++) {
      const complex ur = u[r];
      if (ur == complex(0.0, 0.0)) continue;
      const double *vr = &f.V[(size_t)r * n];
      for (int k = 0; k < n; k++) w[k] += vr[k] * ur;
    }
  vector_c &y = st.w4;
  y.resize(n);
  this->FormalApply(j, w.data(), y.data());
  if (f.identity)
    for (int k = 0; k < n; k++) out[k] = y[k];
  else
    for (int r = 0; r < n; r++) {
      const double *vr = &f.V[(size_t)r * n];
      complex s(0.0, 0.0);
      for (int k = 0; k < n; k++) s += vr[k] * y[k];
      out[r] = s;
    }
}

/*!
 * \f$u^T A v\f$ in channel space (the channel-capture pathways, whose level
 * vectors are masked and so are not columns of \f$G\f$).
 */
complex AMatrixFunc::ChannelBilinear(int j, const complex *u, const complex *v) {
  const ChannelSpace::FormalLevels &f = (*formal_[j - 1]);
  const int n = f.n;
  vector_c av(n);
  this->ApplyA(j, v, av.data());
  complex s(0.0, 0.0);
  for (int k = 0; k < n; k++) s += u[k] * av[k];
  return s;
}

/*!
 * \f$u_c = A_f \Gamma_c\f$ over the formal levels (zero for a channel without
 * widths), cached per energy.
 */
const vector_c &AMatrixFunc::FormalSolve(int j, int ch0) {
  const ChannelSpace::FormalLevels &f = (*formal_[j - 1]);
  ChannelState &st = chan_[j - 1];
  vector_c &u = st.usolve[ch0];
  if (!st.uvalid[ch0]) {
    AZUREProfile::Timer profileTimer(AZUREProfile::kChannelSolve);
    const int n = f.n, m = (int)f.act.size();
    u.assign(n, complex(0.0, 0.0));
    const int b = f.pos[ch0];
    if (b >= 0) {
      vector_c &w = st.w2;
      w.resize(n);
      for (int k = 0; k < n; k++) w[k] = f.gamma[(size_t)k * m + b];
      this->FormalApply(j, w.data(), u.data());
    }
    st.uvalid[ch0] = 1;
  }
  return u;
}

/*!
 * Adjoint of sum_t beta_t F_{ch,chp} in the formal representation
 * \f$M_f = D - \Gamma\Lambda\Gamma^T\f$, \f$F_{ab} = \Gamma_a^T M_f^{-1}\Gamma_b\f$.  With
 * \f$u_c = M_f^{-1}\Gamma_c\f$ and \f$Y_{db} = \delta_{db} + \Lambda_d F_{db}\f$:
 *
 *   \f$\bar e_\lambda \mathrel{+}= -\mathrm{Re}(\bar\beta\, u_{a\lambda} u_{b\lambda})\f$,
 *   \f$\bar\Gamma_{\lambda d} \mathrel{+}= \mathrm{Re}(\bar\beta\,(u_{a\lambda} Y_{db} + u_{b\lambda} Y_{da}))\f$,
 *
 * the first from \f$dD\f$, the second from the explicit widths and from
 * \f$-\Gamma\Lambda\Gamma^T\f$ (using \f$\Gamma^T u_b = F e_b\f$).  O(n c) per term.
 */
void AMatrixFunc::FormalBilinearAdjoint(int j, const std::vector<BilinearTerm> &terms, GradAccum &accum) {
  const ChannelSpace::FormalLevels &f = (*formal_[j - 1]);
  ChannelState &st = chan_[j - 1];
  const int n = f.n, nch = f.nch;
  const bool toFormal = !f.identity;
  if (toFormal) accum.formalUsed = true;
  std::vector<double> &eOut = toFormal ? accum.fe[j - 1] : accum.e[j - 1];
  std::vector<std::vector<double>> &gOut = toFormal ? accum.fgamma[j - 1] : accum.gamma[j - 1];
  const std::vector<int> &lev = active_levels_[j - 1];
  vector_c ya(nch), yb(nch);
  auto column = [&](int c0, vector_c &y) {
    std::fill(y.begin(), y.end(), complex(0.0, 0.0));
    y[c0] = 1.0;
    const int b = f.pos[c0];
    if (b < 0) return;
    const complex *Fb = this->FColumn(j, b);
    const int m = (int)f.act.size();
    for (int a = 0; a < m; a++) y[f.act[a]] += st.lamAct[a] * Fb[a];
  };
  for (const BilinearTerm &t : terms) {
    const int a = t.ch - 1, b = t.chp - 1;
    const complex cb = conj(t.beta);
    const vector_c &ua = this->FormalSolve(j, a);
    const vector_c &ub = this->FormalSolve(j, b);
    column(a, ya);
    column(b, yb);
    for (int k = 0; k < n; k++) {
      const int row = toFormal ? k : lev[k] - 1;
      const complex pa = cb * ua[k], pb = cb * ub[k];
      eOut[row] -= std::real(pa * ub[k]);
      if (pa == complex(0.0, 0.0) && pb == complex(0.0, 0.0)) continue;
      double *g = &gOut[row][0];
      for (int d = 0; d < nch; d++) g[d] += std::real(pa * yb[d] + pb * ya[d]);
    }
  }
}

/*!
 * Formal cotangents (ebar, Gammabar) -> level parameters.  With K V = N V diag(e),
 * V^T N V = 1 and Gamma = V^T G, a variation gives
 *   de_l = v_l^T (dK - e_l dN) v_l,  dV = V Omega,
 *   Omega_kl = v_k^T (dK - e_l dN) v_l / (e_l - e_k) (k != l),  Omega_ll = -v_l^T dN v_l / 2,
 * so with H = V^T Vbar = Gamma Gammabar^T (Vbar = G Gammabar^T)
 *   Kbar = V Phi V^T,  Phi_ll = ebar_l,  Phi_kl = H_kl / (e_l - e_k),
 *   Nbar = V Psi V^T,  Psi_ll = -e_l ebar_l - H_ll / 2,  Psi_kl = -e_l H_kl / (e_l - e_k),
 * taken symmetric (K and N are), and Gbar = V Gammabar.  Kbar and Nbar go
 * through the structural level-matrix adjoint: with the channel kernel set to
 * zero, LevelMatrixAdjoint(X, E) contracts X with dK - E dN, so
 *   <Kbar, dK> + <Nbar, dN> = LMA(Kbar + Nbar, 0) + LMA(-Nbar, 1).
 */
void AMatrixFunc::ResolveFormalCotangents(GradAccum &accum, const vector_matrix_r *shiftDeriv) {
  if (!accum.formalUsed) return;
  this->ClearMatrices();
  // Deferred level-matrix cotangents: with X_K = sum Re Mbar, X_N = sum E Re Mbar,
  // <X_K, dK> - <X_N, dN> = LMA(X_K - X_N, 0) + LMA(X_N, 1) with a zero kernel.
  for (int j = 1; j <= compound()->NumJGroups() && j <= (int)accum.lmUsed.size(); j++) {
    if (!accum.lmUsed[j - 1]) continue;
    JGroup *jg = compound()->GetJGroup(j);
    const int n = (int)std::lround(std::sqrt((double)accum.lmK[j - 1].size()));
    matrix_c X1(n, vector_c(n)), X2(n, vector_c(n));
    for (int r = 0; r < n; r++)
      for (int c = 0; c < n; c++) {
        const double xk = accum.lmK[j - 1][(size_t)r * n + c], xn = accum.lmN[j - 1][(size_t)r * n + c];
        X1[r][c] = xk - xn;
        X2[r][c] = xn;
      }
    const std::vector<complex> zero(jg->NumChannels(), complex(0.0, 0.0));
    this->LevelMatrixAdjoint(nullptr, j, X1, 0.0, accum, shiftDeriv, zero.data());
    this->LevelMatrixAdjoint(nullptr, j, X2, 1.0, accum, shiftDeriv, zero.data());
    std::fill(accum.lmK[j - 1].begin(), accum.lmK[j - 1].end(), 0.0);
    std::fill(accum.lmN[j - 1].begin(), accum.lmN[j - 1].end(), 0.0);
    accum.lmUsed[j - 1] = 0;
  }
  for (int j = 1; j <= compound()->NumJGroups(); j++) {
    JGroup *jg = compound()->GetJGroup(j);
    if (!jg->IsInRMatrix() || j > (int)accum.fe.size()) continue;
    bool any = false;
    for (size_t k = 0; k < accum.fe[j - 1].size() && !any; k++) {
      if (accum.fe[j - 1][k] != 0.0) any = true;
      for (double g : accum.fgamma[j - 1][k])
        if (g != 0.0) any = true;
    }
    if (!any) continue;
    formal_[j - 1] = ChannelSpace::Cached(jg);
    formal_jgroup_[j - 1] = jg;
    ChannelSpace::FormalLevels &f = (*formal_[j - 1]);
    ChannelSpace::Update(jg, active_levels_[j - 1], configure(), f);
    const int n = f.n, nch = f.nch, m = (int)f.act.size();
    const std::vector<double> &ebar = accum.fe[j - 1];
    const std::vector<std::vector<double>> &gbar = accum.fgamma[j - 1];
    std::vector<double> e(n);
    for (int k = 0; k < n; k++) e[k] = std::real(f.e[k]);
    // H = Gamma Gammabar^T over the coupled channels.
    std::vector<double> H((size_t)n * n, 0.0);
    for (int k = 0; k < n; k++)
      for (int l = 0; l < n; l++) {
        double s = 0.0;
        for (int a = 0; a < m; a++) s += f.gamma[(size_t)k * m + a] * gbar[l][f.act[a]];
        H[(size_t)k * n + l] = s;
      }
    std::vector<double> Phi((size_t)n * n, 0.0), Psi((size_t)n * n, 0.0);
    for (int l = 0; l < n; l++) {
      Phi[(size_t)l * n + l] = ebar[l];
      Psi[(size_t)l * n + l] = -e[l] * ebar[l] - 0.5 * H[(size_t)l * n + l];
      for (int k = 0; k < n; k++) {
        if (k == l) continue;
        const double num = H[(size_t)k * n + l] - H[(size_t)l * n + k];
        const double den = e[l] - e[k];
        const double ratio = (num == 0.0 || den == 0.0) ? 0.0 : num / (2.0 * den);
        Phi[(size_t)k * n + l] = ratio;
        Psi[(size_t)k * n + l] = -e[l] * ratio - 0.5 * H[(size_t)l * n + k];
      }
    }
    // Psi as built is symmetric only up to round-off; use its symmetric part.
    for (int k = 0; k < n; k++)
      for (int l = k + 1; l < n; l++) {
        const double s = 0.5 * (Psi[(size_t)k * n + l] + Psi[(size_t)l * n + k]);
        Psi[(size_t)k * n + l] = Psi[(size_t)l * n + k] = s;
      }
    // V X V^T for X = Phi, Psi.
    auto congruence = [&](const std::vector<double> &X) {
      std::vector<double> T((size_t)n * n, 0.0), out((size_t)n * n, 0.0);
      for (int r = 0; r < n; r++)
        for (int l = 0; l < n; l++) {
          double s = 0.0;
          for (int k = 0; k < n; k++) s += f.V[(size_t)r * n + k] * X[(size_t)k * n + l];
          T[(size_t)r * n + l] = s;
        }
      for (int r = 0; r < n; r++)
        for (int c = 0; c < n; c++) {
          double s = 0.0;
          for (int l = 0; l < n; l++) s += T[(size_t)r * n + l] * f.V[(size_t)c * n + l];
          out[(size_t)r * n + c] = s;
        }
      return out;
    };
    const std::vector<double> Kbar = congruence(Phi), Nbar = congruence(Psi);
    matrix_c X1(n, vector_c(n)), X2(n, vector_c(n));
    for (int r = 0; r < n; r++)
      for (int c = 0; c < n; c++) {
        X1[r][c] = Kbar[(size_t)r * n + c] + Nbar[(size_t)r * n + c];
        X2[r][c] = -Nbar[(size_t)r * n + c];
      }
    // The kernel: zero, but for the reference shifts that Update moved into K
    // (constants, like the kernel: F does not depend on them).
    std::vector<complex> kernel(nch, complex(0.0, 0.0));
    for (int c = 0; c < nch; c++)
      if (jg->GetChannel(c + 1)->GetRadType() == 'P')
        kernel[c] = jg->GetChannel(c + 1)->GetBoundaryCondition() - f.lamShift[c];
    this->LevelMatrixAdjoint(nullptr, j, X1, 0.0, accum, shiftDeriv, kernel.data());
    this->LevelMatrixAdjoint(nullptr, j, X2, 1.0, accum, shiftDeriv, kernel.data());
    // Gbar = V Gammabar.
    const std::vector<int> &lev = active_levels_[j - 1];
    for (int r = 0; r < n; r++)
      for (int c = 0; c < nch; c++) {
        double s = 0.0;
        for (int k = 0; k < n; k++) s += f.V[(size_t)r * n + k] * gbar[k][c];
        accum.gamma[j - 1][lev[r] - 1][c] += s;
      }
    for (size_t k = 0; k < accum.fe[j - 1].size(); k++) {
      accum.fe[j - 1][k] = 0.0;
      std::fill(accum.fgamma[j - 1][k].begin(), accum.fgamma[j - 1][k].end(), 0.0);
    }
  }
  accum.formalUsed = false;
}

/*!
 * AZURE_CHANNEL_VERIFY=1: compares the channel-space bilinear forms
 * \f$G^T A G\f$ of J-group j with those of the level matrix, filled here and
 * LU-factored in double and in extended precision, and reports each new worst
 * difference (relative to the largest form of the J-group) on std::cerr.
 */
void AMatrixFunc::VerifyChannelSpace(int j, EPoint *point, double inEnergy) {
  typedef std::complex<long double> cld;
  static double worst = 0.0;
  this->FillLevelMatrix(j, point, inEnergy);
  LevelMatrixSolver lu;
  if (!lu.Decompose(a_inv_matrices_[j - 1])) return;
  const ChannelSpace::FormalLevels &f = (*formal_[j - 1]);
  const int n = f.n, nch = f.nch, m = (int)f.act.size();
  // Extended-precision LU of the same level matrix, as the referee.
  const matrix_c &M = a_inv_matrices_[j - 1];
  std::vector<cld> a((size_t)n * n);
  for (int r = 0; r < n; r++)
    for (int c = 0; c < n; c++) a[(size_t)r * n + c] = cld(std::real(M[r][c]), std::imag(M[r][c]));
  std::vector<int> pv(n);
  for (int k = 0; k < n; k++) {
    int p = k;
    for (int i = k + 1; i < n; i++)
      if (std::abs(a[(size_t)i * n + k]) > std::abs(a[(size_t)p * n + k])) p = i;
    pv[k] = p;
    if (p != k)
      for (int c = 0; c < n; c++) std::swap(a[(size_t)k * n + c], a[(size_t)p * n + c]);
    for (int i = k + 1; i < n; i++) {
      const cld fct = a[(size_t)i * n + k] / a[(size_t)k * n + k];
      a[(size_t)i * n + k] = fct;
      for (int c = k + 1; c < n; c++) a[(size_t)i * n + c] -= fct * a[(size_t)k * n + c];
    }
  }
  std::vector<complex> ref((size_t)m * m), lev((size_t)m * m);
  std::vector<cld> x(n);
  vector_c ga(n), gb(n);
  for (int b = 0; b < m; b++) {
    for (int k = 0; k < n; k++) x[k] = f.gcut[(size_t)k * nch + f.act[b]];
    for (int k = 0; k < n; k++)
      if (pv[k] != k) std::swap(x[k], x[pv[k]]);
    for (int i = 1; i < n; i++)
      for (int c = 0; c < i; c++) x[i] -= a[(size_t)i * n + c] * x[c];
    for (int i = n - 1; i >= 0; i--) {
      for (int c = i + 1; c < n; c++) x[i] -= a[(size_t)i * n + c] * x[c];
      x[i] /= a[(size_t)i * n + i];
    }
    for (int k = 0; k < n; k++) gb[k] = f.gcut[(size_t)k * nch + f.act[b]];
    for (int aa = 0; aa < m; aa++) {
      cld sum = 0;
      for (int k = 0; k < n; k++) sum += (long double)f.gcut[(size_t)k * nch + f.act[aa]] * x[k];
      ref[(size_t)aa * m + b] = complex((double)sum.real(), (double)sum.imag());
      for (int k = 0; k < n; k++) ga[k] = f.gcut[(size_t)k * nch + f.act[aa]];
      lev[(size_t)aa * m + b] = lu.Bilinear(ga, gb);
    }
  }
  double fmax = 0.0, dChannel = 0.0, dLevel = 0.0, dRoutes = 0.0;
  for (int b = 0; b < m; b++)
    for (int aa = 0; aa < m; aa++) {
      const complex ch = this->FColumn(j, b)[aa];
      const size_t k = (size_t)aa * m + b;
      fmax = std::max(fmax, std::abs(ref[k]));
      dChannel = std::max(dChannel, std::abs(ch - ref[k]));
      dLevel = std::max(dLevel, std::abs(lev[k] - ref[k]));
      dRoutes = std::max(dRoutes, std::abs(ch - lev[k]));
    }
  if (!(fmax > 0.0)) return;
  const double rel = dRoutes / fmax;
#pragma omp critical(channel_verify)
  if (rel > worst) {
    worst = rel;
    std::cerr << "[channel-verify] J-group " << j << " (" << n << " levels, " << m << " channels, "
              << (f.pencil ? "pencil" : f.identity ? "standard" : "formal") << ") at " << inEnergy
              << " MeV: channel vs level route " << rel << " of the largest form; against an extended-precision LU: channel "
              << dChannel / fmax << ", level " << dLevel / fmax << std::endl;
  }
}

/*!
 * The cotangent coefficients of the differential sum of `theDecay` for the
 * T-matrix held now, by pathway and Legendre order: model += cR Re(sum_i
 * Z1Z2_i T1_i conj(T2_i) P_{L_i}) gives T1bar += cR Z1Z2 P_L T2 and
 * T2bar += cR Z1Z2 P_L T1, so each pathway's cotangent is cR sum_L q_L P_L with
 * q_L summed over the interferences once.  `particleExit` keeps only the
 * k-groups of equal channel spins, as the forward sum.  False (nothing built)
 * for an interference type the adjoint does not know.
 */
bool AMatrixFunc::TBarProducts(Decay *theDecay, bool particleExit) {
  if (tbar_decay_ == theDecay && tbar_particle_ == particleExit) return true;
  if (!theDecay->TablesBuilt()) return false;
  const std::vector<FlatInterference> &flat = theDecay->FlatInterferences();
  const int nl = theDecay->NumInterferenceOrders();
  const size_t np = theDecay->InterferencePathways().size();
  tbar_q_.assign(np * nl, complex(0.0, 0.0));
  const std::vector<vector_c> *T = tmatrix_.empty() ? nullptr : &tmatrix_[0];
  for (const FlatInterference &f : flat) {
    if (particleExit && !f.spEqual) continue;
    if (!f.r1 || !f.r2) {  // an interference type the adjoint does not know
      tbar_decay_ = nullptr;
      return false;
    }
    const complex T1 = (f.r1 == 'E') ? ec_tmatrix_[f.k - 1][f.m1 - 1] : (T ? (*T)[f.k - 1][f.m1 - 1] : complex(0.0, 0.0));
    const complex T2 = (f.r2 == 'E') ? ec_tmatrix_[f.k - 1][f.m2 - 1] : (T ? (*T)[f.k - 1][f.m2 - 1] : complex(0.0, 0.0));
    tbar_q_[(size_t)f.s1 * nl + f.l] += f.z * T2;
    tbar_q_[(size_t)f.s2 * nl + f.l] += f.z * T1;
  }
  tbar_decay_ = theDecay;
  tbar_particle_ = particleExit;
  return true;
}

/*!
 * See the declaration.  With conj(Mbar) = -sum_t conj(beta_t) u_ch u_chp^T
 * (u = A g, the widths as they are), a level parameter's derivative is
 *   -Re(sum_t conj(beta_t) u_ch^T dM u_chp),  dM = dK - E dN - (dG Lambda G^T + G Lambda dG^T),
 * and for dK = e_l k^T + k e_l^T:  u^T dK v = u_l (k . v) + (k . u) v_l.  The
 * widths also enter g_ch^T A g_chp explicitly: conj(beta) u_chp and
 * conj(beta) u_ch for the entrance and exit width.  `add(level, ch, value)`
 * receives the derivatives (ch = -1: the level energy).
 */
template <class Add>
bool AMatrixFunc::LevelAdjoint(EPoint *energyPoint, int j, const std::vector<BilinearTerm> &terms,
                               const FreeMask &mask, const vector_matrix_r *shiftDeriv, Add add) {
  if (j < 1 || j > (int)route_.size() || j > (int)mask.e.size()) return false;
  const bool channel = (route_[j - 1] == kRouteChannel);
  if (!channel && !solver_valid_[j - 1]) return false;
  // The standard level matrix in channel space has the cheaper formal adjoint.
  if (channel && (*formal_[j - 1]).identity) return false;
  JGroup *jg = compound()->GetJGroup(j);
  const int n = (int)active_levels_[j - 1].size(), nch = jg->NumChannels();
  if (n == 0 || n != jg->NumLevels()) return false;
  const bool park = (configure().paramMask & Config::USE_PARK_FORMALISM);
  const bool brune = !park && (configure().paramMask & Config::USE_BRUNE_FORMALISM);
  // Brune takes S' from shiftDeriv, Park from the levels (CNuc::CalcShiftFunctions).
  if (brune && !shiftDeriv) return false;
  // The channels with the shift terms, as LevelMatrixAdjoint takes them.
  std::vector<char> &particle = la_particle_;
  particle.assign(nch, 0);
  for (int c = 0; c < nch; c++) {
    AChannel *ch = jg->GetChannel(c + 1);
    particle[c] = park ? (ch->GetRadType() == 'P' && compound()->GetPair(ch->GetPairNum())->GetPType() == 0)
                       : (brune && ch->GetRadType() == 'P');
  }
  const ChannelSpace::LevelDerivatives &ld = ChannelSpace::CurrentLevelDerivatives(
      jg, park ? 2 : brune ? 1 : 0, particle, mask.e[j - 1], mask.g[j - 1],
      (brune && shiftDeriv) ? &(*shiftDeriv)[j - 1] : nullptr);
  const ChannelSpace::LevelArrays &arrays = ChannelSpace::CurrentLevelArrays(jg);
  const std::vector<complex> &lam = this->KernelLambda(energyPoint, j);
  const double E = this->InEnergy(energyPoint);
  // The channels the terms involve, their solves u_c = A g_c and g_d . u_c.
  std::vector<int> &used = la_used_;
  used.clear();
  std::vector<int> &slot = la_slot_;
  slot.assign(nch, -1);
  for (const BilinearTerm &t : terms)
    for (int c : {t.ch - 1, t.chp - 1})
      if (slot[c] < 0) {
        slot[c] = (int)used.size();
        used.push_back(c);
      }
  const int nu = (int)used.size();
  std::vector<const vector_c *> &u = la_u_;
  u.resize(nu);
  for (int s = 0; s < nu; s++) u[s] = &this->GetRawChannelSolve(j, used[s] + 1);
  vector_c &gu = la_gu_;  // [s][d] = g_d . u_s
  gu.assign((size_t)nu * nch, complex(0.0, 0.0));
  for (int i = 0; i < n; i++) {
    const double *g = &arrays.G[(size_t)i * nch];
    for (int s = 0; s < nu; s++) {
      const complex ui = (*u[s])[i];
      complex *row = &gu[(size_t)s * nch];
      for (int d = 0; d < nch; d++) row[d] += g[d] * ui;
    }
  }
  vector_c &ku = la_ku_;
  ku.resize(nu);
  for (const ChannelSpace::LevelDerivatives::Entry &e : ld.entries) {
    const int l = e.level;
    if (!e.zero) {
      std::fill(ku.begin(), ku.end(), complex(0.0, 0.0));
      for (int i = 0; i < n; i++) {
        const double kt = e.k[i] - E * e.q[i];
        if (kt == 0.0) continue;
        for (int s = 0; s < nu; s++) ku[s] += kt * (*u[s])[i];
      }
    }
    double v = 0.0;
    for (const BilinearTerm &t : terms) {
      const int sa = slot[t.ch - 1], sb = slot[t.chp - 1];
      const complex ua = (*u[sa])[l], ub = (*u[sb])[l];
      complex dm(0.0, 0.0);
      if (!e.zero) dm += ua * ku[sb] + ku[sa] * ub;
      if (e.ch >= 0) dm -= lam[e.ch] * (ua * gu[(size_t)sb * nch + e.ch] + gu[(size_t)sa * nch + e.ch] * ub);
      const complex cb = conj(t.beta);
      v -= std::real(cb * dm);
      if (e.ch == t.ch - 1) v += std::real(cb * ub);
      if (e.ch == t.chp - 1) v += std::real(cb * ua);
    }
    add(l, e.ch, v);
  }
  return true;
}

bool AMatrixFunc::LevelBasis(EPoint *energyPoint, int j, const BilinearTerm &term, const FreeMask &mask,
                             const vector_matrix_r *shiftDeriv, std::vector<double> &block) {
  const int n = (int)active_levels_[j - 1].size(), nch = compound()->GetJGroup(j)->NumChannels();
  block.assign((size_t)n * (1 + nch), 0.0);
  return this->LevelAdjoint(energyPoint, j, std::vector<BilinearTerm>{term}, mask, shiftDeriv,
                            [&](int l, int ch, double v) {
                              if (ch < 0)
                                block[l] += v;
                              else
                                block[n + (size_t)l * nch + ch] += v;
                            });
}
