// AMatrixFunc, the forward T matrix: the level matrices of a point's energy (FillMatrices, the level route and the
// channel-space route), their solves and bilinear forms, and the T matrix (CalculateTMatrix).  Reverse mode: the
// T cotangents of an observable in AMatrixFuncTBar.cpp, T -> level parameters in AMatrixFuncAdjoint.cpp.
#include "AMatrixFunc.h"
#include "RunDiagnostics.h"
#include <mutex>
#include "azure2/support/Environment.h"
#include "Profile.h"
#include "AngCoeff.h"
#include "PolarizationFunc.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "Unit.h"
#include "AZUREGrad.h"
#include "azure2/physics/ExternalCapture.h"
#include <stdexcept>
#include "JGroup.h"
#include "ALevel.h"
#include "AChannel.h"
#include "PPair.h"
#include <assert.h>
#include <iomanip>
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
 * The AMatrixFunc object is created with reference to a CompoundSystem object.
 */

AMatrixFunc::AMatrixFunc(const CompoundSystem *compound, const LevelState &state, const Config &configure) :
  configure_(&configure),
  compound_(compound),
  state_(&state) {}

/*!
 * Returns the A-Matrix of a J-Group, materializing it from the stored
 * factorization on first use.  Returns an empty matrix if the J-Group has no
 * usable factorization.
 */

const matrix_c &AMatrixFunc::GetAMatrix(int jGroupNum) {
  static const matrix_c emptyMatrix;
  if (jGroupNum < 1 || jGroupNum > (int)solvers_.size()) return emptyMatrix;
  if (route_[jGroupNum - 1] == kRouteChannel) {
    // A = V A_f V^T, one column per unit vector, kept with the object's per-energy state, like the factorization's
    // own inverse.
    ChannelState &st = chan_[jGroupNum - 1];
    if (!st.inverseValid) {
      const int n = (*formal_[jGroupNum - 1]).n;
      st.inverse.assign(n, vector_c(n, complex(0.0, 0.0)));
      vector_c unit(n), col(n);
      for (int c = 0; c < n; c++) {
        std::fill(unit.begin(), unit.end(), complex(0.0, 0.0));
        unit[c] = 1.0;
        this->ApplyA(jGroupNum, unit.data(), col.data());
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

complex AMatrixFunc::GetAMatrixElement(int jGroupNum, int lambdaNum, int muNum) {
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

    const JGroup *jg = compound()->GetJGroup(j);
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
 * This function creates the inverted A-Matrix from the parameters in the CompoundSystem object.
 */

void AMatrixFunc::FillMatrices(Site *point) {
  AZUREProfile::Timer profileTimer(configure().profile, AZUREProfile::kFill);
  const double inEnergy = this->InEnergy(point);
  lastInEnergy_ = inEnergy;
  // Groups of another compound nucleus are not part of this point's matrices
  // (Unit::CalcEDependentValues left their rows empty).
  activeNucleus_ = compound()->GetPair(compound()->GetPairNumFromKey(point->GetEntranceKey()))->GetNucleus();
  const bool channelSpace = !level_only_ && ChannelSpace::Enabled(configure());
  for (int j = 1; j <= compound()->NumJGroups(); j++) {
    const JGroup *jGroup = compound()->GetJGroup(j);
    if (!jGroup->IsInRMatrix() || jGroup->GetNucleus() != activeNucleus_) continue;
    route_[j - 1] = kRouteLevel;
    if (channelSpace) {
      if (!formal_[j - 1] || formal_jgroup_[j - 1] != jGroup) {
        formal_[j - 1] = ChannelSpace::Cached(ChannelMemo(), jGroup);
        formal_jgroup_[j - 1] = jGroup;
      }
      ChannelSpace::FormalLevels &f = (*formal_[j - 1]);
      ChannelSpace::Update(jGroup, state(), j, active_levels_[j - 1], configure(), f);
      if (f.usable) {
        if (this->SolveChannelSpace(j, point, inEnergy)) {
          route_[j - 1] = kRouteChannel;
          AZUREProfile::Count(configure().profile, AZUREProfile::kChannelRoute);
          static const bool verify = azure2::support::ProcessEnvironment().Set(azure2::support::Switch::kChannelVerify);
          if (verify) this->VerifyChannelSpace(j, point, inEnergy);
          continue;
        }
      } else {
        AZUREProfile::Count(configure().profile, f.status == 1 ? AZUREProfile::kLevelRouteSize : AZUREProfile::kLevelRouteSetup);
      }
    } else {
      AZUREProfile::Count(configure().profile, AZUREProfile::kLevelRouteForced);
    }
    this->FillLevelMatrix(j, point, inEnergy);
  }
}

/*!
 * The level matrix of J-group j at one energy (the level route): the level parameters of the group's active levels
 * and the channel functions of the point are collected into arrays, and azure2::physics::fill_level_matrix
 * (azure2/physics/CollisionMatrix.h) has the arithmetic.
 */

void AMatrixFunc::FillLevelMatrix(int j, Site *point, double inEnergy) {
  const bool park = (configure().paramMask & Config::USE_PARK_FORMALISM);
  const bool brune = (configure().paramMask & Config::USE_BRUNE_FORMALISM);
  const JGroup *jGroup = compound()->GetJGroup(j);
  const std::vector<int> &levels = active_levels_[j - 1];
  const int numLevels = (int)levels.size();
  const int numChannels = jGroup->NumChannels();
  const LevelState &state = this->state();
  LevelMatrixGather &g = level_gather_;
  const size_t cells = (size_t)numLevels * numChannels;
  g.E.resize(numLevels);
  g.gamma.resize(cells);
  g.shift.resize(cells);
  g.lo.resize(numChannels);
  g.radType.resize(numChannels);
  g.boundary.assign(numChannels, 0.0);
  if (park || brune) {
    g.shiftD1.resize(cells);
    g.shiftD2.resize(cells);
  }
  if (park) g.parkNorm.resize(numLevels);
  for (int i = 0; i < numLevels; i++) {
    const int la = levels[i];
    g.E[i] = state.E(j, la);
    const size_t row = (size_t)i * numChannels;
    for (int ch = 1; ch <= numChannels; ch++) {
      g.gamma[row + ch - 1] = state.Gamma(j, la, ch);
      g.shift[row + ch - 1] = state.Shift(j, la, ch);
    }
    if (park || brune)
      for (int ch = 1; ch <= numChannels; ch++) {
        g.shiftD1[row + ch - 1] = state.ShiftDerivative(j, la, ch);
        g.shiftD2[row + ch - 1] = state.ShiftSecondDerivative(j, la, ch);
      }
    if (park) g.parkNorm[i] = state.ParkNorm(j, la);
  }
  for (int ch = 1; ch <= numChannels; ch++) {
    const AChannel *channel = jGroup->GetChannel(ch);
    g.lo[ch - 1] = point->GetLoElement(j, ch);
    g.radType[ch - 1] = channel->GetRadType();
    // B_c is read where the level matrix uses it: a particle channel of Brune or Park that some level couples to
    // (AChannel::GetBoundaryCondition refuses before B_c is set).
    if (!(park || brune) || g.radType[ch - 1] != 'P') continue;
    for (int i = 0; i < numLevels; i++)
      if (fabs(g.gamma[(size_t)i * numChannels + ch - 1]) >= azure2::physics::kWidthCutoff) {
        g.boundary[ch - 1] = channel->GetBoundaryCondition();
        break;
      }
  }
  azure2::physics::LevelMatrixTables tables;
  tables.n = numLevels;
  tables.nc = numChannels;
  tables.E = g.E.data();
  tables.gamma = g.gamma.data();
  tables.shift = g.shift.data();
  tables.shift_d1 = g.shiftD1.data();
  tables.shift_d2 = g.shiftD2.data();
  tables.park_norm = g.parkNorm.data();
  tables.lo = g.lo.data();
  tables.boundary = g.boundary.data();
  tables.rad_type = g.radType.data();
  tables.brune = brune;
  tables.park = park;
  tables.rmc = (configure().paramMask & Config::USE_RMC_FORMALISM);
  azure2::physics::fill_level_matrix(tables, inEnergy, a_inv_matrices_[j - 1]);
}

/*!
 * This function factors the inverse A-Matrix of each J-Group.  The A-Matrix
 * itself is not formed: CalculateTMatrix needs only bilinear forms of it, which
 * the factorization supplies with one triangular solve per channel, and the
 * callers that do need individual elements get them from GetAMatrixElement,
 * which materializes the inverse on demand.
 */

void AMatrixFunc::InvertMatrices() {
  AZUREProfile::Timer profileTimer(configure().profile, AZUREProfile::kFactor);
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
      if (!solver_valid_[j - 1]) {
        // The group then contributes nothing at this energy.  A level without
        // widths no longer gets here (FillLevelMatrix); say so once if
        // anything else does.
        if (configure().FirstTime("singular level matrix")) {
          const JGroup *jg = compound()->GetJGroup(j);
          std::cerr << "Warning: singular level matrix of J-group " << jg->GetJ() << (jg->GetPi() < 0 ? "-" : "+")
                    << " at E_x = " << std::setprecision(10) << lastInEnergy_
                    << " MeV; the group is left out at that energy (reported once)." << std::endl;
        }
      }
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
    const JGroup *jGroup = compound()->GetJGroup(jGroupNum);
    const std::vector<int> &map = level_active_index_[jGroupNum - 1];
    gammas.assign(active_levels_[jGroupNum - 1].size(), 0.0);
    for (int la = 1; la < (int)map.size(); la++) {
      const int row = map[la];
      if (row == 0) continue;
      double gamma = state().Gamma(jGroupNum, la, chNum);
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
    AZUREProfile::Timer profileTimer(configure().profile, AZUREProfile::kChannelSolve);
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
  return azure2::physics::width_bilinear(gammas.data(), solve.data(), solve.size());
}

/*!
 * Channel capture with --ignore-externals excludes the levels whose input width in the internal (photon) channel is
 * zero from both indices (external_capture.md section 3), which is a mask on the level set rather than the
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
    const JGroup *jGroup = compound()->GetJGroup(jGroupNum);
    const std::vector<int> &map = level_active_index_[jGroupNum - 1];
    const int n = (int)active_levels_[jGroupNum - 1].size();
    chan_cap_entrance_.assign(n, complex(0.0, 0.0));
    chan_cap_exit_.assign(n, complex(0.0, 0.0));
    // The --ignore-externals mask (external_capture.md section 3): in the formal basis the projector T^-1 P T of
    // the Brune levels whose input width in the photon channel is zero; otherwise those levels are left out.
    const int maskSlot = maskChannel ? compound()->FormalCaptureMaskSlot(jGroupNum, maskChannel) : -1;
    const std::vector<double> *projector = maskSlot >= 0 ? &state().CaptureProjector(maskSlot) : nullptr;
    if (projector && !projector->empty()) {
      const std::vector<int> &levels = compound()->FormalCaptureMaskOf(maskSlot).levels;
      const std::size_t nl = levels.size();
      for (std::size_t i = 0; i < nl; i++) {
        const int row = map[levels[i]];
        if (row == 0) continue;
        double en = 0.0, ex = 0.0;
        for (std::size_t j = 0; j < nl; j++) {
          en += (*projector)[i * nl + j] * state().Gamma(jGroupNum, levels[j], chNum);
          ex += (*projector)[i * nl + j] * state().Gamma(jGroupNum, levels[j], chpNum);
        }
        chan_cap_entrance_[row - 1] = en;
        chan_cap_exit_[row - 1] = ex;
      }
    } else
      for (int la = 1; la < (int)map.size(); la++) {
        const int row = map[la];
        if (row == 0) continue;
        if (maskChannel && state().ZeroInputWidth(jGroupNum, la, maskChannel)) continue;
        chan_cap_entrance_[row - 1] = state().Gamma(jGroupNum, la, chNum);
        chan_cap_exit_[row - 1] = state().Gamma(jGroupNum, la, chpNum);
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

void AMatrixFunc::CalculateTMatrix(Site *point) {
  AZUREProfile::Timer profileTimer(configure().profile, AZUREProfile::kTMatrix);
  tbar_decay_ = nullptr;
  int entranceKey = point->GetEntranceKey();
  int exitKey = point->GetExitKey();
  int aa = compound()->GetPairNumFromKey(entranceKey);
  int exitPairNum = compound()->GetPairNumFromKey(exitKey);
  const azure2::EntrancePathways &P = compound()->Pathways().Entrance(aa);

  int irEnd;
  int irStart;
  bool isRMC = false;
  if ((configure().paramMask & Config::USE_RMC_FORMALISM) &&
      compound()->GetPair(exitPairNum)->GetPType() == 10) {
    irStart = 1;
    irEnd = (int)P.exits.size();
    isRMC = true;
  } else {
    irStart = P.ExitOf(exitPairNum);
    if (!irStart) throw std::runtime_error(NoPathwayMessage(aa, exitPairNum));
    irEnd = irStart;
  }
  for (int ir = irStart; ir <= irEnd; ir++) {
    const azure2::ExitRow &theDecay = P.exits[ir - 1];
    for (int k = 1; k <= theDecay.num_groups; k++) {
      const azure2::SpinGroupRow &theKGroup = P.Group(theDecay, k);
      for (int m = 1; m <= theKGroup.num_internal; m++) {
        const azure2::InternalRow &theMGroup = P.Internal(theKGroup, m);
        const int jNum = theMGroup.group;
        const int chNum = theMGroup.ch;
        const int chpNum = theMGroup.chp;

        const azure2::physics::ChannelPhases entrance{point->GetExpCoulombPhase(jNum, chNum),
                                                      point->GetExpHardSpherePhase(jNum, chNum)};
        const azure2::physics::ChannelPhases exit{point->GetExpCoulombPhase(jNum, chpNum),
                                                  point->GetExpHardSpherePhase(jNum, chpNum)};
        const complex umatrix = azure2::physics::pathway_amplitude(point->GetSqrtPenetrability(jNum, chNum),
                                                                   point->GetSqrtPenetrability(jNum, chpNum),
                                                                   this->GetUBilinear(jNum, chNum, chpNum));
        if (isRMC)
          this->AddTMatrixElement(k, m, azure2::physics::rmc_t_matrix_element(umatrix), ir);
        else
          this->AddTMatrixElement(k, m,
                                  azure2::physics::t_matrix_element(chNum == chpNum, entrance, exit, umatrix));
      }
      for (int m = 1; m <= theKGroup.num_external; m++) {
        const azure2::ExternalRow &theECMGroup = P.External(theKGroup, m);
        const ALevel *finalLevel = compound()->GetJGroup(theECMGroup.final_group)->GetLevel(theECMGroup.final_level);
        const LevelState &state = this->state();
        // Use energy-shift aware EC amplitude calculation if energy shifts are active
        complex ecAmplitude;
        if (configure().paramMask & Config::USE_EXTERNAL_CAPTURE) {
          // Use the new method that accounts for energy shifts through interpolation
          ecAmplitude = point->GetECAmplitudeWithShift(k, m, compound(), configure());
        } else {
          // Use cached amplitude for backward compatibility
          ecAmplitude = point->GetECAmplitude(k, m);
        }
        complex tmatrix = azure2::physics::external_capture_element(
            state.ECGamma(theECMGroup.final_group, theECMGroup.final_level, theECMGroup.final_channel),
            state.SqrtNF(theECMGroup.final_group, theECMGroup.final_level),
            finalLevel->GetECConversionFactor(theECMGroup.final_channel), ecAmplitude);
        if (theECMGroup.channel_capture) {
          int internalChannel = theECMGroup.photon_channel;
          const azure2::InternalRow &chanMGroup = P.internal[theECMGroup.cc_row];
          int maskChannel = (configure().paramMask & Config::IGNORE_ZERO_WIDTHS) ? internalChannel : 0;
          complex umatrix = azure2::physics::channel_capture_factor(
              point->GetSqrtPenetrability(chanMGroup.group, chanMGroup.ch),
              this->GetChannelCaptureBilinear(chanMGroup.group, chanMGroup.ch, chanMGroup.chp, maskChannel));
          tmatrix = tmatrix * umatrix;
        }
        this->AddECTMatrixElement(k, m, tmatrix);
      }
    }
  }
}

/*!
 * Beam energy in the compound system, as FillMatrices takes it.
 */
double AMatrixFunc::InEnergy(const CompoundSystem *compound, Site *point) {
  if (compound->GetPair(compound->GetPairNumFromKey(point->GetEntranceKey()))->GetPType() == 20) {
    int ep = compound->GetPairNumFromKey(point->GetExitKey());
    return point->GetCMEnergy() + compound->GetPair(ep)->GetSepE() + compound->GetPair(ep)->GetExE();
  }
  int aa = compound->GetPairNumFromKey(point->GetEntranceKey());
  return point->GetCMEnergy() + compound->GetPair(aa)->GetSepE() + compound->GetPair(aa)->GetExE();
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
    AZUREProfile::Timer profileTimer(configure().profile, AZUREProfile::kChannelSolve);
    const JGroup *jGroup = compound()->GetJGroup(jGroupNum);
    const std::vector<int> &map = level_active_index_[jGroupNum - 1];
    const int n = (int)active_levels_[jGroupNum - 1].size();
    solve.assign(n, complex(0.0, 0.0));
    for (int la = 1; la < (int)map.size(); la++) {
      const int row = map[la];
      if (row == 0) continue;
      solve[row - 1] = state().Gamma(jGroupNum, la, chNum);
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
          azure2::physics::channel_space::level_solve_from_formal(f, this->FormalSolve(jGroupNum, chNum - 1), solve);
        } else if (same) {
          // Pencil: A g_c = Z X (e_c + z),  (1 - Lambda R) z = Lambda R e_c.
          azure2::physics::channel_space::pencil_channel_solve(f, chan_[jGroupNum - 1], f.pos[chNum - 1], solve);
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
const std::vector<complex> &AMatrixFunc::KernelLambda(Site *point, int j) {
  const JGroup *jg = compound()->GetJGroup(j);
  const bool park = (configure().paramMask & Config::USE_PARK_FORMALISM);
  const bool brune = (configure().paramMask & Config::USE_BRUNE_FORMALISM);
  const int nch = jg->NumChannels();
  kernel_scratch_.resize(nch);
  for (int d = 1; d <= nch; d++) {
    const AChannel *ch = jg->GetChannel(d);
    complex l = point->GetLoElement(j, d);
    const bool particle = park ? (ch->GetRadType() == 'P' && compound()->GetPair(ch->GetPairNum())->GetPType() == 0)
                               : (brune && ch->GetRadType() == 'P');
    if (particle) l += ch->GetBoundaryCondition();
    kernel_scratch_[d - 1] = l;
  }
  return kernel_scratch_;
}

/*!
 * The per-energy part of the channel route for J-group j: the kernel
 * \f$\Lambda\f$, \f$D^{-1}\f$, \f$R = \Gamma^T D^{-1}\Gamma\f$ over the coupled channels and
 * the factors of \f$1 - R\Lambda\f$.  Returns false (the level route is then
 * taken) near a formal pole, where \f$R\f$ is dominated by one level and the
 * Woodbury form loses digits the level matrix keeps, or when \f$1 - R\Lambda\f$
 * is ill-conditioned.
 */
bool AMatrixFunc::SolveChannelSpace(int j, Site *point, double inEnergy) {
  const ChannelSpace::FormalLevels &f = (*formal_[j - 1]);
  ChannelState &st = chan_[j - 1];
  vector_c &lo = st.w3;
  lo.resize(f.nch);
  for (int c = 0; c < f.nch; c++) lo[c] = point->GetLoElement(j, c + 1);
  const azure2::physics::channel_space::Route route =
      azure2::physics::channel_space::solve_channel_space(f, lo.data(), inEnergy, st);
  if (route == azure2::physics::channel_space::Route::kPole) {
    AZUREProfile::Count(configure().profile, AZUREProfile::kLevelRoutePole);
    return false;
  }
  if (!st.defl.empty()) AZUREProfile::Count(configure().profile, AZUREProfile::kDeflated);
  if (route == azure2::physics::channel_space::Route::kPivot) {
    AZUREProfile::Count(configure().profile, AZUREProfile::kLevelRoutePivot);
    return false;
  }
  return true;
}

/*!
 * Solves the bordered system of SolveChannelSpace for the right-hand side
 * (wp over the deflated levels, wz over the channels), in place: on return wp
 * holds y_p and wz holds z = Lambda Gamma^T y.
 */
void AMatrixFunc::SolveBordered(int j, complex *wp, complex *wz) {
  azure2::physics::channel_space::solve_bordered((*formal_[j - 1]), chan_[j - 1], wp, wz);
}

/*!
 * Column b of \f$F = \Gamma^T A_f \Gamma = G^T A G\f$: with the deflated levels p,
 * \f$F_b = \Gamma_p^T y_p + R_b + R z\f$ from the bordered system with
 * right-hand side \f$(\Gamma_{p b}, \Lambda R_b)\f$ (without deflated levels
 * \f$(1 - R\Lambda)^{-1} R_b\f$).
 */
const complex *AMatrixFunc::FColumn(int j, int b) {
  return azure2::physics::channel_space::f_column((*formal_[j - 1]), chan_[j - 1], b);
}

/*!
 * \f$y = A_f w\f$ on the formal route (formal-level coordinates): with
 * \f$t = \Gamma_r^T D_r^{-1} w_r\f$ the bordered system with right-hand side
 * \f$(w_p, \Lambda t)\f$ gives \f$y_p\f$ and \f$z\f$, and
 * \f$y_r = D_r^{-1}(w_r + \Gamma_r z)\f$.
 */
void AMatrixFunc::FormalApply(int j, const complex *w, complex *y) {
  azure2::physics::channel_space::formal_apply((*formal_[j - 1]), chan_[j - 1], w, y);
}

/*!
 * \f$A u\f$ in channel space: \f$V A_f V^T u\f$ on the formal route
 * (FormalApply), \f$P^{-1}u + P^{-1}G z\f$ with \f$(1 - \Lambda R) z = \Lambda G^T P^{-1} u\f$
 * on the pencil route (P = K - E N).
 */
void AMatrixFunc::ApplyA(int j, const complex *u, complex *out) {
  azure2::physics::channel_space::apply_a((*formal_[j - 1]), chan_[j - 1], u, out);
}

/*!
 * \f$u^T A v\f$ in channel space (the channel-capture pathways, whose level
 * vectors are masked and so are not columns of \f$G\f$).
 */
complex AMatrixFunc::ChannelBilinear(int j, const complex *u, const complex *v) {
  return azure2::physics::channel_space::channel_bilinear((*formal_[j - 1]), chan_[j - 1], u, v);
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
    AZUREProfile::Timer profileTimer(configure().profile, AZUREProfile::kChannelSolve);
    azure2::physics::channel_space::formal_channel_solve(f, st, ch0, u);
    st.uvalid[ch0] = 1;
  }
  return u;
}

/*!
 * AZURE_CHANNEL_VERIFY=1: compares the channel-space bilinear forms
 * \f$G^T A G\f$ of J-group j with those of the level matrix, filled here and
 * LU-factored in double and in extended precision, and reports each new worst
 * difference (relative to the largest form of the J-group) on std::cerr.
 */
void AMatrixFunc::VerifyChannelSpace(int j, Site *point, double inEnergy) {
  typedef std::complex<long double> cld;
  this->FillLevelMatrix(j, point, inEnergy);
  azure2::physics::LevelMatrixSolver lu;
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
  azure2::RunDiagnostics &diagnostics = *configure().diagnostics;
  std::lock_guard<std::mutex> lock(diagnostics.channelVerifyMutex);
  if (rel > diagnostics.channelVerifyWorst) {
    diagnostics.channelVerifyWorst = rel;
    std::cerr << "[channel-verify] J-group " << j << " (" << n << " levels, " << m << " channels, "
              << (f.pencil ? "pencil" : f.identity ? "standard" : "formal") << ") at " << inEnergy
              << " MeV: channel vs level route " << rel << " of the largest form; against an extended-precision LU: channel "
              << dChannel / fmax << ", level " << dLevel / fmax << std::endl;
  }
}
