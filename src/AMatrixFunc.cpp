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

namespace {
/// The sink of the adjoints of azure2/physics/CollisionMatrixAdjoint.h and ChannelSolveAdjoint.h over the cells of
/// a gradient accumulator: level r and channel d of a J-group (from 0) are e[r] and gamma[r][d], or, where
/// `levelNumber` is given, those of level number levelNumber[r] (1-based: a block that leaves levels out).
struct LevelSink {
  LevelSink(std::vector<double> &e, std::vector<std::vector<double>> &gamma, const int *levelNumber = nullptr)
      : e(e), gamma(gamma), levelNumber(levelNumber) {}
  /// The level energies and widths of J-group j; `levelOf` maps an active level (1-based) to its level number.
  LevelSink(GradAccum &accum, int j, const std::vector<int> *levelOf = nullptr)
      : LevelSink(accum.e[j - 1], accum.gamma[j - 1], levelOf ? levelOf->data() + 1 : nullptr) {}
  int row(int r) const { return levelNumber ? levelNumber[r] - 1 : r; }
  void energy(int r, double v) { e[row(r)] += v; }
  void width(int r, int d, double v) { gamma[row(r)][d] += v; }
  double *widths(int r) { return &gamma[row(r)][0]; }
  std::vector<double> &e;
  std::vector<std::vector<double>> &gamma;
  const int *levelNumber;
};
}  // namespace

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
 * Reverse mode of one unit's observable: PointTBar (cross section -> T) and TBarAdjoint (T -> level energies and
 * reduced widths); the gradient over the energy tasks (GradEnergyTasks.cpp) sums the T cotangents of the units of a
 * T-matrix group between the two.
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
 * Forward chain (per point, standard A-matrix):
 *   T(k,m) = [tphase] - uphase ( [1] + umatrix )         ([..] only if chNum==chpNum)
 *   umatrix = sum_{l,l'} 2i sqrtP_en sqrtP_ex gamma_en,l gamma_ex,l' A_j[l][l']
 *   A_j = M_j^{-1},   M_j[r][c] = (E_r - inE) delta_rc - sum_ch gamma_r,ch gamma_c,ch L_ch
 * with the cross-section -> T step depending on whether the point is
 * angle-integrated (quadratic |T_temp|^2 sum) or differential (Legendre
 * interference sum + elastic Coulomb interference).
 */

/// A gate of PointTBar no project reaches (why in the message): a broken invariant, never a fallback point.
[[noreturn]] static void UnreachableGate(const char *w) { throw std::logic_error(std::string("PointTBar: ") + w); }
/*!
 * Cross section -> T: the T-matrix cotangents of one point's observable
 * (AMatrixFunc.h).  Reads the T-matrix held by this object and the
 * angle-dependent quantities of `point`.  False: the finite-difference fallback (park_chi2_fallback_test.cpp).
 */
bool AMatrixFunc::PointTBar(Site *point, const azure2::Results &results, double fitBar, TMatrixCotangent &bar,
                            int xsComponent, bool ayCrossSection, const Site &rows, int rowSlot) {
  AZUREProfile::Timer profileTimer(configure().profile, AZUREProfile::kGradTBar);
  const complex I(0.0, 1.0);
  const Site *angles = &rows;  // the angular row (Legendre rows, coefficient order) is the reader's; the rest the site's

  // --- Common gating. ---
  if (configure().paramMask & Config::USE_RMC_FORMALISM) return false;
  if (point->IsAngularDist() && !angles->IsCoefficient()) UnreachableGate("an extrapolation's computed coefficients");
  // xsComponent < 0 asks for A_0 = sigma_int / 4 pi of an angular-distribution point: the angle-integrated
  // cotangent scaled by 1 / 4 pi (GradTargetEffectAdjoint, normalised coefficients under a target effect).
  const bool a0Mode = (xsComponent < 0);
  if (a0Mode) xsComponent = 0;
  if (xsComponent != 0 && (point->IsDifferential() || point->IsPhase()))  // field multipole: Quantities.cpp
    UnreachableGate("an E1/E2 part exists for sigma, sfactor and modified_sfactor only, not differential or phase");

  int aa = compound()->GetPairNumFromKey(point->GetEntranceKey());
  int exitPairNum = compound()->GetPairNumFromKey(point->GetExitKey());
  // Beta-delayed particle emission is supported below; its observable is a plain
  // incoherent sum over pathways, so only the T cotangents differ from the
  // cross-section case. E1/E2 component selection has no meaning for it.
  const bool isBetaDelayed = (compound()->GetPair(aa)->GetPType() == 20);
  if (isBetaDelayed && xsComponent != 0) return false;

  // Active-level compaction unsupported here (original idx == active; ClearMatrices records it).
  if (!all_levels_active_) UnreachableGate("a level outside the R matrix: every level read from <levels> is in it");

  // The exit of the point's exit pair (CalculateTMatrix, called for the point's group first, throws without one).
  const azure2::EntrancePathways &P = compound()->Pathways().Entrance(aa);
  const int ir = P.ExitOf(exitPairNum);
  if (!ir) UnreachableGate("no exit for the point's pair, after CalculateTMatrix of its group");
  const azure2::ExitRow &theDecay = P.exits[ir - 1];

  const double geom = point->GetGeometricalFactor();
  const double i1i2 = compound()->GetPair(aa)->GetI1I2Factor();
  const bool isDiff = point->IsDifferential();
  const bool isPhase = point->IsPhase();
  if (isPhase && aa != exitPairNum) UnreachableGate("an inelastic phase shift (validator and Build refuse it)");
  // The cross section carried by an analyzing-power point is the differential
  // one at the point's angle (GenMatrixFunc computes it before A_y and a
  // target-effect sub-point keeps it); anything else is not supported.
  if (ayCrossSection && (!point->IsAnalyzingPower() || !isDiff || isPhase || xsComponent != 0)) return false;

  // Electron screening (GenMatrixFunc::CalculateCrossSection): a screened observable is f times the bare
  // one, so the bare one's cotangent is f fitBar, and dY/dU_e = k Y per MeV, 1e-6 k Y per eV, with Y the
  // screened value this point holds from its forward calculation.  A phase shift and an A_y value are
  // not screened (an A_y sub-point's cross section is: ayCrossSection).
  bar.screen.clear();
  if (point->GetScreeningParam() > 0 && point->GetScreeningK() != 0.0 && !isPhase && !angles->IsCoefficient() &&
      !(point->IsAnalyzingPower() && !ayCrossSection)) {
    double y = point->ValueFor(*angles, rowSlot, results);
    if (xsComponent == 1) y = point->GetFitE1CrossSection(results);
    else if (xsComponent == 2) y = point->GetFitE2CrossSection(results);
    bar.screen.push_back({point->GetScreeningParam(), fitBar * 1.0e-6 * point->GetScreeningK() * y});
    fitBar *= point->ScreeningFactor();
  }

  // === Cross-section -> T:  build the T-matrix cotangents. ===
  // tBar  : cotangents of the internal T(k,m).
  // ecBar : cotangents of the external-capture T(k,ecm).
  int nK = theDecay.num_groups;
  std::vector<vector_c> &tBar = bar.t, &ecBar = bar.ec;
  tBar.assign(nK, vector_c());
  ecBar.assign(nK, vector_c());
  for (int k = 1; k <= nK; k++) {
    tBar[k - 1].assign(P.Group(theDecay, k).num_internal, complex(0.0, 0.0));
    ecBar[k - 1].assign(P.Group(theDecay, k).num_external, complex(0.0, 0.0));
  }

  // Keep only the contributions feeding the requested cross-section component
  // (0 = full, 1 = E1, 2 = E2); E1/E2 select electric channels of that L.
  auto includeInternal = [&](const azure2::InternalRow *mg) -> bool {
    if (xsComponent == 0) return true;
    const AChannel *ex = compound()->GetJGroup(mg->group)->GetChannel(mg->chp);
    return ex->GetRadType() == 'E' && ex->GetL() == xsComponent;
  };
  auto includeEC = [&](const azure2::ExternalRow *ecg) -> bool {
    if (xsComponent == 0) return true;
    return ecg->rad_type == 'E' && ecg->mult == xsComponent;
  };

  // The angle-integrated cotangent as a function of its own cotangent `bar`: the chain below calls it with
  // fitBar; a normalised Legendre coefficient a_k = N_k / sigma_int calls it with -fitBar a_k / sigma_int.
  // symmetryFactor: the (1+delta_12) of an identical entrance pair where the forward value carries it
  // (the reported cross section), 1 where it does not (sigma_int of the angular coefficients).
  auto integratedBar = [&](double bar, double symmetryFactor) {
    // ---- Angle-integrated: model = (1/100) sum_k sum_temp w_temp |T_temp|^2,
    //      w_temp = geom (2J+1) I1I2 (1+delta_12),  T_temp = sum_m T(k,m) sharing (J,l,l'). ----
    for (int k = 1; k <= nK; k++) {
      // Mirror the angle-integrated UPOS over-counting guard.
      if (aa != exitPairNum && compound()->GetPair(exitPairNum)->GetPType() == 0) {
        if (P.Group(theDecay, k).sp != P.Group(theDecay, k).sp2) continue;
      }
      int numM = P.Group(theDecay, k).num_internal;
      int numEC = P.Group(theDecay, k).num_external;
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
        const azure2::InternalRow *mg = &P.Internal(P.Group(theDecay, k), m);
        if (!includeInternal(mg)) continue;
        const JGroup *jg = compound()->GetJGroup(mg->group);
        int gi = findGroup(jg->GetJ(),
                           jg->GetChannel(mg->ch)->GetL(),
                           jg->GetChannel(mg->chp)->GetL());
        groups[gi].tsum += this->GetTMatrixElement(k, m);
        mGroupIdx[m] = gi;
      }
      // External-capture pathways T(k,ecm), grouped by (J, L, mult) and merged
      // with the matching internal group (forward: same temp T-matrix list).
      for (int m = 1; m <= numEC; m++) {
        const azure2::ExternalRow *ecg = &P.External(P.Group(theDecay, k), m);
        if (!includeEC(ecg)) continue;
        int gi = findGroup(ecg->J, ecg->l, ecg->mult);
        groups[gi].tsum += this->GetECTMatrixElement(k, m);
        ecGroupIdx[m] = gi;
      }
      std::vector<complex> tBarGroup(groups.size());
      for (int g = 0; g < (int)groups.size(); g++) {
        double w = geom * (2.0 * groups[g].j + 1.0) * i1i2 * symmetryFactor;
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
    if (xsComponent != 0 || isDiff) UnreachableGate("sigma_total is angle-integrated and takes no multipole");
    for (int k = 1; k <= nK; k++) {
      if (P.Group(theDecay, k).CrossSpin()) continue;
      for (int m = 1; m <= P.Group(theDecay, k).num_internal; m++) {
        const azure2::InternalRow *mg = &P.Internal(P.Group(theDecay, k), m);
        if (mg->ch != mg->chp) continue;
        const double jValue = compound()->GetJGroup(mg->group)->GetJ();
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
      for (int m = 1; m <= P.Group(theDecay, k).num_internal; m++)
        tBar[k - 1][m - 1] += fitBar * (2.0 * 25.0 / 100.0) * this->GetTMatrixElement(k, m);
  } else if (point->IsAnalyzingPower() && !ayCrossSection) {
    // ---- Vector analyzing power. A_y is built from the channel-spin amplitude
    //      matrix, which is linear in T(k,m) with the coefficients AddPathway
    //      applies, so the whole derivative is a contraction against those same
    //      coefficients. The Coulomb amplitude carries no parameter dependence
    //      and drops out. ----
    if (compound()->GetPair(aa)->GetPType() != 0) UnreachableGate("an A_y needs a particle entrance (Build refuses)");

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
      const int maxL = angles->GetMaxLOrder();
      if (!theDecay.has_capture_ay) return false;  // terms not requested: the forward pass refused the point
      const int nTerms = theDecay.num_capture_ay;
      if (nTerms == 0) return false;

      const double xCos = cos(angles->SlotAngle(rowSlot) * pi / 180.0);  // the reader's angle at the slot
      std::vector<double> assocP(maxL + 1, 0.0);
      for (int k = 1; k <= maxL; k++) assocP[k] = AngCoeff::LegendreP1(k, xCos);

      auto element = [&](int kg, int p, bool isEC) -> complex {
        return isEC ? this->GetECTMatrixElement(kg, p)
                    : this->GetTMatrixElement(kg, p);
      };
      double num = 0.0, den = 0.0;
      for (int i = 0; i < nTerms; i++) {
        const azure2::CaptureAyTerm *t = &P.capture_ay[theDecay.first_capture_ay + i];
        if (t->kOrder > maxL) continue;
        complex prod = element(t->kGroup1, t->path1, t->isEC1) * conj(element(t->kGroup2, t->path2, t->isEC2));
        if (t->ak != 0.0) den += t->ak * real(prod) * angles->RowP(rowSlot, t->kOrder);
        if (t->bk != 0.0) num += t->bk * (-imag(prod)) * assocP[t->kOrder];
      }
      if (fabs(den) < 1.e-300) return false;
      const double ay = num / den;

      const complex I2(0.0, 1.0);
      for (int i = 0; i < nTerms; i++) {
        const azure2::CaptureAyTerm *t = &P.capture_ay[theDecay.first_capture_ay + i];
        if (t->kOrder > maxL) continue;
        const complex T1 = element(t->kGroup1, t->path1, t->isEC1);
        const complex T2 = element(t->kGroup2, t->path2, t->isEC2);
        const double pk = angles->RowP(rowSlot, t->kOrder);
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
      Polarization::AmplitudeMatrix M(compound(), angles->SlotAngle(rowSlot), aa, exitPairNum);  // the reader's, at the slot
      for (int k = 1; k <= nK; k++) {
        if (P.Group(theDecay, k).CrossSpin()) continue;  // as the forward sum
        for (int m = 1; m <= P.Group(theDecay, k).num_internal; m++) {
          const azure2::InternalRow *mg = &P.Internal(P.Group(theDecay, k), m);
          M.AddPathway(mg->group, mg->ch, mg->chp,
                       this->GetTMatrixElement(k, m));
        }
      }
      complex fdAy(0.0, 0.0), fxAy(0.0, 0.0);  // the reader's Coulomb amplitudes at the slot's energy
      if (aa == exitPairNum) angles->CoulombAt(rowSlot, fdAy, fxAy), M.AddCoulomb(fdAy, fxAy);
      if (M.size() == 0) return false;

      const std::vector<complex> bar = M.AnalyzingPowerBar();
      for (int k = 1; k <= nK; k++) {
        if (P.Group(theDecay, k).CrossSpin()) continue;
        for (int m = 1; m <= P.Group(theDecay, k).num_internal; m++) {
          const azure2::InternalRow *mg = &P.Internal(P.Group(theDecay, k), m);
          if (!includeInternal(mg)) continue;
          tBar[k - 1][m - 1] += fitBar * M.PathwayAdjoint(mg->group, mg->ch, mg->chp, bar);
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
      for (int m = 1; m <= P.Group(theDecay, k).num_internal; m++) {
        const azure2::InternalRow *mg = &P.Internal(P.Group(theDecay, k), m);
        const JGroup *jg = compound()->GetJGroup(mg->group);
        const AChannel *en = jg->GetChannel(mg->ch);
        const AChannel *ex = jg->GetChannel(mg->chp);
        if (jg->GetJ() == segJ && en->GetL() == segL && en == ex) {
          complex expCP = point->GetExpCoulombPhase(mg->group, mg->ch);
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
  } else if (!isDiff && (!angles->IsCoefficient() || a0Mode)) {
    // The reported angle-integrated cross section carries (1+delta_12) of an identical entrance pair
    // (GenMatrixFunc); the A_0 of an angular-coefficient point is the unscaled integral.
    integratedBar(a0Mode ? fitBar / (4.0 * pi) : fitBar,
                  a0Mode ? 1.0 : compound()->GetPair(aa)->GetEntranceSymmetryFactor());
  } else {
    // ---- Differential: model = (Re(CT)+Re(RT)+Re(IT))/100. ----
    //   RT = (geom I1I2 rtFactor/pi) Re( sum_{kL,inter} weight T1 conj(T2) P_L )
    //   IT = (geom itFactor/sqrt(pi)) Re( i sum_{k,m: en==ex} statSpin coulombAmp conj(T) P_l )  (elastic)
    //   CT = const (no T dependence).
    // The standard case (entrance==exit pair, or gamma exit) with external-capture interferences, and the
    // inelastic-particle UPOS / normal-angular-distribution case.
    const PPair *exitPairDiff = compound()->GetPair(exitPairNum);
    // Identical entrance pair: x4 / x2 for elastic, (1+delta_12) otherwise,
    // mirroring GenMatrixFunc::CalculateCrossSection.
    double rtFactor = 1.0, itFactor = 1.0;
    if (compound()->GetPair(aa)->IsIdentical()) {
      if (aa == exitPairNum) {
        rtFactor = 4.0;
        itFactor = 2.0;
      } else {
        rtFactor = compound()->GetPair(aa)->GetEntranceSymmetryFactor();
      }
    }
    double cR = geom * i1i2 * rtFactor / (pi * 100.0);
    // An angular-coefficient point (Segment isDiff 9): the forward value is
    //   a_k = (4 geom i1i2 / (100 sigma_int)) sum_{L=k} z T1 T2*   (normalised)   or   A_k = a_k sigma_int / 4 pi,
    // i.e. the differential model's Legendre coefficient of order k, coeff_k = cR sum_{L=k} z T1 T2*, times
    // 4 pi / (rtFactor sigma_int) (normalised) or 1 / rtFactor (absolute).  So the differential machinery below
    // with a Kronecker delta in place of P_l(cos theta) gives the N_k part; the normalised one then adds the
    // sigma_int cotangent -fitBar a_k / sigma_int (integratedBar).  No Coulomb terms: elastic pairs are refused.
    const bool coef = angles->IsCoefficient();
    const int kOrd = coef ? angles->GetCoefficientOrder() : -1;
    const double sigmaInt = coef ? 4.0 * pi * point->GetAngDistA0(results) : 0.0;
    if (coef && point->IsUPOS()) UnreachableGate("UPOS coefficients: quantity legendre takes no secondary_gamma");
    if (coef && !(sigmaInt > 0.0)) return false;  // a vanishing integrated cross section (parameter values)
    const double dfac = coef ? (angles->GetCoefficientAbsolute() ? 1.0 / rtFactor : 4.0 * pi / (rtFactor * sigmaInt)) : 0.0;
    auto pl = [&](int l) -> double { return coef ? ((l == kOrd) ? dfac : 0.0) : angles->RowP(rowSlot, l); };

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
    if (!(inelasticParticle && point->IsUPOS()) && this->TBarProducts(P, theDecay, inelasticParticle)) {
      const complex scale = fitBar * cR;
      const azure2::InterferencePathway *paths = P.interference_pathways.data() + theDecay.first_interference_pathway;
      const size_t numPaths = theDecay.num_interference_pathways;
      const int nl = theDecay.num_orders;
      tbar_legendre_.assign(nl, 0.0);
      const int nP = std::min(nl, angles->GetMaxLOrder() + 1);
      for (int l = 0; l < nP; l++) tbar_legendre_[l] = pl(l);
      for (size_t p = 0; p < numPaths; p++) {
        const complex *q = &tbar_q_[p * nl];
        complex s(0.0, 0.0);
        for (int l = 0; l < nl; l++) s += q[l] * tbar_legendre_[l];
        if (s != complex(0.0, 0.0)) Tadd(paths[p].k, paths[p].m, paths[p].rad, scale * s);
      }
    } else
    for (int inter = 0; inter < theDecay.num_interferences; inter++) {
      const azure2::FlatInterference &itf = P.interferences[theDecay.first_interference + inter];
      int K = itf.k;
      int lOrder = itf.l;
      double legP = pl(lOrder);
      const azure2::SpinGroupRow &kg = P.Group(theDecay, K);
      {
        const std::string type{itf.r1, itf.r2};  // "RR","ER","RE","EE": BuildPathwayTable writes 'R' or 'E' only
        char rad1 = type[0], rad2 = type[1];  // 'R' internal, 'E' external
        int m1 = itf.m1, m2 = itf.m2;

        if (inelasticParticle) {
          // Inelastic particle exit: pick the normal or UPOS sub-case.
          double sp1 = kg.sp, sp2 = kg.sp2;
          int lp1 = compound()->GetJGroup(P.Internal(kg, m1).group)->GetChannel(P.Internal(kg, m1).chp)->GetL();
          int lp2 = compound()->GetJGroup(P.Internal(kg, m2).group)->GetChannel(P.Internal(kg, m2).chp)->GetL();
          if (sp1 == sp2 && !point->IsUPOS()) {
            addQuad(K, m1, rad1, m2, rad2, cR * itf.z * legP);
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
                // unconditionally, as in GenMatrixFunc (checklist t3-24)
                R_L = (R_L + 2.0 * delta * R_LLp + delta * delta * R_LpLp) / (1.0 + delta * delta);
              }
            }
            double W = cR * itf.z_upos * std::pow(2.0 * lOrder + 1.0, 0.5) / 4.0 * R_L * legP;
            addQuad(K, m1, 'R', m2, 'R', W);
          }
        } else {
          // Standard case (entrance==exit pair, or gamma exit).
          addQuad(K, m1, rad1, m2, rad2, cR * itf.z * legP);
        }
      }
    }

    if (coef && !angles->GetCoefficientAbsolute()) {
      // a_k = N_k / sigma_int: the sigma_int cotangent through the angle-integrated form
      integratedBar(-fitBar * point->ValueFor(*angles, rowSlot, results) / sigmaInt, 1.0);
    }

    // f_C of channel spin s, the reader's at its slot (Site::CoulombAt): f_d + (-1)^s f_x (f_x = 0 but with spin exchange)
    complex fd(0.0, 0.0), fx(0.0, 0.0);
    if (!coef && aa == exitPairNum) angles->CoulombAt(rowSlot, fd, fx);
    auto readerCoulomb = [&](double s) { return (((std::lround(2.0 * s) / 2) % 2 == 0) ? 1.0 : -1.0) * fx + fd; };
    if (!coef && aa == exitPairNum && theDecay.has_interferences) {  // elastic Coulomb-nuclear interference
      const complex cI = I * geom * itFactor / (std::sqrt(pi) * 100.0);
      for (int n = 0; n < theDecay.num_diagonal; n++) {
        const azure2::DiagonalPathway &d = P.diagonal[theDecay.first_diagonal + n];
        const complex B = cI * d.statSpin * readerCoulomb(d.s) * angles->RowP(rowSlot, d.l);
        tBar[d.k - 1][d.m - 1] += fitBar * B;
      }
    } else if (!coef && aa == exitPairNum) {  // elastic Coulomb-nuclear interference
      complex cI = I * geom * itFactor / (std::sqrt(pi) * 100.0);
      for (int k = 1; k <= nK; k++) {
        for (int m = 1; m <= P.Group(theDecay, k).num_internal; m++) {
          const azure2::InternalRow *mg = &P.Internal(P.Group(theDecay, k), m);
          const AChannel *en = compound()->GetJGroup(mg->group)->GetChannel(mg->ch);
          const AChannel *ex = compound()->GetJGroup(mg->group)->GetChannel(mg->chp);
          if (en == ex) {
            int l = en->GetL();
            complex B = cI * mg->stat_spin_factor * readerCoulomb(en->GetS()) * angles->RowP(rowSlot, l);
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
 * T -> level energies and reduced widths (AMatrixFunc.h), at the energy of
 * `energyPoint`, whose T-matrix this object holds.
 */
void AMatrixFunc::TBarAdjoint(Site *energyPoint, const TMatrixCotangent &bar, GradAccum &accum,
                              const vector_matrix_r *shiftDeriv) {
  AZUREProfile::Timer profileTimer(configure().profile, AZUREProfile::kGradBackprop);
  const complex I(0.0, 1.0);
  for (const auto &s : bar.screen) accum.AddScreening(s.first, s.second);
  Site *point = energyPoint;
  // Brune formalism: the boundary/shift terms change the M-construction kernel
  // (handled in the structural gamma step via cached shift values) and add
  // explicit E_lambda dependence through S(E_lambda) (handled in the structural
  // energy step using shiftDeriv = dS/dE).
  // Park: the level matrix is (E_r - E) J[r][c] - sum_d g_r g_c (L_d + B_d - S_r,d),
  // with the overlap J of Park's Eqs. (9) and (28) depending on the widths and,
  // through dS/dE at E_lambda, on the level energies.  Its derivatives come from
  // the per-level values CompoundSystem::CalcShiftFunctions stores (S, S', S'').
  const bool park = (configure().paramMask & Config::USE_PARK_FORMALISM);
  int aa = compound()->GetPairNumFromKey(point->GetEntranceKey());
  int exitPairNum = compound()->GetPairNumFromKey(point->GetExitKey());
  // The exit of the point's exit pair (CalculateTMatrix has refused a point without one).
  const azure2::EntrancePathways &P = compound()->Pathways().Entrance(aa);
  const int ir = P.ExitOf(exitPairNum);
  if (!ir) return;
  const azure2::ExitRow &theDecay = P.exits[ir - 1];
  bool hasEC = false;
  for (int k = 1; k <= theDecay.num_groups; k++)
    if (P.Group(theDecay, k).num_external > 0) {
      hasEC = true;
      break;
    }
  const int nK = theDecay.num_groups;
  const std::vector<vector_c> &tBar = bar.t;
  const std::vector<vector_c> &ecBar = bar.ec;

  // Without external capture the level matrix enters only through bilinear
  // forms, whose adjoint is formed from rank-one terms (BilinearAdjoint).
  // AZURE_GRAD_DENSE=1 (the run's Config::environment) keeps the dense A^H Abar A^H products below.
  if (!hasEC && !configure().environment.denseAdjoint) {
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
      int numEC = P.Group(theDecay, k).num_external;
      for (int m = 1; m <= numEC; m++) {
        complex eb = ecBar[k - 1][m - 1];
        if (eb == complex(0.0, 0.0)) continue;
        const azure2::ExternalRow *ecg = &P.External(P.Group(theDecay, k), m);
        int jf = ecg->final_group;
        int lf = ecg->final_level;
        int fc = ecg->final_channel;
        const ALevel *finalLevel = compound()->GetJGroup(jf)->GetLevel(lf);
        const JGroup *jgf = compound()->GetJGroup(jf);
        double sqrtNF = state().SqrtNF(jf, lf);
        double conv = finalLevel->GetECConversionFactor(fc);
        complex ecAmp = useEC ? point->GetECAmplitudeWithShift(k, m, compound(), configure())
                              : point->GetECAmplitude(k, m);

        // base = gamma_fc * sqrtNF * conv * ecAmp  (the part without umatrix_cc).
        double gamma_fc = state().ECGamma(jf, lf, fc);  // the final state's Brune amplitude (cnuc_split.md section 8)
        complex base = gamma_fc * sqrtNF * conv * ecAmp;

        // "effective amplitude" carrying the umatrix_cc factor for the
        // gamma_fc / sqrtNF derivatives (= ecAmp for direct capture).
        complex effAmp = ecAmp;

        if (ecg->channel_capture) {
          int internalChannel = ecg->photon_channel;
          const azure2::InternalRow *ccMg = &P.internal[ecg->cc_row];
          int jcc = ccMg->group;
          int ccChN = ccMg->ch;
          int ccChpN = ccMg->chp;
          const JGroup *jgcc = compound()->GetJGroup(jcc);
          int ccLevels = jgcc->NumLevels();
          complex sqrtPcc = point->GetSqrtPenetrability(jcc, ccChN);  // single sqrtP

          // Forward umatrix_cc (mirrors the IGNORE_ZERO_WIDTHS skips).
          complex umatrix_cc(0.0, 0.0);
          for (int la = 1; la <= ccLevels; la++) {
            if (!jgcc->GetLevel(la)->IsInRMatrix()) continue;
            if (internalChannel && ignoreZero &&
                state().ZeroInputWidth(jcc, la, internalChannel)) continue;
            for (int lap = 1; lap <= ccLevels; lap++) {
              if (!jgcc->GetLevel(lap)->IsInRMatrix()) continue;
              if (internalChannel && ignoreZero &&
                  state().ZeroInputWidth(jcc, lap, internalChannel)) continue;
              umatrix_cc += I2 * sqrtPcc *
                  state().Gamma(jcc, la, ccChN) *
                  state().Gamma(jcc, lap, ccChpN) *
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
                state().ZeroInputWidth(jcc, la, internalChannel)) continue;
            double gEn = state().Gamma(jcc, la, ccChN);
            for (int lap = 1; lap <= ccLevels; lap++) {
              if (!jgcc->GetLevel(lap)->IsInRMatrix()) continue;
              if (internalChannel && ignoreZero &&
                  state().ZeroInputWidth(jcc, lap, internalChannel)) continue;
              double gEx = state().Gamma(jcc, lap, ccChpN);
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
        //   nFSum = 1 + sum_d A_d gamma_d^2,  A_d = 2 chRad redMass uconv/hbarc^2 NF_d (physics/ExternalCapture.h),
        //   sqrtNF = nFSum^{-1/2},  d sqrtNF / d gamma_d = -A_d gamma_d sqrtNF^3.
        // Park: sqrtNF = (J + sum_d A_d gamma_d^2)^{-1/2} with J = 1 - sum_d gamma_d^2 S'_d
        // over the particle channels (CompoundSystem::CalcShiftFunctions), so every particle
        // channel contributes  d sqrtNF / d gamma_d = (S'_d - A_d) gamma_d sqrtNF^3.
        complex pre = gamma_fc * conv * effAmp;  // d T_ec / d sqrtNF
        int numNF = finalLevel->NumNFIntegrals();
        double sqrtNF3 = sqrtNF * sqrtNF * sqrtNF;
        const bool parkFinal = park && finalLevel->IsInRMatrix();
        int numD = parkFinal ? jgf->NumChannels() : numNF;
        for (int d = 1; d <= numD; d++) {
          const PPair *pd = compound()->GetPair(jgf->GetChannel(d)->GetPairNum());
          double Ad = (d <= numNF) ? azure2::physics::bound_state_coefficient(pd->GetChRad(), pd->GetRedMass()) *
                                         finalLevel->GetNFIntegral(d)
                                   : 0.0;
          if (parkFinal && jgf->GetChannel(d)->GetRadType() == 'P' && pd->GetPType() == 0)
            Ad -= state().ShiftDerivative(jf, lf, d);
          if (Ad == 0.0) continue;
          double gd = state().ECGamma(jf, lf, d);
          double dSqrtNF = -Ad * gd * sqrtNF3;
          accum.AddGamma(jf, lf, d, std::real(conj(eb) * (pre * dSqrtNF)));
        }
      }
    }
  }

  // === T -> U -> A (+ explicit gamma), back-propagating tBar over all pathways. ===
  for (int k = 1; k <= nK; k++) {
    int numM = P.Group(theDecay, k).num_internal;
    for (int m = 1; m <= numM; m++) {
      complex tb = tBar[k - 1][m - 1];
      if (tb == complex(0.0, 0.0)) continue;
      const azure2::InternalRow *mg = &P.Internal(P.Group(theDecay, k), m);
      int jNum = mg->group;
      int chNum = mg->ch;
      int chpNum = mg->chp;
      const JGroup *jg = compound()->GetJGroup(jNum);

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
        double gammaEn = state().Gamma(jNum, la, chNum);
        for (int lap = 1; lap <= numLevels; lap++) {
          if (!jg->GetLevel(lap)->IsInRMatrix()) continue;
          double gammaEx = state().Gamma(jNum, lap, chpNum);
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
    const JGroup *jg = compound()->GetJGroup(jNum);
    const matrix_c &A = this->GetAMatrix(jNum);
    const matrix_c &Ab = aBar[jNum - 1];
    int n = (int)A.size();
    if (n == 0) continue;

    // Mbar = -A^H * Abar * A^H.
    const matrix_c Mbar = azure2::physics::inverse_adjoint(A, Ab);
    this->LevelMatrixAdjoint(point, jNum, Mbar, inEnergy, accum, shiftDeriv);
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
 * The level-matrix part of the adjoint of J-group jNum: given the cotangent
 * Mbar of its level matrix, the structural energy and reduced-width gradients
 * of the standard, Brune and Park level matrices (exactly mirroring
 * FillMatrices).
 */
void AMatrixFunc::LevelMatrixAdjoint(Site *energyPoint, int jNum, const matrix_c &Mbar, double inEnergy,
                                     GradAccum &accum, const vector_matrix_r *shiftDeriv, const complex *lambda,
                                     const std::vector<int> *levelOf) {
  AZUREProfile::Timer profileTimer(configure().profile, AZUREProfile::kGradLevelMatrix);
  const bool brune = (configure().paramMask & Config::USE_BRUNE_FORMALISM);
  const bool park = (configure().paramMask & Config::USE_PARK_FORMALISM);
  const JGroup *jg = compound()->GetJGroup(jNum);
  const int n = (int)Mbar.size();
  if (n == 0) return;
  const int nch = jg->NumChannels();
  const LevelState &state = this->state();
  LevelAdjointGather &g = lma_gather_;
  azure2::physics::LevelMatrixAdjointTables t;
  t.n = n;
  t.nc = nch;
  t.brune = brune;
  t.park = park;
  // The channels whose kernel carries the level's shift: the particle channels of Park, the 'P' channels of Brune.
  g.shiftChannel.resize(nch);
  for (int d = 0; d < nch; d++) {
    const AChannel *ch = (park || brune) ? jg->GetChannel(d + 1) : nullptr;
    g.shiftChannel[d] = !ch ? false
                        : park ? (ch->GetRadType() == 'P' && compound()->GetPair(ch->GetPairNum())->GetPType() == 0)
                               : (ch->GetRadType() == 'P');
  }
  t.shift_channel = g.shiftChannel.data();
  // The channel kernel: the caller's, or the Lo-matrix element of the point with B_c of the shift channels.
  if (lambda) {
    t.kernel = lambda;
  } else {
    g.kernel.resize(nch);
    g.boundary.resize(nch);
    for (int d = 0; d < nch; d++) {
      g.kernel[d] = energyPoint->GetLoElement(jNum, d + 1);
      g.boundary[d] = g.shiftChannel[d] ? jg->GetChannel(d + 1)->GetBoundaryCondition() : 0.0;
    }
    t.kernel = g.kernel.data();
    t.boundary = g.boundary.data();
  }
  // The level data.  Every level: the state's own arrays, and for the standard and Brune matrices E, G, S, S G and
  // 1/(E_r - E_c) as kept per parameter set.  A THM block that leaves levels out of the R matrix (levelOf): gathered
  // here over its active levels.
  if (levelOf) {
    ChannelSpace::LevelArrays &la = lma_active_;
    const size_t cells = (size_t)n * nch;
    la.n = n;
    la.nch = nch;
    la.E.resize(n);
    la.G.resize(cells);
    la.S.resize(cells);
    g.d1.resize(cells);
    g.d2.resize(cells);
    g.parkNorm.resize(n);
    if (!park) la.SG.resize(cells);
    for (int r = 0; r < n; r++) {
      const int level = (*levelOf)[r + 1];
      la.E[r] = state.E(jNum, level);
      g.parkNorm[r] = state.ParkNorm(jNum, level);
      for (int ch = 0; ch < nch; ch++) {
        const size_t i = (size_t)r * nch + ch;
        la.G[i] = state.Gamma(jNum, level, ch + 1);
        la.S[i] = state.Shift(jNum, level, ch + 1);
        if (!park) la.SG[i] = la.S[i] * la.G[i];
        g.d1[i] = state.ShiftDerivative(jNum, level, ch + 1);
        g.d2[i] = state.ShiftSecondDerivative(jNum, level, ch + 1);
      }
    }
    if (!park) {
      la.inv.assign((size_t)n * n, 0.0);
      for (int r = 0; r < n; r++)
        for (int c = 0; c < n; c++)
          if (r != c) la.inv[(size_t)r * n + c] = azure2::physics::inverse_level_gap(la.E[r], la.E[c]);
    }
    t.E = la.E.data();
    t.gamma = la.G.data();
    t.shift = la.S.data();
    t.shift_gamma = la.SG.data();
    t.inv_gap = la.inv.data();
    t.shift_d1 = g.d1.data();
    t.shift_d2 = g.d2.data();
    t.park_norm = g.parkNorm.data();
  } else {
    t.E = state.Energies(jNum);
    t.gamma = state.Gammas(jNum);
    t.shift = state.Shifts(jNum);
    t.shift_d1 = state.ShiftDerivatives(jNum);
    t.shift_d2 = state.ShiftSecondDerivatives(jNum);
    t.park_norm = state.ParkNorms(jNum);
    if (!park) {
      const ChannelSpace::LevelArrays &la = ChannelSpace::CurrentLevelArrays(ChannelMemo(), jg, state, jNum);
      t.E = la.E.data();
      t.gamma = la.G.data();
      t.shift = la.S.data();
      t.shift_gamma = la.SG.data();
      t.inv_gap = la.inv.data();
    }
  }
  if (brune && shiftDeriv) t.brune_shift_d1 = &(*shiftDeriv)[jNum - 1];
  azure2::physics::level_matrix_adjoint(t, Mbar, inEnergy, lma_scratch_, LevelSink(accum, jNum, levelOf));
}

/*!
 * Back-propagates the cotangent Abar of the level matrix A of J group jNum
 * (indexed by the active, in-R-matrix levels) to the level energies and
 * reduced widths: the inverse lemma Mbar = -A^H Abar A^H (dense), then the
 * structural part (LevelMatrixAdjoint).  levelOf[r] is the level number of
 * active level r (1-based); null when every level is in the R matrix.  Used
 * by the THM (HOES) adjoint, THMMatrixFunc::THMPointAdjoint.
 */
void AMatrixFunc::LevelMatrixABarAdjoint(Site *point, int jNum, const matrix_c &aBar, double inEnergy,
                                         const vector_matrix_r *shiftDeriv, GradAccum &accum,
                                         const std::vector<int> *levelOf) {
  const matrix_c &A = this->GetAMatrix(jNum);
  const matrix_c &Ab = aBar;
  int n = (int)A.size();
  if (n == 0) return;

  // Mbar = -A^H * Abar * A^H.
  const matrix_c Mbar = azure2::physics::inverse_adjoint(A, Ab);
  this->LevelMatrixAdjoint(point, jNum, Mbar, inEnergy, accum, shiftDeriv, nullptr, levelOf);
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
 * T -> bilinear forms.  Each internal pathway T(k,m) of J-group j is
 *   T = [tphase] - uphase (1 + 2i sqrtP_en sqrtP_ex b),   b = g_en^T A g_ex,
 * so its cotangent tbar gives the bilinear form the cotangent
 *   beta = conj(2i sqrtP_en sqrtP_ex) conj(-uphase) tbar.
 * Pathways of one J-group sharing a channel pair add their betas.
 */
bool AMatrixFunc::CollectBilinearTerms(Site *energyPoint, const TMatrixCotangent &bar,
                                       std::vector<std::vector<BilinearTerm>> &terms) {
  Site *point = energyPoint;
  const complex I(0.0, 1.0);
  int aa = compound()->GetPairNumFromKey(point->GetEntranceKey());
  int exitPairNum = compound()->GetPairNumFromKey(point->GetExitKey());
  // The exit of the point's exit pair (CalculateTMatrix has refused a point without one).
  const azure2::EntrancePathways &P = compound()->Pathways().Entrance(aa);
  const int ir = P.ExitOf(exitPairNum);
  if (!ir) return false;
  const azure2::ExitRow &theDecay = P.exits[ir - 1];
  terms.assign(compound()->NumJGroups(), std::vector<BilinearTerm>());
  for (int k = 1; k <= theDecay.num_groups; k++)
    if (P.Group(theDecay, k).num_external > 0) return false;
  for (int k = 1; k <= theDecay.num_groups && k <= (int)bar.t.size(); k++) {
    const int numM = P.Group(theDecay, k).num_internal;
    for (int m = 1; m <= numM && m <= (int)bar.t[k - 1].size(); m++) {
      const complex tb = bar.t[k - 1][m - 1];
      if (tb == complex(0.0, 0.0)) continue;
      const azure2::InternalRow *mg = &P.Internal(P.Group(theDecay, k), m);
      const int jNum = mg->group;
      const int chNum = mg->ch;
      const int chpNum = mg->chp;
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
void AMatrixFunc::BilinearAdjoint(Site *energyPoint, int jNum, const std::vector<BilinearTerm> &terms,
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
    const JGroup *jg = compound()->GetJGroup(jNum);
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
    const std::vector<double> &G = ChannelSpace::CurrentLevelArrays(ChannelMemo(), jg, state(), jNum).G;
    vector_c ga(nch), gb(nch);
    LevelSink sink(accum, jNum);
    for (const BilinearTerm &t : terms) {
      const vector_c &uEn = this->GetRawChannelSolve(jNum, t.ch);
      const vector_c &uEx = this->GetRawChannelSolve(jNum, t.chp);
      azure2::physics::add_explicit_widths(t, uEn.data(), uEx.data(), n, sink);
      azure2::physics::deferred_bilinear_adjoint(t, uEn.data(), uEx.data(), G.data(), lam.data(), inEnergy, n, nch,
                                                 K.data(), N.data(), ga, gb, sink);
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
  LevelSink sink(accum, jNum);
  for (const BilinearTerm &t : terms) {
    const vector_c &uEn = this->GetRawChannelSolve(jNum, t.ch);
    const vector_c &uEx = this->GetRawChannelSolve(jNum, t.chp);
    azure2::physics::add_explicit_widths(t, uEn.data(), uEx.data(), n, sink);
    azure2::physics::add_rank_one_cotangent(t.beta, uEn.data(), uEx.data(), n, Mbar);
  }
  this->LevelMatrixAdjoint(energyPoint, jNum, Mbar, this->InEnergy(energyPoint), accum, shiftDeriv);
}

/*!
 * One diagonal bilinear form g_ch^T A g_ch with cotangent mBar (the scattering length's M_cc).  The
 * general BilinearAdjoint takes every level of the group as active, so groups with a level left out
 * of the R-matrix are refused here, as is the RMC formalism.
 */
bool AMatrixFunc::BilinearAdjoint(Site *point, int jGroupNum, int chNum, complex mBar, GradAccum &accum,
                                  const vector_matrix_r *shiftDeriv) {
  if (configure().paramMask & Config::USE_RMC_FORMALISM) return false;
  if (jGroupNum < 1 || jGroupNum > (int)solvers_.size()) return false;
  // Solved in channel space (FillMatrices) or as a level matrix (InvertMatrices).
  if (route_[jGroupNum - 1] != kRouteChannel && !solver_valid_[jGroupNum - 1]) return false;
  const JGroup *jg = compound()->GetJGroup(jGroupNum);
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
  LevelSink sink(eOut, gOut, toFormal ? nullptr : lev.data());
  for (const BilinearTerm &t : terms) {
    const int a = t.ch - 1, b = t.chp - 1;
    const vector_c &ua = this->FormalSolve(j, a);
    const vector_c &ub = this->FormalSolve(j, b);
    azure2::physics::channel_space::kernel_column(f, st, a, ya);
    azure2::physics::channel_space::kernel_column(f, st, b, yb);
    azure2::physics::channel_space::formal_bilinear_adjoint(t, ua.data(), ub.data(), ya.data(), yb.data(), n, nch, sink);
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
    const JGroup *jg = compound()->GetJGroup(j);
    const int n = (int)std::lround(std::sqrt((double)accum.lmK[j - 1].size()));
    matrix_c X1, X2;
    azure2::physics::deferred_level_cotangents(accum.lmK[j - 1].data(), accum.lmN[j - 1].data(), n, X1, X2);
    const std::vector<complex> zero(jg->NumChannels(), complex(0.0, 0.0));
    this->LevelMatrixAdjoint(nullptr, j, X1, 0.0, accum, shiftDeriv, zero.data());
    this->LevelMatrixAdjoint(nullptr, j, X2, 1.0, accum, shiftDeriv, zero.data());
    std::fill(accum.lmK[j - 1].begin(), accum.lmK[j - 1].end(), 0.0);
    std::fill(accum.lmN[j - 1].begin(), accum.lmN[j - 1].end(), 0.0);
    accum.lmUsed[j - 1] = 0;
  }
  for (int j = 1; j <= compound()->NumJGroups(); j++) {
    const JGroup *jg = compound()->GetJGroup(j);
    if (!jg->IsInRMatrix() || j > (int)accum.fe.size()) continue;
    bool any = false;
    for (size_t k = 0; k < accum.fe[j - 1].size() && !any; k++) {
      if (accum.fe[j - 1][k] != 0.0) any = true;
      for (double g : accum.fgamma[j - 1][k])
        if (g != 0.0) any = true;
    }
    if (!any) continue;
    formal_[j - 1] = ChannelSpace::Cached(ChannelMemo(), jg);
    formal_jgroup_[j - 1] = jg;
    ChannelSpace::FormalLevels &f = (*formal_[j - 1]);
    ChannelSpace::Update(jg, state(), j, active_levels_[j - 1], configure(), f);
    const int nch = f.nch;
    azure2::physics::channel_space::FormalResolve resolved;
    azure2::physics::channel_space::resolve_formal_cotangents(f, accum.fe[j - 1], accum.fgamma[j - 1], resolved);
    // The kernel: zero, but for the reference shifts that Update moved into K
    // (constants, like the kernel: F does not depend on them).
    std::vector<complex> kernel(nch, complex(0.0, 0.0));
    for (int c = 0; c < nch; c++)
      if (jg->GetChannel(c + 1)->GetRadType() == 'P')
        kernel[c] = jg->GetChannel(c + 1)->GetBoundaryCondition() - f.lamShift[c];
    this->LevelMatrixAdjoint(nullptr, j, resolved.X1, 0.0, accum, shiftDeriv, kernel.data());
    this->LevelMatrixAdjoint(nullptr, j, resolved.X2, 1.0, accum, shiftDeriv, kernel.data());
    // Gbar = V Gammabar.
    azure2::physics::channel_space::add_level_widths_from_formal(
        f, accum.fgamma[j - 1], LevelSink(accum.e[j - 1], accum.gamma[j - 1], active_levels_[j - 1].data()));
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

/*!
 * The cotangent coefficients of the differential sum of `theDecay` for the
 * T-matrix held now, by pathway and Legendre order: model += cR Re(sum_i
 * Z1Z2_i T1_i conj(T2_i) P_{L_i}) gives T1bar += cR Z1Z2 P_L T2 and
 * T2bar += cR Z1Z2 P_L T1, so each pathway's cotangent is cR sum_L q_L P_L with
 * q_L summed over the interferences once.  `particleExit` keeps only the
 * k-groups of equal channel spins, as the forward sum.  False (nothing built)
 * for an interference type the adjoint does not know.
 */
bool AMatrixFunc::TBarProducts(const azure2::EntrancePathways &P, const azure2::ExitRow &theDecay, bool particleExit) {
  if (tbar_decay_ == &theDecay && tbar_particle_ == particleExit) return true;
  if (!theDecay.has_interferences) return false;
  const azure2::FlatInterference *flat = P.interferences.data() + theDecay.first_interference;
  const int nl = theDecay.num_orders;
  const size_t np = theDecay.num_interference_pathways;
  tbar_q_.assign(np * nl, complex(0.0, 0.0));
  const std::vector<vector_c> *T = tmatrix_.empty() ? nullptr : &tmatrix_[0];
  for (const azure2::FlatInterference *fp = flat; fp != flat + theDecay.num_interferences; fp++) {
    const azure2::FlatInterference &f = *fp;
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
  tbar_decay_ = &theDecay;
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
bool AMatrixFunc::LevelAdjoint(Site *energyPoint, int j, const std::vector<BilinearTerm> &terms,
                               const FreeMask &mask, const vector_matrix_r *shiftDeriv, Add add) {
  if (j < 1 || j > (int)route_.size() || j > (int)mask.e.size()) return false;
  const bool channel = (route_[j - 1] == kRouteChannel);
  if (!channel && !solver_valid_[j - 1]) return false;
  // The standard level matrix in channel space has the cheaper formal adjoint.
  if (channel && (*formal_[j - 1]).identity) return false;
  const JGroup *jg = compound()->GetJGroup(j);
  const int n = (int)active_levels_[j - 1].size(), nch = jg->NumChannels();
  if (n == 0 || n != jg->NumLevels()) return false;
  const bool park = (configure().paramMask & Config::USE_PARK_FORMALISM);
  const bool brune = !park && (configure().paramMask & Config::USE_BRUNE_FORMALISM);
  // Brune takes S' from shiftDeriv, Park from the levels (CompoundSystem::CalcShiftFunctions).
  if (brune && !shiftDeriv) return false;
  // The channels with the shift terms, as LevelMatrixAdjoint takes them.
  std::vector<char> &particle = la_particle_;
  particle.assign(nch, 0);
  for (int c = 0; c < nch; c++) {
    const AChannel *ch = jg->GetChannel(c + 1);
    particle[c] = park ? (ch->GetRadType() == 'P' && compound()->GetPair(ch->GetPairNum())->GetPType() == 0)
                       : (brune && ch->GetRadType() == 'P');
  }
  const ChannelSpace::LevelDerivatives &ld = ChannelSpace::CurrentLevelDerivatives(
      ChannelMemo(), jg, state(), j, park ? 2 : brune ? 1 : 0, particle, mask.e[j - 1], mask.g[j - 1],
      (brune && shiftDeriv) ? &(*shiftDeriv)[j - 1] : nullptr);
  const ChannelSpace::LevelArrays &arrays = ChannelSpace::CurrentLevelArrays(ChannelMemo(), jg, state(), j);
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
  azure2::physics::width_solve_products(arrays.G.data(), u, n, nch, la_gu_);  // [s][d] = g_d . u_s
  azure2::physics::level_parameter_adjoint(ld.entries, terms, slot, u, la_gu_, lam.data(), E, n, nch, la_ku_, add);
  return true;
}

bool AMatrixFunc::LevelBasis(Site *energyPoint, int j, const BilinearTerm &term, const FreeMask &mask,
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
