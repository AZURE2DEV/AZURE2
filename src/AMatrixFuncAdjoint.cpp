// AMatrixFunc, T -> level parameters: the adjoint of the T matrix (TBarAdjoint) through the level matrix
// (LevelMatrixAdjoint, LevelMatrixABarAdjoint) and the channel-space bilinear forms (CollectBilinearTerms,
// BilinearAdjoint, FormalBilinearAdjoint, ResolveFormalCotangents), and the derivatives by free level parameter
// (LevelAdjoint, LevelBasis).  The cotangent convention is stated above PointTBar (AMatrixFuncTBar.cpp); the
// forward solve these mirror is AMatrixFunc.cpp.
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
  // The dense path takes the level index as the active one (PointTBar refuses anything else).
  if (!all_levels_active_) throw std::logic_error("TBarAdjoint: a level outside the R matrix");
  const double inEnergy = this->InEnergy(point);
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
        complex base = azure2::physics::external_capture_element(gamma_fc, sqrtNF, conv, ecAmp);

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

          // Forward umatrix_cc on the level set of the forward pass (the IGNORE_ZERO_WIDTHS mask), and its adjoint:
          // T_ec = base * umatrix_cc -> ubar = conj(base) * eb.
          const complex factor = I2 * sqrtPcc;
          const int maskChannel = (internalChannel && ignoreZero) ? internalChannel - 1 : -1;
          const matrix_c &Acc = this->GetAMatrix(jcc);
          const double *gcc = state().Gammas(jcc);
          const unsigned char *zero = state().ZeroInputWidths(jcc);
          const int nchcc = jgcc->NumChannels();
          const complex umatrix_cc =
              azure2::physics::masked_amplitude(factor, gcc, nchcc, ccChN - 1, ccChpN - 1, Acc, zero, maskChannel);
          effAmp = ecAmp * umatrix_cc;
          complex uBarcc = conj(base) * eb;
          ensureABar(jcc, ccLevels);
          azure2::physics::dense_bilinear_adjoint(factor, gcc, nchcc, ccChN - 1, ccChpN - 1, Acc, zero, maskChannel,
                                                  uBarcc, aBar[jcc - 1], LevelSink(accum, jcc));
        }

        // Explicit gamma_fc:  d T_ec / d gamma_fc = sqrtNF * conv * effAmp.
        accum.AddGamma(jf, lf, fc, azure2::physics::external_capture_width_cotangent(eb, sqrtNF, conv, effAmp));

        // sqrtNF couples to every NF channel d of the final level:
        //   nFSum = 1 + sum_d A_d gamma_d^2,  A_d = 2 chRad redMass uconv/hbarc^2 NF_d (physics/ExternalCapture.h),
        //   sqrtNF = nFSum^{-1/2},  d sqrtNF / d gamma_d = -A_d gamma_d sqrtNF^3.
        // Park: sqrtNF = (J + sum_d A_d gamma_d^2)^{-1/2} with J = 1 - sum_d gamma_d^2 S'_d
        // over the particle channels (CompoundSystem::CalcShiftFunctions), so every particle
        // channel contributes  d sqrtNF / d gamma_d = (S'_d - A_d) gamma_d sqrtNF^3.
        complex pre = azure2::physics::external_capture_norm_derivative(gamma_fc, conv, effAmp);  // d T_ec / d sqrtNF
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
          accum.AddGamma(jf, lf, d, azure2::physics::bound_norm_width_cotangent(eb, pre, Ad, gd, sqrtNF3));
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

      const azure2::physics::ChannelPhases entrance{point->GetExpCoulombPhase(jNum, chNum),
                                                    point->GetExpHardSpherePhase(jNum, chNum)};
      const azure2::physics::ChannelPhases exit{point->GetExpCoulombPhase(jNum, chpNum),
                                                point->GetExpHardSpherePhase(jNum, chpNum)};
      complex sqrtPenEn = point->GetSqrtPenetrability(jNum, chNum);
      complex sqrtPenEx = point->GetSqrtPenetrability(jNum, chpNum);

      // d T / d umatrix = -uphase  ->  ubar = conj(-uphase) * tBar.
      const complex uBar = azure2::physics::amplitude_cotangent(entrance, exit, tb);

      ensureABar(jNum, jg->NumLevels());
      // A cotangent Abar += conj(w) ubar with w = pen2 gamma_en gamma_ex, and the widths as they stand in umatrix.
      const complex pen2 = azure2::physics::penetrability_factor(sqrtPenEn, sqrtPenEx);
      azure2::physics::dense_bilinear_adjoint(pen2, state().Gammas(jNum), jg->NumChannels(), chNum - 1, chpNum - 1,
                                              this->GetAMatrix(jNum), nullptr, -1, uBar, aBar[jNum - 1],
                                              LevelSink(accum, jNum));
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
 * T -> bilinear forms.  Each internal pathway T(k,m) of J-group j is
 *   T = [tphase] - uphase (1 + 2i sqrtP_en sqrtP_ex b),   b = g_en^T A g_ex,
 * so its cotangent tbar gives the bilinear form the cotangent
 *   beta = conj(2i sqrtP_en sqrtP_ex) conj(-uphase) tbar.
 * Pathways of one J-group sharing a channel pair add their betas.
 */
bool AMatrixFunc::CollectBilinearTerms(Site *energyPoint, const TMatrixCotangent &bar,
                                       std::vector<std::vector<BilinearTerm>> &terms) {
  Site *point = energyPoint;
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
      const azure2::physics::ChannelPhases entrance{point->GetExpCoulombPhase(jNum, chNum),
                                                    point->GetExpHardSpherePhase(jNum, chNum)};
      const azure2::physics::ChannelPhases exit{point->GetExpCoulombPhase(jNum, chpNum),
                                                point->GetExpHardSpherePhase(jNum, chpNum)};
      const complex uBar = azure2::physics::amplitude_cotangent(entrance, exit, tb);
      const complex beta = azure2::physics::bilinear_cotangent(
          azure2::physics::penetrability_factor(point->GetSqrtPenetrability(jNum, chNum),
                                                point->GetSqrtPenetrability(jNum, chpNum)),
          uBar);
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
