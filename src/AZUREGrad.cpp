#include "AZUREGrad.h"
#include "Profile.h"
#include <atomic>
#include <chrono>
#include <iostream>
#include <cstdlib>
#include "CNuc.h"
#include "EData.h"
#include "JGroup.h"
#include "ALevel.h"
#include "AChannel.h"
#include "PPair.h"
#include "ESegment.h"
#include "EPoint.h"
#include "TargetEffect.h"
#include "AMatrixFunc.h"
#include "Config.h"
#include "CoulFunc.h"
#include "ShftFunc.h"
#include <cmath>
#include <algorithm>
#include <string>
#include <unordered_map>
#ifdef _OPENMP
#include <omp.h>
#endif

/*!
 * Precompute dS_ch/dE_level for every (jGroup, level, channel).  Mirrors the
 * bound/unbound branch of CNuc::CalcShiftFunctions.  In the Park formalism
 * CalcShiftFunctions (which the callers run first) already holds the same
 * derivative on each level, from the same function at the same energy: it is
 * read from there instead of being differentiated again.
 */
vector_matrix_r BuildShiftDerivTable(CNuc *compound, const Config &configure) {
  const bool park = (configure.paramMask & Config::USE_PARK_FORMALISM);
  vector_matrix_r tab(compound->NumJGroups());
  for (int j = 1; j <= compound->NumJGroups(); j++) {
    JGroup *jg = compound->GetJGroup(j);
    tab[j - 1].assign(jg->NumLevels(), std::vector<double>());
    if (!jg->IsInRMatrix()) continue;
    for (int la = 1; la <= jg->NumLevels(); la++) {
      tab[j - 1][la - 1].assign(jg->NumChannels(), 0.0);
      ALevel *level = jg->GetLevel(la);
      if (!level->IsInRMatrix()) continue;
      for (int ch = 1; ch <= jg->NumChannels(); ch++) {
        AChannel *channel = jg->GetChannel(ch);
        PPair *pair = compound->GetPair(channel->GetPairNum());
        // The Brune boundary/shift block applies only to particle 'P' channels.
        if (pair->GetPType() != 0 || channel->GetRadType() != 'P') continue;
        if (park && level->ShiftFunctionsValidAt(level->GetFitE())) {
          tab[j - 1][la - 1][ch - 1] = level->GetShiftDerivative(ch);
          continue;
        }
        int l = channel->GetL();
        double levelEnergy = level->GetFitE();
        double resonanceEnergy = levelEnergy - (pair->GetSepE() + pair->GetExE());
        double dS;
        if (resonanceEnergy < 0.0) {
          ShftFunc theShiftFunction(pair);
          dS = theShiftFunction.EnergyDerivative(l, levelEnergy);
        } else {
          CoulFunc theCoulombFunction(pair,
                                      !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC));
          double radius = pair->GetChRad();
          dS = theCoulombFunction.PEShift_dE(l, radius, resonanceEnergy);
        }
        tab[j - 1][la - 1][ch - 1] = dS;
      }
    }
  }
  return tab;
}

/*!
 * d/dx of sum_lambda (J_lambda/s)^2 over J < 0, with J = 1 - sum_d g_d^2 S'_d:
 *   d/d g_d  = 2 J/s^2 (-2 g_d S'_d),    d/d E = 2 J/s^2 (-sum_d g_d^2 S''_d).
 */
void AddParkPenaltyGradient(CNuc *compound, const Config &configure, GradAccum &accum) {
  if (!(configure.paramMask & Config::USE_PARK_FORMALISM)) return;
  const double s2 = CNuc::kParkNormScale * CNuc::kParkNormScale;
  for (int j = 1; j <= compound->NumJGroups(); j++) {
    JGroup *jg = compound->GetJGroup(j);
    if (!jg->IsInRMatrix()) continue;
    for (int la = 1; la <= jg->NumLevels(); la++) {
      ALevel *level = jg->GetLevel(la);
      if (!level->IsInRMatrix() || level->GetParkNorm() >= 0.0) continue;
      const double factor = 2.0 * level->GetParkNorm() / s2;
      double dE = 0.0;
      for (int ch = 1; ch <= jg->NumChannels(); ch++) {
        AChannel *channel = jg->GetChannel(ch);
        if (channel->GetRadType() != 'P' || compound->GetPair(channel->GetPairNum())->GetPType() != 0) continue;
        double g = level->GetFitGamma(ch);
        accum.AddGamma(j, la, ch, factor * (-2.0 * g * level->GetShiftDerivative(ch)));
        dE -= g * g * level->GetShiftSecondDerivative(ch);
      }
      accum.AddE(j, la, factor * dE);
    }
  }
}

/*!
 * Gradient of the abundance penalty (CNuc::AbundancePenalty): the prior
 * ((w - w0)/sigma)^2 of every nucleus with a sigma, and the sum-to-one term
 * ((sum w) - 1)^2 / s^2, which couples every abundance through the sum.
 */
void AddAbundancePenaltyGradient(CNuc *compound, const Config &configure, GradAccum &accum) {
  if (!compound->HasAbundanceBlock()) return;
  for (int k = 1; k <= compound->NumAbundanceParams(); k++) {
    int sample, nucleus;
    compound->AbundanceOf(k, sample, nucleus);
    const double s = compound->SampleSumSigma(sample, configure.abundanceSumSigma);
    double g = (s > 0.0) ? 2.0 * (compound->AbundanceSum(sample) - 1.0) / (s * s) : 0.0;
    const SampleAbundance &e = compound->GetSample(sample)->entries[nucleus - 1];
    if (e.sigma > 0.0) g += 2.0 * (e.abundance - e.nominal) / (e.sigma * e.sigma);
    accum.AddAbundance(k, g);
  }
}

/*!
 * Allocates and zeroes the physics-coordinate gradient accumulators.
 */
void GradAccum::Init(CNuc *compound) {
  const int nJ = compound->NumJGroups();
  e.assign(nJ, {});
  gamma.assign(nJ, {});
  for (int j = 1; j <= nJ; j++) {
    JGroup *jg = compound->GetJGroup(j);
    const int nLevels = jg->NumLevels();
    const int nCh = jg->NumChannels();
    e[j - 1].assign(nLevels, 0.0);
    gamma[j - 1].assign(nLevels, std::vector<double>(nCh, 0.0));
  }
  fe = e;
  fgamma = gamma;
  formalUsed = false;
  lmK.assign(nJ, {});
  lmN.assign(nJ, {});
  lmUsed.assign(nJ, 0);
  abundance.assign(compound->NumAbundanceParams(), 0.0);
  screening.clear();
}

/*!
 * Resets all accumulators to zero without reallocating.
 */
void GradAccum::Zero() {
  for (auto &row : e) std::fill(row.begin(), row.end(), 0.0);
  for (auto &lev : gamma)
    for (auto &row : lev) std::fill(row.begin(), row.end(), 0.0);
  std::fill(abundance.begin(), abundance.end(), 0.0);
  std::fill(screening.begin(), screening.end(), 0.0);
  if (formalUsed) {
    for (auto &row : fe) std::fill(row.begin(), row.end(), 0.0);
    for (auto &lev : fgamma)
      for (auto &row : lev) std::fill(row.begin(), row.end(), 0.0);
    formalUsed = false;
  }
  for (size_t j = 0; j < lmUsed.size(); j++)
    if (lmUsed[j]) {
      std::fill(lmK[j].begin(), lmK[j].end(), 0.0);
      std::fill(lmN[j].begin(), lmN[j].end(), 0.0);
      lmUsed[j] = 0;
    }
}

/*!
 * Adds another (same-dimension) accumulator into this one, term by term.
 */
void GradAccum::Add(const GradAccum &o) {
  for (int j = 0; j < (int)e.size(); j++) {
    for (int la = 0; la < (int)e[j].size(); la++) {
      e[j][la] += o.e[j][la];
      for (int ch = 0; ch < (int)gamma[j][la].size(); ch++)
        gamma[j][la][ch] += o.gamma[j][la][ch];
    }
  }
  for (int k = 0; k < (int)abundance.size() && k < (int)o.abundance.size(); k++) abundance[k] += o.abundance[k];
  if (o.screening.size() > screening.size()) screening.resize(o.screening.size(), 0.0);
  for (int k = 0; k < (int)o.screening.size(); k++) screening[k] += o.screening[k];
  for (size_t j = 0; j < o.lmUsed.size() && j < lmUsed.size(); j++)
    if (o.lmUsed[j]) {
      if (lmK[j].size() != o.lmK[j].size()) {
        lmK[j].assign(o.lmK[j].size(), 0.0);
        lmN[j].assign(o.lmN[j].size(), 0.0);
      }
      for (size_t k = 0; k < lmK[j].size(); k++) {
        lmK[j][k] += o.lmK[j][k];
        lmN[j][k] += o.lmN[j][k];
      }
      lmUsed[j] = 1;
      formalUsed = true;
    }
  if (o.formalUsed) {
    for (int j = 0; j < (int)fe.size() && j < (int)o.fe.size(); j++)
      for (int k = 0; k < (int)fe[j].size() && k < (int)o.fe[j].size(); k++) {
        fe[j][k] += o.fe[j][k];
        for (int ch = 0; ch < (int)fgamma[j][k].size(); ch++) fgamma[j][k][ch] += o.fgamma[j][k][ch];
      }
    formalUsed = true;
  }
}

/*!
 * Scatters the physics-coordinate gradients into the flat parameter vector.
 */
void GradAccum::Scatter(const ParamIndexMap &map, vector_r &gradFull) const {
  for (int j = 0; j < (int)e.size(); j++) {
    for (int la = 0; la < (int)e[j].size(); la++) {
      int ei = map.EIndexFlat(j + 1, la + 1);
      if (ei >= 0 && ei < (int)gradFull.size()) gradFull[ei] += e[j][la];
      for (int ch = 0; ch < (int)gamma[j][la].size(); ch++) {
        int gi = map.GammaIndexFlat(j + 1, la + 1, ch + 1);
        if (gi >= 0 && gi < (int)gradFull.size()) gradFull[gi] += gamma[j][la][ch];
      }
    }
  }
  for (int k = 0; k < (int)abundance.size(); k++) {
    int ai = map.AbundanceIndex(k + 1);
    if (ai >= 0 && ai < (int)gradFull.size()) gradFull[ai] += abundance[k];
  }
  for (int k = 0; k < (int)screening.size(); k++) {
    int si = map.ScreeningIndex(k + 1);
    if (si >= 0 && si < (int)gradFull.size()) gradFull[si] += screening[k];
  }
}

/*!
 * Returns the full-vector index of E_{jGroup,level}, or -1 if absent.
 */
int ParamIndexMap::EIndex(int jGroup, int level) const {
  auto it = eIndex_.find({jGroup, level});
  return (it == eIndex_.end()) ? -1 : it->second;
}

/*!
 * Returns the full-vector index of gamma_{jGroup,level,channel}, or -1 if absent.
 */
int ParamIndexMap::GammaIndex(int jGroup, int level, int channel) const {
  auto it = gammaIndex_.find({jGroup, level, channel});
  return (it == gammaIndex_.end()) ? -1 : it->second;
}

/*!
 * Returns the full-vector index of the norm parameter for a segment, or -1.
 */
int ParamIndexMap::NormIndex(int segment) const {
  auto it = normIndex_.find(segment);
  return (it == normIndex_.end()) ? -1 : it->second;
}

/*!
 * Returns the full-vector index of the energy-shift parameter for a segment.
 */
int ParamIndexMap::EnergyShiftIndex(int segment) const {
  auto it = shiftIndex_.find(segment);
  return (it == shiftIndex_.end()) ? -1 : it->second;
}

/*!
 * Returns the full-vector index of the sqrt(E) energy-shift coefficient for a segment.
 */
int ParamIndexMap::EnergyShiftSqrtIndex(int segment) const {
  auto it = sqrtShiftIndex_.find(segment);
  return (it == sqrtShiftIndex_.end()) ? -1 : it->second;
}

int ParamIndexMap::AbundanceIndex(int ordinal) const {
  auto it = abundanceIndex_.find(ordinal);
  return (it == abundanceIndex_.end()) ? -1 : it->second;
}

/*!
 * Returns the full-vector index of parameter aj of a segment's shift function, or -1.
 */
int ParamIndexMap::EnergyShiftFnIndex(int segment, int j) const {
  auto it = fnShiftIndex_.find({segment, j});
  return (it == fnShiftIndex_.end()) ? -1 : it->second;
}

/*!
 * Builds the parameter-index map by mirroring, in the exact same order, the
 * assignment loops in:
 *   - CNuc::FillCompoundFromParams      (level energies + gammas)
 *   - EData::FillNormsFromParams        (norms, only for IsVaryNorm segments)
 *   - EData::FillEnergyShiftsFromParams (one energy shift per segment, then
 *                                        one sqrt(E) coefficient per segment,
 *                                        then the shift-function parameters of
 *                                        each segment that declares them)
 *
 * Keeping this in lock-step with those routines is what guarantees the returned
 * gradient vector is a permutation-correct match for the sampler's parameter
 * ordering (see PLAN.md, "Parameter-vector order must exactly match").
 */
ParamIndexMap BuildParamIndexMap(CNuc *compound, EData *data,
                                 const std::vector<bool> &fixed) {
  ParamIndexMap map;
  int i = 0;

  // --- R-matrix parameters: mirror FillCompoundFromParams exactly. ---
  map.eFlat_.assign(compound->NumJGroups(), {});
  map.gammaFlat_.assign(compound->NumJGroups(), {});
  for (int j = 1; j <= compound->NumJGroups(); j++) {
    JGroup *jGroup = compound->GetJGroup(j);
    map.eFlat_[j - 1].assign(jGroup->NumLevels(), -1);
    map.gammaFlat_[j - 1].assign(jGroup->NumLevels(), std::vector<int>(jGroup->NumChannels(), -1));
    for (int la = 1; la <= jGroup->NumLevels(); la++) {
      // Energy parameter.
      map.eFlat_[j - 1][la - 1] = i;
      map.eIndex_[{j, la}] = i;
      map.desc_.push_back(ParamDesc{ParamKind::LevelEnergy, j, la, -1, -1, false});
      i++;
      // One gamma per channel.
      for (int ch = 1; ch <= jGroup->NumChannels(); ch++) {
        map.gammaFlat_[j - 1][la - 1][ch - 1] = i;
        map.gammaIndex_[{j, la, ch}] = i;
        map.desc_.push_back(ParamDesc{ParamKind::Gamma, j, la, ch, -1, false});
        i++;
      }
    }
  }

  // --- Isotopic abundances: one per (sample, nucleus), only with an <abundances> block. ---
  for (int k = 1; k <= compound->NumAbundanceParams(); k++) {
    map.abundanceIndex_[k] = i;
    ParamDesc d{ParamKind::Abundance, -1, -1, -1, -1, false};
    compound->AbundanceOf(k, d.sample, d.nucleus);
    map.desc_.push_back(d);
    i++;
  }

  // --- Norm parameters: one per varying-norm segment. ---
  map.normOffset_ = i;
  if (data) {
    for (int s = 1; s <= data->NumSegments(); s++) {
      ESegment *segment = data->GetSegment(s);
      if (segment && segment->IsVaryNorm()) {
        map.normIndex_[s] = i;
        map.desc_.push_back(ParamDesc{ParamKind::Norm, -1, -1, -1, s, false});
        i++;
      }
    }
  }

  // --- Energy-shift parameters: one per segment (all segments). ---
  map.energyShiftOffset_ = i;
  if (data) {
    for (int s = 1; s <= data->NumSegments(); s++) {
      ESegment *segment = data->GetSegment(s);
      if (segment) {
        map.shiftIndex_[s] = i;
        map.desc_.push_back(ParamDesc{ParamKind::EnergyShift, -1, -1, -1, s, false});
        i++;
      }
    }
  }

  // --- sqrt(E) energy-shift coefficients: one per segment (all segments). ---
  map.energyShiftSqrtOffset_ = i;
  if (data) {
    for (int s = 1; s <= data->NumSegments(); s++) {
      ESegment *segment = data->GetSegment(s);
      if (segment) {
        map.sqrtShiftIndex_[s] = i;
        map.desc_.push_back(ParamDesc{ParamKind::EnergyShiftSqrt, -1, -1, -1, s, false});
        i++;
      }
    }
  }

  // --- Shift-function parameters: those of each segment that declares one (EData::FillMnParams). ---
  map.energyShiftFnOffset_ = i;
  if (data) {
    for (int s = 1; s <= data->NumSegments(); s++) {
      ESegment *segment = data->GetSegment(s);
      if (!segment || !segment->OwnsShiftFn() || !segment->HasShiftFn()) continue;
      for (int j = 0; j < segment->NumShiftFnParams(); j++) {
        map.fnShiftIndex_[{s, j}] = i;
        ParamDesc d{ParamKind::EnergyShiftFn, -1, -1, -1, s, false};
        d.fnTerm = j;
        map.desc_.push_back(d);
        i++;
      }
    }
  }

  // --- Electron-screening potentials: one per <targetInt> line with a screening block (EData::FillMnParams). ---
  if (data) {
    for (int k = 1; k <= data->NumScreeningParams(); k++) {
      map.screeningIndex_[k] = i;
      ParamDesc d{ParamKind::Screening, -1, -1, -1, -1, false};
      d.screening = k;
      map.desc_.push_back(d);
      i++;
    }
  }

  // --- Apply the fixed mask and build packed<->full lookups. ---
  const int nFull = (int)map.desc_.size();
  map.fullToPacked_.assign(nFull, -1);
  map.packedToFull_.clear();
  for (int f = 0; f < nFull; f++) {
    bool isFixed = (f < (int)fixed.size()) ? fixed[f] : false;
    map.desc_[f].fixed = isFixed;
    if (!isFixed) {
      map.fullToPacked_[f] = (int)map.packedToFull_.size();
      map.packedToFull_.push_back(f);
    }
  }
  map.numPacked_ = (int)map.packedToFull_.size();

  return map;
}

// ===========================================================================
// Shared adjoint orchestration (used by AZUREAPI and AZURECalc).
// ===========================================================================

namespace {

// Forward + adjoint for one target-effect (sub-point integrated) point.  The
// integrator IntegrateTargetEffect is a fixed linear functional of the values
// in the sub-points' cross-section slots, Y = sum_i c_i s_i, so the weights
// c_i = dY/ds_i come from central differences on that cheap combiner (exact up
// to rounding); each sub-point is then back-propagated with PointAdjoint.
//
// A cross section (or E1/E2 component) is Y itself.  A target-averaged
// analyzing power is the ratio of two such functionals with the SAME weights
// (EPoint::IntegrateTargetEffectForObservable runs the one integrator on s_i
// and on s_i A_i, s_i being the sub-point's differential cross section at the
// point's angle):
//     <A> = N / D,   N = sum_i c_i s_i A_i,   D = sum_i c_i s_i,
//     d<A> = sum_i c_i [ s_i dA_i + (A_i - <A>) ds_i ] / D,
// so each sub-point contributes an A_y cotangent fitBar c_i s_i / D and a
// differential-cross-section cotangent fitBar c_i (A_i - <A>) / D.
// AZURE_GRAD_TIMING: thread-summed time in the target-effect combiner weights
// (central differences on IntegrateTargetEffect) and in the sub-point adjoints.
std::atomic<long long> g_teCombinerNs{0}, g_teAdjointNs{0};

bool GradTargetEffectAdjoint(EPoint *point, double fitBar, GradAccum &accum,
                             const vector_matrix_r *shiftDeriv, int xsComponent,
                             CNuc *compound, const Config &config) {
  int nsp = point->NumSubPoints();
  if (nsp <= 0) return true;

  // EPoint::Calculate blends the integrated value with the bare one at the
  // edge of an energy-restricted effect (blend < 1) and, with an automatic
  // tolerance, may replace it by the bare value altogether.  Neither is the
  // pure integrator differentiated below, so leave such points to finite
  // differences.
  {
    TargetEffect *te = point->GetParentData()->GetTargetEffect(point->GetTargetEffectNum());
    if (point->GetTargetBlendWeight() < 1.0) return false;
    if (te && te->GetAutoTolerance() > 0.) return false;
  }

  const bool isAy = point->IsAnalyzingPower();
  if (isAy) {
    // Supported: particle-exit A_y, full observable, no point mapping in
    // effect (a mapped point's sub-points would be written by another point's
    // calculation).  GradOnePoint holds a ScopedPointMappingOff, so there every
    // point, mapped or not, carries and calculates its own sub-points.  Capture
    // A_y keeps its refusal: the attenuation folded into P_k has no
    // counterpart for the P_k^1 of the numerator.
    if (xsComponent != 0) return false;
    if (!ScopedPointMappingOff::Active() && (point->IsMapped() || point->NumLocalMappedPoints() > 0)) return false;
    const int exitPairNum = compound->GetPairNumFromKey(point->GetExitKey());
    if (compound->GetPair(exitPairNum)->GetPType() != 0) return false;
  }

  std::vector<double> sigma(nsp);
  for (int i = 1; i <= nsp; i++) sigma[i - 1] = point->GetSubPoint(i)->GetFitCrossSection();
  double model = point->GetFitCrossSection();

  const auto tc0 = std::chrono::steady_clock::now();
  // Exact weights in one pass of the integrator on a linear form
  // (EPoint::TargetEffectWeights).  The central differences used before cost
  // 2 nsp integrator calls, each O(nsp) or O(nsp^2) (straggling, convolution
  // + target integration): O(nsp^3) per point, which made one Jacobian of the
  // 11B+a base grid (CASPAR thick targets, ~2200 sub-points per point) take
  // longer than 30 min.  They remain as the fallback and as a check
  // (AZURE_TE_WEIGHTS_CHECK=N compares the first N points).
  auto fdWeights = [&](std::vector<double> &w) {
    w.assign(nsp, 0.0);
    for (int i = 1; i <= nsp; i++) {
      double s0 = sigma[i - 1];
      double h = 1.0e-3 * (std::fabs(s0) + 1.0);
      point->GetSubPoint(i)->SetFitCrossSection(s0 + h);
      point->IntegrateTargetEffect(config);
      double yp = point->GetFitCrossSection();
      point->GetSubPoint(i)->SetFitCrossSection(s0 - h);
      point->IntegrateTargetEffect(config);
      double ym = point->GetFitCrossSection();
      point->GetSubPoint(i)->SetFitCrossSection(s0);
      w[i - 1] = (yp - ym) / (2.0 * h);
    }
  };
  std::vector<double> c;
  if (std::getenv("AZURE_TE_WEIGHTS_FD") || !point->TargetEffectWeights(config, c)) {
    fdWeights(c);
  } else if (const char *wc = std::getenv("AZURE_TE_WEIGHTS_CHECK")) {
    static std::atomic<int> nChecked{0};
    if (nChecked++ < std::atoi(wc)) {
      std::vector<double> cfd;
      fdWeights(cfd);
      double cmax = 0.0, dmax = 0.0, yW = 0.0;
      for (int i = 0; i < nsp; i++) {
        cmax = std::max(cmax, std::fabs(cfd[i]));
        dmax = std::max(dmax, std::fabs(c[i] - cfd[i]));
        yW += c[i] * sigma[i];
      }
      point->IntegrateTargetEffect(config);
      const double yI = point->GetFitCrossSection();
#pragma omp critical
      std::cerr << "[te-weights] entrance " << point->GetEntranceKey() << " exit " << point->GetExitKey() << " E " << point->GetCMEnergy() << " nsp " << nsp
                << ": max|c - c_fd|/max|c_fd| " << (cmax > 0 ? dmax / cmax : dmax) << ", sum c*s " << yW
                << " vs integrator " << yI << std::endl;
    }
  }
  // For A_y: the denominator D and the average <A>, from the same integrator
  // calls the forward pass makes.
  double den = 0.0, avg = 0.0;
  if (isAy) {
    point->IntegrateTargetEffect(config);
    den = point->GetFitCrossSection();
    for (int i = 1; i <= nsp; i++)
      point->GetSubPoint(i)->SetFitCrossSection(sigma[i - 1] * point->GetSubPoint(i)->GetAnalyzingPower());
    point->IntegrateTargetEffect(config);
    const double num = point->GetFitCrossSection();
    for (int i = 1; i <= nsp; i++) point->GetSubPoint(i)->SetFitCrossSection(sigma[i - 1]);
    // The forward pass returns 0 for D = 0, a constant: no gradient.
    if (!(std::fabs(den) > 0.0)) {
      point->SetFitCrossSection(model);
      return true;
    }
    avg = num / den;
  }
  point->SetFitCrossSection(model);
  const auto tc1 = std::chrono::steady_clock::now();
  g_teCombinerNs += std::chrono::duration_cast<std::chrono::nanoseconds>(tc1 - tc0).count();

  bool ok = true;
  for (int i = 1; i <= nsp && ok; i++) {
    if (c[i - 1] == 0.0) continue;
    EPoint *sp = point->GetSubPoint(i);
    AMatrixFunc amf(compound, config);
    amf.ClearMatrices();
    amf.FillMatrices(sp);
    amf.InvertMatrices();
    amf.CalculateTMatrix(sp);
    if (!isAy) {
      if (!amf.PointAdjoint(sp, fitBar * c[i - 1], accum, shiftDeriv, xsComponent)) ok = false;
      continue;
    }
    const double ayBar = fitBar * c[i - 1] * sigma[i - 1] / den;
    const double xsBar = fitBar * c[i - 1] * (sp->GetAnalyzingPower() - avg) / den;
    if (ayBar != 0.0 && !amf.PointAdjoint(sp, ayBar, accum, shiftDeriv, 0, false)) ok = false;
    if (ok && xsBar != 0.0 && !amf.PointAdjoint(sp, xsBar, accum, shiftDeriv, 0, true)) ok = false;
  }
  g_teAdjointNs += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - tc1).count();
  return ok;
}

// Forward + adjoint for one data point (target-effect or direct A-matrix).
bool GradAdjointPoint(EPoint *point, int xsComponent, double fitBar, GradAccum &accum,
                      const vector_matrix_r *shiftDeriv, CNuc *compound, const Config &config) {
  TargetEffect *te = point->IsTargetEffect()
      ? point->GetParentData()->GetTargetEffect(point->GetTargetEffectNum())
      : nullptr;
  bool subpointTE = te && te->IsSubPointEffect();
  if (subpointTE)
    return GradTargetEffectAdjoint(point, fitBar, accum, shiftDeriv, xsComponent, compound, config);

  AMatrixFunc amf(compound, config);
  amf.ClearMatrices();
  amf.FillMatrices(point);
  amf.InvertMatrices();
  amf.CalculateTMatrix(point);
  return amf.PointAdjoint(point, fitBar, accum, shiftDeriv, xsComponent);
}

// Forward + adjoint for one (segment, point): returns the combined model in
// `outModel` and accumulates fitBar * d(model)/d(E,gamma) into `accum`, where
// fitBar = fitBarFn(segment, i, pointIdx, model).  Handles plain points (single
// forward, reused for model + adjoint), component (SUM/RATIO) and target-effect
// points.  Returns false on an unsupported configuration.
bool GradOnePoint(ESegment *segment, EData *data, int i, int pointIdx,
                  CNuc *compound, const Config &config,
                  const vector_matrix_r *shiftDeriv, const FitBarFn &fitBarFn,
                  GradAccum &accum, double &outModel) {
  // The points of a segment are processed in parallel here, so each one is
  // calculated on its own, without point mapping (see ScopedPointMappingOff).
  ScopedPointMappingOff noPointMapping;
  EPoint *point = segment->GetPoint(pointIdx + 1);
  outModel = 0.0;
  if (!point) return true;

  int baseComp = segment->GetCrossSectionComponent();
  bool hasComp = segment->HasComponents();
  // A natural-target segment weights each term by its nucleus's abundance;
  // d(model)/d(abundance_k) is the unweighted sum of that nucleus's terms.
  const bool natural = segment->IsNatural();
  const int sample = segment->GetNaturalSample();
  auto abundanceOf = [&](int entranceKey) {
    return compound->GetPair(compound->GetPairNumFromKey(entranceKey))->GetNucleus();
  };

  TargetEffect *te = point->IsTargetEffect()
      ? point->GetParentData()->GetTargetEffect(point->GetTargetEffectNum())
      : nullptr;
  bool subpointTE = te && te->IsSubPointEffect();

  // Fast path: plain point, single forward reused for model and adjoint.
  if (!hasComp && !subpointTE) {
    AMatrixFunc amf(compound, config);
    amf.ClearMatrices();
    amf.FillMatrices(point);
    amf.InvertMatrices();
    amf.CalculateTMatrix(point);
    amf.GenMatrixFunc::CalculateCrossSection(point);
    double model = (baseComp == 1) ? point->GetFitE1CrossSection()
        : (baseComp == 2)          ? point->GetFitE2CrossSection()
                                   : point->GetFitCrossSection();
    double weight = 1.0;
    if (natural) {
      const int nucleus = abundanceOf(segment->GetEntranceKey());
      weight = compound->GetAbundance(sample, nucleus);
      outModel = weight * model;
      double fitBar = fitBarFn(segment, i, pointIdx, outModel);
      if (fitBar == 0.0) return true;
      accum.AddAbundance(compound->AbundanceOrdinal(sample, nucleus), fitBar * model);
      return amf.PointAdjoint(point, fitBar * weight, accum, shiftDeriv, baseComp);
    }
    outModel = model;
    double fitBar = fitBarFn(segment, i, pointIdx, model);
    if (fitBar != 0.0)
      return amf.PointAdjoint(point, fitBar, accum, shiftDeriv, baseComp);
    return true;
  }

  // A SUM or RATIO segment whose components share the base point's energy (EData::MapComponentSegments): one
  // level-matrix solution for every term, each term's T-matrix and cotangent in turn, and one back-propagation of
  // the merged cotangents (TBarAdjointMany).  Before 2026-10-24 each term was calculated twice on its own (forward,
  // then forward + adjoint).  AZURE_NO_COMPONENT_MAPPING=1 keeps the per-term route below.
  static const bool componentMappingOff = [] {
    const char *e = std::getenv("AZURE_NO_COMPONENT_MAPPING");
    return e && *e && std::string(e) != "0";
  }();
  if (hasComp && !subpointTE && !componentMappingOff) {
    std::vector<EPoint *> pts{point};
    std::vector<int> xsComp{baseComp};
    std::vector<ESegment *> termSeg{segment};
    bool shared = true;
    for (ESegment *comp : segment->GetComponentSegments()) {
      if (!comp || pointIdx >= comp->NumPoints()) continue;
      EPoint *cp = comp->GetPoint(pointIdx + 1);
      if (!cp) continue;
      if (cp->GetEntranceKey() != point->GetEntranceKey() || cp->GetCMEnergy() != point->GetCMEnergy() ||
          cp->NumSubPoints() != 0)
        shared = false;
      pts.push_back(cp);
      xsComp.push_back(comp->GetCrossSectionComponent());
      termSeg.push_back(comp);
    }
    if (segment->GetOperationType() != SUM && pts.size() < 2) shared = false;
    if (shared) {
      AMatrixFunc amf(compound, config);
      amf.ClearMatrices();
      amf.FillMatrices(point);
      amf.InvertMatrices();
      int tExit = point->GetExitKey();
      amf.CalculateTMatrix(point);
      auto useExitOf = [&](EPoint *p) {
        if (p->GetExitKey() == tExit) return;
        amf.RecalculateTMatrix(p);
        tExit = p->GetExitKey();
      };
      for (EPoint *p : pts) {
        useExitOf(p);
        amf.GenMatrixFunc::CalculateCrossSection(p);
      }
      // The combined model from the values the points now hold, exactly as the forward calculation forms it.
      double model = segment->CalculateTheoreticalCrossSection(pointIdx, compound, config, data, true);
      outModel = model;
      double fitBar = fitBarFn(segment, i, pointIdx, model);
      if (fitBar == 0.0) return true;
      auto valueOf = [](EPoint *p, int c) {
        return (c == 1) ? p->GetFitE1CrossSection() : (c == 2) ? p->GetFitE2CrossSection() : p->GetFitCrossSection();
      };
      // d(model)/d(term q)
      std::vector<double> coef(pts.size(), 0.0);
      if (segment->GetOperationType() == SUM) {
        for (size_t q = 0; q < pts.size(); q++) {
          const double v = valueOf(pts[q], xsComp[q]);
          if (q > 0 && (std::isnan(v) || std::isinf(v))) continue;
          double w = 1.0;
          if (natural) {
            const int nucleus = abundanceOf(termSeg[q]->GetEntranceKey());
            w = compound->GetAbundance(sample, nucleus);
            accum.AddAbundance(compound->AbundanceOrdinal(sample, nucleus),
                               fitBar * (q > 0 ? termSeg[q]->GetComponentScaling() : 1.0) * v);
          }
          coef[q] = (q > 0 ? termSeg[q]->GetComponentScaling() : 1.0) * w;
        }
      } else {  // RATIO: base / first component
        const double baseVal = valueOf(pts[0], xsComp[0]);
        const double denomVal = pts[1]->GetFitCrossSection();
        if (natural) {
          const int nucleus = abundanceOf(segment->GetEntranceKey());
          accum.AddAbundance(compound->AbundanceOrdinal(sample, nucleus), fitBar * baseVal);
        }
        if (baseVal != 0.0 && denomVal != 0.0) {
          coef[0] = model / baseVal;
          coef[1] = -model / denomVal;
          xsComp[1] = 0;
        }
      }
      std::vector<EPoint *> barPts;
      std::vector<TMatrixCotangent> bars;
      for (size_t q = 0; q < pts.size(); q++) {
        if (coef[q] == 0.0) continue;
        useExitOf(pts[q]);
        TMatrixCotangent bar;
        if (!amf.PointTBar(pts[q], fitBar * coef[q], bar, xsComp[q])) return false;
        barPts.push_back(pts[q]);
        bars.push_back(std::move(bar));
      }
      amf.TBarAdjointMany(barPts, bars, accum, shiftDeriv);
      return true;
    }
  }

  // Combined model (components and/or target-effect integration).
  double model = segment->CalculateTheoreticalCrossSection(pointIdx, compound, config, data);
  outModel = model;
  double fitBar = fitBarFn(segment, i, pointIdx, model);
  if (fitBar == 0.0) return true;

  // The base term's own (unweighted) value, as the point still holds it.
  double baseVal = (baseComp == 1) ? point->GetFitE1CrossSection()
      : (baseComp == 2)            ? point->GetFitE2CrossSection()
                                   : point->GetFitCrossSection();
  double baseWeight = 1.0;
  if (natural) {
    const int nucleus = abundanceOf(segment->GetEntranceKey());
    baseWeight = compound->GetAbundance(sample, nucleus);
    accum.AddAbundance(compound->AbundanceOrdinal(sample, nucleus), fitBar * baseVal);
  }
  if (!hasComp) {
    return GradAdjointPoint(point, baseComp, fitBar * baseWeight, accum, shiftDeriv, compound, config);
  } else if (segment->GetOperationType() == SUM) {
    if (!GradAdjointPoint(point, baseComp, fitBar * baseWeight, accum, shiftDeriv, compound, config))
      return false;
    for (ESegment *comp : segment->GetComponentSegments()) {
      if (!comp || pointIdx >= comp->NumPoints()) continue;
      EPoint *cp = comp->GetPoint(pointIdx + 1);
      if (!cp) continue;
      double cv = cp->GetFitCrossSection();
      if (std::isnan(cv) || std::isinf(cv)) continue;
      double compWeight = 1.0;
      if (natural) {
        const int nucleus = abundanceOf(comp->GetEntranceKey());
        compWeight = compound->GetAbundance(sample, nucleus);
        accum.AddAbundance(compound->AbundanceOrdinal(sample, nucleus), fitBar * comp->GetComponentScaling() * cv);
      }
      if (!GradAdjointPoint(cp, comp->GetCrossSectionComponent(),
                            fitBar * comp->GetComponentScaling() * compWeight, accum,
                            shiftDeriv, compound, config))
        return false;
    }
    return true;
  } else {  // RATIO
    const std::vector<ESegment *> &comps = segment->GetComponentSegments();
    if (comps.empty() || !comps[0] || pointIdx >= comps[0]->NumPoints()) return false;
    EPoint *denomPoint = comps[0]->GetPoint(pointIdx + 1);
    double baseVal = (baseComp == 1) ? point->GetFitE1CrossSection()
        : (baseComp == 2)            ? point->GetFitE2CrossSection()
                                     : point->GetFitCrossSection();
    double denomVal = denomPoint ? denomPoint->GetFitCrossSection() : 0.0;
    if (baseVal != 0.0 && denomVal != 0.0 && denomPoint) {
      if (!GradAdjointPoint(point, baseComp, fitBar * model / baseVal,
                            accum, shiftDeriv, compound, config))
        return false;
      if (!GradAdjointPoint(denomPoint, 0, -fitBar * model / denomVal,
                            accum, shiftDeriv, compound, config))
        return false;
    }
    return true;
  }
}

// ---------------------------------------------------------------------------
// Energy groups.  Point mapping (EData::MapData) collects the points that share
// one energy-dependent calculation -- one level matrix and T-matrix -- within a
// segment and across segments.  The gradient works on the same groups: one
// forward solve per group, each member's observable from the shared T-matrix,
// and the back-propagation through the level matrix once per group.
// ---------------------------------------------------------------------------

struct GroupMember {
  ESegment *seg;
  int segNum;    // 1-based
  int pointIdx;  // 0-based
  EPoint *pt;
};

struct GradItem {
  std::vector<GroupMember> members;  // members[0] is the host
  bool grouped;                      // false: every member on its own (GradOnePoint)
};

// A point whose gradient can be taken from a shared T-matrix: no components
// (SUM/RATIO) and no sub-point target effect.  The others keep GradOnePoint.
bool Groupable(ESegment *seg, EPoint *pt) {
  if (seg->HasComponents()) return false;
  if (pt->IsTargetEffect()) {
    TargetEffect *te = pt->GetParentData()->GetTargetEffect(pt->GetTargetEffectNum());
    if (te && te->IsSubPointEffect()) return false;
  }
  return true;
}

bool GroupedGradientOff() {
  static const bool off = [] {
    const char *e = std::getenv("AZURE_NO_GROUPED_GRADIENT");
    return e && *e && std::string(e) != "0";
  }();
  return off;
}

// The work items of one gradient pass: each groupable host with the groupable
// points mapped onto it, and every other point on its own.  Largest first, so
// a big group does not end up last on one thread.
std::vector<GradItem> BuildGradItems(EData *data) {
  const bool off = GroupedGradientOff();
  std::unordered_map<const EPoint *, std::pair<int, int>> where;
  for (int s = 1; s <= data->NumSegments(); s++) {
    ESegment *seg = data->GetSegment(s);
    if (!seg) continue;
    for (int p = 1; p <= seg->NumPoints(); p++) where[seg->GetPoint(p)] = {s, p - 1};
  }
  auto hostOf = [&](EPoint *pt) -> EPoint * {
    if (!pt->IsMapped()) return nullptr;
    EnergyMap m = pt->GetMap();
    if (m.segment < 1 || m.segment > data->NumSegments()) return nullptr;
    ESegment *hs = data->GetSegment(m.segment);
    if (m.point < 1 || m.point > hs->NumPoints()) return nullptr;
    EPoint *h = hs->GetPoint(m.point);
    if (!Groupable(hs, h) || h->GetCMEnergy() != pt->GetCMEnergy()) return nullptr;
    for (int k = 1; k <= h->NumLocalMappedPoints(); k++)
      if (h->GetLocalMappedPoint(k) == pt) return h;
    return nullptr;
  };
  std::vector<GradItem> items;
  for (int s = 1; s <= data->NumSegments(); s++) {
    ESegment *seg = data->GetSegment(s);
    if (!seg) continue;
    for (int p = 1; p <= seg->NumPoints(); p++) {
      EPoint *pt = seg->GetPoint(p);
      if (!pt) continue;
      const bool g = !off && Groupable(seg, pt);
      if (g && hostOf(pt)) continue;  // its host's item calculates it
      GradItem it;
      it.grouped = g;
      it.members.push_back({seg, s, p - 1, pt});
      if (g)
        for (int k = 1; k <= pt->NumLocalMappedPoints(); k++) {
          EPoint *m = pt->GetLocalMappedPoint(k);
          auto w = where.find(m);
          if (w == where.end()) continue;
          ESegment *ms = data->GetSegment(w->second.first);
          if (Groupable(ms, m) && hostOf(m) == pt) it.members.push_back({ms, w->second.first, w->second.second, m});
        }
      items.push_back(std::move(it));
    }
  }
  std::stable_sort(items.begin(), items.end(),
                   [](const GradItem &a, const GradItem &b) { return a.members.size() > b.members.size(); });
  return items;
}

void AddCotangent(TMatrixCotangent &dst, const TMatrixCotangent &src) {
  if (dst.t.empty() && dst.ec.empty()) {
    dst = src;
    return;
  }
  for (size_t k = 0; k < src.t.size() && k < dst.t.size(); k++)
    for (size_t m = 0; m < src.t[k].size() && m < dst.t[k].size(); m++) dst.t[k][m] += src.t[k][m];
  for (size_t k = 0; k < src.ec.size() && k < dst.ec.size(); k++)
    for (size_t m = 0; m < src.ec[k].size() && m < dst.ec[k].size(); m++) dst.ec[k][m] += src.ec[k][m];
  dst.screen.insert(dst.screen.end(), src.screen.begin(), src.screen.end());
}

// The forward solve of a group, at its host's energy, and each member's model
// from the shared T-matrix (the plain-point forward of GradOnePoint).
// `rawModel` is the member's unweighted observable, `weight` its natural-target
// abundance (1 otherwise).
void GroupForward(AMatrixFunc &amf, const GradItem &it, CNuc *compound, std::vector<double> &rawModel,
                  std::vector<double> &weight) {
  EPoint *host = it.members[0].pt;
  AZUREProfile::Timer profileTimer(AZUREProfile::kGradForward);
  amf.ClearMatrices();
  amf.FillMatrices(host);
  amf.InvertMatrices();
  amf.CalculateTMatrix(host);
  const size_t g = it.members.size();
  rawModel.assign(g, 0.0);
  weight.assign(g, 1.0);
  for (size_t q = 0; q < g; q++) {
    const GroupMember &m = it.members[q];
    amf.GenMatrixFunc::CalculateCrossSection(m.pt);
    const int comp = m.seg->GetCrossSectionComponent();
    rawModel[q] = (comp == 1) ? m.pt->GetFitE1CrossSection() : (comp == 2) ? m.pt->GetFitE2CrossSection()
                                                                          : m.pt->GetFitCrossSection();
    if (m.seg->IsNatural()) {
      const int nucleus = compound->GetPair(compound->GetPairNumFromKey(m.seg->GetEntranceKey()))->GetNucleus();
      weight[q] = compound->GetAbundance(m.seg->GetNaturalSample(), nucleus);
    }
  }
}

int AbundanceOrdinalOf(CNuc *compound, ESegment *seg) {
  const int nucleus = compound->GetPair(compound->GetPairNumFromKey(seg->GetEntranceKey()))->GetNucleus();
  return compound->AbundanceOrdinal(seg->GetNaturalSample(), nucleus);
}

// The scalar-objective adjoint of a group: the members' T-matrix cotangents are
// summed (the back-propagation is linear in them) and back-propagated once.
bool GroupGradient(AMatrixFunc &amf, const GradItem &it, CNuc *compound, const vector_matrix_r *shiftDeriv,
                   const FitBarFn &fitBarFn, GradAccum &accum) {
  std::vector<double> rawModel, weight;
  GroupForward(amf, it, compound, rawModel, weight);
  TMatrixCotangent sum, bar;
  bool any = false;
  for (size_t q = 0; q < it.members.size(); q++) {
    const GroupMember &m = it.members[q];
    const double model = weight[q] * rawModel[q];
    const double fitBar = fitBarFn(m.seg, m.segNum, m.pointIdx, model);
    if (fitBar == 0.0) continue;
    if (m.seg->IsNatural()) accum.AddAbundance(AbundanceOrdinalOf(compound, m.seg), fitBar * rawModel[q]);
    if (!amf.PointTBar(m.pt, fitBar * weight[q], bar, m.seg->GetCrossSectionComponent())) return false;
    AddCotangent(sum, bar);
    any = true;
  }
  if (any) amf.TBarAdjoint(it.members[0].pt, sum, accum, shiftDeriv);
  return true;
}

// d(member model)/d(E, gamma) for every member of a group, as rows.  For each
// member `rowFn(q, model, accum)` receives the member's combined model and its
// gradient (scaled by fitBarOf(q, model)) in `accum`.
//
// The derivatives of the level matrix and T-matrix are taken once per group:
// each member's row is a combination, with its own angle-dependent
// coefficients, of the derivatives of the group's bilinear forms g_c^T A g_c'
// (AMatrixFunc::BilinearAdjoint with the cotangents 1 and i, the real and
// imaginary parts of the derivative).  A member then costs one contraction
// with that basis instead of a back-propagation through the level matrix.
// Where the basis would cost more than the members (few members, many channel
// pairs) or the decay has external capture, each member is back-propagated on
// its own -- still from the one shared forward solve.
template <class FitBarOf, class RowFn>
bool GroupRows(AMatrixFunc &amf, const GradItem &it, CNuc *compound, const vector_matrix_r *shiftDeriv,
               FitBarOf fitBarOf, RowFn rowFn, GradAccum &accum, GradAccum &scratch,
               const AMatrixFunc::FreeMask *mask) {
  std::vector<double> rawModel, weight;
  GroupForward(amf, it, compound, rawModel, weight);
  EPoint *host = it.members[0].pt;
  const size_t g = it.members.size();
  std::vector<double> model(g), fitBar(g);
  std::vector<TMatrixCotangent> bars(g);
  std::vector<std::vector<std::vector<BilinearTerm>>> terms(g);
  bool bilinear = true;
  for (size_t q = 0; q < g; q++) {
    const GroupMember &m = it.members[q];
    model[q] = weight[q] * rawModel[q];
    fitBar[q] = fitBarOf(q, model[q]);
    if (fitBar[q] == 0.0) continue;
    if (!amf.PointTBar(m.pt, fitBar[q] * weight[q], bars[q], m.seg->GetCrossSectionComponent())) return false;
    if (bilinear && !amf.CollectBilinearTerms(host, bars[q], terms[q])) bilinear = false;
  }

  // Cost of the two routes, in units of level-matrix elements.
  struct Key { int j, ch, chp; };
  std::vector<Key> keys;
  std::vector<std::vector<int>> keyOf(g);
  double costMembers = 0.0, costBasis = 0.0;
  const bool levelBasis = (mask != nullptr);
  std::vector<double> freeCount;
  if (levelBasis) {
    freeCount.assign(mask->e.size(), 0.0);
    for (size_t j = 0; j < mask->e.size(); j++)
      for (size_t la = 0; la < mask->e[j].size(); la++) {
        freeCount[j] += mask->e[j][la];
        for (char f : mask->g[j][la]) freeCount[j] += f;
      }
  }
  if (bilinear) {
    for (size_t q = 0; q < g; q++) {
      for (int j = 1; j <= (int)terms[q].size(); j++) {
        if (terms[q][j - 1].empty()) continue;
        JGroup *jg = compound->GetJGroup(j);
        const double n2 = (double)jg->NumLevels() * jg->NumLevels();
        costMembers += n2 * (terms[q][j - 1].size() + jg->NumChannels());
        for (const BilinearTerm &t : terms[q][j - 1]) {
          int k = 0;
          while (k < (int)keys.size() && !(keys[k].j == j && keys[k].ch == t.ch && keys[k].chp == t.chp)) k++;
          if (k == (int)keys.size()) {
            keys.push_back({j, t.ch, t.chp});
            // LevelBasis: O(n) per free parameter of the group; else the
            // structural adjoint, O(n^2 c).
            costBasis += (levelBasis ? 2.0 * freeCount[j - 1] * jg->NumLevels() : 2.0 * n2 * (1.0 + jg->NumChannels()));
          }
          keyOf[q].push_back(k);
          costBasis += 2.0 * jg->NumLevels() * (1.0 + jg->NumChannels());
        }
      }
    }
  }
  const bool useBasis = bilinear && costBasis < costMembers;

  if (!useBasis) {
    for (size_t q = 0; q < g; q++) {
      const GroupMember &m = it.members[q];
      accum.Zero();
      if (fitBar[q] != 0.0) {
        if (m.seg->IsNatural()) accum.AddAbundance(AbundanceOrdinalOf(compound, m.seg), fitBar[q] * rawModel[q]);
        amf.TBarAdjoint(host, bars[q], accum, shiftDeriv);
      }
      rowFn(q, model[q], accum);
    }
    return true;
  }

  // The basis: for each channel pair, the J-group block (E then gamma, level by
  // level) of the adjoint with cotangent 1 and with cotangent i.
  std::vector<std::vector<double>> basisRe(keys.size()), basisIm(keys.size());
  AZUREProfile::Timer profileBasis(AZUREProfile::kGradBackprop);
  for (size_t k = 0; k < keys.size(); k++) {
    const int j = keys[k].j;
    JGroup *jg = compound->GetJGroup(j);
    const int n = jg->NumLevels(), nch = jg->NumChannels();
    for (int part = 0; part < 2; part++) {
      std::vector<double> &block = part ? basisIm[k] : basisRe[k];
      // The rank-two derivatives of the level matrix by free parameter
      // (LevelBasis), O(n) each, where they apply.
      const BilinearTerm term{keys[k].ch, keys[k].chp, part ? complex(0.0, 1.0) : complex(1.0, 0.0)};
      if (mask && amf.LevelBasis(host, j, term, *mask, shiftDeriv, block)) continue;
      for (int la = 0; la < n; la++) {
        scratch.e[j - 1][la] = 0.0;
        std::fill(scratch.gamma[j - 1][la].begin(), scratch.gamma[j - 1][la].end(), 0.0);
      }
      std::vector<BilinearTerm> one{{keys[k].ch, keys[k].chp, part ? complex(0.0, 1.0) : complex(1.0, 0.0)}};
      amf.BilinearAdjoint(host, j, one, scratch, shiftDeriv);
      block.assign((size_t)n * (1 + nch), 0.0);
      for (int la = 0; la < n; la++) {
        block[la] = scratch.e[j - 1][la];
        for (int d = 0; d < nch; d++) block[n + (size_t)la * nch + d] = scratch.gamma[j - 1][la][d];
      }
    }
  }
  profileBasis.Stop();
  AZUREProfile::Timer profileRows(AZUREProfile::kGradRows);
  for (size_t q = 0; q < g; q++) {
    const GroupMember &m = it.members[q];
    accum.Zero();
    if (fitBar[q] != 0.0) {
      if (m.seg->IsNatural()) accum.AddAbundance(AbundanceOrdinalOf(compound, m.seg), fitBar[q] * rawModel[q]);
      // The electron-screening derivative is not in the bilinear forms (TBarAdjoint adds it on the other route).
      for (const auto &s : bars[q].screen) accum.AddScreening(s.first, s.second);
      size_t t = 0;
      for (int j = 1; j <= (int)terms[q].size(); j++)
        for (const BilinearTerm &term : terms[q][j - 1]) {
          const int k = keyOf[q][t++];
          const double a = std::real(term.beta), b = std::imag(term.beta);
          const std::vector<double> &re = basisRe[k], &im = basisIm[k];
          JGroup *jg = compound->GetJGroup(j);
          const int n = jg->NumLevels(), nch = jg->NumChannels();
          for (int la = 0; la < n; la++) {
            accum.e[j - 1][la] += a * re[la] + b * im[la];
            double *gr = &accum.gamma[j - 1][la][0];
            const double *r0 = &re[n + (size_t)la * nch], *i0 = &im[n + (size_t)la * nch];
            for (int d = 0; d < nch; d++) gr[d] += a * r0[d] + b * i0[d];
          }
        }
    }
    rowFn(q, model[q], accum);
  }
  return true;
}

void ReportBail(const char *what, ESegment *segment, int i) {
  if (!std::getenv("AZURE_GRAD_DEBUG") || !segment) return;
  EPoint *pt = segment->GetPoint(1);
  std::cerr << "[grad] " << what << " adjoint bailed on segment " << i << " (key " << segment->GetSegmentKey()
            << ", entrance " << segment->GetEntranceKey() << ", exit " << segment->GetExitKey()
            << ", components " << segment->HasComponents() << ", angular-dist " << (pt && pt->IsAngularDist())
            << ", A_y " << (pt && pt->IsAnalyzingPower()) << ", target effect " << (pt && pt->IsTargetEffect())
            << ")" << std::endl;
}

}  // namespace

bool AccumulateEGammaGradient(CNuc *compound, EData *data, const Config &config,
                              const ParamIndexMap &pmap,
                              const vector_matrix_r *shiftDeriv,
                              const FitBarFn &fitBarFn,
                              GradAccum &accum) {
  // The work items (energy groups and single points, BuildGradItems) are
  // independent and only read the shared `compound`, so they run in parallel.
  // The one shared write is the accumulation into `accum`, so each thread
  // accumulates into its OWN GradAccum and they are summed once at the end (a
  // reduction).  The norm gradient is accumulated through `fitBarFn`'s side
  // effect, which the callers guard with `#pragma omp atomic`.
#ifdef _OPENMP
  const int nThreads = std::max(1, omp_get_max_threads());
#else
  const int nThreads = 1;
#endif
  std::vector<GradAccum> threadAccum(nThreads);
  for (GradAccum &ta : threadAccum) {
    ta.Init(compound);
    // A scalar pass: the Brune and Park level-matrix adjoints are summed in the
    // formal parameters and resolved once below.
    ta.formal = true;
  }

  const std::vector<GradItem> items = BuildGradItems(data);
  bool bail = false;
  int bailSeg = 0;
#pragma omp parallel shared(bail, bailSeg)
  {
    AMatrixFunc amf(compound, config);
#pragma omp for schedule(dynamic, 1)
    for (int k = 0; k < (int)items.size(); k++) {
      if (bail) continue;
#ifdef _OPENMP
      GradAccum &local = threadAccum[omp_get_thread_num()];
#else
      GradAccum &local = threadAccum[0];
#endif
      const GradItem &it = items[k];
      bool ok = true;
      if (it.grouped && it.members.size() > 1) {
        ok = GroupGradient(amf, it, compound, shiftDeriv, fitBarFn, local);
      } else {
        for (const GroupMember &m : it.members) {
          double model;
          if (!GradOnePoint(m.seg, data, m.segNum, m.pointIdx, compound, config, shiftDeriv, fitBarFn, local, model)) {
            ok = false;
            break;
          }
        }
      }
      if (!ok) {
#pragma omp critical
        {
          if (!bail) bailSeg = it.members[0].segNum;
          bail = true;
        }
      }
    }
  }
  if (bail) {
    if (std::getenv("AZURE_GRAD_DEBUG")) {
      static std::atomic<bool> told{false};
      if (!told.exchange(true)) {
        bool inactiveLevel = false;
        for (int j = 1; j <= compound->NumJGroups(); j++) {
          JGroup *jg = compound->GetJGroup(j);
          if (!jg->IsInRMatrix()) continue;
          for (int la = 1; la <= jg->NumLevels(); la++)
            if (!jg->GetLevel(la)->IsInRMatrix()) inactiveLevel = true;
        }
        ESegment *segment = data->GetSegment(bailSeg);
        EPoint *pt = segment ? segment->GetPoint(1) : nullptr;
        std::cerr << "[grad] adjoint bailed on segment " << bailSeg << " (entrance "
                  << (segment ? segment->GetEntranceKey() : 0) << ", exit " << (segment ? segment->GetExitKey() : 0)
                  << ", components " << (segment && segment->HasComponents()) << ", angular-dist "
                  << (pt && pt->IsAngularDist()) << ", A_y " << (pt && pt->IsAnalyzingPower()) << ", target effect "
                  << (pt && pt->IsTargetEffect()) << "); a level outside the R matrix in an "
                  << "R-matrix J group: " << (inactiveLevel ? "YES" : "no") << std::endl;
      }
    }
    // On bail the caller ignores `accum` and falls back to finite differences,
    // so leave it untouched.
    return false;
  }
  for (const GradAccum &ta : threadAccum) accum.Add(ta);
  if (accum.formalUsed) {
    AMatrixFunc amf(compound, config);
    amf.ResolveFormalCotangents(accum, shiftDeriv);
  }
  return true;
}

namespace {

// Rows of d(model)/d(theta) for every point, through the energy groups:
// writeRow(member, model, accum) receives each point's combined model and
// fitBarOf(member) times its gradient.  Shared by the residual Jacobian and the
// model gradients.  Returns false, with `bailSeg` set, if a point is outside
// the analytic path.
// The free level parameters of the packed vector, per J-group (FormalBasis).
AMatrixFunc::FreeMask BuildFreeMask(CNuc *compound, const ParamIndexMap &pmap) {
  AMatrixFunc::FreeMask mask;
  const int nJ = compound->NumJGroups();
  mask.e.assign(nJ, {});
  mask.g.assign(nJ, {});
  auto isFree = [&pmap](int full) { return full >= 0 && full < pmap.NumFull() && pmap.FullToPacked(full) >= 0; };
  for (int j = 1; j <= nJ; j++) {
    JGroup *jg = compound->GetJGroup(j);
    const int nl = jg->NumLevels(), nch = jg->NumChannels();
    mask.e[j - 1].assign(nl, 0);
    mask.g[j - 1].assign(nl, std::vector<char>(nch, 0));
    for (int la = 1; la <= nl; la++) {
      mask.e[j - 1][la - 1] = isFree(pmap.EIndexFlat(j, la));
      for (int ch = 1; ch <= nch; ch++) mask.g[j - 1][la - 1][ch - 1] = isFree(pmap.GammaIndexFlat(j, la, ch));
    }
  }
  return mask;
}

template <class FitBarOf, class WriteRow>
bool ForEachRow(CNuc *compound, EData *data, const Config &config, const vector_matrix_r *shiftDeriv,
                FitBarOf fitBarOf, WriteRow writeRow, int &bailSeg, int &numGroups,
                const AMatrixFunc::FreeMask *mask) {
  const std::vector<GradItem> items = BuildGradItems(data);
  numGroups = 0;
  for (const GradItem &it : items)
    if (it.grouped && it.members.size() > 1) numGroups++;
  bool bail = false;
  bailSeg = 0;
#pragma omp parallel shared(bail, bailSeg)
  {
    GradAccum accum, scratch;  // per-thread, reused across its rows
    accum.Init(compound);
    scratch.Init(compound);
    accum.mask = mask;
    AMatrixFunc amf(compound, config);
#pragma omp for schedule(dynamic, 1)
    for (int k = 0; k < (int)items.size(); k++) {
      if (bail) continue;
      const GradItem &it = items[k];
      bool ok = true;
      if (it.grouped && it.members.size() > 1) {
        auto fb = [&](size_t q, double) { return fitBarOf(it.members[q]); };
        auto row = [&](size_t q, double model, GradAccum &acc) { writeRow(it.members[q], model, acc); };
        ok = GroupRows(amf, it, compound, shiftDeriv, fb, row, accum, scratch, mask);
      } else {
        for (const GroupMember &m : it.members) {
          const double fitBar = fitBarOf(m);
          FitBarFn oneFitBar = [fitBar](ESegment *, int, int, double) -> double { return fitBar; };
          accum.Zero();
          double model = 0.0;
          if (!GradOnePoint(m.seg, data, m.segNum, m.pointIdx, compound, config, shiftDeriv, oneFitBar, accum, model)) {
            ok = false;
            break;
          }
          writeRow(m, model, accum);
        }
      }
      if (!ok) {
#pragma omp critical
        {
          if (!bail) bailSeg = it.members[0].segNum;
          bail = true;
        }
      }
    }
  }
  return !bail;
}

}  // namespace

bool ComputeResidualJacobian(CNuc *compound, EData *data, const Config &config,
                             const ParamIndexMap &pmap,
                             const vector_matrix_r *shiftDeriv,
                             vector_r &residuals, vector_r &jacobian, int &nCols) {
  // Columns of J are the non-fixed (packed) parameters; rows are data points.
  nCols = pmap.NumPacked();
  residuals.clear();
  jacobian.clear();

  // Each (non-null) point is one independent row, so pre-assign every point a
  // global row index (mirroring the serial null-skip) and pre-size the outputs.
  // Rows can then be filled in parallel with no contention -- each work item
  // writes only the rows of its own points.
  const int nSeg = data->NumSegments();
  std::vector<std::vector<int>> rowOf(nSeg + 1);
  std::vector<double> norm(nSeg + 1, 0.0);
  std::vector<int> normFull(nSeg + 1, -1);
  int totalRows = 0;
  for (int i = 1; i <= nSeg; i++) {
    ESegment *segment = data->GetSegment(i);
    if (!segment) continue;
    const int nPoints = segment->NumPoints();
    rowOf[i].assign(nPoints, -1);
    for (int pointIdx = 0; pointIdx < nPoints; pointIdx++)
      if (segment->GetPoint(pointIdx + 1)) rowOf[i][pointIdx] = totalRows++;
    norm[i] = segment->GetNorm();
    // A norm-group follower's rows depend on its leader's norm parameter.
    const int normOwner = data->NormOwner(i);
    normFull[i] = data->GetSegment(normOwner)->IsVaryNorm() ? pmap.NormIndex(normOwner) : -1;
  }
  residuals.assign(totalRows, 0.0);
  jacobian.assign((size_t)totalRows * nCols, 0.0);

  const bool timing = std::getenv("AZURE_GRAD_TIMING") != nullptr;
  const long long comb0 = g_teCombinerNs.load(), adj0 = g_teAdjointNs.load();
  const auto ts0 = std::chrono::steady_clock::now();

  // Residual r = (model - data*n)/(cmErr*n);  d r / d model = 1/(cmErr*n), the
  // cotangent that yields d r / d(E,gamma) directly.
  auto fitBarOf = [&](const GroupMember &m) {
    const double denom = m.pt->GetCMCrossSectionError() * norm[m.segNum];
    return (denom != 0.0) ? 1.0 / denom : 0.0;
  };
  static thread_local vector_r fullRow;
  auto writeRow = [&](const GroupMember &m, double model, const GradAccum &accum) {
    const int row = rowOf[m.segNum][m.pointIdx];
    if (row < 0) return;
    const double n = norm[m.segNum];
    const double cmErr = m.pt->GetCMCrossSectionError();
    const double denom = cmErr * n;
    residuals[row] = (denom != 0.0) ? (model - m.pt->GetCMCrossSection() * n) / denom : 0.0;
    // Scatter this row's d r / d(E,gamma) into the packed columns.
    fullRow.assign(pmap.NumFull(), 0.0);
    accum.Scatter(pmap, fullRow);
    // Norm column:  d r / d n = -model/(cmErr*n^2).
    if (normFull[m.segNum] >= 0 && cmErr != 0.0 && n != 0.0) fullRow[normFull[m.segNum]] += -model / (cmErr * n * n);
    const size_t base = (size_t)row * nCols;
    for (int f = 0; f < pmap.NumFull(); f++) {
      int packed = pmap.FullToPacked(f);
      if (packed >= 0 && packed < nCols) jacobian[base + packed] = fullRow[f];
    }
  };
  int bailSeg = 0, numGroups = 0;
  const AMatrixFunc::FreeMask mask = BuildFreeMask(compound, pmap);
  const bool ok = ForEachRow(compound, data, config, shiftDeriv, fitBarOf, writeRow, bailSeg, numGroups, &mask);
  if (!ok) ReportBail("Jacobian/gradient", data->GetSegment(bailSeg), bailSeg);

  if (timing) {
    static std::atomic<int> calls{0};
    const int c = ++calls;
    if (c <= 3 || c % 50 == 0)
      std::cerr << "[grad-timing] Jacobian rows " << c << ": " << totalRows << " rows, " << numGroups
                << " energy groups, "
                << std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ts0).count()
                << " ms; target-effect combiner weights " << (g_teCombinerNs.load() - comb0) / 1.0e6
                << " ms, sub-point adjoints " << (g_teAdjointNs.load() - adj0) / 1.0e6 << " ms (thread-summed)"
                << std::endl;
  }

  if (!ok) {
    residuals.clear();
    jacobian.clear();
  }
  return ok;
}

bool ComputeModelGradients(CNuc *compound, EData *data, const Config &config,
                           const ParamIndexMap &pmap,
                           const vector_matrix_r *shiftDeriv,
                           std::map<EPoint *, vector_r> &gradByPoint) {
  gradByPoint.clear();
  const int nCols = pmap.NumPacked();
  const int nSeg = data->NumSegments();

  // Pre-create every point's (zero-initialized) packed row so the parallel loop
  // only writes into slots that already exist -- std::map insertion is not
  // thread-safe, so all keys must be present before the parallel region.
  for (int i = 1; i <= nSeg; i++) {
    ESegment *segment = data->GetSegment(i);
    if (!segment) continue;
    for (int pointIdx = 0; pointIdx < segment->NumPoints(); pointIdx++) {
      EPoint *pt = segment->GetPoint(pointIdx + 1);
      if (pt) gradByPoint[pt] = vector_r(nCols, 0.0);
    }
  }

  // fitBar = d(model)/d(model) = 1 gives d(model)/d(E,gamma) directly.
  auto fitBarOf = [](const GroupMember &) { return 1.0; };
  static thread_local vector_r fullRow;
  auto writeRow = [&](const GroupMember &m, double, const GradAccum &accum) {
    fullRow.assign(pmap.NumFull(), 0.0);
    accum.Scatter(pmap, fullRow);
    vector_r &row = gradByPoint[m.pt];  // key pre-inserted above
    for (int f = 0; f < pmap.NumFull(); f++) {
      int packed = pmap.FullToPacked(f);
      if (packed >= 0 && packed < nCols) row[packed] = fullRow[f];
    }
  };
  int bailSeg = 0, numGroups = 0;
  const AMatrixFunc::FreeMask mask = BuildFreeMask(compound, pmap);
  const bool ok = ForEachRow(compound, data, config, shiftDeriv, fitBarOf, writeRow, bailSeg, numGroups, &mask);
  if (!ok) {
    ReportBail("Jacobian/gradient", data->GetSegment(bailSeg), bailSeg);
    gradByPoint.clear();
  }
  return ok;
}
