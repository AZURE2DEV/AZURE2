#include "AZUREGrad.h"
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
#ifdef _OPENMP
#include <omp.h>
#endif

/*!
 * Precompute dS_ch/dE_level for every (jGroup, level, channel).  Mirrors the
 * bound/unbound branch of CNuc::CalcShiftFunctions.
 */
vector_matrix_r BuildShiftDerivTable(CNuc *compound, const Config &configure) {
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
  abundance.assign(compound->NumAbundanceParams(), 0.0);
}

/*!
 * Resets all accumulators to zero without reallocating.
 */
void GradAccum::Zero() {
  for (auto &row : e) std::fill(row.begin(), row.end(), 0.0);
  for (auto &lev : gamma)
    for (auto &row : lev) std::fill(row.begin(), row.end(), 0.0);
  std::fill(abundance.begin(), abundance.end(), 0.0);
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
}

/*!
 * Scatters the physics-coordinate gradients into the flat parameter vector.
 */
void GradAccum::Scatter(const ParamIndexMap &map, vector_r &gradFull) const {
  for (int j = 0; j < (int)e.size(); j++) {
    for (int la = 0; la < (int)e[j].size(); la++) {
      int ei = map.EIndex(j + 1, la + 1);
      if (ei >= 0 && ei < (int)gradFull.size()) gradFull[ei] += e[j][la];
      for (int ch = 0; ch < (int)gamma[j][la].size(); ch++) {
        int gi = map.GammaIndex(j + 1, la + 1, ch + 1);
        if (gi >= 0 && gi < (int)gradFull.size()) gradFull[gi] += gamma[j][la][ch];
      }
    }
  }
  for (int k = 0; k < (int)abundance.size(); k++) {
    int ai = map.AbundanceIndex(k + 1);
    if (ai >= 0 && ai < (int)gradFull.size()) gradFull[ai] += abundance[k];
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
  for (int j = 1; j <= compound->NumJGroups(); j++) {
    JGroup *jGroup = compound->GetJGroup(j);
    for (int la = 1; la <= jGroup->NumLevels(); la++) {
      // Energy parameter.
      map.eIndex_[{j, la}] = i;
      map.desc_.push_back(ParamDesc{ParamKind::LevelEnergy, j, la, -1, -1, false});
      i++;
      // One gamma per channel.
      for (int ch = 1; ch <= jGroup->NumChannels(); ch++) {
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

}  // namespace

bool AccumulateEGammaGradient(CNuc *compound, EData *data, const Config &config,
                              const ParamIndexMap &pmap,
                              const vector_matrix_r *shiftDeriv,
                              const FitBarFn &fitBarFn,
                              GradAccum &accum) {
  // Each point's forward+adjoint solve is independent and only reads the shared
  // `compound`, exactly like the forward cross-section loops in EData.cpp, so we
  // parallelize each segment's point loop the same way.  The one shared write is
  // the per-point accumulation into `accum`, so each thread accumulates into its
  // OWN GradAccum and they are summed once at the end (a reduction).  The norm
  // gradient is accumulated through `fitBarFn`'s side effect, which the callers
  // guard with `#pragma omp atomic`.
#ifdef _OPENMP
  const int nThreads = std::max(1, omp_get_max_threads());
#else
  const int nThreads = 1;
#endif
  std::vector<GradAccum> threadAccum(nThreads);
  for (GradAccum &ta : threadAccum) ta.Init(compound);

  bool ok = true;
  for (int i = 1; i <= data->NumSegments() && ok; i++) {
    ESegment *segment = data->GetSegment(i);
    if (!segment) continue;
    const int nPoints = segment->NumPoints();
    bool bail = false;
#pragma omp parallel for shared(bail)
    for (int pointIdx = 0; pointIdx < nPoints; pointIdx++) {
      if (bail) continue;
#ifdef _OPENMP
      GradAccum &local = threadAccum[omp_get_thread_num()];
#else
      GradAccum &local = threadAccum[0];
#endif
      double model;
      if (!GradOnePoint(segment, data, i, pointIdx, compound, config, shiftDeriv,
                        fitBarFn, local, model)) {
#pragma omp critical
        bail = true;
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
          EPoint *pt = segment->GetPoint(1);
          std::cerr << "[grad] adjoint bailed on segment " << i << " (entrance " << segment->GetEntranceKey()
                    << ", exit " << segment->GetExitKey() << ", components " << segment->HasComponents()
                    << ", angular-dist " << (pt && pt->IsAngularDist()) << ", A_y " << (pt && pt->IsAnalyzingPower())
                    << ", target effect " << (pt && pt->IsTargetEffect()) << "); a level outside the R matrix in an "
                    << "R-matrix J group: " << (inactiveLevel ? "YES" : "no") << std::endl;
        }
      }
      ok = false;
      break;
    }
  }

  // Reduce the per-thread partials only on success; on bail the caller ignores
  // `accum` and falls back to finite differences, so leave it untouched.
  if (ok)
    for (const GradAccum &ta : threadAccum) accum.Add(ta);
  return ok;
}

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
  // Rows can then be filled in parallel with no contention -- each thread owns
  // its own accumulator and scatter buffer and writes only its own rows.
  const int nSeg = data->NumSegments();
  std::vector<std::vector<int>> rowOf(nSeg + 1);
  int totalRows = 0;
  for (int i = 1; i <= nSeg; i++) {
    ESegment *segment = data->GetSegment(i);
    if (!segment) continue;
    const int nPoints = segment->NumPoints();
    rowOf[i].assign(nPoints, -1);
    for (int pointIdx = 0; pointIdx < nPoints; pointIdx++)
      if (segment->GetPoint(pointIdx + 1)) rowOf[i][pointIdx] = totalRows++;
  }
  residuals.assign(totalRows, 0.0);
  jacobian.assign((size_t)totalRows * nCols, 0.0);

  bool ok = true;
  const bool timing = std::getenv("AZURE_GRAD_TIMING") != nullptr;
  const long long comb0 = g_teCombinerNs.load(), adj0 = g_teAdjointNs.load();
  std::vector<std::pair<double, int>> segMs;
  for (int i = 1; i <= nSeg && ok; i++) {
    ESegment *segment = data->GetSegment(i);
    if (!segment) continue;
    const auto ts0 = std::chrono::steady_clock::now();
    const double norm = segment->GetNorm();
    // A norm-group follower's rows depend on its leader's norm parameter.
    const int normOwner = data->NormOwner(i);
    const int normFull = data->GetSegment(normOwner)->IsVaryNorm() ? pmap.NormIndex(normOwner) : -1;
    const int nPoints = segment->NumPoints();
    bool bail = false;

#pragma omp parallel
    {
      GradAccum accum;  // per-thread, reused across its rows
      accum.Init(compound);
      vector_r fullRow(pmap.NumFull(), 0.0);
#pragma omp for
      for (int pointIdx = 0; pointIdx < nPoints; pointIdx++) {
        if (bail) continue;
        const int row = rowOf[i][pointIdx];
        if (row < 0) continue;
        EPoint *pt = segment->GetPoint(pointIdx + 1);
        double cmErr = pt->GetCMCrossSectionError();
        double dataval = pt->GetCMCrossSection();

        // Residual r = (model - data*n)/(cmErr*n);  d r / d model = 1/(cmErr*n).
        // PointAdjoint with this cotangent yields d r / d(E,gamma) directly.
        double denom = cmErr * norm;
        FitBarFn resFitBar = [denom](ESegment *, int, int, double) -> double {
          return (denom != 0.0) ? 1.0 / denom : 0.0;
        };

        accum.Zero();
        double model = 0.0;
        if (!GradOnePoint(segment, data, i, pointIdx, compound, config, shiftDeriv,
                          resFitBar, accum, model)) {
#pragma omp critical
          bail = true;
          continue;
        }

        residuals[row] = (denom != 0.0) ? (model - dataval * norm) / denom : 0.0;

        // Scatter this row's d r / d(E,gamma) into the packed columns.
        std::fill(fullRow.begin(), fullRow.end(), 0.0);
        accum.Scatter(pmap, fullRow);
        // Norm column:  d r / d n = -model/(cmErr*n^2).
        if (normFull >= 0 && cmErr != 0.0 && norm != 0.0)
          fullRow[normFull] += -model / (cmErr * norm * norm);

        const size_t base = (size_t)row * nCols;
        for (int f = 0; f < pmap.NumFull(); f++) {
          int packed = pmap.FullToPacked(f);
          if (packed >= 0 && packed < nCols) jacobian[base + packed] = fullRow[f];
        }
      }
    }
    if (bail) {
      ok = false;
      if (std::getenv("AZURE_GRAD_DEBUG")) {
        EPoint *pt = segment->GetPoint(1);
        std::cerr << "[grad] Jacobian/gradient adjoint bailed on segment " << i << " (key " << segment->GetSegmentKey()
                  << ", entrance " << segment->GetEntranceKey() << ", exit " << segment->GetExitKey()
                  << ", components " << segment->HasComponents() << ", angular-dist " << (pt && pt->IsAngularDist())
                  << ", A_y " << (pt && pt->IsAnalyzingPower()) << ", target effect " << (pt && pt->IsTargetEffect())
                  << ")" << std::endl;
      }
    }
    if (timing)
      segMs.push_back({std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ts0).count(),
                       segment->GetSegmentKey()});
  }
  if (timing) {
    static std::atomic<int> calls{0};
    const int c = ++calls;
    if (c <= 3 || c % 50 == 0) {
      std::sort(segMs.begin(), segMs.end(), [](const std::pair<double, int> &a, const std::pair<double, int> &b) { return a.first > b.first; });
      std::cerr << "[grad-timing] Jacobian rows " << c << ": target-effect combiner weights "
                << (g_teCombinerNs.load() - comb0) / 1.0e6 << " ms, sub-point adjoints " << (g_teAdjointNs.load() - adj0) / 1.0e6
                << " ms (thread-summed); slowest segments (key:ms)";
      for (size_t k = 0; k < segMs.size() && k < 6; k++) std::cerr << " " << segMs[k].second << ":" << (long)segMs[k].first;
      std::cerr << std::endl;
    }
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
  std::vector<std::vector<EPoint *>> ptsOf(nSeg + 1);
  for (int i = 1; i <= nSeg; i++) {
    ESegment *segment = data->GetSegment(i);
    if (!segment) continue;
    const int nPoints = segment->NumPoints();
    ptsOf[i].assign(nPoints, nullptr);
    for (int pointIdx = 0; pointIdx < nPoints; pointIdx++) {
      EPoint *pt = segment->GetPoint(pointIdx + 1);
      ptsOf[i][pointIdx] = pt;
      if (pt) gradByPoint[pt] = vector_r(nCols, 0.0);
    }
  }

  bool ok = true;
  for (int i = 1; i <= nSeg && ok; i++) {
    ESegment *segment = data->GetSegment(i);
    if (!segment) continue;
    const int nPoints = segment->NumPoints();
    bool bail = false;

#pragma omp parallel
    {
      GradAccum accum;  // per-thread, reused across its rows
      accum.Init(compound);
      vector_r fullRow(pmap.NumFull(), 0.0);
#pragma omp for
      for (int pointIdx = 0; pointIdx < nPoints; pointIdx++) {
        if (bail) continue;
        EPoint *pt = ptsOf[i][pointIdx];
        if (!pt) continue;

        // fitBar = d(model)/d(model) = 1 gives d(model)/d(E,gamma) directly.
        FitBarFn oneFitBar = [](ESegment *, int, int, double) -> double { return 1.0; };

        accum.Zero();
        double model = 0.0;
        if (!GradOnePoint(segment, data, i, pointIdx, compound, config, shiftDeriv,
                          oneFitBar, accum, model)) {
#pragma omp critical
          bail = true;
          continue;
        }

        std::fill(fullRow.begin(), fullRow.end(), 0.0);
        accum.Scatter(pmap, fullRow);

        vector_r &row = gradByPoint[pt];  // key pre-inserted above
        for (int f = 0; f < pmap.NumFull(); f++) {
          int packed = pmap.FullToPacked(f);
          if (packed >= 0 && packed < nCols) row[packed] = fullRow[f];
        }
      }
    }
    if (bail) {
      ok = false;
      if (std::getenv("AZURE_GRAD_DEBUG")) {
        EPoint *pt = segment->GetPoint(1);
        std::cerr << "[grad] Jacobian/gradient adjoint bailed on segment " << i << " (key " << segment->GetSegmentKey()
                  << ", entrance " << segment->GetEntranceKey() << ", exit " << segment->GetExitKey()
                  << ", components " << segment->HasComponents() << ", angular-dist " << (pt && pt->IsAngularDist())
                  << ", A_y " << (pt && pt->IsAnalyzingPower()) << ", target effect " << (pt && pt->IsTargetEffect())
                  << ")" << std::endl;
      }
    }
  }

  if (!ok) gradByPoint.clear();
  return ok;
}
