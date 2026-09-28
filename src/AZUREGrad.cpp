#include "AZUREGrad.h"
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
#include "ChannelFunc.h"
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
        double dS = ChannelFunc(pair, !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC))
                        .ShiftDerivative(l, resonanceEnergy);
        tab[j - 1][la - 1][ch - 1] = dS;
      }
    }
  }
  return tab;
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
}

/*!
 * Resets all accumulators to zero without reallocating.
 */
void GradAccum::Zero() {
  for (auto &row : e) std::fill(row.begin(), row.end(), 0.0);
  for (auto &lev : gamma)
    for (auto &row : lev) std::fill(row.begin(), row.end(), 0.0);
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
 * Builds the parameter-index map by mirroring, in the exact same order, the
 * assignment loops in:
 *   - CNuc::FillCompoundFromParams      (level energies + gammas)
 *   - EData::FillNormsFromParams        (norms, only for IsVaryNorm segments)
 *   - EData::FillEnergyShiftsFromParams (one energy shift per segment)
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

  // --- Norm parameters: one per varying-norm segment. ---
  map.normOffset_ = i;
  if (data) {
    for (int s = 1; s <= data->NumSegments(); s++) {
      ESegment *segment = data->GetSegment(s);
      if (segment && segment->IsVaryNorm() && !segment->IsProfiledNorm()) {
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

// Forward + adjoint for one target-effect (sub-point integrated) point.  Its
// yield is linear in the sub-point cross sections, so the combination weights
// c_i = d(yield)/d(sigma_i) come from central differences on the cheap
// IntegrateTargetEffect combiner; each sub-point is then back-propagated.
bool GradTargetEffectAdjoint(EPoint *point, double fitBar, GradAccum &accum,
                             const vector_matrix_r *shiftDeriv, int xsComponent,
                             CNuc *compound, const Config &config) {
  int nsp = point->NumSubPoints();
  if (nsp <= 0) return true;
  // A target-averaged analyzing power is a ratio of two integrals,
  // Int A_y sigma dE / Int sigma dE (EPoint::IntegrateTargetEffectForObservable),
  // not a fixed linear combination of the sub-point values.  The weights below
  // would treat it as one, and each sub-point's own A_y adjoint would then be
  // combined without the cross-section weighting -- a wrong Jacobian with no
  // diagnostic.  PointAdjoint refuses such a point, but it is only ever called
  // here on the sub-points, so the refusal has to be made at this level.
  if (point->IsAnalyzingPower()) return false;

  std::vector<double> sigma(nsp);
  for (int i = 1; i <= nsp; i++) sigma[i - 1] = point->GetSubPoint(i)->GetFitCrossSection();
  double model = point->GetFitCrossSection();

  std::vector<double> c(nsp, 0.0);
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
    c[i - 1] = (yp - ym) / (2.0 * h);
  }
  point->SetFitCrossSection(model);

  bool ok = true;
  for (int i = 1; i <= nsp && ok; i++) {
    if (c[i - 1] == 0.0) continue;
    EPoint *sp = point->GetSubPoint(i);
    AMatrixFunc amf(compound, config);
    amf.ClearMatrices();
    amf.FillMatrices(sp);
    amf.InvertMatrices();
    amf.CalculateTMatrix(sp);
    if (!amf.PointAdjoint(sp, fitBar * c[i - 1], accum, shiftDeriv, xsComponent)) ok = false;
  }
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
  EPoint *point = segment->GetPoint(pointIdx + 1);
  outModel = 0.0;
  if (!point) return true;

  int baseComp = segment->GetCrossSectionComponent();
  bool hasComp = segment->HasComponents();

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

  if (!hasComp) {
    return GradAdjointPoint(point, baseComp, fitBar, accum, shiftDeriv, compound, config);
  } else if (segment->GetOperationType() == SUM) {
    if (!GradAdjointPoint(point, baseComp, fitBar, accum, shiftDeriv, compound, config))
      return false;
    for (ESegment *comp : segment->GetComponentSegments()) {
      if (!comp || pointIdx >= comp->NumPoints()) continue;
      EPoint *cp = comp->GetPoint(pointIdx + 1);
      if (!cp) continue;
      double cv = cp->GetFitCrossSection();
      if (std::isnan(cv) || std::isinf(cv)) continue;
      if (!GradAdjointPoint(cp, comp->GetCrossSectionComponent(),
                            fitBar * comp->GetComponentScaling(), accum,
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
    // THM (HOES) segments are computed by THMMatrixFunc, not the T-matrix
    // observable this adjoint differentiates; skip them here. Their
    // contribution to the gradient is added in AZURECalc::Gradient by finite
    // differences of a chi^2 restricted to THM segments, so mixed THM+direct
    // datasets keep the analytic speed-up on the direct part.
    if (segment->IsTHM()) continue;
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
                             vector_r &residuals, vector_r &jacobian, int &nCols,
                             bool skipTHM) {
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
  for (int i = 1; i <= nSeg && ok; i++) {
    ESegment *segment = data->GetSegment(i);
    if (!segment) continue;
    // THM (HOES) segments are not supported here: the adjoint differentiates
    // the T-matrix observable and the fast forward path computes the standard
    // cross section, both wrong for HOES. Bail so the LM caller falls back to
    // MIGRAD (which handles THM through Gradient()'s hybrid path).  A caller
    // that fills THM rows itself (AZUREAPI) asks for them to be left zero.
    if (segment->IsTHM()) {
      if (skipTHM) continue;
      return false;
    }
    const double norm = segment->GetNorm();
    const int normFull = segment->IsVaryNorm() ? pmap.NormIndex(i) : -1;
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
    if (bail) ok = false;
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
                           std::map<EPoint *, vector_r> &gradByPoint,
                           bool skipTHM) {
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
    if (skipTHM && segment->IsTHM()) continue;  // rows left zero for the caller
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
    if (bail) ok = false;
  }

  if (!ok) gradByPoint.clear();
  return ok;
}

// THM (HOES) rows: see the declaration in AZUREGrad.h.  Shared by the API
// (residual Jacobian, chi2 gradient, model gradients) and the CLI band.
std::vector<THMRows> ComputeTHMRows(CNuc *lc, EData *ld, const Config &config,
                                    const vector_r &full, const ParamIndexMap &pmap) {
  std::vector<THMRows> rows;
  const bool brune = (config.paramMask & Config::USE_BRUNE_FORMALISM);
  const int nCols = pmap.NumPacked();

  // Row bookkeeping mirrors ComputeResidualJacobian / EvaluateFilledChi2: one
  // row per non-null point, segments in order.
  int row = 0;
  for (int s = 1; s <= ld->NumSegments(); s++) {
    ESegment *seg = ld->GetSegment(s);
    if (!seg) continue;
    if (seg->IsTHM()) {
      THMRows tr;
      tr.segment = s;
      tr.firstRow = row;
      rows.push_back(tr);
    }
    for (int pid = 0; pid < seg->NumPoints(); pid++)
      if (seg->GetPoint(pid + 1)) row++;
  }
  if (rows.empty()) return rows;

  // Models of every THM point at the currently filled compound.  The points
  // are independent (as in every other forward pass, the compound is only
  // read), so they are evaluated in parallel: this runs 2 x (free E and gamma
  // columns) + 1 times per Jacobian, and serially it left all but one thread
  // idle -- 280 s for the 230 columns of the 12C+12C THM example against a
  // 0.67 s forward pass.  A mapped point is filled by its host's calculation,
  // so the mapped ones are read afterwards.
  std::vector<std::pair<int, int>> order;  // (row block t, point index pid)
  for (size_t t = 0; t < rows.size(); t++) {
    ESegment *seg = ld->GetSegment(rows[t].segment);
    for (int pid = 0; pid < seg->NumPoints(); pid++)
      if (seg->GetPoint(pid + 1)) order.push_back(std::make_pair((int)t, pid));
  }
  auto models = [&](std::vector<vector_r> &out) {
    out.assign(rows.size(), vector_r());
    std::vector<size_t> slot(order.size());
    for (size_t k = 0; k < order.size(); k++) {
      slot[k] = out[order[k].first].size();
      out[order[k].first].push_back(0.0);
    }
    const long n = (long)order.size();
#pragma omp parallel for schedule(dynamic)
    for (long k = 0; k < n; k++) {
      ESegment *seg = ld->GetSegment(rows[order[k].first].segment);
      if (seg->GetPoint(order[k].second + 1)->IsMapped()) continue;
      out[order[k].first][slot[k]] = seg->CalculateTheoreticalCrossSection(order[k].second, lc, config, ld);
    }
    for (long k = 0; k < n; k++) {
      ESegment *seg = ld->GetSegment(rows[order[k].first].segment);
      if (seg->GetPoint(order[k].second + 1)->IsMapped())
        out[order[k].first][slot[k]] = seg->CalculateTheoreticalCrossSection(order[k].second, lc, config, ld);
    }
  };
  auto fillCompound = [&](const vector_r &p) {
    lc->FillCompoundFromParams(p);
    if (brune) lc->CalcShiftFunctions(config);
  };

  // Norms and shifts are the caller's; the THM model depends on the R-matrix
  // parameters (and the shift already applied to the point energies).
  ld->FillNormsFromParams(full);
  ld->FillEnergyShiftsFromParams(full, ld, lc, &config);
  fillCompound(full);
  std::vector<vector_r> m0, mp, mm;
  models(m0);
  for (size_t t = 0; t < rows.size(); t++) {
    rows[t].m = m0[t];
    rows[t].Jm.assign(m0[t].size() * (size_t)nCols, 0.0);
    rows[t].J.assign(m0[t].size() * (size_t)nCols, 0.0);
  }

  // d m / d p by central differences, E and gamma columns only.
  for (int c = 0; c < nCols; c++) {
    const int f = pmap.PackedToFull(c);
    const ParamKind kind = pmap.Desc(f).kind;
    if (kind != ParamKind::LevelEnergy && kind != ParamKind::Gamma) continue;
    const double x0 = full[f];
    const double h = 1.0e-6 * (std::fabs(x0) + 1.0);
    vector_r pp = full;
    pp[f] = x0 + h;
    fillCompound(pp);
    models(mp);
    pp[f] = x0 - h;
    fillCompound(pp);
    models(mm);
    for (size_t t = 0; t < rows.size(); t++)
      for (size_t i = 0; i < m0[t].size(); i++)
        rows[t].Jm[i * nCols + c] = (mp[t][i] - mm[t][i]) / (2.0 * h);
  }
  fillCompound(full);

  // Residuals and the chain rule through the profiled scale.
  for (size_t t = 0; t < rows.size(); t++) {
    THMRows &tr = rows[t];
    ESegment *seg = ld->GetSegment(tr.segment);
    std::vector<double> d, e;
    for (int pid = 0; pid < seg->NumPoints(); pid++) {
      EPoint *pt = seg->GetPoint(pid + 1);
      if (!pt) continue;
      d.push_back(pt->GetCMCrossSection());
      e.push_back(pt->GetCMCrossSectionError());
      pt->SetFitCrossSection(tr.m[d.size() - 1]);
    }
    const size_t n = d.size();
    double Smm = 0.0, Smd = 0.0;
    for (size_t i = 0; i < n; i++) {
      if (e[i] == 0.0) continue;
      double w = 1.0 / (e[i] * e[i]);
      Smm += tr.m[i] * tr.m[i] * w;
      Smd += tr.m[i] * d[i] * w;
    }
    // s multiplies the model (s = 1/n, n the norm on the data).
    bool varyScale = false;
    double s;
    if (seg->IsProfiledNorm()) {
      seg->ProfileNormChiSquared();  // sets n* (or 1 if degenerate), as the CLI
      varyScale = (Smd > 0.0 && Smm > 0.0);
      s = varyScale ? Smd / Smm : 1.0;
    } else {
      double norm = seg->GetNorm();
      s = (norm != 0.0) ? 1.0 / norm : 0.0;
    }
    tr.r.assign(n, 0.0);
    for (size_t i = 0; i < n; i++)
      if (e[i] != 0.0) tr.r[i] = (s * tr.m[i] - d[i]) / e[i];

    vector_r ds(nCols, 0.0);
    if (varyScale) {
      for (int c = 0; c < nCols; c++) {
        double sdJ = 0.0, smJ = 0.0;
        for (size_t i = 0; i < n; i++) {
          if (e[i] == 0.0) continue;
          double w = 1.0 / (e[i] * e[i]);
          sdJ += d[i] * tr.Jm[i * nCols + c] * w;
          smJ += tr.m[i] * tr.Jm[i * nCols + c] * w;
        }
        ds[c] = (sdJ - 2.0 * s * smJ) / Smm;
      }
    }
    for (size_t i = 0; i < n; i++) {
      if (e[i] == 0.0) continue;
      for (int c = 0; c < nCols; c++)
        tr.J[i * nCols + c] = (s * tr.Jm[i * nCols + c] + tr.m[i] * ds[c]) / e[i];
    }
    tr.s = s;
    tr.profiled = varyScale;
    tr.ds = ds;
  }
  return rows;
}
