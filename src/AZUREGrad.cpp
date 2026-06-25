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
vector_matrix_r BuildShiftDerivTable(CNuc* compound, const Config& configure) {
  vector_matrix_r tab(compound->NumJGroups());
  for (int j = 1; j <= compound->NumJGroups(); j++) {
    JGroup* jg = compound->GetJGroup(j);
    tab[j-1].assign(jg->NumLevels(), std::vector<double>());
    if (!jg->IsInRMatrix()) continue;
    for (int la = 1; la <= jg->NumLevels(); la++) {
      tab[j-1][la-1].assign(jg->NumChannels(), 0.0);
      ALevel* level = jg->GetLevel(la);
      if (!level->IsInRMatrix()) continue;
      for (int ch = 1; ch <= jg->NumChannels(); ch++) {
        AChannel* channel = jg->GetChannel(ch);
        PPair* pair = compound->GetPair(channel->GetPairNum());
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
        tab[j-1][la-1][ch-1] = dS;
      }
    }
  }
  return tab;
}

/*!
 * Allocates and zeroes the physics-coordinate gradient accumulators.
 */
void GradAccum::Init(CNuc* compound) {
  const int nJ = compound->NumJGroups();
  e.assign(nJ, {});
  gamma.assign(nJ, {});
  for (int j = 1; j <= nJ; j++) {
    JGroup* jg = compound->GetJGroup(j);
    const int nLevels = jg->NumLevels();
    const int nCh = jg->NumChannels();
    e[j-1].assign(nLevels, 0.0);
    gamma[j-1].assign(nLevels, std::vector<double>(nCh, 0.0));
  }
}

/*!
 * Resets all accumulators to zero without reallocating.
 */
void GradAccum::Zero() {
  for (auto& row : e) std::fill(row.begin(), row.end(), 0.0);
  for (auto& lev : gamma)
    for (auto& row : lev) std::fill(row.begin(), row.end(), 0.0);
}

/*!
 * Adds another (same-dimension) accumulator into this one, term by term.
 */
void GradAccum::Add(const GradAccum& o) {
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
void GradAccum::Scatter(const ParamIndexMap& map, vector_r& gradFull) const {
  for (int j = 0; j < (int)e.size(); j++) {
    for (int la = 0; la < (int)e[j].size(); la++) {
      int ei = map.EIndex(j+1, la+1);
      if (ei >= 0 && ei < (int)gradFull.size()) gradFull[ei] += e[j][la];
      for (int ch = 0; ch < (int)gamma[j][la].size(); ch++) {
        int gi = map.GammaIndex(j+1, la+1, ch+1);
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
ParamIndexMap BuildParamIndexMap(CNuc* compound, EData* data,
                                 const std::vector<bool>& fixed) {
  ParamIndexMap map;
  int i = 0;

  // --- R-matrix parameters: mirror FillCompoundFromParams exactly. ---
  for (int j = 1; j <= compound->NumJGroups(); j++) {
    JGroup* jGroup = compound->GetJGroup(j);
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
      ESegment* segment = data->GetSegment(s);
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
      ESegment* segment = data->GetSegment(s);
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
bool GradTargetEffectAdjoint(EPoint* point, double fitBar, GradAccum& accum,
                             const vector_matrix_r* shiftDeriv, int xsComponent,
                             CNuc* compound, const Config& config) {
  int nsp = point->NumSubPoints();
  if(nsp <= 0) return true;

  std::vector<double> sigma(nsp);
  for(int i = 1; i <= nsp; i++) sigma[i-1] = point->GetSubPoint(i)->GetFitCrossSection();
  double model = point->GetFitCrossSection();

  std::vector<double> c(nsp, 0.0);
  for(int i = 1; i <= nsp; i++) {
    double s0 = sigma[i-1];
    double h = 1.0e-3 * (std::fabs(s0) + 1.0);
    point->GetSubPoint(i)->SetFitCrossSection(s0 + h);
    point->IntegrateTargetEffect(config);
    double yp = point->GetFitCrossSection();
    point->GetSubPoint(i)->SetFitCrossSection(s0 - h);
    point->IntegrateTargetEffect(config);
    double ym = point->GetFitCrossSection();
    point->GetSubPoint(i)->SetFitCrossSection(s0);
    c[i-1] = (yp - ym) / (2.0 * h);
  }
  point->SetFitCrossSection(model);

  bool ok = true;
  for(int i = 1; i <= nsp && ok; i++) {
    if(c[i-1] == 0.0) continue;
    EPoint* sp = point->GetSubPoint(i);
    AMatrixFunc amf(compound, config);
    amf.ClearMatrices();
    amf.FillMatrices(sp);
    amf.InvertMatrices();
    amf.CalculateTMatrix(sp);
    if(!amf.PointAdjoint(sp, fitBar * c[i-1], accum, shiftDeriv, xsComponent)) ok = false;
  }
  return ok;
}

// Forward + adjoint for one data point (target-effect or direct A-matrix).
bool GradAdjointPoint(EPoint* point, int xsComponent, double fitBar, GradAccum& accum,
                      const vector_matrix_r* shiftDeriv, CNuc* compound, const Config& config) {
  TargetEffect* te = point->IsTargetEffect()
    ? point->GetParentData()->GetTargetEffect(point->GetTargetEffectNum()) : nullptr;
  bool subpointTE = te && (te->IsConvolution() || te->IsTargetIntegration() ||
                           te->IsConvCoefficients());
  if(subpointTE)
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
bool GradOnePoint(ESegment* segment, EData* data, int i, int pointIdx,
                  CNuc* compound, const Config& config,
                  const vector_matrix_r* shiftDeriv, const FitBarFn& fitBarFn,
                  GradAccum& accum, double& outModel) {
  EPoint* point = segment->GetPoint(pointIdx + 1);
  outModel = 0.0;
  if(!point) return true;

  int baseComp = segment->GetCrossSectionComponent();
  bool hasComp = segment->HasComponents();

  TargetEffect* te = point->IsTargetEffect()
    ? point->GetParentData()->GetTargetEffect(point->GetTargetEffectNum()) : nullptr;
  bool subpointTE = te && (te->IsConvolution() || te->IsTargetIntegration() ||
                           te->IsConvCoefficients());

  // Fast path: plain point, single forward reused for model and adjoint.
  if(!hasComp && !subpointTE) {
    AMatrixFunc amf(compound, config);
    amf.ClearMatrices();
    amf.FillMatrices(point);
    amf.InvertMatrices();
    amf.CalculateTMatrix(point);
    amf.GenMatrixFunc::CalculateCrossSection(point);
    double model = (baseComp == 1) ? point->GetFitE1CrossSection()
                 : (baseComp == 2) ? point->GetFitE2CrossSection()
                                   : point->GetFitCrossSection();
    outModel = model;
    double fitBar = fitBarFn(segment, i, pointIdx, model);
    if(fitBar != 0.0)
      return amf.PointAdjoint(point, fitBar, accum, shiftDeriv, baseComp);
    return true;
  }

  // Combined model (components and/or target-effect integration).
  double model = segment->CalculateTheoreticalCrossSection(pointIdx, compound, config, data);
  outModel = model;
  double fitBar = fitBarFn(segment, i, pointIdx, model);
  if(fitBar == 0.0) return true;

  if(!hasComp) {
    return GradAdjointPoint(point, baseComp, fitBar, accum, shiftDeriv, compound, config);
  } else if(segment->GetOperationType() == SUM) {
    if(!GradAdjointPoint(point, baseComp, fitBar, accum, shiftDeriv, compound, config))
      return false;
    for(ESegment* comp : segment->GetComponentSegments()) {
      if(!comp || pointIdx >= comp->NumPoints()) continue;
      EPoint* cp = comp->GetPoint(pointIdx + 1);
      if(!cp) continue;
      double cv = cp->GetFitCrossSection();
      if(std::isnan(cv) || std::isinf(cv)) continue;
      if(!GradAdjointPoint(cp, comp->GetCrossSectionComponent(),
                           fitBar * comp->GetComponentScaling(), accum,
                           shiftDeriv, compound, config))
        return false;
    }
    return true;
  } else {  // RATIO
    const std::vector<ESegment*>& comps = segment->GetComponentSegments();
    if(comps.empty() || !comps[0] || pointIdx >= comps[0]->NumPoints()) return false;
    EPoint* denomPoint = comps[0]->GetPoint(pointIdx + 1);
    double baseVal = (baseComp == 1) ? point->GetFitE1CrossSection()
                   : (baseComp == 2) ? point->GetFitE2CrossSection()
                                     : point->GetFitCrossSection();
    double denomVal = denomPoint ? denomPoint->GetFitCrossSection() : 0.0;
    if(baseVal != 0.0 && denomVal != 0.0 && denomPoint) {
      if(!GradAdjointPoint(point, baseComp, fitBar * model / baseVal,
                           accum, shiftDeriv, compound, config))
        return false;
      if(!GradAdjointPoint(denomPoint, 0, -fitBar * model / denomVal,
                           accum, shiftDeriv, compound, config))
        return false;
    }
    return true;
  }
}

}  // namespace

bool AccumulateEGammaGradient(CNuc* compound, EData* data, const Config& config,
                              const ParamIndexMap& pmap,
                              const vector_matrix_r* shiftDeriv,
                              const FitBarFn& fitBarFn,
                              GradAccum& accum) {
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
  for(GradAccum& ta : threadAccum) ta.Init(compound);

  bool ok = true;
  for(int i = 1; i <= data->NumSegments() && ok; i++) {
    ESegment* segment = data->GetSegment(i);
    if(!segment) continue;
    const int nPoints = segment->NumPoints();
    bool bail = false;
#pragma omp parallel for shared(bail)
    for(int pointIdx = 0; pointIdx < nPoints; pointIdx++) {
      if(bail) continue;
#ifdef _OPENMP
      GradAccum& local = threadAccum[omp_get_thread_num()];
#else
      GradAccum& local = threadAccum[0];
#endif
      double model;
      if(!GradOnePoint(segment, data, i, pointIdx, compound, config, shiftDeriv,
                       fitBarFn, local, model)) {
#pragma omp critical
        bail = true;
      }
    }
    if(bail) { ok = false; break; }
  }

  // Reduce the per-thread partials only on success; on bail the caller ignores
  // `accum` and falls back to finite differences, so leave it untouched.
  if(ok)
    for(const GradAccum& ta : threadAccum) accum.Add(ta);
  return ok;
}

bool ComputeResidualJacobian(CNuc* compound, EData* data, const Config& config,
                             const ParamIndexMap& pmap,
                             const vector_matrix_r* shiftDeriv,
                             vector_r& residuals, vector_r& jacobian, int& nCols) {
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
  for(int i = 1; i <= nSeg; i++) {
    ESegment* segment = data->GetSegment(i);
    if(!segment) continue;
    const int nPoints = segment->NumPoints();
    rowOf[i].assign(nPoints, -1);
    for(int pointIdx = 0; pointIdx < nPoints; pointIdx++)
      if(segment->GetPoint(pointIdx + 1)) rowOf[i][pointIdx] = totalRows++;
  }
  residuals.assign(totalRows, 0.0);
  jacobian.assign((size_t)totalRows * nCols, 0.0);

  bool ok = true;
  for(int i = 1; i <= nSeg && ok; i++) {
    ESegment* segment = data->GetSegment(i);
    if(!segment) continue;
    const double norm = segment->GetNorm();
    const int normFull = segment->IsVaryNorm() ? pmap.NormIndex(i) : -1;
    const int nPoints = segment->NumPoints();
    bool bail = false;

#pragma omp parallel
    {
      GradAccum accum;                       // per-thread, reused across its rows
      accum.Init(compound);
      vector_r fullRow(pmap.NumFull(), 0.0);
#pragma omp for
      for(int pointIdx = 0; pointIdx < nPoints; pointIdx++) {
        if(bail) continue;
        const int row = rowOf[i][pointIdx];
        if(row < 0) continue;
        EPoint* pt = segment->GetPoint(pointIdx + 1);
        double cmErr = pt->GetCMCrossSectionError();
        double dataval = pt->GetCMCrossSection();

        // Residual r = (model - data*n)/(cmErr*n);  d r / d model = 1/(cmErr*n).
        // PointAdjoint with this cotangent yields d r / d(E,gamma) directly.
        double denom = cmErr * norm;
        FitBarFn resFitBar = [denom](ESegment*, int, int, double) -> double {
          return (denom != 0.0) ? 1.0 / denom : 0.0;
        };

        accum.Zero();
        double model = 0.0;
        if(!GradOnePoint(segment, data, i, pointIdx, compound, config, shiftDeriv,
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
        if(normFull >= 0 && cmErr != 0.0 && norm != 0.0)
          fullRow[normFull] += -model / (cmErr * norm * norm);

        const size_t base = (size_t)row * nCols;
        for(int f = 0; f < pmap.NumFull(); f++) {
          int packed = pmap.FullToPacked(f);
          if(packed >= 0 && packed < nCols) jacobian[base + packed] = fullRow[f];
        }
      }
    }
    if(bail) ok = false;
  }

  if(!ok) { residuals.clear(); jacobian.clear(); }
  return ok;
}
