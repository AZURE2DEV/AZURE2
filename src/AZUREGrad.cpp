#include "AZUREGrad.h"
#include "RunDiagnostics.h"
#include "azure2/support/Environment.h"
#include "GradEnergyTasks.h"
#include "ModelPass.h"
#include "Profile.h"
#include <atomic>
#include <chrono>
#include <iostream>
#include <cstdlib>
#include "CNuc.h"
#include "CalculationPlan.h"
#include "WorkingCopy.h"
#include "JGroup.h"
#include "ALevel.h"
#include "AChannel.h"
#include "PPair.h"
#include "Segment.h"
#include "Unit.h"
#include "TargetEffect.h"
#include "AMatrixFunc.h"
#include "Config.h"
#include "ChannelFunc.h"
#include "CoulFunc.h"
#include "ShftFunc.h"
#include "THMMatrixFunc.h"
#include "PointFold.h"
#include <cstring>
#include <iostream>
#include <cmath>
#include <stdexcept>
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
    const JGroup *jg = compound->GetJGroup(j);
    tab[j - 1].assign(jg->NumLevels(), std::vector<double>());
    if (!jg->IsInRMatrix()) continue;
    for (int la = 1; la <= jg->NumLevels(); la++) {
      tab[j - 1][la - 1].assign(jg->NumChannels(), 0.0);
      const ALevel *level = jg->GetLevel(la);
      if (!level->IsInRMatrix()) continue;
      for (int ch = 1; ch <= jg->NumChannels(); ch++) {
        const AChannel *channel = jg->GetChannel(ch);
        const PPair *pair = compound->GetPair(channel->GetPairNum());
        // The Brune boundary/shift block applies only to particle 'P' channels.
        if (pair->GetPType() != 0 || channel->GetRadType() != 'P') continue;
        if (park && level->ShiftFunctionsValidAt(level->GetFitE())) {
          tab[j - 1][la - 1][ch - 1] = level->GetShiftDerivative(ch);
          continue;
        }
        int l = channel->GetL();
        double levelEnergy = level->GetFitE();
        double resonanceEnergy = levelEnergy - (pair->GetSepE() + pair->GetExE());
        double dS = ChannelFunc(pair, !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC))
                        .ShiftDerivative(l, pair->SpinOrbitTwoJ(jg->GetJ()), resonanceEnergy);
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
    const JGroup *jg = compound->GetJGroup(j);
    if (!jg->IsInRMatrix()) continue;
    for (int la = 1; la <= jg->NumLevels(); la++) {
      const ALevel *level = jg->GetLevel(la);
      if (!level->IsInRMatrix() || level->GetParkNorm() >= 0.0) continue;
      const double factor = 2.0 * level->GetParkNorm() / s2;
      double dE = 0.0;
      for (int ch = 1; ch <= jg->NumChannels(); ch++) {
        const AChannel *channel = jg->GetChannel(ch);
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
 * ((sum w) - c)^2 / s^2 (c its centre, 1 unless written), which couples every abundance through the sum.
 */
void AddAbundancePenaltyGradient(CNuc *compound, const Config &configure, GradAccum &accum) {
  if (!compound->HasAbundanceBlock()) return;
  for (int k = 1; k <= compound->NumAbundanceParams(); k++) {
    int sample, nucleus;
    compound->AbundanceOf(k, sample, nucleus);
    const double s = compound->SampleSumSigma(sample, configure.abundanceSumSigma);
    const double centre = compound->GetSample(sample)->sumCentre;
    double g = (s > 0.0) ? 2.0 * (compound->AbundanceSum(sample) - centre) / (s * s) : 0.0;
    const SampleAbundance &e = compound->GetSample(sample)->entries[nucleus - 1];
    const double w = compound->AbundanceValues()[k - 1];
    if (e.sigma > 0.0 && !e.priorReplaced) g += 2.0 * (w - e.nominal) / (e.sigma * e.sigma);
    accum.AddAbundance(k, g);
  }
}

/*!
 * Allocates and zeroes the physics-coordinate gradient accumulators.
 */
void GradAccum::Init(CNuc *compound) {
  source = compound;
  const int nJ = compound->NumJGroups();
  e.assign(nJ, {});
  gamma.assign(nJ, {});
  for (int j = 1; j <= nJ; j++) {
    const JGroup *jg = compound->GetJGroup(j);
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
  doppler.clear();
  thickness.clear();
}

/*!
 * Resets all accumulators to zero without reallocating.
 */
void GradAccum::Zero() {
  direct.clear();
  for (auto &row : e) std::fill(row.begin(), row.end(), 0.0);
  for (auto &lev : gamma)
    for (auto &row : lev) std::fill(row.begin(), row.end(), 0.0);
  std::fill(abundance.begin(), abundance.end(), 0.0);
  std::fill(screening.begin(), screening.end(), 0.0);
  std::fill(doppler.begin(), doppler.end(), 0.0);
  std::fill(thickness.begin(), thickness.end(), 0.0);
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
  direct.insert(direct.end(), o.direct.begin(), o.direct.end());
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
  if (o.doppler.size() > doppler.size()) doppler.resize(o.doppler.size(), 0.0);
  for (int k = 0; k < (int)o.doppler.size(); k++) doppler[k] += o.doppler[k];
  if (o.thickness.size() > thickness.size()) thickness.resize(o.thickness.size(), 0.0);
  for (int k = 0; k < (int)o.thickness.size(); k++) thickness[k] += o.thickness[k];
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
  // The photon amplitudes of a gamma_decay level are functions of the record's numbers, which hold their slots: their
  // cotangents go through the map's Jacobian (GammaDecayMap.h), not to the slot.
  const GammaDecaySet *decays = (source && !source->GetGammaDecays().Empty()) ? &source->GetGammaDecays() : nullptr;
  for (int j = 0; j < (int)e.size(); j++) {
    for (int la = 0; la < (int)e[j].size(); la++) {
      int ei = map.EIndexFlat(j + 1, la + 1);
      if (ei >= 0 && ei < (int)gradFull.size()) gradFull[ei] += e[j][la];
      for (int ch = 0; ch < (int)gamma[j][la].size(); ch++) {
        if (decays && decays->OwnsSlot(j + 1, la + 1, ch + 1)) continue;
        int gi = map.GammaIndexFlat(j + 1, la + 1, ch + 1);
        if (gi >= 0 && gi < (int)gradFull.size()) gradFull[gi] += gamma[j][la][ch];
      }
    }
  }
  if (decays) decays->PullBack(*source, *this, map, gradFull);
  for (const auto &d : direct)
    if (d.first >= 0 && d.first < (int)gradFull.size()) gradFull[d.first] += d.second;
  for (int k = 0; k < (int)abundance.size(); k++) {
    int ai = map.AbundanceIndex(k + 1);
    if (ai >= 0 && ai < (int)gradFull.size()) gradFull[ai] += abundance[k];
  }
  for (int k = 0; k < (int)screening.size(); k++) {
    int si = map.ScreeningIndex(k + 1);
    if (si >= 0 && si < (int)gradFull.size()) gradFull[si] += screening[k];
  }
  for (int k = 0; k < (int)doppler.size(); k++) {
    int di = map.DopplerIndex(k + 1);
    if (di >= 0 && di < (int)gradFull.size()) gradFull[di] += doppler[k];
  }
  for (int k = 0; k < (int)thickness.size(); k++) {
    int ti = map.ThicknessIndex(k + 1);
    if (ti >= 0 && ti < (int)gradFull.size()) gradFull[ti] += thickness[k];
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
 *   - WorkingCopy::FillNormsFromParams        (norms, only for IsVaryNorm segments)
 *   - WorkingCopy::FillEnergyShiftsFromParams (one energy shift per segment, then
 *                                        one sqrt(E) coefficient per segment,
 *                                        then the shift-function parameters of
 *                                        each segment that declares them)
 *
 * Keeping this in lock-step with those routines is what guarantees the returned
 * gradient vector is a permutation-correct match for the sampler's parameter
 * ordering (see PLAN.md, "Parameter-vector order must exactly match").
 */
ParamIndexMap BuildParamIndexMap(CNuc *compound, WorkingCopy *data,
                                 const std::vector<bool> &fixed) {
  ParamIndexMap map;
  int i = 0;

  // --- R-matrix parameters: mirror FillCompoundFromParams exactly. ---
  map.eFlat_.assign(compound->NumJGroups(), {});
  map.gammaFlat_.assign(compound->NumJGroups(), {});
  for (int j = 1; j <= compound->NumJGroups(); j++) {
    const JGroup *jGroup = compound->GetJGroup(j);
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
        const bool radiative = compound->GetGammaDecays().OwnsSlot(j, la, ch);
        map.desc_.push_back(ParamDesc{radiative ? ParamKind::RadiativeDecay : ParamKind::Gamma, j, la, ch, -1, false});
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
    for (int s = 1; s <= data->Data().NumSegments(); s++) {
      const Segment *segment = data->Data().GetSegment(s);
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
    for (int s = 1; s <= data->Data().NumSegments(); s++) {
      const Segment *segment = data->Data().GetSegment(s);
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
    for (int s = 1; s <= data->Data().NumSegments(); s++) {
      const Segment *segment = data->Data().GetSegment(s);
      if (segment) {
        map.sqrtShiftIndex_[s] = i;
        map.desc_.push_back(ParamDesc{ParamKind::EnergyShiftSqrt, -1, -1, -1, s, false});
        i++;
      }
    }
  }

  // --- Shift-function parameters: those of each segment that declares one (WorkingCopy::FillMnParams). ---
  map.energyShiftFnOffset_ = i;
  if (data) {
    for (int s = 1; s <= data->Data().NumSegments(); s++) {
      const Segment *segment = data->Data().GetSegment(s);
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

  // --- Electron-screening potentials: one per target with a screening field (WorkingCopy::FillMnParams). ---
  if (data) {
    for (int k = 1; k <= data->Data().NumScreeningParams(); k++) {
      map.screeningIndex_[k] = i;
      ParamDesc d{ParamKind::Screening, -1, -1, -1, -1, false};
      d.screening = k;
      map.desc_.push_back(d);
      i++;
    }
  }

  // --- Fitted target densities: one per densityfit line or density group (WorkingCopy::FillMnParams). ---
  // No analytic derivative: the target thickness enters through the depth window, finite differences.
  if (data) {
    for (int k = 1; k <= data->Data().NumDensityParams(); k++) {
      map.densityIndex_[k] = i;
      ParamDesc d{ParamKind::Density, -1, -1, -1, -1, false};
      d.density = k;
      map.desc_.push_back(d);
      i++;
    }
  }

  // --- Fitted Doppler temperatures: one per 'doppler T vary sigma' line (WorkingCopy::FillMnParams). ---
  // Analytic: d/dT through the Doppler weights (Unit::DopplerWeights), GradAccum::AddDoppler.
  if (data) {
    for (int k = 1; k <= data->Data().NumDopplerParams(); k++) {
      map.dopplerIndex_[k] = i;
      ParamDesc d{ParamKind::Doppler, -1, -1, -1, -1, false};
      d.doppler = k;
      map.desc_.push_back(d);
      i++;
    }
  }

  // --- Fitted sample thicknesses: one per transmission / yield line with the fit form (WorkingCopy::FillMnParams). ---
  // Analytic: d/dn of the transform (the tree's backward, AddParamCotangents), GradAccum::AddThickness.
  if (data) {
    for (int k = 1; k <= data->Data().NumThicknessParams(); k++) {
      map.thicknessIndex_[k] = i;
      ParamDesc d{ParamKind::Thickness, -1, -1, -1, data->Data().GetThicknessSegment(k), false};
      d.thickness = k;
      map.desc_.push_back(d);
      i++;
    }
  }

  // --- THM coherent backgrounds (cbackground=): the last block, if any. ---
  if (data)
    for (int k = 0; k < data->Data().NumThmCoherentParams(); k++) {
      map.desc_.push_back(ParamDesc{ParamKind::ThmCoherent, -1, -1, k + 1, -1, false});
      i++;
    }
  AppendApparatusParams(map, data, i);  // the apparatus parameters, last (ApparatusParameters.cpp)
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

// The cotangents of the sub-points of one target-effect (sub-point integrated)
// point.  The integrator IntegrateTargetEffect is a fixed linear functional of
// the values in the sub-points' cross-section slots, Y = sum_i c_i s_i, so the
// weights c_i = dY/ds_i come from one pass of the integrator on a linear form
// (or central differences on it); each sub-point's cotangent is then fitBar c_i,
// back-propagated with the other units at its energy (GradEnergyTasks.cpp).
//
// A cross section (or E1/E2 component) is Y itself.  A target-averaged
// analyzing power is the ratio of two such functionals with the SAME weights
// (Unit::FoldSubPoints runs the one integrator on s_i
// and on s_i A_i, s_i being the sub-point's differential cross section at the
// point's angle):
//     <A> = N / D,   N = sum_i c_i s_i A_i,   D = sum_i c_i s_i,
//     d<A> = sum_i c_i [ s_i dA_i + (A_i - <A>) ds_i ] / D,
// so each sub-point contributes an A_y cotangent fitBar c_i s_i / D and a
// differential-cross-section cotangent fitBar c_i (A_i - <A>) / D.
// AZURE_GRAD_TIMING: thread-summed time in the target-effect combiner weights.

bool GradTargetEffectAdjoint(WorkingCopy &data, Unit *point, azure2::Results &results, double fitBar, GradAccum &accum,
                             int xsComponent, CNuc *compound, const Config &config, int row,
                             std::vector<azure2::UnitCotangent> &out) {
  int nsp = point->NumSubPoints();
  if (nsp <= 0) return true;

  // A ratio (Unit::FoldsRatio: A_y over its cross section, a normalised coefficient a_k over A_0) is
  //   <r> = N / D,  N = sum c_i w_i r_i,  D = sum c_i w_i,   d<r> = sum c_i [ w_i dr_i + (r_i - <r>) dw_i ] / D:
  // an element cotangent fitBar c_i w_i / D (PointTBar on the observable) and a weight cotangent
  // fitBar c_i (r_i - <r>) / D on the weight's route (Unit::RatioWeightRoute).  Anything else is linear in the slots.
  const bool ratio = point->FoldsRatio();
  // An E1 or E2 part exists for sigma, sfactor and modified_sfactor only (the field multipole,
  // experiment/src/Quantities.cpp): never for a ratio (a coefficient, an analyzing power).
  if (ratio && xsComponent != 0)
    throw std::logic_error("GradTargetEffectAdjoint: an E1 or E2 part of a coefficient or an analyzing power "
                           "(an invariant of the set-up is broken)");
  if (point->IsAnalyzingPower()) {
    // Supported: particle-exit A_y, full observable (every point's sub-points are units of the point pass, with their
    // own records).  Capture A_y keeps its refusal: the attenuation folded into P_k has no counterpart for the P_k^1
    // of the numerator.
    const int exitPairNum = compound->GetPairNumFromKey(point->GetExitKey());
    if (compound->GetPair(exitPairNum)->GetPType() != 0) return false;
  }

  std::vector<double> sigma(nsp);
  for (int i = 1; i <= nsp; i++) sigma[i - 1] = point->SlotValue(*point->GetSubPoint(data, i), i, results);
  double model = point->GetFitCrossSection(results);

  const auto tc0 = std::chrono::steady_clock::now();
  // Exact weights in one pass of the integrator on a linear form
  // (Unit::TargetEffectWeights).  The central differences used before cost
  // 2 nsp integrator calls, each O(nsp) or O(nsp^2) (straggling, convolution
  // + target integration): O(nsp^3) per point, which made one Jacobian of the
  // 11B+a base grid (CASPAR thick targets, ~2200 sub-points per point) take
  // longer than 30 min.  They remain as the fallback and as a check
  // (AZURE_TE_WEIGHTS_CHECK=N compares the first N points).
  // The slots' records are read, never written: other points of the group read them too (SlotGroup.h).
  auto fdWeights = [&](std::vector<double> &w) {
    w.assign(nsp, 0.0);
    std::vector<double> v(sigma);
    for (int i = 1; i <= nsp; i++) {
      double s0 = sigma[i - 1];
      double h = 1.0e-3 * (std::fabs(s0) + 1.0);
      v[i - 1] = s0 + h;
      double yp = FoldValues(*point, data, config, v);
      v[i - 1] = s0 - h;
      double ym = FoldValues(*point, data, config, v);
      v[i - 1] = s0;
      w[i - 1] = (yp - ym) / (2.0 * h);
    }
  };
  std::vector<double> c;
  if (config.environment.foldWeightsByDifferences || !FoldWeights(*point, data, config, c)) {
    fdWeights(c);
  } else if (azure2::support::ProcessEnvironment().Set(azure2::support::Switch::kFoldWeightsCheck)) {
    if (config.diagnostics->foldWeightChecks++ < azure2::support::ProcessEnvironment().Integer(azure2::support::Switch::kFoldWeightsCheck, 0)) {
      std::vector<double> cfd;
      fdWeights(cfd);
      double cmax = 0.0, dmax = 0.0, yW = 0.0;
      for (int i = 0; i < nsp; i++) {
        cmax = std::max(cmax, std::fabs(cfd[i]));
        dmax = std::max(dmax, std::fabs(c[i] - cfd[i]));
        yW += c[i] * sigma[i];
      }
      const double yI = FoldValues(*point, data, config, sigma);
#pragma omp critical
      std::cerr << "[te-weights] entrance " << point->GetEntranceKey() << " exit " << point->GetExitKey() << " E " << point->GetCMEnergy() << " nsp " << nsp
                << ": max|c - c_fd|/max|c_fd| " << (cmax > 0 ? dmax / cmax : dmax) << ", sum c*s " << yW
                << " vs integrator " << yI << std::endl;
    }
  }
  // For a ratio: its weights, elements, the denominator D and the average <r>, from the same integrator calls the
  // forward pass makes.
  double den = 0.0, avg = 0.0;
  std::vector<double> weight, element, product;
  if (ratio) {
    weight.resize(nsp), element.resize(nsp), product.resize(nsp);
    for (int i = 1; i <= nsp; i++) {
      point->RatioPartsAt(data, i, results, weight[i - 1], element[i - 1]);
      product[i - 1] = weight[i - 1] * element[i - 1];
    }
    den = FoldValues(*point, data, config, weight);
    const double num = FoldValues(*point, data, config, product);
    // The forward pass returns 0 for D = 0, a constant: no gradient.
    if (!(std::fabs(den) > 0.0)) {
      point->SetFitCrossSection(results, model);
      return true;
    }
    avg = num / den;
  }
  point->SetFitCrossSection(results, model);
  // A fitted Doppler temperature (doppler_<n>): d/dT through the Doppler weights, the sub-point values fixed.  The
  // ratios (A_y, a normalised coefficient) take d(N/D) = (dN - <A> dD)/D on the same integrator.
  if (const int dk = data.Data().DopplerParamOfEffect(point->GetTargetEffectNum())) {
    std::vector<double> vals(nsp);
    double dT = 0.0;
    if (ratio) {
      const double dN = FoldTemperatureDerivative(*point, data, config, product);
      const double dD = FoldTemperatureDerivative(*point, data, config, weight);
      dT = (dN - avg * dD) / den;
    } else {
      for (int i = 1; i <= nsp; i++) {
        Site *sp = point->GetSubPoint(data, i);
        vals[i - 1] = (xsComponent == 1) ? sp->GetFitE1CrossSection(results)
                      : (xsComponent == 2) ? sp->GetFitE2CrossSection(results) : sigma[i - 1];
      }
      dT = point->LabFoldFactor() * FoldTemperatureDerivative(*point, data, config, vals);
    }
    accum.AddDoppler(dk, fitBar * dT);
  }
  const auto tc1 = std::chrono::steady_clock::now();
  config.diagnostics->foldWeightNs += std::chrono::duration_cast<std::chrono::nanoseconds>(tc1 - tc0).count();

  auto emit = [&](int i, double value, int component, bool ayCrossSection) {
    if (value != 0.0)
      out.push_back({data.SlotSite(*point, i), row, value, component, ayCrossSection, point->PlanTask(), point, i});
  };
  // Read through a lab-frame row set (RowSet.h) a slot's value is its observable over J_i and the fold is multiplied
  // by the point's own Jacobian: the cotangent of the observable is the slot's times LabFoldFactor / J_i (both 1
  // otherwise); a ratio's element (A_y) is not divided, its weight (the cross section) is.
  const double labFactor = point->LabFoldFactor();
  int weightComponent = 0;
  bool weightAy = false;
  if (ratio) point->RatioWeightRoute(weightComponent, weightAy);
  for (int i = 1; i <= nsp; i++) {
    if (c[i - 1] == 0.0) continue;
    if (!ratio) {
      emit(i, fitBar * c[i - 1] * labFactor / point->SlotJacobian(i), xsComponent, false);
    } else {
      emit(i, fitBar * c[i - 1] * weight[i - 1] / den, 0, false);
      emit(i, fitBar * c[i - 1] * (element[i - 1] - avg) / den / point->SlotJacobian(i), weightComponent, weightAy);
    }
  }
  return true;
}

// The parameters' share of a point's model cotangents (ModelPass.h): each weighted term's abundance, a transform's
// sample thickness (segment_<key>_thickness, its group leader's for a follower) and a fitted Doppler temperature.
void AddParamCotangents(const azure2::ModelCotangents &c, const Segment *segment, Unit *base, const WorkingCopy &data, int i,
                        CNuc *compound, GradAccum &accum) {
  using azure2::experiment::ParamRef;
  for (std::size_t q = 0; q < c.refs->size(); q++) {
    if (c.params[q] == 0.0) continue;
    const ParamRef &r = (*c.refs)[q];
    if (r.kind == ParamRef::Kind::Abundance) {
      accum.AddAbundance(compound->AbundanceOrdinal(segment->GetNaturalSample(), r.nucleus), c.params[q]);
    } else {
      const int k = data.Data().ThicknessParamOfSegment(i);
      if (k > 0) accum.AddThickness(k, c.params[q]);
    }
  }
  if (c.params.back() != 0.0) {
    const int dk = data.Data().DopplerParamOfEffect(base->GetTargetEffectNum());
    if (dk > 0) accum.AddDoppler(dk, c.params.back());
  }
}

}  // namespace

bool azure2::PointCotangents(WorkingCopy &data, const GradPoint &p, CNuc *compound, const Config &config,
                             azure2::Results &results, double fitBar, int row, GradAccum &params,
                             std::vector<UnitCotangent> &out) {
  if (fitBar == 0.0) return true;
  // The cotangents of the point's expression tree: its leaves (points, or sub-points of a transform's terms) and its
  // parameters.
  ModelCotangents c;
  if (!ModelPointBackward(data, *p.segment, p.pointIdx, compound, config, results, fitBar, c)) return false;
  AddParamCotangents(c, p.segment, p.point, data, p.segNum, compound, params);
  auto emit = [&](const Unit *reader, const Site *u, int slot, double value, int component) {
    if (value != 0.0) out.push_back({data.SiteIndex(u), row, value, component, false, reader->PlanTask(), reader, slot});
  };
  if (c.subPoints) {
    for (const ModelCotangents::Leaf &l : c.leaves)
      for (std::size_t k = 0; k < l.sub.size(); k++) emit(l.point, l.point->GetSubPoint(data, (int)k + 1), (int)k + 1, l.sub[k], l.xsComponent);
    return true;
  }
  for (const ModelCotangents::Leaf &l : c.leaves) {
    if (l.cotangent == 0.0) continue;
    if (!HasSubPointFold(data, *l.point))
      emit(l.point, l.point, 0, l.cotangent, l.xsComponent);
    else if (!GradTargetEffectAdjoint(data, l.point, results, l.cotangent, params, l.xsComponent, compound, config, row,
                                      out))
      return false;
  }
  return true;
}

namespace {

// ---------------------------------------------------------------------------
// THM (HOES) model derivative of one point -- a point or a mapped point, with
// its folding if it has one -- by the HOES adjoint (THMMatrixFunc::
// THMPointAdjoint): d m/d(E, gamma) into `accum`, d m/d(coherent value k)
// into coherentBar[k].  A folded point is linear in its sub-point cross
// sections; the weights come from the combiner as in GradTargetEffectAdjoint.
// A folded analyzing power (a ratio of two folds): not covered (false).
bool ThmPointModelGradient(WorkingCopy &data, Unit *point, azure2::Results &results, CNuc *compound, const Config &config,
                           const ThmAdjointContext &ctx, THMMatrixFunc &f,
                           GradAccum &accum, std::vector<double> &coherentBar) {
  f.Reset(compound, config);
  auto adjoint = [&](Site *pt, double weight) -> bool {
    f.ClearMatrices();
    f.FillMatrices(pt);
    f.InvertMatrices();
    return f.THMPointAdjoint(pt, weight, accum, ctx, &coherentBar);
  };
  const TargetEffect *te = point->IsTargetEffect() ? data.Data().GetTargetEffect(point->GetTargetEffectNum()) : nullptr;
  if (!(te && te->IsSubPointEffect())) return adjoint(point, 1.0);
  if (point->IsAnalyzingPower()) return false;
  const int nsp = point->NumSubPoints();
  if (nsp <= 0) return true;
  std::vector<double> sigma(nsp), c(nsp, 0.0);
  for (int i = 1; i <= nsp; i++) sigma[i - 1] = point->GetSubPoint(data, i)->GetFitCrossSection(results);
  std::vector<double> v(sigma);  // the slots' records are read, never written (SlotGroup.h)
  for (int i = 1; i <= nsp; i++) {
    const double s0 = sigma[i - 1], h = 1.0e-3 * (std::fabs(s0) + 1.0);
    v[i - 1] = s0 + h;
    const double yp = FoldValues(*point, data, config, v);
    v[i - 1] = s0 - h;
    const double ym = FoldValues(*point, data, config, v);
    v[i - 1] = s0;
    c[i - 1] = (yp - ym) / (2.0 * h);
  }
  for (int i = 1; i <= nsp; i++)
    if (c[i - 1] != 0.0 && !adjoint(point->GetSubPoint(data, i), c[i - 1])) return false;
  return true;
}

void ReportBail(WorkingCopy &data, const char *what, const Segment *segment, int i) {
  if (!azure2::support::ProcessEnvironment().Set(azure2::support::Switch::kGradDebug) || !segment) return;
  Unit *pt = segment->GetPoint(data, 1);
  std::cerr << "[grad] " << what << " adjoint bailed on segment " << i << " (key " << segment->GetSegmentKey()
            << ", entrance " << segment->GetEntranceKey() << ", exit " << segment->GetExitKey()
            << ", components " << segment->HasComponents() << ", angular-dist " << (pt && pt->IsAngularDist())
            << ", A_y " << (pt && pt->IsAnalyzingPower()) << ", target effect " << (pt && pt->IsTargetEffect())
            << ")" << std::endl;
}

// The formal-basis parts that the adjoint does not differentiate: the Brune amplitude of an external-capture final
// state that shares its block (CNuc::ECFinalStatesFromFormal, an iteration on the block's formal parameters) and the
// --ignore-externals projector of a channel capture (CNuc::FormalCaptureMasksOutsideAdjoint, external_capture.md
// section 7).  Every adjoint entry point asks this and falls back to finite differences.
bool FormalPartsOutsideAdjoint(CNuc *compound, const ParamIndexMap &pmap) {
  if (compound->HasFormalECFinalBlocks()) return true;
  return compound->FormalCaptureMasksOutsideAdjoint([&pmap](int j, int la, int ch) {
    const int full = pmap.GammaIndexFlat(j, la, ch);
    return full >= 0 && full < pmap.NumFull() && pmap.FullToPacked(full) >= 0;
  });
}

}  // namespace

/// AZURE2_THM_FD=1 (any value but empty or 0): THM rows by central
/// differences, the path before the HOES adjoint (kept for testing).
bool ThmFiniteDifferencesRequested(const Config &config) { return config.environment.thmByDifferences; }


bool AccumulateEGammaGradient(CNuc *compound, WorkingCopy *data, azure2::Results &results, const Config &config,
                              const ParamIndexMap &pmap,
                              const vector_matrix_r *shiftDeriv,
                              const FitBarFn &fitBarFn,
                              GradAccum &accum, const ThmAdjointContext *hoes) {
  // THM (HOES) segments are computed by THMMatrixFunc, not the T-matrix
  // observable this adjoint differentiates; they are left out here. Their
  // contribution to the gradient is added in AZURECalc::Gradient by the HOES
  // adjoint or finite differences of a chi^2 restricted to THM segments, so
  // mixed THM+direct datasets keep the analytic speed-up on the direct part.
  if (FormalPartsOutsideAdjoint(compound, pmap)) return false;
  int bailSeg = 0;
  if (!azure2::EnergyTaskGradient(compound, data, results, config, shiftDeriv, fitBarFn, accum, bailSeg, hoes)) {
    if (azure2::support::ProcessEnvironment().Set(azure2::support::Switch::kGradDebug)) {
      if (config.FirstTime("AZURE_GRAD_DEBUG energy-task gradient declined")) {
        bool inactiveLevel = false;
        for (int j = 1; j <= compound->NumJGroups(); j++) {
          const JGroup *jg = compound->GetJGroup(j);
          if (!jg->IsInRMatrix()) continue;
          for (int la = 1; la <= jg->NumLevels(); la++)
            if (!jg->GetLevel(la)->IsInRMatrix()) inactiveLevel = true;
        }
        const Segment *segment = data->Data().GetSegment(bailSeg);
        Unit *pt = segment ? segment->GetPoint(*data, 1) : nullptr;
        std::cerr << "[grad] adjoint bailed on segment " << bailSeg << " (entrance "
                  << (segment ? segment->GetEntranceKey() : 0) << ", exit " << (segment ? segment->GetExitKey() : 0)
                  << ", components " << (segment && segment->HasComponents()) << ", angular-dist "
                  << (pt && pt->IsAngularDist()) << ", A_y " << (pt && pt->IsAnalyzingPower()) << ", target effect "
                  << (pt && pt->IsTargetEffect()) << "); a level outside the R matrix in an "
                  << "R-matrix J group: " << (inactiveLevel ? "YES" : "no") << std::endl;
      }
    }
    // On bail the caller ignores `accum` and falls back to finite differences.
    return false;
  }
  if (accum.formalUsed) {
    AMatrixFunc amf(compound, config);
    amf.ResolveFormalCotangents(accum, shiftDeriv);
  }
  return true;
}

namespace {

// The free level parameters of the packed vector, per J-group (FormalBasis).
AMatrixFunc::FreeMask BuildFreeMask(CNuc *compound, const ParamIndexMap &pmap) {
  AMatrixFunc::FreeMask mask;
  const int nJ = compound->NumJGroups();
  mask.e.assign(nJ, {});
  mask.g.assign(nJ, {});
  auto isFree = [&pmap](int full) { return full >= 0 && full < pmap.NumFull() && pmap.FullToPacked(full) >= 0; };
  for (int j = 1; j <= nJ; j++) {
    const JGroup *jg = compound->GetJGroup(j);
    const int nl = jg->NumLevels(), nch = jg->NumChannels();
    mask.e[j - 1].assign(nl, 0);
    mask.g[j - 1].assign(nl, std::vector<char>(nch, 0));
    for (int la = 1; la <= nl; la++) {
      mask.e[j - 1][la - 1] = isFree(pmap.EIndexFlat(j, la));
      for (int ch = 1; ch <= nch; ch++) mask.g[j - 1][la - 1][ch - 1] = isFree(pmap.GammaIndexFlat(j, la, ch));
      // A gamma_decay level: its photon amplitudes move with any free number of the level (the record's, its energy,
      // its particle amplitudes), so their cotangents are needed for the pull-back (GammaDecayMap.h).
      // With external capture also with the final states' amplitudes (GammaDecaySet::InputIndices).
      const GammaDecaySet &decays = compound->GetGammaDecays();
      const int record = decays.RecordOf(j, la);
      if (record < 0) continue;
      bool anyFree = false;
      for (const int index : decays.InputIndices(pmap, record)) anyFree = anyFree || isFree(index);
      if (anyFree)
        for (const int ch : decays.Levels()[record].photonChannels) mask.g[j - 1][la - 1][ch - 1] = 1;
    }
  }
  return mask;
}

}  // namespace

bool ComputeResidualJacobian(CNuc *compound, WorkingCopy *data, azure2::Results &results, const Config &config,
                             const ParamIndexMap &pmap,
                             const vector_matrix_r *shiftDeriv,
                             vector_r &residuals, vector_r &jacobian, int &nCols,
                             bool skipTHM) {
  // Columns of J are the non-fixed (packed) parameters; rows are data points.
  nCols = pmap.NumPacked();
  residuals.clear();
  jacobian.clear();

  // THM (HOES) segments are not supported here: the adjoint differentiates
  // the T-matrix observable and the fast forward path computes the standard
  // cross section, both wrong for HOES. Bail so the LM caller falls back to
  // MIGRAD (which handles THM through Gradient()'s hybrid path).  A caller
  // that fills THM rows itself (AZUREAPI) asks for them to be left zero.
  if (!skipTHM)
    for (int i = 1; i <= data->Data().NumSegments(); i++)
      if (data->Data().GetSegment(i) && data->Data().GetSegment(i)->IsTHM()) return false;

  // One row per (non-null) point, segments in order (THM rows left zero).
  const int nSeg = data->Data().NumSegments();
  std::vector<std::vector<int>> rowOf(nSeg + 1);
  std::vector<double> norm(nSeg + 1, 0.0);
  std::vector<int> normPacked(nSeg + 1, -1);
  int totalRows = 0;
  for (int i = 1; i <= nSeg; i++) {
    const Segment *segment = data->Data().GetSegment(i);
    if (!segment) continue;
    const int nPoints = segment->NumPoints();
    rowOf[i].assign(nPoints, -1);
    for (int pointIdx = 0; pointIdx < nPoints; pointIdx++)
      if (segment->GetPoint(*data, pointIdx + 1)) rowOf[i][pointIdx] = totalRows++;
    norm[i] = data->State(*segment).norm;
    // A norm-group follower's rows depend on its leader's norm parameter.
    const int normOwner = data->Data().NormOwner(i);
    const int full = data->Data().GetSegment(normOwner)->IsVaryNorm() ? pmap.NormIndex(normOwner) : -1;
    normPacked[i] = (full >= 0 && full < pmap.NumFull()) ? pmap.FullToPacked(full) : -1;
  }
  residuals.assign(totalRows, 0.0);
  jacobian.assign((size_t)totalRows * nCols, 0.0);
  if (FormalPartsOutsideAdjoint(compound, pmap) && nSeg > 0) {
    residuals.clear();
    jacobian.clear();
    return false;
  }

  const bool timing = azure2::support::ProcessEnvironment().Set(azure2::support::Switch::kGradTiming);
  const long long comb0 = config.diagnostics->foldWeightNs.load();
  const auto ts0 = std::chrono::steady_clock::now();

  // r = (model/J - data*n)/(err*n) (DataResidual); d r/d model = 1/(err*n*J) yields d r/d(E,gamma) directly.
  auto fitBarOf = [&](const azure2::GradPoint &m, double) {
    const double denom = m.segment->LabError(m.pointIdx) * norm[m.segNum];
    return (denom != 0.0) ? 1.0 / denom / Segment::Jacobian(data->State(*m.segment), 0, m.pointIdx) : 0.0;
  };
  const AMatrixFunc::FreeMask mask = BuildFreeMask(compound, pmap);
  std::vector<azure2::GradPoint> points;
  std::vector<double> models, rows;
  int bailSeg = 0, numTasks = 0;
  const bool ok = azure2::EnergyTaskRows(compound, data, results, config, shiftDeriv, pmap, &mask, fitBarOf, points,
                                         models, rows, bailSeg, numTasks);
  if (!ok) {
    ReportBail(*data, "Jacobian/gradient", data->Data().GetSegment(bailSeg), bailSeg);
  } else {
    for (std::size_t k = 0; k < points.size(); k++) {
      const azure2::GradPoint &m = points[k];
      const int row = rowOf[m.segNum][m.pointIdx];
      if (row < 0) continue;
      const double n = norm[m.segNum];
      const double labErr = m.segment->LabError(m.pointIdx), jac = Segment::Jacobian(data->State(*m.segment), 0, m.pointIdx);
      const azure2::Residual res = azure2::DataResidual(models[k], jac, m.segment->LabValue(m.pointIdx), labErr, n);
      residuals[row] = res.Pull();
      double *out = jacobian.data() + (size_t)row * nCols;
      std::copy(rows.begin() + k * nCols, rows.begin() + (k + 1) * nCols, out);
      // Norm column:  d r / d n = -(model/J)/(err*n^2).
      const int np = normPacked[m.segNum];
      if (np >= 0 && np < nCols && labErr != 0.0 && n != 0.0) out[np] += -res.model / (labErr * n * n);
    }
  }

  if (timing) {
    const int c = ++config.diagnostics->jacobianRowTimings;
    if (c <= 3 || c % 50 == 0)
      std::cerr << "[grad-timing] Jacobian rows " << c << ": " << totalRows << " rows, " << numTasks
                << " energy tasks, "
                << std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ts0).count()
                << " ms; target-effect combiner weights " << (config.diagnostics->foldWeightNs.load() - comb0) / 1.0e6
                << " ms (thread-summed)" << std::endl;
  }

  if (!ok) {
    residuals.clear();
    jacobian.clear();
  }
  return ok;
}

bool ComputeModelGradients(CNuc *compound, WorkingCopy *data, azure2::Results &results, const Config &config,
                           const ParamIndexMap &pmap,
                           const vector_matrix_r *shiftDeriv,
                           std::map<Unit *, vector_r> &gradByPoint,
                           bool skipTHM) {
  gradByPoint.clear();
  const int nCols = pmap.NumPacked();
  const int nSeg = data->Data().NumSegments();
  // The rows of THM (HOES) points come from the HOES adjoint (ComputeTHMRows): without skipTHM their presence takes
  // the caller off the analytic path; with it they are left zero for the caller to fill.
  for (int i = 1; i <= nSeg; i++) {
    const Segment *segment = data->Data().GetSegment(i);
    if (!segment) continue;
    if (segment->IsTHM() && !skipTHM) return false;
    for (int pointIdx = 0; pointIdx < segment->NumPoints(); pointIdx++) {
      Unit *pt = segment->GetPoint(*data, pointIdx + 1);
      if (pt) gradByPoint[pt] = vector_r(nCols, 0.0);
    }
  }
  if (FormalPartsOutsideAdjoint(compound, pmap) && nSeg > 0) {
    gradByPoint.clear();
    return false;
  }
  // fitBar = d(model)/d(model) = 1 gives d(model)/d(E,gamma) directly.
  auto fitBarOf = [](const azure2::GradPoint &, double) { return 1.0; };
  const AMatrixFunc::FreeMask mask = BuildFreeMask(compound, pmap);
  std::vector<azure2::GradPoint> points;
  std::vector<double> models, rows;
  int bailSeg = 0, numTasks = 0;
  if (!azure2::EnergyTaskRows(compound, data, results, config, shiftDeriv, pmap, &mask, fitBarOf, points, models, rows,
                              bailSeg, numTasks)) {
    ReportBail(*data, "Jacobian/gradient", data->Data().GetSegment(bailSeg), bailSeg);
    gradByPoint.clear();
    return false;
  }
  for (std::size_t k = 0; k < points.size(); k++)
    std::copy(rows.begin() + k * nCols, rows.begin() + (k + 1) * nCols, gradByPoint[points[k].point].begin());
  return true;
}

// THM (HOES) rows: see the declaration in AZUREGrad.h.  Shared by the API
// (residual Jacobian, chi2 gradient, model gradients) and the CLI band.
std::vector<THMRows> ComputeTHMRows(CNuc *lc, WorkingCopy *ld, azure2::Results &results, const Config &config,
                                    const vector_r &full, const ParamIndexMap &pmap) {
  std::vector<THMRows> rows;
  const bool brune = (config.paramMask & Config::USE_BRUNE_FORMALISM);
  const int nCols = pmap.NumPacked();

  // Row bookkeeping mirrors ComputeResidualJacobian / EvaluateFilledChi2: one
  // row per non-null point, segments in order.
  int row = 0;
  for (int s = 1; s <= ld->Data().NumSegments(); s++) {
    const Segment *seg = ld->Data().GetSegment(s);
    if (!seg) continue;
    if (seg->IsTHM()) {
      THMRows tr;
      tr.segment = s;
      tr.firstRow = row;
      rows.push_back(tr);
    }
    for (int pid = 0; pid < seg->NumPoints(); pid++)
      if (seg->GetPoint(*ld, pid + 1)) row++;
  }
  if (rows.empty()) return rows;

  // Models of every THM point at the currently filled compound.  The points
  // are independent (as in every other forward pass, the compound is only
  // read), so they are evaluated in parallel: this runs 2 x (free E and gamma
  // columns) + 1 times per Jacobian, and serially it left all but one thread
  // idle -- 280 s for the 230 columns of the 12C+12C THM example against a
  // 0.67 s forward pass.
  std::vector<std::pair<int, int>> order;  // (row block t, point index pid)
  for (size_t t = 0; t < rows.size(); t++) {
    const Segment *seg = ld->Data().GetSegment(rows[t].segment);
    for (int pid = 0; pid < seg->NumPoints(); pid++)
      if (seg->GetPoint(*ld, pid + 1)) order.push_back(std::make_pair((int)t, pid));
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
      const Segment *seg = ld->Data().GetSegment(rows[order[k].first].segment);
      out[order[k].first][slot[k]] =
          azure2::EvaluateModelPoint(*ld, *seg, order[k].second, lc, config, results, azure2::LeafValues::Calculate);
    }
  };
  auto fillCompound = [&](const vector_r &p) {
    lc->FillCompoundFromParams(p);
    if (brune) lc->CalcShiftFunctions(config);
  };

  // Norms and shifts are the caller's; the THM model depends on the R-matrix
  // parameters (and the shift already applied to the point energies).
  ld->FillNormsFromParams(full);
  ld->FillEnergyShiftsFromParams(full, results, lc, &config);
  fillCompound(full);
  std::vector<vector_r> m0, mp, mm;
  models(m0);
  for (size_t t = 0; t < rows.size(); t++) {
    rows[t].m = m0[t];
    rows[t].Jm.assign(m0[t].size() * (size_t)nCols, 0.0);
    rows[t].J.assign(m0[t].size() * (size_t)nCols, 0.0);
  }

  // d m / d p, E and gamma columns -- and the THM coherent backgrounds
  // (cbackground=), which only THM points see -- by the HOES adjoint: one
  // reverse sweep per point (ThmPointModelGradient).  A configuration it does
  // not cover, or AZURE2_THM_FD=1, takes the central differences below.
  bool analytic = !ThmFiniteDifferencesRequested(config);
  for (size_t t = 0; t < rows.size() && analytic; t++)
    if (ld->Data().GetSegment(rows[t].segment)->HasComponents()) analytic = false;
  if (analytic) {
    vector_matrix_r shiftDeriv;
    const vector_matrix_r *sdp = nullptr;
    if (brune) {
      shiftDeriv = BuildShiftDerivTable(lc, config);
      sdp = &shiftDeriv;
    }
    bool lineshape = false;
    for (size_t t = 0; t < rows.size(); t++)
      if (ld->Data().GetSegment(rows[t].segment)->GetThmLineshape()) lineshape = true;
    const ThmAdjointContext ctx = BuildThmAdjointContext(lc, config, sdp, lineshape);
    const int nCoherent = ld->Data().NumThmCoherentParams();
    std::vector<size_t> slot(order.size());
    {
      std::vector<size_t> next(rows.size(), 0);
      for (size_t k = 0; k < order.size(); k++) slot[k] = next[order[k].first]++;
    }
    bool bail = false;
    const long n = (long)order.size();
#pragma omp parallel
    {
      GradAccum accum;
      accum.Init(lc);
      vector_r fullRow(pmap.NumFull(), 0.0);
      std::vector<double> coherentBar(nCoherent, 0.0);
      THMMatrixFunc thm(lc, config);  // this thread's solver for the region (re-pointed per point)
#pragma omp for schedule(dynamic)
      for (long k = 0; k < n; k++) {
        if (bail) continue;
        const Segment *seg = ld->Data().GetSegment(rows[order[k].first].segment);
        Unit *pt = seg->GetPoint(*ld, order[k].second + 1);
        accum.Zero();
        std::fill(coherentBar.begin(), coherentBar.end(), 0.0);
        if (!ThmPointModelGradient(*ld, pt, results, lc, config, ctx, thm, accum, coherentBar)) {
#pragma omp critical
          bail = true;
          continue;
        }
        std::fill(fullRow.begin(), fullRow.end(), 0.0);
        accum.Scatter(pmap, fullRow);
        double *row = &rows[order[k].first].Jm[slot[k] * nCols];
        for (int c = 0; c < nCols; c++) {
          const int f = pmap.PackedToFull(c);
          const ParamDesc &desc = pmap.Desc(f);
          if (IsLevelSlot(desc.kind))
            row[c] = fullRow[f];
          else if (desc.kind == ParamKind::ThmCoherent && desc.channel >= 1 && desc.channel <= nCoherent)
            row[c] = coherentBar[desc.channel - 1];
        }
      }
    }
    if (bail) {
      analytic = false;
      for (size_t t = 0; t < rows.size(); t++) std::fill(rows[t].Jm.begin(), rows[t].Jm.end(), 0.0);
      if (azure2::support::ProcessEnvironment().Set(azure2::support::Switch::kGradDebug))
        std::cerr << "[grad] THM adjoint: a point outside its path (a target effect with an automatic "
                     "tolerance, RMC, ...) -> central differences for the THM rows" << std::endl;
    }
  }

  // Central differences (AZURE2_THM_FD=1, or the adjoint declined).
  for (int c = 0; c < nCols && !analytic; c++) {
    const int f = pmap.PackedToFull(c);
    const ParamKind kind = pmap.Desc(f).kind;
    if (kind == ParamKind::ThmCoherent) {
      const double x0 = full[f];
      const double h = 1.0e-6 * (std::fabs(x0) + 1.0);
      vector_r pp = full;
      pp[f] = x0 + h;
      ld->FillThmCoherentFromParams(pp, lc);
      models(mp);
      pp[f] = x0 - h;
      ld->FillThmCoherentFromParams(pp, lc);
      models(mm);
      ld->FillThmCoherentFromParams(full, lc);
      for (size_t t = 0; t < rows.size(); t++)
        for (size_t i = 0; i < m0[t].size(); i++)
          rows[t].Jm[i * nCols + c] = (mp[t][i] - mm[t][i]) / (2.0 * h);
      continue;
    }
    if (!IsLevelSlot(kind)) continue;
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
    const Segment *seg = ld->Data().GetSegment(tr.segment);
    std::vector<double> d, e, ml, Jml;
    for (int pid = 0; pid < seg->NumPoints(); pid++) {
      Unit *pt = seg->GetPoint(*ld, pid + 1);
      if (!pt) continue;
      d.push_back(seg->LabValue(pid));
      e.push_back(seg->LabError(pid));
      pt->SetModelValue(results, tr.m[d.size() - 1]);  // the residuals below take the model in the datum's frame
      azure2::AppendInDataFrame(*ld, *seg, pid, tr.m[d.size() - 1], tr.Jm.data() + (d.size() - 1) * nCols, nCols, ml, Jml);
    }
    const size_t n = d.size();
    if (ld->Data().ThmGroupOf(tr.segment) >= 0) continue;  // a THM experiment: below
    double Smm = 0.0, Smd = 0.0;
    for (size_t i = 0; i < n; i++) {
      if (e[i] == 0.0) continue;
      double w = 1.0 / (e[i] * e[i]);
      Smm += ml[i] * ml[i] * w;
      Smd += ml[i] * d[i] * w;
    }
    // s multiplies the model (s = 1/n, n the norm on the data).
    bool varyScale = false;
    double s;
    if (seg->IsProfiledNorm()) {
      varyScale = (Smd > 0.0 && Smm > 0.0);
      s = varyScale ? Smd / Smm : 1.0;
    } else {
      double norm = ld->State(*seg).norm;
      s = (norm != 0.0) ? 1.0 / norm : 0.0;
    }
    tr.r.assign(n, 0.0);
    for (size_t i = 0; i < n; i++)
      if (e[i] != 0.0) tr.r[i] = (s * ml[i] - d[i]) / e[i];

    vector_r ds(nCols, 0.0);
    if (varyScale) {
      for (int c = 0; c < nCols; c++) {
        double sdJ = 0.0, smJ = 0.0;
        for (size_t i = 0; i < n; i++) {
          if (e[i] == 0.0) continue;
          double w = 1.0 / (e[i] * e[i]);
          sdJ += d[i] * Jml[i * nCols + c] * w;
          smJ += ml[i] * Jml[i * nCols + c] * w;
        }
        ds[c] = (sdJ - 2.0 * s * smJ) / Smm;
      }
    }
    for (size_t i = 0; i < n; i++) {
      if (e[i] == 0.0) continue;
      for (int c = 0; c < nCols; c++)
        tr.J[i * nCols + c] = (s * Jml[i * nCols + c] + ml[i] * ds[c]) / e[i];
    }
    tr.s = s;
    tr.profiled = varyScale;
    tr.ds = ds;
  }

  // THM experiments: one linear profile (shared norm, background) over the
  // points of all their segments, differentiated through the profile.
  for (int g = 0; g < ld->Data().NumThmGroups(); g++) {
    std::vector<size_t> members;  // indices into rows
    for (size_t t = 0; t < rows.size(); t++)
      if (ld->Data().ThmGroupOf(rows[t].segment) == g) members.push_back(t);
    if (members.empty()) continue;
    ld->ProfileThmGroup(g, results);  // norms and profile from the models set above
    const ThmProfile &prof = results.objective.thmProfiles[g];
    std::vector<double> m, d, e, energy, Jm;
    for (size_t t : members) {
      const Segment *seg = ld->Data().GetSegment(rows[t].segment);
      size_t i = 0;
      for (int pid = 0; pid < seg->NumPoints(); pid++) {
        Unit *pt = seg->GetPoint(*ld, pid + 1);
        if (!pt) continue;
        d.push_back(seg->LabValue(pid));
        e.push_back(seg->LabError(pid));
        energy.push_back(pt->GetCMEnergy());
        azure2::AppendInDataFrame(*ld, *seg, pid, rows[t].m[i], rows[t].Jm.data() + i * nCols, nCols, m, Jm);
        i++;
      }
    }
    std::vector<double> J, G;
    ThmProfileDerivative(prof, m, d, e, energy, Jm, nCols, J, G);
    size_t k = 0;
    for (size_t t : members) {
      THMRows &tr = rows[t];
      const size_t n = tr.m.size();
      tr.r.assign(n, 0.0);
      for (size_t i = 0; i < n; i++) tr.r[i] = prof.Residual(m[k + i], d[k + i], e[k + i], energy[k + i]);
      tr.J.assign(J.begin() + k * nCols, J.begin() + (k + n) * nCols);
      tr.G.assign(G.begin() + k * nCols, G.begin() + (k + n) * nCols);
      tr.s = prof.s;
      tr.profiled = false;  // the whole derivative is in J and G
      tr.ds.assign(nCols, 0.0);
      k += n;
    }
  }
  return rows;
}
