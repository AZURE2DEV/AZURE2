#include "AZUREGrad.h"
#include "GradFreeMask.h"
#include "Objective.h"
#include "RunDiagnostics.h"
#include "azure2/support/Environment.h"
#include "GradEnergyTasks.h"
#include "ModelPass.h"
#include "Profile.h"
#include <atomic>
#include <chrono>
#include <iostream>
#include <cstdlib>
#include "CompoundSystem.h"
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
#include "azure2/support/FloatingPointPrecision.h"
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
 * bound/unbound branch of CompoundSystem::CalcShiftFunctions.  In the Park formalism
 * CalcShiftFunctions (which the callers run first) already holds the same
 * derivative on each level, from the same function at the same energy: it is
 * read from there instead of being differentiated again.
 */
vector_matrix_r BuildShiftDerivTable(const CompoundSystem *compound, const LevelState &state, const Config &configure) {
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
        if (park && state.ShiftsValidAt(j, la, state.E(j, la))) {
          tab[j - 1][la - 1][ch - 1] = state.ShiftDerivative(j, la, ch);
          continue;
        }
        int l = channel->GetL();
        double levelEnergy = state.E(j, la);
        double resonanceEnergy = levelEnergy - (pair->GetSepE() + pair->GetExE());
        double dS = ChannelFunc(pair, channel->GetRadius())
                        .ShiftDerivative(l, pair->SpinOrbitTwoJ(jg->GetJ()), resonanceEnergy);
        tab[j - 1][la - 1][ch - 1] = dS;
      }
    }
  }
  return tab;
}

/*!
 * Allocates and zeroes the physics-coordinate gradient accumulators.
 */
void GradAccum::Init(const CompoundSystem *compound, const LevelState &state) {
  source = compound;
  sourceState = &state;
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
  abundance.clear();
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
  if (o.abundance.size() > abundance.size()) abundance.resize(o.abundance.size(), 0.0);
  for (int k = 0; k < (int)o.abundance.size(); k++) abundance[k] += o.abundance[k];
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
  if (decays) decays->PullBack(*source, *sourceState, *this, map, gradFull);
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
  // Isospin levels (IsospinTies.h): a derived slot is a linear function of free slots, so its gradient passes on to
  // them (chain rule); its own entry is cleared (it is fixed).  Last, after everything else has reached the slots.
  for (const IsospinTie &t : map.IsospinTies()) {
    if (t.slot < 0 || t.slot >= (int)gradFull.size()) continue;
    const double g = gradFull[t.slot];
    for (const auto &src : t.sources)
      if (src.first >= 0 && src.first < (int)gradFull.size()) gradFull[src.first] += src.second * g;
    gradFull[t.slot] = 0.0;
  }
}

/*!
 * Builds the parameter-index map: the level block from the compound's layout
 * (CompoundSystem::ParameterLayout), then by mirroring, in the exact same order, the
 * assignment loops in:
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
ParamIndexMap BuildParamIndexMap(const CompoundSystem *compound, const WorkingCopy *data,
                                 const std::vector<bool> &fixed) {
  ParamIndexMap map;

  // --- R-matrix parameters: the level block, as the compound's layout places it. ---
  std::vector<int> numLevels, numChannels;
  for (int j = 1; j <= compound->NumJGroups(); j++) {
    const JGroup *jGroup = compound->GetJGroup(j);
    numLevels.push_back(jGroup->NumLevels());
    numChannels.push_back(jGroup->NumChannels());
  }
  map.SizeLevelBlock(numLevels, numChannels);
  for (const LevelParameterLayout::Slot &slot : compound->ParameterLayout().Slots()) {
    const int j = slot.group, la = slot.level, ch = slot.channel;
    if (slot.kind == LevelParameterLayout::Slot::Kind::kEnergy)
      map.AppendLevelEnergy(j, la);
    else
      map.AppendAmplitude(j, la, ch, compound->GetGammaDecays().OwnsSlot(j, la, ch));
  }

  // --- Isospin shifts, then mix amplitudes, only with <isospin> records (CompoundSystem::FillMnParams). ---
  for (int k = 0; k < compound->GetIsospin().NumParams(); k++) {
    const ParameterRecord record = compound->GetIsospin().RecordOf(k);
    map.AppendIsospin(record.kind, record.jGroup, record.level, record.channel, record.nucleus);
  }
  map.SetIsospinTies(compound->GetIsospin().Ties());

  if (data) {
    const DataSet &set = data->Data();
    // --- Isotopic abundances: one per (sample, nucleus), only with an <abundances> block. ---
    for (int k = 1; k <= set.Abundances().NumParams(); k++) {
      int sample = -1, nucleus = -1;
      set.Abundances().Of(k, sample, nucleus);
      map.AppendAbundance(k, sample, nucleus);
    }
    // --- Norm parameters: one per varying-norm segment. ---
    for (int s = 1; s <= set.NumSegments(); s++) {
      const Segment *segment = set.GetSegment(s);
      if (segment && segment->IsVaryNorm() && !segment->IsProfiledNorm()) map.AppendNorm(s);
    }
    // --- Energy-shift parameters: one per segment (all segments). ---
    for (int s = 1; s <= set.NumSegments(); s++)
      if (set.GetSegment(s)) map.AppendEnergyShift(s);
    // --- sqrt(E) energy-shift coefficients: one per segment (all segments). ---
    for (int s = 1; s <= set.NumSegments(); s++)
      if (set.GetSegment(s)) map.AppendEnergyShiftSqrt(s);
    // --- Shift-function parameters: those of each segment that declares one (WorkingCopy::FillMnParams). ---
    for (int s = 1; s <= set.NumSegments(); s++) {
      const Segment *segment = set.GetSegment(s);
      if (!segment || !segment->OwnsShiftFn() || !segment->HasShiftFn()) continue;
      for (int j = 0; j < segment->NumShiftFnParams(); j++) map.AppendEnergyShiftFn(s, j);
    }
    // --- Electron-screening potentials: one per target with a screening field (WorkingCopy::FillMnParams). ---
    for (int k = 1; k <= set.NumScreeningParams(); k++) map.AppendScreening(k);
    // --- Fitted target densities: one per densityfit line or density group (WorkingCopy::FillMnParams). ---
    // No analytic derivative: the target thickness enters through the depth window, finite differences.
    for (int k = 1; k <= set.NumDensityParams(); k++) map.AppendDensity(k);
    // --- Fitted Doppler temperatures: one per 'doppler T vary sigma' line (WorkingCopy::FillMnParams). ---
    // Analytic: d/dT through the Doppler weights (Unit::DopplerWeights), GradAccum::AddDoppler.
    for (int k = 1; k <= set.NumDopplerParams(); k++) map.AppendDoppler(k);
    // --- Fitted sample thicknesses: one per transmission / yield line with the fit form (WorkingCopy::FillMnParams). ---
    // Analytic: d/dn of the transform (the tree's backward, AddParamCotangents), GradAccum::AddThickness.
    for (int k = 1; k <= set.NumThicknessParams(); k++) map.AppendThickness(k, set.GetThicknessSegment(k));
    // --- THM coherent backgrounds (cbackground=): the last block but one, if any. ---
    for (int k = 0; k < set.NumThmCoherentParams(); k++) map.AppendThmCoherent(k + 1);
  }
  AppendApparatusParams(map, data);  // the apparatus parameters, last (ApparatusParameters.cpp)
  // --- Apply the fixed mask and build packed<->full lookups. ---
  map.Pack(fixed);
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
                             int xsComponent, const CompoundSystem *compound, const Config &config, int row,
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
  // The weights the forward fold of this gradient call kept (PointFold.h, KeepFoldWeights), or computed here.
  std::vector<double> c;
  bool exact = false;
  if (const std::vector<double> *kept = KeptWeights(*point, results)) {
    c = *kept;
    exact = true;
  } else if (!config.environment.foldWeightsByDifferences) {
    exact = FoldWeights(*point, data, config, c);
  }
  if (!exact) {
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
    den = FoldKept(*point, data, config, weight, results);
    const double num = FoldKept(*point, data, config, product, results);
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
                        const CompoundSystem *compound, GradAccum &accum) {
  using azure2::experiment::ParamRef;
  for (std::size_t q = 0; q < c.refs->size(); q++) {
    if (c.params[q] == 0.0) continue;
    const ParamRef &r = (*c.refs)[q];
    if (r.kind == ParamRef::Kind::Abundance) {
      accum.AddAbundance(data.Data().Abundances().Ordinal(segment->GetNaturalSample(), r.nucleus), c.params[q]);
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

bool azure2::PointCotangents(WorkingCopy &data, const GradPoint &p, const CompoundSystem *compound, const LevelState &state, const Config &config,
                             azure2::Results &results, double fitBar, int row, GradAccum &params,
                             std::vector<UnitCotangent> &out) {
  if (fitBar == 0.0) return true;
  // The cotangents of the point's expression tree: its leaves (points, or sub-points of a transform's terms) and its
  // parameters.
  ModelCotangents c;
  if (!ModelPointBackward(data, *p.segment, p.pointIdx, compound, state, config, results, fitBar, c)) return false;
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

// The formal-basis parts that the adjoint does not differentiate: the Brune amplitude of an external-capture final
// state that shares its block (CompoundSystem::ECFinalStatesFromFormal, an iteration on the block's formal parameters) and the
// --ignore-externals projector of a channel capture (CompoundSystem::FormalCaptureMasksOutsideAdjoint, external_capture.md
// section 7).  Every adjoint entry point asks this and falls back to finite differences.
bool FormalPartsOutsideAdjoint(const CompoundSystem *compound, const LevelState &state, const ParamIndexMap &pmap) {
  if (compound->HasFormalECFinalBlocks()) return true;
  return compound->FormalCaptureMasksOutsideAdjoint(state, [&pmap](int j, int la, int ch) {
    const int full = pmap.GammaIndexFlat(j, la, ch);
    return full >= 0 && full < pmap.NumFull() && pmap.FullToPacked(full) >= 0;
  });
}

}  // namespace

/// AZURE2_THM_FD=1 (any value but empty or 0): THM rows by central
/// differences, the path before the HOES adjoint (kept for testing).
bool ThmFiniteDifferencesRequested(const Config &config) { return config.environment.thmByDifferences; }


bool AccumulateEGammaGradient(const CompoundSystem *compound, const LevelState &state, WorkingCopy *data, azure2::Results &results, const Config &config,
                              const ParamIndexMap &pmap,
                              const vector_matrix_r *shiftDeriv,
                              const FitBarFn &fitBarFn,
                              GradAccum &accum, const ThmAdjointContext *hoes) {
  // THM (HOES) segments are computed by THMMatrixFunc, not the T-matrix
  // observable this adjoint differentiates; they are left out here. Their
  // contribution to the gradient is added in AZURECalc::Gradient by the HOES
  // adjoint or finite differences of a chi^2 restricted to THM segments, so
  // mixed THM+direct datasets keep the analytic speed-up on the direct part.
  if (FormalPartsOutsideAdjoint(compound, state, pmap)) return false;
  int bailSeg = 0;
  if (!azure2::EnergyTaskGradient(compound, state, data, results, config, shiftDeriv, fitBarFn, accum, bailSeg, hoes)) {
    // On bail the caller ignores `accum` and falls back to finite differences.
    return false;
  }
  if (accum.formalUsed) {
    AMatrixFunc amf(compound, state, config);
    amf.ResolveFormalCotangents(accum, shiftDeriv);
  }
  return true;
}

bool ComputeResidualJacobian(const CompoundSystem *compound, const LevelState &state, WorkingCopy *data, azure2::Results &results, const Config &config,
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
  if (FormalPartsOutsideAdjoint(compound, state, pmap) && nSeg > 0) {
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
  const bool ok = azure2::EnergyTaskRows(compound, state, data, results, config, shiftDeriv, pmap, &mask, fitBarOf, points,
                                         models, rows, bailSeg, numTasks);
  if (ok) {
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

bool ComputeModelGradients(const CompoundSystem *compound, const LevelState &state, WorkingCopy *data, azure2::Results &results, const Config &config,
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
  if (FormalPartsOutsideAdjoint(compound, state, pmap) && nSeg > 0) {
    gradByPoint.clear();
    return false;
  }
  // fitBar = d(model)/d(model) = 1 gives d(model)/d(E,gamma) directly.
  auto fitBarOf = [](const azure2::GradPoint &, double) { return 1.0; };
  const AMatrixFunc::FreeMask mask = BuildFreeMask(compound, pmap);
  std::vector<azure2::GradPoint> points;
  std::vector<double> models, rows;
  int bailSeg = 0, numTasks = 0;
  if (!azure2::EnergyTaskRows(compound, state, data, results, config, shiftDeriv, pmap, &mask, fitBarOf, points, models, rows,
                              bailSeg, numTasks)) {
    gradByPoint.clear();
    return false;
  }
  for (std::size_t k = 0; k < points.size(); k++)
    std::copy(rows.begin() + k * nCols, rows.begin() + (k + 1) * nCols, gradByPoint[points[k].point].begin());
  return true;
}
