// The gathers of the objective's diagonal priors (DiagonalPriors.h; docs/upgrade/objective_arrays.md section 2.1): the
// one place where the rule of eligibility of a normalization, energy-shift, target, apparatus or abundance prior (T3,
// T4, T8, T7's per-nucleus priors) is written, and the value it reads.  The value, the gradient, the curvature, the
// least-squares rows, the sampler and the chiSquared.out writer all take their priors from here.
#include "DiagonalPriors.h"

#include <algorithm>
#include <limits>

#include "AbundanceSamples.h"
#include "ApparatusParameters.h"
#include "ParameterLimitsManager.h"
#include "Segment.h"
#include "TargetEffect.h"
#include "WorkingCopy.h"

namespace {

void Add(DiagonalPriors &out, int index, double value, double centre, double sigma, PriorSource source, int group) {
  DiagonalPrior e;
  e.index = index;
  e.value = value;
  e.centre = centre;
  e.sigma = sigma;
  e.source = source;
  e.group = group;
  out.push_back(e);
}

}  // namespace

void AppendNormPrior(const WorkingCopy &data, int s, const ParamIndexMap *pmap, DiagonalPriors &out) {
  const Segment *seg = data.Data().GetSegment(s);
  // A segment of a THM experiment has its norm profiled with the experiment's (or fixed); objective.md D4, D5.
  if (!seg || data.Data().ThmGroupOf(s) >= 0 || seg->IsProfiledNorm() || !seg->IsVaryNorm() || seg->IsNormPriorReplaced())
    return;
  const double n0 = seg->GetNominalNorm(), sigma = n0 / 100. * seg->GetNormError();
  if (sigma != 0.) Add(out, pmap ? pmap->NormIndex(s) : -1, data.State(*seg).norm, n0, sigma, PriorSource::kNorm, s);
}

void AppendShiftPriors(const WorkingCopy &data, int s, const ParamIndexMap *pmap, DiagonalPriors &out) {
  const Segment *seg = data.Data().GetSegment(s);
  if (!seg) return;
  const SegmentState &state = data.State(*seg);
  // Each term is penalised from the centre of its prior, 0 where none is written (an old line; DeBoer: "It should
  // measure the penalty from 0"); the value on the line is only the start (parameter_kinds.md section 4).  An old
  // line's prior_centre row moves the constant term's centre (GetNominalEnergyShift).
  if (seg->IsVaryEnergyShift() && seg->GetEnergyShiftError() != 0.0 && !seg->IsShiftPriorReplaced())
    Add(out, pmap ? pmap->EnergyShiftIndex(s) : -1, state.energyShift, seg->GetNominalEnergyShift(),
        seg->GetEnergyShiftError(), PriorSource::kShift, s);
  if (seg->IsVaryEnergyShiftSqrt() && seg->GetEnergyShiftSqrtError() != 0.0 && !seg->IsShiftSqrtPriorReplaced())
    Add(out, pmap ? pmap->EnergyShiftSqrtIndex(s) : -1, state.energyShiftSqrt, seg->GetEnergyShiftSqrtCentre(),
        seg->GetEnergyShiftSqrtError(), PriorSource::kShiftSqrt, s);
  // The shift function's parameters are penalised once, on the segment that owns them.
  if (seg->OwnsShiftFn() && seg->HasShiftFn())
    for (int j = 0; j < (int)state.shiftFnValues.size(); j++)
      if (seg->IsVaryShiftFn(j) && seg->GetShiftFnErrors()[j] != 0.0 && !seg->IsShiftFnPriorReplaced(j))
        Add(out, pmap ? pmap->EnergyShiftFnIndex(s, j) : -1, state.shiftFnValues[j], seg->GetShiftFnCentres()[j],
            seg->GetShiftFnErrors()[j], PriorSource::kShiftFn, s);
}

void AppendSegmentPriors(const WorkingCopy &data, int s, const ParamIndexMap *pmap, DiagonalPriors &out) {
  AppendNormPrior(data, s, pmap, out);
  AppendShiftPriors(data, s, pmap, out);
}

void AppendAbundancePriors(const WorkingCopy &data, const ParamIndexMap *pmap, DiagonalPriors &out) {
  const AbundanceSamples &samples = data.Data().Abundances();
  const std::vector<double> &values = data.AbundanceValues();
  for (int k = 1; k <= samples.NumParams(); k++) {
    int sample, nucleus;
    samples.Of(k, sample, nucleus);
    const SampleAbundance &e = samples.Sample(sample).entries[nucleus - 1];
    if (!(e.sigma > 0.0) || e.priorReplaced) continue;
    Add(out, pmap ? pmap->AbundanceIndex(k) : -1, (k - 1 < (int)values.size()) ? values[k - 1] : std::numeric_limits<double>::quiet_NaN(), e.nominal,
        e.sigma, PriorSource::kAbundance, sample);  // no value before the first hand-over
  }
}

DiagonalPriors DiagonalPriorsOf(const WorkingCopy &data, const ParamIndexMap &pmap) {
  DiagonalPriors all;
  const int numSegments = data.Data().NumSegments();
  for (int s = 1; s <= numSegments; s++) AppendNormPrior(data, s, &pmap, all);
  for (int s = 1; s <= numSegments; s++) AppendShiftPriors(data, s, &pmap, all);
  for (PriorSource source : {PriorSource::kScreening, PriorSource::kDensity, PriorSource::kDoppler,
                             PriorSource::kThickness, PriorSource::kApparatus})
    AppendTargetPriors(data, source, &pmap, all);
  AppendAbundancePriors(data, &pmap, all);
  DiagonalPriors indexed;
  indexed.reserve(all.size());
  for (const DiagonalPrior &e : all)
    if (e.index >= 0) indexed.push_back(e);
  return indexed;
}

void AppendTargetPriors(const WorkingCopy &data, PriorSource source, const ParamIndexMap *pmap, DiagonalPriors &out) {
  const DataSet &set = data.Data();
  auto add = [&](int index, double value, double centre, double sigma, int group) {
    Add(out, index, value, centre, sigma, source, group);
  };
  switch (source) {
    case PriorSource::kScreening:  // U_e in eV, the value of the parameter's first copy
      for (int k = 1; k <= set.NumScreeningParams(); k++) {
        const TargetEffect *te = set.GetScreeningTargetEffect(k);
        if (!te->HasScreening() || !te->IsScreeningVaried() || te->IsScreeningPriorReplaced() || te->GetScreeningSigma() <= 0.0)
          continue;
        add(pmap ? pmap->ScreeningIndex(k) : -1, data.EffectValues(set.ScreeningCopies(k).front()).screeningUe,
            te->GetScreeningNominal(), te->GetScreeningSigma(), 0);
      }
      break;
    case PriorSource::kDensity:  // N in atoms/cm^2, the lead copy's
      for (int k = 1; k <= set.NumDensityParams(); k++) {
        const TargetEffect *te = set.GetDensityTargetEffect(k);
        if (!te->HasDensityFit() || !te->IsDensityVaried() || te->IsDensityPriorReplaced() || te->GetDensitySigma() <= 0.0)
          continue;
        add(pmap ? pmap->DensityIndex(k) : -1, data.EffectValues(set.GetDensityLeadCopy(k)).density, te->GetDensityCentre(),
            te->GetDensitySigma(), 0);
      }
      break;
    case PriorSource::kDoppler:  // T in K, the value of the parameter's first copy
      for (int k = 1; k <= set.NumDopplerParams(); k++) {
        const TargetEffect *te = set.GetDopplerTargetEffect(k);
        if (!te->HasDopplerFit() || !te->IsDopplerVaried() || te->IsDopplerPriorReplaced() || te->GetDopplerSigma() <= 0.0)
          continue;
        add(pmap ? pmap->DopplerIndex(k) : -1, data.EffectValues(set.DopplerCopies(k).front()).dopplerTemperature,
            te->GetDopplerCentre(), te->GetDopplerSigma(), 0);
      }
      break;
    case PriorSource::kThickness:  // n in atoms/b; a thickness-group follower is not varied (IsThicknessVaried)
      for (int k = 1; k <= set.NumThicknessParams(); k++) {
        const Segment &seg = *set.GetSegment(set.GetThicknessSegment(k));
        if (!seg.IsThicknessVaried() || seg.IsThicknessPriorReplaced() || !(seg.GetThicknessSigma() > 0.0)) continue;
        add(pmap ? pmap->ThicknessIndex(k) : -1, data.State(seg).thickness, seg.GetThicknessCentre(), seg.GetThicknessSigma(), 0);
      }
      break;
    case PriorSource::kApparatus:  // the kernel parameter's engine unit; group: the row of ApparatusRows()
      for (int k = 1; k <= set.NumApparatusParams(); k++) {
        const ApparatusParam &a = set.GetApparatusParam(k);
        if (!a.declaration.vary || a.priorReplaced || !(a.declaration.sigma > 0.0)) continue;
        add(pmap ? pmap->ApparatusIndex(k) : -1,
            a.Row().read(*set.GetTargetEffect(a.lead), a.declaration.slot, data.EffectValues(a.lead)), a.declaration.centre,
            a.declaration.sigma, a.row);
      }
      break;
    default:
      break;
  }
}

bool NuisancePriorOf(const ParameterLimitsManager *limits, const PriorCovariance *pc, int f, double &centre, double &sigma) {
  if (!limits) return false;
  // A parameter of the prior covariance has no diagonal prior of its own: the multivariate Gaussian replaces it.
  if (pc && std::find(pc->FullIndex().begin(), pc->FullIndex().end(), f) != pc->FullIndex().end()) return false;
  return limits->NuisancePrior(f, centre, sigma);
}

DiagonalPriors NuisancePriorsOf(const ParameterLimitsManager *limits, const PriorCovariance *pc, const vector_r &p) {
  DiagonalPriors out;
  if (!limits) return out;
  for (int f = 0; f < (int)p.size(); f++) {
    double centre, sigma;
    if (NuisancePriorOf(limits, pc, f, centre, sigma)) Add(out, f, p[f], centre, sigma, PriorSource::kNuisance, 0);
  }
  return out;
}
