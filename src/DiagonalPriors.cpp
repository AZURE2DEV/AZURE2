// The gathers of the objective's diagonal priors (DiagonalPriors.h; docs/upgrade/objective_arrays.md section 2.1): the
// one place where the value path's rule of eligibility of a target or apparatus prior (T8) is written, and the value
// it reads.
#include "DiagonalPriors.h"

#include "ApparatusParameters.h"
#include "Segment.h"
#include "TargetEffect.h"
#include "WorkingCopy.h"

void AppendTargetPriors(const WorkingCopy &data, PriorSource source, const ParamIndexMap *pmap, DiagonalPriors &out) {
  const DataSet &set = data.Data();
  auto add = [&](int index, double value, double centre, double sigma, int group) {
    DiagonalPrior e;
    e.index = index;
    e.value = value;
    e.centre = centre;
    e.sigma = sigma;
    e.source = source;
    e.group = group;
    out.push_back(e);
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
