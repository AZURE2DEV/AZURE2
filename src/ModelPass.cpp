// The model pass and the residual (ModelPass.h; docs/upgrade/data_model.md step S4).  The combination below is the
// one ESegment computed before S4, in the same arithmetic order, so every model value keeps its bits; what changed is
// where it goes: Results::model, never the base point's leaf slot.
#include "ModelPass.h"

#include <atomic>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "CNuc.h"
#include "Config.h"
#include "EData.h"
#include "EPoint.h"
#include "ESegment.h"
#include "PPair.h"
#include "TargetEffect.h"

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();

/*!
 * Calculate a point; if the calculation throws (a singular matrix, a Coulomb function out of range, ...), mark it
 * failed (EPoint::MarkFailed: NaN, reported once) and return false, so that neither the point nor chi^2 can pass for
 * a result.  Before 2026-11-01 a failed base point gave 0, a failed or non-finite component was left out of a SUM and
 * a failed RATIO denominator gave 0, all without a word: the natural SUM segments of tests/natural_neon on the
 * R-matrix route came back 0 at 84 points.
 */
bool CalculatePoint(EPoint *point, const char *role, CNuc *cnuc, const Config &configure, azure2::Results &results) {
  try {
    point->Calculate(cnuc, configure, results);
    return true;
  } catch (const std::exception &e) {
    point->MarkFailed(results, role, configure, e.what());
  } catch (...) {
    point->MarkFailed(results, role, configure, "unknown exception.");
  }
  return false;
}

// Whether the combination calculates `point` itself: the leaves are not in yet, and the point is its own host.
bool CalculatesItself(const EPoint *point, azure2::LeafValues leaves) {
  return leaves == azure2::LeafValues::Calculate && (!point->IsMapped() || ScopedPointMappingOff::Active());
}

// A sample segment's point whose combination fails: its per-nucleus terms are NaN too (written by this evaluation).
double FailedSample(const ESegment &segment, EPoint *basePoint, CNuc *cnuc, azure2::Results &results) {
  if (segment.IsNatural()) basePoint->SetNaturalTerms(results, std::vector<double>(cnuc->NumNuclei(), kNaN));
  return kNaN;
}

/*!
 * The value of a transmission or yield segment at a point (ESegment::IsTransformed).  The points, base and
 * components, hold plain cross sections; they are summed with their weights into s (and s_t for a yield), at every
 * sub-point when the point carries a sub-point target effect, transformed there, and only then integrated with the
 * base point's kernel (EPoint::IntegrateSubPointValues, on a local vector of the transformed values).
 * EData::CheckTransformedSegments makes sure every term shares the base point's sub-point grid.
 */
double TransformedValue(ESegment &segment, int pointIndex, CNuc *cnuc, const Config &configure,
                        azure2::Results &results, azure2::LeafValues leaves) {
  EPoint *basePoint = segment.GetPoint(pointIndex + 1);
  for (ESegment *comp : segment.GetComponentSegments()) {
    if (!comp || pointIndex >= comp->NumPoints()) continue;
    EPoint *cp = comp->GetPoint(pointIndex + 1);
    if (cp && CalculatesItself(cp, leaves) &&
        !CalculatePoint(cp, comp->IsYieldTotalTerm() ? "yield total" : "SUM component", cnuc, configure, results))
      return FailedSample(segment, basePoint, cnuc, results);
  }
  std::vector<ESegment::TransformTerm> terms;
  if (!segment.TransformTerms(pointIndex, cnuc, terms)) return FailedSample(segment, basePoint, cnuc, results);
  const double n = segment.GetSampleThickness();
  const bool yield = segment.IsYield();
  auto valueOf = [&results](EPoint *p, int c) {
    return (c == 1) ? p->GetFitE1CrossSection(results) : (c == 2) ? p->GetFitE2CrossSection(results)
                                                                  : p->GetFitCrossSection(results);
  };
  // the observable from the term values v_q (a point, or one sub-point of every term)
  auto transform = [&](const std::vector<double> &v) {
    double sum = 0.0, total = 0.0;
    for (size_t q = 0; q < terms.size(); q++) (terms[q].total ? total : sum) += terms[q].weight * v[q];
    return yield ? ESegment::YieldY0(n, sum, total) : ESegment::Transmission(n, sum);
  };

  if (segment.IsNatural()) {
    std::vector<double> naturalTerms(cnuc->NumNuclei(), 0.0);
    for (const ESegment::TransformTerm &t : terms)
      if (!t.total) naturalTerms[t.nucleus - 1] += t.weight * valueOf(t.point, t.xsComponent);
    basePoint->SetNaturalTerms(results, naturalTerms);
  }

  std::vector<double> v(terms.size());
  if (!ESegment::HasSubPointIntegration(basePoint)) {
    for (size_t q = 0; q < terms.size(); q++) v[q] = valueOf(terms[q].point, terms[q].xsComponent);
    return transform(v);
  }

  const int nsp = basePoint->NumSubPoints();
  for (const ESegment::TransformTerm &t : terms)
    if (t.point->NumSubPoints() != nsp) return kNaN;  // refused at load (EData::CheckTransformedSegments)
  std::vector<double> f(nsp);
  // Free-gas Doppler broadening comes before the transform (SAMMY's order): every term's sub-point values are
  // broadened on the common grid, the broadened values transformed, and the kernel then applied without its own
  // Doppler step.  A line with no kernel after the Doppler step gives the transform of sigma_D at the point.
  const TargetEffect *effect = basePoint->GetParentData()->GetTargetEffect(basePoint->GetTargetEffectNum());
  const bool doppler = effect && effect->IsDoppler();
  std::vector<std::vector<double>> raw;
  std::vector<std::pair<int, double>> dw;
  auto broadenAt = [&](double e) {
    basePoint->DopplerWeights(e, dw);
    for (size_t q = 0; q < terms.size(); q++) {
      double s = 0.0;
      for (const auto &t : dw) s += t.second * raw[q][t.first - 1];
      v[q] = s;
    }
  };
  if (doppler) {
    raw.assign(terms.size(), std::vector<double>(nsp));
    for (size_t q = 0; q < terms.size(); q++)
      for (int j = 1; j <= nsp; j++) raw[q][j - 1] = valueOf(terms[q].point->GetSubPoint(j), terms[q].xsComponent);
    if (!effect->HasResolutionStage()) {
      broadenAt(basePoint->GetCMEnergy());
      return transform(v);
    }
  }
  for (int j = 1; j <= nsp; j++) {
    if (doppler) {
      broadenAt(basePoint->GetSubPoint(j)->GetCMEnergy());
    } else {
      for (size_t q = 0; q < terms.size(); q++) v[q] = valueOf(terms[q].point->GetSubPoint(j), terms[q].xsComponent);
    }
    f[j - 1] = transform(v);
  }
  EPoint::DopplerAppliedScope dopplerDone(basePoint);
  return basePoint->IntegrateSubPointValues(configure, f);
}

/*!
 * The combination of a segment at a point (ModelPass.h): the base point's value (its E1 / E2 part for an E1 / E2
 * segment), weighted by its abundance in a sample segment, then the SUM of the scaled components in the base point's
 * frame, or the RATIO over the first component, or the transform.  The per-nucleus terms of a sample segment go into
 * the base point's leaf record (the decomposition output).
 */
double ComposeModelValue(ESegment &segment, int pointIndex, CNuc *cnuc, const Config &configure,
                         azure2::Results &results, azure2::LeafValues leaves) {
  EPoint *basePoint = segment.GetPoint(pointIndex + 1);
  if (CalculatesItself(basePoint, leaves) && !CalculatePoint(basePoint, "base term", cnuc, configure, results))
    return FailedSample(segment, basePoint, cnuc, results);

  if (segment.IsTransformed()) return TransformedValue(segment, pointIndex, cnuc, configure, results, leaves);

  // For E1/E2-only segments compare against the corresponding component of the
  // (angle-integrated capture) cross section instead of the full cross section.
  double baseValue;
  if (segment.GetCrossSectionComponent() == 1)
    baseValue = basePoint->GetFitE1CrossSection(results);
  else if (segment.GetCrossSectionComponent() == 2)
    baseValue = basePoint->GetFitE2CrossSection(results);
  else
    baseValue = basePoint->GetFitCrossSection(results);

  // A sample segment: every term is weighted by the isotopic abundance of its entrance pair's compound nucleus, and
  // the weighted terms are kept per nucleus on the base point for the decomposition output.
  std::vector<double> naturalTerms;
  const bool natural = segment.IsNatural();
  const int sample = segment.GetNaturalSample();
  if (natural) {
    naturalTerms.assign(cnuc->NumNuclei(), 0.0);
    const int nucleus = cnuc->GetPair(cnuc->GetPairNumFromKey(segment.GetEntranceKey()))->GetNucleus();
    baseValue *= cnuc->GetAbundance(sample, nucleus);
    naturalTerms[nucleus - 1] += baseValue;
  }

  if (!segment.HasComponents()) {
    if (natural) basePoint->SetNaturalTerms(results, naturalTerms);
    return baseValue;
  }

  if (segment.GetOperationType() == SUM) {
    double sum = baseValue;
    for (ESegment *componentSegment : segment.GetComponentSegments()) {
      if (!componentSegment || pointIndex >= componentSegment->NumPoints()) continue;
      EPoint *componentPoint = componentSegment->GetPoint(pointIndex + 1);
      if (!componentPoint) continue;
      if (CalculatesItself(componentPoint, leaves) &&
          !CalculatePoint(componentPoint, "SUM component", cnuc, configure, results))
        return FailedSample(segment, basePoint, cnuc, results);
      // A non-finite term makes the sum non-finite (it was left out before 2026-11-01).
      const double compValue = componentPoint->GetFitCrossSection(results);
      // the term in the base point's frame (ESegment::TermFrameFactor)
      double term = componentSegment->GetComponentScaling() * ESegment::TermFrameFactor(basePoint, componentPoint) * compValue;
      if (natural) {
        const int nucleus = cnuc->GetPair(cnuc->GetPairNumFromKey(componentSegment->GetEntranceKey()))->GetNucleus();
        term *= cnuc->GetAbundance(sample, nucleus);
        naturalTerms[nucleus - 1] += term;
      }
      sum += term;
    }
    if (natural) basePoint->SetNaturalTerms(results, naturalTerms);
    return sum;
  }
  if (segment.GetOperationType() != RATIO) return baseValue;

  // RATIO: base / first component (never a sample segment: refused at load).
  ESegment *firstComponent = segment.GetComponentSegments()[0];
  EPoint *denominatorPoint = (firstComponent && pointIndex < firstComponent->NumPoints())
                                 ? firstComponent->GetPoint(pointIndex + 1) : nullptr;
  if (!denominatorPoint) return kNaN;  // component segments copy the base's points: not expected
  if (CalculatesItself(denominatorPoint, leaves) &&
      !CalculatePoint(denominatorPoint, "RATIO denominator", cnuc, configure, results))
    return kNaN;
  const double denominatorValue = denominatorPoint->GetFitCrossSection(results);
  if (denominatorValue == 0.0) {
    // No ratio exists (e.g. the denominator's exit channel is closed at this energy): NaN and a warning, once.
    // Before 2026-11-01 the point was 0, which a fit compared with the datum as if it were a calculated ratio.
    static std::atomic<bool> reported(false);
    if (!reported.exchange(true)) {
#pragma omp critical(azure_failed_point)
      configure.outStream << "WARNING: the RATIO denominator of segment " << basePoint->GetSegmentKey()
                          << " is 0 at E_lab " << basePoint->GetLabEnergy() << " MeV, angle "
                          << basePoint->GetLabAngle() << " deg: the ratio is NaN there (further zero"
                          << " denominators are not reported)." << std::endl;
    }
    return kNaN;
  }
  // A lab-frame differential ratio is the ratio of the two LAB cross sections, sigma_lab = sigma_cm / J: each term
  // with its own pair's Jacobian (ESegment::RatioTermFrameFactor; audit item 15).
  const double baseConversionFactor = ESegment::RatioTermFrameFactor(basePoint, &segment, cnuc);
  const double denominatorConversionFactor = ESegment::RatioTermFrameFactor(denominatorPoint, firstComponent, cnuc);
  return baseValue / denominatorValue * (denominatorConversionFactor / baseConversionFactor);
}

}  // namespace

namespace azure2 {

double EvaluateModelPoint(ESegment &segment, int pointIndex, CNuc *compound, const Config &config, Results &results,
                          LeafValues leaves) {
  EPoint *point = (pointIndex >= 0 && pointIndex < segment.NumPoints()) ? segment.GetPoint(pointIndex + 1) : nullptr;
  if (!point) throw std::logic_error("EvaluateModelPoint: segment " + std::to_string(segment.GetSegmentKey()) +
                                     " has no point " + std::to_string(pointIndex));
  const double m = ComposeModelValue(segment, pointIndex, compound, config, results, leaves);
  point->SetModelValue(results, m);
  point->CopyAngularDistsToModel(results);
  return m;
}

void EvaluateModelSegment(ESegment &segment, CNuc *compound, const Config &config, Results &results,
                          LeafValues leaves) {
  for (int pid = 0; pid < segment.NumPoints(); pid++)
    EvaluateModelPoint(segment, pid, compound, config, results, leaves);
}

void EvaluateModel(EData &data, CNuc *compound, const Config &config, Results &results, LeafValues leaves) {
  for (int i = 1; i <= data.NumSegments(); i++) EvaluateModelSegment(*data.GetSegment(i), compound, config, results, leaves);
}

Residual PointResidual(const EPoint &point, const Results &results, double norm) {
  return DataResidual(point.GetModelValue(results), point.GetCMCrossSection(), point.GetCMCrossSectionError(), norm);
}

}  // namespace azure2

void EPoint::NoModelRecord(const azure2::Results &r) const {
  throw std::logic_error(
      "EPoint: a model value was read or written on a point without a model record in the buffer it was given "
      "(segment " + std::to_string(segment_key_) + ", E_lab " + std::to_string(lab_energy_) + " MeV, offset " +
      (model_slot_ == kNoRecord ? std::string("none") : std::to_string(model_slot_)) + ", " +
      std::to_string(r.model.size()) + " model values): only the points of a data set's own segments have one "
      "(EData::LayoutResults), and EData::BeginEvaluation sizes the buffer");
}
