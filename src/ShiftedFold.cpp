// Folded data with an energy shift on the fixed nodes of the shared mesh (include/ShiftedFold.h).
#include "ShiftedFold.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "CalculationPlan.h"
#include "Config.h"
#include "Segment.h"
#include "ModelPass.h"
#include "SegmentState.h"
#include "TargetEffect.h"
#include "WorkingCopy.h"

namespace azure2 {

namespace {
/// A move of the support beyond the planned nodes smaller than this (MeV) is rounding of the sums that place them.
constexpr double kNodeRounding = 1.0e-9;
/// The reach in prior widths (question 1).
constexpr double kReachWidths = 3.0;
/// The halvings of a shift column's difference step before the column is left to the caller's path (question 5):
/// 40 take 1e-6 MeV to 1e-18 MeV.
constexpr int kStepHalvings = 40;

/// The part of a segment's state TotalEnergyShift reads, copied.
SegmentState ShiftPart(const SegmentState &state) {
  SegmentState s;
  s.energyShift = state.energyShift;
  s.energyShiftSqrt = state.energyShiftSqrt;
  s.shiftFn = state.shiftFn;
  s.shiftFnValues = state.shiftFnValues;
  return s;
}
}  // namespace

bool HasEnergyShift(const Segment &segment) {
  return segment.IsVaryEnergyShift() || segment.StartEnergyShift() != 0.0 || segment.GetShiftLeaderPos() >= 0 ||
         segment.HasShiftFn() || segment.IsVaryEnergyShiftSqrt() || segment.StartEnergyShiftSqrt() != 0.0;
}

double PointShift(const Segment &segment, const SegmentState &state, double labEnergy, double cmFactor) {
  return cmFactor * segment.TotalEnergyShift(state, labEnergy);
}

double ShiftReach(const std::vector<Segment> &segments, const Segment &segment, const SegmentState &start, double labEnergy,
                  double cmFactor, double supportWidth) {
  const Segment &constant = segment.GetShiftLeaderPos() >= 0 ? segments[(size_t)segment.GetShiftLeaderPos()] : segment;
  const Segment &function = (segment.HasShiftFn() && segment.GetShiftFnSourcePos() >= 0)
                                ? segments[(size_t)segment.GetShiftFnSourcePos()]
                                : segment;
  const double at = PointShift(segment, start, labEnergy, cmFactor);
  double sum = 0.0;
  // The move of the point's c.m. shift when one parameter moves by its prior width; the support's width without one.
  auto add = [&](double width, void (*apply)(const Segment &, SegmentState &, int, double), int j) {
    if (width == 0.0) {
      sum += supportWidth;
      return;
    }
    SegmentState s = ShiftPart(start);
    apply(segment, s, j, width);
    sum += std::fabs(PointShift(segment, s, labEnergy, cmFactor) - at);
  };
  if (constant.IsVaryEnergyShift())
    add(constant.GetEnergyShiftError(), [](const Segment &, SegmentState &s, int, double w) { s.energyShift += w; }, 0);
  if (segment.IsVaryEnergyShiftSqrt())
    add(segment.GetEnergyShiftSqrtError(), [](const Segment &, SegmentState &s, int, double w) { s.energyShiftSqrt += w; }, 0);
  if (segment.HasShiftFn())
    for (int j = 0; j < segment.NumShiftFnParams() && j < (int)start.shiftFnValues.size(); j++)
      if (function.IsVaryShiftFn(j))
        add(j < (int)function.GetShiftFnErrors().size() ? function.GetShiftFnErrors()[(size_t)j] : 0.0,
            [](const Segment &seg, SegmentState &s, int k, double w) { seg.SetShiftFnValue(s, k, s.shiftFnValues[(size_t)k] + w); },
            j);
  return kReachWidths * sum;
}

double ShiftCentre(double delta, double startDelta, double reach) {
  if (!(reach > 0.0)) return delta;
  return startDelta + 2.0 * reach * std::round((delta - startDelta) / (2.0 * reach));
}

double SupportThicknessChange(const Unit &point, const Unit::SubGridSpec &spec, const TargetEffect &effect,
                              const Config &configure, double delta) {
  if (!spec.targetIntegration || delta == 0.0) return 0.0;
  const double f = spec.stoppingConversion;
  return effect.SupportThickness(point.GetOriginalEnergy() + delta / f, spec.centreEnergy + delta, f, configure) -
         spec.supportThickness;
}

MeshSupport PlannedShiftSupport(const MeshSupport &nominal, const Unit::SubGridSpec &spec, double centre,
                                const std::function<double(double)> &thicknessChange) {
  // The low edge a_i + Delta - dT(Delta) over the reach |Delta - c_i| <= R_i lies above a_i + c_i - R_i - G (question 1):
  // G the largest dT of the three samples c_i - R_i, c_i, c_i + R_i plus their spread.
  const double R = spec.shiftReach;
  double largest = thicknessChange(centre), smallest = largest;
  if (R > 0.0)
    for (double d : {centre - R, centre + R}) {
      const double t = thicknessChange(d);
      largest = std::max(largest, t);
      smallest = std::min(smallest, t);
    }
  const double growth = largest + (largest - smallest);
  MeshSupport s = nominal;
  s.low = nominal.low + centre - spec.shiftReach - growth;
  s.high = nominal.high + centre + spec.shiftReach;
  const double floor = spec.floorEnergy + spec.frameShift;
  if (s.low <= floor) {
    s.low = floor;
    s.lowHard = true;
  }
  return s;
}

bool IsShiftedOnMesh(const Unit &point, const WorkingCopy &data) {
  const Unit::SubGridSpec &spec = data.Data().GridSpecOf(point);
  return spec.mesh >= 0 && (spec.shifted || data.MeshCentre(spec) != 0.0);
}

namespace {
/// The moved support of a shifted point (r.delta, r.low, r.high, clipped to its planned nodes within rounding); the
/// edge ("upper", "lower") that leaves the planned nodes by more than that, nullptr when none does.
const char *MovedSupport(const Unit &point, const WorkingCopy &data, ShiftedNodes &r) {
  const Unit::SubGridSpec &spec = data.Data().GridSpecOf(point);
  r.delta = point.GetCMEnergy() - spec.centreEnergy;
  r.high = spec.startEnergy + r.delta;
  r.low = std::max(spec.endEnergy + r.delta - point.SupportThicknessChange(), spec.floorEnergy);
  const int n = point.NumSubPoints();
  if (n == 0) return nullptr;
  const double top = point.GetSubPoint(data, 1)->GetCMEnergy(), bottom = point.GetSubPoint(data, n)->GetCMEnergy();
  if (n == 1) {
    // A single planned node, the exact upper edge of a planned support narrower than MeshZeroWidth(): the fold is the
    // thin limit at the node whatever the shift (header, "A support narrower than ...").  The moved support lies inside
    // the planned one, so at most MeshZeroWidth() below the node.
    if (r.high - top > kNodeRounding) return "upper";
    return top - r.low > MeshZeroWidth() + kNodeRounding ? "lower" : nullptr;
  }
  if (r.high - top > kNodeRounding) return "upper";
  if (bottom - r.low > kNodeRounding) return "lower";
  r.high = std::min(r.high, top);
  r.low = std::max(r.low, bottom);
  return nullptr;
}
}  // namespace

ShiftedNodes ShiftedFoldNodes(const Unit &point, const WorkingCopy &data) {
  ShiftedNodes r;
  const int n = point.NumSubPoints();
  auto energy = [&](int j) { return point.GetSubPoint(data, j)->GetCMEnergy(); };  // j = 1 .. n, descending
  if (const char *edge = MovedSupport(point, data, r))
    throw std::logic_error(std::string("shifted fold: the ") + edge + " edge of the moved support of segment " +
                           std::to_string(point.GetSegmentKey()) + " at E_cm " +
                           std::to_string(edge[0] == 'u' ? r.high : r.low) + " MeV is outside the point's planned nodes [" +
                           std::to_string(energy(n)) + ", " + std::to_string(energy(1)) + "] (ShiftedFold.h, question 2)");
  if (n == 0) return r;
  if (n == 1) {
    r.count = 1;
    return r;
  }
  int first = 1, last = n;
  while (first < n && energy(first + 1) >= r.high) first++;
  while (last > first && energy(last - 1) <= r.low) last--;
  if (last == first && n >= 2) (last < n) ? last++ : first--;  // at least two nodes: a kernel needs an interval
  r.first = first - 1;
  r.count = last - first + 1;
  return r;
}

bool ShiftedSupportsInside(const WorkingCopy &data, const std::vector<const Segment *> &segments) {
  ShiftedNodes r;
  for (const Segment *segment : segments)
    for (int t = 0; t < segment->NumTerms(); t++)
      for (const Unit &point : segment->TermPoints(data, t))
        if (IsShiftedOnMesh(point, data) && MovedSupport(point, data, r)) return false;
  return true;
}

double InsideShiftStep(const WorkingCopy &data, const std::vector<int> &segmentNumbers, double x0, double h,
                       const std::function<void(double)> &move) {
  std::vector<const Segment *> segments;
  for (int s : segmentNumbers) segments.push_back(data.Data().GetSegment(s));
  double step = 0.0;
  for (int k = 0; k < kStepHalvings && step == 0.0; k++, h *= 0.5) {
    move(x0 + h);
    if (!ShiftedSupportsInside(data, segments)) continue;
    move(x0 - h);
    if (ShiftedSupportsInside(data, segments)) step = h;
  }
  move(x0);
  return step;
}

void MoveTargetToShift(Unit &point, const WorkingCopy &data, const Config &configure) {
  const Unit::SubGridSpec &spec = data.Data().GridSpecOf(point);
  if (!spec.valid || !spec.targetIntegration) return;
  const TargetEffect &effect = *data.Data().GetTargetEffect(point.GetTargetEffectNum());
  const TargetEffectValues &values = data.EffectValues(point.GetTargetEffectNum());
  const double f = spec.stoppingConversion, lab = point.GetLabEnergy(), cm = point.GetCMEnergy();
  // the lab energy a fitted density or stopping scale recomputes the thickness at (WorkingCopy::FillDensitiesFromParams,
  // FollowStoppingScales)
  point.SetTargetLabEnergy(lab);
  point.SetTargetThickness(effect.KernelThickness(lab, cm, f, configure, values.density, values.stoppingScale));
  point.SetSupportThicknessChange(
      cm == spec.centreEnergy ? 0.0 : effect.SupportThickness(lab, cm, f, configure) - spec.supportThickness);
}

bool FoldOnlyShift(const WorkingCopy &data, const Segment &segment) {
  if (!segment.IsTargetEffect()) return false;
  // The node values do not depend on the shift; neither does the row set of a point given at a lab angle (its c.m.
  // angle at a node is a function of the node's energy and the lab angle, RowSet.h), only the point's own Jacobian.
  for (int t = 0; t < segment.NumTerms(); t++)
    for (const Unit &point : segment.TermPoints(data, t))
      if (!IsShiftedOnMesh(point, data) || !HasSubPointFold(data, point)) return false;
  return true;
}

bool RefoldShiftedSegments(WorkingCopy &data, const std::vector<const Segment *> &segments, const Config &configure,
                           Results &results) {
  for (const Segment *segment : segments)
    if (!FoldOnlyShift(data, *segment)) return false;
  for (const Segment *segment : segments)
    for (int t = 0; t < segment->NumTerms(); t++)
      for (Unit &point : segment->TermPoints(data, t)) point.FoldSubPoints(data, configure, results);
  return true;
}

double SegmentLocalChi2(const Segment &segment, const CNuc *compound, const LevelState &state, const Config &config, WorkingCopy &data,
                        Results &results, double norm) {
  double chi = 0.0;
  const bool refolded = RefoldShiftedSegments(data, {&segment}, config, results);
  for (int pid = 0; pid < segment.NumPoints(); pid++) {
    EvaluateModelPoint(data, segment, pid, compound, state, config, results,
                       refolded ? LeafValues::Ready : LeafValues::Calculate);
    chi += PointResidual(data, segment, pid, results, norm).Chi2();
  }
  segment.AddShiftPenalty(data.State(segment), chi);
  return chi;
}

std::vector<double> HeldNorms(WorkingCopy &data, const std::vector<int> &segments, const CNuc *compound, const LevelState &state, const Config &config,
                              Results &results) {
  std::vector<double> norms;
  for (int s : segments) {
    const Segment &segment = *data.Data().GetSegment(s);
    double norm = data.State(segment).norm;
    if (segment.IsProfiledNorm()) {
      SegmentLocalChi2(segment, compound, state, config, data, results, norm);  // the models at the present shift
      segment.ProfileNormChiSquared(data, results, norm);
    }
    norms.push_back(norm);
  }
  return norms;
}

}  // namespace azure2
