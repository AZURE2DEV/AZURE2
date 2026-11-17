#ifndef SHIFTPARAMTARGET_H
#define SHIFTPARAMTARGET_H

// What one energy-shift parameter moves, for the two places that difference a shift on its own segments: the
// energy-shift columns of the residual Jacobian (FillShiftColumns, src/AZURECalcJacobian.cpp) and the
// segment-local shift derivative of the gradient (AZURECalc::GradientBy, src/AZURECalcGradient.cpp).  A private
// header of those two files (it was in the unnamed namespace of AZURECalc.cpp while they were one file).

#include <vector>
#include "AZUREGrad.h"
#include "Segment.h"
#include "WorkingCopy.h"

// The segments (1-based) whose points one energy-shift parameter moves, and how to
// set it on them: its own segment, plus every follower of its shift group for a
// group's constant shift and for each shift-function parameter.
struct ShiftParamTarget {
  ParamKind kind;
  int fnTerm;
  std::vector<int> segs;
  double Get(WorkingCopy *ld) const {
    const Segment *seg = ld->Data().GetSegment(segs.front());
    if (kind == ParamKind::EnergyShiftSqrt) return ld->State(*seg).energyShiftSqrt;
    if (kind == ParamKind::EnergyShiftFn) return ld->State(*seg).shiftFnValues[fnTerm];
    return ld->State(*seg).energyShift;
  }
  void Set(WorkingCopy *ld, double v) const {
    for (int s : segs) {
      const Segment *seg = ld->Data().GetSegment(s);
      if (kind == ParamKind::EnergyShiftSqrt) ld->State(*seg).energyShiftSqrt = v;
      else if (kind == ParamKind::EnergyShiftFn) seg->SetShiftFnValue(ld->State(*seg), fnTerm, v);
      else ld->State(*seg).energyShift = v;
    }
  }
  // Every member self-contained (no components, not a total-capture sum), so the
  // parameter moves these segments' rows and nothing else.
  bool Local(WorkingCopy *ld) const {
    for (int s : segs) {
      const Segment *seg = ld->Data().GetSegment(s);
      if (seg->HasComponents() || seg->IsTotalCapture()) return false;
    }
    return true;
  }
};

/// The target of the shift parameter `d` on the working copy `ld` (src/AZURECalcJacobian.cpp).
ShiftParamTarget ShiftTargets(WorkingCopy *ld, const ParamDesc &d);

#endif
