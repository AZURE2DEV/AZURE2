// The sub-point accessors of Unit, inline: a point's sub-points are a range of the units of its data set
// (WorkingCopy::Units).  Included at the end of WorkingCopy.h, where both classes are complete.
#ifndef AZURE2_SUBPOINTACCESS_H
#define AZURE2_SUBPOINTACCESS_H

inline Unit *Unit::GetSubPoint(WorkingCopy &data, int i) {
  return i <= sub_count_ ? &data.UnitAt(sub_first_ + i - 1) : nullptr;
}
inline const Unit *Unit::GetSubPoint(const WorkingCopy &data, int i) const {
  return i <= sub_count_ ? &data.UnitAt(sub_first_ + i - 1) : nullptr;
}
inline UnitRange<Unit> Unit::GetSubPoints(WorkingCopy &data) { return data.Units(sub_first_, sub_count_); }

#endif
