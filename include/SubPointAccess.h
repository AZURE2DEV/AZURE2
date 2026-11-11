// The sub-point accessors of EPoint, inline: a point's sub-points are a range of the units of its data set
// (EData::Units).  Included at the end of EData.h, where both classes are complete.
#ifndef AZURE2_SUBPOINTACCESS_H
#define AZURE2_SUBPOINTACCESS_H

inline EPoint *EPoint::GetSubPoint(EData &data, int i) {
  return i <= sub_count_ ? &data.Unit(sub_first_ + i - 1) : nullptr;
}
inline const EPoint *EPoint::GetSubPoint(const EData &data, int i) const {
  return i <= sub_count_ ? &data.Unit(sub_first_ + i - 1) : nullptr;
}
inline UnitRange<EPoint> EPoint::GetSubPoints(EData &data) { return data.Units(sub_first_, sub_count_); }

#endif
