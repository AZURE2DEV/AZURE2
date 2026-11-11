// The data set of a project after the set-up (DataSet.h): the lookups that are more than one line.
#include "DataSet.h"

#include <algorithm>

const ESegment *DataSet::GetSegmentFromKey(int key) const {
  for (const ESegment &segment : segments_)
    if (segment.GetSegmentKey() == key) return &segment;
  return nullptr;
}

int DataSet::NormOwner(int s) const {
  if (s < 1 || s > (int)segments_.size()) return s;
  const int lead = segments_[s - 1].GetNormLeaderPos();
  return (lead >= 0) ? lead + 1 : s;
}

int DataSet::ThicknessParamOfSegment(int pos) const {
  if (pos < 1 || pos > (int)segments_.size()) return 0;
  const int lead = segments_[pos - 1].GetThicknessLeaderPos();
  const int owner = (lead >= 0) ? lead + 1 : pos;
  for (int k = 1; k <= this->NumThicknessParams(); k++)
    if (thicknessSegments_[k - 1] == owner) return k;
  return 0;
}

int DataSet::DopplerParamOfEffect(int n) const {
  for (int k = 1; k <= this->NumDopplerParams(); k++)
    for (int c : dopplerCopies_[k - 1])
      if (c == n) return k;
  return 0;
}

bool DataSet::HasRelativeCovariance() const {
  for (const std::shared_ptr<DataCovariance::Group> &gp : covGroups_)
    if (gp->HasRelative()) return true;
  return false;
}

int DataSet::ThmGroupOf(int i) const {
  for (size_t g = 0; g < thmGroups_.size(); g++) {
    if (thmGroups_[g].trivial) continue;
    const std::vector<int> &s = thmGroups_[g].segments;
    if (std::find(s.begin(), s.end(), i) != s.end()) return GetSegment(i)->IsProfiledNorm() ? (int)g : -1;
  }
  return -1;
}

bool DataSet::IsLastOfThmGroup(int g, int i) const {
  const std::vector<int> &s = thmGroups_[g].segments;
  for (size_t k = s.size(); k-- > 0;)
    if (GetSegment(s[k])->IsProfiledNorm()) return s[k] == i;
  return false;
}
