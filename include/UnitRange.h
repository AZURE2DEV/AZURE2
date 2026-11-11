#ifndef AZURE2_UNITRANGE_H
#define AZURE2_UNITRANGE_H

/*!
 * A view of a run of calculation units in the working copy's unit vector (EData::Units, docs/upgrade/
 * data_structure.md section 1.4): the points of a term, or the sub-points of a point.  Made on demand from the
 * owner's vector and a (first, count) range, never stored: a copy of the working copy holds ranges, not views.
 */
template <class P>
class UnitRange {
 public:
  UnitRange() = default;
  UnitRange(P *first, int count) : first_(first), count_(count) {}
  P *begin() const { return first_; }
  P *end() const { return first_ + count_; }
  int size() const { return count_; }
  bool empty() const { return count_ == 0; }
  P &operator[](int i) const { return first_[i]; }
  P &front() const { return first_[0]; }
  P &back() const { return first_[count_ - 1]; }

 private:
  P *first_ = nullptr;
  int count_ = 0;
};

#endif
