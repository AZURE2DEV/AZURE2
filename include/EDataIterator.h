#ifndef EDATAITERATOR_H
#define EDATAITERATOR_H

#include <vector>

class EData;  // EData.h
class EPoint;
class ESegment;

typedef EPoint *EPointIterator;
typedef std::vector<ESegment>::iterator ESegmentIterator;

/// Loops through the segments of an EData object and their points.

class EDataIterator {
 public:
  explicit EDataIterator(EData *data);  // the points of the segments of `data` that have no components
  EDataIterator(const EDataIterator &it);
  EDataIterator &operator++();
  EDataIterator operator++(int);
  bool operator==(const EDataIterator &);
  bool operator!=(const EDataIterator &);
  EDataIterator &SetEnd();
  ESegmentIterator &segment();
  EPointIterator &point();

 private:
  EData *data_;
  std::vector<ESegment> *segments_;
  ESegmentIterator segmentIterator_;
  EPointIterator pointIterator_;
};

#endif
