#ifndef EDATAITERATOR_H
#define EDATAITERATOR_H

#include <vector>

class EData;
class EPoint;
class ESegment;

typedef EPoint *EPointIterator;
typedef std::vector<ESegment>::iterator ESegmentIterator;

/// An iterator class for an EData object.

/*!
 * The EDataIterator class is used to loop through all data segments and points contained in an
 * EData object.
 */

class EDataIterator {
 public:
  /// Over the points of the segments of `data` (the segments without components).
  explicit EDataIterator(EData *data);
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
