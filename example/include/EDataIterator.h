#ifndef EDATAITERATOR_H
#define EDATAITERATOR_H

#include <vector>

class EPoint;
class ESegment;

typedef std::vector<EPoint>::iterator EPointIterator;
typedef std::vector<EPoint*>::iterator EPointMapIterator;
typedef std::vector<ESegment>::iterator ESegmentIterator;

///An iterator class for an EData object.

/*!
 * The EDataIterator class is used to loop through all data segments and points contained in an
 * EData object.
 */

class EDataIterator {
 public:
  EDataIterator(std::vector<ESegment>*);
  EDataIterator(const EDataIterator& it);
  EDataIterator& operator++();
  EDataIterator operator++(int);
  bool operator==(const EDataIterator&);
  bool operator!=(const EDataIterator&);
  EDataIterator& SetEnd();
  ESegmentIterator& segment();
  EPointIterator& point();
  
 private:
  std::vector<ESegment>* segments_;
  ESegmentIterator segmentIterator_;
  EPointIterator pointIterator_;
};

#endif
