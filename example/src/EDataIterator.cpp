#include "EDataIterator.h"
#include "ESegment.h"


/*!
 * The EDataIterator object is created with reference to a vector of ESegment objects.
 */

EDataIterator::EDataIterator(std::vector<ESegment>* segments) :
  segments_(segments) {
  segmentIterator_=segments_->begin();
  pointIterator_=(*segmentIterator_).GetPoints().begin();
}

/*!
 * The copy constructor creates a new EDataIterator object with reference to an already
 * existing object.
 */

EDataIterator::EDataIterator(const EDataIterator& it) :
  segments_(it.segments_), segmentIterator_(it.segmentIterator_), pointIterator_(it.pointIterator_) {
}

/*!
 * This function defines the prefix version of the iterating operator.
 */

EDataIterator& EDataIterator::operator++() {
  if(pointIterator_ < ((*segmentIterator_).GetPoints().end())-1) {
    pointIterator_++;
  } else if(segmentIterator_ < segments_->end()-1) {
    segmentIterator_++;
    pointIterator_=(*segmentIterator_).GetPoints().begin();
  } else if ((segmentIterator_ == (segments_->end()-1)) && 
	     (pointIterator_ == ((*segmentIterator_).GetPoints().end()-1)))
    pointIterator_++;
  return *this;
}

/*!
 * This function defines the postfix version of the iterating operator.
 */

EDataIterator EDataIterator::operator++(int) {
  EDataIterator temp(*this);
  operator++();
  return temp;
}

/*!
 * This function defines the boolean equal operator.
 */

bool EDataIterator::operator==(const EDataIterator& rhs) {
  return (rhs.segmentIterator_==segmentIterator_ && rhs.pointIterator_==pointIterator_) ? (true) : (false);
}

/*!
 * This function defines the boolean not-equal operator.
 */

bool EDataIterator::operator!=(const EDataIterator& rhs) {
  return (rhs.segmentIterator_==segmentIterator_ && rhs.pointIterator_==pointIterator_) ? (false) : (true);
}

/*!
 * Sets the contained std::vector<ESegment>::iterator object to the refer to the last segment, and 
 * the contained std::vector<EData>::iterator object to refer to the default end() position (one-past).
 */

EDataIterator& EDataIterator::SetEnd() {
  segmentIterator_=segments_->end()-1;
  pointIterator_=(*segmentIterator_).GetPoints().end();
  return *this;
}

/*!
 * Returns a reference to the contained std::vector<ESegment>::iterator object.
 */

ESegmentIterator& EDataIterator::segment() {
  return segmentIterator_;
}

/*!
 * Returns a reference to the contained std::vector<EPoint>::iterator object.
 */

EPointIterator& EDataIterator::point() {
  return pointIterator_;
}
