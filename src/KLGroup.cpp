#include "KLGroup.h"

/*!
 * The object is created with reference to a specfic KGroup number as well as Legendre polynomial order.
 */

KLGroup::KLGroup(int kGroupNum, int lOrder) :
  k_(kGroupNum),
  lorder_(lOrder) {};

/*!
 * Returns the position of the \f$ s,s' \f$ combination in the KGroup vector.
 */

int KLGroup::GetK() const {
  return k_;
}

/*!
 * Returns the Legendre polynomial order.
 */

int KLGroup::GetLOrder() const {
  return lorder_;
}

/*!
 *  Returns the number of interference combinations in the Interference vector.
 */

int KLGroup::NumInterferences() const {
  return interferences_.size();
}

/*!
 * Tests an interference combination to determine if it exists in the Interference vector.
 * If the combination exists, its position in the vector is returned.  Otherwise, the function returns 0.
 */

int KLGroup::IsInterference(Interference interference) {
  std::tuple<int, int, std::string> key(interference.GetM1(), interference.GetM2(),
                                        interference.GetInterferenceType());
  std::map<std::tuple<int, int, std::string>, int>::const_iterator it = interferenceIndex_.find(key);
  if (it != interferenceIndex_.end()) return it->second;
  return 0;
}

/*!
 * Adds an interference combination to the Interference vector.
 */

void KLGroup::AddInterference(Interference interference) {
  interferences_.push_back(interference);
  std::tuple<int, int, std::string> key(interference.GetM1(), interference.GetM2(),
                                        interference.GetInterferenceType());
  // Store the 1-based position, matching the previous IsInterference() return value.
  interferenceIndex_[key] = (int)interferences_.size();
}

/*!
 * Returns a pointer to an interference combination specified by a position in the Interference vector.
 */

Interference *KLGroup::GetInterference(int interferenceNum) {
  Interference *b = &interferences_[interferenceNum - 1];
  return b;
}
