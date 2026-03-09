#include "KLGroup.h"

/*!
 * The object is created with reference to a specfic KGroup number as well as Legendre polynomial order.
 */

KLGroup::KLGroup(int kGroupNum,int lOrder) :
  k_(kGroupNum), lorder_(lOrder) {};

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
  bool b=false;
  int c=0;
  while(!b&&c<this->NumInterferences())
    {
      if(interference.GetM1()==this->GetInterference(c+1)->GetM1()&&
	 interference.GetM2()==this->GetInterference(c+1)->GetM2()&&
	 interference.GetInterferenceType()==this->GetInterference(c+1)->GetInterferenceType()) b=true;
	 c++;
    }
  if(b) return c;
  else return 0;
}

/*!
 * Adds an interference combination to the Interference vector. 
 */

void KLGroup::AddInterference(Interference interference) {
  interferences_.push_back(interference);
}

/*!
 * Returns a pointer to an interference combination specified by a position in the Interference vector.
 */

Interference *KLGroup::GetInterference(int interferenceNum) {
  Interference *b=&interferences_[interferenceNum-1];
  return b;
}
