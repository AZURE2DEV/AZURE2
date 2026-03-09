#include "Decay.h"

/*!
 * The Decay object is created with a reference to a particle pair number, which represents
 * a position in the PPair vector.  
 */

Decay::Decay(int pairNum) :
  pair_(pairNum) {};

/*!
 * Returns the pair number of the decay.
 */

int Decay::GetPairNum() const {
  return pair_;
}

/*!
 * Returns the number of \f$ s,s' \f$ combinations in the KGroup vector.
 */

int Decay::NumKGroups() const {
  return kgroups_.size();
}

/*!
 * Returns the number of \f$ k,L \f$ combinations in the KLGroup vector.
 */

int Decay::NumKLGroups() const {
  return klgroups_.size();
}

/*!
 * Tests a specific \f$ s,s' \f$ combination to determine if it exists in the KGroup vector.
 * If the combination exists, the position of the combination in the vector is returned.  Otherwise, the function returns 0.
 */

int Decay::IsKGroup(KGroup a) {
  bool b=false;
  int c=0;
  while(!b&&c<this->NumKGroups())
    {
      if(a.GetS()==this->GetKGroup(c+1)->GetS()&&
	 a.GetSp()==this->GetKGroup(c+1)->GetSp()) b=true;
	 c++;
    }
  if(b) return c;
  else return 0;
}

/*!
 * Tests a specific \f$ s,s',s'2 \f$ combination to determine if it exists in the KGroup vector.
 * If the combination exists, the position of the combination in the vector is returned.  Otherwise, the function returns 0.
 */

int Decay::IsKGroup(KGroup a, bool flag) {
  bool b=false;
  int c=0;
  while(!b&&c<this->NumKGroups())
    {
      if(a.GetS()==this->GetKGroup(c+1)->GetS()&&
	 a.GetSp()==this->GetKGroup(c+1)->GetSp()&&
         a.GetSp2()==this->GetKGroup(c+1)->GetSp2()) b=true;
	 c++;
    }
  if(b) return c;
  else return 0;
}

/*!
 * Tests a specific \f$ k,L \f$ combination to determine if it exists in the KLGroup vector.
 * If the combination exists, the position of the combination in the vector is returned.  Otherwise, the function returns 0.
 */

int Decay::IsKLGroup(KLGroup a) {
  bool b=false;
  int c=0;
  while(!b&&c<this->NumKLGroups())
    {
      if(a.GetK()==this->GetKLGroup(c+1)->GetK()&&
	 a.GetLOrder()==this->GetKLGroup(c+1)->GetLOrder()) b=true;
	 c++;
    }
  if(b) return c;
  else return 0;
}

/*!
 * Adds a \f$ s,s' \f$ combination to the KGroup vector.
 */

void Decay::AddKGroup(KGroup kGroup) {
  kgroups_.push_back(kGroup);
}

/*!
 * Adds a \f$ k,L \f$ combination to the KLGroup vector.
 */

void Decay::AddKLGroup(KLGroup klGroup) {
  klgroups_.push_back(klGroup);
}

/*!
 * Returns a pointer to a \f$ s,s' \f$ combination specified by a position in the KGroup vector.
 */

KGroup *Decay::GetKGroup(int kGroupNum) {
  KGroup *b=&kgroups_[kGroupNum-1];
  return b;
}

/*!
 * Returns a pointer to a \f$ k,L \f$ combination specified by a position in the KLGroup vector.
 */

KLGroup *Decay::GetKLGroup(int klGroupNum) {
  KLGroup *b=&klgroups_[klGroupNum-1];
  return b;
}
