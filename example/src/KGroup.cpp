#include "KGroup.h"

/*!
 * The KGroup is created from a specific combination of entrance and exit channel spin values.
 */

KGroup::KGroup(double s, double sPrime) :
    s_(s), sp_(sPrime) {};

/*!
 * The version of KGroup is overloaded to accomidate the multiple final spin combinations for unobserved primary, observed secondary reactions.
 */

KGroup::KGroup(double s, double sPrime, double sPrime2) :
    s_(s), sp_(sPrime), sp2_(sPrime2) {};

/*!
 * Returns the number of internal reaction pathways in the MGroup vector.
 */

int KGroup::NumMGroups() const {
  return mgroups_.size();
}

/*!
 * Returns the number of external reaction pathways in the ECMGroup vector;
 */

int KGroup::NumECMGroups() const {
  return ec_mgroups_.size();
}

/*!
 * Tests a specific internal reaction pathway to see if it already exists in the MGroup vector. 
 * If the pathway exists, its position in the vector is returned.  Otherwise, the function
 * returns 0.
 */

int KGroup::IsMGroup(MGroup mGroup) {
  bool b=false;
  int c=0;
  while(!b&&c<this->NumMGroups())
    {
      if(mGroup.GetChNum()==this->GetMGroup(c+1)->GetChNum()&&
	 mGroup.GetChpNum()==this->GetMGroup(c+1)->GetChpNum()&&
	 mGroup.GetJNum()==this->GetMGroup(c+1)->GetJNum()) b=true;
      c++;
    }
  if(b) return c;
  else return 0;
}

/*!
 * Returns the value of the entrance channel spin.
 */

double KGroup::GetS() const {
  return s_;
}

/*!
 * Returns the value of the exit channel spin.
 */

double KGroup::GetSp() const {
  return sp_;
}

/*!
 * Returns the value of the second exit channel spin.
 */

double KGroup::GetSp2() const {
  return sp2_;
}

/*!
 * Adds a new internal reaction pathway to the MGroup vector. 
 */

void KGroup::AddMGroup(MGroup mGroup) {
  mgroups_.push_back(mGroup);
}

/*!
 * Adds a new external reaction pathway to the ECMGroup vector.
 */

void KGroup::AddECMGroup(ECMGroup ecMGroup) {
  ec_mgroups_.push_back(ecMGroup);
}

/*!
 * Returns a pointer to the internal reaction pathway specified by a position in the MGroup vector. 
 */

MGroup *KGroup::GetMGroup(int mGroupNum) {
  MGroup *b=&mgroups_[mGroupNum-1];
  return b;
}
/*!
 * Returns a pointer to the external reaction pathway specified by a position in the ECMGroup vector. 
 */

ECMGroup *KGroup::GetECMGroup(int ecMGroupNum) {
  ECMGroup *b=&ec_mgroups_[ecMGroupNum-1];
  return b;
}
