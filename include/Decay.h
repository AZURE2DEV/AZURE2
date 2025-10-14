#ifndef DECAY_H
#define DECAY_H

#include "KGroup.h"
#include "KLGroup.h"

///An AZURE decay pair.

/*!
 * In AZURE, a Decay object represents a decay pair of the compound nucleus.  The Decay object is keyed to a particle pair in the PPair vector,
 * and serves as a container class for the KGroup and the KLGroup vectors and their subsequent reaction pathways.  
 */

class Decay {
 public:
  Decay(int);
  int GetPairNum() const;
  int NumKGroups() const;
  int NumKLGroups() const;
  int IsKGroup(KGroup);
  int IsKGroup(KGroup,bool);
  int IsKLGroup(KLGroup) ; 
  void AddKGroup(KGroup);
  void AddKLGroup(KLGroup);
  KGroup *GetKGroup(int);
  KLGroup *GetKLGroup(int);
 private:
  int pair_;
  std::vector<KGroup> kgroups_;
  std::vector<KLGroup> klgroups_;
};


#endif
