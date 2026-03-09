#ifndef ECMGROUP_H
#define ECMGROUP_H

#include "Constants.h"

///An AZURE external reaction pathway.

/*!
 * An external reaction pathways in AZURE is one of two types: a hard sphere pathway or a resonant pathway.
 * The hard sphere pathways refers to the portion of the initial incoming plus outgoing scattering 
 * wavefunction that is hard-shpere scattered and captured directly to a final state.  The resonant pathways
 * refer to the portion of the outgoing scattering wavefunction that is scattered/transformed by the R-Matrix.
 * This type of pathway is linked to the internal resonant pathways, and can be thought of as first passing through
 * a the resonant T-matrix before being captured directly to a final state.
 */

class ECMGroup {
 public:
  ECMGroup(char, int, int, double, int, int, int);
  ECMGroup(char, int, int, double, int, int, int, int, int, int, int);
  bool IsChannelCapture() const;
  char GetRadType() const;
  int GetMult() const;
  int GetL() const;
  int GetFinalChannel() const;
  int GetJGroupNum() const;
  int GetLevelNum() const;
  int GetChanCapDecay() const;
  int GetChanCapKGroup() const;
  int GetChanCapMGroup() const;
  int GetIntChannelNum() const;
  double GetJ() const;
  double GetStatSpinFactor() const;
  void SetStatSpinFactor(double);
 private:
  char radtype_;
  int mult_;
  int li_;
  int chf_;
  int jGroupNum_;
  int levelNum_;
  bool ischancap_;
  int chdecay_;
  int chkgroup_;
  int chmgroup_;
  int internalChannel_;
  double ji_;
  double statspinfactor_;
  complex tmatrix_;
};

#endif
