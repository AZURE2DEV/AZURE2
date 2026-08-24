#ifndef ECMGROUP_H
#define ECMGROUP_H

#include "Constants.h"

/// An AZURE external reaction pathway.

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
  /// Hard-sphere external pathway: radiation type, multipolarity, l, J, final
  /// channel, and the J-group and level of the final state.
  ECMGroup(char, int, int, double, int, int, int);
  /// Resonant external (channel-capture) pathway: as above, plus the decay,
  /// KGroup and MGroup of the internal pathway it passes through.
  ECMGroup(char, int, int, double, int, int, int, int, int, int, int);
  /// Is this a resonant external pathway (channel capture) rather than hard sphere?
  bool IsChannelCapture() const;
  /// Radiation type of the capture gamma: 'E' or 'M'.
  char GetRadType() const;
  /// Multipolarity of the capture gamma.
  int GetMult() const;
  /// Entrance orbital angular momentum of the pathway.
  int GetL() const;
  /// 1-based final channel number.
  int GetFinalChannel() const;
  /// 1-based J-group of the final state.
  int GetJGroupNum() const;
  /// 1-based level within that J-group.
  int GetLevelNum() const;
  /// 1-based exit pair in the Decay vector. Channel capture only.
  int GetChanCapDecay() const;
  /// 1-based KGroup of the internal pathway. Channel capture only.
  int GetChanCapKGroup() const;
  /// 1-based MGroup of the internal pathway. Channel capture only.
  int GetChanCapMGroup() const;
  /// 1-based internal channel this external pathway corresponds to.
  int GetIntChannelNum() const;
  /// Entrance spin of the pathway.
  double GetJ() const;
  /// Statistical spin factor \f$g_J\f$ for the pathway.
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
