#ifndef DECAY_H
#define DECAY_H

#include "KGroup.h"
#include "KLGroup.h"
#include "CaptureAyTerm.h"

/// An AZURE decay pair.

/*!
 * In AZURE, a Decay object represents a decay pair of the compound nucleus.  The Decay object is keyed to a particle pair in the PPair vector,
 * and serves as a container class for the KGroup and the KLGroup vectors and their subsequent reaction pathways.
 */

class Decay {
 public:
  /// Build for the exit pair with this 1-based pair number.
  Decay(int);
  /// 1-based number of the exit particle pair this decay leads to.
  int GetPairNum() const;
  /// Number of \f$s,s'\f$ channel-spin combinations for this decay.
  int NumKGroups() const;
  /// Number of \f$k,L\f$ combinations contributing to the angular distribution.
  int NumKLGroups() const;
  /// 1-based position of the \f$s,s'\f$ combination, or 0 if absent.
  int IsKGroup(KGroup);
  /// As above, but also matching sp2 -- needed for the UPOS groups the
  /// analyzing power adds, where \f$s'\f$ alone does not identify the group.
  int IsKGroup(KGroup, bool);
  /// 1-based position of the \f$k,L\f$ combination, or 0 if absent.
  int IsKLGroup(KLGroup);
  void AddKGroup(KGroup);
  void AddKLGroup(KLGroup);
  /// \f$s,s'\f$ group \p i, 1-based.
  KGroup *GetKGroup(int);
  /// \f$k,L\f$ group \p i, 1-based.
  KLGroup *GetKLGroup(int);

  //! Capture analyzing-power terms, filled by CNuc::CalcCaptureAnalyzingPower.
  int NumCaptureAyTerms() const;
  /// Add one \f$(t,t')\f$ pathway pair to the capture analyzing-power table.
  void AddCaptureAyTerm(const CaptureAyTerm &);
  /// Term \p i, 1-based.
  const CaptureAyTerm *GetCaptureAyTerm(int) const;
  /// Has the table been built? It is filled on first use, not at
  /// initialization: it costs 9-j symbols a run without polarization segments
  /// would never look at.
  bool IsCaptureAyBuilt() const;
  void SetCaptureAyBuilt();

 private:
  int pair_;
  std::vector<KGroup> kgroups_;
  std::vector<KLGroup> klgroups_;
  std::vector<CaptureAyTerm> captureAy_;
  bool captureAyBuilt_ = false;
};


#endif
