#ifndef JGROUP_H
#define JGROUP_H

#include "Constants.h"
#include "ALevel.h"
#include "AChannel.h"

class NucLine;

/// An AZURE \f$ J^\pi \f$ group.

/*!
 * In R-Matrix theory, levels are grouped according to their \f$ J^\pi \f$ values.  There is one R-/A-Matrix,
 * and thus one T-Matrix, for each \f$ J^\pi \f$ group. A JGroup object holds vectors of ALevel and AChannel objects.
 */

class JGroup {
 public:
  /// Build from a line of the nuclear input file.
  JGroup(NucLine);
  /// Build directly from a spin and a parity (\f$\pm1\f$).
  JGroup(double, int);
  /// Is this group part of the A-/R-Matrix calculation?
  /// False for a group that exists only to carry a bound state for external
  /// capture -- a subthreshold state has no R-Matrix level of its own.
  bool IsInRMatrix() const;
  /// Puts a group built from a spin and parity into the R-matrix calculation: a hard-sphere group
  /// (channels, no levels) scatters through its hard-sphere phases alone, since its level sum is empty.
  void SetInRMatrix(bool inRMatrix) { isinrmatrix_ = inRMatrix; }
  /// 1-based position of the first level in the group whose energy lies within tol (MeV) of the given
  /// level's, or 0 if there is none.  The default 1 keV is for matching a state typed separately in the
  /// particle-pairs tab (a capture final state) to a level; the rows of one level in the <levels> block
  /// carry the same energy and are grouped with tol = 0 (CNuc::Fill).
  int IsLevel(ALevel, double tol = 1e-3);
  /// Parity of the group, \f$\pm1\f$.
  int GetPi() const;
  /// Number of levels in the group.
  int NumLevels() const;
  /// Number of channels shared by every level of the group.
  int NumChannels();
  /// 1-based position of the channel in the group, or 0 if it is not there.
  int IsChannel(AChannel);
  /// Total angular momentum of the group.
  double GetJ() const;
  /// Append a level.
  void AddLevel(ALevel);
  /// Append a channel; it applies to every level in the group.
  void AddChannel(AChannel);
  /// Channel \p i, 1-based. Every level of the group shares this channel set.
  AChannel *GetChannel(int);
  /// Level \p i, 1-based, numbered within this group rather than globally.
  ALevel *GetLevel(int);

 private:
  bool isinrmatrix_;
  int pi_;
  double j_;
  std::vector<ALevel> levels_;
  std::vector<AChannel> channels_;
};

#endif
