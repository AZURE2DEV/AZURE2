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
  /// 1-based position of the level in the group, or 0 if it is not there.
  int IsLevel(ALevel);
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
  /// The compound nucleus the group belongs to: its 1-based index in the
  /// CNuc's nucleus list (0 until set) and its charge and mass numbers.  A
  /// group holds channels of one nucleus only, so the level matrix of a
  /// project with several nuclei is block-diagonal by construction.
  int GetNucleus() const { return nucleus_; }
  int GetNucleusZ() const { return nucleusZ_; }
  int GetNucleusA() const { return nucleusA_; }
  void SetNucleus(int index, int z, int a) { nucleus_ = index; nucleusZ_ = z; nucleusA_ = a; }
  /// Whether labels of this group's levels should name the nucleus: set by
  /// CNuc::Fill on every group when the project holds more than one nucleus.
  bool NucleusLabelled() const { return nucleusLabelled_; }
  void SetNucleusLabelled(bool on) { nucleusLabelled_ = on; }
  /// B_c = -l_c in every particle channel of the group (Hale's fixed boundary
  /// values), instead of the shift function at the first level's energy.
  /// Set for groups that hold isospin levels, so that isospin partners share B.
  bool BoundaryMinusL() const { return boundaryMinusL_; }
  void SetBoundaryMinusL(bool on) { boundaryMinusL_ = on; }
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
  int nucleus_ = 0;
  int nucleusZ_ = 0;
  int nucleusA_ = 0;
  bool nucleusLabelled_ = false;
  bool boundaryMinusL_ = false;
  double j_;
  std::vector<ALevel> levels_;
  std::vector<AChannel> channels_;
};

#endif
