#ifndef ALEVEL_H
#define ALEVEL_H

#include "Constants.h"

class NucLine;

/// An AZURE level object.

/*!
 * An R-matrix level  represents a specific eigenstate of the compound nucleus.
 *
 * A level carries its energy and one reduced width amplitude per channel of
 * its \f$J^\pi\f$ group, in three parallel sets that are easy to confuse:
 *
 *  - **input / formal** -- GetE, GetGamma: what the nuclear input file said,
 *    in the formal R-matrix parameterization.
 *  - **fitted** -- GetFitE, GetFitGamma: the working values the minimizer
 *    moves. This is the set a calculation actually reads.
 *  - **physical / transformed** -- GetTransformE, GetTransformGamma and
 *    GetBigGamma: the observed energy, reduced width amplitude and partial
 *    width \f$\Gamma\f$ that CNuc::TransformOut produces for reporting.
 *
 * Channel numbers are 1-based and index the channels of the owning JGroup, so
 * the same index means the same channel for every level in the group.
 */

class ALevel {
 public:
  /// Build from a line of the nuclear input file.
  ALevel(NucLine);
  /// Build with only an energy; widths are added afterwards with AddGamma.
  ALevel(double);
  /// Is this level part of the A-/R-Matrix calculation?
  /// False for a level that exists only as an external-capture final state,
  /// such as a subthreshold bound state.
  bool IsInRMatrix() const;
  /// Is the level energy held fixed during fitting?
  bool EnergyFixed() const;
  /// Is this channel's reduced width amplitude held fixed during fitting?
  bool ChannelFixed(int) const;
  /// Is the level a final state for external capture?
  bool IsECLevel() const;
  int NumNFIntegrals() const;
  /// Iterations CNuc::TransformOut needed to converge for this level.
  int GetTransformIterations() const;
  /// 1-based pair number this level is the bound state of; 0 if it is not one.
  int GetECPairNum() const;
  /// Bit mask of the external-capture multipolarities feeding this level.
  unsigned char GetECMultMask() const;

  // -- input / formal parameters ---------------------------------------------
  /// Level energy as read from the input file (MeV).
  double GetE() const;
  /// Formal internal reduced width amplitude for a channel.
  double GetGamma(int) const;

  // -- fitted parameters -----------------------------------------------------
  /// Working reduced width amplitude the minimizer moves; what a calculation reads.
  double GetFitGamma(int) const;
  /// Working level energy the minimizer moves.
  double GetFitE() const;

  // -- physical parameters, from TransformOut --------------------------------
  /// Observed reduced width amplitude after the formal-to-physical transformation.
  double GetTransformGamma(int) const;
  /// Observed level energy after the formal-to-physical transformation.
  double GetTransformE() const;
  /// Observed partial width \f$\Gamma\f$ for a channel (Breit-Wigner).
  double GetBigGamma(int) const;

  // -- external capture ------------------------------------------------------
  /// Channel integral in the denominator of \f$N_f^{1/2}\f$,
  /// \f$\int_a^\infty [W_c(kr)/W_c(ka_c)]^2\f$.
  double GetNFIntegral(int) const;
  /// The \f$N_f^{1/2}\f$ normalization of the bound-state wave function.
  double GetSqrtNFFactor() const;
  /// Factor converting this channel's reduced width amplitude to an ANC.
  double GetECConversionFactor(int) const;
  /// External part of the reduced width amplitude for a channel.
  complex GetExternalGamma(int) const;

  /// Shift function for a channel, evaluated at the resonance energy.
  double GetShiftFunction(int) const;

  /// Append a channel, taking its initial width from the input file line.
  void AddGamma(NucLine);
  /// Append a channel with the given initial reduced width amplitude.
  void AddGamma(double);
  void SetGamma(int, double);
  void SetE(double);
  void SetFitGamma(int, double);
  void SetFitE(double);
  void AddNFIntegral(double);
  void SetSqrtNFFactor(double);
  void AddECConversionFactor(double);
  void SetTransformGamma(int, double);
  void SetTransformE(double);
  void SetBigGamma(int, double);
  void SetTransformIterations(int);
  void SetExternalGamma(int, complex);
  void SetShiftFunction(int, double);
  /// Mark the level as an external-capture final state of a pair, with its
  /// multipolarity mask.
  void SetECParams(int, unsigned char);

 private:
  bool isinrmatrix_;
  bool energyfixed_;
  bool isECLevel_;
  int transform_iter_;
  int ecPairNum_;
  unsigned char ecMultMask_;
  double level_e_;
  double fitlevel_e_;
  double sqrt_nf_factor_;
  double transform_e_;
  std::vector<bool> channelfixed_;
  vector_r gammas_;
  vector_r fitgammas_;
  vector_r nf_integrals_;
  vector_r ec_conv_factors_;
  vector_r transform_gammas_;
  vector_r big_gammas_;
  vector_r shifts_;
  vector_c external_gammas_;
};

#endif
