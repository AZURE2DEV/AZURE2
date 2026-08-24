#ifndef RMATRIXFUNC_H
#define RMATRIXFUNC_H

#include "GenMatrixFunc.h"

/// A function class to calculate the T-Matrix using the R-Matrix

/*!
 * The RMatrixFunc function class calculates the T-Matrix for a given energy point using the compound
 * nucleus object.  The RMatrixFunc class is a child class of GenMatrixFunc, where the cross section is
 * calculated from the T-Matrix.
 */

class RMatrixFunc : public GenMatrixFunc {
 public:
  /// Build against a compound nucleus.
  RMatrixFunc(CNuc *, const Config &);
  /*!
   * Returns a pointer to the compound nucleus object.
   */
  CNuc *compound() const { return compound_; };
  const Config &configure() const { return configure_; };

  /// Drop every stored matrix before the next point.
  void ClearMatrices();
  /// Build \f$R\f$ and \f$[1-RL]\f$ from the level parameters at this point.
  void FillMatrices(EPoint *);
  /// Invert \f$[1-RL]\f$ and form \f$[1-RL]^{-1}R\f$.
  void InvertMatrices();
  /// T-matrix per reaction pathway from \f$[1-RL]^{-1}R\f$.
  void CalculateTMatrix(EPoint *);
  /*!
   * Instantiated in the parent class.
   */
  void CalculateCrossSection();

  /// \f$R\f$ element at (J-group, channel, channel'), all 1-based.
  complex GetRMatrixElement(int, int, int) const;
  /// \f$[1-RL]\f$ element at (J-group, channel, channel').
  complex GetRLMatrixElement(int, int, int) const;
  /// \f$[1-RL]^{-1}\f$ element at (J-group, channel, channel').
  complex GetRLInvMatrixElement(int, int, int) const;
  /// \f$[1-RL]^{-1}R\f$ element at (J-group, channel, channel').
  complex GetRLInvRMatrixElement(int, int, int) const;
  /// The whole \f$[1-RL]\f$ matrix for one J-group.
  matrix_c *GetJSpecRLMatrix(int);
  /// Store an \f$R\f$ element.
  void AddRMatrixElement(int, int, int, complex);
  /// Store a \f$[1-RL]\f$ element.
  void AddRLMatrixElement(int, int, int, complex);
  /// Append a whole inverted matrix.
  void AddRLInvMatrix(matrix_c);
  /// Store a \f$[1-RL]^{-1}R\f$ element.
  void AddRLInvRMatrixElement(int, int, int, complex);

 private:
  CNuc *compound_;
  const Config &configure_;
  vector_matrix_c r_matrices_;
  vector_matrix_c rl_matrices_;
  vector_matrix_c rl_inv_matrices_;
  vector_matrix_c rl_inv_r_matrices_;
};

#endif
