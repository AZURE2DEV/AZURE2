#ifndef AMATRIXFUNC_H
#define AMATRIXFUNC_H

#include "GenMatrixFunc.h"

///A function class to calculate the T-Matrix using the A-Matrix

/*!
 * The AMatrixFunc function class calculates the T-Matrix for a given energy point using the compound
 * nucleus object.  The AMatrixFunc class is a child class of GenMatrixFunc, where the cross section is
 * calculated from the T-Matrix.
 */

class AMatrixFunc : public GenMatrixFunc {
 public:
  AMatrixFunc(CNuc*, const Config &configure);
  /*!
   * Returns a pointer to the compound nucleus object.
   */
  CNuc *compound() const {return compound_;};
  /*!
   * Returns a reference to the Config structure.
   */
  const Config &configure() const {return configure_;};
  
  void ClearMatrices();
  void FillMatrices(EPoint*);
  void InvertMatrices();
  void CalculateTMatrix(EPoint*);
  /*!
   * Instantiated in the parent class.
   */
  void CalculateCrossSection();

  complex GetAMatrixElement(int,int,int) const;
  matrix_c *GetJSpecAInvMatrix(int);
  void AddAInvMatrixElement(int,int,int,complex);
  void AddAMatrix(matrix_c);
 private:
  const Config &configure_;
  CNuc *compound_;
  vector_matrix_c a_inv_matrices_;
  vector_matrix_c a_matrices_;
};

#endif
