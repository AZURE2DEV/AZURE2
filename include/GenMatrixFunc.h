#ifndef GENMATRIXFUNC_H
#define GENMATRIXFUNC_H

#include "Constants.h"
#include "Config.h"

class EPoint;
class CNuc;

/// A temporaray T-Matrix structure

/*!
 * The TempTMatrix structure is used to coherently add T-matrix elements from
 * pathways  with like \f$ J,l,l' \f$ values for the calculation of angle integrated cross section.
 * This is primarly used to facilitate the interference between internal and external pathways.
 */

struct TempTMatrix {
  /// Total spin value of temporary matrix element
  double jValue;
  /// Entrance orbital angular momentum for temporary matrix element
  int lValue;
  /// Exit orbital angular momentum for temporary matrix element
  int lpValue;
  /// Value of temporary matrix element
  complex TMatrix;
};


/// A generalized function class to calculate cross sections

/*!
 * The GenMatrixFunc function class is the general form of the function used
 * to calculate cross section from R-Matrix parameters.  It is the parent class of
 * AMatrixFunc and RMatrixFunc.
 */

class GenMatrixFunc {
 public:
  GenMatrixFunc() {};
  virtual ~GenMatrixFunc() {};
  /*!
   *This virtual function in implemented in the child class.
   */
  virtual void ClearMatrices() = 0;
  /*!
   *This virtual function in implemented in the child class.
   */
  virtual void FillMatrices(EPoint *) = 0;
  /*!
   *This virtual function in implemented in the child class.
   */
  virtual void InvertMatrices() = 0;
  /*!
   *This virtual function in implemented in the child class.
   */
  virtual void CalculateTMatrix(EPoint *) = 0;
  /// Cross section at a point from the T-matrix the subclass built.
  void CalculateCrossSection(EPoint *);

  /*!
   * Spin-summed |M|^2 from the channel-spin amplitude matrix, and the vector
   * analyzing power, for one point.
   *
   * Built from the same T-matrix elements CalculateCrossSection uses, by the
   * independent route of Seyler Eq. (4) rather than the Blatt-Biedenharn
   * contraction. Returns false when the point has no pathway, or when the
   * entrance is not a particle channel.
   *
   * The spin sum is returned without AZURE2's geometrical and unit factors, so
   * it is proportional to the differential cross section rather than equal to
   * it. That is enough for the check that matters: at fixed energy the ratio to
   * the existing cross section has to be constant in angle, which is what pins
   * down the coupling order and the choice of l versus l'.
   */
  bool CalculateAmplitudeMatrix(EPoint *point, double *spinSum, double *analyzingPower);
  /*!
   * Vector analyzing power of a capture reaction, from the Legendre
   * coefficients of R. G. Seyler and H. R. Weller, Phys. Rev. C \b 20 (1979)
   * 453, Eqs. (20) and (21):
   * \f$A_y = \sum_k b_k P_k^1(\cos\theta) / \sum_k a_k P_k(\cos\theta)\f$.
   * The coefficient table is built on first use by
   * CNuc::CalcCaptureAnalyzingPower; the T-matrix elements are the ones this
   * object has already computed, external capture included.
   * Returns false for anything that is not a photon exit channel.
   */
  bool CalculateCaptureAnalyzingPower(EPoint *point, double *unpolarized,
                                      double *analyzingPower);
  /// Start a temporary T-matrix element for a \f$J,l,l'\f$ combination.
  void NewTempTMatrix(TempTMatrix);
  /// As NewTempTMatrix, for the E1 component kept separately.
  void NewTempTMatrixE1(TempTMatrix);
  /// As NewTempTMatrix, for the E2 component.
  void NewTempTMatrixE2(TempTMatrix);
  /// Accumulate into temporary element \p i, 1-based.
  void AddToTempTMatrix(int, complex);
  /// Accumulate into an E1 temporary element.
  void AddToTempTMatrixE1(int, complex);
  /// Accumulate into an E2 temporary element.
  void AddToTempTMatrixE2(int, complex);
  /// Drop the temporary elements.
  void ClearTempTMatrices();
  /// Drop the E1 temporaries.
  void ClearTempTMatricesE1();
  /// Drop the E2 temporaries.
  void ClearTempTMatricesE2();
  /// Store an internal T-matrix element for a resonant pathway.
  void AddTMatrixElement(int, int, complex, int decayNum = 1);
  /// Store an external T-matrix element for an external-capture pathway.
  void AddECTMatrixElement(int, int, complex);
  /// 1-based position of the \f$J,l,l'\f$ temporary element, or 0 if absent.
  int IsTempTMatrix(double, int, int);
  /// 1-based position of the E1 temporary, or 0.
  int IsTempTMatrixE1(double, int, int);
  /// 1-based position of the E2 temporary, or 0.
  int IsTempTMatrixE2(double, int, int);
  /// Number of temporary elements.
  int NumTempTMatrices() const;
  /// Number of E1 temporaries.
  int NumTempTMatricesE1() const;
  /// Number of E2 temporaries.
  int NumTempTMatricesE2() const;
  TempTMatrix *GetTempTMatrix(int);
  TempTMatrix *GetTempTMatrixE1(int);
  TempTMatrix *GetTempTMatrixE2(int);
  /// Internal T-matrix element for a resonant pathway.
  complex GetTMatrixElement(int, int, int decayNum = 1) const;
  /// External T-matrix element for an external-capture pathway.
  complex GetECTMatrixElement(int, int) const;
  /// Racah coefficient \f$R_k\f$ of the angular-distribution algebra.
  double GetRk(double, double, double, double, int);

  /*!
   *This virtual function in implemented in the child class.
   */
  virtual CNuc *compound() const = 0;
  /*!
   *This virtual function in implemented in the child class.
   */
  virtual const Config &configure() const = 0;

 protected:
  /// Vector of internal T-matrix elements accessable to child class.
  std::vector<matrix_c> tmatrix_;
  /// Vector of external T-matrix elements accessable to child class.
  matrix_c ec_tmatrix_;

 private:
  std::vector<TempTMatrix> temp_t_matrices_;
  std::vector<TempTMatrix> temp_t_matrices_E1_;
  std::vector<TempTMatrix> temp_t_matrices_E2_;
};

#endif
