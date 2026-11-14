#ifndef GENMATRIXFUNC_H
#define GENMATRIXFUNC_H

#include "Constants.h"
#include "Config.h"
#include "Site.h"
#include "LevelState.h"

#include <cstring>
#include <string>
#include <vector>

class CompoundSystem;
namespace azure2 {
struct Results;
struct EntrancePathways;
struct ExitRow;
}

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
  virtual void FillMatrices(Site *) = 0;
  /*!
   *This virtual function in implemented in the child class.
   */
  virtual void InvertMatrices() = 0;
  /*!
   *This virtual function in implemented in the child class.
   */
  virtual void CalculateTMatrix(Site *) = 0;
  /// The T-matrix of another exit pair at the same energy: drops the elements held now (and what was derived from
  /// them) and builds those of `point`'s exit from the solved level matrix, which does not depend on the exit.  Used
  /// for the T-matrix groups of an energy task after the first (SolveAtOneEnergy).
  void RecalculateTMatrix(Site *point) {
    ResetTMatrix();
    CalculateTMatrix(point);
  }
  /// True when `other` holds the same T-matrix bit for bit, internal and external elements (the shared-solve check
  /// of AZURE_PLAN_CHECK, CalculationPlan.h).
  bool SameTMatrix(const GenMatrixFunc &other) const {
    return SameBits(tmatrix_, other.tmatrix_) && SameBits(ec_tmatrix_, other.ec_tmatrix_);
  }
  /// The observable of the point from the T-matrix the subclass built, times the electron-screening factor where the
  /// point carries one (CalculateBareCrossSection, then Unit::ScreeningFactor).
  void CalculateCrossSection(Site *, azure2::Results &results);
  /// The observable without electron screening.
  void CalculateBareCrossSection(Site *, azure2::Results &results);

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
  /// The refusal of a point whose exit pair no pathway of its entrance pair reaches (pair numbers, 1-based).
  std::string NoPathwayMessage(int entrancePair, int exitPair) const;
  /// The refusal of a capture-A_y point whose terms were not requested at set-up (CompoundSetUpRequests::AddCaptureAy).
  std::string CaptureAyNotBuiltMessage(const Site *point) const;

  /*!
   *This virtual function in implemented in the child class.
   */
  virtual const CompoundSystem *compound() const = 0;
  /// The levels at the evaluated parameter vector (LevelState.h), held by the evaluation that made this object.
  virtual const LevelState &state() const = 0;
  /*!
   *This virtual function in implemented in the child class.
   */
  virtual const Config &configure() const = 0;

 protected:
  /// Vector of internal T-matrix elements accessable to child class.
  std::vector<matrix_c> tmatrix_;
  /// Vector of external T-matrix elements accessable to child class.
  matrix_c ec_tmatrix_;
  static bool SameBits(const vector_c &a, const vector_c &b) {
    return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(complex)) == 0);
  }
  template <class V>
  static bool SameBits(const std::vector<V> &a, const std::vector<V> &b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); i++)
      if (!SameBits(a[i], b[i])) return false;
    return true;
  }
  /// Forget everything derived from the T-matrix (ClearMatrices).
  void InvalidateTMatrixProducts() { rt_decay_ = nullptr; }
  /// Drop the T-matrix elements and their products, keeping the solved matrices (RecalculateTMatrix).
  virtual void ResetTMatrix() {
    tmatrix_.clear();
    ec_tmatrix_.clear();
    InvalidateTMatrixProducts();
  }

 private:
  /// The interference products Z1Z2 T1 conj(T2) of the differential sum,
  /// summed by Legendre order (rt_wl_[L]), for the T-matrix held now: they do
  /// not depend on the angle, so the points that share one T-matrix
  /// (the units of one T-matrix group of an energy task) compute them once, and each point's sum is one term per
  /// order.  Entrance = exit pair, capture, and particle exits other than
  /// UPOS (whose sum keeps only the k-groups of equal channel spins).
  const azure2::ExitRow *rt_decay_ = nullptr;
  std::vector<complex> rt_wl_;
  /// The Coulomb-nuclear interference per (channel spin, l), BuildCoulombInterference.
  struct CoulombTerm {
    double s;
    int l;
    complex c;
  };
  std::vector<CoulombTerm> rt_it_;
  void BuildCoulombInterference(const azure2::EntrancePathways &P, const azure2::ExitRow &exit, bool elastic);
  /// rt_wl_ (and for elastic scattering rt_it_) for the T-matrix held now, once per T-matrix.
  void BuildOrderProducts(const azure2::EntrancePathways &P, const azure2::ExitRow &exit, bool particleExit, bool elastic);
  /// The record of an analysing power at the T-matrix held now (SiteObservable.h, ay_record): the cross section's
  /// coefficients and A_y's -- the amplitude matrix's per exit orbital momentum (Polarization::AmplitudeMatrix::
  /// AddPathwayCoefficients) for a particle exit, Seyler and Weller's per order for capture.
  Site::AnalyzingPowerCoefficients AnalyzingPowerCoefficientsOf(const Site *point, const azure2::EntrancePathways &P,
                                                                const azure2::ExitRow &theDecay, int aa,
                                                                int exitPairNum);
  /// The observables of elastic scattering of `point`'s pair `aa` at its energy (SiteObservable::kElastic).
  Site::ElasticCoefficients ElasticCoefficientsOf(const Site *point, const azure2::EntrancePathways &P,
                                                  const azure2::ExitRow &exit, int aa);
  /// The Legendre coefficients B_L (b/sr) of dsigma/dOmega_c.m. of `point`'s reaction at its energy (for elastic
  /// scattering the nuclear amplitude's part).
  std::vector<double> LegendreCoefficientsOf(const Site *point, const azure2::EntrancePathways &P,
                                             const azure2::ExitRow &exit, int aa, int exitPairNum);
  std::vector<TempTMatrix> temp_t_matrices_;
  std::vector<TempTMatrix> temp_t_matrices_E1_;
  std::vector<TempTMatrix> temp_t_matrices_E2_;
};

#endif
