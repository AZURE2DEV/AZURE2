#ifndef AZUREAPI_H
#define AZUREAPI_H

#include "AZUREMain.h"

#include "Constants.h"
#include <vector>

class Config;
class EData;
class CNuc;
class EPoint;
struct GradAccum;

///A function class to perform the calculation of the chi-squared value

/*!
 * The AZUREAPI function class calculates the cross section based on a 
 * parameter set for all available data, and returns a chi-squared value.
 * This function class is what Minuit calls repeatedly during the fitting
 * process to perform the minimization.
 */

class AZUREAPI {
 public:
  /*!
   * The AZUREAPI object is created with reference to an EData and CNuc object.
   *. The runtime configurations are also passed through a Config structure.
   */
  AZUREAPI(Config& configure) : configure_(configure) { };
  
  ~AZUREAPI() {
  };

  bool Initialize( );

  // Update data objects, returns number of segments
  int UpdateData( );
  // Update segments values
  int UpdateSegments(vector_r& p);
  // Update segments values for RWA
  int UpdateSegmentsRWA(vector_r& p);
  // Update segments values for all RWA parameters
  int UpdateSegmentsAllRWA(vector_r& p);
  // Calculate the external capture for data
  bool CalculateExternalCapture( );
  // Reads the parameters values
  bool UpdateParameters( );
  // Reads the norms values
  void UpdateNorms( );
  // Set AZURE2 to calculate data points
  void SetData( );
  // Set AZURE2 to calculate extrapolations
  void SetExtrap( );
  // Set radius to a fixed value
  void SetRadius( int idx, double r );
  // Get indeces of normalization parameters
  vector_r GetNormalizationIndices( );
  // Get indeces of energy shift parameters
  vector_r GetEnergyShiftIndices( );
  // Number of numeric fields packed per parameter by GetParameterInfo().
  // Keep this in sync with pyazr/parameters.py (Parameter._NFIELDS).
  static const int kParamInfoFields = 15;
  /*!
   * Returns structured metadata describing every parameter, in the same order
   * as params_names() / params_all() / params_fixed().
   *
   * The result is a flat vector of kParamInfoFields doubles per parameter; the
   * fields (and -1 for "not applicable") are:
   *   0  type          0=energy, 1=width, 2=norm, 3=energy-shift
   *   1  jgroup        1-based J-group index (R-matrix params only)
   *   2  J             total spin of the J-group / level
   *   3  parity        +1 / -1
   *   4  level         1-based level index within the J-group
   *   5  level_energy  level energy (MeV)
   *   6  channel       1-based channel index (width params only)
   *   7  L             channel orbital angular momentum
   *   8  S             channel spin
   *   9  pair          channel particle-pair number
   *   10 radtype       ASCII code of the channel radiation type ('P','E','M')
   *   11 fixed         1 if the parameter is fixed, 0 if free
   *   12 value         current (physical) parameter value
   *   13 segment_key   data-segment key (norm / energy-shift params only)
   *   14 wigner_limit  Wigner limit of the channel reduced width, i.e. the
   *                    bound on |reduced-width amplitude| (width params only)
   */
  vector_r GetParameterInfo( ) const;
  
  /*!
   * Returns a reference to the Config structure.
   */
  Config &configure() const {return configure_;};
  /*!
   * Returns a pointer to the EData object.
   */
  EData *data() const {return data_;};
  /*!
   * Returns a pointer to the CNuc object.
   */
  CNuc *compound() const {return compound_;};
  /*!
   * Returns a pointer to the parameter values object.
   */
  vector_r params_values() const {return values_;};
  /*!
   * Returns a pointer to the parameter values object for RWA.
   */
  vector_r params_values_rwa() const {return values_rwa_;};
  /*!
   * Returns a pointer to the parameter names object.
   */
  std::string params_names(int i) const {return names_[i];};
  /*!
   * Returns a pointer to the parameter names object.
   */
  vector_r params_all() const {return all_;};
  /*!
   * Returns a pointer to the parameter names object for RWA.
   */
  vector_r params_all_rwa() const {return all_rwa_;};
  /*!
   * Returns a pointer to the fixed parameters object.
   */
  std::vector<bool> params_fixed() const {return fixed_;};
  /*!
   * Returns a pointer to the calculated segments object.
   */
  vector_r data_energies(int i) const {return dataEnergies_[i];};
  /*!
   * Returns a pointer to the calculated segments object.
   */
  vector_r data_angles(int i) const {return dataAngles_[i];};
  /*!
   * Returns a pointer to the calculated segments object.
   */
  vector_r data_segments(int i) const {return dataSegments_[i];};
  /*!
   * Returns a pointer to the calculated segments object.
   */
  vector_r data_segments_errors(int i) const {return dataSegmentsErrors_[i];};
  /*!
   * Returns a pointer to the calculated segments object.
   */
  vector_r calculated_segments(int i) const {return calculatedSegments_[i];};
  /*!
   * Returns a pointer to the calculated E1 segments object.
   */
  vector_r calculated_segments_e1(int i) const {return calculatedSegmentsE1_[i];};
  /*!
   * Returns a pointer to the calculated E2 segments object.
   */
  vector_r calculated_segments_e2(int i) const {return calculatedSegmentsE2_[i];};
  /*!
   * Returns a pointer to the calculated energies object.
   */
  vector_r calculated_energies(int i) const {return calculatedEnergies_[i];};
  /*!
   * Returns a pointer to the calculated energies object.
   */
  vector_r calculated_angles(int i) const {return calculatedAngles_[i];};
  /*!
   * Returns a pointer to the data excitation energies.
   */
  vector_r data_excitation_energies(int i) const {return dataExcitationEnergies_[i];};
  /*!
   * Returns a pointer to the calculated excitation energies.
   */
  vector_r calculated_excitation_energies(int i) const {return calculatedExcitationEnergies_[i];};
  /*!
   * Returns a pointer to the segments norms.
   */
  vector_r norms( ) const {return norms_;};
  /*!
   * Returns a pointer to the segments norms errors.
   */
  vector_r norms_errors( ) const {return normsErrors_;};
  /*!
   * Returns a pointer to the data sfactor conversion.
   */
  vector_r data_conv( int i ) const {return dataConv_[i];};
  /*!
   * Returns a pointer to the calculated sfactor conversion.
   */
  vector_r calculated_conv( int i ) const {return calculatedConv_[i];};
  /*!
   * Transform RWA parameters to physical values
   */
  vector_r TransformRWAParameters(const vector_r& p) const;
  /*!
   * Transform all RWA parameters to physical values
   */
  vector_r TransformAllRWAParameters(const vector_r& p) const;

  /*!
   * Calculate chi-squared from RWA parameters
   */
  double CalculateChi2RWA(const vector_r& rwaParams) const;
  
  /*!
   * Calculate chi-squared from physical parameters
   */
  double CalculateChi2Physical(const vector_r& physicalParams) const;

  /*!
   * Value and analytic gradient of the (data) chi-squared with respect to the
   * non-fixed RWA parameters.  Input: the non-fixed RWA parameters (energies,
   * reduced widths, normalizations), as for CalculateChi2RWA.  Returns
   *   [ chi2, d(chi2)/dp_0, ..., d(chi2)/dp_{n-1} ].
   * Energies / reduced widths / normalizations are analytic; energy shifts are
   * finite-differenced.  (For a log-likelihood use lnL = -0.5*(chi2 + const),
   * grad lnL = -0.5*grad chi2.)
   */
  vector_r CalculateChi2GradRWA(const vector_r& params) const;

  /*!
   * Standardized residuals r_i = (fit_i - data_i*n)/(cmErr_i*n) (so sum r_i^2 =
   * chi2) and their analytic Jacobian J_{ij} = d r_i / d theta_j, for
   * Gauss-Newton / Levenberg-Marquardt.  Columns are the non-fixed parameters
   * (input order).  Returns
   *   [ nRes, nCols, r_0..r_{nRes-1}, J row-major (nRes x nCols) ],
   * or [ -1 ] if a point is outside the supported analytic path.  Energy-shift
   * columns are left zero.
   */
  vector_r CalculateResidualJacobianRWA(const vector_r& params) const;


 private:

  /*!
   * Analytic reverse-mode gradient of the data chi-squared w.r.t. the energy and
   * reduced-width (gamma) parameters, accumulated into the energy/gamma entries
   * of gradFull; the data-term normalization gradient is accumulated into the
   * norm entries.  Returns false (touching nothing) if any data point is outside
   * the supported analytic path, so the caller falls back to finite differences.
   */
  bool Chi2GradEGammaNorm(const vector_r& fullParams, vector_r& gradFull) const;

  // Configuration
  Config &configure_;
  EData *data_;
  CNuc *compound_;

  // Parameters
  std::vector<bool> fixed_;
  std::vector<std::string> names_;
  vector_r all_, all_rwa_, values_, values_rwa_, norms_, normsErrors_;

  // Data
  std::vector<vector_r> dataConv_;
  std::vector<vector_r> dataEnergies_;
  std::vector<vector_r> dataAngles_;
  std::vector<vector_r> dataSegments_;
  std::vector<vector_r> dataSegmentsErrors_;
  std::vector<vector_r> dataExcitationEnergies_;

  std::vector<vector_r> calculatedConv_;
  std::vector<vector_r> calculatedEnergies_;
  std::vector<vector_r> calculatedAngles_;
  std::vector<vector_r> calculatedSegments_;
  std::vector<vector_r> calculatedSegmentsE1_;
  std::vector<vector_r> calculatedSegmentsE2_;
  std::vector<vector_r> calculatedExcitationEnergies_;

};

#endif