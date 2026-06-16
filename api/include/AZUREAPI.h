#ifndef AZUREAPI_H
#define AZUREAPI_H

#include "AZUREMain.h"

#include "Constants.h"
#include <vector>

class Config;
class EData;
class CNuc;

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
   * Calculate the Gaussian log-likelihood from RWA parameters with per-segment
   * error inflation.
   *
   * The input vector is the concatenation of the non-fixed RWA parameters
   * (which include the normalizations) followed by one error-inflation factor
   * per segment, in the same segment order as norms() / UpdateData(). For each
   * data point the variance is inflated as
   *   var = (dataErr * norm)^2 + (f * model)^2
   * where model is the theoretical (fit) cross section and f is the segment's
   * inflation factor. The returned value is
   *   lnL = -0.5 * sum_i [ (fit - data*norm)^2 / var_i + ln(2*pi*var_i) ]
   * i.e. it includes the error-normalization term so the inflation factors are
   * self-regulating.
   */
  double CalculateLnLRWA(const vector_r& params) const;

  /*!
   * Calculate the Gaussian log-likelihood from RWA parameters using a full
   * per-segment covariance matrix with correlated error inflation.
   *
   * The input vector is packed identically to CalculateLnLRWA: the non-fixed
   * RWA parameters (including normalizations) followed by one error-inflation
   * factor per segment. Within each segment the data share a single
   * normalization, so the inflation is treated as 100% correlated. The
   * covariance block for a segment is
   *   C_ij = (dataErr_i * norm)^2 * delta_ij + f^2 * model_i * model_j
   * i.e. the diagonal carries statistical-plus-inflation variance while the
   * off-diagonal carries only the (fully-correlated) inflation term. Different
   * segments are uncorrelated (block-diagonal C). The returned value is the
   * multivariate Gaussian log-likelihood
   *   lnL = -0.5 * [ r^T C^{-1} r + ln det(2*pi*C) ]
   * with r_i = fit_i - data_i*norm, summed over the segment blocks.
   */
  double CalculateLnLCovRWA(const vector_r& params) const;


 private:

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