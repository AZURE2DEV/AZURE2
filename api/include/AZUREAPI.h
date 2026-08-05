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
  // Set the channel radius of particle pair idx (1-based) to r fm and rebuild
  // the compound nucleus, data and parameters.  Returns false if the rebuild
  // failed, in which case the instance is no longer usable.
  bool SetRadius( int idx, double r );
  // Get indeces of normalization parameters
  vector_r GetNormalizationIndices( );
  // Get indeces of energy shift parameters
  vector_r GetEnergyShiftIndices( );
  // Number of numeric fields packed per parameter by GetParameterInfo().
  // Keep this in sync with pyazr/parameters.py (Parameter._NFIELDS).
  static const int kParamInfoFields = 16;
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
   *   15 input_is_rwa  1 if the width's .azr input value was declared as a
   *                    reduced width amplitude (MeV^(1/2)) rather than a
   *                    physical partial width/ANC, 0 otherwise; the parameter
   *                    transformations pass such widths through unconverted
   *                    (width params only, -1 otherwise)
   */
  vector_r GetParameterInfo( ) const;

  // Number of numeric fields packed per particle pair by GetPairsInfo().
  // Keep this in sync with pyazr/parameters.py (Pair._NFIELDS).
  static const int kPairInfoFields = 16;
  /*!
   * Returns structured metadata describing every particle pair, in 1-based
   * pair-number order (the same number stored in field 9 -- "pair" -- of
   * GetParameterInfo(), so a width parameter can be matched to its pair).
   *
   * The result is a flat vector of kPairInfoFields doubles per pair; the fields
   * are:
   *   0  pair          1-based pair number (matches Parameter "pair")
   *   1  pair_key      user pair key from the .azr file
   *   2  ptype         particle type (0 = particle channel, otherwise photon)
   *   3  is_entrance   1 if this is the entrance pair, else 0
   *   4  J1            intrinsic spin of particle 1
   *   5  parity1       parity of particle 1 (+1 / -1)
   *   6  Z1            charge number of particle 1
   *   7  M1            mass of particle 1 (amu)
   *   8  J2            intrinsic spin of particle 2
   *   9  parity2       parity of particle 2 (+1 / -1)
   *   10 Z2            charge number of particle 2
   *   11 M2            mass of particle 2 (amu)
   *   12 sepE          separation energy (MeV)
   *   13 exE           excitation energy of the pair (MeV)
   *   14 chRad         channel radius (fm)
   *   15 i1i2factor    1 / ((2 J1 + 1)(2 J2 + 1)); the entrance pair's value is
   *                    the denominator of the statistical spin factor
   *                    omega = (2 J + 1) * i1i2factor
   *   16 bindingE      Trojan Horse binding energy of the pair (MeV); 0 unless
   *                    the pair is used as a THM entrance channel
   */
  vector_r GetPairsInfo( ) const;
  
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
   * Legendre coefficients of the angular distribution for segment \p i, one
   * group per point and each group self-describing, so a segment whose points
   * carry different orders still round-trips:
   *
   *   [ n_0, c_0_0 ... c_0_(n_0-1), n_1, c_1_0 ... c_1_(n_1-1), ... ]
   *
   * Read a count, consume that many coefficients, repeat until the array is
   * exhausted; the number of groups is the number of points.
   *
   * Only points belonging to an angular-distribution segment carry any; every
   * other point contributes a count of zero.
   */
  vector_r calculated_angular_dists(int i) const {return calculatedAngularDists_[i];};
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

  /*!
   * Per-point sensitivities d(model)/d(theta) of the calculated segments, for
   * covariance uncertainty bands: sigma^2 = g^T C g (SAMMY Eq. IV E4.2).
   *
   * Columns are the free *R-matrix* parameters -- level energies and reduced
   * width amplitudes, in packed order -- which is exactly what
   * output/covariance.dat spans; normalizations and energy shifts are omitted
   * because no calculated observable depends on them.  Rows follow the segment
   * and point order of UpdateSegments / GET_CALCULATED_SEGMENT, so row k of
   * segment s belongs to calculated point k of segment s.
   *
   * Input: the non-fixed RWA parameters, as for UpdateSegmentsRWA.  Returns
   *   [ nSegments, nCols, nPoints_0 .. nPoints_{nSeg-1}, G row-major ],
   * with G holding sum(nPoints) rows of nCols, or [ -1 ] if a point is outside
   * the supported analytic path.
   *
   * One reverse-mode adjoint per point gives that point's whole row, so this
   * costs about two forward evaluations regardless of the parameter count --
   * against the 2*nCols forward passes a finite-difference band would need.
   */
  vector_r CalculateModelGradientsRWA(const vector_r& params) const;

  /*!
   * Coulomb wave functions on a requested energy grid.
   *
   * Request: [pairKey, l, radius, nE, E_1 ... E_nE], energies in MeV (centre of
   * mass), radius in fm.  A radius of zero means "use the pair's own channel
   * radius", which is where a penetrability or a hard-sphere phase is wanted.
   *
   * Response: [nE, then per energy: F, dF, G, dG, P, S, deltaHS], with P the
   * penetrability, S the shift function and deltaHS the hard-sphere phase shift
   * in radians.  Whether the values come from the accurate Coulomb routine,
   * from GSL, or from Numerov integration through a nuclear potential follows
   * the run's own configuration -- so this is also how one sees what the hybrid
   * model does to the external region.
   */
  vector_r GetCoulombFunctions(const vector_r& request) const;

  /*!
   * External-capture integrals on a requested energy grid.
   *
   * Request: [pairKey, nE, E_1 ... E_nE].  Every external-capture pathway the
   * compound nucleus generates from that entrance pair is evaluated at every
   * energy.
   *
   * Response: [nPathways, nE, then per pathway six descriptors
   * (li, lf, 2*si, 2*sf, multipolarity, radiationType) followed by 2*nE numbers
   * (real, imaginary part of the integral at each energy).
   *
   * These are the integrals the capture cross section is built from, and they
   * are the most expensive thing in a capture calculation --- which is why they
   * are cached.  Exposing them makes both facts checkable from a script.
   */
  vector_r GetECIntegrals(const vector_r& request) const;

  /*!
   * Coulomb-function cache counters, aggregated over threads.
   *
   * Response: [queries, hits, entries, keys, disabledKeys, threads].
   */
  vector_r GetCacheStats( ) const;


 private:

  /*!
   * Analytic reverse-mode gradient of the data chi-squared w.r.t. the energy and
   * reduced-width (gamma) parameters, accumulated into the energy/gamma entries
   * of gradFull; the data-term normalization gradient is accumulated into the
   * norm entries.  Returns false (touching nothing) if any data point is outside
   * the supported analytic path, so the caller falls back to finite differences.
   *
   * On success also returns, via chi2Out, the data chi-squared -- a free
   * byproduct of the forward model the adjoint already evaluates at every point,
   * consistent with the gradient (same per-point residual). Left untouched when
   * the analytic path bails.
   */
  bool Chi2GradEGammaNorm(const vector_r& fullParams, vector_r& gradFull,
                          double& chi2Out) const;

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
  std::vector<vector_r> calculatedAngularDists_;
  std::vector<vector_r> calculatedSegments_;
  std::vector<vector_r> calculatedSegmentsE1_;
  std::vector<vector_r> calculatedSegmentsE2_;
  std::vector<vector_r> calculatedExcitationEnergies_;

};

#endif