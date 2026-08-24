#ifndef EPOINT_H
#define EPOINT_H

#include <cstring>
#include "Constants.h"

/// A container structure for a reference to a data point.

/*!
 * If a point is mapped back to another in the calculation, this structure
 * hold the relevant indices of the map point.
 */

struct EnergyMap {
  /// The segment index for the map
  int segment;
  /// The point index for the map
  int point;
};

class ESegment;
class EData;
class CNuc;
class PPair;
class TargetEffect;
class DataLine;
class Config;

/// An AZURE data point

/*!
 * A data point object in AZURE consists of a defined entrance and exit pair, an energy, an angle,
 * measured cross section and uncertainty, s-factor conversions, and
 * several flags that determine the type of data (angle integrated or differential) to be analysed.
 */

class EPoint {
 public:
  /// Build from a line of a data file, taking its conventions from the parent segment.
  EPoint(DataLine, ESegment *);
  /// Build at a given energy and angle, for a segment generated from a grid.
  EPoint(double, double, ESegment *);
  /// Build with every convention given explicitly, rather than inherited from a segment.
  EPoint(double, double, int, int, bool, bool, bool, double, int, int);
  /// Differential cross section?
  bool IsDifferential() const;
  /// Phase shift?
  bool IsPhase() const;
  /// THM (modified R-matrix, half-off-shell) point?
  bool IsTHM() const;
  /// Is this point calculated by another? Points at equal energies are mapped onto one so the energy-dependent work is done once.
  bool IsMapped() const;
  /// Does the point carry target effects?
  bool IsTargetEffect() const;
  /// Angular distribution?
  bool IsAngularDist() const;
  //! Vector analyzing power point; the fit value is A_y, not a cross section.
  bool IsAnalyzingPower() const { return is_analyzing_power_; };
  void SetIsAnalyzingPower(bool v) { is_analyzing_power_ = v; };
  //! A_y is kept beside the cross section rather than replacing it, because
  //! target-effect integration needs the cross section as the weight.
  double GetAnalyzingPower() const { return analyzing_power_; };
  void SetAnalyzingPower(double v) { analyzing_power_ = v; };
  /// Is this one of the sub-points a target-effect integral is built from?
  bool IsSubPoint() const { return is_sub_point_; };
  /// Unobserved-primary, observed-secondary point?
  bool IsUPOS() const;
  /// Entrance pair key, as written in the input files; not necessarily its PPair position.
  int GetEntranceKey() const;
  /// Exit pair key, as written in the input files.
  int GetExitKey() const;
  /// Highest Legendre order stored for this point.
  int GetMaxLOrder() const;
  /// Orbital angular momentum. Phase-shift points only.
  int GetL() const;
  /// How many points are mapped onto this one.
  int NumLocalMappedPoints() const;
  /// Number of sub-points, the samples a target-effect integral is evaluated on.
  int NumSubPoints() const;
  /// 1-based position of the TargetEffect in the parent EData.
  int GetTargetEffectNum() const;
  /// Highest Legendre order. Angular-distribution points only.
  int GetMaxAngDistOrder() const;
  /// Number of stored angular-distribution coefficients.
  int GetNumAngularDists() const;
  /// Angular momentum of the secondary decay. UPOS points only.
  int GetSecondaryDecayL() const;
  /// Angle in the lab frame.
  double GetLabAngle() const;
  /// Angle in the centre-of-mass frame.
  double GetCMAngle() const;
  /// Energy in the lab frame.
  double GetLabEnergy() const;
  /// Energy in the centre-of-mass frame.
  double GetCMEnergy() const;
  /// Excitation energy of the compound nucleus -- the axis shared by every entrance pair.
  double GetExcitationEnergy() const;
  /// Energy as read, before any energy shift was applied.
  double GetOriginalEnergy() const;
  /// Legendre polynomial of the given order at this point's angle.
  double GetLegendreP(int) const;
  /// Measured cross section, lab frame.
  double GetLabCrossSection() const;
  /// Measured cross section, centre-of-mass frame.
  double GetCMCrossSection() const;
  /// Measured uncertainty, lab frame.
  double GetLabCrossSectionError() const;
  /// Measured uncertainty, centre-of-mass frame.
  double GetCMCrossSectionError() const;
  /// Geometrical factor \f$\pi/k^2\f$.
  double GetGeometricalFactor() const;
  /// Cross section AZURE2 calculated here.
  double GetFitCrossSection() const;
  /// Calculated cross section, E1 component only.
  double GetFitE1CrossSection() const;
  /// Calculated cross section, E2 component only.
  double GetFitE2CrossSection() const;
  /// Multiply a cross section by this to get the astrophysical S-factor.
  double GetSFactorConversion() const;
  /// \f$\sqrt{P_c}\f$ for the channel at (J-group, channel), both 1-based.
  double GetSqrtPenetrability(int, int) const;
  /// THM entrance transfer form factor M_l = (b-1) j_l - rho dj_l/drho for
  /// JGroup/channel, assembled with the given boundary b (0 if not stored).
  double GetThmFormFactor(int, int, double) const;
  /// Total spin. Phase-shift points only.
  double GetJ() const;
  /// Stopping cross section at this sub-point, for a yield-curve target integration.
  double GetStoppingPower() const;
  /// Beam energy loss across the target at this point.
  double GetTargetThickness() const;
  /// Angular-distribution coefficient of the given order.
  double GetAngularDist(int) const;
  /// Jacobian converting the angle from lab to centre of mass.
  double GetAngleKinFactor() const;
  /// Jacobian converting the cross section from lab to centre of mass.
  double GetCrossSectionKinFactor() const;
  /// Final-state spin. UPOS points only.
  double GetIc() const;
  /// Multipole mixing ratio. UPOS points only.
  double GetDelta() const;
  /// Diagonal \f$L_o\f$ element for the channel at (J-group, channel), both 1-based.
  complex GetLoElement(int, int) const;
  /// \f$\exp(i\omega_c)\f$, the Coulomb phase, for the channel at (J-group, channel).
  complex GetExpCoulombPhase(int, int) const;
  /// \f$\exp(i\delta_c)\f$, the hard-sphere phase, for the channel at (J-group, channel).
  complex GetExpHardSpherePhase(int, int) const;
  /// Coulomb (Rutherford) amplitude \f$C_\alpha\f$.
  complex GetCoulombAmplitude() const;
  /// External-capture amplitude for the pathway at (KGroup, ECMGroup), both 1-based.
  complex GetECAmplitude(int, int) const;
  /// As GetECAmplitude, but interpolated to the shifted energy through the amplitude cache.
  complex GetECAmplitudeWithShift(int, int, CNuc *, const Config &) const;
  /// Where this point is mapped, if it is.
  EnergyMap GetMap() const;
  /// Compute everything that depends on energy but not on the fit parameters. Run once before fitting.
  void Initialize(CNuc *, const Config &);
  /// Lab energy to centre-of-mass energy.
  double ConvertLabValue(double, PPair *);
  /// Centre-of-mass energy to lab energy.
  double ConvertCMValue(double, PPair *);
  /// Fill the centre-of-mass energy from the lab energy.
  void ConvertLabEnergy(PPair *);
  /// Fill the compound excitation energy.
  void ConvertExcitationEnergy(PPair *);
  /// Total decay energy from the light particle's, with the parent at rest.
  void ConvertDecayEnergy(PPair *);
  /// Fill the centre-of-mass angle from the lab angle.
  void ConvertLabAngle(PPair *);
  /// As above for a reaction with different entrance and exit pairs.
  void ConvertLabAngle(PPair *, PPair *, const Config &);
  /// Centre-of-mass angle by the relativistic form, Iliadis C.37-C.38.
  void ConvertCMAngle(PPair *, PPair *, const Config &);
  /// Fill the centre-of-mass cross section and uncertainty from the lab ones.
  void ConvertCrossSection(PPair *, PPair *);
  /// Centre-of-mass gamma angle, with the relativistic correction.
  void ConvertLabAngleGammas(PPair *);
  /// Centre-of-mass gamma cross section.
  void ConvertCrossSectionGammas(PPair *);
  double CalculateCrossSectionConversionFactor(PPair *, PPair *);
  double CalculateCrossSectionGammaConversionFactor(PPair *);
  /// Append the next Legendre polynomial; they are stored in order L = 0, 1, 2, ...
  void AddLegendreP(double);
  /// Drop the stored polynomials, so they can be recomputed after an energy shift.
  void ClearLegendrePolynomials();
  /// Set \f$\pi/k^2\f$.
  void SetGeometricalFactor(double);
  /// Record the calculated cross section.
  void SetFitCrossSection(double);
  void SetFitE1CrossSection(double);
  void SetFitE2CrossSection(double);
  /// Set the cross-section to S-factor conversion.
  void SetSFactorConversion(double);
  void SetLabEnergy(double);
  void SetCMEnergy(double);
  void SetExcitationEnergy(double);
  void SetLabAngle(double);
  void SetCMAngle(double);
  void SetExitKey(int);
  void SetEntranceKey(int);
  /// Compute the Legendre polynomials up to a maximum order.
  void CalcLegendreP(int, CNuc *, TargetEffect *);
  /// Compute the geometrical factor, S-factor conversion, \f$L_o\f$ elements, penetrabilities and phases together.
  void CalcEDependentValues(CNuc *, const Config &);
  /// Recompute those at the current, possibly shifted, energy.
  void RecalcEDependentValues(CNuc *, const Config &);
  /// Store an \f$L_o\f$ element at (J-group, channel).
  void AddLoElement(int, int, complex);
  /// Store a \f$\sqrt{P_c}\f$ at (J-group, channel).
  void AddSqrtPenetrability(int, int, double);
  /// Store the pieces of the THM form factor at (J-group, channel).
  void AddThmFormFactor(int, int, double, double);
  /// Store a Coulomb phase factor at (J-group, channel).
  void AddExpCoulombPhase(int, int, complex);
  /// Store a hard-sphere phase factor at (J-group, channel).
  void AddExpHardSpherePhase(int, int, complex);
  /// Compute the Coulomb amplitude \f$C_\alpha\f$.
  void CalcCoulombAmplitude(CNuc *);
  /// Set the Coulomb amplitude directly.
  void SetCoulombAmplitude(complex);
  /// Compute the external-capture amplitudes for every matching pathway.
  void CalculateECAmplitudes(CNuc *, const Config &);
  void AddECAmplitude(int, int, complex);
  void AddECAmplitude(int, int, complex, double);
  void ClearECAmplitudes();
  void Calculate(CNuc *, const Config &configure, EPoint *parent = NULL, int subPointNum = 0);
  void SetMap(int, int);
  void ClearMapping();
  void AddLocalMappedPoint(EPoint *);
  void ClearLocalMappedPoints();
  void SetTargetEffectNum(int);
  void AddSubPoint(EPoint);
  void IntegrateTargetEffect(const Config &);
  void IntegrateTargetEffectForObservable(const Config &);
  void SetParentData(EData *);
  void SetStoppingPower(double);
  void SetTargetThickness(double);
  void SetAngularDists(vector_r);
  void SetAngleKinFactor(double);
  void SetCrossSectionKinFactor(double);
  EData *GetParentData() const;
  EPoint *GetLocalMappedPoint(int) const;
  EPoint *GetSubPoint(int);
  std::vector<EPoint> &GetSubPoints();
  std::vector<EPoint *> &GetMappedPoints();
  void StoreSubpointOffsets();                                                          // Store offsets for adaptive grid preservation
  void ApplySubpointShift(double energyShift, CNuc *theCNuc, const Config &configure);  // Intelligent shift preserving resonance structure
 private:
  bool is_differential_;
  bool is_phase_;
  bool is_thm_;
  bool is_mapped_;
  bool is_ang_dist_;
  bool is_analyzing_power_ = false;
  bool is_sub_point_ = false;
  double analyzing_power_ = 0.0;
  int entrance_key_;
  int exit_key_;
  int segment_key_;
  int l_value_;
  int targetEffectNum_;
  int max_ang_dist_order_;
  double cm_angle_;
  double lab_angle_;
  double original_energy_;
  double cm_energy_;
  double lab_energy_;
  double excitation_energy_;
  double cm_crosssection_;
  double cm_dcrosssection_;
  double lab_crosssection_;
  double lab_dcrosssection_;
  double geofactor_;
  double fitcrosssection_;
  double fitE1crosssection_;
  double fitE2crosssection_;
  double sfactorconv_;
  double j_value_;
  double stoppingPower_;
  double targetThickness_;
  double angleKinFactor_;
  double crossSectionKinFactor_;
  bool isUPOS_;
  int secondaryDecayL_;
  double Ic_;
  double delta_;
  struct EnergyMap energy_map_;
  complex coulombamplitude_;
  vector_r legendreP_;
  vector_r angularDists_;
  matrix_c lo_elements_;
  matrix_r penetrabilities_;
  /// Boundary-independent THM form-factor pieces j_l(rho) and rho dj_l/drho,
  /// indexed [jGroup-1][channel-1]. M_l is assembled at the entrance vertex
  /// with the per-level boundary (Brune) or the channel boundary constant.
  matrix_r thm_jl_;
  matrix_r thm_rhodjl_;
  matrix_c coulombphase_;
  matrix_c hardspherephase_;
  matrix_c ec_amplitudes_;
  matrix_r ec_energies_;  // Energies at which EC amplitudes were calculated
  std::vector<EPoint *> local_mapped_points_;
  std::vector<EPoint> integrationPoints_;
  EData *parentData_;
  ESegment *parentSegment_;
};

#endif
