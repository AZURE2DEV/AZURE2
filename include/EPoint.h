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
  EPoint(DataLine, ESegment *);
  EPoint(double, double, ESegment *);
  EPoint(double, double, int, int, bool, bool, bool, double, int, int);
  bool IsDifferential() const;
  bool IsPhase() const;
  bool IsTHM() const;
  bool IsMapped() const;
  bool IsTargetEffect() const;
  bool IsAngularDist() const;
  //! Vector analyzing power point; the fit value is A_y, not a cross section.
  bool IsAnalyzingPower() const { return is_analyzing_power_; };
  void SetIsAnalyzingPower(bool v) { is_analyzing_power_ = v; };
  //! A_y is kept beside the cross section rather than replacing it, because
  //! target-effect integration needs the cross section as the weight.
  double GetAnalyzingPower() const { return analyzing_power_; };
  void SetAnalyzingPower(double v) { analyzing_power_ = v; };
  bool IsSubPoint() const { return is_sub_point_; };
  bool IsUPOS() const;
  int GetEntranceKey() const;
  int GetExitKey() const;
  int GetMaxLOrder() const;
  int GetL() const;
  int NumLocalMappedPoints() const;
  int NumSubPoints() const;
  int GetTargetEffectNum() const;
  int GetMaxAngDistOrder() const;
  int GetNumAngularDists() const;
  int GetSecondaryDecayL() const;
  double GetLabAngle() const;
  double GetCMAngle() const;
  double GetLabEnergy() const;
  double GetCMEnergy() const;
  double GetExcitationEnergy() const;
  double GetOriginalEnergy() const;
  double GetLegendreP(int) const;
  double GetLabCrossSection() const;
  double GetCMCrossSection() const;
  double GetLabCrossSectionError() const;
  double GetCMCrossSectionError() const;
  double GetGeometricalFactor() const;
  double GetFitCrossSection() const;
  double GetFitE1CrossSection() const;
  double GetFitE2CrossSection() const;
  double GetSFactorConversion() const;
  double GetSqrtPenetrability(int, int) const;
  /// THM entrance transfer form factor M_l = (b-1) j_l - rho dj_l/drho for
  /// JGroup/channel, assembled with the given boundary b (0 if not stored).
  double GetThmFormFactor(int, int, double) const;
  double GetJ() const;
  double GetStoppingPower() const;
  double GetTargetThickness() const;
  double GetAngularDist(int) const;
  double GetAngleKinFactor() const;
  double GetCrossSectionKinFactor() const;
  double GetIc() const;
  double GetDelta() const;
  complex GetLoElement(int, int) const;
  complex GetExpCoulombPhase(int, int) const;
  complex GetExpHardSpherePhase(int, int) const;
  complex GetCoulombAmplitude() const;
  complex GetECAmplitude(int, int) const;
  complex GetECAmplitudeWithShift(int, int, CNuc *, const Config &) const;
  EnergyMap GetMap() const;
  void Initialize(CNuc *, const Config &);
  double ConvertLabValue(double, PPair *);
  double ConvertCMValue(double, PPair *);
  void ConvertLabEnergy(PPair *);
  void ConvertExcitationEnergy(PPair *);
  void ConvertDecayEnergy(PPair *);
  void ConvertLabAngle(PPair *);
  void ConvertLabAngle(PPair *, PPair *, const Config &);
  void ConvertCMAngle(PPair *, PPair *, const Config &);
  void ConvertCrossSection(PPair *, PPair *);
  void ConvertLabAngleGammas(PPair *);
  void ConvertCrossSectionGammas(PPair *);
  double CalculateCrossSectionConversionFactor(PPair *, PPair *);
  double CalculateCrossSectionGammaConversionFactor(PPair *);
  void AddLegendreP(double);
  void ClearLegendrePolynomials();
  void SetGeometricalFactor(double);
  void SetFitCrossSection(double);
  void SetFitE1CrossSection(double);
  void SetFitE2CrossSection(double);
  void SetSFactorConversion(double);
  void SetLabEnergy(double);
  void SetCMEnergy(double);
  void SetExcitationEnergy(double);
  void SetLabAngle(double);
  void SetCMAngle(double);
  void SetExitKey(int);
  void SetEntranceKey(int);
  void CalcLegendreP(int, CNuc *, TargetEffect *);
  void CalcEDependentValues(CNuc *, const Config &);
  void RecalcEDependentValues(CNuc *, const Config &);
  void AddLoElement(int, int, complex);
  void AddSqrtPenetrability(int, int, double);
  void AddThmFormFactor(int, int, double, double);
  void AddExpCoulombPhase(int, int, complex);
  void AddExpHardSpherePhase(int, int, complex);
  void CalcCoulombAmplitude(CNuc *);
  void SetCoulombAmplitude(complex);
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
