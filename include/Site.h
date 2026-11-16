#ifndef AZURE2_SITE_H
#define AZURE2_SITE_H
#include "CalculationPlan.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>
#include <string>
#include <tuple>
#include <utility>
#include <cassert>
#include <memory>

#include <cstring>
#include <cstdint>
#include "Constants.h"
#include "SharedTable.h"
#include "ChannelTable.h"
#include "UnitTables.h"
#include "Results.h"
#include "TermSettings.h"
#include "UnitRange.h"
#include "SiteObservable.h"
#include "RowSet.h"

class Segment;
class WorkingCopy;
class CompoundSystem;
class LevelState;
class PPair;
class TargetEffect;
namespace azure2 {
namespace experiment {
struct AngleLimit;
struct FrameRow;
struct MeasurementRows;
}
}  // namespace azure2
class Config;
struct ThmWeightTable;
class ThmDistortion;
class ThmDwVertex;
struct ThmAngleWindow;
struct ThmCoherentBackground;
struct ThmLineshape;
struct ThmSpectatorWindow;


/*!
 * One calculation site: an energy and an angle at which the point pass solves the level matrix and writes an
 * observable -- an unfolded data point, or a slot of a folded one (a sample of its target-effect integral).  It holds
 * what depends on the point at that energy: the c.m., lab and compound energies, the angle and its Legendre rows,
 * the channel rows (shared with the node of the table it is at: NodeTable.h), the Coulomb amplitude, the screening
 * factor, the stopping, the external-capture amplitudes and where its result record lies.  The observables
 * (GenMatrixFunc, THMMatrixFunc, the external capture) read a site, whatever it is.  A data point (Unit) is a site
 * with its row, its fold and its records; the slots of a folded point are sites held by its working copy
 * (WorkingCopy::SubSites), one per node of its integral, in energy order.
 */
/// The identical-particle factor of dsigma/dOmega in the Legendre coefficients B_L of a reaction from `entrance`: 4 for
/// elastic scattering of identical particles (|f(theta) + f(pi - theta)|^2 of the nuclear amplitude), (1 + delta_12)
/// = 2 for another exit of an identical entrance pair, 1 otherwise (GenMatrixFunc::CalculateBareCrossSection).
double LegendreSymmetryFactor(const PPair &entrance, bool elastic);

class Site {
 public:
  /// A slot at a given c.m. energy and angle that is no row of a measurement (a sample of a target-effect integral),
  /// with the conventions of term `term` of the segment (Segment::TermEntranceKey ...).
  Site(double angle, double energy, const Segment *segment, int term);
 protected:
  using ThmPsTable = ::ThmPsTable;  // UnitTables.h
 public:

  using EDependentTables = ::EDependentTables;  // UnitTables.h
  using EDependentShare = ::EDependentShare;    // UnitTables.h
  /// The point's energy-dependent values, for another point at the same energy (requires them to be valid).
  EDependentShare ShareEDependentValues() const;
  /// Takes a share made at exactly this point's energy, angle and pairs; returns false and changes nothing otherwise.
  bool AdoptEDependentValues(const EDependentShare &);
  /// Differential cross section?
  bool IsDifferential() const;
  /// Phase shift?
  bool IsPhase() const;
  /// THM (modified R-matrix, half-off-shell) point?
  bool IsTHM() const { return settings_->thm; }
  /// Weight table multiplying the THM model cross section (<thm> weight[k]=),
  /// or null.  Owned by the parent Segment.
  const ThmWeightTable *GetThmWeight() const { return settings_->thmWeight; }
  /// Coulomb line shape of the point's THM experiment (lineshape=on), or
  /// null.  Owned by the parent Segment.
  const ThmLineshape *GetThmLineshape() const { return settings_->thmLineshape; }
  /// Spectator-momentum window of the point's THM experiment (ps=...), or
  /// null.  Owned by the parent Segment.
  const ThmSpectatorWindow *GetThmSpectatorWindow() const { return settings_->thmWindow; }
  /// Distortion factor R(E) of the point's THM experiment (distortion=...),
  /// or null.  Owned by the parent Segment.
  const ThmDistortion *GetThmDistortion() const { return settings_->thmDistortion; }
  /// Distorted-wave entrance vertex of the point's THM experiment
  /// (vertexModel=dw), or null (the plane-wave vertex).
  const ThmDwVertex *GetThmDwVertex() const { return settings_->thmDw; }
  /// Angular window of the point's THM experiment (theta=thmin-thmax), or
  /// null (the angle-integrated observable).  Owned by the parent Segment.
  const ThmAngleWindow *GetThmAngleWindow() const { return settings_->thmAngle; }
  /// THM coherent background of the point's experiment (cbackground=), or null.
  const ThmCoherentBackground *GetThmCoherent() const { return settings_->thmCoherent; }
  /// Angular distribution?
  bool IsAngularDist() const;
  //! Vector analyzing power point; the fit value is A_y, not a cross section.
  bool IsAnalyzingPower() const { return settings_->analyzingPower; };
  //! Half-width (keV) of the final-state line window of a capture into a broad unbound final state; 0: none.
  double GetLineWindow() const { return settings_->lineWindow; };
  //! The unit's whole record and its length: the line-shape integration sums the records of its nodes.
  double *Record(azure2::Results &r) const { return Rec(r); }
  std::size_t RecordLength() const { return azure2::Results::RecordSize(num_ang_slots_, num_natural_slots_); }
  //! Whether the record has an angular part (its count at Results::AngularOffset()).
  bool HasAngularPart() const { return num_ang_slots_ != 0; }
  /// The total cross section of the entrance pair (optical theorem), Segment::IsTotalTerm.
  bool IsTotal() const { return settings_->total; };
  //! A_y is kept beside the cross section rather than replacing it, because
  //! target-effect integration needs the cross section as the weight.
  double GetAnalyzingPower(const azure2::Results &r) const { return Rec(r)[azure2::Results::kAy]; };
  void SetAnalyzingPower(azure2::Results &r, double v) const { Rec(r)[azure2::Results::kAy] = v; };
  /// Is this one of the sub-points a target-effect integral is built from?
  bool IsSubPoint() const { return is_sub_point_; };
  /*!
   * Electron screening (a `screening` block of the point's target effect): the calculated
   * cross section is the bare one times f = exp(k U_e), k = pi eta / E_cm of the entrance pair
   * (CalcEDependentValues), U_e in MeV.  param is the 1-based number of the screening
   * parameter (WorkingCopy::ScreeningEffects), 0 = not screened.  Set by WorkingCopy::ApplyScreening on every
   * point and sub-point of the effect's segments; a new sub-point inherits it (WorkingCopy::AddSubPoint).
   */
  void SetScreening(double ueMeV, int param) { values_.screen_ue = ueMeV; screen_param_ = param; }
  int GetScreeningParam() const { return screen_param_; }
  double GetScreeningUe() const { return values_.screen_ue; }
  /// d ln f / d U_e (1/MeV); 0 for a neutral or non-particle entrance pair.
  double GetScreeningK() const { return values_.screen_k; }
  /// The factor f applied to this point's cross section (1 when not screened).
  double ScreeningFactor() const { return (screen_param_ > 0 && values_.screen_k != 0.0) ? std::exp(values_.screen_k * values_.screen_ue) : 1.0; }
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
  /// Highest Legendre order. Angular-distribution points only.
  int GetMaxAngDistOrder() const;
  /// An angular-coefficient data point (segment isDiff 9): the fit value is the Legendre
  /// coefficient of this order, a_k = A_k / A_0 (or A_k itself, b/sr, when absolute).
  bool IsCoefficient() const;
  int GetCoefficientOrder() const;
  bool GetCoefficientAbsolute() const;
  /// A_0 = sigma_int / 4 pi (b/sr) of an angular-distribution point, set with its coefficients.
  double GetAngDistA0(const azure2::Results &r) const { return Rec(r)[azure2::Results::kA0]; };
  void SetAngDistA0(azure2::Results &r, double a0) const { Rec(r)[azure2::Results::kA0] = a0; };
  /// Angular momentum of the secondary decay. UPOS points only.
  int GetSecondaryDecayL() const;
  /// The kind of observable the site's record holds (SiteObservable.h), from its term.
  SiteObservable Observable() const;
  /// Whether the record holds observables that its readers contract with their own angular rows (ValueFor): B_L, for
  /// elastic scattering with the Coulomb terms, A_y (ReadThroughRow).  A scalar kind (kValue, kHoes) is read as it is.
  bool AngleFree() const { return ReadThroughRow(Observable()); }
  /*!
   * The observable of this site as `reader` sees it at the reader's slot `slot` (1-based; 0: the reader is the site):
   * the reader's angular row contracted with the record (the row of the site's kind in kSiteObservables); for kValue
   * the value the record holds (the reader's angle is the site's).
   */
  double ValueFor(const Site &reader, int slot, const azure2::Results &r) const;
  /// A_y of this site as `reader` sees it (kAnalyzingPower: the second value of the contraction), else the record's.
  double AnalyzingPowerFor(const Site &reader, int slot, const azure2::Results &r) const;
  /// The values of this site at slot k (0-based) of the row set `set` (RowSet.h) into values[0 .. set.valuesPerSlot):
  /// the set's row contracted with the record, as ValueFor (and AnalyzingPowerFor) give them to each reader of the set,
  /// the first divided by the slot's Jacobian (RowSet::JacobianAt: 1 but for data given at a lab angle).
  void ValuesWith(const azure2::RowSet &set, int k, const azure2::Results &r, double *values) const;
  /// The reader's Coulomb amplitudes at its slot `slot` (0: its own energy): f_d, the (symmetrized) f_C(theta), or
  /// the direct one for identical particles with spin, and f_x = f_C(pi - theta) there, else 0 (a reader of elastic
  /// slots reads them from its row set: RowSet.h).
  void CoulombAt(int slot, complex &direct, complex &exchange) const;
  /// f_C at c.m. energy `energy` (MeV) and c.m. angle `angle` (deg) of `pair` scattering elastically: `direct`
  /// (symmetrized for an identical pair without spin-dependent exchange) and `exchange` (0 unless spin-dependent).
  static void CoulombAmplitudesAt(const PPair *pair, double energy, double angle, complex &direct, complex &exchange,
                                  bool &spinDependent);
  /// The reader's angular row at its slot `slot` (1-based; 0: its own energy): P_L(cos theta) with the Q coefficients
  /// of its effect, at the slot's energy where the detector's lab transform makes it depend on the energy
  /// (RowDependsOnEnergy; from its row set, RowSet.h), else at the point's.
  double RowP(int slot, int l) const { return Row(slot)[l]; }
  /// The whole row at slot `slot` (RowP(slot, l) = Row(slot)[l]).
  const double *Row(int slot) const {
    return (slot > 0 && values_.rowset) ? values_.rowset->Row(values_.rowset_at + slot - 1) : values_.legendreP->data();
  }
  /// The row set the reader's slots are read with (RowSet.h), null where it reads no angle-free slots.
  const azure2::RowSet *RowSetOf() const { return values_.rowset ? &*values_.rowset : nullptr; }
  /// The position of the reader's slot `slot` (1-based) in its row set.
  int RowSetIndex(int slot) const { return values_.rowset_at + slot - 1; }
  /// Read the slots with the row set `set` from its position `at` (WorkingCopy::BuildRowSets).
  void ReadWithRowSet(SharedTable<azure2::RowSet> set, int at) {
    values_.rowset = std::move(set);
    values_.rowset_at = at;
  }
  /*!
   * The value of this reader's slot `slot` (1-based), the site `site`: the contraction of the evaluation where the
   * reader has a row set (Results::contracted, RowSet.h), else site.ValueFor(*this, slot, r) -- the same number.
   */
  double SlotValue(const Site &site, int slot, const azure2::Results &r) const {
    if (values_.rowset) return r.contracted[values_.rowset->ValueIndex(values_.rowset_at + slot - 1, 0)];
    return site.ValueFor(*this, slot, r) / SlotJacobian(slot);
  }
  /// The reader's c.m. angle (deg) at its slot `slot` (1-based; 0: its own): its row set's (RowSet::AngleAt; for data
  /// given at a lab angle theta_cm(E_n; theta_lab)), else its own.
  double SlotAngle(int slot) const {
    return (slot > 0 && values_.rowset) ? values_.rowset->AngleAt(values_.rowset_at + slot - 1) : GetCMAngle();
  }
  /// dOmega_lab/dOmega_cm at the reader's slot `slot` (RowSet::JacobianAt): SlotValue is the contraction divided by
  /// it, 1 but for data given at a lab angle.
  double SlotJacobian(int slot) const { return (slot > 0 && values_.rowset) ? values_.rowset->JacobianAt(values_.rowset_at + slot - 1) : 1.0; }
  /// A_y of this reader's slot `slot` (an analysing power; SlotValue is the cross section that weights it).
  double SlotAnalyzingPower(const Site &site, int slot, const azure2::Results &r) const {
    if (values_.rowset && values_.rowset->valuesPerSlot > 1)
      return r.contracted[values_.rowset->ValueIndex(values_.rowset_at + slot - 1, 1)];
    return site.AnalyzingPowerFor(*this, slot, r);
  }
  /// Whether the point's angular row depends on the energy: Q coefficients on a particle-to-particle reaction that is
  /// not charged elastic scattering, UPOS or A_y (the detector averages the lab distribution: CalcLegendreP).
  bool RowDependsOnEnergy(const CompoundSystem *, const TargetEffect *) const;
  /// The row P_L Q-transformed at c.m. energy `cmEnergy` (MeV) and c.m. angle `cmAngle` (deg; the point's own where
  /// not given), L = 0 .. maxL, into p.
  void LegendreRowAt(double cmEnergy, int maxL, const CompoundSystem *, const TargetEffect *, double *p) const;
  void LegendreRowAt(double cmEnergy, double cmAngle, int maxL, const CompoundSystem *, const TargetEffect *, double *p) const;
  /// The number K of Legendre coefficients the record holds (AngleFree) and the first of them.
  int NumLegendreCoefficients(const azure2::Results &r) const {
    return num_ang_slots_ ? (int)Rec(r)[azure2::Results::AngularOffset()] : 0;
  }
  const double *LegendreCoefficients(const azure2::Results &r) const {
    return Rec(r) + azure2::Results::AngularOffset() + 1;
  }
  /// a_k A_0 = B_k / f of order k as `reader`'s angular distribution has it (0 beyond its max_ang_dist_order or the
  /// coefficients held).
  double CoefficientOfOrder(const Site &reader, const azure2::Results &r, int k) const;
  /// Multiply every observable the angular part of the record holds by `f` (electron screening).
  void ScaleObservable(azure2::Results &r, double f) const;
  /// NaN in every observable the angular part holds, with counts that make every reader NaN (MarkFailed).
  void PoisonObservable(azure2::Results &r) const;
  /// The identical-particle factor of the point's B_L (LegendreSymmetryFactor), by which they exceed the a_k A_0 of
  /// an angular distribution (sigma_int without it); set with the point's Legendre rows (CalcLegendreP).
  double EntranceSymmetry() const { return values_.entrance_symmetry; }
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
  /// The point's row in its segment's measurement (0-based), -1 for a point that is no row (a sub-point).
  int Row() const { return is_sub_point_ ? -1 : row_; }
  /// Give the point its term's settings (Segment::SetTermSettings, at set-up).
  void UseSettings(const TermSettings *settings) { settings_ = settings; }
  /// The settings of the site's term (shared by its term's units).
  const TermSettings &Settings() const { return *settings_; }
  /// Make this point a sub-point of `parent`'s row: it takes the row's Legendre order (an angular-coefficient row)
  /// but none of its data.
  /// Make this site a slot of `point` (Unit::MakeSubPoint): it takes the point's screening.
  void TakeSlotSettingsOf(const Site &point);
  void BelongToRowOf(const Site &parent) {
    // a slot holding Legendre coefficients serves every order: the reader's row picks its own (ValueFor)
    legendre_order_ = parent.AngleFree() ? -1 : parent.legendre_order_;
    row_ = parent.row_;
    is_sub_point_ = true;
  }
  /// Geometrical factor \f$\pi/k^2\f$.
  double GetGeometricalFactor() const;
  /// Cross section AZURE2 calculated here.
  double GetFitCrossSection(const azure2::Results &r) const { return Rec(r)[azure2::Results::kSigma]; }
  /// Calculated cross section, E1 component only.
  double GetFitE1CrossSection(const azure2::Results &r) const { return Rec(r)[azure2::Results::kSigmaE1]; }
  /// Calculated cross section, E2 component only.
  double GetFitE2CrossSection(const azure2::Results &r) const { return Rec(r)[azure2::Results::kSigmaE2]; }
  /// Multiply a cross section by this to get the astrophysical S-factor.
  double GetSFactorConversion() const;
  /// \f$\sqrt{P_c}\f$ for the channel at (J-group, channel), both 1-based.
  double GetSqrtPenetrability(int, int) const;
  /// THM entrance transfer form factor M_l = (b-1) j_l - rho dj_l/drho + C_l for
  /// JGroup/channel, assembled with the given (possibly complex) boundary b;
  /// C_l is the external Coulomb term, 0 unless requested (0 if not stored).
  complex GetThmFormFactor(int, int, complex) const;
  /// The same at node `node` of the point's spectator-momentum window
  /// (0-based; ps=...).  0 if not stored.
  complex GetThmFormFactor(int, int, complex, int node) const;
  /// Number of nodes of the point's spectator-momentum window (0: none).
  int NumThmPsNodes() const { return values_.thm_ps ? values_.thm_ps->nodes : 0; }
  /// Normalized weight of node `node` of that window at the point's energy.
  double GetThmPsWeight(int node) const { return values_.thm_ps->weight[node]; }
  /// Bytes held by the point's spectator-window table (0 without one; shared
  /// tables are counted by every point that refers to them).
  size_t ThmPsTableBytes() const { return values_.thm_ps ? values_.thm_ps->Bytes() : 0; }
  /// Total spin. Phase-shift points only.
  double GetJ() const;
  /// Channel spin of a phase-shift point (< 0: the segment gave none, so a single spin must match).
  double GetPhaseS() const { return settings_->phaseS; }
  /// Stopping cross section at this sub-point, for a yield-curve target integration.
  double GetStoppingPower() const;
  /// Beam energy loss across the target at this point.
  double GetTargetThickness() const;
  /// The lab energy the target thickness was computed at: the nominal one at the set-up (WorkingCopy::
  /// SetUpTargetEffects), the shifted, true one after an energy shift (azure2::MoveTargetToShift), so a fitted density
  /// or stopping scale recomputes it at the same energy (WorkingCopy::FillDensitiesFromParams, FollowStoppingScales).
  double GetTargetLabEnergy() const { return values_.targetLabEnergy; }
  void SetTargetLabEnergy(double e) { values_.targetLabEnergy = e; }
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
  /// Coulomb (Rutherford) amplitude \f$C_\alpha\f$.  For an identical spin-0
  /// pair this is already the symmetrized f_C(theta) + f_C(pi - theta); for an
  /// identical pair with spin it is only the direct term -- use the
  /// channel-spin overload there.
  complex GetCoulombAmplitude() const;
  /// Coulomb amplitude seen by channel spin \p s.  Identical to
  /// GetCoulombAmplitude() except for an identical pair with spin j != 0,
  /// where exchange enters with (-1)^s: f_C(theta) + (-1)^s f_C(pi - theta).
  complex GetCoulombAmplitude(double s) const;
  /// External-capture amplitude of external pathway m of spin group k of the point's exit (PathwayTable.h), both
  /// 1-based.
  complex GetECAmplitude(int, int) const;
  /// As GetECAmplitude, at the point's current (shifted) energy: from the energy tables
  /// (ECTable), or interpolated through the amplitude cache with --no-ec-table.  A value computed from the tables
  /// is held in the site until its energy moves, so the caller's site is written, as by RecalcEDependentValues.
  complex GetECAmplitudeWithShift(int, int, const CompoundSystem *, const Config &);
  /// Fill the compound excitation energy from the c.m. energy (a sub-point moved by an energy shift).
  void ConvertExcitationEnergy(const PPair *);
  /// Relativistic lab -> c.m. angle and dOmega_lab/dOmega_cm of the detected particle (relativistic=1).
  bool RelativisticLabToCM(const PPair *entrancePair, const PPair *exitPair, double nrCMAngle, double &cmAngle, double &factor) const;
  /// dOmega_lab/dOmega_cm of a particle at the point, not stored (the RATIO frame factors).
  double CalculateCrossSectionConversionFactor(const PPair *, const PPair *);
  /// Lab motion (g, boost gamma) of the detected exit particle for the attenuation coefficients at c.m. energy
  /// `cmEnergy` (MeV); false: Q as entered.
  bool QLabBoostAt(double cmEnergy, const PPair *entrancePair, const PPair *exitPair, double &g, double &gb) const;
  /// dOmega_lab/dOmega_cm of a photon at the point (the RATIO frame factors).
  double CalculateCrossSectionGammaConversionFactor(const PPair *);
  /// Append the next Legendre polynomial; they are stored in order L = 0, 1, 2, ...
  void AddLegendreP(double);
  /// Drop the stored polynomials, so they can be recomputed after an energy shift.
  void ClearLegendrePolynomials();
  /// Set \f$\pi/k^2\f$.
  void SetGeometricalFactor(double);
  /// Record the calculated cross section.
  void SetFitCrossSection(azure2::Results &r, double v) const { Rec(r)[azure2::Results::kSigma] = v; }
  /// The abundance-weighted term of each compound nucleus at this point of a
  /// natural-target segment, indexed by nucleus (azure2::EvaluateModelPoint).
  std::vector<double> GetNaturalTerms(const azure2::Results &r) const;
  void SetNaturalTerms(azure2::Results &r, const std::vector<double> &terms);
  void SetFitE1CrossSection(azure2::Results &r, double v) const { Rec(r)[azure2::Results::kSigmaE1] = v; }
  void SetFitE2CrossSection(azure2::Results &r, double v) const { Rec(r)[azure2::Results::kSigmaE2] = v; }
  /// Where the point's result record lies in the buffer of an evaluation (Results.h, docs/upgrade/data_model.md S3):
  /// its offset, assigned by the layout of the data set that holds the point (WorkingCopy::LayoutResults) or by the caller
  /// of a calculation on a point outside a data set (UsePrivateRecord); kNoRecord before either.  The point holds
  /// no buffer: the values are in the caller's.
  static constexpr std::size_t kNoRecord = static_cast<std::size_t>(-1);
  std::size_t ResultSlot() const { return slot_; }
  /// Slots of the angular part a site's record needs where the reactions' Legendre orders are at most `orders`
  /// and the amplitudes of an analysing power are at most `capacity` (WorkingCopy::Capacity): its kind's
  /// (SiteObservable.h), 0 for kValue.
  int ResultAngularCapacity(const RecordCapacity &capacity) const;
  /// Slots for the coefficients a_k of an angular-distribution or coefficient point in its model record:
  /// max_ang_dist_order + 1, else 0.
  int ModelAngularCapacity() const {
    return (settings_->angularDist && settings_->maxAngDistOrder >= 0) ? settings_->maxAngDistOrder + 1 : 0;
  }
  int ResultAngularSlots() const { return num_ang_slots_; }
  int ResultNaturalSlots() const { return num_natural_slots_; }
  /// Place the point's record at offset `slot`, with nAng angular and nNat natural slots.
  void SetResultRecord(std::size_t slot, int nAng, int nNat);
  /// Give the point (a point outside any data set) the only record of `scratch`, at offset 0, holding zeros (NaN under
  /// AZURE_POISON_RESULTS=1: the record of a new evaluation, as WorkingCopy::BeginEvaluation gives).
  void UsePrivateRecord(azure2::Results &scratch);
  /// A calculation that threw: NaN cross sections (total, E1, E2) on this point, and a warning naming the segment, the
  /// role of the point, its energy and the error, the first time in a run.  Used by azure2::EvaluateModelPoint and
  /// WorkingCopy::CalculatePointsInParallel (before 2026-11-01 both gave 0 or left the term out, silently).
  void MarkFailed(azure2::Results &results, const char *role, const Config &configure, const char *why);
  int GetSegmentKey() const { return segment_key_; }
  /// Set the cross-section to S-factor conversion.
  void SetSFactorConversion(double);
  void SetLabEnergy(double);
  void SetCMEnergy(double);
  void SetLabAngle(double);
  void SetCMAngle(double);
  /// Compute the Legendre polynomials up to a maximum order.  `data`: the working copy of the site (null for a point
  /// outside one).  These functions write this site only; Unit's functions of the same names reach its slots too.
  void CalcLegendreP(WorkingCopy *data, int, const CompoundSystem *, const TargetEffect *);
  /// Compute the geometrical factor, S-factor conversion, \f$L_o\f$ elements, penetrabilities and phases together.
  /// The rows of the open and closed particle channels are read from the tables of `rows` (docs/upgrade/
  /// channel_table.md) when one holds them, and otherwise computed (ChannelTable::Compute), once per (pair, l, E_c)
  /// of the point; rows->used, when set, receives every row used.  `rowsOf`, when given, is a site with the same
  /// ChannelRowsKey whose rows are computed: this site holds that table instead of computing its own (the rest --
  /// the conversion factors, the THM vertex pieces -- is its own).
  void CalcEDependentValues(WorkingCopy *data, const CompoundSystem *, const Config &,
                            const ChannelRowSource *rows = nullptr, const Site *rowsOf = nullptr);
  /// What the channel rows (EDependentTables) depend on: the entrance pair's compound nucleus, the bits of the compound
  /// energy E_x (CompoundEnergy), and the exit pair's charge Z_1 + Z_2 where a beta-decay channel of that nucleus
  /// reads it (its Fermi function; -1 elsewhere).  Two sites with one key compute the same rows bit for bit.
  std::tuple<int, std::uint64_t, int> ChannelRowsKey(const CompoundSystem *);
  /// Hold the channel rows of `other`, a site with the same ChannelRowsKey whose rows are computed, in place of this
  /// site's own (equal bit for bit, dropped).
  void TakeChannelRowsOf(const Site &other) { values_.edep = other.values_.edep; }
  /// The keys (pair, l, E_c) of the channel-table rows CalcEDependentValues reads for this site.
  void ChannelKeys(WorkingCopy *data, const CompoundSystem *, std::vector<ChannelTable::Key> &);
  /// The compound-nucleus energy E_x (MeV) at which CalcEDependentValues evaluates the channels.
  double CompoundEnergy(const CompoundSystem *);
  /// Recompute those at the current, possibly shifted, energy; nothing while they are at that energy.
  void RecalcEDependentValues(WorkingCopy *data, const CompoundSystem *, const Config &, const ChannelRowSource *rows = nullptr);
  /// Are the energy-dependent values those of the site's present energy?
  bool EDependentCurrent() const { return values_.eDependentValid && this->GetCMEnergy() == values_.eDependentEnergy; }
  /// The channel rows the point holds (L_o, sqrt(P), phases per J group and channel; null before they are computed),
  /// and a holder of them (the node table shares them: NodeTable.h).
  const EDependentTables *ChannelRows() const { return values_.edep ? &*values_.edep : nullptr; }
  SharedTable<EDependentTables> ChannelRowsShared() const { return values_.edep; }
  /// Empties the per-JGroup rows of the energy-dependent arrays, keeping their storage.
  void ClearEDependentRows();
  /// Store an \f$L_o\f$ element at (J-group, channel).
  void AddLoElement(int, int, complex);
  /// Store a \f$\sqrt{P_c}\f$ at (J-group, channel).
  void AddSqrtPenetrability(int, int, double);
  /// Store a Coulomb phase factor at (J-group, channel).
  void AddExpCoulombPhase(int, int, complex);
  /// Store a hard-sphere phase factor at (J-group, channel).
  void AddExpHardSpherePhase(int, int, complex);
  /// Compute the Coulomb amplitude \f$C_\alpha\f$.
  void CalcCoulombAmplitude(WorkingCopy *data, const CompoundSystem *);
  /// Set the Coulomb amplitude directly.
  void SetCoulombAmplitude(complex);
  /// Compute the external-capture amplitudes for every matching pathway.
  /// `table`: the stage's channel table, for sqrt(P) and exp(-i phi) at a_c (null: computed).
  void CalculateECAmplitudes(WorkingCopy *data, const CompoundSystem *, const Config &, const ChannelTable *table = nullptr);
  /// Keep an external-capture amplitude of pathway (k, m) at the point's energy (or at `energy`, MeV c.m.); with
  /// --no-ec-table also in the run's cache (RunECAmplitudes).
  void AddECAmplitude(int, int, complex, const Config &);
  void AddECAmplitude(int, int, complex, double, const Config &);
  void ClearECAmplitudes();
  /// The external-capture amplitudes the point holds at c.m. energy `energy` (MeV) -- the set-up's, which
  /// GetECAmplitudeWithShift returns there instead of computing them; they need not agree to the last bit between
  /// points at one energy -- as a 64-bit hash of their positions and bits (FNV-1a), 0 when it holds none there.
  std::uint64_t ECAmplitudesHeldAt(double energy) const {
    std::uint64_t h = 1469598103934665603ULL;
    bool any = false;
    auto mix = [&h](const void *p, std::size_t n) {
      for (std::size_t i = 0; i < n; i++) h = (h ^ static_cast<const unsigned char *>(p)[i]) * 1099511628211ULL;
    };
    for (std::size_t k = 0; k < ec_energies_.size(); k++)
      for (std::size_t m = 0; m < ec_energies_[k].size(); m++) {
        if (ec_energies_[k][m] != energy || k >= ec_amplitudes_.size() || m >= ec_amplitudes_[k].size()) continue;
        any = true;
        mix(&k, sizeof k), mix(&m, sizeof m), mix(&ec_amplitudes_[k][m], sizeof(complex));
      }
    return any ? (h ? h : 1) : 0;
  }
  /// Whether `other` holds at its c.m. energy bit for bit the amplitudes this point holds at its own
  /// (ECAmplitudesHeldAt).
  bool SameECAmplitudesAs(const Site &other) const {
    auto held = [](const Site &u) {
      const double energy = u.GetCMEnergy();
      std::vector<std::pair<std::pair<std::size_t, std::size_t>, complex>> v;
      for (std::size_t k = 0; k < u.ec_energies_.size(); k++)
        for (std::size_t m = 0; m < u.ec_energies_[k].size(); m++)
          if (u.ec_energies_[k][m] == energy && k < u.ec_amplitudes_.size() && m < u.ec_amplitudes_[k].size())
            v.push_back({{k, m}, u.ec_amplitudes_[k][m]});
      return v;
    };
    const auto a = held(*this), b = held(other);
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); i++)
      if (a[i].first != b[i].first || std::memcmp(&a[i].second, &b[i].second, sizeof(complex)) != 0) return false;
    return true;
  }
  /// The test hook of failed points: true when AZURE_TEST_FAIL_E (and AZURE_TEST_FAIL_EXIT, if set) name this point.
  bool InjectedFailure() const;
  void SetStoppingPower(double);
  void SetTargetThickness(double);
  /// Record the Legendre coefficients B_0 .. B_{K-1} (b/sr) of the site (kLegendre, and the nuclear part of kElastic).
  void SetLegendreCoefficients(azure2::Results &r, const vector_r &coefficients) const;
  /// The observables of elastic scattering at the site's energy (SiteObservable::kElastic).
  struct ElasticCoefficients {
    vector_r nuclear;     ///< B_L^NN (b/sr)
    double coulomb = 0.0; ///< c_0
    double cross = 0.0;   ///< c_x
    vector_c direct;      ///< D_l (b/sr)
    vector_c exchange;    ///< X_l (b/sr), empty without spin-dependent exchange
  };
  void SetElasticCoefficients(azure2::Results &r, const ElasticCoefficients &e) const;
  /// The coefficients of an analysing power (SiteObservable.h, ay_record), written by
  /// GenMatrixFunc::AnalyzingPowerCoefficientsOf.
  struct AnalyzingPowerCoefficients {
    bool elastic = false;  ///< the cross section's coefficients are weightElastic's, else weightLegendre's
    ElasticCoefficients weightElastic;
    vector_r weightLegendre;
    bool capture = true;
    std::array<double, 4> spins{};  ///< particle: j1, j2 (entrance; j1 the projectile's), j1', j2'
    int amplitudes = 0, orders = 0;
    vector_c coefficients;          ///< particle: C(a, l') at a orders + l'
    vector_r denominator, numerator;  ///< capture: d_k, n_k
  };
  void SetAnalyzingPowerCoefficients(azure2::Results &r, const AnalyzingPowerCoefficients &a) const;
  /// A value without an angular dependence (B_0 alone) in the record of an angle-free site, laid out as its kind's.
  void SetIsotropicObservable(azure2::Results &r, double value) const;
  /// Write coefficients into the angular part `angular` of a record (its count first) laid out as kLegendre's or
  /// kElastic's (SiteObservable.h); the sizes are the caller's to check.
  static void WriteLegendre(double *angular, const vector_r &coefficients);
  static void WriteElastic(double *angular, const ElasticCoefficients &e);
 protected:
  // The members are in three groups, by who writes them (docs/upgrade/evaluation.md section 7.8): what the set-up
  // fixes, the values an evaluation's hand-over writes (Values, below), and the place of the result record.
  //
  // Fixed at set-up.  No evaluation writes these of a point; a slot's are fixed when the slot is made (at set-up,
  // or by the evaluation that re-places its grid: the slot is a new site then).
  bool is_sub_point_ = false;
  int segment_key_;
  double lab_angle_;       ///< the lab angle the point is calculated at: its row's, a fixed angle, or a parent's
  double nominal_energy_;  ///< a point that is no row: its energy before any shift (MeV)
  /// What the point's term computes (shared, immutable once set up; TermSettings.h).
  /// Not owning: the term (segment or leaf) holds the record by a shared pointer that every copy of the data set
  /// copies, so a copy of a unit touches no reference count.  A point outside a data set owns its own (own_settings_).
  const TermSettings *settings_ = nullptr;
  int row_ = -1;            ///< the point's row in its measurement (0-based), -1 for a point that is no row
  int legendre_order_ = -1; ///< the Legendre order k of an angular-coefficient row (and its sub-points), else -1
  int screen_param_ = 0;    ///< the screening parameter of the point's target effect (SetScreening: the same each time)
  // The external-capture amplitudes of the set-up and their energies (CalculateECAmplitudes).  A working copy made
  // by WorkingCopy::Clone starts without them and evaluates at its points' energies (Values::ec_shifted).
  matrix_c ec_amplitudes_;
  matrix_r ec_energies_;  // Energies at which EC amplitudes were calculated
  /// The amplitude at the last shifted energy of each pathway (GetECAmplitudeWithShift).
  struct ShiftedECAmplitude {
    double energy = 0.0;
    complex amplitude = complex(0.0, 0.0);
    bool valid = false;
  };
  /// What an evaluation writes into a site: the values of the hand-over (the energies and the c.m. angle after an
  /// energy shift, the channel rows and Coulomb amplitude at that energy, the screening energy, the stopping power
  /// and target thickness of a fitted density or stopping scale, the rows a re-placed grid is read with) and the
  /// external-capture amplitude an evaluation computes at a shifted energy.  Each is first written at set-up; the
  /// set-up's value is the evaluation's until a parameter of the hand-over moves it.  One struct, so that the part of
  /// a site that belongs to an evaluation is one named thing (section 7.8 has the writer of each member).
  struct Values {
    double cm_angle;
    double cm_energy;
    double lab_energy;
    double excitation_energy;
    double geofactor;
    double sfactorconv;
    double entrance_symmetry = 1.0;
    double screen_k = 0.0;
    double screen_ue = 0.0;
    double stoppingPower;
    double targetThickness;
    double targetLabEnergy = 0.0;
    complex coulombamplitude;
    // Direct and exchange Coulomb amplitudes, kept apart only for an identical
    // pair with spin, whose exchange sign depends on the channel spin.
    bool spinDependentCoulomb = false;
    complex coulombdirect = complex(0., 0.);
    complex coulombexchange = complex(0., 0.);
    // The Legendre polynomials P_l(cos theta) of the point, by order, and the energy-dependent channel tables (one
    // row per J group: L_o, sqrt(P), Coulomb and hard-sphere phase factors).  Filled at set-up and rewritten only
    // when the point's energy or grid changes, so the working copies share them (a copy of a folded project with
    // ~10^6 sub-points held ~7 GB of these rows); every writer goes through MutableEDep() / MutableLegendre(), which
    // give the point its own copy first if the tables are shared (copy-on-write).
    SharedTable<vector_r> legendreP;
    /// A reader of angle-free slots: the rows it reads them with (RowSet.h) and the position of its first slot there.
    SharedTable<azure2::RowSet> rowset;
    int rowset_at = 0;
    SharedTable<EDependentTables> edep;
    std::shared_ptr<const ThmPsTable> thm_ps;
    /// The THM vertex pieces of the entrance pair's channels (a THM point only; UnitTables.h).
    std::shared_ptr<const ThmVertexRows> thm_vertex;
    // Energy at which CalcEDependentValues last ran: RecalcEDependentValues is
    // a no-op while the point's energy is unchanged (its values depend only on
    // the energy, the angle and the channel radii).
    double eDependentEnergy = 0.0;
    bool eDependentValid = false;
    std::vector<std::vector<ShiftedECAmplitude>> ec_shifted;
  };
  Values values_;
  /// This point's own, writable tables (copied first if shared with another point).
  EDependentTables &MutableEDep();
  /// This point's own, writable Legendre polynomials (copied first if shared).
  vector_r &MutableLegendre();
  // Where the result record lies (see ResultSlot): no value of an evaluation is held in the point itself.  Laid out at
  // set-up (WorkingCopy::LayoutResults) and again by an evaluation that re-places a grid.
  std::size_t slot_ = kNoRecord;
  std::uint16_t num_ang_slots_ = 0;
  std::uint16_t num_natural_slots_ = 0;
  // The point's record in the buffer `r`; a point without a record, or a buffer too short for it (not sized for the
  // layout of the point's data set: WorkingCopy::BeginEvaluation), throws.
  const double *Rec(const azure2::Results &r) const {
    if (slot_ == kNoRecord || slot_ + azure2::Results::RecordSize(num_ang_slots_, num_natural_slots_) > r.slot.size())
      NoResultRecord(r);
    return r.slot.data() + slot_;
  }
  double *Rec(azure2::Results &r) const {
    if (slot_ == kNoRecord || slot_ + azure2::Results::RecordSize(num_ang_slots_, num_natural_slots_) > r.slot.size())
      NoResultRecord(r);
    return r.slot.data() + slot_;
  }
  [[noreturn]] void NoResultRecord(const azure2::Results &r) const;
 protected:
  Site() = default;
};

#endif