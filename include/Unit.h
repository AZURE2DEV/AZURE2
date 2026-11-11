#ifndef UNIT_H
#define UNIT_H
#include "CalculationPlan.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <string>
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

class Segment;
class WorkingCopy;
class Unit;
typedef Unit *UnitIterator;
class CNuc;
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


/// An AZURE data point

/*!
 * A data point object in AZURE consists of a defined entrance and exit pair, an energy, an angle,
 * measured cross section and uncertainty, s-factor conversions, and
 * several flags that determine the type of data (angle integrated or differential) to be analysed.
 */

class Unit {
 private:
  using ThmPsTable = ::ThmPsTable;  // UnitTables.h

 public:
  using EDependentTables = ::EDependentTables;  // UnitTables.h
  using EDependentShare = ::EDependentShare;    // UnitTables.h
  /// The point's energy-dependent values, for another point at the same energy (requires them to be valid).
  EDependentShare ShareEDependentValues() const;
  /// Takes a share made at exactly this point's energy, angle and pairs; returns false and changes nothing otherwise.
  bool AdoptEDependentValues(const EDependentShare &);
  /// The point of row `row` (0-based) of a measurement, taking its conventions from the segment.  The point reads
  /// what was measured there (energy, value, uncertainty, bins, file record, Legendre order) from the rows, which the
  /// segment's measurement keeps alive and never changes.
  Unit(const azure2::experiment::MeasurementRows &rows, int row, const Segment *segment, int term);
  /// A point at a given c.m. energy and angle that is no row of a measurement (a sub-point of a target effect), with
  /// the conventions of term `term` of the segment (Segment::TermEntranceKey ...).
  Unit(double angle, double energy, const Segment *segment, int term);
  /// Build with every convention given explicitly, rather than inherited from a segment.
  Unit(double, double, int, int, bool, bool, bool, double, int, int);
  /// Differential cross section?
  bool IsDifferential() const;
  /// Phase shift?
  bool IsPhase() const;
  /// THM (modified R-matrix, half-off-shell) point?
  bool IsTHM() const;
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
  /// Does the point carry target effects?
  bool IsTargetEffect() const;
  /// Angular distribution?
  bool IsAngularDist() const;
  //! Vector analyzing power point; the fit value is A_y, not a cross section.
  bool IsAnalyzingPower() const { return settings_->analyzingPower; };
  /// The total cross section of the entrance pair (optical theorem), Segment::IsTotalTerm.
  bool IsTotal() const { return settings_->total; };
  //! A_y is kept beside the cross section rather than replacing it, because
  //! target-effect integration needs the cross section as the weight.
  double GetAnalyzingPower(const azure2::Results &r) const { return Rec(r)[azure2::Results::kAy]; };
  void SetAnalyzingPower(azure2::Results &r, double v) const { Rec(r)[azure2::Results::kAy] = v; };
  /// How the adaptive sub-point grid of a point was laid out (WorkingCopy::SetUpTargetEffects), so WorkingCopy::
  /// ReplaceSubPointGrids can lay it out again for other level parameters; held by the data set (DataSet::PointSetUp).
  /// Energies are c.m., unshifted.
  struct SubGridSpec {
    bool valid = false;
    double startEnergy = 0.0, endEnergy = 0.0;
    int maxPoints = 0, entranceKey = 0;
    double baseEnergyStep = 0.0, resonanceWidthMultiplier = 0.0, pointsPerWidth = 0.0;
    bool targetIntegration = false;
    double stoppingConversion = 1.0;  ///< c.m./lab factor of the stopping power
  };
  /// What the set-up worked out for this point that its evaluations read again (its sub-point grid's layout, its
  /// user-defined resolution kernel): record SetUpIndex of the data set's DataSet::PointSetUps, -1 for none.
  int SetUpIndex() const { return set_up_; }
  void SetSetUpIndex(int i) { set_up_ = i; }
  /// The point's own lab -> c.m. factor and entrance-pair charges, with which every target-effect kernel
  /// converts the line's lab energies (beam sigma, beam profile, Bohr straggling).  Set by
  /// WorkingCopy::SetUpTargetEffects for each point and component point; a line shared by several pairs
  /// used to be converted once with the first pair's factor (2026-10-27).
  void SetTargetKinematics(double cmFactor, int z1, int z2) { cm_factor_ = cmFactor; bohr_z1_ = z1; bohr_z2_ = z2; }
  double GetCMFactor() const { return cm_factor_; }
  int GetBeamCharge() const { return bohr_z1_; }
  int GetTargetCharge() const { return bohr_z2_; }
  /// Is this one of the sub-points a target-effect integral is built from?
  bool IsSubPoint() const { return is_sub_point_; };
  /*!
   * Electron screening (a `screening` block on the point's <targetInt> line): the calculated
   * cross section is the bare one times f = exp(k U_e), k = pi eta / E_cm of the entrance pair
   * (CalcEDependentValues), U_e in MeV.  param is the 1-based number of the screening
   * parameter (WorkingCopy::ScreeningEffects), 0 = not screened.  Set by WorkingCopy::ApplyScreening on every
   * point and sub-point of the effect's segments; a new sub-point inherits it (WorkingCopy::AddSubPoint).
   */
  void SetScreening(double ueMeV, int param) { screen_ue_ = ueMeV; screen_param_ = param; }
  int GetScreeningParam() const { return screen_param_; }
  double GetScreeningUe() const { return screen_ue_; }
  /// d ln f / d U_e (1/MeV); 0 for a neutral or non-particle entrance pair.
  double GetScreeningK() const { return screen_k_; }
  /// The factor f applied to this point's cross section (1 when not screened).
  double ScreeningFactor() const { return (screen_param_ > 0 && screen_k_ != 0.0) ? std::exp(screen_k_ * screen_ue_) : 1.0; }
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
  /// Number of sub-points, the samples a target-effect integral is evaluated on.
  int NumSubPoints() const { return sub_count_; }
  /// 1-based position of the TargetEffect in the parent WorkingCopy.
  int GetTargetEffectNum() const;
  /// Weight in [0,1] with which the target effect applies to this point;
  /// 1 everywhere unless the effect declares energy ranges.
  double GetTargetBlendWeight() const { return target_blend_weight_; };
  void SetTargetBlendWeight(double w) { target_blend_weight_ = w; };
  /// Thin-target integration restricted to energy windows (2026-11-05): outside the windows (blend weight 0) the
  /// point is evaluated at the single depth-averaged energy E_bar = E - D of its own target instead of being a
  /// plain point at its data (beam) energy, which was off by the centroid shift ~Delta/2.  Such a point carries
  /// one sub-point at E_bar (fixed density) or three spanning the range of E_bar a fitted density can reach,
  /// and SetTargetFallbackOnly(true) marks it.  D (c.m. MeV, TargetEffect::MeanEnergyLoss) is set on every thin-target
  /// point: edge blending and the automatic decision fall back to sigma(E_bar) as well.
  void SetTargetFallbackOnly(bool f) { target_fallback_only_ = f; }
  bool IsTargetFallbackOnly() const { return target_fallback_only_; }
  /// Inside a TargetFallbackScope.
  bool IsTargetFallbackNow() const { return target_fallback_now_; }
  double GetTargetFallbackDepth() const { return target_fallback_depth_; }
  void SetTargetFallbackDepth(double d) { target_fallback_depth_ = d; }
  /// While alive, the target integration of the point returns sigma(E_bar) from its sub-points (an automatic
  /// tolerance decided the integration is not needed at this point).
  class TargetFallbackScope {
   public:
    explicit TargetFallbackScope(Unit *p) : p_(p), prev_(p->target_fallback_now_) { p_->target_fallback_now_ = true; }
    ~TargetFallbackScope() { p_->target_fallback_now_ = prev_; }
    TargetFallbackScope(const TargetFallbackScope &) = delete;
    TargetFallbackScope &operator=(const TargetFallbackScope &) = delete;
   private:
    Unit *p_;
    bool prev_;
  };
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
  /// Make this point a sub-point of `parent`'s row: it takes the row's Legendre order (an angular-coefficient row)
  /// but none of its data.
  void BelongToRowOf(const Unit &parent) {
    legendre_order_ = parent.legendre_order_;
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
  int NumThmPsNodes() const { return thm_ps_ ? thm_ps_->nodes : 0; }
  /// Normalized weight of node `node` of that window at the point's energy.
  double GetThmPsWeight(int node) const { return thm_ps_->weight[node]; }
  /// Bytes held by the point's spectator-window table (0 without one; shared
  /// tables are counted by every point that refers to them).
  size_t ThmPsTableBytes() const { return thm_ps_ ? thm_ps_->Bytes() : 0; }
  /// Total spin. Phase-shift points only.
  double GetJ() const;
  /// Channel spin of a phase-shift point (< 0: the segment gave none, so a single spin must match).
  double GetPhaseS() const { return settings_->phaseS; }
  /// Stopping cross section at this sub-point, for a yield-curve target integration.
  double GetStoppingPower() const;
  /// Beam energy loss across the target at this point.
  double GetTargetThickness() const;
  /// The lab energy the target thickness was computed at (WorkingCopy::SetUpTargetEffects, unshifted), so a
  /// fitted density recomputes it at the same energy (WorkingCopy::FillDensitiesFromParams).
  double GetTargetLabEnergy() const { return targetLabEnergy_; }
  void SetTargetLabEnergy(double e) { targetLabEnergy_ = e; }
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
  /// (ECTable), or interpolated through the amplitude cache with --no-ec-table.
  complex GetECAmplitudeWithShift(int, int, CNuc *, const Config &) const;
  /// Compute everything that depends on energy but not on the fit parameters. Run once before fitting.
  /// `data`: the data set that holds the point (its target effect may carry Q coefficients), null outside one.
  void Initialize(WorkingCopy *data, CNuc *, const Config &);
  /// Fill the compound excitation energy from the c.m. energy (a sub-point moved by an energy shift).
  void ConvertExcitationEnergy(PPair *);
  /// Take the energies, c.m. angle and c.m. bin edges of the point's row from its segment's frame
  /// (Segment::ConvertToCM, experiment::ConvertRow).
  void SetFrame(const azure2::experiment::FrameRow &row);
  /// Relativistic lab -> c.m. angle and dOmega_lab/dOmega_cm of the detected particle (relativistic=1).
  bool RelativisticLabToCM(PPair *entrancePair, PPair *exitPair, double nrCMAngle, double &cmAngle, double &factor) const;
  /// dOmega_lab/dOmega_cm of a particle at the point, not stored (the RATIO frame factors).
  double CalculateCrossSectionConversionFactor(PPair *, PPair *);
  /// Lab motion (g, boost gamma) of the detected exit particle for the attenuation coefficients; false: Q as entered.
  bool QLabBoost(PPair *entrancePair, PPair *exitPair, double &g, double &gb) const;
  /// dOmega_lab/dOmega_cm of a photon at the point (the RATIO frame factors).
  double CalculateCrossSectionGammaConversionFactor(PPair *);
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
  /// Slots for angular-distribution coefficients a unit needs: max_ang_dist_order + 1 for an angular-distribution
  /// or coefficient point (GenMatrixFunc writes min(maxLOrder, max_ang_dist_order) + 1 of them), else 0.
  int ResultAngularCapacity() const {
    return (settings_->angularDist && settings_->maxAngDistOrder >= 0) ? settings_->maxAngDistOrder + 1 : 0;
  }
  int ResultAngularSlots() const { return num_ang_slots_; }
  int ResultNaturalSlots() const { return num_natural_slots_; }
  /// Place the point's record at offset `slot`, with nAng angular and nNat natural slots.
  void SetResultRecord(std::size_t slot, int nAng, int nNat);
  /// Give the point (a probe copy, a point outside any data set) the only record of `scratch`, at offset 0, holding
  /// the values of its record in `values` (a probe copy reads its original's record); without `values`, zeros (NaN
  /// under AZURE_POISON_RESULTS=1: the record of a new evaluation, as WorkingCopy::BeginEvaluation gives).
  void UsePrivateRecord(azure2::Results &scratch, const azure2::Results &values);
  void UsePrivateRecord(azure2::Results &scratch);
  /// Where the point's model record lies in Results::model (WorkingCopy::LayoutResults; kNoRecord for a component point, a
  /// sub-point or a point outside a data set), and how many angular-distribution coefficients it holds.  The model
  /// value is what the data set predicts at the point (ModelPass.h); the leaf values above are what the point itself
  /// calculates.
  void SetModelRecord(std::size_t slot, int nAng) { model_slot_ = slot, num_model_ang_ = (std::uint16_t)nAng; }
  double GetModelValue(const azure2::Results &r) const { return r.model[ModelRecord(r)]; }
  void SetModelValue(azure2::Results &r, double v) const { r.model[ModelRecord(r)] = v; }
  /// The point's covariance seed (M^-1 r)_i (Results::ObjectiveValues::covarianceSeed), at its model record.
  double GetCovarianceSeed(const azure2::Results &r) const { return r.objective.covarianceSeed.at(ModelRecord(r)); }
  void SetCovarianceSeed(azure2::Results &r, double z) const { r.objective.covarianceSeed.at(ModelRecord(r)) = z; }
  /// The coefficients a_k of an angular-distribution point as the model record holds them (count, then values).
  int GetNumModelAngularDists(const azure2::Results &r) const {
    if (!num_model_ang_) return 0;
    const double n = r.model[ModelRecord(r) + 1];  // NaN (poisoned, not written): the capacity, as the leaf's count
    return std::isnan(n) ? (int)num_model_ang_ : (int)n;
  }
  double GetModelAngularDist(const azure2::Results &r, int k) const { return r.model[ModelRecord(r) + 2 + k]; }
  /// Copy the leaf record's angular-distribution count and coefficients into the model record (EvaluateModelPoint).
  void CopyAngularDistsToModel(azure2::Results &r) const {
    if (!num_model_ang_) return;
    const double *from = Rec(r) + azure2::Results::AngularOffset();
    std::copy(from, from + 1 + std::min<int>(num_model_ang_, num_ang_slots_), r.model.begin() + ModelRecord(r) + 1);
  }
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
  /// Compute the Legendre polynomials up to a maximum order.  The functions below that reach the point's sub-points
  /// take `data`, the data set that holds them (null for a point outside one, which has none).
  void CalcLegendreP(WorkingCopy *data, int, CNuc *, const TargetEffect *);
  /// Compute the geometrical factor, S-factor conversion, \f$L_o\f$ elements, penetrabilities and phases together.
  /// With includeSubPoints false the sub-points are left alone (their energies are updated, and their values recomputed, separately).
  /// The rows of the open and closed particle channels are read from the tables of `rows` (docs/upgrade/
  /// channel_table.md) when one holds them, and otherwise computed (ChannelTable::Compute), once per (pair, l, E_c)
  /// of the point; rows->used, when set, receives every row used.
  void CalcEDependentValues(WorkingCopy *data, CNuc *, const Config &, bool includeSubPoints,
                            const ChannelRowSource *rows = nullptr);
  /// The keys (pair, l, E_c) of the channel-table rows CalcEDependentValues reads for this point and, with
  /// includeSubPoints, its sub-points: the same blocks, channels and channel energies.
  void ChannelKeys(WorkingCopy *data, CNuc *, std::vector<ChannelTable::Key> &, bool includeSubPoints = true);
  /// The compound-nucleus energy E_x (MeV) at which CalcEDependentValues evaluates the channels.
  double CompoundEnergy(CNuc *);
  /// Recompute those at the current, possibly shifted, energy.  Writes this point (and its sub-points) only.
  void RecalcEDependentValues(WorkingCopy *data, CNuc *, const Config &, bool includeSubPoints = true,
                              const ChannelRowSource *rows = nullptr);
  /// The channel rows the point holds (L_o, sqrt(P), phases per J group and channel; null before they are computed),
  /// and a holder of them (the node table shares them: NodeTable.h).
  const EDependentTables *ChannelRows() const { return edep_ ? &*edep_ : nullptr; }
  SharedTable<EDependentTables> ChannelRowsShared() const { return edep_; }
  /// Empties the per-JGroup rows of the energy-dependent arrays, keeping their storage.
  void ClearEDependentRows();
  /// Store an \f$L_o\f$ element at (J-group, channel).
  void AddLoElement(int, int, complex);
  /// Store a \f$\sqrt{P_c}\f$ at (J-group, channel).
  void AddSqrtPenetrability(int, int, double);
  /// Store the pieces of the THM form factor at (J-group, channel).
  void AddThmFormFactor(int, int, double, double, complex = complex(0.0, 0.0));
  /// Store a Coulomb phase factor at (J-group, channel).
  void AddExpCoulombPhase(int, int, complex);
  /// Store a hard-sphere phase factor at (J-group, channel).
  void AddExpHardSpherePhase(int, int, complex);
  /// Compute the Coulomb amplitude \f$C_\alpha\f$.
  void CalcCoulombAmplitude(WorkingCopy *data, CNuc *, bool includeSubPoints = true);
  /// Set the Coulomb amplitude directly.
  void SetCoulombAmplitude(complex);
  /// Compute the external-capture amplitudes for every matching pathway.
  /// `table`: the stage's channel table, for sqrt(P) and exp(-i phi) at a_c (null: computed).
  void CalculateECAmplitudes(WorkingCopy *data, CNuc *, const Config &, const ChannelTable *table = nullptr);
  void AddECAmplitude(int, int, complex);
  void AddECAmplitude(int, int, complex, double);
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
  bool SameECAmplitudesAs(const Unit &other) const {
    auto held = [](const Unit &u) {
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
  /// Calculate the point (its sub-points and its fold) for the current parameters, outside the energy pass
  /// (WorkingCopy::CalculatePointsInParallel), writing the values into `results`, the buffer of the evaluation laid out
  /// by `data`, the data set that holds the point (its target effects give the fold).  One solve at the point's energy,
  /// or one per sub-point (azure2::SolveAtOneEnergy), then FoldSubPoints.  Null `data` for a point outside any data
  /// set (no target effect: a reaction-rate point).
  void Calculate(WorkingCopy *data, CNuc *, const Config &configure, azure2::Results &results);
  /// The test hook of failed points: true when AZURE_TEST_FAIL_E (and AZURE_TEST_FAIL_EXIT, if set) name this point.
  bool InjectedFailure() const;
  /// Before its sub-points are calculated: with external capture, their energy-dependent values at their energies.
  void PrepareSubPoints(WorkingCopy *data, CNuc *, const Config &configure);
  /// The fold of this point from its sub-points' values already in `results`: the bare probe, the automatic tolerance
  /// and edge blending, the thin-target fallback, the integrations.
  void FoldSubPoints(WorkingCopy *data, CNuc *, const Config &configure, azure2::Results &results);
  /// The fold of a thin-target integration, with its sigma(E_bar) fallback (see the .cpp).
  void FoldThinTarget(WorkingCopy &data, azure2::Results &results, const Config &configure, const TargetEffect *effect);
  /// The task of the data set's calculation plan that calculates this point (CalculationPlan.h; -1: none, as for a
  /// sub-point, a probe copy or a point outside a data set).  Set when the plan is built (WorkingCopy::BuildPlan).
  int PlanTask() const { return plan_task_; }
  void SetPlanTask(int task) { plan_task_ = task; }
  /// Take a copy of the point out of the plan: it calculates itself alone (a probe).
  void DetachFromPlan() { plan_task_ = -1; }
  void SetTargetEffectNum(int);
  /// Make `sub` a sub-point of this point (WorkingCopy::AddSubPoint places it): it shares the term's settings and the
  /// point's screening.
  void MakeSubPoint(Unit &sub) const;
  /// The point's sub-points, a range of the units of `data` (WorkingCopy::Units); set by the data set.
  int SubPointFirst() const { return sub_first_; }
  void SetSubPointRange(int first, int count) { sub_first_ = first, sub_count_ = count; }
  /// The point's fold of its observable slots (FoldPoint, include/PointFold.h, on each slot the data set compares).
  void IntegrateTargetEffectForObservable(const WorkingCopy &data, const Config &, azure2::Results &results);
  /// Integrate the E1 and E2 components of the cross section over the same
  /// sub-points as the total, for a segment compared against one component.
  void IntegrateTargetEffectComponents(const WorkingCopy &data, const Config &, azure2::Results &results);
  /// Per-point energy window (columns 5-6 of the data file) of a beam-profile
  /// effect: the detector-reconstructed energy slice the point was built from.
  bool HasBinWindow() const { return bin_low_cm_ == bin_low_cm_ && bin_high_cm_ == bin_high_cm_; };
  /// Set the per-point energy window directly, in the c.m. frame.  A row of a measurement takes it from columns 5-6
  /// of the data file through its frame; this is for a point built outside that path, as the beam-profile
  /// reference test does.
  void SetBinWindow(double low, double high) {
    bin_low_cm_ = low;
    bin_high_cm_ = high;
  };
  double GetBinLowCM() const { return bin_low_cm_; };
  double GetBinHighCM() const { return bin_high_cm_; };
  /// Kinematics of the inverse photodissociation reaction for the
  /// detailed-balance weight: Q_gamma = S_entrance - Ex_final (MeV) and the
  /// compound-nucleus mass (MeV).
  void SetPhotoKinematics(double qGamma, double compoundMassMeV) { photo_qgamma_ = qGamma; photo_mcn_ = compoundMassMeV; };
  double GetPhotoQ() const { return photo_qgamma_; }
  double GetPhotoCompoundMass() const { return photo_mcn_; }
  void SetStoppingPower(double);
  void SetTargetThickness(double);
  void SetAngularDists(azure2::Results &r, vector_r dists);
  /// Sub-point i (1-based) in the units of `data`, the data set that holds the point; null beyond the last.
  Unit *GetSubPoint(WorkingCopy &data, int i);
  const Unit *GetSubPoint(const WorkingCopy &data, int i) const;
  UnitRange<Unit> GetSubPoints(WorkingCopy &data);
 private:
  bool is_sub_point_ = false;
  double cm_factor_ = 1.0;   ///< c.m./lab energy factor of the point's entrance pair (SetTargetKinematics)
  int bohr_z1_ = 0;          ///< beam charge of the entrance pair, for Bohr straggling
  int bohr_z2_ = 0;          ///< target charge of the entrance pair
  int set_up_ = -1;  ///< SetUpIndex
  int segment_key_;
  int targetEffectNum_;
  double target_blend_weight_ = 1.0;
  bool target_fallback_only_ = false;   ///< see SetTargetFallbackOnly
  bool target_fallback_now_ = false;    ///< see TargetFallbackScope
  double target_fallback_depth_ = 0.0;  ///< D = E - E_bar, c.m. MeV (see SetTargetFallbackOnly)
  double cm_angle_;
  double lab_angle_;       ///< the lab angle the point is calculated at: its row's, a fixed angle, or a parent's
  double nominal_energy_;  ///< a point that is no row: its energy before any shift (MeV)
  /// What the point's term computes (shared, immutable once set up; TermSettings.h).
  /// Not owning: the term (segment or leaf) holds the record by a shared pointer that every copy of the data set
  /// copies, so a copy of a unit touches no reference count.  A point outside a data set owns its own (own_settings_).
  const TermSettings *settings_ = nullptr;
  std::shared_ptr<const TermSettings> own_settings_;
  int row_ = -1;            ///< the point's row in its measurement (0-based), -1 for a point that is no row
  int legendre_order_ = -1; ///< the Legendre order k of an angular-coefficient row (and its sub-points), else -1
  double cm_energy_;
  double lab_energy_;
  double excitation_energy_;
  double geofactor_;
  double sfactorconv_;
  double screen_k_ = 0.0;
  double screen_ue_ = 0.0;
  int screen_param_ = 0;
  double stoppingPower_;
  double targetThickness_;
  double targetLabEnergy_ = 0.0;
  // Beam-profile effect: per-point energy window (c.m.; NaN when the data file has no such columns) and the
  // inverse reaction kinematics for the detailed-balance weight.
  double bin_low_cm_ = std::numeric_limits<double>::quiet_NaN();
  double bin_high_cm_ = std::numeric_limits<double>::quiet_NaN();
  double photo_qgamma_ = 0.0;
  double photo_mcn_ = 0.0;
  complex coulombamplitude_;
  // Direct and exchange Coulomb amplitudes, kept apart only for an identical
  // pair with spin, whose exchange sign depends on the channel spin.
  bool spinDependentCoulomb_ = false;
  complex coulombdirect_ = complex(0., 0.);
  complex coulombexchange_ = complex(0., 0.);
  // The Legendre polynomials P_l(cos theta) of the point, by order, and the energy-dependent channel tables (one row
  // per J group: L_o, sqrt(P), Coulomb and hard-sphere phase factors).  Filled at set-up and rewritten only when the
  // point's energy or grid changes, so the working copies share them (a copy of a folded project with ~10^6
  // sub-points held ~7 GB of these rows); every writer goes through MutableEDep() / MutableLegendre(), which give
  // the point its own copy first if the tables are shared (copy-on-write).
  SharedTable<vector_r> legendreP_;
  SharedTable<EDependentTables> edep_;
  /// This point's own, writable tables (copied first if shared with another point).
  EDependentTables &MutableEDep();
  /// This point's own, writable Legendre polynomials (copied first if shared).
  vector_r &MutableLegendre();
  std::shared_ptr<const ThmPsTable> thm_ps_;
  // Energy at which CalcEDependentValues last ran: RecalcEDependentValues is
  // a no-op while the point's energy is unchanged (its values depend only on
  // the energy, the angle and the channel radii).
  double eDependentEnergy_ = 0.0;
  bool eDependentValid_ = false;
  matrix_c ec_amplitudes_;
  matrix_r ec_energies_;  // Energies at which EC amplitudes were calculated
  /// The amplitude at the last shifted energy of each pathway (GetECAmplitudeWithShift).
  struct ShiftedECAmplitude {
    double energy = 0.0;
    complex amplitude = complex(0.0, 0.0);
    bool valid = false;
  };
  mutable std::vector<std::vector<ShiftedECAmplitude>> ec_shifted_;
  int plan_task_ = -1;  // see PlanTask
  int sub_first_ = 0, sub_count_ = 0;  ///< the sub-points: units [first, first + count) of the data set
  // Where the result record lies (see ResultSlot): no value of an evaluation is held in the point itself.
  std::size_t slot_ = kNoRecord;
  std::uint16_t num_ang_slots_ = 0;
  std::uint16_t num_natural_slots_ = 0;
  std::uint16_t num_model_ang_ = 0;  // see SetModelRecord
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
  // Where the model record lies (see SetModelRecord).
  std::size_t model_slot_ = kNoRecord;
  std::size_t ModelRecord(const azure2::Results &r) const {
    if (model_slot_ == kNoRecord || model_slot_ + azure2::Results::ModelRecordSize(num_model_ang_) > r.model.size())
      NoModelRecord(r);
    return model_slot_;
  }
  [[noreturn]] void NoModelRecord(const azure2::Results &r) const;
};

#endif
