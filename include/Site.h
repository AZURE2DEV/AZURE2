#ifndef AZURE2_SITE_H
#define AZURE2_SITE_H
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
  void BelongToRowOf(const Site &parent) {
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
  /// Fill the compound excitation energy from the c.m. energy (a sub-point moved by an energy shift).
  void ConvertExcitationEnergy(PPair *);
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
  void CalcLegendreP(WorkingCopy *data, int, CNuc *, const TargetEffect *);
  /// Compute the geometrical factor, S-factor conversion, \f$L_o\f$ elements, penetrabilities and phases together.
  /// The rows of the open and closed particle channels are read from the tables of `rows` (docs/upgrade/
  /// channel_table.md) when one holds them, and otherwise computed (ChannelTable::Compute), once per (pair, l, E_c)
  /// of the point; rows->used, when set, receives every row used.
  void CalcEDependentValues(WorkingCopy *data, CNuc *, const Config &, const ChannelRowSource *rows = nullptr);
  /// The keys (pair, l, E_c) of the channel-table rows CalcEDependentValues reads for this site.
  void ChannelKeys(WorkingCopy *data, CNuc *, std::vector<ChannelTable::Key> &);
  /// The compound-nucleus energy E_x (MeV) at which CalcEDependentValues evaluates the channels.
  double CompoundEnergy(CNuc *);
  /// Recompute those at the current, possibly shifted, energy; nothing while they are at that energy.
  void RecalcEDependentValues(WorkingCopy *data, CNuc *, const Config &, const ChannelRowSource *rows = nullptr);
  /// Are the energy-dependent values those of the site's present energy?
  bool EDependentCurrent() const { return eDependentValid_ && this->GetCMEnergy() == eDependentEnergy_; }
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
  void CalcCoulombAmplitude(WorkingCopy *data, CNuc *);
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
  void SetAngularDists(azure2::Results &r, vector_r dists);
 protected:
  bool is_sub_point_ = false;
  int segment_key_;
  double cm_angle_;
  double lab_angle_;       ///< the lab angle the point is calculated at: its row's, a fixed angle, or a parent's
  double nominal_energy_;  ///< a point that is no row: its energy before any shift (MeV)
  /// What the point's term computes (shared, immutable once set up; TermSettings.h).
  /// Not owning: the term (segment or leaf) holds the record by a shared pointer that every copy of the data set
  /// copies, so a copy of a unit touches no reference count.  A point outside a data set owns its own (own_settings_).
  const TermSettings *settings_ = nullptr;
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
  // Where the result record lies (see ResultSlot): no value of an evaluation is held in the point itself.
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