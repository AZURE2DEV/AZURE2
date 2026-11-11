#ifndef AZURE2_DATASET_H
#define AZURE2_DATASET_H

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "CalculationPlan.h"
#include "DataCovariance.h"
#include "Segment.h"
#include "TargetEffect.h"
#include "ThmGroup.h"

class PriorCovariance;

/*!
 * The data of a project as the set-up leaves it (docs/upgrade/data_structure.md section 1.3): the data sets (segments)
 * with their measurements, settings and the terms of their predictions; the target effects (lines, kernels, priors);
 * the covariance groups; the THM experiments; the parameter layout.  Built once by the set-up
 * and then never written: every working copy (WorkingCopy) holds it through a std::shared_ptr<const DataSet>, so a copy of
 * a working copy shares it and an evaluation can only read it.  The set-up writes the object before it is frozen
 * (WorkingCopy::Draft, WorkingCopy::Freeze); a later change of a setting (a prior covariance given by the caller) makes a new
 * DataSet (WorkingCopy::ReviseDataSet), which copies made before do not see.
 *
 * Shared and written after the set-up, on purpose (named here so that nobody adds another):
 *  - a covariance group's frozen model of its relvector terms (DataCovariance::Group, held through a shared pointer
 *    to a non-const group): every copy must use the one model frozen at the first evaluation of a minimization.  It
 *    is frozen by Group::EnsureFactored under the group's mutex (double-checked through an atomic flag, so the
 *    evaluations that run at once on several copies freeze it once and then only read it); it is forgotten by
 *    WorkingCopy::RefreezeCovariance, which the CLI calls between two minimizations, when no evaluation runs.
 */
class DataSet {
 public:
  using LeafRef = std::pair<int, int>;  ///< (0-based position of a segment, its term)

  // ---- the data sets
  int NumSegments() const { return (int)segments_.size(); }
  /// Segment i, 1-based.
  const Segment *GetSegment(int i) const { return &segments_[i - 1]; }
  Segment *GetSegment(int i) { return &segments_[i - 1]; }
  const std::vector<Segment> &GetSegments() const { return segments_; }
  std::vector<Segment> &GetSegments() { return segments_; }
  /// The segment with this key (its data set number in the file), or null.
  const Segment *GetSegmentFromKey(int key) const;
  Segment *GetSegmentFromKey(int key) { return const_cast<Segment *>(static_cast<const DataSet *>(this)->GetSegmentFromKey(key)); }
  bool IsSegmentKey(int key) const { return GetSegmentFromKey(key) != nullptr; }
  /// Every leaf as (0-based position of its segment, its term), in segment order and term order: the order the leaves
  /// are set up, laid out and written to the intEC file in.
  const std::vector<LeafRef> &LeafRefs() const { return leafRefs_; }

  /*!
   * What the set-up worked out for one point (a row of a term) that its evaluations read again: the layout of its
   * adaptive sub-point grid (WorkingCopy::ReplaceSubPointGrids lays it out again for other levels) and its user-defined
   * resolution kernel with the projectile rest energy (MeV) its time-of-flight conversion uses.  A unit names its
   * record by Unit::SetUpIndex; the record is shared by every copy.
   */
  struct PointSetUp {
    Unit::SubGridSpec grid;
    std::shared_ptr<const UdrKernel> udr;
    double udrMass = 0.0;
  };
  /// The record of point `p`, null when the set-up gave it none.
  const PointSetUp *SetUpOf(const Unit &p) const {
    return p.SetUpIndex() >= 0 ? &pointSetUps_[p.SetUpIndex()] : nullptr;
  }
  /// The layout of the adaptive grid of point `p` (an invalid one without).
  const Unit::SubGridSpec &GridSpecOf(const Unit &p) const {
    static const Unit::SubGridSpec none;
    return p.SetUpIndex() >= 0 ? pointSetUps_[p.SetUpIndex()].grid : none;
  }
  /// The set-up: the record of point `p`, made on first use.
  PointSetUp &SetUpFor(Unit &p) {
    if (p.SetUpIndex() < 0) {
      p.SetSetUpIndex((int)pointSetUps_.size());
      pointSetUps_.emplace_back();
    }
    return pointSetUps_[p.SetUpIndex()];
  }

  // ---- the target effects (1-based)
  int NumTargetEffects() const { return (int)targetEffects_.size(); }
  const TargetEffect *GetTargetEffect(int n) const {
    return (n >= 1 && n <= (int)targetEffects_.size()) ? &targetEffects_[n - 1] : nullptr;
  }
  TargetEffect *GetTargetEffect(int n) {
    return (n >= 1 && n <= (int)targetEffects_.size()) ? &targetEffects_[n - 1] : nullptr;
  }

  // ---- the parameter layout (WorkingCopy::FillMnParams)
  int GetNormParamOffset() const { return normParamOffset_; }
  int GetEnergyShiftParamOffset() const { return energyShiftParamOffset_; }
  int GetEnergyShiftSqrtParamOffset() const { return energyShiftSqrtParamOffset_; }
  int GetEnergyShiftFnParamOffset() const { return energyShiftFnParamOffset_; }
  /*!
   * Electron screening.  The <targetInt> lines with a `screening` block, in line order, each with one parameter
   * `screening_<n>` (n = the effect number: its place among the active <targetInt> lines; format 3:
   * `<target>.screening`), U_e in eV, after the shift-function parameters.  A line listing several segments has one
   * copy of its effect per segment (WorkingCopy::SetUpTargetEffects); the parameter is the line's, shared by every copy.
   */
  int NumScreeningParams() const { return (int)screeningEffects_.size(); }
  int GetScreeningEffect(int k) const { return screeningEffects_[k - 1]; }
  /// The target effects (copies) of screening parameter k's line.
  const std::vector<int> &ScreeningCopies(int k) const { return screeningCopies_[k - 1]; }
  /// The first copy of screening parameter k's line (its prior is the line's).
  const TargetEffect *GetScreeningTargetEffect(int k) const { return GetTargetEffect(screeningCopies_[k - 1].front()); }
  TargetEffect *GetScreeningTargetEffect(int k) { return GetTargetEffect(screeningCopies_[k - 1].front()); }
  /*!
   * Fitted target densities (TargetEffect::HasDensityFit).  One parameter `density_<n>` (format 3:
   * `<target>.areal_density`; atoms/cm^2) per ungrouped `densityfit` line and per `densitygroup`, n = the effect
   * number of the line (of the group's first line), after the screening parameters.
   */
  int NumDensityParams() const { return (int)densityEffects_.size(); }
  int GetDensityEffect(int k) const { return densityEffects_[k - 1]; }
  int GetDensityLeadCopy(int k) const { return densityLeadCopies_[k - 1]; }
  /// The density parameter (1-based) of target effect (copy) n (1-based), 0 = none.
  int DensityParamOfEffect(int n) const {
    return (n >= 1 && n <= (int)densityParamOfEffect_.size()) ? densityParamOfEffect_[n - 1] : 0;
  }
  const TargetEffect *GetDensityTargetEffect(int k) const { return GetTargetEffect(densityLeadCopies_[k - 1]); }
  TargetEffect *GetDensityTargetEffect(int k) { return GetTargetEffect(densityLeadCopies_[k - 1]); }
  /*!
   * Fitted Doppler temperatures (TargetEffect::HasDopplerFit, `doppler T vary sigma [Tmax]`).  One parameter
   * `doppler_<n>` (format 3: `<target>.doppler`; K) per such line, after the density parameters.
   */
  int NumDopplerParams() const { return (int)dopplerEffects_.size(); }
  int GetDopplerEffect(int k) const { return dopplerEffects_[k - 1]; }
  const std::vector<int> &DopplerCopies(int k) const { return dopplerCopies_[k - 1]; }
  const TargetEffect *GetDopplerTargetEffect(int k) const { return GetTargetEffect(dopplerCopies_[k - 1].front()); }
  TargetEffect *GetDopplerTargetEffect(int k) { return GetTargetEffect(dopplerCopies_[k - 1].front()); }
  /// The Doppler parameter (1-based) of effect \p n, 0 if its temperature is not a fit parameter.
  int DopplerParamOfEffect(int n) const;
  /*!
   * Fitted sample thicknesses (Segment::HasThicknessFit).  One parameter (atoms/barn) per data segment whose line has
   * the fit form, in segment order, after the Doppler temperatures; a thickness-group follower's is fixed at its
   * leader's value.
   */
  int NumThicknessParams() const { return (int)thicknessSegments_.size(); }
  /// Position (1-based, as GetSegment) of the segment owning thickness parameter k (1-based).
  int GetThicknessSegment(int k) const { return thicknessSegments_[k - 1]; }
  /// The thickness parameter (1-based) whose value segment \p pos (1-based) uses; 0 when that thickness is no
  /// parameter.
  int ThicknessParamOfSegment(int pos) const;
  /// The segment (1-based) whose norm segment s uses: its norm-group leader, or s itself.
  int NormOwner(int s) const;
  /*!
   * Stopping scales (TargetEffect::HasStoppingScale; energy_grid.md section 6.2, decision D6).  One parameter
   * `stopping_scale.<material>.<projectile>` (dimensionless) per distinct TargetEffect::GetStoppingScaleKey of the
   * effects, in the order of their first effect, after the THM coherent backgrounds (the last kind).  The parameter
   * holds the first effect's settings (Schema::Build refuses targets that share it and write it differently).
   */
  int NumStoppingScaleParams() const { return (int)stoppingScaleLeads_.size(); }
  /// The first target effect (copy) of stopping-scale parameter k (1-based).
  int GetStoppingScaleLeadCopy(int k) const { return stoppingScaleLeads_[k - 1]; }
  const TargetEffect *GetStoppingScaleTargetEffect(int k) const { return GetTargetEffect(stoppingScaleLeads_[k - 1]); }
  TargetEffect *GetStoppingScaleTargetEffect(int k) { return GetTargetEffect(stoppingScaleLeads_[k - 1]); }
  /// The stopping-scale parameter (1-based) of target effect (copy) n (1-based), 0 = none.
  int StoppingScaleParamOfEffect(int n) const {
    return (n >= 1 && n <= (int)stoppingScaleParamOfEffect_.size()) ? stoppingScaleParamOfEffect_[n - 1] : 0;
  }
  int StoppingScaleParamOffset() const { return stoppingScaleParamOffset_; }
  int ScreeningParamOffset() const { return screeningParamOffset_; }
  int DensityParamOffset() const { return densityParamOffset_; }
  int DopplerParamOffset() const { return dopplerParamOffset_; }
  int ThicknessParamOffset() const { return thicknessParamOffset_; }

  // ---- covariance (DataCovariance.h); the groups' frozen model is shared and mutable (see above)
  bool HasCovariance() const { return !covGroups_.empty(); }
  int NumCovarianceGroups() const { return (int)covGroups_.size(); }
  DataCovariance::Group &GetCovarianceGroup(int g) const { return *covGroups_[g]; }
  /// True when any group has relvector terms (a model frozen at the first evaluation).
  bool HasRelativeCovariance() const;

  // ---- THM experiments (ThmGroup.h)
  int NumThmGroups() const { return (int)thmGroups_.size(); }
  const ThmGroup &GetThmGroup(int g) const { return thmGroups_[g]; }
  int NumThmCoherentParams() const { return (int)thmCoherentParams_.size(); }
  /// Coherent-background parameter k (0-based); its value is the start value (the working copy holds the present one).
  const ThmCoherentParam &GetThmCoherentParam(int k) const { return thmCoherentParams_[k]; }
  int ThmCoherentParamOffset() const { return thmCoherentParamOffset_; }
  /// The non-trivial THM experiment that profiles segment i (1-based) jointly with others or with a background, or -1.
  int ThmGroupOf(int i) const;
  /// Is i the last profiled segment of group g?
  bool IsLastOfThmGroup(int g, int i) const;

  // ---- the run
  /// The prior parameter covariance (<config> priorCovariance=, PriorCovariance.h); null without one.
  const PriorCovariance *GetPriorCovariance() const { return priorCov_.get(); }
  /// The labels of the compound nuclei (the sample decomposition files) and their number.
  const std::vector<std::string> &NucleusLabels() const { return nucleusLabels_; }
  int NumNuclei() const { return numNuclei_; }

 private:
  friend class WorkingCopy;  // the set-up writes the draft (WorkingCopy::Draft)
  std::vector<Segment> segments_;
  std::vector<TargetEffect> targetEffects_;
  std::vector<LeafRef> leafRefs_;
  std::vector<PointSetUp> pointSetUps_;  // SetUpOf
  /// The key of the next leaf (WorkingCopy::AddLeaf): -1000, -1001, ... in creation order, one counter per data set, so a
  /// leaf's key depends on its project alone and no two leaves share one (data_model.md section 1.3, H10).
  int nextComponentKey_ = -1000;
  std::vector<std::shared_ptr<DataCovariance::Group>> covGroups_;
  int normParamOffset_ = 0;
  int energyShiftParamOffset_ = 0;
  int energyShiftSqrtParamOffset_ = 0;
  int energyShiftFnParamOffset_ = 0;
  int screeningParamOffset_ = -1;  // -1 until FillMnParams has laid the vector out
  std::vector<int> screeningEffects_;  // line numbers (effect numbers) with a screening block, in order
  std::vector<std::vector<int>> screeningCopies_;  // per screening parameter: the target effects of its line
  int densityParamOffset_ = -1;
  std::vector<int> densityEffects_;        // line number (effect number) holding each density parameter, in order
  std::vector<int> densityLeadCopies_;     // per density parameter: the first target effect (copy) of that line
  std::vector<int> densityParamOfEffect_;  // per target effect (copy, 0-based): its density parameter, 0 = none
  int stoppingScaleParamOffset_ = -1;
  std::vector<int> stoppingScaleLeads_;          // per stopping-scale parameter: its first target effect (copy)
  std::vector<int> stoppingScaleParamOfEffect_;  // per target effect (copy, 0-based): its scale parameter, 0 = none
  int dopplerParamOffset_ = -1;
  std::vector<int> dopplerEffects_;               // line number holding each Doppler parameter, in order
  std::vector<std::vector<int>> dopplerCopies_;   // per Doppler parameter: the target effects (copies) of its line
  int thicknessParamOffset_ = -1;
  std::vector<int> thicknessSegments_;  // positions (1-based) of the segments with a thickness parameter, in order
  std::vector<ThmGroup> thmGroups_;
  std::vector<ThmCoherentParam> thmCoherentParams_;
  int thmCoherentParamOffset_ = -1;
  std::shared_ptr<const PriorCovariance> priorCov_;
  std::vector<std::string> nucleusLabels_;
  int numNuclei_ = 0;
};

#endif
