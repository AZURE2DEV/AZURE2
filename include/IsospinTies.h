#ifndef ISOSPINTIES_H
#define ISOSPINTIES_H

#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "Constants.h"
#include "ParameterKind.h"
#include "azure2/experiment/Isospin.h"

class AZUREParams;
class CompoundSystem;
class Config;
class LevelState;
namespace ROOT {
namespace Minuit2 {
class MnUserParameters;
}
}  // namespace ROOT

/*!
 * Isospin-labelled levels across compound nuclei (the records of <isospin>, azure2/experiment/Isospin.h; Hale's
 * charge-independent R-matrix, as in his A = 4 analysis; dev's IsospinLevels of 2026-10-17, re-written here).
 *
 * A level whose <levels> lines span several compound nuclei under one level ID is one isospin level of isospin T.  In
 * every nucleus it lives in, its reduced-width amplitudes are the Clebsch-Gordan projection of one amplitude per
 * isospin channel (arrangement, l, s), and its energy above the nucleus's threshold pair is the reference member's
 * above the reference's, plus the difference of the nuclei's internal Coulomb shifts.  The member in the reference
 * nucleus carries the free parameters; the other members' are derived from them.  An isospin-forbidden channel
 * (factor 0) is zero unless an isospin_mix record gives it an amplitude of its own.  A soft level's derived slots are
 * free instead, each with a penalty ((value - tie)/sigma)^2.
 *
 * The ties act on the fitted parameters of the run's basis:
 *  - formal (--no-brune; Hale's formulation): formal amplitudes and eigenvalues, and every J^pi group that holds an
 *    isospin level has B_c = -l_c (CompoundSystem::CalcBoundaryConditions), so that partners share their boundary values;
 *  - Brune (the default): Brune amplitudes and observed energies; no boundary value enters.
 * Park's basis and the RMC formalism are refused.
 *
 * A set-up fact of the compound system: built once by CompoundSystem::Fill (Read), constant afterwards, shared by every
 * evaluation.  Everything is a function of the parameter vector it is given: Apply writes the derived slots of a
 * vector from its free ones, and CompoundSystem::FillCompoundFromParams calls it first, so every fill sees consistent
 * members; the penalties and the report read the vector.  The isospin parameters (the shifts, then the mixes) follow
 * the level block in the vector (CompoundSystem::FillMnParams).
 */

/// A derived slot of the parameter vector as a linear function of free slots: d(slot)/d(source) = coefficient.
struct IsospinTie {
  int slot = -1;
  std::vector<std::pair<int, double>> sources;
};

class IsospinSet {
 public:
  /// The records of <isospin> resolved against the compound system after its levels, channels and parameter layout
  /// are final: members, factors, owners, slots.  Reports what it built.  nullptr after a refusal (reported); an
  /// inactive set without records.
  static std::shared_ptr<const IsospinSet> Read(const CompoundSystem &compound, const Config &configure);

  bool Active() const { return active_; }
  /// The J^pi groups (1-based) that hold a member of an isospin level: B_c = -l_c there in the formal basis.
  const std::set<int> &Groups() const { return groups_; }
  /// Number of isospin parameters (the shifts, then the mixes) and the index of the first in the vector.
  int NumParams() const { return (int)shifts_.size() + (int)mixes_.size(); }
  int FirstParam() const { return firstParam_; }
  /// The record of isospin parameter k (0-based: the shifts, then the mixes).
  ParameterRecord RecordOf(int k) const;
  /// CompoundSystem::FillMnParams: the isospin parameters after the level block, the derived slots written and fixed, the
  /// soft slots set on their targets and freed.
  void AddParameters(AZUREParams &params) const;
  /// The values of the isospin parameters as the records write them (the shifts, then the mixes): a level state's
  /// until a fill from a vector sets others (LevelState::IsospinValues).
  std::vector<double> FileValues() const;
  /// The values of the isospin parameters in a parameter vector that holds them; empty for a shorter vector.
  std::vector<double> ValuesIn(const vector_r &q) const;
  /// Overwrites the derived slots of a parameter vector from the free ones: the level slots of q and the isospin
  /// parameters `values` (the shifts, then the mixes).  The one-argument form takes them from q itself.
  void Apply(vector_r &q, const std::vector<double> &values) const;
  void Apply(vector_r &q) const { Apply(q, ValuesIn(q)); }
  /// Writes the derived values into the Minuit parameters and keeps those slots fixed (after a parameter file).
  void Sync(ROOT::Minuit2::MnUserParameters &p) const;
  /// Every derived slot as a linear function of the free ones, for the chain rule of the gradient.
  const std::vector<IsospinTie> &Ties() const { return ties_; }
  bool IsDerivedSlot(int slot) const { return derived_.count(slot) > 0; }
  bool HasSoftTies() const { return !softTerms_.empty(); }

  /// The priors of the fitted shifts and mixes, then the soft ties, as residuals at the vector q: r and dr/dq as
  /// (index, value) pairs.  chi^2 = sum r^2.
  struct PenaltyRow {
    double r = 0.0;
    std::vector<std::pair<int, double>> d;
  };
  /// At the levels and isospin values of `state` (what the last fill left: nothing here reads a vector).
  std::vector<PenaltyRow> PenaltyRows(const CompoundSystem &compound, const LevelState &state) const;
  bool HasPenalty() const;
  double Penalty(const CompoundSystem &compound, const LevelState &state) const;
  /// The soft ties' part of Penalty alone.
  double SoftPenalty(const CompoundSystem &compound, const LevelState &state) const;
  /// d chi^2 / d q as (index, value) pairs (GradAccum::direct).
  void AddPenaltyGradient(const CompoundSystem &compound, const LevelState &state, std::vector<std::pair<int, double>> &grad) const;
  /// f(index, prior centre, prior width) for every fitted shift and mix with a prior.
  template <class F>
  void ForEachPrior(F f) const {
    for (const Shift &e : shifts_)
      if (e.spec.shift.fit && e.spec.shift.hasPrior) f(e.param, e.spec.shift.priorCentre, e.spec.shift.priorWidth);
    for (const Mix &m : mixes_)
      if (m.spec.amplitude.fit && m.spec.amplitude.hasPrior) f(m.param, m.spec.amplitude.priorCentre, m.spec.amplitude.priorWidth);
  }

  /// The label of parameter `index` if it is an isospin parameter, else "".
  std::string ParamLabel(int index) const;
  /// "ISOSPIN LEVEL <id> T=<t> (reference|derived|soft)" for a level of a J^pi group, or "".
  std::string MemberLabel(int jGroup, int level) const;
  /// isospin.out: the shifts, mixes and every member's energy and amplitudes at the values of `state`.
  void Write(const Config &configure, const CompoundSystem &compound, const LevelState &state) const;

 private:
  struct Shift {
    azure2::experiment::IsospinNucleusEntry spec;
    int nucleus = 0;         ///< 1-based nucleus of the compound system
    double threshold = 0.0;  ///< E_x of the threshold pair (S + excitations)
    int param = -1;          ///< index of its parameter; -1 for the reference
  };
  struct Mix {
    azure2::experiment::IsospinMixEntry spec;
    int param = -1;  ///< index of its parameter
    int slot = -1;   ///< index of the amplitude it sets
  };
  /// One amplitude of one member: the owner of its isospin channel (free), a multiple of the owner, or forbidden.
  struct Width {
    int channel = 0;     ///< 1-based channel of the member's J^pi group
    int slot = -1;       ///< index of this amplitude
    int ownerSlot = -1;  ///< the owner's index; -1 for the owner itself, -2 when forbidden
    double cg = 0.0, cgComputed = 0.0;
    bool overridden = false;
    double ratio = 1.0;  ///< cg / cg(owner)
    int mix = -1;        ///< index into mixes_, or -1
    std::string key;     ///< the isospin channel (arrangement, l, s), for the report
  };
  struct Member {
    int nucleus = 0, jGroup = 0, level = 0;
    double t3 = 0.0;
    double threshold = 0.0;
    int energySlot = -1;
    std::vector<Width> widths;  ///< particle channels only
    int untied = 0;             ///< photon and beta channels of the group: outside the ties
  };
  struct Level {
    azure2::experiment::IsospinLevelEntry spec;
    double t = 0.0;
    int ref = 0;  ///< index of the reference member
    std::vector<Member> members;
  };
  /// One soft-tie penalty: ((q[slot] - target)/sigma)^2, target = constant + sum coefficient q[source].
  struct SoftTerm {
    int slot = -1;
    std::vector<std::pair<int, double>> sources;
    double constant = 0.0, sigma = 1.0;
  };
  const Shift *ShiftOf(int nucleus) const;
  int ShiftParamOf(int nucleus) const;
  double ShiftValue(int nucleus, const std::vector<double> &values) const;
  /// The value at index `index` of the parameter vector as `state` holds it: a level's energy or amplitude, or an
  /// isospin parameter.
  double At(const CompoundSystem &compound, const LevelState &state, int index) const;

  bool active_ = false;
  bool brune_ = false;
  std::vector<Shift> shifts_;  ///< the nuclei with a parameter (not the reference), in record order
  Shift reference_;
  std::vector<Shift> nuclei_;  ///< the reference and the shifts, in record order (the report)
  std::vector<Mix> mixes_;
  std::vector<Level> levels_;
  std::set<int> derived_, groups_;
  std::vector<IsospinTie> ties_;
  std::vector<SoftTerm> softTerms_;
  int firstParam_ = -1;
};

#endif
