// Internal to Schema::Build (experiment/src/Build.cpp, BuildEffects.cpp): the state of one build and the helpers both
// files use.  Not for other callers: Schema::Build (Schema.h) is the interface.
#ifndef AZURE2_EXPERIMENT_BUILDCONTEXT_H
#define AZURE2_EXPERIMENT_BUILDCONTEXT_H

#include <map>
#include <string>
#include <vector>

#include "azure2/experiment/ProjectPaths.h"
#include "azure2/experiment/Schema.h"

namespace azure2 {
namespace experiment {
namespace build {

/// A data set or an extrapolation as Build resolved it: its record, its key and the pair keys of its first term.
struct Member {
  const ValidatedRecord *record = nullptr;
  int key = 0;
  bool active = true;
  bool isData = true;
  int entranceKey = 0;               ///< the pair key of the first term's entrance
  std::vector<int> termEntranceKeys; ///< every term's entrance pair key, in order
};

class Context {
 public:
  Context(const ValidationResult &v, const std::vector<PairInfo> &pairs, const BuildOptions &options, BuiltProject &out)
      : v_(v), pairs_(pairs), opt_(options), out_(out) {}

  const ValidationResult &v_;
  const std::vector<PairInfo> &pairs_;
  const BuildOptions &opt_;
  BuiltProject &out_;
  std::vector<Member> data_, extrapolations_;

  /// An error of `record`, at `field` (its first line when given, else the record's header).
  void Error(const char *code, const ValidatedRecord &record, const std::string &field, const std::string &message,
             const std::vector<std::string> &allowed = {});
  /// A warning of `record` at `field` (as Error, severity Warning).
  void Warn(const char *code, const ValidatedRecord &record, const std::string &field, const std::string &message,
            const std::string &found = "");
  const PairInfo *Pair(int key) const;
  /// A file's path as the engine opens it: relative to the project directory ("" = the working directory;
  /// ProjectPaths.h).  For opening only: messages and specs keep the name the project writes.
  std::string ProjectPath(const std::string &file) const { return PathToOpen(opt_.projectDir, file); }
  /// The record named `name` of kind `kind`, or nullptr.
  const ValidatedRecord *Record(const std::string &name, const std::string &kind) const;
  /// The position of a record in the project (records in input order), for ordering.
  int Position(const ValidatedRecord *r) const;
  /// The value in the engine's unit of its dimension (Units.h ToBase: the unit table's factor; bit-exact for factor 1).
  static double InBase(const Number &n, const Unit *unit) { return ToBase(n.value, unit); }
};

/// The centre of a nuisance's prior as the specs carry it (PriorCentre, docs/upgrade/parameter_kinds.md section 4):
/// written when the nuisance has a prior; `centre` is its value in the engine's unit of the kind, `what` names the
/// field in messages.
inline PriorCentre PriorCentreOf(const Nuisance &n, double centre, const std::string &what) {
  PriorCentre c;
  if (!n.hasPrior) return c;
  c.written = true;
  c.value = centre;
  c.source = what + " prior " + n.priorCentre.text;
  return c;
}

/// The pair key `labelField` names by its label (light, heavy, state; `text` as written, `what` for the messages:
/// "the entrance of 15N(p,g)16O", "the pair p+12C"), with the pair-key field `keyField` where the label names several
/// pairs (section 3.5, D2/D3): accepted only there and required there; "" where the record cannot carry a key.
/// Returns false after an error (E_UNKNOWN_PAIR, E_AMBIGUOUS_PAIR, E_PAIR_KEY).
bool ResolvePair(Context &c, const ValidatedRecord &r, const std::string &labelField, const std::string &what,
                 const Particle &light, const Particle &heavy, const ResidualState &state, const std::string &text,
                 const std::string &keyField, int &key);

/// The unit of a data set's or an extrapolation's values (BuildUnits.cpp): for a data set (`columns` given) the factor
/// of its value column; for an extrapolation a constant on the model; for a yield the factor per unit charge and the
/// efficiency (ReactionSpec::scale), and fold windows refused on a layer.  Returns false after an error.
bool ValueUnit(Context &c, const ValidatedRecord &r, ReactionSpec &spec, ColumnConversion *columns);
/// The column layout of a data set's file and the conversion of its rows (DatasetSpec::columns): units of the energy
/// and angle roles, an E_cm column taken to E_lab with the entrance pair's masses.  Returns false after an error.
bool DataColumns(Context &c, const ValidatedRecord &r, const Columns &cols, const std::string &secondRole, DatasetSpec &d);
/// The data set's energy window in lab MeV (the frame and unit of its energy column), after DataColumns.
void EnergyWindow(const ValidatedRecord &r, DatasetSpec &d);

/// The `potential` records of <potentials> (BuildPotentials.cpp): one setting per pair into out_.potentials, and the
/// refusals that need the pairs, the THM data sets and the relativistic flag.  After the data sets are built.
void BuildPotentials(Context &c);

/// The effects of targets, beams and detectors on the data sets and extrapolations, the abundances of the targets'
/// isotopes and the sample thicknesses (BuildEffects.cpp).
void BuildEffects(Context &c);

}  // namespace build
}  // namespace experiment
}  // namespace azure2

#endif
