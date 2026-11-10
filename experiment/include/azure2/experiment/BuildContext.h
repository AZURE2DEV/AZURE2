// Internal to Schema::Build (experiment/src/Build.cpp, BuildEffects.cpp): the state of one build and the helpers both
// files use.  Not for other callers: Schema::Build (Schema.h) is the interface.
#ifndef AZURE2_EXPERIMENT_BUILDCONTEXT_H
#define AZURE2_EXPERIMENT_BUILDCONTEXT_H

#include <map>
#include <string>
#include <vector>

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
  Context(const ValidationResult &v, const std::vector<PairRow> &pairs, const BuildOptions &options, BuiltProject &out)
      : v_(v), pairs_(pairs), opt_(options), out_(out) {}

  const ValidationResult &v_;
  const std::vector<PairRow> &pairs_;
  const BuildOptions &opt_;
  BuiltProject &out_;
  std::vector<Member> data_, extrapolations_;

  /// An error of `record`, at `field` (its first line when given, else the record's header).
  void Error(const char *code, const ValidatedRecord &record, const std::string &field, const std::string &message,
             const std::vector<std::string> &allowed = {});
  const PairRow *Pair(int key) const;
  /// The record named `name` of kind `kind`, or nullptr.
  const ValidatedRecord *Record(const std::string &name, const std::string &kind) const;
  /// The position of a record in the project (records in input order), for ordering.
  int Position(const ValidatedRecord *r) const;
  /// "MeV": the value in the engine's unit of its dimension (the unit factor); bit-exact for factor 1.
  static double InBase(const Number &n, const Unit *unit) { return unit && unit->factor != 1.0 ? n.value * unit->factor : n.value; }
};

/// The effects of targets, beams and detectors on the data sets and extrapolations, the abundances of the targets'
/// isotopes and the sample thicknesses (BuildEffects.cpp).
void BuildEffects(Context &c);

}  // namespace build
}  // namespace experiment
}  // namespace azure2

#endif
