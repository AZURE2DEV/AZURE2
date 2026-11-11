#ifndef AZURE2_EXPERIMENT_DIAGNOSTICS_H
#define AZURE2_EXPERIMENT_DIAGNOSTICS_H

#include <string>
#include <vector>

#include "azure2/io/Record.h"

namespace azure2 {
namespace experiment {

/*!
 * Diagnostics of a format-3 project (docs/upgrade/project_format3.md section 6.2): every problem the reader and the
 * validator find, as data -- a stable code, where (line, column, section, record, field), what was found, what is
 * allowed, a message that names the record and field, and optionally a fix that changes syntax only (a misspelled key,
 * a unit's spelling, a missing field as a '?' placeholder; never a number, decision F18).
 *
 * The codes are listed with their severity and meaning by DiagnosticCodes(); the schema exports that list
 * (Schema::ToJson), and every code has a test that triggers it alone (repo_consistency check 24).
 */
enum class Severity { Error, Warning, Ask };

struct Fix {
  std::string action;  ///< "replace_key", "replace_token" (at line, col) or "insert_field" (after afterLine)
  int line = 0;
  int col = 0;
  int afterLine = 0;
  std::string text;  ///< the new key or token, or the line to insert
};

struct Diagnostic {
  Severity severity = Severity::Error;
  std::string code;
  std::string file;  ///< set by the caller that knows it (the validator works on records)
  int line = 0;
  int column = 0;
  std::string section;
  std::string recordKind;
  std::string record;
  std::string field;
  bool hasFound = false;
  std::string found;
  std::vector<std::string> allowed;
  std::string message;
  bool hasFix = false;
  Fix fix;

  /// One JSON object, keys as in section 6.2 (absent members omitted, "found": null when nothing was found).
  std::string ToJson() const;
  /// "file:line: record.field: message" (section 6.2's text form; parts that are empty are left out).
  std::string ToText() const;
};

struct DiagnosticCode {
  const char *code;
  Severity severity;
  const char *meaning;
};

/// Every code the reader (io::record_codes) and the validator emit.
const std::vector<DiagnosticCode> &DiagnosticCodes();

/// The reader's diagnostic as a project diagnostic of `section`.
Diagnostic FromRecordDiagnostic(const io::RecordDiagnostic &d, const std::string &section);

const char *SeverityName(Severity s);

/// The validator's codes (the reader's are io::record_codes).
namespace codes {
// Records and names
constexpr const char *kWrongSection = "E_WRONG_SECTION";
constexpr const char *kBadName = "E_BAD_NAME";
constexpr const char *kDuplicateName = "E_DUPLICATE_NAME";
// Fields
constexpr const char *kUnknownField = "E_UNKNOWN_FIELD";
constexpr const char *kFieldTwice = "E_FIELD_TWICE";
constexpr const char *kFieldRequired = "E_FIELD_REQUIRED";
constexpr const char *kFieldNotApplicable = "E_FIELD_NOT_APPLICABLE";
constexpr const char *kValueMissing = "E_VALUE_MISSING";
constexpr const char *kPlaceholder = "E_PLACEHOLDER";
// Values
constexpr const char *kTooFewTokens = "E_TOO_FEW_TOKENS";
constexpr const char *kTooManyTokens = "E_TOO_MANY_TOKENS";
constexpr const char *kNotANumber = "E_NOT_A_NUMBER";
constexpr const char *kNotAnInteger = "E_NOT_AN_INTEGER";
constexpr const char *kOutOfRange = "E_OUT_OF_RANGE";
constexpr const char *kEmptyRange = "E_EMPTY_RANGE";
constexpr const char *kNotInList = "E_NOT_IN_LIST";
constexpr const char *kUnknownUnit = "E_UNKNOWN_UNIT";
constexpr const char *kUnitMissing = "E_UNIT_MISSING";
constexpr const char *kUnitNotAllowed = "E_UNIT_NOT_ALLOWED";
constexpr const char *kNuisanceMode = "E_NUISANCE_MODE";
constexpr const char *kModeNotAllowed = "E_MODE_NOT_ALLOWED";
constexpr const char *kPriorSyntax = "E_PRIOR_SYNTAX";
constexpr const char *kPartMissing = "E_PART_MISSING";
constexpr const char *kPartTwice = "E_PART_TWICE";
constexpr const char *kFormulaSyntax = "E_FORMULA_SYNTAX";
constexpr const char *kSpinSyntax = "E_SPIN_SYNTAX";
constexpr const char *kNucleusSyntax = "E_NUCLEUS_SYNTAX";
constexpr const char *kUnknownElement = "E_UNKNOWN_ELEMENT";
constexpr const char *kPairSyntax = "E_PAIR_SYNTAX";
constexpr const char *kWindowSyntax = "E_WINDOW_SYNTAX";
// Reactions
constexpr const char *kReactionSyntax = "E_REACTION_SYNTAX";
constexpr const char *kReactionOrdinal = "E_REACTION_ORDINAL";
constexpr const char *kReactionBalance = "E_REACTION_BALANCE";
constexpr const char *kReactionNotAllowed = "E_REACTION_NOT_ALLOWED";
// Column roles
constexpr const char *kUnknownRole = "E_UNKNOWN_ROLE";
constexpr const char *kRoleUnit = "E_ROLE_UNIT";
constexpr const char *kRoleTwice = "E_ROLE_TWICE";
constexpr const char *kRoleRequired = "E_ROLE_REQUIRED";
constexpr const char *kRoleNotAllowed = "E_ROLE_NOT_ALLOWED";
// Between records
constexpr const char *kUnknownReference = "E_UNKNOWN_REFERENCE";
constexpr const char *kReferenceKind = "E_REFERENCE_KIND";
constexpr const char *kUnknownEpoch = "E_UNKNOWN_EPOCH";
constexpr const char *kTargetKind = "E_TARGET_KIND";
constexpr const char *kDegenerateNuisance = "E_DEGENERATE_NUISANCE";
// Model data against <levels> (ModelDatum.h, BuildModelData)
constexpr const char *kUnknownLevel = "E_UNKNOWN_LEVEL";
constexpr const char *kLevelJpi = "E_LEVEL_JPI";
constexpr const char *kUnknownPair = "E_UNKNOWN_PAIR";
constexpr const char *kAmbiguousPair = "E_AMBIGUOUS_PAIR";
constexpr const char *kNotAChannel = "E_NOT_A_CHANNEL";
// Build checks (Schema::Build: they need <levels> or the data files)
constexpr const char *kPairKey = "E_PAIR_KEY";
constexpr const char *kDataFile = "E_DATA_FILE";
constexpr const char *kDataFileNotRead = "W_DATA_FILE_NOT_READ";
constexpr const char *kDataRow = "E_DATA_ROW";
constexpr const char *kDataError = "W_DATA_ERROR";
constexpr const char *kCovarianceFile = "E_COVARIANCE_FILE";
constexpr const char *kFinalStateNotLevel = "E_FINAL_STATE_NOT_LEVEL";
constexpr const char *kFinalStateShared = "E_FINAL_STATE_SHARED";
constexpr const char *kNotBuilt = "E_NOT_BUILT";
constexpr const char *kThmRelativistic = "E_THM_RELATIVISTIC";
constexpr const char *kThmCovariance = "E_THM_COVARIANCE";
constexpr const char *kFoldConflict = "E_FOLD_CONFLICT";
// Hybrid potentials (`potential` records: BuildPotentials; the old block in a format-3 file: ReadFormat3Project)
constexpr const char *kPotentialTwice = "E_POTENTIAL_TWICE";
constexpr const char *kPotentialPairKind = "E_POTENTIAL_PAIR_KIND";
constexpr const char *kPotentialRelativistic = "E_POTENTIAL_RELATIVISTIC";
constexpr const char *kThmHybridPotential = "E_THM_HYBRID_POTENTIAL";
constexpr const char *kOldPotentialBlock = "E_OLD_POTENTIAL_BLOCK";
// Spin-orbit terms (a potential record's spin_orbit, pyazr, the API: SpinOrbitChecks.h)
constexpr const char *kSoTwoSpins = "E_SO_TWO_SPINS";
constexpr const char *kSoZeroSpin = "E_SO_ZERO_SPIN";
constexpr const char *kSoRelativistic = "E_SO_RELATIVISTIC";
constexpr const char *kSoThm = "E_SO_THM";
constexpr const char *kSoNeedsJ = "E_SO_NEEDS_J";
constexpr const char *kSoNoTerm = "E_SO_NO_TERM";
// Warnings
constexpr const char *kNotMeasurable = "W_NOT_MEASURABLE";
constexpr const char *kBroadLevel = "W_BROAD_LEVEL";
constexpr const char *kNearbyLevel = "W_NEARBY_LEVEL";
}  // namespace codes

}  // namespace experiment
}  // namespace azure2

#endif
