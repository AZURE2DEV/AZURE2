// Diagnostics of format 3 (Diagnostics.h; docs/upgrade/project_format3.md section 6.2): the list of codes and the
// JSON and text forms of one diagnostic.
#include "azure2/experiment/Diagnostics.h"

#include "azure2/experiment/Json.h"

namespace azure2 {
namespace experiment {

const char *SeverityName(Severity s) {
  switch (s) {
    case Severity::Error: return "error";
    case Severity::Warning: return "warning";
    case Severity::Ask: return "ask";
  }
  return "?";
}

const std::vector<DiagnosticCode> &DiagnosticCodes() {
  using namespace codes;
  namespace rc = io::record_codes;
  const Severity E = Severity::Error;
  // {code, severity, meaning}
  static const std::vector<DiagnosticCode> table = {
      // The reader (io/RecordReader)
      {rc::kUnterminatedQuote, E, "a quoted token has no closing quote on its line"},
      {rc::kQuoteInToken, E, "a quote inside a token, or text right after a closing quote"},
      {rc::kUnknownKind, E, "a line that should start a record names no record kind of its section"},
      {rc::kHeaderName, E, "a header without a name, or with more than '<kind> <name>' on its line"},
      {rc::kKeySyntax, E, "a field key that is not a letter followed by letters, digits and '_'"},
      // Records and names
      {kWrongSection, E, "a record in a section that does not hold its kind"},
      {kBadName, E, "a name that is not [A-Za-z][A-Za-z0-9_-]*"},
      {kDuplicateName, E, "a name used by two records (one namespace, F3)"},
      // Fields
      {kUnknownField, E, "a key that is no field of the record's kind (fix: the nearest key)"},
      {kFieldTwice, E, "a field that is not repeatable, given twice"},
      {kFieldRequired, E, "a required or conditionally required field is missing (fix: the field with '?')"},
      {kFieldNotApplicable, E, "a field that does not apply to the record's kind, quantity or other fields"},
      {kValueMissing, E, "a key without a value"},
      {kPlaceholder, E, "a value that is '?': a value is needed (F18)"},
      // Values
      {kTooFewTokens, E, "the value stops before its grammar is complete"},
      {kTooManyTokens, E, "tokens after a complete value"},
      {kNotANumber, E, "a token that must be a number"},
      {kNotAnInteger, E, "a token that must be an integer"},
      {kOutOfRange, E, "a number or integer outside its allowed range"},
      {kEmptyRange, E, "a range or window whose low end is above its high end"},
      {kNotInList, E, "a word that is not one of the closed list allowed here"},
      {kUnknownUnit, E, "a token that is no unit (fix: the unit whose spelling differs only in case)"},
      {kUnitMissing, E, "a dimensioned value without its unit (there is no default unit)"},
      {kUnitNotAllowed, E, "a unit of another dimension, or one the quantity does not allow"},
      {kNuisanceMode, E, "a nuisance without its mode fixed or fit"},
      {kModeNotAllowed, E, "a nuisance mode or prior the quantity or field does not allow"},
      {kPriorSyntax, E, "a prior that is not 'prior <centre> <width>[%]' with a width > 0"},
      {kPartMissing, E, "a structured value without one of its required keywords"},
      {kPartTwice, E, "a keyword given twice in a structured value"},
      {kFormulaSyntax, E, "a formula that is not quoted, or a parameter that is not <name>=<number>"},
      {kSpinSyntax, E, "a spin that is not 0, 1, 3/2, ... (with + or - where a parity is needed)"},
      {kNucleusSyntax, E, "a nuclide that is not <A><Symbol>, or with A < Z"},
      {kUnknownElement, E, "an element symbol that is not in the table (physics/Elements.h)"},
      {kPairSyntax, E, "a particle pair that is not <particle>+<nucleus>[<E_x>]"},
      {kWindowSyntax, E, "fold windows that are not <lo>-<hi> ... <unit> or auto[:tol[:margin]]"},
      // Reactions
      {kReactionSyntax, E, "a reaction expression that does not follow section 3.5"},
      {kReactionOrdinal, E, "an ordinal particle label (a1, p2): name the state by its excitation energy (F4)"},
      {kReactionBalance, E, "a reaction that does not conserve Z and A"},
      {kReactionNotAllowed, E, "a reaction (or pair) the quantity is not defined for"},
      // Column roles
      {kUnknownRole, E, "a column role that is not in the role table"},
      {kRoleUnit, E, "an energy or angle role without its unit in brackets, or a unit on a role without one"},
      {kRoleTwice, E, "a column role given twice"},
      {kRoleRequired, E, "a column the quantity requires is missing"},
      {kRoleNotAllowed, E, "a column role the quantity does not take"},
      // Between records
      {kUnknownReference, E, "a name that no record has (fix: the nearest name of the right kind)"},
      {kReferenceKind, E, "a name of a record of another kind than the field needs"},
      {kUnknownEpoch, E, "target@epoch where the target has no such epoch"},
      {kTargetKind, E, "a target whose kind the quantity does not take"},
      {kDegenerateNuisance, E, "efficiency and norm both free without a prior (F9)"},
      // Model data against <levels>
      {kUnknownLevel, E, "a datum's level ID that no active level of <levels> has, or two levels of its J^pi have"},
      {kLevelJpi, E, "a datum's J^pi that is not the J^pi of the level with its ID (F13)"},
      {kUnknownPair, E, "a reaction or pair that no pair of <levels> has (particles and states)"},
      {kAmbiguousPair, E, "a reaction or pair that two pairs of <levels> have (a data set says which with a pair key)"},
      {kNotAChannel, E, "a datum's pair (or pair, l, s) that is no channel of its level, or no s-wave channel spin of the pair"},
      {kPairKey, E, "entrance_pair / exit_pair where the label names one pair, or a key the label does not name"},
      {kDataFile, E, "a data file that cannot be read, or has fewer columns than its roles"},
      {kNotBuilt, E, "a value the table takes and this engine cannot yet set up (a unit, frame or form of a later step)"},
      {kThmRelativistic, E, "a THM data set or extrapolation in a project with relativistic kinematics"},
      {kThmCovariance, E, "a covariance file on a THM data set (its chi^2 is its experiment's profiled one)"},
      {kFoldConflict, E, "data sets sharing a target with a free parameter differ in beam, detector or quadrature"},
      // Warnings
      {kNotMeasurable, Severity::Warning, "a data set of a quantity that is computed for tests and comparison, not measured"},
      {kBroadLevel, Severity::Warning, "a resonance strength of a level with Gamma / E_r > 1 %: eq. (8) and the area of the resonance differ"},
      {kNearbyLevel, Severity::Warning, "a resonance strength of a level with another level of its J^pi within 10 Gamma"},
  };
  return table;
}

Diagnostic FromRecordDiagnostic(const io::RecordDiagnostic &r, const std::string &section) {
  Diagnostic d;
  d.code = r.code;
  d.line = r.line;
  d.column = r.col;
  d.section = section;
  d.hasFound = !r.found.empty();
  d.found = r.found;
  d.message = r.message;
  return d;
}

std::string Diagnostic::ToJson() const {
  JsonWriter j(0);
  j.BeginObject();
  j.Key("severity").String(SeverityName(severity));
  j.Key("code").String(code);
  if (!file.empty()) j.Key("file").String(file);
  j.Key("line").Integer(line);
  j.Key("column").Integer(column);
  if (!section.empty()) j.Key("section").String(section);
  if (!recordKind.empty()) j.Key("record_kind").String(recordKind);
  if (!record.empty()) j.Key("record").String(record);
  if (!field.empty()) j.Key("field").String(field);
  j.Key("found");
  if (hasFound)
    j.String(found);
  else
    j.Null();
  if (!allowed.empty()) j.Key("allowed").StringArray(allowed);
  j.Key("message").String(message);
  if (hasFix) {
    j.Key("fix").BeginObject();
    j.Key("action").String(fix.action);
    if (fix.action == "insert_field") {
      j.Key("after_line").Integer(fix.afterLine);
    } else {
      j.Key("line").Integer(fix.line);
      j.Key("column").Integer(fix.col);
    }
    j.Key("text").String(fix.text);
    j.EndObject();
  }
  j.EndObject();
  return j.Str();
}

std::string Diagnostic::ToText() const {
  std::string s = file;
  if (line > 0) s += (s.empty() ? "" : ":") + std::to_string(line);
  if (!s.empty()) s += ": ";
  if (!record.empty()) s += record + (field.empty() ? "" : "." + field) + ": ";
  return s + SeverityName(severity) + " " + code + ": " + message;
}

}  // namespace experiment
}  // namespace azure2
