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
      {kFieldDeleted, E, "a field an earlier format 3 had and this one does not (fix: remove it, as the message says)"},
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
      {kPriorCentre, E, "a fitted norm whose prior is centred at 0 or below (its width is relative to the centre)"},
      {kPartMissing, E, "a structured value without one of its required keywords"},
      {kPartTwice, E, "a keyword given twice in a structured value"},
      {kFormulaSyntax, E, "a formula that is not quoted, or a parameter that is not <name>=<number>"},
      {kSpinSyntax, E, "a spin that is not 0, 1, 3/2, ... (with + or - where a parity is needed)"},
      {kSpinOperator, E, "the operator of a spin observable that does not follow its grammar (b:y, b:y*t:y, b:T20, b:y>e1:y)"},
      {kNucleusSyntax, E, "a nuclide that is not <A><Symbol>, or with A < Z"},
      {kUnknownElement, E, "an element symbol that is not in the table (physics/Elements.h)"},
      {kPairSyntax, E, "a particle pair that is not <particle>+<nucleus>[<E_x>]"},
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
      {kAmbiguousPair, E, "a reaction or pair that two pairs of <levels> have (a data set or a potential says which with a pair key)"},
      {kNotAChannel, E, "a datum's pair (or pair, l, s) that is no channel of its level, or no s-wave channel spin of the pair; a gamma_decay's final state that is no photon channel of its level"},
      {kGammaMixture, E, "a gamma_decay level that also has a non-zero photon width in <levels> (write 0 there or delete the record)"},
      {kBranchingSum, E, "gamma_decay branchings that are all fixed and do not sum to 1, or a fitted one without a remainder"},
      {kBranchingRemainder, E, "gamma_decay branchings whose sum exceeds 1, leaving a negative remainder"},
      {kMixingOneMultipole, E, "a gamma_decay mixing ratio for a final state with one multipolarity in the level's channels"},
      {kGammaRecord, E, "a gamma_decay record that names a level twice, a final state twice, a level without photon channels, or a malformed remainder prior"},
      {kPairKey, E, "entrance_pair / exit_pair / pair_key where the label names one pair, or a key the label does not name"},
      {kDataFile, E, "a data file that cannot be read, holds no data row, or has a row with fewer columns than its roles"},
      {kDataRow, E, "a row of a data file that does not parse (four numbers), or an angular-coefficient row without a Legendre order 0-20"},
      {kDataError, Severity::Warning, "a row inside a data set's window with an uncertainty <= 0 and no covariance file: the chi^2 skips a 0 and squares a negative one"},
      {kCovarianceFile, E, "a covariance file that cannot be read, does not list the data set's file with its row count, or has a variance <= 0 on a row in the window"},
      {kFinalStateNotLevel, E, "the final state of a photon pair with external capture that a data set or extrapolation exits into is not one level of <levels> (none, or several within 1 keV)"},
      {kFinalStateShared, E, "two photon pairs with external capture whose final states are the same level"},
      {kFinalStateUnbound, E, "a capture final state above a particle threshold: with external capture, at or above the threshold of a pair with a channel in its J^pi group; or above the entrance threshold of a capture asked into it (not computed by this version)"},
      {kNotBuilt, E, "a value the table takes and this engine cannot yet set up (a unit, frame or form of a later step)"},
      {kThmRelativistic, E, "a THM data set or extrapolation in a project with relativistic kinematics"},
      {kThmCovariance, E, "a covariance file on a THM data set (its chi^2 is its experiment's profiled one)"},
      {kThmExperimentUnknown, E, "a thm field naming an experiment that <thm> has no experiment[<name>] line for"},
      {kThmOldList, E, "an old member or weight list in <thm> (segments=, weight[<k>]=, weightTest[<k>]=): in format 3 the data set's thm and thm_weight fields say it"},
      {kFoldConflict, E, "data sets sharing a target with a free parameter differ in beam, detector or quadrature"},
      {kLayerStopping, E, "a yield on a layer whose stopping is not positive and finite at an energy of its sub-point grid (set-up)"},
      {kLayerExitBelowModel, E, "a yield on a layer the beam does not leave above the lowest energy the model is evaluated at (set-up)"},
      {kStoppingProjectile, E, "a target whose stopping depends on the projectile (evaluation, table, scale) crossed by a neutral particle, an ion above Z = 92, or several projectiles in one data set"},
      {kStoppingActive, E, "a stopping evaluation or table whose composition does not contain the active nucleus's element"},
      {kStoppingScaleConflict, E, "two targets that share a stopping scale (one material and projectile) and write it differently"},
      {kStoppingData, E, "a stopping evaluation in an engine built without SRIM's data files where neither --srim-data-dir=DIR nor the switch AZURE2_SRIM_DATA_DIR names a directory with them, or a stopping table that cannot be read (set-up)"},
      {kReferenceTable, E, "the reference table of a ratio_to_reference record: unreadable, not increasing in energy, a value <= 0, or not covering an energy of the data set or grid (build)"},
      {kStoppingRange, E, "a target's stopping asked for at an energy outside its range: SRIM's above 1e7 keV/u, or a table's ends (set-up)"},
      {kLabJacobianAtLimit, E, "a laboratory-frame differential datum at the kinematic maximum lab angle, where dOmega_lab/dOmega_cm is 0 and the lab-frame model diverges (set-up)"},
      {kPotentialTwice, E, "two potential records for one pair"},
      {kPotentialPairKind, E, "a potential record on a photon or beta pair (a potential belongs to a particle pair)"},
      {kChannelRadiusTwice, E, "two channel_radius records for one selection of channels (pair, l, s, J)"},
      {kChannelRadiusPairKind, E, "a channel_radius record on a photon or beta pair (a radius belongs to a particle channel)"},
      {kChannelRadiusThm, E, "a channel_radius record on the entrance pair of a THM data set (the vertex is built at the pair's radius)"},
      {kIsospinRecord, E, "a record of <isospin> that another record or <levels> contradicts: a nucleus or a level named twice, no or a second reference, a level ID no level has, a mix or a factor on a level without an isospin_level record, a soft tie with one width, a profiled number"},
      {kIsospinPair, E, "a pair of an <isospin> record that is not a particle pair, or a threshold pair that does not form the record's nucleus"},
      {kPotentialRelativistic, E, "a potential record in a project with relativistic kinematics"},
      {kThmHybridPotential, E, "a THM data set or extrapolation whose entrance pair has a potential record"},
      {kOldPotentialBlock, E, "the old <potential> block in a format-3 project (AZURE2 --convert writes potential records)"},
      {kRemovedConfigKey, E, "a <config> key that format 3 no longer has (useAdaptiveGrid: every fold is on the shared mesh; use_gsl_coul: GSL's Coulomb functions are gone)"},
      {kUnknownSection, E, "a section format 3 does not have (a newer AZURE2's section, a misspelt name), or a section inside one that does not hold it: nothing would read it"},
      {kUnreadLine, E, "a <config> key format 3 does not have, a <config> line after the eleven entries that is not key=value, or a line of <mcmc> outside <parameters>: nothing would read it"},
      {kConfigOptionValue, E, "an option of the calculation in <config> (use_brune, transform, ...) with a value other than 0, 1, true or false"},
      {kConfigOptionsContradict, E, "two options of the calculation that cannot hold together, whatever gave them: a flag, a keyword, a <config> key or the default (use_park=1 with use_brune=0, use_rmc=1 with use_brune=1, ec_convection_current=0 with use_long_wavelength=1)"},
      {kParameterFileNoMatch, E, "a parameter file (param.sav as a run's start) none of whose names is a parameter of the run"},
      {kParameterFileNames, Severity::Warning, "a parameter file with names the project has no parameter for, or without names of parameters it has"},
      {kParameterSettingUnknown, Severity::Warning, "a <parameterSettings> row that names no parameter of the project (ignored)"},
      {kSoTwoSpins, E, "a spin-orbit term on a pair whose two particles both have spin (the coupled treatment is not implemented)"},
      {kSoZeroSpin, E, "a spin-orbit term acting on a particle of spin 0 (it would vanish)"},
      {kSoRelativistic, E, "a spin-orbit term in a project with relativistic kinematics"},
      {kSoThm, E, "a THM data set or extrapolation whose entrance pair's potential has a spin-orbit term"},
      {kSoNeedsJ, E, "channel functions of a pair with a spin-orbit term asked for without J (pyazr, the API)"},
      {kSoNoTerm, E, "channel functions of a pair without a spin-orbit term asked for with a J (pyazr, the API)"},
      // Warnings
      {kNotMeasurable, Severity::Warning, "a data set of a quantity that is computed for tests and comparison, not measured"},
      {kBroadLevel, Severity::Warning, "a resonance strength of a level with Gamma / E_r > 1 %: eq. (8) and the area of the resonance differ"},
      {kNearbyLevel, Severity::Warning, "a resonance strength of a level with another level of its J^pi within 10 Gamma"},
      {kDataFileNotRead, Severity::Warning, "a data file that cannot be read in a run that does not use the data (an extrapolation): not read"},
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
