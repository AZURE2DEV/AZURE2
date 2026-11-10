#ifndef AZURE2_EXPERIMENT_SCHEMA_H
#define AZURE2_EXPERIMENT_SCHEMA_H

#include <climits>
#include <string>
#include <utility>
#include <vector>

#include "azure2/experiment/Conditions.h"
#include "azure2/experiment/Diagnostics.h"
#include "azure2/experiment/Quantities.h"
#include "azure2/experiment/Values.h"
#include "azure2/io/Record.h"

namespace azure2 {
namespace experiment {

/*!
 * The table of format 3 (docs/upgrade/project_format3.md sections 3, 6, 7): the record kinds, their sections and
 * fields with value types (Schema.cpp), the quantities (Quantities.cpp) and the units (Units.cpp); one table-driven
 * validator over records read by io::ReadRecords; and the table exported as JSON (section 6.1).
 *
 * A new field is one FieldDef row of its kind (and a quantity's QuantityField where it depends on the quantity); a new
 * quantity one Quantity row; a new unit one Unit row; a new value type (rare) a parser in Parsers.h and a case in
 * ParseValue.  Each comes with the test named in CLAUDE.md's recipes.
 */

struct FormSpec;

/// The value type of a field or of a part of a structured field, with its parameters.
struct ValueSpec {
  ValueSpec() = default;
  explicit ValueSpec(ValueType t, std::string dim = "") : type(t), dimension(std::move(dim)) {}

  ValueType type = ValueType::Text;
  /// Unit, NumberWithUnit, Range, Nuisance, Grid, Windows: the dimension of the unit ("" for a dimensionless
  /// nuisance, which takes no unit token); "@quantity": a unit the record's quantity allows.
  std::string dimension;
  std::vector<std::string> values;  ///< Enum: the words allowed; Reference: the record kinds it may name
  std::vector<std::pair<std::string, std::string>> aliases;  ///< Enum: (alias, word)
  long min = LONG_MIN;      ///< Integer: the smallest value allowed
  bool positive = false;    ///< Number, NumberWithUnit: > 0
  bool epoch = false;       ///< Reference: <name>@<epoch> allowed
  bool profiled = false;    ///< Nuisance: 'profiled' allowed (where, the validator decides)
  bool none = false;        ///< Grid: 'none' allowed
  std::vector<std::string> refKinds;  ///< Nuisance: a name of a record of these kinds may stand instead
  const FormSpec *form = nullptr;     ///< Form
};

struct PartSpec {
  std::string name;     ///< the part's name in Value::parts (the keyword, for a keyword part)
  bool keyword = false; ///< written as "<name> <value>"; positional parts come first, in order
  ValueSpec spec;
  bool optional = false;
  std::string defaultText;  ///< the value used when an optional part is absent (logged by the engine)
  std::string doc;
};

struct Alternative {
  std::string keyword;  ///< the leading word; "" for the alternative chosen when no keyword matches
  std::vector<PartSpec> parts;
  std::string doc;
};

struct FormSpec {
  std::vector<Alternative> alternatives;
};

/// The grammar of a value spec in words ("<value> [<unit>] fixed|fit ...").
std::string Grammar(const ValueSpec &spec);

struct FieldDef {
  std::string key;
  std::string title;
  ValueSpec spec;
  bool byQuantity = false;  ///< the record's quantity row decides where it applies (Quantity::fields)
  Condition applies;        ///< where the field may be given (Always for most)
  Condition required;       ///< where it must be given, among the records it applies to (Always, or never: negated)
  bool repeatable = false;
  std::string defaultText;  ///< an optional field's value when absent (logged by the engine when used)
  std::string doc;
};

struct RecordKind {
  std::string kind;
  std::vector<std::string> sections;
  std::string title;
  std::string doc;
  std::vector<FieldDef> fields;
  std::vector<std::vector<std::string>> oneOf;  ///< at least one field of each group
  /// Keys of other kinds that a user may write here by mistake, with the reason they do not belong.
  std::vector<std::pair<std::string, std::string>> refused;
  const FieldDef *Find(const std::string &key) const;
};

/// The records of one section (as io::ReadRecords returns them), with the section's name.
struct SectionRecords {
  std::string section;
  std::vector<io::Record> records;
};

struct ValidatedField {
  std::string key;
  int line = 0;
  int col = 0;
  Value value;
};

/// A record with its fields parsed into values; fields that did not parse are absent (their diagnostics say why).
struct ValidatedRecord {
  std::string kind;
  std::string name;
  std::string section;
  int line = 0;
  int lastLine = 0;                     ///< the record's last line (where a missing field is inserted)
  const Quantity *quantity = nullptr;   ///< the row of its `quantity` field, when it has one that parsed
  std::vector<ValidatedField> fields;   ///< in file order; a repeatable field may appear more than once
  bool Has(const std::string &key) const;
  const Value *Get(const std::string &key) const;  ///< the first value of `key`, or nullptr
  std::vector<const Value *> GetAll(const std::string &key) const;
};

struct ValidationResult {
  std::vector<Diagnostic> diagnostics;   ///< every problem, ordered by line
  std::vector<ValidatedRecord> records;  ///< every record of a known kind, in input order
  bool Ok() const;                       ///< no diagnostic of severity error
  const ValidatedRecord *Find(const std::string &name) const;
};

class Schema {
 public:
  static constexpr const char *kVersion = "3.0";

  /// The record kinds, in schema order.
  static const std::vector<RecordKind> &Kinds();
  static const RecordKind *FindKind(const std::string &kind);
  /// The format-3 sections and the record kinds each holds (section 3.2).
  static const std::vector<std::pair<std::string, std::vector<std::string>>> &Sections();
  /// The record kinds of a section, empty for a section that holds no records: the `kinds` argument of
  /// io::ReadRecords.  Within a section, kinds and the field keys of its records are disjoint (schema_test checks it),
  /// so a line is a header or a field by its first token (section 3.1); across sections they are not (a data set's
  /// field `target` names a record of kind `target` in <targets>), so a section is read with its own kinds only.
  static std::vector<std::string> KindsOfSection(const std::string &section);

  /*!
   * Validate the records of a project's format-3 sections against the table and return every problem
   * (section 6.2): unknown kinds and fields, fields given twice, missing required and conditionally required fields,
   * value grammars, units, column roles, reaction expressions, the fields and units the quantity allows, '?'
   * placeholders, duplicate names, references to missing records or records of the wrong kind.  What needs
   * <levels> or the data files is not checked here (Schema::Build, step 3).  Never throws for bad input.
   */
  static ValidationResult Validate(const std::vector<SectionRecords> &sections);

  /// The table as the JSON document of section 6.1 (schema_version "3.0"), deterministic, ending without a newline.
  static std::string ToJson();
};

}  // namespace experiment
}  // namespace azure2

#endif
