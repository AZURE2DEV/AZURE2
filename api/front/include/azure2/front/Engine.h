#ifndef AZURE2_FRONT_ENGINE_H
#define AZURE2_FRONT_ENGINE_H

// The engine as a front end in its own process sees it (docs/upgrade/qt_gui.md): a project file read, edited,
// checked and written with the engine's own grammars, and the program's commands run with progress, messages, a stop
// request and results returned as data.  Plain C++ only: no Qt type, no engine header.  The Qt setup program (gui/)
// includes this header and no other of the engine; its test (api/front/tests/front_test.cpp) runs without Qt.
//
// Every function below calls what the command line and pyazr call: io::Document and its record, <config> and
// <levels> grammars (the lossless document pyazr's Project edits), ValidateProjectBytes (AZURE2 --validate),
// convert::ConvertToFormat3 (AZURE2 --convert), RunProgram (main's run, src/ProgramMain.cpp), the schema's JSON
// (AZURE2 --schema).  Nothing here holds state between calls but the Project object a caller owns.

#include <atomic>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace azure2 {
namespace front {

/// One diagnostic of the engine's reader and validator (the codes AZURE2 --validate prints; the schema's
/// "diagnostic_codes" lists them with their meaning).
struct Diagnostic {
  std::string severity;  ///< "error", "warning" or "ask"
  std::string code;      ///< E_..., W_..., A_...
  std::string file;
  int line = 0;          ///< 1-based line of the file, 0 when the diagnostic has none
  std::string section;
  std::string record;
  std::string field;
  std::string message;
  std::string text;      ///< the one line --validate prints: "file:line: record.field: message"
};

/// What --validate says of a project.
struct Validation {
  bool ok = false;             ///< no error
  bool format3 = true;         ///< false: the file is no format-3 project; `refusal` names the conversion
  std::string refusal;
  std::vector<Diagnostic> diagnostics;
  std::string stoppingReport;  ///< the targets' stopping sources (or why a target has none)
};

/// One field line of a record: its key and its value as written (tokens with their quotes).
struct Field {
  std::string key;
  std::string value;
  std::string comment;  ///< the "# ..." after the value, without the '#'; empty for none
};

/// One record of a format-3 section: `kind name` and its field lines in file order.
struct Record {
  std::string kind;
  std::string name;
  std::vector<Field> fields;
  int line = 0;  ///< the header's line in the file as last read or written (0 for a record not yet written)
  /// The value of the first field `key`, or "".
  std::string Value(const std::string &key) const;
};

/// One entry of <config>: a positional line (the eleven at the top) or a key=value line.
struct ConfigEntry {
  std::string name;  ///< the positional entry's name (io::ConfigPositionalNames) or the key
  std::string value;
  bool positional = false;
};

/// One row of <levels> as the engine reads it (io::LevelRow): the 31 to 34 columns of io::LevelColumnNames() as
/// written, in that order ("" for an optional column the row does not write).
struct LevelRow {
  std::vector<std::string> tokens;
};

/// A format-3 project in memory: the file's lossless document.  An edit replaces only the lines of what it changes;
/// a project opened and saved without an edit is written back byte for byte.
class Project {
 public:
  Project();
  ~Project();
  Project(Project &&) noexcept;
  Project &operator=(Project &&) noexcept;
  Project(const Project &) = delete;
  Project &operator=(const Project &) = delete;

  /// The project file at `path`, or false with `why`: the file cannot be read, or it is not format 3 (the engine's
  /// refusal, io::NotFormat3Message, which names AZURE2 --convert; ConvertProject below does it).
  static bool Open(const std::string &path, Project &out, std::string &why);
  /// An empty project to be saved at `path`: the eleven positional <config> entries (first value of each; output/
  /// and checks/), format=3, empty <levels>, <data>, <extrapolations> and <targets> (as pyazr's empty project).
  static Project Empty(const std::string &path);
  /// The project of `text` at `path` (a file not read from disk; used by tests and by New).
  static Project FromText(const std::string &text, const std::string &path);

  const std::string &Path() const;
  /// The directory relative paths of the project are read from.
  std::string Directory() const;
  /// The whole text as it would be saved.
  std::string Text() const;
  /// The text differs from what was last read from or written to Path().
  bool Modified() const;
  /// Write Text() to `path` (empty: Path(); a new path becomes the project's).  False with `why`.
  bool Save(const std::string &path, std::string &why);

  /// The names of the document's top-level sections, in file order.
  std::vector<std::string> Sections() const;
  bool HasSection(const std::string &name) const;
  /// The lines between a section's tags (empty for a section the file has not).
  std::vector<std::string> SectionBody(const std::string &name) const;
  /// Replace a section's body (the section is added at the end of the file when it has none; an empty body of a
  /// section the file has not adds nothing).
  void SetSectionBody(const std::string &name, const std::vector<std::string> &lines);

  /// <config>, read as the engine reads it (io::ReadConfigSection).
  std::vector<ConfigEntry> Config() const;
  /// The value of a <config> entry, or `absent`.
  std::string ConfigValue(const std::string &name, const std::string &absent = "") const;
  /// Set (or, `remove`, delete) one <config> entry with the engine's edit (io::ConfigEdits: the line's comment kept).
  /// False with `why` for a name the engine has not or a value it would not read.
  bool SetConfig(const std::string &name, const std::string &value, std::string &why, bool remove = false);

  /// The rows of <levels> in file order; a blank line between two rows (a new level) is an empty row.
  std::vector<LevelRow> Levels() const;
  /// Replace the rows of <levels>.  A row equal to the one at its place keeps its line; an edited row is written by
  /// io::EditLevelRow over the old line (its layout kept); a new row is its tokens joined by blanks.  Every row is
  /// read back with the engine's reader: false with `why` (and nothing changed) for one it refuses.
  bool SetLevels(const std::vector<LevelRow> &rows, std::string &why);

  /// The records of a section of records (`data`, `extrapolations`, `targets`, `modelData`, `potentials`,
  /// `gammaDecays`, `channelRadii`, `isospin`; the schema's "sections").  Empty for a section the file has not.
  std::vector<Record> Records(const std::string &section) const;
  /// Replace the record named record.name of `section` (its comments kept), or append it at the end of the section
  /// (the section is made when the file has none).  Written by io::WriteRecords.  False with `why` for a kind the
  /// section does not take or a name or key that is no token.
  bool SetRecord(const std::string &section, const Record &record, std::string &why);
  /// Rename a record (its references elsewhere are not renamed: --validate names them).
  bool RenameRecord(const std::string &section, const std::string &from, const std::string &to, std::string &why);
  /// Remove the record `name` of `section` with its comment lines.
  bool RemoveRecord(const std::string &section, const std::string &name, std::string &why);
  /// Move the record `name` of `section` one place up (`delta` -1) or down (+1) among the section's records: its
  /// lines (comments above it included) change places with its neighbour's, every other line stays.  False with
  /// `why` at either end or for a name the section has not.
  bool MoveRecord(const std::string &section, const std::string &name, int delta, std::string &why);
  /// The records of `section` switched off (SetRecordOn), each as Records() gives a record (its fields as written; line
  /// = the line of its header in the file).  Not seen by the engine, Records() or Validate().
  std::vector<Record> RecordsOff(const std::string &section) const;
  /// Switch record `name` of `section` off or on.  Off: each of its lines (the comments above it, the header, the fields)
  /// is prefixed with RecordOffMark(), so that the engine reads them as comments and the record does nothing; on: the
  /// mark is taken away again, which gives back the lines byte for byte.  No other line moves.  False with `why` for a
  /// record that is not there (on, or off), or when an active record of that name exists (on).
  bool SetRecordOn(const std::string &section, const std::string &name, bool on, std::string &why);

  /// AZURE2 --validate of the text as it stands (data files are read relative to Directory()).
  Validation Validate() const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

/// The table of project format 3 as JSON (AZURE2 --schema): sections, record kinds and their fields with titles,
/// docs, types, closed choices and units, quantities, units, diagnostic codes.
std::string SchemaJson();
/// The mark a record switched off carries at the start of each of its lines ("#off "): a comment of format 3.
std::string RecordOffMark();

// ---- The shape of a value: what a form offers to edit it (docs/upgrade/qt_gui.md section 11) ----

struct PartInfo;
/// One alternative of a structured value (the schema's FormSpec): led by `keyword`, or "" for the one taken when no
/// keyword matches.
struct AlternativeInfo {
  std::string keyword;
  std::string doc;
  std::vector<PartInfo> parts;  ///< positional parts first, in order, then the keyword parts
};
/// The value type of a field or of a part of a structured field, with what its grammar lets a form choose (the
/// schema's ValueSpec).  Every word a form offers comes from here: units, enum words, alternatives, grid kinds.
struct ValueInfo {
  std::string type;       ///< the value type's name ("number", "nuisance", "form", "grid", ...: ValueTypeName)
  std::string grammar;    ///< the grammar in words, for a placeholder
  std::string dimension;  ///< of the unit ("" none; "@quantity": a unit of the record's quantity)
  std::vector<std::string> units;  ///< the units of `dimension` (empty for "@quantity": QuantityUnits)
  bool unitOfQuantity = false;
  std::vector<std::string> values;  ///< enum: the words; reference: the record kinds it may name
  bool positive = false;            ///< number: > 0
  bool hasMin = false;              ///< integer: a smallest value
  long min = 0;
  bool epoch = false;      ///< reference: <name>@<epoch> allowed
  bool profiled = false;   ///< nuisance: 'profiled' allowed
  bool none = false;       ///< grid: 'none' allowed
  std::vector<std::string> refKinds;          ///< nuisance: the kinds of record whose name may stand instead
  std::vector<AlternativeInfo> alternatives;  ///< form: its alternatives
};
/// One part of an alternative: positional (in order) or written "<name> <value>".
struct PartInfo {
  std::string name;
  bool keyword = false;
  bool optional = false;
  std::string defaultText;  ///< the value used when an optional part is absent
  std::string doc;
  ValueInfo value;
};

/// One field of a record kind (the schema's row, in the kind's order): what a form needs to offer it.
struct FieldInfo {
  std::string key;
  std::string title;
  std::string type;      ///< the value type ("enum", "number_with_unit", "nuisance", "path", "structured", ...)
  std::string grammar;   ///< the value's grammar, for a placeholder
  std::string doc;
  bool required = false; ///< required wherever the field applies (a conditional requirement is "requiredWhen")
  std::string requiredWhen;  ///< in words, "" when unconditional
  std::string appliesWhen;   ///< in words, "" when the field applies to every record of the kind
  bool byQuantity = false;   ///< the record's quantity decides whether it applies
  bool rare = false;     ///< set rarely: a form shows it one click away unless the record sets it
  bool repeatable = false;
  std::string defaultText;
  std::vector<std::string> values;  ///< enum: the words; reference: the kinds it names
  std::vector<std::string> units;   ///< the units of its dimension; for "@quantity" those of the record's quantity
  bool unitOfQuantity = false;      ///< units come from the record's quantity (QuantityUnits)
  ValueInfo value;                  ///< the value's whole shape (parts of a structured value, grid kinds, ...)
  /// A field a table of records shows as a column (KindInfo::summaryFields gives their order).
  bool summary = false;
  /// The key of the field this one belongs with on a form (the field its "applies"/"required" condition names, as
  /// a stopping's phase with the stopping, a profile's truncation with the profile); "" for a field of its own.
  std::string group;
};
/// One record kind: its sections, title, doc and fields in the schema's order.
struct KindInfo {
  std::string kind;
  std::string title;
  std::string doc;
  std::vector<std::string> sections;
  std::vector<FieldInfo> fields;
  /// The columns a table of these records shows, in order: the required fields first (the quantity's and the
  /// reaction's), then the references to other records, the windows and grids, the main numbers; never free text or
  /// the column roles.  Every FieldInfo with `summary`.
  std::vector<std::string> summaryFields;
  /// The field that switches a record on and off (an optional yes/no field whose default is the first word: a data
  /// set's or an extrapolation's "active"), or "".
  std::string includeField;
  std::string includeOn, includeOff;  ///< its two words (on first)
  /// The groups at least one field of which a record must give (a beam: spread or profile).
  std::vector<std::vector<std::string>> oneOf;
};
std::vector<KindInfo> RecordKinds();
/// The sections of records and the kinds each holds, in the schema's order.
std::vector<std::pair<std::string, std::vector<std::string>>> RecordSections();
/// The quantities a record kind may use (names, schema order), and the units of one quantity.
std::vector<std::string> QuantityNames(const std::string &kind);
std::vector<std::string> QuantityUnits(const std::string &quantity);
/// One quantity of the schema (section 7.2): its name, the title and doc a form shows, its units.
struct QuantityInfo {
  std::string name;
  std::string title;
  std::string doc;
  std::vector<std::string> units;
};
/// The quantities a record kind may use, in the schema's order.
std::vector<QuantityInfo> Quantities(const std::string &kind);
/// The diagnostic codes and their meaning.
std::vector<std::pair<std::string, std::string>> DiagnosticCodes();

/// A value split into the pieces a form edits, each as written (numbers keep their text).  Which members hold it
/// depends on `type`; FieldValueText writes it back.  Made by ParseFieldValue with the engine's own parsers.
struct NamedValue;
struct ValueParts {
  std::string type;    ///< the value type (ValueInfo::type)
  /// Number, integer, unit, enum, name, text, spin, nucleus, pair, composition-free scalars, reaction: the value's
  /// text as written; reference: the record's name; path: the path (without quotes); nuisance with `reference`: the
  /// name of the record that stands for it.
  std::string text;
  bool quoted = false;  ///< path: written in quotes
  std::string epoch;    ///< reference: the epoch after '@' ("" none)
  std::string number;   ///< number with unit, nuisance: the number
  std::string unit;     ///< number with unit, range, nuisance, grid: the unit token ("" none)
  std::string lo, hi;   ///< range; a grid's span
  // nuisance
  bool reference = false;  ///< written as the name of a record (`text`)
  bool fit = false;        ///< mode fit (else fixed)
  bool profiled = false;   ///< 'profiled' (no number)
  bool hasPrior = false;
  std::string priorCentre, priorWidth;
  bool priorPercent = false;  ///< the width is a percentage of the centre
  std::string gridMax;        ///< "" none
  // formula: the expression in x and its parameters (name, number)
  std::string expression;
  std::vector<std::pair<std::string, std::string>> params;
  std::vector<std::string> numbers;  ///< number list; a list grid's values
  std::vector<std::pair<std::string, std::string>> composition;  ///< (element symbol, atoms)
  // columns
  bool fromFile = false;  ///< the roles are in the data file's header
  std::vector<std::pair<std::string, std::string>> columns;  ///< (role, unit or "")
  // grid
  std::string gridKind;  ///< a word of GridKinds()
  std::string frame;     ///< a word of Frames()
  bool byPoints = false; ///< a range grid by its number of points (else by its step)
  std::string step, points;
  bool log = false;
  std::string tolerance, maxPoints, startStep;
  // form
  std::string alternative;        ///< the alternative's keyword
  std::vector<NamedValue> parts;  ///< the parts given, by their name, in the alternative's order
  const ValueParts *Part(const std::string &name) const;
};
struct NamedValue {
  std::string name;
  ValueParts value;
};

/// The value `text` of field `key` of a record of `kind` (with the record's `quantity` for the units "@quantity", ""
/// when unknown) split into its pieces by the engine's parser (the one --validate reads it with).  False with `why`,
/// the validator's message for that field, when the parser refuses the text.
bool ParseFieldValue(const std::string &kind, const std::string &key, const std::string &text,
                     const std::string &quantity, ValueParts &out, std::string &why);
/// The text of a value from its pieces, as a writer of format 3 puts it (one blank between tokens, quotes where the
/// reader needs them; a nuisance by the engine's NuisanceText).  Parsing that text gives the pieces back.
std::string FieldValueText(const std::string &kind, const std::string &key, const ValueParts &parts);
/// The same for a value of shape `info` (a part of a structured value).
std::string ValueText(const ValueInfo &info, const ValueParts &parts);
/// The field's value as a cell of a table shows it, short: a nuisance "1 (fit)", a form its alternative and numbers.
std::string ShortValueText(const std::string &kind, const std::string &key, const std::string &text);

/// One kind of grid (section 3.8) and what it is given by.
struct GridKind {
  std::string word;
  std::string doc;
  bool span = false;      ///< lo, hi
  bool spacing = false;   ///< step or points [log]
  bool list = false;      ///< a list of values
  bool file = false;      ///< a file of values
  bool adaptive = false;  ///< tolerance, max_points, start_step
  bool unit = true;       ///< unit and frame follow ('none' has none)
};
/// The kinds of grid ('none' last, when `withNone`).
std::vector<GridKind> GridKinds(bool withNone);
/// The frames of an energy or angle ("lab", "cm").
std::vector<std::string> Frames();
/// The units of a dimension ("energy", "angle", ...), in the unit table's order.
std::vector<std::string> DimensionUnits(const std::string &dimension);
/// The column roles of a data file (section 3.4.1): the role and the dimension of its unit ("" for a role without).
std::vector<std::pair<std::string, std::string>> ColumnRoles();

/// The columns a new data set of `quantity` starts with: its first energy role and, for a quantity with an angle, its
/// first angle role (each in the first unit of its dimension), then the value and its error.
std::vector<std::pair<std::string, std::string>> DefaultColumns(const std::string &quantity);

/// The shape of an apparatus record's kernel for a small preview: x and y of a curve, its title and axis label.
struct ShapeCurve {
  std::string title;
  std::string xLabel;
  std::vector<double> x, y;
};
/// The shape a record's kernel gives where the engine can evaluate it from the record alone: a beam's energy profile
/// (the sum of its skewed-Gaussian components, truncated as the record says: the kernel's own weight,
/// experiment::BeamProfileWeight, over its support, lab MeV).  False (`why` says why) for a record without one, or
/// whose values do not read; a Doppler broadening needs the pair and the energy of a point and has none here.
bool PreviewShape(const Record &record, int points, ShapeCurve &out, std::string &why);

/// One pair of the project's <levels> (its first row), as a reaction picker lists it.
struct PairChoice {
  int key = 0;
  std::string label;       ///< "p+15N", "g+16O[6.13]" (experiment::PairLabel)
  int kind = 0;            ///< 0 particle pair, 10 photon, 20 beta
  std::string nucleus;     ///< the compound nucleus it forms ("16O")
  double threshold = 0;    ///< MeV
};
std::vector<PairChoice> Pairs(const Project &project);
/// The label of the reaction from pair `entrance` to pair `exit` ("15N(p,a)12C"), `exit` 0 the total capture of
/// the entrance ("12C(p,g)13N[total]"); "" for a key the project has not.  The converter's ReactionLabel.
std::string ReactionLabel(const Project &project, int entrance, int exit);
/// The pair keys a one-term reaction names (0 for a side that names no pair or several: then the record's
/// entrance_pair / exit_pair decides); exit -1 for a total capture.  False for a sum, a ratio or a text the parser
/// refuses.
bool ReactionPairs(const Project &project, const std::string &reaction, int &entrance, int &exit);

/// The <config> keys and positional entries as JSON: [{"name", "kind", "title", "group", "doc", "values",
/// "positional", "option", "flag", "negation", "default", "rare", "absent"}] from io::ConfigPositionals() and
/// io::ConfigKeys() in their order.
std::string ConfigTableJson();

/// One command of the program (RunCommands()): `name` as on the command line, whether this build has it.
struct Command {
  std::string name;
  std::string help;
  bool built = true;
};
std::vector<Command> Commands();
/// One option of a command (RunQuestions()) or of the program (a fitter's switch, --max-calls, --threads ...).
struct Option {
  std::string option;     ///< "--covariance-band", "--temperature-min=" (with '=' when it takes a value)
  std::string kind;       ///< "flag", "number", "whole", "file", "word"
  std::vector<std::string> commands;  ///< the commands that take it; empty: every command
  bool required = false;
  double low = 0, high = 0, defaultValue = 0;
  std::string help;
  std::vector<std::string> words;  ///< "word": the words it takes ("rwa", "physical")
  bool lowOpen = false;            ///< the range excludes `low` ("--delta-chi2": above 0)
  /// The unit of a number, as the row's help names it last ("the rate's lowest temperature, GK": "GK"); "" for none.
  std::string unit;
  /// When the option is asked, in words ("with --covariance-band"); "" for always.
  std::string when;
  /// The flag it applies with ("--covariance-band" for --scale-covariance); "" for none.
  std::string needs;
  /// For a required option that is one of two ways: the other way ("--temperature-file"); "" for none.
  std::string instead;
};
std::vector<Option> Options();

/// The fitters of this build (the rows of Fitters()): name, its switch ("" for the default MIGRAD), its doc.
struct Fitter {
  std::string name;
  std::string option;
  std::string doc;
};
std::vector<Fitter> Fitters();

/// What a run is asked: the command line `AZURE2 <command> <options...> <project>` with `threads` OpenMP threads.
struct RunRequest {
  std::string command;               ///< "calculate", "fit", "extrapolate", "minos", "rate", "mcmc"
  std::vector<std::string> options;  ///< as on the command line ("--covariance-band", "--temperature-min=0.01")
  std::string project;               ///< the project file (saved: the run reads the file, as the program does)
  int threads = 1;
};

/// One progress record of a run (the console's "Iteration:" lines): the objective after an evaluation.
struct RunProgress {
  std::string stage;   ///< "fit", "hesse", "minos", "mcmc" (a step of the sampler), or "evaluation" outside a fit
  long evaluation = 0; ///< the run's evaluation number, the fitter's iteration, or the sampler's step (from 1)
  long calls = 0;      ///< objective evaluations of the fit so far; "mcmc": the steps asked
  double objective = 0;  ///< "mcmc": the walkers' mean log-posterior after the step
  double best = 0;     ///< lowest objective of the fit so far (NaN before the first); "mcmc": the mean log-likelihood
};

/// The receiver of a run's messages and progress.  Called from the run's threads (Progress from several at once):
/// an implementation must be thread-safe and must not call back into the engine.
class RunObserver {
 public:
  virtual ~RunObserver();
  /// One line the program prints (a line ending in '\r' -- a progress line the console overwrites -- is passed
  /// with `transient` set).
  virtual void Message(const std::string &line, bool transient) = 0;
  virtual void Progress(const RunProgress &) {}
};

/// What a run gives back as data, beside the files it writes in the project's output directory.
struct DataSetChiSquared {
  int segment = 0;  ///< the console's "Segment #" (the data set's place among the active ones, from 1)
  double chiSquared = 0;
  int points = 0;
};
struct RunResult {
  int exitCode = 0;  ///< the program's: 0, 1 aborted or refused, 2 a refused command line, 3 not format 3, -1 failed
  bool stopped = false;               ///< the stop request was set during the run
  bool haveChiSquared = false;
  double totalChiSquared = 0;
  std::vector<DataSetChiSquared> dataSets;
  std::string outputDirectory;        ///< the project's output directory (as the run opened it)
  std::vector<std::string> outputFiles;  ///< files of that directory written by the run (names, sorted)
};

/// Run one command as the program does (RunProgram, src/ProgramMain.cpp) in the calling thread: the same lines, the
/// same files.  `stop` (may be null) is read by the run: when set, a fit ends at its best point and writes its outputs
/// there.  Not reentrant: one run at a time per process (the OpenMP thread count is the calling thread's).
RunResult Run(const RunRequest &request, RunObserver &observer, const std::atomic<bool> *stop);

/// AZURE2 --convert OLD NEW (format 1 or 2 to format 3, with the side files): false with `report` saying why.
/// `report` holds the converter's notes and questions either way.
bool ConvertProject(const std::string &from, const std::string &to, std::string &report);

// ---- What the last run wrote, as data (the Plot tab; docs/upgrade/qt_gui.md section 5) ----

/// One point of a curve: one row of an AZUREOut_* file (src/DataOutputFiles.cpp writes them; c.m. frame, or lab
/// columns when the project's directory holds kinflag.dat with 1).  NaN where the row has no such column.
struct CurvePoint {
  double energy = 0;      ///< E (MeV), c.m. of the entrance pair
  double excitation = 0;  ///< E_x of the compound nucleus (MeV)
  double angle = 0;       ///< theta (deg), c.m.
  double model = 0;       ///< the calculated value in the cross-section column (b, b/sr, or the observable's unit)
  double modelS = 0;      ///< the calculated value in the S-factor column (MeV b)
  double data = 0, dataError = 0, dataS = 0, dataSError = 0;  ///< data sets only (NaN for an extrapolation)
  double band = 0, bandS = 0;  ///< 1 sigma of `model`, `modelS` from the .band file of the run (NaN without one)
};
/// One data set or extrapolation of the project as the last run wrote it.
struct Curve {
  std::string section;   ///< "data" or "extrapolations"
  std::string record;    ///< the record's name
  int key = 0;           ///< its segment key (place among the records of its kind, from 1, inactive ones counted)
  std::string reaction, quantity, unit;  ///< the record's fields as written
  std::string observable;  ///< the engine's observable: "cross_section", "differential", "analyzing_power",
                           ///< "phase_shift", "legendre" (rows of coefficients), "total_capture"
  bool sFactor = false;  ///< the S-factor column means an S factor (an angle-integrated or differential cross section)
  bool positive = true;  ///< the value is a cross section (a log axis applies); false for A_y, phases, coefficients
  std::string file;      ///< the output file, a name in the output directory
  int block = 0;         ///< its block in that file, from 0
  bool haveBand = false; ///< the file's .band sibling holds this block
  std::vector<CurvePoint> points;  ///< empty for "legendre"
  std::vector<std::vector<double>> coefficients;  ///< "legendre": E then the coefficients, one row per point
};
struct Curves {
  std::string outputDirectory;
  bool labFrame = false;            ///< kinflag.dat asked for lab columns
  std::vector<Curve> curves;        ///< data sets first, then extrapolations, each in record order
  std::vector<std::string> notes;   ///< what could not be read or matched, in words (a missing file: not run yet)
};
/// The curves of the project's last run (the files of its output directory as they stand), matched to its records
/// by the engine's own build of the project (each active data set's and extrapolation's pair keys and observable:
/// the file a segment is written to and its block among the segments of that file).  A data set the run left out
/// (no points) is matched through chiSquared.out's keys; a block whose rows differ in number from chiSquared.out's N
/// is refused with a note (the file is older than the project).
Curves ReadCurves(const Project &project);

// ---- The Overview page: the project and its last run at a glance (docs/upgrade/qt_gui.md section 15) ----

/// The levels of one J^pi group of a compound nucleus (the engine's grouping, experiment::LevelTablesOfRows).
struct LevelGroupCount {
  std::string nucleus;  ///< "16O"
  int twoJ = 0;
  int parity = 1;
  std::string label;    ///< "3/2-"
  int levels = 0;
};
/// One option of the calculation (a row of io::ConfigKeys() that is an option) as the project sets it.
struct OptionState {
  std::string name, title, effect;  ///< effect: "physics", "convention", "numerics"
  bool on = false;                  ///< the value the run takes from the project (the key, else the default)
  bool defaultOn = false;
  bool written = false;             ///< the project's <config> writes the key
};
/// The records of one kind of a section, and how many of them an active data set names.
struct RecordUse {
  std::string section, kind;
  int records = 0;
  int used = 0;
};
/// What the project holds, as the Overview page's first card shows it.
struct ProjectSummary {
  std::string file;
  bool built = false;  ///< the engine's build of the project succeeded (else the counts the text gives, no activity)
  std::vector<std::string> nuclei;           ///< the compound nuclei of the levels, in pair order
  std::vector<PairChoice> pairs;             ///< as Pairs() gives them
  std::vector<LevelGroupCount> levelGroups;  ///< in the engine's group order
  int levels = 0;                            ///< levels of the active <levels> rows
  int dataSets = 0, activeDataSets = 0;
  int extrapolations = 0, activeExtrapolations = 0;
  std::vector<RecordUse> effects;            ///< the records of <targets> by kind
  std::string basis;   ///< "Brune", "formal (Lane-Thomas)" or "Park" (use_brune, use_park)
  std::string matrix;  ///< "A-matrix" or "R-matrix" (the positional amatrix)
  std::vector<int> externalCapturePairs;     ///< pairs with external capture on an active row (empty: off)
  std::vector<OptionState> options;          ///< every option row, in the table's order
  std::vector<std::string> notes;            ///< what could not be read, in words
};
/// The project as it stands (its text, not the saved file): pairs, levels per J^pi, records, options.  No set-up of
/// the levels: milliseconds.
ProjectSummary SummarizeProject(const Project &project);

/// One data set of the last run (a line of chiSquared.out), matched to its record by key.
struct DataSetResult {
  int key = 0;
  std::string record, reaction, quantity;  ///< "" when the key matches no record of the project as it stands
  int points = 0;
  double chiSquared = 0;
  double chiSquaredPerPoint = 0;  ///< chiSquared / points (NaN for no point)
  double norm = 1;                ///< the norm the run used
  double normChiSquared = 0;      ///< the norm's prior chi^2 as the file gives it (0 without a prior)
  bool normPrior = false;         ///< the record's norm has a prior (the build's normError > 0)
  double normPull = 0;            ///< sign(norm - centre) sqrt(normChiSquared); NaN without a prior
};
/// One output file and its modification time.
struct OutputFileTime {
  std::string name;
  long long time = 0;  ///< seconds since 1970-01-01 UTC
};
/// The last run as its files in the project's output directory tell it (nothing recomputed).
struct RunSummary {
  std::string outputDirectory;
  bool present = false;      ///< a run's outputs are there (chiSquared.out, or an extrapolation's files)
  bool haveChiSquared = false;  ///< chiSquared.out was read
  /// "fit", "calculate" or "extrapolate", inferred from the files (docs/upgrade/qt_gui.md section 15); "" none.
  std::string command;
  long long runTime = 0;      ///< chiSquared.out's time (an extrapolation: its newest .extrap), seconds since 1970
  long long projectTime = 0;  ///< the project file's time on disk (0 when it is not on disk)
  bool stale = false;         ///< projectNewer, or an input file is newer (changedInputs)
  bool projectNewer = false;  ///< the project file is newer than the run's outputs
  // fitStatus.out (a fit's; empty for a calculation)
  std::vector<std::pair<std::string, std::string>> fitStatus;  ///< its key / value lines, in order
  std::string status;     ///< "converged", "call_limit", "stopped", ... (the engine's FitOutcome words)
  std::string minimizer;  ///< the fitter that ran
  long calls = 0;
  // chiSquared.out
  double totalChiSquared = 0;  ///< Total-Chi-Squared: the data's chi^2
  int points = 0;              ///< Total-N
  std::vector<std::pair<std::string, double>> totals;  ///< every field of the total line, in the file's order
  std::vector<DataSetResult> dataSets;                 ///< in the file's order
  /// The active data sets of the project that the run did not use (no point in their window), by record name.
  std::vector<std::string> leftOut;
  /// The records whose input file (InputFiles) is newer than the run's outputs, by name; each makes the run stale.
  std::vector<std::string> changedInputs;
  std::string parameterFile;  ///< the run's parameter file: param.sav after a fit, else param.par ("" none)
  std::vector<OutputFileTime> files;  ///< the files of the output directory, by name
  std::vector<std::string> notes;
};
/// The files of the project's output directory as they stand.  Milliseconds.
RunSummary ReadLastRun(const Project &project);

/// A file a record of the project reads (the engine's build of the project: a data set's data and covariance files,
/// an extrapolation's grid file), resolved against the project's directory as the engine opens it.
struct InputFile {
  std::string path;
  std::string section, record;
  std::string role;  ///< "data", "covariance", "grid"
};
/// The input files of the active records (none when the project does not build).  Milliseconds.
std::vector<InputFile> InputFiles(const Project &project);

/// A parameter of the run's parameter file that ended at a limit of the fit or, after a fit, without an error.
struct ParameterFlag {
  std::string name, kind, unit;
  double value = 0, error = 0, lower = 0, upper = 0;
  bool atLimit = false;
  bool noError = false;
};
/// The parameters of the saved project as the engine sets them up, against the run's parameter file.
struct RunParameters {
  int parameters = 0;
  int free = 0;
  std::vector<std::pair<std::string, int>> freeByKind;  ///< kind key (ParamKindTable) -> free ones, in vector order
  std::string file;                    ///< the parameter file read ("" none)
  std::vector<ParameterFlag> flagged;  ///< at a limit or without an error, in vector order
};
/// ReadParameters on `projectFile`, then `parameterFile` ("" none) read by name: a free parameter is at a limit when
/// |value - limit| <= 1e-6 max(1, |limit|) (fit units), without an error when `fit` and its error is 0 or missing.
/// Sets up a session (seconds for a project with external capture): call it off a window's thread.
bool ReadRunParameters(const std::string &projectFile, const std::string &parameterFile, bool fit,
                       RunParameters &out, std::string &why);
/// chi^2 / (N - p) with chi^2 and N of the run and p its free parameters; NaN when N <= p.
double ChiSquaredPerDegreeOfFreedom(const RunSummary &run, const RunParameters &parameters);

// ---- Fit parameters (the Fitting Settings tab) ----

/// One parameter of the fit as the engine sets it up (a session of the project, as pyazr's Session): the vector of
/// param.sav, in its order.
struct Parameter {
  std::string name;       ///< as param.sav names it (energy_1, width_1_2, el_lab090.norm, ...)
  std::string kind;       ///< the kind's key (ParamKindTable: "energy", "width", "norm", "shift", ...)
  std::string unit;       ///< the unit of `value` (the fit's: MeV^1/2 for a reduced-width amplitude)
  double value = 0;       ///< the start value (fit units, as param.sav writes it)
  bool fixed = false;
  double lower = 0, upper = 0;  ///< the fit's limits in fit units (-inf / +inf: none), <parameterSettings> applied
  int levelNumber = -1;   ///< level counted over all J groups (level parameters), else -1
  int channel = -1;       ///< channel within the J group (widths), else -1
};
/// The parameters of a fit of the SAVED project file (the levels' and the data sets' nuisances): the session's set-up (levels, data, external capture), nothing written but what a session
/// writes (an intEC.dat it computes).  False with `why` when the project cannot be set up.
bool ReadParameters(const std::string &projectFile, std::vector<Parameter> &out, std::string &why);
/// The same with what the set-up writes (an intEC.dat it computes) written into `scratchDirectory` (an existing
/// directory) instead of the project's output directory, which is left as it is.
bool ReadParameters(const std::string &projectFile, const std::string &scratchDirectory, std::vector<Parameter> &out,
                    std::string &why);

/// One line of a parameter file (param.sav, param.par, param.errors' first three columns): name, value, error.
struct ParameterValue {
  std::string name;
  double value = 0, error = 0;
};
/// The lines of a parameter file as the engine reads it (AZUREParams: "name value error", "#parametrization" first);
/// false with `why` when it cannot be read.
bool ReadParameterFile(const std::string &path, std::vector<ParameterValue> &out, std::string &why);

/// One row of <parameterSettings> as the engine reads it (ParameterLimitsManager::ReadParameterSettings): the one-token
/// form "name nominal lower upper error fitError nuisance category index" or the GUI's "Level N Energy (MeV) ..." and
/// "Level N Channel M Width (eV) ..." forms.  lower == upper == 0: no limits.  A width's numbers are partial widths
/// (eV), converted by the engine.
struct ParameterSetting {
  std::string name;
  double nominal = 0, lower = 0, upper = 0, error = 0, fitError = 0;
  bool nuisance = false;
  std::string category = "level";
  int index = 0;
  int line = 0;  ///< 1-based line in the file (0 for a row not yet written)
};
std::vector<ParameterSetting> ParameterSettings(const Project &project);
/// Replace (by name) or append one row, written in the one-token form when `name` is one token; `remove` deletes it.
/// False with `why` for a name with a '#' or an empty one.
bool SetParameterSetting(Project &project, const ParameterSetting &row, std::string &why, bool remove = false);

/// The sampler's prior rows (<mcmc> <parameters>, src/SamplerRun.cpp): name (as param.sav or "label (name)"), value,
/// prior mean, prior sigma, use.
struct SamplerPrior {
  std::string name;
  double value = 0, mean = 0, sigma = 0;
  bool use = false;
};
std::vector<SamplerPrior> SamplerPriors(const Project &project);
/// Write the <parameters> block of <mcmc> with `rows` (the section's other lines kept; the section made when the file
/// has none and `rows` is not empty).
void SetSamplerPriors(Project &project, const std::vector<SamplerPrior> &rows);

// ---- <thm>: the engine's tables of its keys and its edits of the section's lines (the THM tab) ----

/// The sampler's chain as the last mcmc run wrote it (samples.mcmc in the project's output directory), read by the
/// engine's one reader of that file (ReadSampleFile, api/include/AZUREAPISampler.h).  Every array is row-major:
/// chain [step][walker][parameter] in the order of the run's free parameters and in its basis (rwa: the fit vector;
/// physical: the project's parameters); the log densities [step][walker].
struct SamplerChain {
  std::string file;  ///< the file read
  int steps = 0, walkers = 0, parameters = 0;
  std::vector<double> chain;
  std::vector<double> logProbability, logLikelihood, logPrior;
};
/// False with `why` when the project's output directory has no samples.mcmc or the file cannot be read.
bool ReadSamplerChain(const Project &project, SamplerChain &out, std::string &why);

/// The keys of <thm> as two record kinds in the schema's shape (ThmOptionKeys, ThmExperimentKeys, the tables the
/// engine's parser reads): kind "<thm>" (the global options) and kind "experiment" (the keys of an experiment[<name>]
/// line).  FieldInfo::appliesWhen names the key that must be given for one to apply ("vertexModel=ff" ...).
std::vector<KindInfo> ThmKinds();
/// One key=value of <thm> as the reader takes it (ThmBlockEntries): a global option (experiment "") or a key of an
/// experiment; `line` 1-based in the file.
struct ThmEntry {
  std::string experiment, key, value;
  int line = 0;
};
std::vector<ThmEntry> ThmEntries(const Project &project);
/// Set (or remove) a global option of <thm> (ThmOptionEdits; the section is made when the file has none).
bool SetThmOption(Project &project, const std::string &key, const std::string &value, std::string &why,
                  bool remove = false);
/// Set (or remove) a key of experiment[name] (ThmExperimentEdits); with `remove` and key "" the experiment's lines.
bool SetThmExperimentKey(Project &project, const std::string &name, const std::string &key, const std::string &value,
                         std::string &why, bool remove = false);

/// Whether this build's program has `command` (Commands(): a command gone from the engine is gone here).
bool HasCommand(const std::string &command);

/// The single-particle (Wigner-limit) width Gamma_sp = 2 P_l(E) gamma_W^2 of a particle channel, in MeV: charges
/// z1, z2, reduced mass mu (u), orbital momentum l, channel radius a (fm), c.m. energy E > 0 (MeV).  The engine's
/// function (physics::single_particle_width), not a copy; NaN where it is not defined (E <= 0).
double SingleParticleWidth(int z1, int z2, double mu_u, int l, double a_fm, double E_MeV);

// ---- The Levels and Pairs tabs (docs/upgrade/qt_gui.md section 11) ----

/// The columns of a <levels> row, in the order of LevelRow::tokens (io::LevelColumnNames: "levelJ", "levelPi",
/// "levelE", ..., "e1"), and the index of one by name (-1 for none), so that a front end names a column, never a place.
const std::vector<std::string> &LevelColumnNames();
int LevelColumnIndex(const std::string &name);

/// What the width column of a <levels> row is and the channel's quantities at the row's level energy E_x (MeV).
/// NaN where a quantity is not defined for the channel (a photon channel has no penetrability of its own here).
struct ChannelQuantities {
  /// "open" (particle channel, E_c > 0: Gamma in eV), "closed" (E_c < 0: an ANC in fm^-1/2), "threshold" (E_c = 0:
  /// no width defined), "amplitude" (gammaIsRWA: a reduced-width amplitude in MeV^1/2), "photon" (Gamma in eV),
  /// "magnetic_moment" (a static M1: mu in mu_N), "quadrupole_moment" (a static E2: Q in b), "beta" (B).
  std::string kind;
  std::string symbol;   ///< UTF-8: "Γ", "C", "γ", "μ", "Q", "B"
  std::string unit;     ///< "eV", "fm^-1/2", "MeV^1/2", "μ_N", "b", ""
  std::string meaning;  ///< in words: "partial width of the open channel", ...
  bool particle = false;      ///< a particle pair's channel (pType 0)
  double value = 0;           ///< the width column as read
  double channelEnergy = 0;   ///< E_c = E_x - S - e (MeV, c.m.), as the engine forms it at a level
  double radius = 0;          ///< the radius used (fm)
  double P = 0, S = 0, dSdE = 0;  ///< penetrability (or mu a / (hbar^2 W^2) below threshold, MeV^-1 fm^-1), shift, dS/dE
  double wignerLimit = 0;         ///< gamma_W^2 = 3 hbar^2 / (2 mu a^2) (MeV), physics::wigner_limit
  double singleParticleWidth = 0; ///< Gamma_sp = 2 P gamma_W^2 (MeV), open channels only
  double reducedWidth = 0;        ///< the observed reduced width gamma^2 = Gamma / (2 P) or C^2 / (2 P) (MeV)
  double wignerRatio = 0;         ///< theta^2 = gamma^2 / gamma_W^2
  std::string error;              ///< why the row could not be read ("" when it was)
};
/// The row's channel at its level energy: P, S, dS/dE from the engine's ParticleChannelValues on the pair the engine
/// builds from the row (PPair(LevelRow)), as the binding's level doors.  radius_fm > 0 replaces the row's radius (a
/// channel's own radius, ChannelRadius); relativistic: the project's <config> relativistic.
ChannelQuantities ChannelAtLevel(const LevelRow &row, double radius_fm = 0, bool relativistic = false);

/// The radius of channel (l, 2s) of a J group 2J of pair `pairKey` as the set-up gives it: the most specific record of
/// <channelRadii> that names it (experiment::ChannelRadiusFor on the engine's build of the records), else the pair's.
struct ChannelRadiusInfo {
  double radius = 0;   ///< fm
  bool own = false;    ///< a record gives it
  std::string record;  ///< that record's name
  bool ambiguous = false;  ///< two records name it and neither is the more specific (the set-up refuses)
};
ChannelRadiusInfo ChannelRadius(const Project &project, int pairKey, int l, int twoS, int twoJ, double pairRadius);

/// What names a pair besides its <levels> rows: the data sets and extrapolations whose reaction resolves to it (the
/// engine's build of the project: entrance and exit keys), as "data: <record>" / "extrapolations: <record>".
std::vector<std::string> PairUses(const Project &project, int pairKey);

/// The engine's version string ("2.1.0").
std::string Version();

// ---- Physics (the Physics tab; docs/upgrade/physics_tab.md section 7 is this interface) ----
//
// Two kinds of call.  Session-free: the rows of the physics table (api/physics_table, the rows pyazr.physics and the
// web page's Physics tab read) and one row evaluated element-wise.  Per project: a PhysicsSession, the engine's
// session of a saved project at the start values or a parameter file (as pyazr's Session(file, parameter_file=)),
// whose calls give the session's channel functions, a J^pi block's matrices on an energy grid, the levels' physical
// energies and widths, and the THM diagnostics.  Units: energies MeV, radii fm, angles deg, phases rad unless a field
// says deg.  Not reentrant on one PhysicsSession (a caller runs its calls in one thread at a time); two sessions are
// independent.

/// One argument of a physics function (its row: api/physics_table/include/azure2/physics_table/PhysicsTable.h).
struct PhysicsArgument {
  std::string name, unit, doc;
  std::string kind;    ///< "real", "integer", "flag" (0 or 1), "choice" (the index of one of `choices`)
  std::string domain;  ///< in words: "finite > 0", "integer in [0, 10000]", ...
  double lower = 0, upper = 0;        ///< the domain's ends (-inf / +inf: unbounded)
  bool lowerOpen = false, upperOpen = false;
  std::vector<std::string> choices;   ///< "choice": the names, in index order
  bool hasDefault = false;
  double defaultValue = 0;            ///< a choice: the index; a flag: 0 or 1
  bool hasExample = false;
  double example = 0;                 ///< a typical value: a form's start value where there is no default
  std::string exampleText;            ///< "text" (a stopping row's formula): a typical value ("Ta2O5")
  bool hasSweep = false;              ///< the argument may be a plot's axis; the range of a first plot over it:
  double sweepLower = 0, sweepUpper = 0;
  bool sweepLog = false;              ///< log spacing (else linear)
};
struct PhysicsReturn {
  std::string name, unit, doc;
};
/// One function of the physics table (pyazr.kinds()["physics"] gives the same rows).
struct PhysicsFunction {
  std::string name, summary, doc;
  std::string library;  ///< the library function it calls (header: function)
  std::string axis;     ///< the argument a plot puts on its horizontal axis
  std::string group;    ///< the table's section ("kinematics", "Coulomb functions", ...)
  std::vector<PhysicsArgument> arguments;
  std::vector<PhysicsReturn> returns;
};
/// The rows, in the table's order.  Fast (no calculation).
std::vector<PhysicsFunction> PhysicsFunctions();
/// What one row gives on n elements: values[k][i] is return k at element i (NaN where status[i] != 0).  status[i]:
/// 0 computed; k > 0 argument k (1-based) outside its domain; -2 the arguments together outside the domain; -1 not
/// converged or failed; a row's own positive code where its doc says so.
struct PhysicsValues {
  std::vector<std::vector<double>> values;
  std::vector<int> status;
};
/// Row `name` on the arguments: one column per argument in the row's order (a choice as its index, a flag as 0 or 1),
/// each of length n or 1 (a length-1 column is repeated: broadcasting).  False with `why` for an unknown name or
/// columns of other lengths; a value outside a domain is a status, not a refusal.  The row's arithmetic is the one
/// pyazr.physics.<name> runs, bit for bit.  About 1-50 microseconds per element (coulomb_waves, the THM form factor
/// about 1 ms).
bool EvaluatePhysics(const std::string &name, const std::vector<std::vector<double>> &arguments, PhysicsValues &out,
                     std::string &why);
/// The words of a status of row `name` ("e_cm outside its domain (finite > 0)", "not converged or failed", ...).
std::string PhysicsStatusText(const std::string &name, int status);

/// The stopping table (stopping/: SRIM-2013 stopping, energy loss, range and straggling; the rows of
/// pyazr.kinds()["stopping"] and the web page's stopping functions): rows of the physics table's shape whose first
/// argument is the material's formula (kind "text", one per call), then the numeric ones (phase, compound_correction,
/// z1, mass_u, then the row's own).  Energies lab MeV of the projectile; stopping in eV/(1e15 atoms/cm^2) per atom of
/// the material; areal densities atoms/cm^2; straggling MeV.  Every use cites SRIM (StoppingCitation()).
std::vector<PhysicsFunction> StoppingFunctions();
std::string StoppingCitation();
/// The words of a status of stopping row `name` (k > 0: the k-th NUMERIC argument outside its domain).
std::string StoppingStatusText(const std::string &name, int status);
/// Row `name` for the material `formula`: one column per NUMERIC argument (the arguments after the formula), each of
/// length n or 1; as EvaluatePhysics (status -2 where a layer is at least the range).  False with `why` for an unknown
/// row, columns of unequal lengths or a formula stopping/ refuses.  Bit for bit pyazr.stopping.evaluate.
bool EvaluateStopping(const std::string &name, const std::string &formula,
                      const std::vector<std::vector<double>> &arguments, PhysicsValues &out, std::string &why);

/// One series of a per-session answer: a field over the answer's axis, for one line (a channel, a pair of channels,
/// a level).  `keys` label the line with numbers ("jgroup", "channel", "channel_out", "level", "l", "s", "J", "pi",
/// "pair"); `label` is the same in words.
struct PhysicsSeries {
  std::string name;   ///< the field: "P", "S", "delta", "U_abs", ...
  std::string unit;
  std::string label;
  std::vector<std::pair<std::string, double>> keys;
  std::vector<double> values;  ///< one per axis value; NaN where not defined there
  double Key(const std::string &key, double absent = -1) const;
};
/// A per-session answer drawn against one axis: the axis values, the series, scalar facts ("E_aA" ...), notes.
struct PhysicsCurves {
  std::string axis, axisUnit;
  std::vector<double> axisValues;
  std::vector<PhysicsSeries> series;
  std::vector<std::pair<std::string, std::string>> facts;  ///< name, value in words (with its unit)
  std::vector<std::string> notes;                           ///< what was left out, and why
  /// The first series named `name` whose keys include every (key, value) of `match`, or nullptr.
  const PhysicsSeries *Find(const std::string &name,
                            const std::vector<std::pair<std::string, double>> &match = {}) const;
};

/// One particle or photon pair of the session (<levels>' pair columns, as the engine set them up).
struct PhysicsPair {
  int key = 0;            ///< the pair key of <levels>
  std::string kind;       ///< "particle", "photon", "beta"
  int z1 = 0, z2 = 0;
  double m1 = 0, m2 = 0;  ///< atomic masses (u)
  double j1 = 0, j2 = 0;
  double separation = 0, excitation = 0;  ///< MeV: the threshold is E_x = separation + excitation
  double radius = 0;      ///< the pair's channel radius (fm)
  bool identical = false;
};
/// One channel of a J^pi group.
struct PhysicsChannel {
  int channel = 0;  ///< 1-based in the group
  int pair = 0;     ///< pair key
  std::string radiation;  ///< "P" particle, "E" or "M" photon (multipolarity l)
  int l = 0;
  double s = 0;
  double radius = 0;  ///< a_c (fm), 0 for a photon channel
};
/// One J^pi group (the J groups of parameters.out, 1-based).
struct PhysicsGroup {
  int jgroup = 0;
  double J = 0;
  int parity = 1;
  int levels = 0;     ///< levels in the R-matrix
  std::vector<PhysicsChannel> channels;
};
/// One level's physical (observed) energy and widths: the engine's transformation (CompoundSystem::TransformOut, the
/// numbers parameters.out writes), at the session's vector.
struct PhysicsPole {
  int jgroup = 0, level = 0;  ///< 1-based, as parameters.out
  double J = 0;
  int parity = 1;
  double energy = 0;          ///< E_x of the level (MeV), the physical energy in every parametrization
  double formalEnergy = 0;    ///< E_lambda of the vector (MeV): the formal pole energy in the standard basis
  double totalWidth = 0;      ///< the sum of the observed widths of the open particle channels and the photon channels (MeV)
  struct Width {
    int channel = 0, pair = 0;
    std::string radiation;
    int l = 0;
    double s = 0;
    std::string kind;   ///< "width" (an open channel, W in MeV) or "anc" (a closed one: `value` in fm^-1/2)
    double width = 0;   ///< observed width W (MeV); 0 for a closed channel
    double value = 0;   ///< the reported value (the width in eV, the ANC in fm^-1/2, as parameters.out)
    double penetrability = 0;  ///< the P the width was formed with (0 closed)
  };
  std::vector<Width> widths;
};

struct PotentialPole;  // the Potential panel's (below)

/// The engine's session of a saved project (docs/upgrade/physics_tab.md section 7.2).  Opening one sets up the
/// levels and the data (or the extrapolations) as AZURE2 calculate would: seconds for a project with external capture
/// (its integrals are read from intEC.dat when current).  Every call after that takes milliseconds per energy.
class PhysicsSession {
 public:
  PhysicsSession();
  ~PhysicsSession();
  PhysicsSession(PhysicsSession &&) noexcept;
  PhysicsSession &operator=(PhysicsSession &&) noexcept;
  PhysicsSession(const PhysicsSession &) = delete;
  PhysicsSession &operator=(const PhysicsSession &) = delete;

  /// The session of the saved project file at its start values, or at the values of `parameterFile` ("" for none),
  /// read as AZURE2 --parameter-file reads it (AZUREParams::ReadUserParameters, every name it gives, fixed ones
  /// included; a file of another parametrization converted or refused as the program does).  `extrapolations`: set
  /// up the extrapolations instead of the data (a project with no data set needs it).  False with `why` (the
  /// engine's words) when the project cannot be set up.
  static bool Open(const std::string &projectFile, const std::string &parameterFile, bool extrapolations,
                   PhysicsSession &out, std::string &why);
  bool IsOpen() const;
  /// "standard", "brune" or "park": the parametrization the levels' numbers are in.
  std::string Parametrization() const;

  std::vector<PhysicsPair> Pairs() const;
  std::vector<PhysicsGroup> Groups() const;

  /// The levels' physical energies and widths (the poles the engine knows on the real axis).
  std::vector<PhysicsPole> Poles();

  /// The session's channel functions of pair `pair` (key), orbital `l`, at c.m. energies `eCm` of the pair (MeV):
  /// series F, Fp, G, Gp, P, S, delta_hs (rad, the hard-sphere phase shift) and open (1 or 0), axis "e_cm".  `radius`
  /// fm, 0 for the channel radius; `twoJ` the channel's 2J for a pair with a spin-orbit term (else -1).  As pyazr's
  /// Session.channel_functions (the Coulomb table when one is in use; the hybrid model's solution for a pair with a
  /// nuclear potential).  False with `why` for an unknown or photon pair or a J the pair cannot take.
  bool ChannelFunctions(int pair, int l, double radius, int twoJ, const std::vector<double> &eCm, PhysicsCurves &out,
                        std::string &why);

  /// One J^pi group on excitation energies `excitation` (E_x, MeV; axis "E_x"), as pyazr's Session.jpi_block.  Series
  /// (keys: jgroup; channel; channel_out; level; level_out):
  ///  per channel   E_c (MeV), P, S, L_re and L_im (S - B_c and P: the level matrix's L_c - B_c);
  ///  per pair of open particle channels c, c' (c <= c'):  U_re, U_im, U_abs (|U_cc'|), T_re, T_im, T_abs (|T_cc'|), and on the diagonal
  ///                delta (deg, the nuclear phase shift arg(U_cc)/2 as the engine's phase_shift data kind: in
  ///                (-90, 90], [0, 180) for identical particles), eta (the inelasticity |U_cc|);
  ///  per pair of levels l <= m:  A_abs (|A_lm|, MeV^-1);
  ///  once:          det_re, det_im (MeV^n), log10_abs_det (log10 |det A^-1|: minima at the resonances).
  /// NaN where an element is not defined (a closed channel).  `parameterValues` empty: the session's vector.
  /// About 0.1 ms per energy for a few levels; the 2-D maps of section 5 are these series as rows.
  bool Block(int jgroup, const std::vector<double> &excitation, PhysicsCurves &out, std::string &why);

  /// The THM experiments of the project (<thm> experiment[<name>]), by name.
  std::vector<std::string> ThmExperiments() const;
  /// THM diagnostics at the session's vector on c.m. energies of x + A (MeV, axis "E"), as pyazr's
  /// Session.thm_vertex, thm_lineshape and thm_distortion (session_thm.py; the same engine tables):
  ///  vertex:     M2, M2_qf (and M2_pw for vertexModel=dw) per channel and level, reached (1 or 0);
  ///  lineshape:  NC2 per exit pair and level, zeta and eta_sb per exit pair, E_sF, eta_0;
  ///  distortion: R, R_model, M2, M2_PW, E_sF, k_sF, eta_sF, theta_cm, x, q, lmax.
  /// False with `why` (the engine's words) for an experiment that is not there or has no such part.
  /// The box states of pair `pair`'s potential in (l, 2J) through the session (AZUREAPI::GetPotentialBoxStates, pyazr's
  /// Session.potential_box_states): `energy`, `gamma2` and `nodes` of each PotentialPole (the rest left 0).  False with
  /// the engine's words.  The Potential panel computes them from the record (PotentialWaveOf); this is the session's.
  bool PotentialBoxStates(int pair, int l, int twoJ, double radius, double boundary, double eMin, double eMax,
                          std::vector<PotentialPole> &out, std::string &why);
  /// The potential's own a and r in (l, 2J) through the session (AZUREAPI::GetPotentialEffectiveRange, pyazr's
  /// Session.potential_scattering_length); false with its reason.
  bool PotentialScatteringLength(int pair, int l, int twoJ, double &a, double &r, std::string &why);

  bool ThmVertex(const std::string &experiment, const std::vector<double> &energies, PhysicsCurves &out,
                 std::string &why);
  bool ThmLineshape(const std::string &experiment, const std::vector<double> &energies, PhysicsCurves &out,
                    std::string &why);
  bool ThmDistortion(const std::string &experiment, const std::vector<double> &energies, PhysicsCurves &out,
                     std::string &why);

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

// ---- What the Physics tab draws from the per-session calls (docs/upgrade/qt_gui.md section 15.2) ----

/// One per-session call of PhysicsSession that answers with a PhysicsCurves.
struct PhysicsCall {
  std::string name;   ///< "ChannelFunctions", "Block", "ThmVertex", "ThmLineshape", "ThmDistortion"
  std::string title;  ///< in words, for a drop-down
  std::string doc;    ///< what it is, one or two sentences
  std::string axisTitle;  ///< the horizontal axis in words with its unit ("E_x of the compound nucleus (MeV)")
  bool slow = false;  ///< milliseconds or more per energy (a DWBA per energy) or a chi^2 evaluation first
};
std::vector<PhysicsCall> PhysicsCalls();

/// One thing a page can draw from a call's answer: the series named `y` (every series of that name, one line each,
/// labelled by its `label`) against the call's axis, or against the series named `x` with the axis as the curve's
/// parameter (the Argand plot).
struct PhysicsPlot {
  std::string call;   ///< the PhysicsCall it reads
  std::string name;   ///< an identifier ("delta", "argand", "map_det", ...)
  std::string title, doc;
  std::string y, yTitle;   ///< the series on the vertical axis; its axis title with unit
  std::string x, xTitle;   ///< "" the call's axis; else the series on the horizontal axis
  bool diagonalOnly = false;  ///< only the series whose channel equals channel_out
  bool logY = false;          ///< a log axis suits it (a value over orders of magnitude, > 0)
  /// A map (Pole map page): "" not a map; "series" rows = the series named `y` of one answer (one J group);
  /// "jgroup" rows = one answer per J group, the series named `y` of each.
  std::string mapRows;
};
/// The plots of one call (`call` = "" for every call's), in the order a drop-down lists them.
std::vector<PhysicsPlot> PhysicsPlots(const std::string &call);

/// The parameter file of the project's last fit, `param.sav` in its output directory (as the program writes it),
/// when that file exists; "" otherwise.
std::string LastParameterFile(const Project &project);

// ---- Round 8: what the quantity decides, kernels drawn, sums of reactions, the channel rule (qt_gui.md section 15) ----

/// One field a quantity uses (experiment::Quantity::fields): the fields of a record kind marked FieldInfo::byQuantity
/// that the quantity does not list are refused for it.
struct QuantityFieldRule {
  std::string field;
  bool required = false;
  bool conditional = false;  ///< the rule holds only when `when` holds
  std::string when;          ///< the condition in words ("every reaction term has a photon exit"), "" for always
};
/// The fields quantity `quantity` uses; false for a quantity (or alias) the table has not.
bool QuantityFieldRules(const std::string &quantity, std::vector<QuantityFieldRule> &out);

/// The shape of one kernel at an energy: x lab MeV of the beam, y the kernel's weight per lab MeV (area 1).
struct KernelCurve {
  std::string kernel;  ///< the kernel table's row: "doppler", "gaussian", "rectangular", "triangular", "beam_profile"
  std::string title;   ///< in words, with its parameter ("free-gas Doppler broadening, T = 300 K")
  std::string xLabel;
  std::vector<double> x, y;
  double fwhm_keV = 0;    ///< full width at half maximum, keV lab (NaN when the curve has no half maximum on both sides)
  double fwhmCm_keV = 0;  ///< the same in the c.m. frame of the pair
};
/// What a record of <targets> does to the points of a data set that names it, drawn from the engine's kernels.
struct KernelPreview {
  std::vector<std::string> users;  ///< the data sets and extrapolations that name the record ("data: s1", ...)
  std::string user;                ///< the one drawn for ("" when none names it: the first particle pair is used)
  std::string pair;                ///< its entrance pair ("p+15N")
  double cmFactor = 0;             ///< E_c.m. / E_lab of the pair (experiment::BeamCMFactor)
  double labEnergy = 0, cmEnergy = 0;  ///< E_0 (MeV)
  double lowEnergy = 0, highEnergy = 0;  ///< the user's data or grid, lab MeV (0, 0 unknown)
  std::vector<KernelCurve> curves;     ///< one per kernel the record selects
  std::vector<std::string> notes;      ///< what was assumed or left out, in words
};
/// The kernels a record of section `section` (a target, beam or detector, as a form holds it: applied or not) selects,
/// for the points of `user` ("data: <name>" or "extrapolations: <name>"; "" the first that names it) at lab energy
/// `labEnergy` (MeV; NaN: the middle of the user's energies), each on `points` nodes.  The record is set into a copy
/// of `project`, built by the engine's reader, and each kernel made by the engine's TargetEffect from the built spec
/// (docs/upgrade/physics_tab.md section 6.2).  False with `why` when the copy does not build or the record selects
/// no kernel.
bool KernelPreviewOf(const Project &project, const std::string &section, const Record &record, const std::string &user,
                     double labEnergy, int points, KernelPreview &out, std::string &why);

/// The columns of a data file of `quantity` with `columns` columns: DefaultColumns(quantity), then columns the engine
/// skips (the role "_") up to that number, or the first `columns` of them.
std::vector<std::pair<std::string, std::string>> ColumnsForFile(const std::string &quantity, int columns);

/// One term of a reaction expression (section 3.5 of the format): [scale*]reaction [@ angle unit frame].
struct ReactionTermInfo {
  std::string scale;     ///< the weight as written ("" for none: 1)
  std::string reaction;  ///< the reaction as written ("20Ne(n,a)17O")
  int entrance = 0;      ///< its pair keys (0: it names no single pair of the project); exit -1 a total capture
  int exit = 0;
  std::string angle, angleUnit, angleFrame;  ///< a fixed angle ("" for none)
};
struct ReactionExpression {
  std::string op = "single";  ///< "single", "sum", "ratio"
  std::vector<ReactionTermInfo> terms;
};
/// The reaction expression `text` split by the engine's parser (experiment::ParseReaction); false with the parser's
/// words otherwise.
bool ParseReactionExpression(const Project &project, const std::string &text, ReactionExpression &out, std::string &why);
/// The text of an expression: its terms as the grammar writes them, joined by " + " or " / ".
std::string ReactionExpressionText(const ReactionExpression &expression);

/// One channel a level may have: its pair key, l (a photon channel: the multipolarity; a beta channel: 1
/// Gamow-Teller, 0 Fermi) and 2s.
struct AllowedChannelInfo {
  int pairKey = 0;
  int l = 0;
  int twoS = 0;
};
/// The channels a level of 2J `twoJ` and parity `parity` has with the pairs of `pairRows` (one <levels> row per
/// pair, the GUI's unsaved pairs; keys as the rows give them), limited to `pairKeys` (empty: all): the engine's rule
/// (experiment::AllowedChannels, pyazr's) with the limits lMax, multipolarityMax and multipolaritiesMax, less the
/// channels the engine's set-up refuses: l + s odd for identical particles, a photon channel between an integer and a
/// half-integer spin.  False with the engine's words for a refused J^pi, limit, key or row.
bool AllowedChannels(const std::vector<LevelRow> &pairRows, int twoJ, int parity, int lMax, int multipolarityMax,
                     int multipolaritiesMax, const std::vector<int> &pairKeys, std::vector<AllowedChannelInfo> &out,
                     std::string &why);

// ---- The Potential panel of the Pairs tab and the Levels tab's partial waves (docs/upgrade/qt_gui.md section 19) ----

/// One `wave` field of a `potential` record in its meaning, split by the engine's parser (ParseFieldValue); every
/// number and spin keeps the text it was written with, so PotentialWaveText of an unchanged wave is the text read.
struct PotentialWave {
  int l = 0;
  std::string J;            ///< "" every J of l; else the spin as written ("1/2")
  int twoJ = -1;            ///< 2J, -1 without J
  std::string depth, width; ///< the wave's own V0 and Gaussian r0 as written ("" the central one's)
  std::string depthUnit = "MeV", widthUnit = "fm";
  bool background = false;  ///< `background interior`: the external + internal mode
  std::string boundary;     ///< B_bg as written ("" none)
  std::string subtract;     ///< the node number of the subtracted box state as written ("" none)
  bool capture = false;     ///< `capture interior`
  bool boundStateOnly = false;  ///< `only bound_state`: point-Coulomb channel functions in this wave
};
/// The value `text` of a `wave` field split by the engine's parser; false with its words.
bool ParsePotentialWave(const std::string &text, PotentialWave &out, std::string &why);
/// The text of a `wave` field, written by the record writer from the parts (FieldValueText).
std::string PotentialWaveText(const PotentialWave &wave);

/// The potential record of a pair: the active one the engine's build gives the pair, else a record switched off whose
/// pair is the pair.  `found` false: the pair has none.
struct PairPotential {
  bool found = false;
  bool on = true;  ///< false: switched off (Project::RecordsOff)
  Record record;
};
PairPotential PotentialOfPair(const Project &project, int pairKey);

/// A box state of a wave's potential: a pole of its interior R-function at the channel radius (physics/
/// InteriorBoxStates, as Session.potential_box_states).
struct PotentialPole {
  int nodes = 0;            ///< what `subtract` names
  double energy = 0;        ///< E_p (MeV, c.m. of the pair)
  double excitation = 0;    ///< E_x = E_p + the pair's threshold (MeV)
  double gamma2 = 0;        ///< gamma_p^2 (MeV)
  double wignerRatio = 0;   ///< gamma_p^2 / gamma_W^2, gamma_W^2 = 3 hbar^2/(2 mu a^2) (physics::wigner_limit)
  bool subtracted = false;  ///< named by the wave's `subtract`
};
/// A level of <levels> in the J^pi of a wave with a channel of the pair and l: what the panel shows beside the poles.
struct PotentialLevel {
  double J = 0;
  int parity = 1;
  double excitation = 0;   ///< E_x (MeV)
  double energy = 0;       ///< E_c = E_x - threshold (MeV)
  bool active = true;
  double s = 0;            ///< channel spin
  std::string kind;        ///< ChannelAtLevel's: "open" (value Gamma, eV), "closed" (ANC, fm^-1/2), "amplitude", ...
  double value = 0;        ///< the channel's value as written
};
/// What one partial wave of a pair's potential is and gives, at the channel radius.
struct PotentialWaveReport {
  int l = 0, twoJ = -1;     ///< twoJ: the 2J computed (a wave without J of a J-dependent potential: its lowest)
  int parity = 1;           ///< of the J^pi the wave forms (pi1 pi2 (-1)^l)
  double radius = 0;        ///< a_c (fm): the pair's, or a <channelRadii> record's
  std::string radiusNote;   ///< where it comes from, or that its channel spins have different ones (the first is used)
  double threshold = 0;     ///< E_x of the pair's threshold (MeV)
  double depth = 0, width = 0;  ///< the wave's V0 (MeV) and Gaussian r0 (fm; 0 for a Woods-Saxon term)
  bool background = false;
  double boundary = 0;
  std::vector<PotentialPole> poles;  ///< E_p in [eMin, eMax], ascending
  std::string subtractError;         ///< the set-up's refusal of `subtract` (SubtractedBoxStates), "" when none
  std::string radiusError;           ///< the set-up's refusal of a background beyond the matching radius, "" when none
  // the potential alone
  int boundStates = 0;      ///< bound states of the potential in the wave (forbidden ones included)
  int zeroEnergyNodes = 0;  ///< nodes of its zero-energy solution inside r_m
  bool effectiveRange = false;  ///< a and r below were computed (else `effectiveRangeWhy`)
  std::string effectiveRangeWhy;
  double a = 0, r = 0, aSpread = 0, rSpread = 0, eLow = 0, eHigh = 0;  ///< physics/ PotentialEffectiveRange
  std::vector<PotentialLevel> levels;  ///< the project's levels of this J^pi with a channel (pair, l)
};
/// The wave (l, twoJ; -1 none) of pair `pairKey` with the potential of `record` (as a form holds it, applied or not;
/// on or off), set into a copy of `project` built by the engine's reader, the pair as the engine makes it from its
/// <levels> row: the box states in [eMin, eMax] (MeV, c.m.) at the channel radius with the wave's boundary B_bg (0
/// without a background), the set-up's checks of the wave, and the potential's own bound states, a and r.  False with
/// `why` when the copy does not build (the engine's words) or the wave cannot be asked (a J the potential needs).
bool PotentialWaveOf(const Project &project, const Record &record, int pairKey, int l, int twoJ, double eMin, double eMax,
                     PotentialWaveReport &out, std::string &why);

/// A depth found for a wave: the tuning's answer.
struct PotentialDepth {
  bool found = false;
  std::string why;          ///< the engine's words when not found
  double depth = 0;         ///< V0 (MeV)
  double a = 0, r = 0;      ///< the wave's a and r at that depth (scattering-length tuning)
  int boundStates = 0;      ///< bound states at that depth (the branch kept)
  int nodes = 0;            ///< nodes of the zero-energy solution (scattering length) or of the bound state (binding)
  double energy = 0, anc = 0;  ///< binding: the eigenvalue (MeV, c.m.) and the potential's own ANC (fm^-1/2) there
  double width = 0;         ///< depth and range: the Gaussian r0 (fm)
  double rLow = 0, rHigh = 0;  ///< depth and range: the r found with the asked a over the widths scanned
};
/// The depth of the wave for a scattering length `a` (fm^(2l+1)) in the branch of the record's depth
/// (physics/ DepthForScatteringLength).
bool PotentialDepthForScatteringLength(const Project &project, const Record &record, int pairKey, int l, int twoJ,
                                       double a, PotentialDepth &out, std::string &why);
/// The depth and Gaussian width of the wave for a scattering length `a` and an effective range `r` in the branch of
/// the record's depth (physics/ DepthWidthForEffectiveRange; potential_background.md section 13.4).  `found` false with
/// the reason (the range of r reached) when the branch has no such pair.
bool PotentialDepthForEffectiveRange(const Project &project, const Record &record, int pairKey, int l, int twoJ, double a,
                                     double r, PotentialDepth &out, std::string &why);
/// The depth of the wave nearest the record's within +-30 % that binds it at `energy` (MeV, c.m., < 0: the level's
/// E_x minus the threshold; physics/ DepthBindingAt, the depth the set-up's refusal of `capture interior` prints), the
/// eigenvalue and the potential's own ANC there (PotentialBoundStateANC).
bool PotentialDepthForBinding(const Project &project, const Record &record, int pairKey, int l, int twoJ, double energy,
                              PotentialDepth &out, std::string &why);

/// Is the partial wave (pair, l, J^pi) in the calculation?  The engine makes a J^pi group only for the J^pi of an active
/// level (CompoundSystem::Fill), and with <config> includeAllJPi=1 also for every J^pi a particle channel with l up to
/// maxOrbitalMomentum (absent: the largest l of a particle channel of <levels>) makes (AddHardSphereGroups); a J^pi
/// without a group is left out, its potential, background and hard-sphere phase with it.
struct PartialWaveState {
  int twoJ = 0, parity = 1;
  bool level = false;   ///< an active level of this J^pi has the channel (pair, l)
  bool key = false;     ///< includeAllJPi would add the J^pi (the key on or not)
  bool included = false;  ///< level, or key with includeAllJPi on
};
/// The J^pi of the wave (pair `pairKey`, l, 2J; -1: every J of l) as the rows `rows` (every row of <levels>, as the
/// Levels tab writes them) and the two keys make them.
std::vector<PartialWaveState> PartialWaveStates(const std::vector<LevelRow> &rows, bool includeAllJPi,
                                                int maxOrbitalMomentum, int pairKey, int l, int twoJ);
/// A J^pi that the particle pairs of `rows` form with l <= lMax and that no active level carries.
struct MissingPartialWave {
  int twoJ = 0, parity = 1;
  int nucleusZ = 0, nucleusA = 0;  ///< the compound nucleus (Z1 + Z2, round(M1 + M2) of a pair that forms it)
  int pairKey = 0;                 ///< a particle pair that forms it, with its lowest l
  int l = 0;
  bool key = false;                ///< includeAllJPi adds it (its lowest l <= the key's maximum l)
};
/// Every J^pi without an active level, nucleus by nucleus, in the order of J then parity; `keyMaxL` the largest l the
/// key's groups take (maxOrbitalMomentum, or else the largest l of a particle channel of the rows; -1 for none).
std::vector<MissingPartialWave> PartialWavesWithoutLevel(const std::vector<LevelRow> &rows, int lMax,
                                                         int maxOrbitalMomentum, int &keyMaxL);

}  // namespace front
}  // namespace azure2

#endif
