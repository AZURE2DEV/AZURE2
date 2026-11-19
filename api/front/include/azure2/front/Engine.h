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

  /// AZURE2 --validate of the text as it stands (data files are read relative to Directory()).
  Validation Validate() const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

/// The table of project format 3 as JSON (AZURE2 --schema): sections, record kinds and their fields with titles,
/// docs, types, closed choices and units, quantities, units, diagnostic codes.
std::string SchemaJson();
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
};
/// One record kind: its sections, title, doc and fields in the schema's order.
struct KindInfo {
  std::string kind;
  std::string title;
  std::string doc;
  std::vector<std::string> sections;
  std::vector<FieldInfo> fields;
};
std::vector<KindInfo> RecordKinds();
/// The sections of records and the kinds each holds, in the schema's order.
std::vector<std::pair<std::string, std::vector<std::string>>> RecordSections();
/// The quantities a record kind may use (names, schema order), and the units of one quantity.
std::vector<std::string> QuantityNames(const std::string &kind);
std::vector<std::string> QuantityUnits(const std::string &quantity);
/// The diagnostic codes and their meaning.
std::vector<std::pair<std::string, std::string>> DiagnosticCodes();

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

/// The engine's version string ("2.1.0").
std::string Version();

}  // namespace front
}  // namespace azure2

#endif
