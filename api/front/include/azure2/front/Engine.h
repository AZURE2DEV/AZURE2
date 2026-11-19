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
  ///  per pair of open particle channels c, c' (c <= c'):  U_re, U_im, U_abs (|U_cc'|), T_re, T_im, and on the diagonal
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

}  // namespace front
}  // namespace azure2

#endif
