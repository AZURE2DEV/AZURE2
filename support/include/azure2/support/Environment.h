#ifndef AZURE2_SUPPORT_ENVIRONMENT_H
#define AZURE2_SUPPORT_ENVIRONMENT_H

#include <array>
#include <functional>
#include <string>
#include <vector>

/*!
 * The environment switches: one table of every variable the engine, its libraries and the binding read, and one
 * reader of the process's environment (docs/upgrade/evaluation.md section 8; the manual's list is
 * docs/source/reference/command_line.rst, "Environment switches", checked against this table by repo_consistency).
 *
 * ProcessEnvironment() is the record of the process's environment, read once on its first call and never written
 * after: an immutable table, like the constant tables of docs/upgrade/ARCHITECTURE.md section 9.  The switches that
 * can move a number are copied from it into the run's Config when the Config is made (Config::Reset), and read there;
 * the printing and checking switches and the test hooks are read from the record where they are used.  Nothing else
 * calls getenv (repo_consistency check 39).
 *
 * Add a switch: a value of Switch, its row in the table of src/Environment.cpp (same order; a switch of kind numbers
 * with its option), its line in the manual's list, and its reader (CLAUDE.md, recipe "Add a diagnostic switch").  The
 * program's parser, --help and --describe take the option from the row.
 */
namespace azure2 {
namespace support {

/// What setting a switch changes.
enum class SwitchKind {
  kNumbers,  ///< an output number can move: other physics, another input, or another route equal to rounding
  kCost,     ///< the same bits in another time
  kPrints,   ///< prints, dumps or checks (and may exit)
  kTest      ///< injects a failure, for a test
};

/// How a switch's value is read.
enum class SwitchType {
  kFlag,     ///< set or not (each reader states whether "" and "0" count as set: Environment::Set or On)
  kText,     ///< a path or a word
  kNumber,   ///< a real number (atof)
  kInteger   ///< a whole number (atoi, or a checked parse where the reader refuses a bad value)
};

/// The switches, in the order of the table's rows.
enum class Switch : int {
  // numbers
  kCoulombTable,
  kSrimDataDir,
  kOptionsFile,
  kRefineQuadrature,
  kEcTableNodes,
  kSwapCouplingOrder,
  kBoundShiftNumerical,
  kNoNeutralClosedForm,
  kShiftDerivative1eV,
  kLevelSpace,
  kThmByDifferences,
  kGradAdjoint,
  kGradDense,
  kFoldWeightsByDifferences,
  // cost
  kSerialChi2,
  // prints
  kGradCheck,
  kGradCheckStep,
  kGradCheckIndices,
  kGradCheckCount,
  kGradCheckRows,
  kGradCheckSegments,
  kShiftGradCheck,
  kFoldWeightsCheck,
  kProfile,
  kGradTiming,
  kGridTiming,
  kPlanCheck,
  kNodeCheck,
  kPoolCheck,
  kChannelVerify,
  kEcTableCheck,
  kTraceFcn,
  kDumpNonFinite,
  kDebugGrid,
  kPoisonResults,
  // test
  kTestFailEnergy,
  kTestFailExit,
  kCount
};

constexpr int kSwitchCount = static_cast<int>(Switch::kCount);

/// One row of the table.
struct SwitchRow {
  Switch id;
  const char *name;          ///< the environment variable
  SwitchType type;
  const char *unset;         ///< what holds when it is not set
  SwitchKind kind;
  const char *what;          ///< one line: what setting it does
  const char *documented;    ///< the page of docs/source that documents it
  /// The program's option of a switch of kind numbers ("" for the other kinds): "--level-space" for a flag (its value
  /// is "1", as the variable set to 1), "--refine-quadrature=K" for one that takes a value (docs/upgrade/evaluation.md
  /// section 9).
  const char *option;
  /// A whole number its reader checks: low <= value <= high, else the run exits 2 (low > high: not checked; the reader
  /// reads it as it always did).
  long low;
  long high;
};

/// Where the value in force of a switch came from, the first that gives one: the option on the command line, the
/// options file (--options-file / AZURE_OPTIONS_FILE), the environment variable; else the default.
enum class SwitchSource { kDefault, kEnvironment, kOptionsFile, kOption };
/// "default", "environment", "options file", "option".
const char *SwitchSourceName(SwitchSource source);

/// The table, one row per Switch in its order.
const std::array<SwitchRow, kSwitchCount> &SwitchTable();
/// The row of `id`.
const SwitchRow &RowOf(Switch id);
/// "numbers", "cost", "prints", "test".
const char *KindName(SwitchKind kind);

/// The row whose option `argument` is ("--level-space", or "--refine-quadrature=4" for "--refine-quadrature=K"),
/// with its value in `value` ("1" for a flag); null when it is no switch option (also the bare "--srim-data-dir").
/// std::invalid_argument for a value row given with an empty value ("--srim-data-dir=").
const SwitchRow *ParseSwitchOption(const std::string &argument, std::string &value);

/// The whole number of `text` with low <= N <= high, or std::invalid_argument "<what> needs a whole number from <low>
/// to <high>, not '<text>'" (the one checked parse of the command line's values and of the checked switches).
long CheckedWholeNumber(const std::string &text, long low, long high, const std::string &what);

/// The values of the switches in an environment: for each, whether it is set, its text and where it came from.
class Environment {
 public:
  /// The value of an environment variable by name, or null when it is not set.
  using Lookup = std::function<const char *(const char *name)>;

  /// The record of the environment `lookup` reads: every switch of the table, the home directory (HOME, else
  /// USERPROFILE; "." without either) for the default Coulomb-table path, and OMP_NUM_THREADS (for what --threads
  /// says it overrides).  `arguments`: the program's command line (argv[1..]); its switch options, then those of the
  /// options file (named by --options-file=FILE, else AZURE_OPTIONS_FILE; its lines starting with "--"), win over the
  /// variables.  std::invalid_argument for a switch option without its value, or an options file naming another.
  static Environment Read(const Lookup &lookup, const std::vector<std::string> &arguments = {});

  /// The variable is set, to anything (also "" or "0").
  bool Set(Switch id) const { return set_[Index(id)]; }
  /// The variable is set, not empty and not "0".
  bool On(Switch id) const;
  /// Its text; "" when it is not set.
  const std::string &Text(Switch id) const { return text_[Index(id)]; }
  /// atof of its text when it is set (so "" is 0), else `unset`.
  double Number(Switch id, double unset) const;
  /// atoi of its text when it is set (so "" is 0), else `unset`.
  long Integer(Switch id, long unset) const;
  /// HOME, else USERPROFILE, else "." (each only when set and not empty).
  const std::string &Home() const { return home_; }
  /// The switches that are set, in the table's order.
  std::vector<Switch> SetSwitches() const;
  /// Where the value of `id` came from.
  SwitchSource Source(Switch id) const { return source_[Index(id)]; }
  /// The options file's lines that start with "--" (spaces and tabs removed), in order; empty without a file.
  const std::vector<std::string> &OptionsFileLines() const { return optionsFileLines_; }
  /// An options file is named but cannot be read.
  bool OptionsFileUnreadable() const { return optionsFileUnreadable_; }
  /// OMP_NUM_THREADS as set ("" when it is not).
  const std::string &OmpNumThreads() const { return ompNumThreads_; }

 private:
  static int Index(Switch id) { return static_cast<int>(id); }
  void Give(Switch id, const std::string &text, SwitchSource source);
  std::array<bool, kSwitchCount> set_{};
  std::array<std::string, kSwitchCount> text_;
  std::array<SwitchSource, kSwitchCount> source_{};
  std::string home_;
  std::vector<std::string> optionsFileLines_;
  bool optionsFileUnreadable_ = false;
  std::string ompNumThreads_;
};

/// The whole number in force of a checked switch (a row with low <= high), else `unset` when it is not set;
/// std::invalid_argument naming the option or the variable its value came from.
long SwitchWholeNumber(const Environment &environment, Switch id, long unset);
/// Every checked switch of `environment` (SwitchWholeNumber): the program calls it at start, so a bad value exits 2
/// whatever gave it.
void CheckSwitchValues(const Environment &environment);

/// The switches set in `environment` that are not of kind numbers, for AZURE2 --describe: "Environment switches set:"
/// and one line "  NAME=value (kind)" per switch in the table's order; "" when none is set.
std::string SetSwitchesText(const Environment &environment);
/// Every switch of kind numbers with its value in force and its source, for AZURE2 --describe: one line
/// "  --option (VARIABLE) = value (source)" each, in the table's order (the row's `unset` text when not set).
std::string NumberSwitchesText(const Environment &environment);
/// The same as a JSON object {"VARIABLE": {"option": ..., "value": ..., "source": ...}, ...} (value null when not set).
std::string NumberSwitchesJson(const Environment &environment);
/// The switch options for AZURE2 --help, from the table: one line "\t--option:" + what it does + the variable each.
std::string SwitchOptionsHelp();

/// The process's environment, made on the first call and never again: a variable set later in the same process changes
/// nothing.  The program's main makes it first, with its command line (StartProcessEnvironment); a pyazr session's
/// first Config reads it from the environment alone.
const Environment &ProcessEnvironment();
/// The program's main, before anything reads the record: the record made with the command line's switch options
/// (Environment::Read).  std::logic_error when the record was made before (its options would be silently lost);
/// std::invalid_argument as Environment::Read.
const Environment &StartProcessEnvironment(const std::vector<std::string> &arguments);

}  // namespace support
}  // namespace azure2

#endif
