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
 * Add a switch: a value of Switch, its row in the table of src/Environment.cpp (same order), its line in the manual's
 * list, and its reader (CLAUDE.md, recipe "Add a diagnostic switch").
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
  kGradDebug,
  kChannelDebug,
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
};

/// The table, one row per Switch in its order.
const std::array<SwitchRow, kSwitchCount> &SwitchTable();
/// The row of `id`.
const SwitchRow &RowOf(Switch id);
/// "numbers", "cost", "prints", "test".
const char *KindName(SwitchKind kind);

/// The values of the switches in an environment: for each, whether it is set and its text.
class Environment {
 public:
  /// The value of an environment variable by name, or null when it is not set.
  using Lookup = std::function<const char *(const char *name)>;

  /// The record of the environment `lookup` reads: every switch of the table, and the home directory (HOME, else
  /// USERPROFILE; "." without either) for the default Coulomb-table path.
  static Environment Read(const Lookup &lookup);

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

 private:
  static int Index(Switch id) { return static_cast<int>(id); }
  std::array<bool, kSwitchCount> set_{};
  std::array<std::string, kSwitchCount> text_;
  std::string home_;
};

/// The process's environment, read on the first call (the CLI's main calls it first; a pyazr session's first Config
/// reads it) and never again: a variable set later in the same process changes nothing.
const Environment &ProcessEnvironment();

}  // namespace support
}  // namespace azure2

#endif
