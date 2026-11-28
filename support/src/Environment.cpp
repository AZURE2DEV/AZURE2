#include "azure2/support/Environment.h"

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <stdexcept>

namespace azure2 {
namespace support {

namespace {

using K = SwitchKind;
using T = SwitchType;
constexpr const char *kManual = "reference/command_line.rst";

// One row per switch, in the order of enum Switch (checked when it compiles).  Columns: id, name, type, value
// when unset, kind, what it does, where it is documented, the program's option (kind numbers only), the range of a
// whole number its reader checks (low > high: none).
constexpr std::array<SwitchRow, kSwitchCount> kTable = {{
    {Switch::kCoulombTable, "AZURE2_COULOMB_TABLE", T::kText, "~/.azure2/coulomb_table.bin", K::kNumbers,
     "the Coulomb table file; the option also gives use_coulomb_table=1", kManual, "--coulomb-table=FILE", 1, 0},
    {Switch::kSrimDataDir, "AZURE2_SRIM_DATA_DIR", T::kText, "the compiled copy", K::kNumbers,
     "the directory of SRIM-2013's SCOEF03.DAT and SNUC03.dat for the stopping powers",
     kManual, "--srim-data-dir=DIR", 1, 0},
    {Switch::kOptionsFile, "AZURE_OPTIONS_FILE", T::kText, "none", K::kNumbers,
     "a file of command-line options, one per line, applied after those on the command line",
     kManual, "--options-file=FILE", 1, 0},
    {Switch::kRefineQuadrature, "AZURE_REFINE_QUADRATURE", T::kInteger, "1", K::kNumbers,
     "every sub-point grid gets K times its points (1..64; a bad value exits 2)",
     kManual, "--refine-quadrature=K", 1, 64},
    {Switch::kFoldDepthSteps, "AZURE_FOLD_DEPTH_STEPS", T::kInteger, "8", K::kNumbers,
     "panels per Gaussian width of a target's depth integral with a spread (0..1000; 0: the node intervals)",
     kManual, "--fold-depth-steps=K", 0, 1000},
    {Switch::kEcTableNodes, "AZURE2_EC_TABLE_NODES", T::kInteger, "ecTableNodesPerDecade (50)", K::kNumbers,
     "nodes per decade of the external-capture tables, over the <config> key", kManual, "--ec-table-nodes=N", 1, 0},
    {Switch::kSwapCouplingOrder, "AZURE2_SWAP_COUPLING_ORDER", T::kFlag, "off", K::kNumbers,
     "couples the two particle spins of a channel spin the other way round (a convention test)",
     kManual, "--swap-coupling-order", 1, 0},
    {Switch::kBoundShiftNumerical, "AZURE_BOUND_SHIFT_NUMERICAL", T::kFlag, "off", K::kNumbers,
     "the bound-state shift function from the numerical Whittaker derivative, not the U ratio",
     kManual, "--bound-shift-numerical", 1, 0},
    {Switch::kNoNeutralClosedForm, "AZURE_NO_NEUTRAL_CLOSED_FORM", T::kFlag, "off", K::kNumbers,
     "the neutral bound-state shift without its closed form", kManual, "--no-neutral-closed-form", 1, 0},
    {Switch::kShiftDerivative1eV, "AZURE_SHIFT_DE_1EV", T::kFlag, "off", K::kNumbers,
     "dS/dE by a 1 eV central difference of the exact S (the route before 2026-10-19)",
     kManual, "--shift-derivative-1ev", 1, 0},
    {Switch::kLevelSpace, "AZURE_LEVEL_SPACE", T::kFlag, "off", K::kNumbers,
     "every J group solved in level space (no channel-space solves); equal to rounding",
     kManual, "--level-space", 1, 0},
    {Switch::kThmByDifferences, "AZURE2_THM_FD", T::kFlag, "off", K::kNumbers,
     "the THM rows of the gradient and Jacobian by finite differences, not the adjoint",
     kManual, "--thm-by-differences", 1, 0},
    {Switch::kGradAdjoint, "AZURE_GRAD_ADJOINT", T::kFlag, "off", K::kNumbers,
     "the adjoint gradient instead of the Jacobian path; equal to rounding", kManual, "--gradient-adjoint", 1, 0},
    {Switch::kGradDense, "AZURE_GRAD_DENSE", T::kFlag, "off", K::kNumbers,
     "the dense level-matrix adjoint instead of the rank-one terms; equal to rounding",
     kManual, "--gradient-dense", 1, 0},
    {Switch::kFoldWeightsByDifferences, "AZURE_TE_WEIGHTS_FD", T::kFlag, "off", K::kNumbers,
     "the target-effect combiner weights by finite differences", kManual, "--fold-weights-by-differences", 1, 0},
    {Switch::kSerialChi2, "AZURE_SERIAL_CHI2", T::kFlag, "off", K::kCost,
     "the points of an evaluation one by one, not in the parallel pass (the same bits)", kManual, "", 1, 0},
    {Switch::kGradCheck, "AZURE_GRAD_CHECK", T::kText, "off", K::kPrints,
     "compares sampled analytic derivatives with finite differences; 'exit' stops after the check", kManual, "", 1, 0},
    {Switch::kGradCheckStep, "AZURE_GRAD_CHECK_H", T::kNumber, "1e-6", K::kPrints,
     "the relative step of the gradient check", kManual, "", 1, 0},
    {Switch::kGradCheckIndices, "AZURE_GRAD_CHECK_IDX", T::kText, "sampled", K::kPrints,
     "the parameters the gradient check differentiates (a list of indices)", kManual, "", 1, 0},
    {Switch::kGradCheckCount, "AZURE_GRAD_CHECK_N", T::kInteger, "sampled", K::kPrints,
     "how many parameters the gradient check samples", kManual, "", 1, 0},
    {Switch::kGradCheckRows, "AZURE_GRAD_CHECK_ROWS", T::kInteger, "60", K::kPrints,
     "how many rows of the gradient check's table are printed", kManual, "", 1, 0},
    {Switch::kGradCheckSegments, "AZURE_GRAD_CHECK_SEGS", T::kText, "all", K::kPrints,
     "the segments the gradient check differentiates", kManual, "", 1, 0},
    {Switch::kShiftGradCheck, "AZURE_SHIFT_GRAD_CHECK", T::kFlag, "off", K::kPrints,
     "checks the energy-shift Jacobian columns against full differences", kManual, "", 1, 0},
    {Switch::kFoldWeightsCheck, "AZURE_TE_WEIGHTS_CHECK", T::kText, "off", K::kPrints,
     "checks the target-effect combiner weights against finite differences", kManual, "", 1, 0},
    {Switch::kProfile, "AZURE_PROFILE", T::kFlag, "off", K::kPrints,
     "per-stage times and counters of the calculation on standard error", kManual, "", 1, 0},
    {Switch::kGradTiming, "AZURE_GRAD_TIMING", T::kFlag, "off", K::kPrints,
     "the time of every gradient call on standard error", kManual, "", 1, 0},
    {Switch::kGridTiming, "AZURE_GRID_TIMING", T::kFlag, "off", K::kPrints,
     "the time and the moves of every re-placement of the sub-point grids", kManual, "", 1, 0},
    {Switch::kPlanCheck, "AZURE_PLAN_CHECK", T::kFlag, "off", K::kPrints,
     "checks the calculation plan's energy tasks and node tables against the units (stops on a difference)",
     kManual, "", 1, 0},
    {Switch::kNodeCheck, "AZURE_NODE_CHECK", T::kFlag, "off", K::kPrints,
     "checks every unit's node and channel rows after each hand-over (stops on a difference)", kManual, "", 1, 0},
    {Switch::kPoolCheck, "AZURE_POOL_CHECK", T::kFlag, "off", K::kPrints,
     "compares the pooled, parallel chi-squared with a fresh serial one (first five calls)", kManual, "", 1, 0},
    {Switch::kChannelVerify, "AZURE_CHANNEL_VERIFY", T::kFlag, "off", K::kPrints,
     "compares every channel-space solve with the level-space solve (slow)", kManual, "", 1, 0},
    {Switch::kEcTableCheck, "AZURE2_EC_TABLE_CHECK", T::kFlag, "off", K::kPrints,
     "computes every external-capture integral also directly and reports the largest difference", kManual, "", 1, 0},
    {Switch::kTraceFcn, "AZURE_TRACE_FCN", T::kText, "off", K::kPrints,
     "prints every chi-squared evaluation of a fit; 'fresh' also on a fresh copy", kManual, "", 1, 0},
    {Switch::kDumpNonFinite, "AZURE_DUMP_NONFINITE", T::kNumber, "off", K::kPrints,
     "dumps an evaluation whose chi-squared is not finite or above the threshold (default 1e10)", kManual, "", 1, 0},
    {Switch::kDebugGrid, "AZR_DEBUG_GRID", T::kFlag, "off", K::kPrints,
     "prints the adaptive integration grids as they are built", kManual, "", 1, 0},
    {Switch::kPoisonResults, "AZURE_POISON_RESULTS", T::kFlag, "off", K::kPrints,
     "fills reused result storage with NaN at the start of every evaluation", kManual, "", 1, 0},
    {Switch::kTestFailEnergy, "AZURE_TEST_FAIL_E", T::kNumber, "off", K::kTest,
     "the calculation of a point at this lab energy (MeV) fails", kManual, "", 1, 0},
    {Switch::kTestFailExit, "AZURE_TEST_FAIL_EXIT", T::kInteger, "every pair", K::kTest,
     "... only for this exit pair key", kManual, "", 1, 0},
}};

constexpr bool InOrder() {
  for (int i = 0; i < kSwitchCount; i++)
    if (static_cast<int>(kTable[i].id) != i) return false;
  return true;
}
static_assert(InOrder(), "the rows of the switch table are in the order of enum Switch");

bool Filled(const char *value) { return value != nullptr && value[0] != '\0'; }

// "--refine-quadrature=" of "--refine-quadrature=K"; the whole option of a flag row.
std::string OptionStem(const SwitchRow &row) {
  const std::string option = row.option;
  const std::size_t equals = option.find('=');
  return equals == std::string::npos ? option : option.substr(0, equals + 1);
}

// The name a value is refused under: the option when an option gave it, else the variable.
std::string GivenAs(const Environment &environment, const SwitchRow &row) {
  const SwitchSource source = environment.Source(row.id);
  if (source != SwitchSource::kOption && source != SwitchSource::kOptionsFile) return row.name;
  const std::string stem = OptionStem(row);
  return stem.back() == '=' ? stem.substr(0, stem.size() - 1) : stem;
}

// The options file's lines that start with "--", spaces and tabs removed (as the program read them before).
std::vector<std::string> OptionsFileLinesOf(std::ifstream &in) {
  std::vector<std::string> lines;
  std::string line;
  while (std::getline(in, line)) {
    std::string trimmed;
    for (char c : line)
      if (c != ' ' && c != '\t' && c != '\n' && c != '\r') trimmed += c;
    if (trimmed.compare(0, 2, "--") == 0) lines.push_back(trimmed);
  }
  return lines;
}

// A JSON string (the values are paths and words).
std::string JsonString(const std::string &text) {
  std::string out = "\"";
  for (char c : text) {
    if (c == '"' || c == '\\')
      out += std::string("\\") + c;
    else if (static_cast<unsigned char>(c) < 0x20)
      out += ' ';
    else
      out += c;
  }
  return out + "\"";
}

}  // namespace

const char *SwitchSourceName(SwitchSource source) {
  switch (source) {
    case SwitchSource::kDefault: return "default";
    case SwitchSource::kEnvironment: return "environment";
    case SwitchSource::kOptionsFile: return "options file";
    case SwitchSource::kOption: return "option";
  }
  return "?";
}

const SwitchRow *ParseSwitchOption(const std::string &argument, std::string &value) {
  for (const SwitchRow &row : kTable) {
    if (!*row.option) continue;
    const std::string stem = OptionStem(row);
    if (stem.back() != '=') {  // a flag: the whole argument
      if (argument != stem) continue;
      value = "1";
      return &row;
    }
    // "--name=VALUE" only: the bare "--name" is another flag or an unknown option (--coulomb-table is a flag of
    // use_coulomb_table)
    if (argument.compare(0, stem.size(), stem) != 0) continue;
    if (argument.size() == stem.size())
      throw std::invalid_argument(stem + " needs a value: " + row.option + " (got '" + argument + "')");
    value = argument.substr(stem.size());
    return &row;
  }
  return nullptr;
}

long CheckedWholeNumber(const std::string &text, long low, long high, const std::string &what) {
  char *end = nullptr;
  errno = 0;
  const long n = text.empty() ? 0 : std::strtol(text.c_str(), &end, 10);
  if (text.empty() || *end != '\0' || errno == ERANGE || n < low || n > high)
    throw std::invalid_argument(what + " needs a whole number from " + std::to_string(low) + " to " +
                                std::to_string(high) + ", not '" + text + "'");
  return n;
}

const std::array<SwitchRow, kSwitchCount> &SwitchTable() { return kTable; }

const SwitchRow &RowOf(Switch id) { return SwitchTable()[static_cast<int>(id)]; }

const char *KindName(SwitchKind kind) {
  switch (kind) {
    case SwitchKind::kNumbers: return "numbers";
    case SwitchKind::kCost: return "cost";
    case SwitchKind::kPrints: return "prints";
    case SwitchKind::kTest: return "test";
  }
  return "?";
}

void Environment::Give(Switch id, const std::string &text, SwitchSource source) {
  set_[Index(id)] = true;
  text_[Index(id)] = text;
  source_[Index(id)] = source;
}

Environment Environment::Read(const Lookup &lookup, const std::vector<std::string> &arguments) {
  Environment e;
  // The variables, then the options file's switch options over them, then the command line's over both.
  for (const SwitchRow &row : SwitchTable())
    if (const char *value = lookup(row.name)) e.Give(row.id, value, SwitchSource::kEnvironment);
  std::vector<std::pair<Switch, std::string>> commandLine;
  for (const std::string &argument : arguments) {
    std::string value;
    if (const SwitchRow *row = ParseSwitchOption(argument, value)) commandLine.emplace_back(row->id, value);
  }
  std::string file = e.Text(Switch::kOptionsFile);
  for (const auto &given : commandLine)
    if (given.first == Switch::kOptionsFile) file = given.second;
  if (e.Set(Switch::kOptionsFile) || !file.empty()) {
    std::ifstream in(file);
    if (in) {
      e.optionsFileLines_ = OptionsFileLinesOf(in);
      for (const std::string &line : e.optionsFileLines_) {
        std::string value;
        const SwitchRow *row = ParseSwitchOption(line, value);
        if (!row) continue;
        if (row->id == Switch::kOptionsFile)
          throw std::invalid_argument("the options file " + file + " names another options file (" + line + ")");
        e.Give(row->id, value, SwitchSource::kOptionsFile);
      }
    } else {
      e.optionsFileUnreadable_ = true;
    }
  }
  for (const auto &given : commandLine) e.Give(given.first, given.second, SwitchSource::kOption);
  const char *home = lookup("HOME");
  if (!Filled(home)) home = lookup("USERPROFILE");
  e.home_ = home ? home : ".";
  if (const char *threads = lookup("OMP_NUM_THREADS")) e.ompNumThreads_ = threads;
  return e;
}

bool Environment::On(Switch id) const {
  const std::string &t = Text(id);
  return Set(id) && !t.empty() && t != "0";
}

double Environment::Number(Switch id, double unset) const { return Set(id) ? std::atof(Text(id).c_str()) : unset; }

long Environment::Integer(Switch id, long unset) const { return Set(id) ? std::atoi(Text(id).c_str()) : unset; }

std::vector<Switch> Environment::SetSwitches() const {
  std::vector<Switch> out;
  for (const SwitchRow &row : SwitchTable())
    if (Set(row.id)) out.push_back(row.id);
  return out;
}

long SwitchWholeNumber(const Environment &environment, Switch id, long unset) {
  const SwitchRow &row = RowOf(id);
  if (!environment.Set(id)) return unset;
  return CheckedWholeNumber(environment.Text(id), row.low, row.high, GivenAs(environment, row));
}

void CheckSwitchValues(const Environment &environment) {
  for (const SwitchRow &row : SwitchTable())
    if (row.low <= row.high) SwitchWholeNumber(environment, row.id, 0);
}

std::string SetSwitchesText(const Environment &environment) {
  std::string text;
  for (Switch id : environment.SetSwitches())
    if (RowOf(id).kind != SwitchKind::kNumbers)
      text += std::string("  ") + RowOf(id).name + "=" + environment.Text(id) + " (" + KindName(RowOf(id).kind) + ")\n";
  return text.empty() ? text : "Environment switches set:\n" + text;
}

std::string NumberSwitchesText(const Environment &environment) {
  std::string text;
  for (const SwitchRow &row : SwitchTable())
    if (row.kind == SwitchKind::kNumbers)
      text += std::string("  ") + row.option + " (" + row.name + ") = " +
              (environment.Set(row.id) ? environment.Text(row.id) : std::string(row.unset)) + " (" +
              SwitchSourceName(environment.Source(row.id)) + ")\n";
  return text;
}

std::string NumberSwitchesJson(const Environment &environment) {
  std::string json = "{";
  for (const SwitchRow &row : SwitchTable()) {
    if (row.kind != SwitchKind::kNumbers) continue;
    json += std::string(json.size() > 1 ? ", " : "") + JsonString(row.name) +
            ": {\"option\": " + JsonString(row.option) +
            ", \"value\": " + (environment.Set(row.id) ? JsonString(environment.Text(row.id)) : std::string("null")) +
            ", \"source\": " + JsonString(SwitchSourceName(environment.Source(row.id))) + "}";
  }
  return json + "}";
}

std::string SwitchOptionsHelp() {
  std::string text;
  for (const SwitchRow &row : SwitchTable()) {
    if (!*row.option) continue;
    std::string flag = std::string("\t") + row.option + ":";
    flag.resize(std::max<std::size_t>(flag.size() + 1, 24), ' ');
    text += flag + row.what + " (the variable " + row.name + ").\n";
  }
  return text;
}

namespace {

// The record, made once, by the first of ProcessEnvironment and StartProcessEnvironment to be called.
const Environment &MakeRecord(const std::vector<std::string> &arguments, bool &madeHere) {
  // The one reading of the process's environment (repo_consistency check 39).
  static const Environment record = [&] {
    madeHere = true;
    return Environment::Read([](const char *name) { return std::getenv(name); }, arguments);
  }();
  return record;
}

}  // namespace

const Environment &ProcessEnvironment() {
  bool madeHere = false;
  return MakeRecord({}, madeHere);
}

const Environment &StartProcessEnvironment(const std::vector<std::string> &arguments) {
  bool madeHere = false;
  const Environment &record = MakeRecord(arguments, madeHere);
  if (!madeHere)
    throw std::logic_error("the environment switches were read before the command line gave its options");
  return record;
}

}  // namespace support
}  // namespace azure2
