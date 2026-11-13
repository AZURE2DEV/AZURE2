#include "azure2/support/Environment.h"

#include <cstdlib>
#include <cstring>

namespace azure2 {
namespace support {

namespace {

using K = SwitchKind;
using T = SwitchType;
constexpr const char *kManual = "reference/command_line.rst";

// One row per switch, in the order of enum Switch (checked when it compiles).  Columns: id, name, type, value
// when unset, kind, what it does, where it is documented.
constexpr std::array<SwitchRow, kSwitchCount> kTable = {{
    {Switch::kCoulombTable, "AZURE2_COULOMB_TABLE", T::kText, "~/.azure2/coulomb_table.bin", K::kNumbers,
     "the Coulomb table file when --coulomb-table names none", kManual},
    {Switch::kSrimDataDir, "AZURE2_SRIM_DATA_DIR", T::kText, "the compiled copy", K::kNumbers,
     "the directory of SRIM-2013's SCOEF03.DAT and SNUC03.dat for the stopping powers", kManual},
    {Switch::kOptionsFile, "AZURE_OPTIONS_FILE", T::kText, "none", K::kNumbers,
     "a file of command-line options, one per line, applied after those on the command line", kManual},
    {Switch::kRefineQuadrature, "AZURE_REFINE_QUADRATURE", T::kInteger, "1", K::kNumbers,
     "every sub-point grid gets K times its points (1..64; a bad value exits 2)", kManual},
    {Switch::kEcTableNodes, "AZURE2_EC_TABLE_NODES", T::kInteger, "ecTableNodesPerDecade (50)", K::kNumbers,
     "nodes per decade of the external-capture tables, over the <config> key", kManual},
    {Switch::kSwapCouplingOrder, "AZURE2_SWAP_COUPLING_ORDER", T::kFlag, "off", K::kNumbers,
     "couples the two particle spins of a channel spin the other way round (a convention test)", kManual},
    {Switch::kBoundShiftNumerical, "AZURE_BOUND_SHIFT_NUMERICAL", T::kFlag, "off", K::kNumbers,
     "the bound-state shift function from the numerical Whittaker derivative, not the U ratio", kManual},
    {Switch::kNoNeutralClosedForm, "AZURE_NO_NEUTRAL_CLOSED_FORM", T::kFlag, "off", K::kNumbers,
     "the neutral bound-state shift without its closed form", kManual},
    {Switch::kShiftDerivative1eV, "AZURE_SHIFT_DE_1EV", T::kFlag, "off", K::kNumbers,
     "dS/dE by a 1 eV central difference of the exact S (the route before 2026-10-19)", kManual},
    {Switch::kLevelSpace, "AZURE_LEVEL_SPACE", T::kFlag, "off", K::kNumbers,
     "every J group solved in level space (no channel-space solves); equal to rounding", kManual},
    {Switch::kThmByDifferences, "AZURE2_THM_FD", T::kFlag, "off", K::kNumbers,
     "the THM rows of the gradient and Jacobian by finite differences, not the adjoint", kManual},
    {Switch::kGradAdjoint, "AZURE_GRAD_ADJOINT", T::kFlag, "off", K::kNumbers,
     "the adjoint gradient instead of the Jacobian path; equal to rounding", kManual},
    {Switch::kGradDense, "AZURE_GRAD_DENSE", T::kFlag, "off", K::kNumbers,
     "the dense level-matrix adjoint instead of the rank-one terms; equal to rounding", kManual},
    {Switch::kFoldWeightsByDifferences, "AZURE_TE_WEIGHTS_FD", T::kFlag, "off", K::kNumbers,
     "the target-effect combiner weights by finite differences", kManual},
    {Switch::kSerialChi2, "AZURE_SERIAL_CHI2", T::kFlag, "off", K::kCost,
     "the points of an evaluation one by one, not in the parallel pass (the same bits)", kManual},
    {Switch::kGradCheck, "AZURE_GRAD_CHECK", T::kText, "off", K::kPrints,
     "compares sampled analytic derivatives with finite differences; 'exit' stops after the check", kManual},
    {Switch::kGradCheckStep, "AZURE_GRAD_CHECK_H", T::kNumber, "1e-6", K::kPrints,
     "the relative step of the gradient check", kManual},
    {Switch::kGradCheckIndices, "AZURE_GRAD_CHECK_IDX", T::kText, "sampled", K::kPrints,
     "the parameters the gradient check differentiates (a list of indices)", kManual},
    {Switch::kGradCheckCount, "AZURE_GRAD_CHECK_N", T::kInteger, "sampled", K::kPrints,
     "how many parameters the gradient check samples", kManual},
    {Switch::kGradCheckRows, "AZURE_GRAD_CHECK_ROWS", T::kInteger, "60", K::kPrints,
     "how many rows of the gradient check's table are printed", kManual},
    {Switch::kGradCheckSegments, "AZURE_GRAD_CHECK_SEGS", T::kText, "all", K::kPrints,
     "the segments the gradient check differentiates", kManual},
    {Switch::kShiftGradCheck, "AZURE_SHIFT_GRAD_CHECK", T::kFlag, "off", K::kPrints,
     "checks the energy-shift Jacobian columns against full differences", kManual},
    {Switch::kFoldWeightsCheck, "AZURE_TE_WEIGHTS_CHECK", T::kText, "off", K::kPrints,
     "checks the target-effect combiner weights against finite differences", kManual},
    {Switch::kProfile, "AZURE_PROFILE", T::kFlag, "off", K::kPrints,
     "per-stage times and counters of the calculation on standard error", kManual},
    {Switch::kGradTiming, "AZURE_GRAD_TIMING", T::kFlag, "off", K::kPrints,
     "the time of every gradient call on standard error", kManual},
    {Switch::kGridTiming, "AZURE_GRID_TIMING", T::kFlag, "off", K::kPrints,
     "the time and the moves of every re-placement of the sub-point grids", kManual},
    {Switch::kPlanCheck, "AZURE_PLAN_CHECK", T::kFlag, "off", K::kPrints,
     "checks the calculation plan's energy tasks and node tables against the units (stops on a difference)", kManual},
    {Switch::kNodeCheck, "AZURE_NODE_CHECK", T::kFlag, "off", K::kPrints,
     "checks every unit's node and channel rows after each hand-over (stops on a difference)", kManual},
    {Switch::kPoolCheck, "AZURE_POOL_CHECK", T::kFlag, "off", K::kPrints,
     "compares the pooled, parallel chi-squared with a fresh serial one (first five calls)", kManual},
    {Switch::kChannelVerify, "AZURE_CHANNEL_VERIFY", T::kFlag, "off", K::kPrints,
     "compares every channel-space solve with the level-space solve (slow)", kManual},
    {Switch::kEcTableCheck, "AZURE2_EC_TABLE_CHECK", T::kFlag, "off", K::kPrints,
     "computes every external-capture integral also directly and reports the largest difference", kManual},
    {Switch::kTraceFcn, "AZURE_TRACE_FCN", T::kText, "off", K::kPrints,
     "prints every chi-squared evaluation of a fit; 'fresh' also on a fresh copy", kManual},
    {Switch::kDumpNonFinite, "AZURE_DUMP_NONFINITE", T::kNumber, "off", K::kPrints,
     "dumps an evaluation whose chi-squared is not finite or above the threshold (default 1e10)", kManual},
    {Switch::kGradDebug, "AZURE_GRAD_DEBUG", T::kFlag, "off", K::kPrints,
     "names the segment and the reason when the gradient falls back to differences", kManual},
    {Switch::kChannelDebug, "AZURE_CHANNEL_DEBUG", T::kFlag, "off", K::kPrints,
     "reports the formal representation of each J group", kManual},
    {Switch::kDebugGrid, "AZR_DEBUG_GRID", T::kFlag, "off", K::kPrints,
     "prints the adaptive integration grids as they are built", kManual},
    {Switch::kPoisonResults, "AZURE_POISON_RESULTS", T::kFlag, "off", K::kPrints,
     "fills reused result storage with NaN at the start of every evaluation", kManual},
    {Switch::kTestFailEnergy, "AZURE_TEST_FAIL_E", T::kNumber, "off", K::kTest,
     "the calculation of a point at this lab energy (MeV) fails", kManual},
    {Switch::kTestFailExit, "AZURE_TEST_FAIL_EXIT", T::kInteger, "every pair", K::kTest,
     "... only for this exit pair key", kManual},
}};

constexpr bool InOrder() {
  for (int i = 0; i < kSwitchCount; i++)
    if (static_cast<int>(kTable[i].id) != i) return false;
  return true;
}
static_assert(InOrder(), "the rows of the switch table are in the order of enum Switch");

bool Filled(const char *value) { return value != nullptr && value[0] != '\0'; }

}  // namespace

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

Environment Environment::Read(const Lookup &lookup) {
  Environment e;
  for (const SwitchRow &row : SwitchTable()) {
    const char *value = lookup(row.name);
    e.set_[Index(row.id)] = value != nullptr;
    if (value) e.text_[Index(row.id)] = value;
  }
  const char *home = lookup("HOME");
  if (!Filled(home)) home = lookup("USERPROFILE");
  e.home_ = home ? home : ".";
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

std::string SetSwitchesText(const Environment &environment) {
  std::string text;
  for (Switch id : environment.SetSwitches())
    text += std::string("  ") + RowOf(id).name + "=" + environment.Text(id) + " (" + KindName(RowOf(id).kind) + ")\n";
  return text.empty() ? text : "Environment switches set:\n" + text;
}

const Environment &ProcessEnvironment() {
  // The one reading of the process's environment (repo_consistency check 39), on the first call.
  static const Environment record = Environment::Read([](const char *name) { return std::getenv(name); });
  return record;
}

}  // namespace support
}  // namespace azure2
