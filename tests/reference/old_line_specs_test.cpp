// The old-line adaptor (convert/'s SpecsFromOldLines.h) and the engine's constructors that take its specs.
//   old_line_specs_test
// 1. Every isDiff code of a <segmentsData> and a <segmentsTest> line gives the Segment flags of the old engine's
//    line constructors (the expected table below is written from those constructors, not from the adaptor).
// 2. Composite components keep entrance, exit, fixed angle and scaling bit for bit.
// 3. The refusals the old engine made while reading a line come with its words, in its order (warning first).
// 4. Prior-centre rows: the old names equal the engine's parameter names; a row is written into its data set
//    (inactive lines counted); a name that matches no line is refused.
// 5. <abundances> rows and <targetInt> data-set lists.
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

#include "WorkingCopy.h"
#include "Segment.h"
#include "azure2/convert/SpecsFromOldLines.h"

Config *g_config = nullptr;  // defined by the executable; the library reads it

namespace {

using namespace azure2::experiment;
int failures = 0;
void check(const std::string &what, bool ok) {
  if (!ok) {
    std::printf("  FAIL  %s\n", what.c_str());
    failures++;
  }
}

SpecList<DatasetSpec> Data(const std::string &body, bool relativistic = false) {
  std::istringstream in("<segmentsData>\n" + body + "</segmentsData>\n");
  return azure2::convert::DatasetSpecsFromOldLines(in, relativistic);
}

SpecList<ExtrapolationSpec> Tests(const std::string &body) {
  std::istringstream in("<segmentsTest>\n" + body + "</segmentsTest>\n");
  return azure2::convert::ExtrapolationSpecsFromOldLines(in, false);
}

std::string OldDataText(int isDiff, const std::string &tail = "") {
  std::string extra = (isDiff % 10 == 2) ? " 1.5 1" : "";
  return "1 1 2 0.1 1.0 10 170 " + std::to_string(isDiff) + extra + " 1.0 0 0.0 0.0 0.0 0 data.dat" + tail + "\n";
}

// The flags of the old Segment(SegLine) / Segment(ExtrapLine): differential, c.m. differential, phase, total
// capture, cross-section component, angular distribution, coefficients, analyzing power, THM.
struct Flags {
  bool diff, cm, phase;
  int totalCapture, component;
  bool angDist, coefficients, analyzing, thm;
};

Flags Of(const Segment &s) {
  return {s.IsDifferential(), s.IsCMDifferential(), s.IsPhase(), s.IsTotalCapture(), s.GetCrossSectionComponent(),
          s.IsAngularDist(), s.IsCoefficients(), s.IsAnalyzingPower(), s.IsTHM()};
}

bool Same(const Flags &a, const Flags &b) {
  return a.diff == b.diff && a.cm == b.cm && a.phase == b.phase && a.totalCapture == b.totalCapture &&
         a.component == b.component && a.angDist == b.angDist && a.coefficients == b.coefficients &&
         a.analyzing == b.analyzing && a.thm == b.thm;
}

}  // namespace

int main() {
  // 1. Data codes (Segment(SegLine)): 1, 4, 7 differential; 4, 7 c.m.; 2 phase; 3 total capture; 5 / 6 the E1 / E2
  // component; 9 coefficients; 7 analyzing power; + 10 THM; any other code angle-integrated.
  const struct {
    int code;
    Flags f;
  } dataCodes[] = {
      {0, {false, false, false, 0, 0, false, false, false, false}}, {1, {true, false, false, 0, 0, false, false, false, false}},
      {2, {false, false, true, 0, 0, false, false, false, false}},  {3, {false, false, false, 1, 0, false, false, false, false}},
      {4, {true, true, false, 0, 0, false, false, false, false}},   {5, {false, false, false, 0, 1, false, false, false, false}},
      {6, {false, false, false, 0, 2, false, false, false, false}}, {7, {true, true, false, 0, 0, false, false, true, false}},
      {9, {false, false, false, 0, 0, false, true, false, false}},  {11, {true, false, false, 0, 0, false, false, false, true}},
      {14, {true, true, false, 0, 0, false, false, false, true}},   {17, {true, true, false, 0, 0, false, false, true, true}},
      {12, {false, false, true, 0, 0, false, false, false, true}},  {42, {false, false, false, 0, 0, false, false, false, true}},
  };
  for (const auto &c : dataCodes) {
    const SpecList<DatasetSpec> d = Data(OldDataText(c.code));
    const bool one = !d.failed && d.items.size() == 1;
    check("data isDiff " + std::to_string(c.code) + " read", one);
    if (one) check("data isDiff " + std::to_string(c.code) + " flags", Same(Of(Segment(d.items[0].spec)), c.f));
  }
  // Extrapolation codes (Segment(ExtrapLine)): 1, 5, 7 differential; 5, 7 c.m.; 3 angular distribution; 4 total.
  const struct {
    int code;
    std::string line;
    Flags f;
  } testCodes[] = {
      {1, "1 1 2 0.1 1.0 0.1 10 170 10 1\n", {true, false, false, 0, 0, false, false, false, false}},
      {3, "1 1 2 0.1 1.0 0.1 10 170 10 3 4\n", {false, false, false, 0, 0, true, false, false, false}},
      {4, "1 1 -1 0.1 1.0 0.1 0 0 0 4\n", {false, false, false, 1, 0, false, false, false, false}},
      {5, "1 1 2 0.1 1.0 0.1 10 170 10 5\n", {true, true, false, 0, 0, false, false, false, false}},
      {6, "1 1 2 0.1 1.0 0.1 0 0 0 6\n", {false, false, false, 0, 0, false, false, false, false}},
      {15, "1 1 2 0.1 1.0 0.1 10 170 10 15\n", {true, true, false, 0, 0, false, false, false, true}},
  };
  for (const auto &c : testCodes) {
    const SpecList<ExtrapolationSpec> t = Tests(c.line);
    const bool one = !t.failed && t.items.size() == 1;
    check("test isDiff " + std::to_string(c.code) + " read", one);
    if (!one) continue;
    const Segment s(t.items[0].spec);
    check("test isDiff " + std::to_string(c.code) + " flags", Same(Of(s), c.f));
    if (c.code == 3) check("test isDiff 3 order", s.GetMaxAngDistOrder() == 4);
  }

  // 2. Components: the -1 (scaling) layout with a fixed angle and a scaling, and one without either.
  {
    const SpecList<DatasetSpec> d = Data(OldDataText(0, " 1 1 -1 2 1 2 30.25 1.4999999999999999e-07 1 3 -999 -999"));
    const bool one = !d.failed && d.items.size() == 1;
    check("components read", one);
    if (one) {
      const ReactionSpec &r = d.items[0].spec.reaction;
      check("ratio", r.combination == ReactionSpec::Combination::Ratio && r.components.size() == 2);
      if (r.components.size() == 2) {
        check("component 1", r.components[0].entranceKey == 1 && r.components[0].exitKey == 2 && r.components[0].hasAngle &&
                                 r.components[0].angle == 30.25 && r.components[0].scaling == 1.4999999999999999e-07);
        check("component 2", r.components[1].exitKey == 3 && !r.components[1].hasAngle && r.components[1].scaling == 1.0);
      }
      const Segment s(d.items[0].spec);
      check("segment ratio", s.IsAdvanced() && s.GetLegacyOperationType() == 1 && s.GetComponents().size() == 2);
    }
    const SpecList<DatasetSpec> bad = Data(OldDataText(0, " 1 2 1 1 2"));
    check("operation type 2 refused", bad.failed && bad.error.find("neither 0 (sum) nor 1 (ratio)") != std::string::npos);
  }

  // 3. Refusals while reading, after the lines before them, with the warning of the refused line first.
  {
    const SpecList<DatasetSpec> d = Data(OldDataText(1) + OldDataText(19));
    check("isDiff 19 after line 1", d.failed && d.items.size() == 1 &&
                                        d.error == "ERROR: Data segment #2 uses isDiff 19 (THM angular-distribution "
                                                   "coefficients), which is not supported.");
    const SpecList<DatasetSpec> r = Data(OldDataText(11), true);
    check("THM relativistic", r.failed && r.items.empty() &&
                                  r.error == "ERROR: Data segment #1 is a THM segment; THM segments are not supported "
                                             "with relativistic=1.");
    const SpecList<DatasetSpec> p = Data(OldDataText(8));
    check("isDiff 8", p.failed && p.error.find("ERROR: Data segment #1 uses isDiff 8 (Polarization x Cross Section)") == 0);
    const SpecList<DatasetSpec> w = Data(OldDataText(0, " 1 0 3 1 2 1 3") + OldDataText(19));
    check("warning kept", w.failed && w.items.size() == 1 && !w.items[0].warning.empty());
    const SpecList<DatasetSpec> open = [] {
      std::istringstream in("<segmentsData>\n" + OldDataText(0));
      return azure2::convert::DatasetSpecsFromOldLines(in, false);
    }();
    check("unclosed section", open.failed && open.error.empty());
    const SpecList<ExtrapolationSpec> cov = Tests("1 1 2 0.1 1.0 0.1 0 0 0 0 covariance c.dat\n");
    check("test covariance", cov.failed && cov.error.find("carries a covariance block") != std::string::npos);
  }

  // 4. Prior centres.
  // The names the adaptor gives a data set's parameters (the engine's names of an old project), and the names its
  // prior_centre rows are matched against.
  {
    const SpecList<DatasetSpec> d = Data(OldDataText(0) + OldDataText(0));
    check("parameter names of data set 2",
          d.items.size() == 2 && d.items[1].spec.names.norm == "segment_2_norm" &&
              d.items[1].spec.names.energyShift == "segment_2_energy_shift" &&
              d.items[1].spec.names.energyShiftSqrt == "segment_2_energy_shift_sqrt" &&
              d.items[1].spec.names.thickness == "segment_2_thickness");
    for (int k = 1; k <= 40; k++) {
      check("norm name " + std::to_string(k), azure2::convert::OldSegmentParameterName(k, true) == "segment_" + std::to_string(k) + "_norm");
      check("shift name " + std::to_string(k),
            azure2::convert::OldSegmentParameterName(k, false) == "segment_" + std::to_string(k) + "_energy_shift");
    }
  }
  {
    SpecList<DatasetSpec> d = Data(OldDataText(0) + "0" + OldDataText(0).substr(1));
    std::istringstream rows("<parameterSettings>\nsegment_2_norm prior_centre 1.25\nsegment_1_energy_shift prior_centre "
                            "0.003\nsegment_1_energy_shift prior_centre 0.004\nwidth_1 fixed\n</parameterSettings>\n");
    const std::string error = azure2::convert::PriorCentresFromOldRows(&rows, d);
    check("centres read", error.empty() && d.items.size() == 2);
    if (d.items.size() == 2) {
      check("norm centre on the inactive line 2", d.items[1].spec.hasNormCentre && d.items[1].spec.normCentre == 1.25 &&
                                                      d.items[1].spec.normCentreSource == "segment_2_norm prior_centre 1.25");
      check("last shift row wins", d.items[0].spec.hasShiftCentre && d.items[0].spec.shiftCentre == 0.004);
    }
    std::istringstream bad("<parameterSettings>\nsegment_3_norm prior_centre 1.0\n</parameterSettings>\n");
    const std::string refused = azure2::convert::PriorCentresFromOldRows(&bad, d);
    check("unknown line refused", refused.find("names no parameter a prior_centre row can move: it takes segment_1_norm "
                                               "... segment_2_norm") != std::string::npos);
  }

  // 5. Abundances and the data sets of a <targetInt> line.
  {
    std::istringstream in("<abundances>\nfoil: 3 6 0.075 1 0.01\n# comment\nfoil: sum 0.001\n</abundances>\n");
    const AbundanceSpecs a = azure2::convert::AbundanceSpecsFromOldLines(&in);
    check("abundances", !a.failed && a.abundances.present && a.abundances.rows.size() == 2 &&
                            a.abundances.rows[0].sample == "foil" && a.abundances.rows[0].z == 3 &&
                            a.abundances.rows[0].value == 0.075 && a.abundances.rows[0].vary && a.abundances.rows[1].isSum);
    std::istringstream badIn("<abundances>\n3 6 x\n</abundances>\n");
    const AbundanceSpecs b = azure2::convert::AbundanceSpecsFromOldLines(&badIn);
    check("bad abundance", b.failed && b.error.find("ERROR: <abundances> line \"3 6 x\" is not") == 0);
    check("no abundances", !azure2::convert::AbundanceSpecsFromOldLines(nullptr).abundances.present);
    std::istringstream te("<targetInt>\n1 \"1,3-5\" 20 1 0.001 0 0 \"\" 0 0 0 0 \"\" 0\n</targetInt>\n");
    const SpecList<TargetEffectSpec> t = azure2::convert::TargetEffectSpecsFromOldLines(te);
    check("target effect", !t.failed && t.items.size() == 1 && t.items[0].spec.active &&
                               t.items[0].spec.datasetKeys == std::vector<int>({1, 3, 4, 5}) && t.items[0].spec.convolution);
  }

  if (failures == 0) std::printf("old_line_specs: all checks passed\n");
  return failures == 0 ? 0 : 1;
}
