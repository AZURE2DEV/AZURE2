// Schema::Build (docs/upgrade/project_format3.md step 3b): one case per check it makes (each project text differs from
// a valid one in one place and gives exactly one diagnostic, with that code; check 24), and the specs of a valid
// project member by member (reaction resolution with states, the pair-key rule, observables, windows, norm and shift,
// a shared normalization, an extrapolation grid by decision F10, the effect of a target and a beam, THM membership).
//   experiment_build_test
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "azure2/experiment/ProjectReader.h"
#include "azure2/io/Document.h"
#include "azure2/io/RecordReader.h"

using namespace azure2::experiment;
using azure2::io::Document;

namespace {

int failures = 0;
int checks = 0;

void check(const std::string &what, bool ok) {
  checks++;
  if (!ok) {
    std::printf("  FAIL  %s\n", what.c_str());
    failures++;
  }
}

// p+15N (1), a+12C (2), a+12C[4.43891] (3), g+16O (4), d+12C twice (5, 6), p+12C (7)
LevelTables Pairs() {
  auto pair = [](int key, int type, int z1, int a1, int z2, int a2, int twoJ1, int twoJ2, double e2) {
    PairInfo p;
    p.key = key;
    p.kind = type;
    p.lightZ = z1;
    p.lightA = a1;
    p.heavyZ = z2;
    p.heavyA = a2;
    p.twoJLight = twoJ1;
    p.twoJHeavy = twoJ2;
    p.heavyExcitation = e2;
    return p;
  };
  LevelTables t;
  t.pairs = {pair(1, 0, 1, 1, 7, 15, 1, 1, 0), pair(2, 0, 2, 4, 6, 12, 0, 0, 0), pair(3, 0, 2, 4, 6, 12, 0, 4, 4.43891),
             pair(4, 10, 0, 0, 8, 16, 2, 0, 0), pair(5, 0, 1, 2, 6, 12, 2, 0, 0), pair(6, 0, 1, 2, 6, 12, 2, 0, 0),
             pair(7, 0, 1, 1, 6, 12, 1, 0, 0)};
  return t;
}

BuiltProject BuildText(const std::string &text, bool relativistic = false, bool readFiles = false) {
  const Document doc = Document::Parse(text);
  std::vector<SectionRecords> sections;
  for (const Document::Section &s : doc.Sections()) {
    if (s.depth != 0 || Schema::KindsOfSection(s.name).empty()) continue;
    sections.push_back({s.name, azure2::io::ReadRecords(doc, s, Schema::KindsOfSection(s.name)).set.records});
  }
  const ValidationResult v = Schema::Validate(sections);
  BuiltProject b;
  if (!v.Ok()) {
    b.diagnostics = v.diagnostics;
    for (Diagnostic &d : b.diagnostics) d.message = "(validation) " + d.message;
    return b;
  }
  BuildOptions o;
  o.relativistic = relativistic;
  o.readDataFiles = readFiles;
  return Schema::Build(v, Pairs(), o);
}

const std::string kDataset =
    "dataset d1\n"
    "    reaction       15N(p,a)12C\n"
    "    quantity       sigma\n"
    "    unit           b\n"
    "    file           build_test_data.dat\n"
    "    columns        E_lab[MeV] _ value error\n"
    "    energy_window  0.2 1.4 MeV\n"
    "    norm           1 fit prior 1 5%\n";

std::string Set(const std::string &base, const std::string &key, const std::string &line) {
  std::string out;
  std::size_t a = 0;
  bool done = false;
  while (a < base.size()) {
    const std::size_t b = base.find('\n', a);
    const std::string l = base.substr(a, b - a);
    if (!done && l.compare(0, 4 + key.size() + 1, "    " + key + " ") == 0) {
      if (!line.empty()) out += line + "\n";
      done = true;
    } else {
      out += l + "\n";
    }
    a = b + 1;
  }
  return out;
}
std::string Data(const std::string &records) { return "<data>\n" + records + "</data>\n"; }
std::string Rename(const std::string &record, const std::string &name) {
  return "dataset " + name + record.substr(record.find('\n'));
}

void ExpectOne(const std::string &code, const BuiltProject &b) {
  const bool ok = b.diagnostics.size() == 1 && b.diagnostics[0].code == code;
  check(code + " alone", ok);
  if (!ok)
    for (const Diagnostic &d : b.diagnostics) std::printf("        got %s: %s\n", d.code.c_str(), d.message.c_str());
}

}  // namespace

int main() {
  std::ofstream("build_test_data.dat") << "0.3 0 1.0 0.1\n";
  const std::string targets =
      "<targets>\n"
      "target t1\n    kind layer\n    active 15N\n    areal_density 1e18 atoms/cm2 fixed\n"
      "    stopping formula \"a0\" a0=1e-20 MeV*cm2/atom\n    straggling none\n    screening 1000 eV fit prior 1000 200\n"
      "beam b1\n    spread 1 keV lab\n"
      "beam b2\n    spread 0.002 MeV lab\n"
      "</targets>\n";

  // a valid project and its specs
  {
    const std::string text =
        Data("normalization D\n    value 1.1 fit prior 1 8%\n" + kDataset +
             Rename(Set(Set(kDataset, "reaction", "    reaction 15N(p,a)12C[4.439] + 2*16O(g,a)12C @ 90 deg lab"), "norm", "    norm D"), "d2") +
             Rename(Set(Set(Set(kDataset, "quantity", "    quantity dsigma_domega"), "columns",
                            "    columns E_lab[MeV] theta_cm[deg] value error"),
                        "unit", "    unit b/sr") +
                        "    shift offset 0.001 MeV fit prior 0.002 0.003\n",
                    "d3") +
             Rename(Set(Set(Set(kDataset, "norm", "    norm D"), "reaction", "    reaction 12C(d,d)12C"), "quantity", "    quantity sigma_nuclear") +
                        "    entrance_pair 6\n    exit_pair 6\n    target t1\n    beam b1\n    quadrature points 20\n",
                    "d4") +
             Rename(Set(Set(kDataset, "norm", "    norm profiled"), "unit", "    unit arbitrary") + "    thm LC\n", "d5")) +
        "<extrapolations>\nextrapolation x1\n    reaction 15N(p,g)16O\n    quantity sigma\n    unit b\n"
        "    energies range 0.1 0.4 MeV lab step 0.1\n</extrapolations>\n" + targets;
    const BuiltProject b = BuildText(text);
    for (const Diagnostic &d : b.diagnostics) std::printf("        %s: %s\n", d.code.c_str(), d.message.c_str());
    check("the valid project builds without a diagnostic", b.Ok() && b.diagnostics.empty());
    if (b.data.datasets.items.size() == 5 && b.extrapolations.extrapolations.items.size() == 1) {
      const DatasetSpec &d1 = b.data.datasets.items[0].spec, &d2 = b.data.datasets.items[1].spec;
      const DatasetSpec &d3 = b.data.datasets.items[2].spec, &d4 = b.data.datasets.items[3].spec;
      check("d1: key 1, pairs 1 -> 2, cross section", d1.key == 1 && d1.reaction.entranceKey == 1 && d1.reaction.exitKey == 2 &&
                                                         d1.reaction.observable == Observable::CrossSection);
      check("d1: energy window 0.2 .. 1.4 MeV", d1.minE == 0.2 && d1.maxE == 1.4);
      check("d1: norm 1 fit, 5 %, centred on the value", d1.norm == 1.0 && d1.varyNorm && d1.normError == 5.0 && !d1.hasNormCentre);
      check("d2: [4.439] matches 4.43891 (pair 3); a sum with a scaled term at 90 deg", d2.reaction.exitKey == 3 &&
                d2.reaction.combination == ReactionSpec::Combination::Sum && d2.reaction.components.size() == 1 &&
                d2.reaction.components[0].entranceKey == 4 && d2.reaction.components[0].exitKey == 2 &&
                d2.reaction.components[0].scaling == 2.0 && d2.reaction.components[0].hasAngle &&
                d2.reaction.components[0].angle == 90.0);
      check("d2, d4: the shared normalization D's values and label", d2.normGroup == "D" && d4.normGroup == "D" && d2.norm == 1.1 &&
                                                                         d2.normError == 8.0 && d2.hasNormCentre && d2.normCentre == 1.0);
      check("d3: theta_cm -> differential c.m.; the shift offset with its prior centre",
            d3.reaction.observable == Observable::DifferentialCM && d3.energyShift == 0.001 && d3.varyEnergyShift &&
                d3.energyShiftError == 0.003 && d3.hasShiftCentre && d3.shiftCentre == 0.002);
      check("d4: entrance_pair 6 picks the second d+12C pair", d4.reaction.entranceKey == 6 && d4.reaction.exitKey == 6);
      const DatasetSpec &d5 = b.data.datasets.items[4].spec;
      check("d5: THM experiment LC, profiled norm", d5.reaction.thm && d5.varyNorm && b.thm.experiments.size() == 1 &&
                                                       b.thm.experiments[0].name == "LC" &&
                                                       b.thm.experiments[0].datasetKeys == std::vector<int>{5});
      const ExtrapolationSpec &x1 = b.extrapolations.extrapolations.items[0].spec;
      check("x1: E_k = lo + k step (F10), 4 points", x1.energies.size() == 4 && x1.energies[3] == 0.1 + 3 * 0.1 && x1.minE == 0.1 &&
                                                         x1.maxE == 0.4 && x1.reaction.exitKey == 4);
      check("one effect: target t1 and beam b1 on data set 4", b.data.targetEffects.items.size() == 1 &&
                                                                  b.data.targetEffects.items[0].spec.datasetKeys == std::vector<int>{4});
      if (b.data.targetEffects.items.size() == 1) {
        const TargetEffectSpec &e = b.data.targetEffects.items[0].spec;
        check("the effect: a layer of 1e18 atoms/cm2, screening 1000 eV fitted with 200 eV, sigma 1 keV = 0.001 MeV, 20 points",
              e.active && e.targetIntegration && e.density == 1e18 && e.stoppingPowerEquation == "a0" && e.screening &&
                  e.screeningUe == 1000 && e.screeningVary && e.screeningSigma == 200 && e.convolution && e.sigma == 1 * 1e-3 &&
                  e.numIntegrationPoints == 20 && e.resonanceWidthMultiplier == 20 && e.pointsPerWidth == 50);
      }
    } else {
      check("5 data sets and 1 extrapolation", false);
    }
  }

  // one case per build check: each gives exactly one diagnostic, with its code
  ExpectOne("E_UNKNOWN_PAIR", BuildText(Data(Set(kDataset, "reaction", "    reaction 14N(p,a)11C"))));
  const std::string dd = Set(Set(kDataset, "reaction", "    reaction 12C(d,d)12C"), "quantity", "    quantity sigma_nuclear");
  ExpectOne("E_AMBIGUOUS_PAIR", BuildText(Data(dd)));
  ExpectOne("E_PAIR_KEY", BuildText(Data(kDataset + "    entrance_pair 1\n")));
  ExpectOne("E_PAIR_KEY", BuildText(Data(dd + "    entrance_pair 1\n    exit_pair 5\n")));
  ExpectOne("E_DATA_FILE", BuildText(Data(Set(kDataset, "file", "    file no_such_file.dat")), false, true));
  check("the data file read in the valid case", BuildText(Data(kDataset), false, true).diagnostics.empty());
  ExpectOne("E_DATA_FILE", BuildText(Data(Set(kDataset, "columns", "    columns E_lab[MeV] _ value error _ _")), false, true));
  ExpectOne("E_NOT_BUILT", BuildText(Data(Set(kDataset, "unit", "    unit mb"))));
  const std::string thm = Set(Set(kDataset, "norm", "    norm profiled"), "unit", "    unit arbitrary") + "    thm global\n";
  check("a THM data set builds", BuildText(Data(thm)).diagnostics.empty());
  ExpectOne("E_THM_RELATIVISTIC", BuildText(Data(thm), true));
  ExpectOne("E_THM_COVARIANCE", BuildText(Data(thm + "    covariance cov.txt\n")));
  ExpectOne("E_FOLD_CONFLICT",
            BuildText(Data(kDataset + "    target t1\n    beam b1\n    quadrature points 20\n" +
                           Rename(kDataset, "d2") + "    target t1\n    beam b2\n    quadrature points 20\n") +
                      targets));
  const std::string phase = Set(Set(Set(Set(kDataset, "quantity", "    quantity phase_shift"), "unit", "    unit deg"), "reaction",
                                    "    reaction 12C(p,p)12C"),
                                "norm", "    norm 1 fixed");
  check("a phase shift (J, l, s) = (1/2, 0, 1/2) of p+12C builds", BuildText(Data(phase + "    J 1/2\n    l 0\n    s 1/2\n")).diagnostics.empty());
  ExpectOne("E_FIELD_NOT_APPLICABLE", BuildText(Data(phase + "    J 1/2\n    l 0\n    s 3/2\n")));
  // the format-3 refusals of the old engine's reading checks that the table makes already (no Build check needed)
  check("a covariance on an extrapolation is refused by the validator (E_UNKNOWN_FIELD)",
        BuildText("<extrapolations>\nextrapolation x1\n    reaction 15N(p,a)12C\n    quantity sigma\n    unit b\n"
                  "    energies list 0.3 MeV lab\n    covariance c.txt\n</extrapolations>\n")
                .diagnostics.size() == 1);
  check("the removed polarization quantity is not in the table (E_NOT_IN_LIST)",
        BuildText(Data(Set(kDataset, "quantity", "    quantity polarization"))).diagnostics.size() >= 1);

  std::printf("experiment_build: %d checks, %d failed\n", checks, failures);
  return failures == 0 ? 0 : 1;
}
