// Schema::Build (docs/upgrade/project_format3.md step 3b): one case per check it makes (each project text differs from
// a valid one in one place and gives exactly one diagnostic, with that code; check 24), and the specs of a valid
// project member by member (reaction resolution with states, the pair-key rule, observables, windows, norm and shift,
// a shared normalization, an extrapolation grid by decision F10, the effect of a target and a beam, THM membership).
//   experiment_build_test
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "azure2/experiment/ProjectReader.h"
#include "azure2/experiment/SpinOrbitChecks.h"
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

// p+15N (1), a+12C (2), a+12C[4.43891] (3), g+16O (4), d+12C twice (5, 6), p+12C (7), a+a (8, identical),
// a+a (9, not identical: a pair whose particles differ in a column the label does not show)
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
             pair(7, 0, 1, 1, 6, 12, 1, 0, 0), pair(8, 0, 2, 4, 2, 4, 0, 0, 0), pair(9, 0, 2, 4, 2, 4, 0, 0, 0)};
  t.pairs[7].identical = true;
  t.pairs[0].lightMass = 1.00783;  // p + 15N: the masses a c.m. grid is converted with
  t.pairs[0].heavyMass = 15.0001;
  return t;
}

BuiltProject BuildText(const std::string &text, bool relativistic = false, bool readFiles = false,
                       DataFiles dataFiles = DataFiles::Required, const LevelTables &pairs = Pairs()) {
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
  o.dataFiles = dataFiles;
  return Schema::Build(v, pairs, o);
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

// Units (project_format3.md sections 3.3, 7.2): the factor of every number Build holds and of every column of a data
// file comes from the unit table; a per-charge yield and an S factor are factors on the model.
void UnitCases(const std::string &targets) {
  auto data = [](const BuiltProject &b) -> const DatasetSpec & { return b.data.datasets.items.at(0).spec; };
  const BuiltProject mb = BuildText(Data(Set(kDataset, "unit", "    unit mb")));
  check("unit mb: the value column times 1e-3 at the points, the model unscaled",
        mb.Ok() && data(mb).columns.value == 1e-3 && data(mb).reaction.scale.Identity() && data(mb).columns.energy == 1.0);
  const BuiltProject keV = BuildText(Data(Set(Set(kDataset, "columns", "    columns E_cm[keV] _ value error"), "energy_window",
                                              "    energy_window 200 1400 keV")));
  {
    const DatasetSpec &d = data(keV);
    const double lab = 16.00793 / 15.0001;  // E_lab / E_cm of pair 1 (masses of Pairs())
    check("E_cm[keV]: energy factor 1e-3, c.m. with pair 1's masses; the window in c.m. keV to lab MeV",
          keV.Ok() && d.columns.energy == 1e-3 && d.columns.energyCM && d.columns.lightMass == 1.00783 &&
              d.columns.heavyMass == 15.0001 && std::fabs(d.minE - 0.2 * lab) < 1e-15 && std::fabs(d.maxE - 1.4 * lab) < 1e-15);
  }
  const BuiltProject rad = BuildText(Data(Set(Set(Set(kDataset, "quantity", "    quantity dsigma_domega"), "columns",
                                                  "    columns E_lab[MeV] theta_lab[rad] value error"),
                                              "unit", "    unit mb/sr")));
  check("theta_lab[rad], mb/sr: angle factor 180/pi, value 1e-3",
        rad.Ok() && data(rad).columns.angle == 180.0 / 3.14159265358979323846 && data(rad).columns.value == 1e-3);
  ExpectOne("E_NOT_BUILT", BuildText(Data(Set(kDataset, "unit", "    unit mb") + "    covariance cov.txt\n")));
  // an extrapolation in another unit: the model in that unit
  const BuiltProject xmb = BuildText("<extrapolations>\nextrapolation x1\n    reaction 15N(p,g)16O\n    quantity sigma\n    unit mb\n"
                                     "    energies list 0.3 MeV lab\n</extrapolations>\n");
  check("an extrapolation in mb: the model x 1000", xmb.Ok() && xmb.extrapolations.extrapolations.items[0].spec.reaction.scale.factor == 1000.0);
  // a wrong dimension is the validator's (the unit is not one of the quantity's)
  ExpectOne("E_UNIT_NOT_ALLOWED", BuildText(Data(Set(kDataset, "unit", "    unit MeV*b"))));
  // the S factor: a factor on the model, defined for a charged entrance
  const BuiltProject sf = BuildText(Data(Set(Set(kDataset, "quantity", "    quantity sfactor"), "unit", "    unit keV*b")));
  check("sfactor keV*b: an S-factor factor on the model, the values times 1e-3 (to MeV b)",
        sf.Ok() && data(sf).reaction.scale.sFactor && data(sf).reaction.scale.factor == 1.0 && data(sf).columns.value == 1e-3 &&
            data(sf).reaction.observable == Observable::CrossSection);
  ExpectOne("E_REACTION_NOT_ALLOWED", BuildText(Data(Set(Set(Set(kDataset, "quantity", "    quantity sfactor"), "unit", "    unit keV*b"),
                                                         "reaction", "    reaction 16O(g,a)12C"))));
  // yields: per particle and per uC on an infinite target; a layer, a fitted efficiency and no charge state refused
  {
    const std::string inf = targets.substr(0, targets.size() - std::string("</targets>\n").size()) +
                            "target t9\n    kind infinite\n    active 15N\n    stopping formula \"a0\" a0=1e-20 MeV*cm2/atom\n"
                            "    straggling none\n</targets>\n";
    const std::string y = Set(Set(kDataset, "quantity", "    quantity yield"), "unit", "    unit 1/uC") + "    quadrature points 20\n";
    const BuiltProject q2 = BuildText(Data(y + "    target t9\n    charge_state 2\n") + inf);
    const double perUC = 1e-6 / (2 * 1.602176634e-19);
    check("yield 1/uC, q = 2: the model x 1e-6 C / (2 e)",
          q2.Ok() && std::fabs(data(q2).reaction.scale.factor / perUC - 1) < 1e-15 && !data(q2).reaction.scale.sFactor);
    const BuiltProject c = BuildText(Data(Set(y, "unit", "    unit counts/mC") + "    target t9\n    charge_state 1\n    efficiency 0.5 fixed\n") + inf);
    check("yield counts/mC, q = 1, efficiency 0.5: the model x 1e-3 C / e x 0.5",
          c.Ok() && std::fabs(data(c).reaction.scale.factor / (0.5e-3 / 1.602176634e-19) - 1) < 1e-15);
    const BuiltProject p = BuildText(Data(Set(y, "unit", "    unit 1/particle") + "    target t9\n") + inf);
    check("yield 1/particle: no factor", p.Ok() && data(p).reaction.scale.Identity());
    ExpectOne("E_NOT_BUILT", BuildText(Data(y + "    target t1\n    charge_state 2\n") + inf));
    ExpectOne("E_NOT_BUILT", BuildText(Data(Set(y, "unit", "    unit counts/uC") + "    target t9\n    charge_state 1\n"
                                                                                      "    efficiency 0.5 fit prior 0.5 3%\n") + inf));
    ExpectOne("E_FIELD_REQUIRED", BuildText(Data(y + "    target t9\n") + inf));
  }
  // targets: areal density in atoms/b, U_e in keV, a stopping formula in keV*cm2/atom
  {
    std::string t = Set(Set(targets, "areal_density", "    areal_density 1e-6 atoms/b fit prior 1e-6 2e-7"), "screening",
                        "    screening 1 keV fit prior 1 0.2");
    t = Set(t, "stopping", "    stopping formula \"a0\" a0=1e-17 keV*cm2/atom");
    const BuiltProject b = BuildText(Data(kDataset + "    target t1\n    quadrature points 20\n") + t);
    const bool one = b.Ok() && b.data.targetEffects.items.size() == 1;
    const TargetEffectSpec *e = one ? &b.data.targetEffects.items[0].spec : nullptr;
    check("target: 1e-6 atoms/b = 1e18 atoms/cm2 (prior 2e17), U_e 1 keV = 1000 eV (prior 200 eV), stopping (a0)*0.001",
          e && e->density == 1e18 && e->densitySigmaText == "2e+17" && e->screeningUe == 1000.0 && e->screeningSigma == 200.0 &&
              e->stoppingPowerEquation == "(a0)*0.001");
    if (!one)
      for (const Diagnostic &d : b.diagnostics) std::printf("        got %s: %s\n", d.code.c_str(), d.message.c_str());
  }
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
        check("the effect: a layer of 1e18 atoms/cm2, screening 1000 eV fitted with 200 eV, sigma 1 keV = 0.001 MeV, 20 points "
              "(multiplier and points per width left out: the engine's defaults)",
              e.active && e.targetIntegration && e.density == 1e18 && e.stoppingPowerEquation == "a0" && e.screening &&
                  e.screeningUe == 1000 && e.screeningVary && e.screeningSigma == 200 && e.convolution && e.sigma == 1 * 1e-3 &&
                  e.numIntegrationPoints == 20 && !e.hasResonanceWidthMultiplier && !e.hasPointsPerWidth);
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
  // a run that does not use the data (an extrapolation): the same file is a warning, and the project builds
  {
    const BuiltProject b = BuildText(Data(Set(kDataset, "file", "    file no_such_file.dat")), false, true, DataFiles::NotUsed);
    ExpectOne("W_DATA_FILE_NOT_READ", b);
    check("W_DATA_FILE_NOT_READ: a warning, the project builds, the path in 'found'",
          b.Ok() && !b.diagnostics.empty() && b.diagnostics[0].severity == Severity::Warning &&
              b.diagnostics[0].found == "no_such_file.dat" && b.data.datasets.items.size() == 1);
  }
  check("the data file read in the valid case", BuildText(Data(kDataset), false, true).diagnostics.empty());
  ExpectOne("E_DATA_FILE", BuildText(Data(Set(kDataset, "columns", "    columns E_lab[MeV] _ value error _ _")), false, true));
  UnitCases(targets);
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
  // the grid forms of an extrapolation (Grids.h): a c.m. list with the entrance pair's masses, a file, an adaptive grid
  {
    std::ofstream("build_test_grid.txt") << "# E_lab (keV)\n100\n200\n";
    const std::string x = "<extrapolations>\nextrapolation x1\n    reaction 15N(p,a)12C\n    quantity sigma\n    unit b\n";
    const BuiltProject cm = BuildText(x + "    energies list 0.3 MeV cm\n</extrapolations>\n");
    check("a c.m. list: E_lab = 0.3 (1.00783 + 15.0001)/15.0001 MeV, with pair 1's masses",
          cm.Ok() && cm.extrapolations.extrapolations.items.size() == 1 &&
              std::fabs(cm.extrapolations.extrapolations.items[0].spec.energies[0] - 0.3 * 16.00793 / 15.0001) < 1e-15);
    const BuiltProject file = BuildText(x + "    energies file \"build_test_grid.txt\" keV lab\n</extrapolations>\n");
    check("a file grid in keV: 0.1 and 0.2 MeV, the path kept for --describe",
          file.Ok() && file.extrapolations.extrapolations.items[0].spec.energies == std::vector<double>({0.1, 0.2}) &&
              file.extrapolations.extrapolations.items[0].spec.gridFile == "build_test_grid.txt");
    ExpectOne("E_DATA_FILE", BuildText(x + "    energies file \"no_such_grid.txt\" MeV lab\n</extrapolations>\n"));
    const BuiltProject adaptive =
        BuildText(x + "    energies adaptive 0.1 0.35 MeV lab tolerance 1e-3 max_points 500 start_step 0.1\n</extrapolations>\n");
    const ExtrapolationSpec &a = adaptive.extrapolations.extrapolations.items[0].spec;
    check("an adaptive grid: its start grid 0.1, 0.2, 0.30000000000000004, 0.35, tolerance and max_points",
          adaptive.Ok() && a.adaptive.on && a.adaptive.tolerance == 1e-3 && a.adaptive.maxPoints == 500 &&
              a.energies == std::vector<double>({0.1, 0.1 + 0.1, 0.1 + 2 * 0.1, 0.35}) && a.minE == 0.1 && a.maxE == 0.35);
  }
  check("the removed polarization quantity is not in the table (E_NOT_IN_LIST)",
        BuildText(Data(Set(kDataset, "quantity", "    quantity polarization"))).diagnostics.size() >= 1);

  // ratio to Rutherford and to Mott (data_kinds.md section 4.3): the spec, and the refusals of Build
  {
    const std::string rth = "dataset r1\n    reaction 12C(a,a)12C\n    quantity ratio_to_rutherford\n    unit 1\n"
                            "    file build_test_data.dat\n    columns E_lab[MeV] theta_lab[deg] value error\n"
                            "    norm 1 fit prior 1 3%\n";
    const std::string mott = Rename(Set(Set(Set(rth, "reaction", "    reaction 4He(a,a)4He\n    entrance_pair 8\n    exit_pair 8"),
                                            "quantity", "    quantity ratio_to_mott"),
                                        "columns", "    columns E_lab[MeV] theta_cm[deg] value error"),
                                    "m1");
    const std::string xr = "<extrapolations>\nextrapolation x1\n    reaction 12C(a,a)12C\n    quantity ratio_to_rutherford\n"
                           "    unit 1\n    energies range 3 4 MeV lab step 0.5\n    angles list 30 60 deg lab\n</extrapolations>\n";
    const BuiltProject b = BuildText(Data(rth + mott) + xr);
    for (const Diagnostic &d : b.diagnostics) std::printf("        %s: %s\n", d.code.c_str(), d.message.c_str());
    check("ratio_to_rutherford, ratio_to_mott and an extrapolation build", b.Ok() && b.diagnostics.empty());
    if (b.data.datasets.items.size() == 2 && b.extrapolations.extrapolations.items.size() == 1) {
      const ReactionSpec &r = b.data.datasets.items[0].spec.reaction, &m = b.data.datasets.items[1].spec.reaction;
      const ReactionSpec &x = b.extrapolations.extrapolations.items[0].spec.reaction;
      check("r1: a+12C elastic (pair 2), lab angles: differential lab over Rutherford, one term",
            r.entranceKey == 2 && r.exitKey == 2 && r.observable == Observable::DifferentialLab &&
                r.coulombReference == CoulombReference::Rutherford && r.combination == ReactionSpec::Combination::None &&
                r.components.empty());
      check("r1: norm 1 fit, 3 %", b.data.datasets.items[0].spec.varyNorm && b.data.datasets.items[0].spec.normError == 3.0);
      check("m1: a+a (pair 8), c.m. angles: differential c.m. over Mott",
            m.entranceKey == 8 && m.exitKey == 8 && m.observable == Observable::DifferentialCM &&
                m.coulombReference == CoulombReference::Mott);
      check("x1: the extrapolation over Rutherford at lab angles",
            x.observable == Observable::DifferentialLab && x.coulombReference == CoulombReference::Rutherford);
    } else {
      check("2 data sets and 1 extrapolation", false);
    }
    check("dsigma_domega has no Coulomb reference",
          BuildText(Data(kDataset)).data.datasets.items[0].spec.reaction.coulombReference == CoulombReference::None);
    // refusals: the validator's (from the labels) and Build's (from the pairs of <levels>), one diagnostic each
    ExpectOne("E_REACTION_NOT_ALLOWED", BuildText(Data(Set(rth, "reaction", "    reaction 4He(a,a)4He\n    entrance_pair 8\n    exit_pair 8"))));
    ExpectOne("E_REACTION_NOT_ALLOWED", BuildText(Data(Set(rth, "reaction", "    reaction 12C(a,a)12C[4.439]"))));
    ExpectOne("E_REACTION_NOT_ALLOWED", BuildText(Data(Set(mott, "reaction", "    reaction 12C(a,a)12C"))));
    {
      const BuiltProject n = BuildText(Data(Set(Set(mott, "entrance_pair", "    entrance_pair 9"), "exit_pair", "    exit_pair 9")));
      ExpectOne("E_REACTION_NOT_ALLOWED", n);
      check("Build: a ratio_to_mott on a pair that is not identical names ratio_to_rutherford",
            n.diagnostics.size() == 1 && n.diagnostics[0].message.find("ratio_to_rutherford") != std::string::npos &&
                n.diagnostics[0].message.find("(validation)") == std::string::npos);
    }
    {
      const std::string id = Set(rth, "reaction", "    reaction 4He(a,a)4He\n    entrance_pair 8\n    exit_pair 8");
      const BuiltProject n = BuildText(Data(id));
      check("ratio_to_rutherford with identical particles names ratio_to_mott",
            n.diagnostics.size() == 1 && n.diagnostics[0].message.find("ratio_to_mott") != std::string::npos);
    }
    // the folds the Coulomb reference does not follow (BuildEffects.cpp, CoulombRatioFolds): one E_NOT_BUILT each
    const std::string t2 = "<targets>\ntarget t2\n    kind layer\n    active 12C\n    areal_density 1e18 atoms/cm2 fixed\n"
                           "    stopping formula \"a0\" a0=1e-20 MeV*cm2/atom\n    straggling none\n";
    const std::string beam = "beam b3\n    spread 0.01 MeV lab\n";
    check("a target layer and a beam spread build",
          BuildText(Data(rth + "    target t2\n    beam b3\n    quadrature points 20\n") + t2 + beam + "</targets>\n")
              .diagnostics.empty());
    ExpectOne("E_NOT_BUILT", BuildText(Data(rth + "    target t2\n    quadrature points 20\n") + t2 +
                                       "    screening 1000 eV fixed\n</targets>\n"));
    ExpectOne("E_NOT_BUILT", BuildText(Data(rth + "    target t2\n    quadrature points 20\n") + t2 +
                                       "    doppler free_gas 300 K fixed\n</targets>\n"));
    ExpectOne("E_NOT_BUILT", BuildText(Data(rth + "    beam b3\n    detector d3\n    quadrature points 20\n") +
                                       "<targets>\n" + beam + "detector d3\n    attenuation 1 0.99 0.98\n</targets>\n"));
    ExpectOne("E_NOT_BUILT", BuildText(Data(rth + "    beam b3\n    quadrature points 20\n"
                                                  "    fold_windows 3-4 MeV auto_tolerance 1e-3\n") +
                                       "<targets>\n" + beam + "</targets>\n"));
  }

  // hybrid potentials (`potential` records, BuildPotentials): the settings member by member, units converted, the
  // defaults noted, and one case per refusal
  {
    const std::string pot =
        "potential p12C\n    pair p+12C\n    coulomb_radius 4.5 fm\n"
        "    central woods_saxon depth 20000 keV radius 3.6 fm diffuseness 0.6 fm\n    matching_radius 16 fm\n";
    const std::string gauss = "potential pN\n    pair p+15N\n    coulomb_radius 3.5 fm\n    central gaussian depth 40 MeV width 2.5 fm\n"
                              "    step 0.02 fm\n";
    auto potentials = [](const std::string &records) { return "<potentials>\n" + records + "</potentials>\n"; };
    const BuiltProject b = BuildText(potentials(pot + gauss));
    for (const Diagnostic &d : b.diagnostics) std::printf("        %s: %s\n", d.code.c_str(), d.message.c_str());
    check("two potential records build", b.Ok() && b.diagnostics.empty());
    check("one setting per pair named (7: p+12C, 1: p+15N), no default", b.potentials.pairs.size() == 2 &&
                                                                            b.potentials.pairs.count(7) && b.potentials.pairs.count(1) &&
                                                                            !b.potentials.defaults.enabled);
    if (b.potentials.pairs.size() == 2) {
      const HybridPotentialSetting &ws = b.potentials.pairs.at(7), &g = b.potentials.pairs.at(1);
      check("p12C: Woods-Saxon, 20000 keV = 20 MeV, R 3.6, a 0.6, R_C 4.5, r_m 16 fm, the default step",
            ws.enabled && ws.shape == azure2::physics::FormShape::WoodsSaxon && ws.V0 == 20000 * 1e-3 && ws.R == 3.6 &&
                ws.a == 0.6 && ws.coulomb_radius_fm == 4.5 && ws.matching_radius_fm == 16.0 && ws.step_fm == 0.0275);
      check("pN: Gaussian 40 MeV, r0 2.5 fm, R_C 3.5 fm, the default r_m, step 0.02 fm",
            g.enabled && g.shape == azure2::physics::FormShape::Gaussian && g.V0 == 40.0 && g.r0 == 2.5 &&
                g.coulomb_radius_fm == 3.5 && g.matching_radius_fm == 15.0 && g.step_fm == 0.02);
      check("the potential values carry R_C, r_m and the step", ws.Potential().coulomb_radius_fm == 4.5 &&
                                                                 ws.Potential().matching_radius_fm == 16.0 &&
                                                                 g.Potential().step_fm == 0.02);
      check("the record of each pair", b.potentialRecords.at(7).first == "p12C" && b.potentialRecords.at(1).first == "pN");
    }
    check("the defaults used are noted with the table's values",
          b.potentialNotes == std::vector<std::string>{"potential p12C: step not given, 0.0275 fm (the default)",
                                                       "potential pN: matching_radius not given, 15 fm (the default)"});
    // the table's default texts are the solver's defaults
    const RecordKind *kind = Schema::FindKind("potential");
    check("matching_radius and step defaults of the table = HybridPotentialSetting's",
          kind && kind->Find("matching_radius")->defaultText == "15 fm" && kind->Find("step")->defaultText == "0.0275 fm" &&
              HybridPotentialSetting().matching_radius_fm == 15.0 && HybridPotentialSetting().step_fm == 0.0275);
    const std::string dd = "potential d\n    pair d+12C\n    coulomb_radius 3.5 fm\n    central gaussian depth 40 MeV width 2.5 fm\n";
    check("pair_key names one of the two d+12C pairs",
          BuildText(potentials(dd + "    pair_key 6\n")).potentials.pairs.count(6) == 1);
    ExpectOne("E_UNKNOWN_PAIR", BuildText(potentials(Set(pot, "pair", "    pair p+13C"))));
    ExpectOne("E_AMBIGUOUS_PAIR", BuildText(potentials(dd)));
    ExpectOne("E_PAIR_KEY", BuildText(potentials(pot + "    pair_key 7\n")));
    ExpectOne("E_PAIR_KEY", BuildText(potentials(dd + "    pair_key 1\n")));
    ExpectOne("E_POTENTIAL_TWICE", BuildText(potentials(pot + "potential again\n" + pot.substr(pot.find('\n') + 1))));
    ExpectOne("E_POTENTIAL_PAIR_KIND", BuildText(potentials(Set(pot, "pair", "    pair g+16O"))));
    ExpectOne("E_POTENTIAL_RELATIVISTIC", BuildText(potentials(pot), true));
    // the THM data set d5 above enters through p+15N: a potential on that pair refuses it
    const std::string thm = Set(Set(kDataset, "norm", "    norm profiled"), "unit", "    unit arbitrary") + "    thm global\n";
    ExpectOne("E_THM_HYBRID_POTENTIAL", BuildText(Data(thm) + potentials(gauss)));
    check("a THM data set with a potential on another pair builds", BuildText(Data(thm) + potentials(pot)).diagnostics.empty());
    // the validator's refusals of the record (no build check needed)
    ExpectOne("E_FIELD_REQUIRED", BuildText(potentials(Set(pot, "coulomb_radius", ""))));
    ExpectOne("E_OUT_OF_RANGE", BuildText(potentials(Set(pot, "coulomb_radius", "    coulomb_radius 0 fm"))));
    ExpectOne("E_OUT_OF_RANGE", BuildText(potentials(Set(gauss, "step", "    step -0.01 fm"))));
    ExpectOne("E_OUT_OF_RANGE", BuildText(potentials(Set(gauss, "central", "    central gaussian depth 40 MeV width 0 fm"))));
    ExpectOne("E_UNIT_NOT_ALLOWED", BuildText(potentials(Set(pot, "matching_radius", "    matching_radius 16 MeV"))));
    ExpectOne("E_PART_MISSING", BuildText(potentials(Set(pot, "central", "    central woods_saxon depth 20 MeV radius 3.6 fm"))));
    // a spin-orbit term (spin_orbit): the term member by member (keV converted), the potential value carries it;
    // one case per refusal
    {
      const std::string so = "    spin_orbit woods_saxon depth 6000 keV radius 3.5 fm diffuseness 0.65 fm on light operator 2l.s\n";
      const BuiltProject s = BuildText(potentials(pot + so));
      for (const Diagnostic &d : s.diagnostics) std::printf("        %s: %s\n", d.code.c_str(), d.message.c_str());
      const bool built = s.Ok() && s.diagnostics.empty() && s.potentials.pairs.count(7) && s.potentials.pairs.at(7).spin_orbit;
      check("p12C with a spin-orbit term builds", built);
      if (built) {
        const azure2::physics::SpinOrbitTerm &t = *s.potentials.pairs.at(7).spin_orbit;
        check("spin_orbit: 6000 keV = 6 MeV, R_so 3.5 fm, a_so 0.65 fm, on the light particle",
              t.depth_MeV == 6000 * 1e-3 && t.form.radius_fm == 3.5 && t.form.diffuseness_fm == 0.65 &&
                  t.acts_on == azure2::physics::SpinOrbitParticle::Light && t.form.shape == azure2::physics::FormShape::WoodsSaxon);
        check("the potential value carries the term", s.potentials.pairs.at(7).Potential().spin_orbit == t);
      }
      check("a record without spin_orbit has none", !b.potentials.pairs.at(7).spin_orbit);
      ExpectOne("E_SO_ZERO_SPIN", BuildText(potentials(pot + Set(so, "spin_orbit",
                                                                 "    spin_orbit woods_saxon depth 6 MeV radius 3.5 fm diffuseness 0.65 fm on heavy operator 2l.s"))));
      const std::string soN = "    spin_orbit woods_saxon depth 6 MeV radius 3.5 fm diffuseness 0.65 fm on light operator 2l.s\n";
      ExpectOne("E_SO_TWO_SPINS", BuildText(potentials(gauss + soN)));
      ExpectOne("E_SO_RELATIVISTIC", BuildText(potentials(pot + so), true));
      // THM through p+15N with the spin of 15N set to 0 (one spinning particle): E_SO_THM in place of E_THM_HYBRID_POTENTIAL
      LevelTables oneSpin = Pairs();
      oneSpin.pairs[0].twoJHeavy = 0;
      ExpectOne("E_SO_THM", BuildText(Data(thm) + potentials(gauss + soN), false, false, DataFiles::Required, oneSpin));
      ExpectOne("E_THM_HYBRID_POTENTIAL", BuildText(Data(thm) + potentials(gauss), false, false, DataFiles::Required, oneSpin));
      // the validator's: the operator is 2l.s only (an l.s strength is half), the shape woods_saxon only
      ExpectOne("E_NOT_IN_LIST", BuildText(potentials(pot + Set(so, "spin_orbit",
                                                         "    spin_orbit woods_saxon depth 6 MeV radius 3.5 fm diffuseness 0.65 fm on light operator l.s"))));
      ExpectOne("E_PART_MISSING", BuildText(potentials(pot + Set(so, "spin_orbit",
                                                                 "    spin_orbit woods_saxon depth 6 MeV radius 3.5 fm diffuseness 0.65 fm on light"))));
      ExpectOne("E_OUT_OF_RANGE", BuildText(potentials(pot + Set(so, "spin_orbit",
                                                                 "    spin_orbit woods_saxon depth 0 MeV radius 3.5 fm diffuseness 0.65 fm on light operator 2l.s"))));
    }
    // the refusals shared with pyazr and the API (SpinOrbitChecks.h)
    {
      using azure2::physics::SpinOrbitParticle;
      check("one spinning particle, the term on it: allowed",
            !SpinOrbitPairRefusal(1, 0, SpinOrbitParticle::Light, "p+12C") && !SpinOrbitPairRefusal(0, 3, SpinOrbitParticle::Heavy, "a+7Li"));
      check("E_SO_TWO_SPINS: both spins non-zero", SpinOrbitPairRefusal(1, 1, SpinOrbitParticle::Light, "p+15N").code == "E_SO_TWO_SPINS");
      check("E_SO_ZERO_SPIN: the term on the spin-0 particle", SpinOrbitPairRefusal(1, 0, SpinOrbitParticle::Heavy, "p+12C").code == "E_SO_ZERO_SPIN");
      check("E_SO_NEEDS_J: channel functions of a spin-orbit pair without J", ChannelFunctionsJRefusal(true, false, "1").code == "E_SO_NEEDS_J");
      check("E_SO_NO_TERM: a J for a pair without a spin-orbit term", ChannelFunctionsJRefusal(false, true, "1").code == "E_SO_NO_TERM");
      check("J with a spin-orbit pair, none without: allowed", !ChannelFunctionsJRefusal(true, true, "1") && !ChannelFunctionsJRefusal(false, false, "1"));
    }
    // the old block in a format-3 project (the reading path: ProjectReader)
    {
      Format3ReadOptions read;
      read.build = false;
      const Format3Project p = ReadFormat3Text("<config>\nformat=3\n</config>\n<potential>\nuseHybridPotential=1\n</potential>\n", read);
      const bool one = p.diagnostics.size() == 1 && p.diagnostics[0].code == "E_OLD_POTENTIAL_BLOCK" &&
                       p.diagnostics[0].line == 4 && p.diagnostics[0].message.find("--convert") != std::string::npos;
      check("E_OLD_POTENTIAL_BLOCK alone, at the block's line, naming the converter", one);
      if (!one)
        for (const Diagnostic &d : p.diagnostics) std::printf("        got %s: %s\n", d.code.c_str(), d.message.c_str());
    }
  }

  std::printf("experiment_build: %d checks, %d failed\n", checks, failures);
  return failures == 0 ? 0 : 1;
}
