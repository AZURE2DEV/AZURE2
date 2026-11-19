// The front ends' facade without Qt (docs/upgrade/qt_gui.md section 10): front_test <source dir> <AZURE2 binary>.
//  (a) the tables: schema, commands, options, fitters, <config> keys;
//  (b) the test projects T1, T2, T3 and 13N opened and written back byte for byte, also after every record, every
//      <levels> row and every <config> entry is set to what it holds; an edit changes only its lines and is undone;
//      a record added, renamed and removed; --validate's diagnostics, a broken field named by its code;
//  (c) an old project refused on opening and converted by the converter of --convert;
//  (d) a run of `calculate` through the facade writes the same output files, byte for byte, as the program on the same
//      project, and returns the chi^2 the program prints; T1 and 13N;
//  (e) a fit stopped through the facade ends at its best point with exit 0 and fitStatus.out saying "stopped".
// Each check prints PASS or FAIL; exit 1 on a failure.
#include "azure2/front/Engine.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;
namespace front = azure2::front;

namespace {

int failures = 0;

void Check(bool ok, const std::string &what) {
  std::cout << (ok ? "PASS " : "FAIL ") << what << std::endl;
  if (!ok) failures++;
}

std::string Bytes(const fs::path &p) {
  std::ifstream in(p, std::ios::binary);
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

/// A copy of a project's directory (the project file, its data and side files) under `to`.
fs::path CopyProject(const fs::path &projectFile, const fs::path &to) {
  fs::remove_all(to);
  fs::create_directories(to);
  fs::copy(projectFile.parent_path(), to, fs::copy_options::recursive | fs::copy_options::skip_symlinks);
  for (const char *d : {"output", "checks"}) {
    fs::remove_all(to / d);
    fs::create_directories(to / d);
  }
  return to / projectFile.filename();
}

class Collect : public front::RunObserver {
 public:
  void Message(const std::string &line, bool transient) override {
    std::lock_guard<std::mutex> lock(mutex);
    if (!transient) lines.push_back(line);
  }
  void Progress(const front::RunProgress &p) override {
    std::lock_guard<std::mutex> lock(mutex);
    progress++;
    if (stopAfter > 0 && progress >= stopAfter && stop) stop->store(true);
    last = p;
  }
  std::mutex mutex;
  std::vector<std::string> lines;
  long progress = 0;
  long stopAfter = 0;
  std::atomic<bool> *stop = nullptr;
  front::RunProgress last;
};

/// The output files of two runs are the same files with the same bytes.
bool SameOutputs(const fs::path &a, const fs::path &b, std::string &why) {
  std::vector<std::string> names;
  for (const auto &e : fs::directory_iterator(a))
    if (e.is_regular_file()) names.push_back(e.path().filename().string());
  int nb = 0;
  for (const auto &e : fs::directory_iterator(b))
    if (e.is_regular_file()) nb++;
  if ((int)names.size() != nb) {
    why = std::to_string(names.size()) + " files against " + std::to_string(nb);
    return false;
  }
  if (names.empty()) {
    why = "no output file";
    return false;
  }
  for (const std::string &n : names)
    if (!fs::exists(b / n) || Bytes(a / n) != Bytes(b / n)) {
      why = n + " differs";
      return false;
    }
  why = std::to_string(names.size()) + " files identical";
  return true;
}

void RoundTrip(const fs::path &file) {
  const std::string name = file.filename().string();
  front::Project p;
  std::string why;
  Check(front::Project::Open(file.string(), p, why), name + ": opens (" + why + ")");
  const std::string original = Bytes(file);
  Check(p.Text() == original && !p.Modified(), name + ": text as read");
  // every record, row and entry set to what it holds: not one byte moves
  for (const std::string &section : p.Sections())
    for (const front::Record &r : p.Records(section))
      if (!p.SetRecord(section, r, why)) Check(false, name + ": SetRecord " + r.name + ": " + why);
  if (!p.SetLevels(p.Levels(), why)) Check(false, name + ": SetLevels: " + why);
  for (const front::ConfigEntry &e : p.Config())
    if (e.name != "format" && !p.SetConfig(e.name, e.value, why)) Check(false, name + ": SetConfig " + e.name + ": " + why);
  Check(p.Text() == original, name + ": every record, row and entry set as it stands: the same bytes");
  const fs::path out = fs::path("front_test_work") / ("roundtrip_" + name);
  fs::create_directories(out.parent_path());
  Check(p.Save(out.string(), why) && Bytes(out) == original, name + ": saved byte for byte");
}


/// The value text with one blank between tokens outside quotes (what a writer of format 3 writes).
std::string OneBlank(const std::string &text) {
  std::string out;
  bool quote = false, blank = false;
  for (char c : text) {
    if (c == '"') quote = !quote;
    if (!quote && (c == ' ' || c == '\t')) {
      blank = true;
      continue;
    }
    if (blank && !out.empty()) out += ' ';
    blank = false;
    out += c;
  }
  return out;
}

/// (f) the forms of values (docs/upgrade/qt_gui.md section 11): every field of every record of the test projects
/// split into its pieces by the engine's parser and written back from them gives the same text; the summary
/// columns, the include switch and the groups come from the schema; pairs, reactions and the move of a record.
void FormsChecks(const fs::path &source) {
  const std::vector<fs::path> files = {source / "tests/upgrade/T1/T1.azr", source / "tests/upgrade/T2/T2.azr",
                                       source / "tests/upgrade/T3/T3.azr", source / "tests/13N/13N.azr",
                                       source / "experiment/tests/data/section_examples.azr"};
  for (const fs::path &file : files) {
    front::Project p;
    std::string why;
    if (!front::Project::Open(file.string(), p, why)) {
      Check(false, file.filename().string() + ": opens");
      continue;
    }
    int fields = 0, same = 0;
    std::vector<std::string> differ;
    for (const auto &section : front::RecordSections())
      for (const front::Record &r : p.Records(section.first))
        for (const front::Field &f : r.fields) {
          front::ValueParts parts;
          fields++;
          if (!front::ParseFieldValue(r.kind, f.key, f.value, r.Value("quantity"), parts, why)) {
            differ.push_back(r.name + "." + f.key + " not read: " + why);
            continue;
          }
          const std::string written = front::FieldValueText(r.kind, f.key, parts);
          if (written == OneBlank(f.value))
            same++;
          else
            differ.push_back(r.name + "." + f.key + ": '" + f.value + "' -> '" + written + "'");
        }
    for (std::size_t i = 0; i < differ.size() && i < 8; i++) std::cout << "    " << differ[i] << std::endl;
    Check(fields > 0 && same == fields, file.filename().string() + ": " + std::to_string(same) + " of " +
                                            std::to_string(fields) + " field values split and written back unchanged");
  }
  // a piece changed: only it moves
  {
    front::ValueParts n;
    std::string why;
    const bool read = front::ParseFieldValue("target", "areal_density", "2e+17 atoms/cm2 fit prior 2e+17 10% grid_max 3e+17", "", n, why);
    Check(read && n.fit && n.hasPrior && n.priorPercent && n.priorWidth == "10" && n.gridMax == "3e+17" && n.unit == "atoms/cm2",
          "a fitted areal density with a prior and grid_max split into its pieces");
    n.fit = false;
    n.hasPrior = false;
    Check(front::FieldValueText("target", "areal_density", n) == "2e+17 atoms/cm2 fixed grid_max 3e+17",
          "the prior dropped and the mode fixed: '" + front::FieldValueText("target", "areal_density", n) + "'");
    front::ValueParts s;
    front::ParseFieldValue("target", "stopping", "formula \"a0*x^a1\" a0=6.0e-20 a1=-0.30 MeV*cm2/atom", "", s, why);
    Check(s.alternative == "formula" && s.Part("expression") && s.Part("expression")->expression == "a0*x^a1" &&
              s.Part("expression")->params.size() == 2 && s.Part("unit") && s.Part("unit")->text == "MeV*cm2/atom",
          "a stopping formula split into the alternative, the expression, its parameters and the unit");
    front::ValueParts b;
    front::ParseFieldValue("beam", "profile", "skewed_gaussian xi 0.34 MeV omega 0.035 MeV alpha 1.5 weight 1.0", "", b, why);
    Check(b.alternative == "skewed_gaussian" && b.Part("omega") && b.Part("omega")->number == "0.035" &&
              b.Part("alpha") && b.Part("alpha")->text == "1.5",
          "a beam profile split into its named parameters");
    const bool unitRefused = !front::ParseFieldValue("dataset", "unit", "barn_per_furlong", "sigma", b, why);
    Check(unitRefused && !why.empty() && why.find("field '") == std::string::npos, "a wrong unit refused in plain words: " + why);
    const bool modeRefused = !front::ParseFieldValue("target", "areal_density", "2e17 atoms/cm2", "", b, why);
    Check(modeRefused, "a nuisance without its mode refused: " + why);
    Check(front::ShortValueText("dataset", "norm", "1 fit prior 1 5%") == "1 (fit)", "the short text of a fitted norm");
    Check(front::ShortValueText("detector", "attenuation", "1 0.995  0.98") == "1 0.995 0.98", "the short text of a number list");
    // the beam profile's preview: the kernel's weight; a skewed Gaussian of alpha > 0 peaks above xi, and its
    // integral is the component's weight (the curve spans xi +- 5 omega)
    front::Record beam;
    beam.kind = "beam";
    beam.name = "b";
    beam.fields = {{"profile", "skewed_gaussian xi 340 keV omega 0.035 MeV alpha 1.5 weight 2", ""}, {"truncation", "0", ""}};
    front::ShapeCurve shape;
    const bool drawn = front::PreviewShape(beam, 401, shape, why);
    double area = 0, peak = 0, at = 0;
    for (std::size_t i = 1; i < shape.x.size(); i++) area += 0.5 * (shape.y[i] + shape.y[i - 1]) * (shape.x[i] - shape.x[i - 1]);
    for (std::size_t i = 0; i < shape.x.size(); i++)
      if (shape.y[i] > peak) {
        peak = shape.y[i];
        at = shape.x[i];
      }
    Check(drawn && std::fabs(area - 2.0) < 2e-3 && at > 0.34 && at < 0.34 + 0.035 && std::fabs(shape.x.front() - (0.34 - 0.175)) < 1e-12,
          "the beam profile's preview: area " + std::to_string(area) + " (weight 2), peak at " + std::to_string(at) + " MeV");
    front::Record target;
    target.kind = "target";
    Check(!front::PreviewShape(target, 10, shape, why), "no preview for a record without a shape (" + why + ")");
  }
  // the shapes from the schema
  {
    const front::KindInfo *dataset = nullptr, *target = nullptr, *beam = nullptr;
    const std::vector<front::KindInfo> kinds = front::RecordKinds();
    for (const front::KindInfo &k : kinds) {
      if (k.kind == "dataset") dataset = &k;
      if (k.kind == "target") target = &k;
      if (k.kind == "beam") beam = &k;
    }
    auto has = [](const front::KindInfo *k, const std::string &key) {
      return k && std::find(k->summaryFields.begin(), k->summaryFields.end(), key) != k->summaryFields.end();
    };
    std::string listed;
    if (dataset) for (const std::string &x : dataset->summaryFields) listed += x + " ";
    Check(dataset && dataset->summaryFields.size() >= 5 && dataset->summaryFields[0] == "reaction" && has(dataset, "norm") &&
              has(dataset, "target") && has(dataset, "energy_window") && !has(dataset, "columns") && !has(dataset, "source") && !has(dataset, "quadrature"),
          "a data set's summary columns: " + listed);
    Check(dataset && dataset->includeField == "active" && dataset->includeOn == "yes", "a data set's include switch is 'active'");
    Check(has(target, "stopping") && has(target, "areal_density") && has(beam, "profile"), "a target's and a beam's summary columns");
    std::string phase, truncation;
    if (target) for (const front::FieldInfo &f : target->fields) if (f.key == "phase") phase = f.group;
    if (beam) for (const front::FieldInfo &f : beam->fields) if (f.key == "truncation") truncation = f.group;
    Check(phase == "stopping" && truncation == "profile", "a stopping's phase grouped with the stopping, a profile's truncation with the profile");
    bool profileShape = false;
    if (beam)
      for (const front::FieldInfo &f : beam->fields)
        if (f.key == "profile")
          profileShape = f.value.alternatives.size() == 1 && f.value.alternatives[0].parts.size() == 4 &&
                         !f.value.alternatives[0].parts[0].value.units.empty();
    Check(profileShape, "the beam profile's shape: one alternative, four named parts with their units");
    const auto columns = front::DefaultColumns("dsigma_domega");
    Check(columns.size() == 4 && !columns[0].second.empty() && !columns[1].second.empty() && columns[3].second.empty(),
          "the columns of a new differential cross section: energy, angle, value, error");
    Check(!front::Quantities("dataset").empty() && !front::Quantities("dataset")[0].title.empty(), "the quantities with their titles");
  }
  // every grid kind the facade offers is one the engine reads
  {
    std::string why;
    int read = 0;
    const std::vector<front::GridKind> kinds = front::GridKinds(true);
    for (const front::GridKind &k : kinds) {
      front::ValueParts g;
      g.type = "grid";
      g.gridKind = k.word;
      g.lo = "0.1";
      g.hi = "1";
      g.unit = "MeV";
      g.frame = front::Frames()[1];
      g.step = "0.01";
      g.numbers = {"0.1", "0.2"};
      g.text = "grid.txt";
      g.tolerance = "1e-3";
      g.maxPoints = "100";
      g.startStep = "0.01";
      g.unit = k.unit ? "deg" : "";
      front::ValueParts back;
      const std::string angles = front::FieldValueText("extrapolation", "angles", g);
      read += front::ParseFieldValue("extrapolation", "angles", angles, "", back, why) && back.gridKind == k.word;
    }
    Check(read == (int)kinds.size(), std::to_string(read) + " of " + std::to_string(kinds.size()) + " grid kinds written and read back");
  }
  // pairs, reactions, the move of a record
  {
    front::Project p;
    std::string why;
    front::Project::Open((source / "tests/upgrade/T1/T1.azr").string(), p, why);
    const std::vector<front::PairChoice> pairs = front::Pairs(p);
    Check(!pairs.empty() && pairs[0].key == 1 && !pairs[0].label.empty(), "T1: the pairs of <levels> (" + (pairs.empty() ? std::string() : pairs[0].label) + ")");
    const std::string reaction = p.Records("data")[0].Value("reaction");
    int in = 0, out = 0;
    const bool resolved = front::ReactionPairs(p, reaction, in, out);
    Check(resolved && in > 0 && out > 0 && front::ReactionLabel(p, in, out) == reaction,
          "T1: the reaction " + reaction + " is the pairs " + std::to_string(in) + " -> " + std::to_string(out) + " and back");
    const std::string original = p.Text();
    const std::vector<front::Record> before = p.Records("data");
    Check(p.MoveRecord("data", before[1].name, -1, why) && p.Records("data")[0].name == before[1].name &&
              p.Records("data")[1].name == before[0].name,
          "T1: a data set moved up");
    Check(p.MoveRecord("data", before[1].name, +1, why) && p.Text() == original, "T1: and down again: the same bytes");
    const bool refused = !p.MoveRecord("data", before[0].name, -1, why);
    Check(refused, "T1: the first record cannot move up (" + why + ")");
  }
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "usage: front_test <source dir> <AZURE2 binary>" << std::endl;
    return 2;
  }
  const fs::path source = argv[1];
  const std::string binary = argv[2];
  fs::create_directories("front_test_work");

  // (a) the tables
  Check(front::SchemaJson().find("\"dataset\"") != std::string::npos, "schema JSON names the record kind dataset");
  bool calculate = false, fit = false;
  for (const front::Command &c : front::Commands()) {
    calculate = calculate || c.name == "calculate";
    fit = fit || c.name == "fit";
  }
  Check(calculate && fit, "the commands calculate and fit");
  Check(!front::Options().empty() && !front::Fitters().empty(), "the options and the fitters");
  Check(front::ConfigTableJson().find("\"use_brune\"") != std::string::npos, "the <config> table names use_brune");
  {
    bool reaction = false;
    for (const front::KindInfo &k : front::RecordKinds())
      if (k.kind == "dataset" && !k.fields.empty() && k.fields[0].key == "reaction") reaction = k.fields[0].required;
    Check(reaction, "the record kinds in the schema's order: dataset's first field the required reaction");
    const std::vector<std::string> units = front::QuantityUnits("sfactor");
    Check(!units.empty() && units[0] == "MeV*b", "the units of the S factor");
  }

  // (b) round trips and edits
  const std::vector<fs::path> projects = {source / "tests/upgrade/T1/T1.azr", source / "tests/upgrade/T2/T2.azr",
                                          source / "tests/upgrade/T3/T3.azr", source / "tests/13N/13N.azr"};
  for (const fs::path &f : projects) RoundTrip(f);
  {
    front::Project p;
    std::string why;
    front::Project::Open(projects[0].string(), p, why);
    const std::string original = p.Text();
    std::vector<front::Record> data = p.Records("data");
    Check(data.size() > 3 && data[0].kind == "dataset", "T1: the data sets");
    front::Record r = data[0];
    for (front::Field &f : r.fields)
      if (f.key == "norm") f.value = "1.05 fit";
    Check(p.SetRecord("data", r, why), "T1: a norm set to 1.05 fit");
    // exactly one line differs
    std::istringstream a(original), b(p.Text());
    std::string la, lb;
    int differ = 0, lines = 0;
    while (std::getline(a, la) && std::getline(b, lb)) {
      lines++;
      if (la != lb) differ++;
    }
    Check(differ == 1, "T1: the edit changed one line of " + std::to_string(lines));
    Check(p.Records("data")[0].Value("norm") == "1.05 fit", "T1: the record reads back the new value");
    r.fields.push_back({"energy_shift", "0 keV fixed", ""});
    Check(p.SetRecord("data", r, why) && p.Records("data")[0].Value("energy_shift") == "0 keV fixed",
          "T1: a field added (" + why + ")");
    front::Record copy = data[1];
    copy.name = "copy_of_" + data[1].name;
    Check(p.SetRecord("data", copy, why) && p.Records("data").size() == data.size() + 1, "T1: a record appended");
    Check(p.RenameRecord("data", copy.name, "renamed", why) && p.Records("data").back().name == "renamed",
          "T1: a record renamed");
    Check(p.RemoveRecord("data", "renamed", why) && p.Records("data").size() == data.size(), "T1: a record removed");
    Check(p.SetRecord("data", data[0], why) && p.Text() == original, "T1: the first record set back: the same bytes");
    front::Validation v = p.Validate();
    Check(v.ok && v.format3, "T1: valid (" + std::to_string(v.diagnostics.size()) + " diagnostics)");
    // a broken field: the engine's code, the record and the field
    front::Record bad = data[0];
    for (front::Field &f : bad.fields)
      if (f.key == "quantity") f.value = "sigmaa";
    p.SetRecord("data", bad, why);
    v = p.Validate();
    bool named = false;
    for (const front::Diagnostic &d : v.diagnostics)
      named = named || (d.severity == "error" && d.record == bad.name && d.field == "quantity" && !d.code.empty());
    Check(!v.ok && named, "T1: quantity sigmaa refused with a code naming the record and the field");
    // <levels>: one width changed through the row writer
    p = front::Project::FromText(original, projects[0].string());
    std::vector<front::LevelRow> rows = p.Levels();
    Check(!rows.empty() && rows[0].tokens.size() >= 31, "T1: the <levels> rows");
    rows[0].tokens[11] = "1.9";
    Check(p.SetLevels(rows, why) && p.Levels()[0].tokens[11] == "1.9", "T1: a level width written by the row writer");
    rows[0].tokens[11] = "x";
    Check(!p.SetLevels(rows, why), "T1: a row the engine cannot read is refused (" + why + ")");
    // <config>
    Check(p.SetConfig("use_brune", "0", why) && p.ConfigValue("use_brune") == "0", "T1: <config> use_brune=0");
    Check(!p.SetConfig("no_such_key", "1", why), "T1: an unknown <config> key refused");
  }
  {
    // an empty project: what New makes, valid once saved
    front::Project e = front::Project::Empty((fs::path("front_test_work") / "empty.azr").string());
    std::string why;
    Check(e.Save("", why), "an empty project saved");
    front::Project back;
    Check(front::Project::Open(e.Path(), back, why), "the empty project opens (" + why + ")");
    Check(back.ConfigValue("format") == "3", "the empty project is format 3");
  }

  // (c) an old project: refused, then converted
  {
    const fs::path old = source / "tests/convert/fixtures/tests/13N/13N.azr";
    front::Project p;
    std::string why;
    Check(!front::Project::Open(old.string(), p, why) && why.find("--convert") != std::string::npos,
          "an old project refused with the conversion command");
    const fs::path dir = fs::path("front_test_work") / "convert";
    fs::remove_all(dir);
    fs::create_directories(dir);
    const fs::path to = dir / "13N.f3.azr";
    std::string report;
    Check(front::ConvertProject(old.string(), to.string(), report) && front::Project::Open(to.string(), p, why),
          "the old project converted and opened");
    Check(!front::ConvertProject(old.string(), to.string(), report), "the converter overwrites nothing");
  }

  // (d) calculate through the facade = the program
  for (const fs::path &f : {projects[0], projects[3]}) {
    const std::string name = f.stem().string();
    const fs::path a = CopyProject(f, fs::path("front_test_work") / ("cli_" + name));
    const fs::path b = CopyProject(f, fs::path("front_test_work") / ("front_" + name));
    std::string command = "\"" + binary + "\" calculate --threads=2 \"" + a.string() + "\" > \"" +
                          (a.parent_path() / "stdout.txt").string() + "\" 2>&1";
#ifdef _WIN32
    // cmd.exe drops the first and the last quote of a line that starts with one: an outer pair keeps the inner ones.
    command = "\"" + command + "\"";
#endif
    const int status = std::system(command.c_str());
    Collect observer;
    front::RunRequest request;
    request.command = "calculate";
    request.project = b.string();
    request.threads = 2;
    const front::RunResult result = front::Run(request, observer, nullptr);
    Check(status == 0 && result.exitCode == 0, name + ": both runs end with exit 0");
    std::string why;
    const bool same = SameOutputs(a.parent_path() / "output", b.parent_path() / "output", why);
    Check(same, name + ": output files of the program and of the facade: " + why);
    Check(!result.outputFiles.empty(), name + ": the run lists " + std::to_string(result.outputFiles.size()) + " files");
    // the chi^2 as data equals the program's line
    const std::string out = Bytes(a.parent_path() / "stdout.txt");
    const std::size_t at = out.find("Total Chi-Squared: ");
    const double printed = at == std::string::npos ? NAN : std::atof(out.c_str() + at + 19);
    Check(result.haveChiSquared && std::fabs(result.totalChiSquared - printed) <= 1e-11 * std::fabs(printed),
          name + ": total chi^2 as data " + std::to_string(result.totalChiSquared) + " = printed " + std::to_string(printed));
    Check(!result.dataSets.empty(), name + ": " + std::to_string(result.dataSets.size()) + " data sets' chi^2");
    bool sawTotal = false;
    for (const std::string &l : observer.lines) sawTotal = sawTotal || l.rfind("Total Chi-Squared: ", 0) == 0;
    Check(sawTotal, name + ": the messages carry the program's lines");

    // (f) the run's curves as data: every row of every AZUREOut_*.out file in one curve, each data set's points
    // as many as chiSquared.out's N
    front::Project opened;
    front::Project::Open(b.string(), opened, why);
    const front::Curves curves = front::ReadCurves(opened);
    std::size_t rows = 0, points = 0;
    for (const auto &e : fs::directory_iterator(b.parent_path() / "output")) {
      const std::string file = e.path().filename().string();
      if (file.rfind("AZUREOut_aa=", 0) != 0 || e.path().extension() != ".out") continue;
      std::ifstream in(e.path());
      std::string line;
      while (std::getline(in, line)) rows += line.find_first_not_of(" \t") != std::string::npos;
    }
    for (const front::Curve &c : curves.curves) points += c.points.size() + c.coefficients.size();
    std::size_t chiN = 0;
    for (const front::DataSetChiSquared &d : result.dataSets) chiN += d.points;
    std::string notes;
    for (const std::string &n : curves.notes) notes += " " + n;
    Check(!curves.curves.empty() && points == rows && points == chiN,
          name + ": the curves hold " + std::to_string(points) + " points, the files " + std::to_string(rows) +
              " rows, chiSquared.out " + std::to_string(chiN) + notes);
    bool named = true;
    for (const front::Curve &c : curves.curves) named = named && !c.record.empty() && c.section == "data";
    Check(named, name + ": every curve names its data set");
  }

  // (g) the fit's parameters, <parameterSettings> and <mcmc> rows through their grammars
  {
    const fs::path b = CopyProject(projects[0], fs::path("front_test_work") / "params_T1");
    std::vector<front::Parameter> params;
    std::string why;
    const bool read = front::ReadParameters(b.string(), params, why);
    Check(read && !params.empty(),
          "T1: " + std::to_string(params.size()) + " parameters of a session (" + why + ")");
    bool energy = false;
    for (const front::Parameter &p : params) energy = energy || (p.name == "energy_1" && p.kind == "energy");
    Check(energy, "T1: energy_1 among them, of kind energy");
    front::Project p;
    front::Project::Open(b.string(), p, why);
    const std::string original = p.Text();
    front::ParameterSetting row;
    row.name = "energy_1";
    row.nominal = 0.5;
    row.lower = 0.4;
    row.upper = 0.6;
    Check(front::SetParameterSetting(p, row, why) && front::ParameterSettings(p).size() == 1 &&
              front::ParameterSettings(p)[0].upper == 0.6,
          "T1: a <parameterSettings> row written and read back");
    Check(p.Validate().ok, "T1: valid with the row");
    Check(front::SetParameterSetting(p, row, why, true) && front::ParameterSettings(p).empty(), "T1: the row removed");
    front::SamplerPrior prior{"energy_1", 0.5, 0.5, 0.01, true};
    front::SetSamplerPriors(p, {prior});
    const std::vector<front::SamplerPrior> back = front::SamplerPriors(p);
    Check(back.size() == 1 && back[0].name == "energy_1" && back[0].sigma == 0.01 && back[0].use,
          "T1: an <mcmc> prior row written and read back");
    front::SetSamplerPriors(p, {});
    Check(front::SamplerPriors(p).empty(), "T1: the <mcmc> prior rows cleared");
    std::vector<front::ParameterValue> values;
    Check(!front::ReadParameterFile((b.parent_path() / "none.sav").string(), values, why), "a missing parameter file");
    (void)original;
  }

  // (e) a fit stopped through the facade
  {
    const fs::path b = CopyProject(source / "tests/upgrade/T1/T1_start.azr", fs::path("front_test_work") / "stop");
    std::atomic<bool> stop{false};
    Collect observer;
    observer.stop = &stop;
    observer.stopAfter = 25;
    front::RunRequest request;
    request.command = "fit";
    request.project = b.string();
    request.threads = 2;
    const front::RunResult result = front::Run(request, observer, &stop);
    const std::string status = Bytes(b.parent_path() / "output" / "fitStatus.out");
    Check(result.stopped && result.exitCode == 0, "T1_start: the stopped fit ends with exit 0");
    Check(status.find("stopped") != std::string::npos, "T1_start: fitStatus.out says stopped");
    Check(observer.progress >= 25 && observer.progress < 2000, "T1_start: " + std::to_string(observer.progress) +
                                                                   " progress records, the fit ended soon after the request");
    Check(result.haveChiSquared, "T1_start: the final evaluation at the best point");
  }

  FormsChecks(source);

  std::cout << (failures ? "FAILED: " + std::to_string(failures) : std::string("all passed")) << std::endl;
  return failures ? 1 : 0;
}
