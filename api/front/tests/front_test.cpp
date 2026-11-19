// The front ends' facade without Qt (docs/upgrade/qt_gui.md section 10): front_test <source dir> <AZURE2 binary>.
//  (a) the tables: schema, commands, options, fitters, <config> keys;
//  (b) the test projects T1, T2, T3 and 13N opened and written back byte for byte, also after every record, every
//      <levels> row and every <config> entry is set to what it holds; an edit changes only its lines and is undone;
//      a record added, renamed and removed; --validate's diagnostics, a broken field named by its code;
//  (c) an old project refused on opening and converted by the converter of --convert;
//  (d) a run of `calculate` through the facade writes the same output files, byte for byte, as the program on the same
//      project, and returns the chi^2 the program prints; T1 and 13N;
//  (e) a fit stopped through the facade ends at its best point with exit 0 and fitStatus.out saying "stopped";
//  (f) physics: the rows of the physics table and one evaluated; a session of T1: poles, channel functions against
//      the row, the 1/2+ block's phase shift against the program's extrapolation, |U| = 1, |det A^-1| least at the
//      level; T2's THM vertex.
// Each check prints PASS or FAIL; exit 1 on a failure.
#include "azure2/front/Engine.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
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

std::string Sci(double x) {
  char b[32];
  std::snprintf(b, sizeof(b), "%.2e", x);
  return b;
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

  // (f) physics (docs/upgrade/physics_tab.md section 7)
  {
    // The rows: the physics table pyazr.physics reads (26 on 2026-11-19), every row with its axis among its arguments.
    const std::vector<front::PhysicsFunction> rows = front::PhysicsFunctions();
    bool axes = !rows.empty(), penetrability = false;
    for (const front::PhysicsFunction &f : rows) {
      bool found = false;
      for (const front::PhysicsArgument &a : f.arguments) found = found || (a.name == f.axis && a.hasSweep);
      axes = axes && found;
      penetrability = penetrability || f.name == "penetrability";
    }
    Check(rows.size() >= 26 && axes && penetrability,
          std::to_string(rows.size()) + " physics rows, each with its axis an argument with a sweep");
    std::string why;
    front::PhysicsValues v;
    // A broadcast column, a value outside its domain (e_cm = -1: argument 7) and its words.
    const bool evaluated = front::EvaluatePhysics(
        "penetrability", {{0}, {1}, {6}, {1.00782503}, {12.0}, {5.0}, {0.5, 1.0, -1.0}, {0}}, v, why);
    Check(evaluated && v.status.size() == 3 && v.status[0] == 0 && v.status[2] == 7 && std::isnan(v.values[0][2]) &&
              v.values[0][1] > v.values[0][0],
          "penetrability on three energies: computed, rising, and e_cm = -1 refused as argument 7 (" + why + ")");
    Check(front::PhysicsStatusText("penetrability", 7).find("e_cm") != std::string::npos,
          "the status words name e_cm: " + front::PhysicsStatusText("penetrability", 7));
    why.clear();
    const bool unknownRow = !front::EvaluatePhysics("no_such_row", {}, v, why);
    Check(unknownRow && !why.empty(), "an unknown row is refused: " + why);

    // A session of T1 (Brune's basis): pairs, groups, poles.
    const fs::path t1 = CopyProject(source / "tests/upgrade/T1/T1.azr", fs::path("front_test_work") / "physics");
    front::PhysicsSession s;
    why.clear();
    const bool opened = front::PhysicsSession::Open(t1.string(), "", false, s, why);
    Check(opened && s.IsOpen(), "T1: a physics session (" + why + ")");
    Check(s.Parametrization() == "brune", "T1: parametrization " + s.Parametrization());
    const std::vector<front::PhysicsPair> pairs = s.Pairs();
    Check(pairs.size() >= 2 && pairs[0].kind == "particle" && pairs[0].z1 * pairs[0].z2 == 6,
          "T1: pair 1 is p + 12C");
    int group = 0;
    for (const front::PhysicsGroup &g : s.Groups())
      if (g.J == 0.5 && g.parity == 1) group = g.jgroup;
    Check(group > 0, "T1: the J group 1/2+ is " + std::to_string(group));
    // The poles: the 1/2+ level at 2.365 MeV with Gamma_p = 33000 eV (T1.azr's <levels>), back from the engine's
    // transformation (Brune's basis: the numbers given are the observed ones).
    bool pole = false;
    for (const front::PhysicsPole &p : s.Poles())
      if (p.jgroup == group && std::fabs(p.energy - 2.365) < 1e-9 && !p.widths.empty())
        pole = std::fabs(p.widths[0].value - 33000.0) < 1e-6 && p.widths[0].kind == "width" &&
               std::fabs(p.widths[0].width - 0.033) < 1e-12 && p.totalWidth >= 0.033;
    Check(pole, "T1: the pole 1/2+ at E_x = 2.365 MeV with Gamma_p = 33000 eV");

    // The session's channel functions against the session-free row of the same pair (the row computes F and G
    // exactly; a session with a Coulomb table interpolates them to about 1e-7).
    front::PhysicsCurves cf;
    const std::vector<double> energies = {0.3, 0.6, 1.2};
    why.clear();
    const bool functions = s.ChannelFunctions(1, 0, 0.0, -1, energies, cf, why);
    Check(functions && cf.Find("P") && cf.Find("open"), "T1: channel functions of pair 1, l = 0 (" + why + ")");
    double radius = 0;
    for (const front::PhysicsGroup &g : s.Groups())
      for (const front::PhysicsChannel &c : g.channels)
        if (c.pair == 1 && c.l == 0) radius = c.radius;
    front::EvaluatePhysics("penetrability",
                           {{0}, {(double)pairs[0].z1}, {(double)pairs[0].z2}, {pairs[0].m1}, {pairs[0].m2}, {radius},
                            energies, {0}},
                           v, why);
    double worst = 0;
    for (std::size_t i = 0; i < energies.size() && cf.Find("P"); i++)
      worst = std::max(worst, std::fabs(cf.Find("P")->values[i] / v.values[0][i] - 1.0));
    Check(cf.Find("P") && worst < 1e-6, "T1: session P = the row's P to " + Sci(worst));

    // The 1/2+ block on the E_x of T1's extrapolation x3 (12C(p,p) phase shift, J = 1/2, l = 0): the facade's delta
    // against the program's output of `extrapolate` (AZUREOut, read by ReadCurves), and |U_pp| = 1 (one open particle
    // channel: the collision matrix is unitary there).
    Collect observer;
    front::RunRequest request;
    request.command = "extrapolate";
    request.project = t1.string();
    request.threads = 2;
    const front::RunResult run = front::Run(request, observer, nullptr);
    front::Project project;
    front::Project::Open(t1.string(), project, why);
    const front::Curves curves = front::ReadCurves(project);
    const front::Curve *x3 = nullptr;
    for (const front::Curve &c : curves.curves)
      if (c.record == "x3") x3 = &c;
    Check(run.exitCode == 0 && x3 && x3->points.size() > 10, "T1: extrapolate wrote the phase shifts of x3");
    if (x3) {
      // E_x of each point from the record's lab energies (the file prints 7 digits; dDelta/dE reaches 5000 deg/MeV at
      // the level): the row lab_to_cm_energy, the engine's conversion of a lab point, plus the pair's threshold.
      std::vector<double> lab;
      for (const front::Record &r : project.Records("extrapolations"))
        if (r.name == "x3") {
          std::istringstream words(r.Value("energies"));
          std::string w;
          while (words >> w)
            if (w != "list" && w.find_first_not_of("0123456789.eE+-") == std::string::npos) lab.push_back(std::stod(w));
        }
      front::EvaluatePhysics("lab_to_cm_energy", {{pairs[0].m1}, {pairs[0].m2}, lab, {0}}, v, why);
      std::vector<double> ex;
      for (std::size_t i = 0; i < lab.size() && v.values.size() == 1; i++)
        ex.push_back(v.values[0][i] + pairs[0].separation + pairs[0].excitation);
      Check(ex.size() == x3->points.size(), "T1: " + std::to_string(ex.size()) + " energies of x3");
      front::PhysicsCurves block;
      why.clear();
      const bool blocked = s.Block(group, ex, block, why);
      Check(blocked, "T1: the 1/2+ block on " + std::to_string(ex.size()) + " energies (" + why + ")");
      const front::PhysicsSeries *delta = block.Find("delta", {{"channel", 1}, {"channel_out", 1}});
      const front::PhysicsSeries *eta = block.Find("eta", {{"channel", 1}, {"channel_out", 1}});
      double dmax = 0, umax = 0;
      for (std::size_t i = 0; delta && eta && i < ex.size() && i < x3->points.size(); i++) {
        // the file's 7 significant digits
        dmax = std::max(dmax, std::fabs(delta->values[i] - x3->points[i].model) / std::max(1.0, std::fabs(x3->points[i].model)));
        umax = std::max(umax, std::fabs(eta->values[i] - 1.0));
      }
      Check(delta && dmax < 1e-6, "T1: delta of the block = the program's phase shift of x3 to " +
                                      Sci(dmax) + " relative");
      Check(eta && umax < 1e-12, "T1: |U_pp| = 1 to " + Sci(umax));
      // log10 |det A^-1| on a grid across the 2.365 MeV level (Gamma = 33 keV): its minimum lies within Gamma.
      std::vector<double> grid;
      for (int i = 0; i <= 400; i++) grid.push_back(2.2 + 0.3 * i / 400.0);
      s.Block(group, grid, block, why);
      const front::PhysicsSeries *det = block.Find("log10_abs_det");
      std::size_t at = 0;
      for (std::size_t i = 0; det && i < grid.size(); i++)
        if (det->values[i] < det->values[at]) at = i;
      Check(det && std::fabs(grid[at] - 2.365) < 0.033,
            "T1: |det A^-1| is least at E_x = " + std::to_string(grid[at]) + " MeV, within Gamma of 2.365");
    }
    why.clear();
    const bool noGroup = !s.Block(99, {2.0}, cf, why);
    Check(noGroup && !why.empty(), "T1: a J group that is not there is refused: " + why);

    // THM: T2's experiment LC09, its vertex at three energies (pyazr's Session.thm_vertex reads the same table).
    const fs::path t2 = CopyProject(source / "tests/upgrade/T2/T2.azr", fs::path("front_test_work") / "physics_t2");
    front::PhysicsSession s2;
    why.clear();
    const bool opened2 = front::PhysicsSession::Open(t2.string(), "", false, s2, why);
    Check(opened2, "T2: a physics session (" + why + ")");
    const std::vector<std::string> experiments = s2.ThmExperiments();
    Check(experiments.size() == 1 && experiments[0] == "LC09", "T2: the THM experiment LC09");
    front::PhysicsCurves vertex;
    why.clear();
    bool positive = s2.ThmVertex("LC09", {0.3, 0.5, 0.8}, vertex, why) && vertex.Find("M2");
    for (const front::PhysicsSeries &series : vertex.series)
      if (series.name == "M2")
        for (double x : series.values) positive = positive && x > 0 && std::isfinite(x);
    Check(positive, "T2: <|M_l|^2> of LC09 finite and positive (" + why + ")");
    why.clear();
    const bool noExperiment = !s2.ThmVertex("none", {0.5}, vertex, why);
    Check(noExperiment && !why.empty(), "T2: an unknown experiment refused: " + why);
  }

  std::cout << (failures ? "FAILED: " + std::to_string(failures) : std::string("all passed")) << std::endl;
  return failures ? 1 : 0;
}
