// The front ends' facade without Qt (docs/upgrade/qt_gui.md section 5): front_test <source dir> <AZURE2 binary>.
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
    const std::string command = "\"" + binary + "\" calculate --threads=2 \"" + a.string() + "\" > \"" +
                                (a.parent_path() / "stdout.txt").string() + "\" 2>&1";
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

  std::cout << (failures ? "FAILED: " + std::to_string(failures) : std::string("all passed")) << std::endl;
  return failures ? 1 : 0;
}
