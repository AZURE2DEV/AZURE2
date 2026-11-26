// The facade does not depend on the caller's numeric locale (docs/upgrade/qt_gui.md section 19.10): a program that
// sets the user's locale (a Qt application does, from the environment) with a decimal-comma language must read,
// check, convert and save projects as the command line does.  Each call is made under "C" and under a decimal-comma
// LC_NUMERIC of this machine, without any window; the results must be equal byte for byte, and the caller's locale
// must be as it was after each call.
//   front_locale_test <repository root> <scratch dir>        exit 77 (skipped) without a decimal-comma locale
#include <algorithm>
#include <clocale>
#include <vector>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "azure2/front/Engine.h"

namespace fs = std::filesystem;
using namespace azure2;

namespace {

int failures = 0;

void check(const std::string &what, bool ok) {
  std::printf("%s  %s\n", ok ? "  ok  " : "  FAIL", what.c_str());
  if (!ok) failures++;
}

std::string Slurp(const fs::path &p) {
  std::ifstream in(p, std::ios::binary);
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

/// What the facade gives for one project: the validation as text, the saved bytes after one edit, the converted bytes.
struct Results {
  bool opened = false;
  std::string validation, saved, converted, convertReport;
};

Results RunAll(const fs::path &project, const fs::path &oldProject, const fs::path &out) {
  Results r;
  fs::remove_all(out);
  fs::create_directories(out);
  front::Project p;
  std::string why;
  r.opened = front::Project::Open(project.string(), p, why);
  const front::Validation v = p.Validate();
  std::ostringstream text;
  text << v.ok << " " << v.diagnostics.size() << "\n";
  for (const front::Diagnostic &d : v.diagnostics) text << d.severity << " " << d.code << " " << d.text << "\n";
  r.validation = text.str();
  // an edit through the engine's grammar (a record's norm), then the file written
  for (const front::Record &rec : p.Records("data")) {
    front::Record edited = rec;
    for (auto &f : edited.fields)
      if (f.key == "norm") f.value = "1.025 fixed";
    p.SetRecord("data", edited, why);
    break;
  }
  p.Save((out / "saved.azr").string(), why);
  r.saved = Slurp(out / "saved.azr");
  front::ConvertProject(oldProject.string(), (out / "converted.azr").string(), r.convertReport);
  r.converted = Slurp(out / "converted.azr");
  return r;
}

struct Quiet : front::RunObserver {
  void Message(const std::string &, bool) override {}
};

/// The output files of a calculation of a copy of `project` (its folder copied to `dir`), name -> bytes.
std::vector<std::pair<std::string, std::string>> Calculate(const fs::path &project, const fs::path &dir, int threads) {
  fs::remove_all(dir);
  fs::copy(project.parent_path(), dir, fs::copy_options::recursive);
  for (const char *d : {"output", "checks"}) {
    fs::remove_all(dir / d);
    fs::create_directories(dir / d);
  }
  front::RunRequest req;
  req.command = "calculate";
  req.project = (dir / project.filename()).string();
  req.threads = threads;
  Quiet quiet;
  front::Run(req, quiet, nullptr);
  std::vector<std::pair<std::string, std::string>> out;
  for (const auto &e : fs::directory_iterator(dir / "output")) out.push_back({e.path().filename().string(), Slurp(e.path())});
  std::sort(out.begin(), out.end());
  return out;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    std::printf("usage: front_locale_test <repository root> <scratch dir>\n");
    return 2;
  }
  const fs::path root = argv[1], scratch = argv[2];
  const char *comma = nullptr;
  for (const char *name : {"de_DE.utf8", "de_DE.UTF-8", "de_DE", "it_IT.utf8", "it_IT.UTF-8", "pl_PL.utf8", "fr_FR.utf8"})
    if (std::setlocale(LC_NUMERIC, name) && std::string(std::localeconv()->decimal_point) == ",") {
      comma = name;
      break;
    }
  std::setlocale(LC_NUMERIC, "C");
  if (!comma) {
    std::printf("skip: no decimal-comma locale on this machine\n");
    return 77;
  }
  struct Case {
    std::string name;
    fs::path project, old;
  };
  const Case cases[] = {
      {"c12c12_tumino2018", root / "examples/c12c12_tumino2018/c12c12_tumino2018.azr",
       root / "tests/convert/fixtures/history/components_scaling_2026_05/components_scaling_2026_05.azr"},
      {"13N", root / "tests/13N/13N.azr", root / "tests/convert/fixtures/tests/13N/13N.azr"},
      {"T2", root / "tests/upgrade/T2/T2.azr",
       root / "tests/convert/fixtures/history/target_2025_01_coefficients/target_2025_01_coefficients.azr"}};
  for (const Case &c : cases) {
    // the converted project goes beside the old one's data: a copy of the old project's folder in scratch
    const fs::path oldDir = scratch / (c.name + "_old");
    fs::remove_all(oldDir);
    fs::copy(c.old.parent_path(), oldDir, fs::copy_options::recursive);
    const fs::path old = oldDir / c.old.filename();
    if (c.name == "13N") fs::copy(root / "tests/13N/data", oldDir / "data", fs::copy_options::recursive | fs::copy_options::skip_existing);
    std::setlocale(LC_NUMERIC, "C");
    const Results a = RunAll(c.project, old, scratch / (c.name + "_C"));
    std::setlocale(LC_NUMERIC, comma);
    const Results b = RunAll(c.project, old, scratch / (c.name + "_comma"));
    const bool still = std::string(std::localeconv()->decimal_point) == ",";
    std::setlocale(LC_NUMERIC, "C");
    check(c.name + ": opened under C and under " + comma, a.opened && b.opened);
    check(c.name + ": the validation under " + comma + " is that under C", a.validation == b.validation);
    if (a.validation != b.validation) std::printf("        C: %s\n        %s: %s\n", a.validation.substr(0, 300).c_str(), comma, b.validation.substr(0, 300).c_str());
    check(c.name + ": the saved file under " + comma + " is that under C, byte for byte", !a.saved.empty() && a.saved == b.saved);
    check(c.name + ": the converted old project under " + comma + " is that under C, byte for byte",
          !a.converted.empty() && a.converted == b.converted);
    if (a.converted != b.converted) std::printf("        %s\n", b.convertReport.substr(0, 400).c_str());
    check(c.name + ": the caller's locale is " + comma + " again after the calls", still);
  }
  // a calculation (T1, two threads): the run's output files under the comma locale are those under C
  {
    const fs::path t1 = root / "tests/upgrade/T1/T1.azr";
    std::setlocale(LC_NUMERIC, "C");
    const auto a = Calculate(t1, scratch / "run_C", 2);
    std::setlocale(LC_NUMERIC, comma);
    const auto b = Calculate(t1, scratch / "run_comma", 2);
    std::setlocale(LC_NUMERIC, "C");
    bool same = !a.empty() && a.size() == b.size();
    for (std::size_t i = 0; same && i < a.size(); i++) {
      same = a[i] == b[i];
      if (!same) std::printf("        differs: %s\n", a[i].first.c_str());
    }
    check("T1 calculated through the facade under " + std::string(comma) + ": " + std::to_string(a.size()) +
              " output files, byte for byte those under C",
          same);
  }
  std::printf("%s\n", failures ? "FAILED" : "all passed");
  return failures ? 1 : 0;
}
