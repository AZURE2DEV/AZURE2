// The format-1 -> format-2 step of the converter (convert/Convert.h; docs/upgrade/project.md section 4), on every .azr under the given
// directories and on synthetic documents:
//  1. the result is format 2, and converting it again changes nothing (the report is empty);
//  2. apart from the lines it adds to <config> ("useAdaptiveGrid=...", "format=2"), the converted file is the
//     original with lines removed, in order (a subsequence), and every removed line lies in <lastRun>, in <mcmc> outside
//     its <parameters>, in a <potential> without a hybrid potential, or is <potential>'s useAdaptiveGrid line;
//  3. <lastRun> is gone, <mcmc> holds only <parameters>, <potential> remains only with useHybridPotential=1, and a
//     useAdaptiveGrid of <potential> reappears in <config> with its value;
//  4. a format-2 file is untouched (byte for byte).
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "azure2/convert/Convert.h"
#include "azure2/io/Document.h"
#include "azure2/io/ProjectFormat.h"

using namespace azure2::io;
using azure2::convert::ConvertToFormat2;
namespace fs = std::filesystem;

namespace {

int failures = 0;
long checks = 0;

void check(const std::string &what, bool ok) {
  checks++;
  if (!ok) {
    if (failures < 20) std::printf("  FAIL  %s\n", what.c_str());
    failures++;
  }
}

std::string Read(const fs::path &p) {
  std::ifstream in(p, std::ios::binary);
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

std::string Trim(const std::string &s) {
  std::size_t a = s.find_first_not_of(" \t"), b = s.find_last_not_of(" \t");
  return a == std::string::npos ? "" : s.substr(a, b - a + 1);
}

// The line indices of the original that a conversion may remove.
std::vector<bool> Removable(const Document &d, std::string &grid) {
  std::vector<bool> r(d.NumLines(), false);
  auto mark = [&](std::size_t a, std::size_t b) { for (std::size_t i = a; i <= b && i < r.size(); i++) r[i] = true; };
  if (const Document::Section *s = d.Find("lastRun")) mark(s->open, s->close);
  if (const Document::Section *m = d.Find("mcmc")) {
    mark(m->open, m->close);
    if (const Document::Section *p = d.FindChild(*m, "parameters")) {
      bool priors = false;
      for (const std::string &l : d.Body(*p)) if (!Trim(l).empty()) priors = true;
      if (priors) {
        for (std::size_t i = p->open; i <= p->close; i++) r[i] = false;
        r[m->open] = r[m->close] = false;
      }
    }
  }
  if (const Document::Section *p = d.Find("potential")) {
    bool hybrid = false;
    for (const std::string &l : d.Body(*p)) if (Trim(l) == "useHybridPotential=1") hybrid = true;
    for (std::size_t i = p->open + 1; i < p->close; i++)
      if (Trim(d.Line(i)).rfind("useAdaptiveGrid=", 0) == 0) {
        r[i] = true;
        grid = Trim(d.Line(i)).substr(16);
      }
    if (!hybrid) mark(p->open, p->close);
  }
  return r;
}

}  // namespace

int main(int argc, char **argv) {
  std::vector<fs::path> files;
  for (int i = 1; i < argc; i++)
    for (const auto &e : fs::recursive_directory_iterator(argv[i]))
      if (e.is_regular_file() && e.path().extension() == ".azr" &&
          e.path().string().find("/output/") == std::string::npos)
        files.push_back(e.path());
  int converted = 0, removedLastRun = 0, removedPotential = 0, movedGrid = 0;
  for (const fs::path &f : files) {
    const std::string name = f.string(), bytes = Read(f);
    const Document orig = Document::Parse(bytes);
    Document d = Document::Parse(bytes);
    const bool wasOne = FormatOf(orig) == 1;
    if (FormatOf(orig) == 3) {  // a format-3 test project (tests/ratio_to_rutherford/models): left as it is
      check(name + ": a format-3 file is not converted to format 2", ConvertToFormat2(d).empty() && d.Write() == bytes);
      continue;
    }
    const std::vector<std::string> report = ConvertToFormat2(d);
    check(name + ": format 2 after conversion", FormatOf(d) == 2);
    check(name + ": a report for a format-1 file", !wasOne || !report.empty());
    Document again = Document::Parse(d.Write());
    check(name + ": converting twice = once", ConvertToFormat2(again).empty() && again.Write() == d.Write());
    if (!wasOne) {
      check(name + ": a format-2 file is untouched", d.Write() == bytes);
      continue;
    }
    converted++;
    std::string grid;
    const std::vector<bool> removable = Removable(orig, grid);
    // subsequence, apart from the added <config> lines
    std::size_t j = 0;
    bool ok = true;
    for (std::size_t i = 0; i < d.NumLines(); i++) {
      const std::string t = Trim(d.Line(i));
      if (t == "format=2" || t == "useAdaptiveGrid=" + grid) continue;
      while (j < orig.NumLines() && orig.Line(j) != d.Line(i)) {
        if (!removable[j]) ok = false;
        j++;
      }
      if (j == orig.NumLines()) { ok = false; break; }
      j++;
    }
    for (; j < orig.NumLines(); j++) if (!removable[j]) ok = false;
    check(name + ": only the designated lines removed, the rest in order", ok);
    check(name + ": <lastRun> gone", d.Find("lastRun") == nullptr);
    if (orig.Find("lastRun")) removedLastRun++;
    if (const Document::Section *m = d.Find("mcmc")) {
      bool onlyParameters = true;
      for (std::size_t i = m->open + 1; i < m->close; i++) {
        const Document::Section *p = d.FindChild(*m, "parameters");
        if (!p || i < p->open || i > p->close) onlyParameters = false;
      }
      check(name + ": <mcmc> holds only <parameters>", onlyParameters);
    }
    if (orig.Find("potential") && !d.Find("potential")) removedPotential++;
    if (!grid.empty()) {
      movedGrid++;
      bool inConfig = false;
      if (const Document::Section *c = d.Find("config"))
        for (const std::string &l : d.Body(*c)) if (Trim(l) == "useAdaptiveGrid=" + grid) inConfig = true;
      check(name + ": useAdaptiveGrid in <config>", inConfig);
    }
  }
  std::printf("  %zu files: %d converted; <lastRun> removed from %d, <potential> from %d, useAdaptiveGrid moved in %d\n",
              files.size(), converted, removedLastRun, removedPotential, movedGrid);
  // synthetic: a hybrid potential keeps its block, loses its useAdaptiveGrid line
  {
    Document d = Document::Parse("<config>\nx\n</config>\n<potential>\nuseHybridPotential=1\nuseAdaptiveGrid=0\nV0=1\n</potential>\n");
    ConvertToFormat2(d);
    check("hybrid potential kept, grid moved",
          d.Write() == "<config>\nx\nuseAdaptiveGrid=0\nformat=2\n</config>\n<potential>\nuseHybridPotential=1\nV0=1\n</potential>\n");
    Document e = Document::Parse("<config>\n</config>\n<potential>\nuseHybridPotential=1\npair=2\nuseHybridPotential=0\n</potential>\n");
    ConvertToFormat2(e);
    check("a hybrid potential of one pair keeps the block", e.Find("potential") != nullptr);
  }
  // format 3 (docs/upgrade/project_format3.md section 8.4): known as such, and left alone by the format-2 conversion
  {
    const std::string text = "<config>\nformat=3\n</config>\n<lastRun>\nx\n</lastRun>\n";
    Document d = Document::Parse(text);
    check("format=3 is format 3", FormatOf(d) == 3);
    check("a format-3 file is not converted to format 2", ConvertToFormat2(d).empty() && d.Write() == text);
  }
  if (files.empty()) {
    std::printf("  FAIL  no .azr found\n");
    failures++;
  }
  std::printf("%s: %ld checks, %d failure(s)\n", failures ? "FAIL" : "PASS", checks, failures);
  return failures ? 1 : 0;
}
