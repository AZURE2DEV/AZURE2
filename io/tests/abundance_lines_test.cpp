// io::ReadAbundanceSection (AbundanceLines.h), the grammar of the <abundances> section: (1) synthetic sections against
// hand values -- the sample prefix, "sum" lines, comments, a prefix alone, the two refusals, the end tag; (2) every
// <abundances> section of the project files under the directories given against the parsing of CNuc::ReadAbundances
// at a16f2a69 (copied below with its actions replaced by records): the same lines, values bit for bit, and the same
// outcome.
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "azure2/io/AbundanceLines.h"
#include "azure2/io/Convert.h"
#include "azure2/io/Document.h"

namespace fs = std::filesystem;
using azure2::io::AbundanceLine;
using azure2::io::AbundanceSection;

namespace {
int failures = 0, checks = 0;
void check(const std::string &what, bool ok) {
  checks++;
  if (!ok) {
    std::printf("  FAIL  %s\n", what.c_str());
    failures++;
  }
}

AbundanceSection Read(const std::string &text) {
  std::istringstream in(text);
  return azure2::io::ReadAbundanceSection(in);
}

// ---- the parsing of CNuc::ReadAbundances at a16f2a69: its loop, with sampleOf(...) and the nucleus search replaced by
//      appending a record, and each "return -1" by the problem it reports ----
AbundanceSection OldRead(std::istream &in) {
  AbundanceSection out;
  std::string line = "";
  getline(in, line);  // <abundances>
  line = "";
  while (!in.eof() && line != "</abundances>") {
    getline(in, line);
    std::string entry = line;
    size_t hash = entry.find('#');
    if (hash != std::string::npos) entry = entry.substr(0, hash);
    if (entry.find_first_not_of(" \t\r") == std::string::npos) continue;
    if (line == "</abundances>") break;
    std::string name;
    std::istringstream probe(entry);
    std::string first;
    probe >> first;
    if (!first.empty() && first.back() == ':') {
      name = first.substr(0, first.size() - 1);
      std::string rest;
      getline(probe, rest);
      entry = rest;
    }
    std::istringstream stm(entry);
    std::string head;
    if (!(stm >> head)) continue;
    if (head == "sum") {
      double s;
      if (!(stm >> s)) {
        out.problem = AbundanceSection::Problem::BadSum;
        out.badLine = line;
        return out;
      }
      AbundanceLine a;
      a.sample = name;
      a.isSum = true;
      a.sumSigma = s;
      a.text = line;
      out.lines.push_back(a);
      continue;
    }
    std::istringstream stm2(entry);
    int z, a, vary;
    double value, sigma;
    if (!(stm2 >> z >> a >> value >> vary >> sigma)) {
      out.problem = AbundanceSection::Problem::BadEntry;
      out.badLine = line;
      return out;
    }
    AbundanceLine r;
    r.sample = name;
    r.Z = z;
    r.A = a;
    r.value = value;
    r.vary = vary;
    r.sigma = sigma;
    r.text = line;
    out.lines.push_back(r);
  }
  out.complete = (line == "</abundances>");
  return out;
}

bool Same(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }
bool SameLine(const AbundanceLine &a, const AbundanceLine &b) {
  return a.sample == b.sample && a.isSum == b.isSum && Same(a.sumSigma, b.sumSigma) && a.Z == b.Z && a.A == b.A &&
         Same(a.value, b.value) && a.vary == b.vary && Same(a.sigma, b.sigma) && a.text == b.text;
}
bool SameSection(const AbundanceSection &a, const AbundanceSection &b) {
  if (a.lines.size() != b.lines.size() || a.problem != b.problem || a.badLine != b.badLine || a.complete != b.complete)
    return false;
  for (std::size_t k = 0; k < a.lines.size(); k++)
    if (!SameLine(a.lines[k], b.lines[k])) return false;
  return true;
}

/// The <abundances> section as Config::SectionStream gives it (trimmed tags), or "" without one.
std::string SectionText(const std::string &bytes) {
  azure2::io::Document doc = azure2::io::Document::Parse(bytes);
  azure2::io::ConvertToFormat2(doc);
  const azure2::io::Document::Section *s = doc.Find("abundances");
  if (!s) return "";
  std::string text;
  const std::size_t end = (s->close == azure2::io::Document::Section::npos) ? doc.NumLines() : s->close + 1;
  for (std::size_t i = s->open; i < end; i++) {
    bool closing = false;
    const std::string tag = azure2::io::TagOf(doc.Line(i), closing);
    text += (i == s->open || (i == s->close && !tag.empty())) ? (closing ? "</" : "<") + tag + ">" : doc.Line(i);
    text += "\n";
  }
  return text;
}
}  // namespace

int main(int argc, char **argv) {
  // (1) Synthetic sections.
  const std::string text =
      "<abundances>\n"
      "# natural neon\n"
      "10 20 0.9048 1 0.001   # 20Ne\n"
      "gas: 10 22 0.0925 0 0\n"
      "\n"
      "gas:\n"
      "gas: sum 0.01\n"
      "sum 0.002 extra\n"
      "</abundances>\n";
  const AbundanceSection s = Read(text);
  check("four lines, complete, no problem", s.lines.size() == 4 && s.complete && s.problem == AbundanceSection::Problem::None);
  if (s.lines.size() == 4) {
    check("an unnamed entry", s.lines[0].sample.empty() && !s.lines[0].isSum && s.lines[0].Z == 10 && s.lines[0].A == 20 &&
                                  s.lines[0].value == 0.9048 && s.lines[0].vary == 1 && s.lines[0].sigma == 0.001 &&
                                  s.lines[0].text == "10 20 0.9048 1 0.001   # 20Ne");
    check("a named entry", s.lines[1].sample == "gas" && s.lines[1].A == 22 && s.lines[1].vary == 0 && s.lines[1].sigma == 0);
    check("a named sum line (a prefix alone is skipped)", s.lines[2].sample == "gas" && s.lines[2].isSum &&
                                                               s.lines[2].sumSigma == 0.01);
    check("an unnamed sum line, extra tokens ignored", s.lines[3].sample.empty() && s.lines[3].isSum && s.lines[3].sumSigma == 0.002);
  }
  const AbundanceSection badSum = Read("<abundances>\n10 20 1 0 0\nsum x\n10 22 1 0 0\n</abundances>\n");
  check("a sum without a number: BadSum, the line before kept", badSum.problem == AbundanceSection::Problem::BadSum &&
                                                                    badSum.badLine == "sum x" && badSum.lines.size() == 1 &&
                                                                    !badSum.complete);
  const AbundanceSection badEntry = Read("<abundances>\n10 20 1 0\n</abundances>\n");
  check("four numbers: BadEntry", badEntry.problem == AbundanceSection::Problem::BadEntry && badEntry.badLine == "10 20 1 0");
  const AbundanceSection open = Read("<abundances>\n10 20 1 0 0\n");
  check("no </abundances>: incomplete", !open.complete && open.lines.size() == 1);
  check("the synthetic sections read as before", SameSection(s, [&] { std::istringstream in(text); return OldRead(in); }()));

  // (2) The project files: the new reader against the old parsing.
  long files = 0;
  for (int i = 1; i < argc; i++)
    for (const auto &e : fs::recursive_directory_iterator(argv[i])) {
      if (!e.is_regular_file() || e.path().extension() != ".azr" || e.path().string().find("/output/") != std::string::npos)
        continue;
      std::ifstream f(e.path(), std::ios::binary);
      std::ostringstream bytes;
      bytes << f.rdbuf();
      const std::string section = SectionText(bytes.str());
      if (section.empty()) continue;
      std::istringstream oldIn(section);
      check(e.path().string() + ": as the parsing before the move", SameSection(OldRead(oldIn), Read(section)));
      files++;
    }
  if (argc > 1) check("project files with <abundances> found", files > 10);
  std::printf("%s: %d checks (%ld project files), %d failure(s)\n", failures ? "FAIL" : "PASS", checks, files, failures);
  return failures ? 1 : 0;
}
