// The <config> section (ConfigSection.h): synthetic sections (the entries, an edit of each kind, the refusals), and
// every project file under the given directories against the stream reading Config::ReadConfigFile did before the
// grammar moved to io/ (copied below as the reference), and written back: each entry set to its own value gives
// the same bytes.
#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "azure2/io/ConfigSection.h"
#include "azure2/io/Document.h"

using namespace azure2::io;
namespace fs = std::filesystem;

namespace {

int failures = 0;
long checks = 0;

void check(const std::string &what, bool ok) {
  checks++;
  if (!ok) {
    if (failures < 30) std::printf("  FAIL  %s\n", what.c_str());
    failures++;
  }
}

std::string Read(const fs::path &p) {
  std::ifstream in(p, std::ios::binary);
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

// The reference: Config::ReadConfigFile's stream reading of the section (before 2026-11-11), the values it kept.
struct Reference {
  std::map<std::string, std::string> values;
  bool closed = false;
};
std::string RefTrim(const std::string &dummy) {
  std::string temp = dummy;
  const std::size_t pound = temp.find_last_of('#');
  if (pound != std::string::npos) temp = temp.substr(0, pound);
  const std::size_t p2 = temp.find_last_not_of(" \n\t\r");
  if (p2 == std::string::npos) return std::string();
  std::size_t p1 = temp.find_first_not_of(" \n\t\r");
  if (p1 == std::string::npos) p1 = 0;
  return temp.substr(p1, p2 - p1 + 1);
}
Reference ReferenceRead(const std::string &text) {
  Reference r;
  std::istringstream in(text);
  std::string line, temp, dummy;
  getline(in, line);
  in >> temp;
  getline(in, dummy);
  r.values["amatrix"] = temp;
  getline(in, dummy);
  r.values["output_dir"] = RefTrim(dummy);
  getline(in, dummy);
  r.values["checks_dir"] = RefTrim(dummy);
  for (const char *name : {"check_compound_nucleus", "check_boundary_conditions", "check_data", "check_energy_dependence",
                           "check_legendre", "check_coulomb_amplitudes", "check_pathways", "check_angular_distributions"}) {
    in >> temp;
    getline(in, dummy);
    r.values[name] = temp;
  }
  line = "";
  while (line != "</config>" && !in.eof()) {
    getline(in, line);
    std::string entry = line;
    const std::size_t hash = entry.find('#');
    if (hash != std::string::npos) entry = entry.substr(0, hash);
    const std::size_t eq = entry.find('=');
    if (eq == std::string::npos) continue;
    std::string key = entry.substr(0, eq), value = entry.substr(eq + 1);
    key.erase(0, key.find_first_not_of(" \t\r"));
    key.erase(key.find_last_not_of(" \t\r") + 1);
    value.erase(0, value.find_first_not_of(" \t\r"));
    value.erase(value.find_last_not_of(" \t\r") + 1);
    r.values[key] = value;
  }
  r.closed = line == "</config>";
  return r;
}

std::vector<std::string> Lines(const std::string &text) {
  std::vector<std::string> out;
  std::istringstream in(text);
  for (std::string l; std::getline(in, l);) out.push_back(l);
  return out;
}

std::vector<std::string> Apply(std::vector<std::string> lines, const std::vector<ConfigLineEdit> &edits) {
  // descending lines, as the binding applies them
  std::vector<ConfigLineEdit> e = edits;
  std::sort(e.begin(), e.end(), [](const ConfigLineEdit &a, const ConfigLineEdit &b) { return a.line > b.line; });
  for (const ConfigLineEdit &x : e) {
    if (x.insert)
      lines.insert(lines.begin() + x.line, x.text);
    else if (x.remove)
      lines.erase(lines.begin() + x.line);
    else
      lines[x.line] = x.text;
  }
  return lines;
}

const std::string kSection =
    "<config>\n"
    "true                #Perform A-Matrix Calculation\n"
    "output/             #Full Path to Output Directory\n"
    "checks/             #Full Path to Checks Directory\n"
    "none                #Compound Nucleus Check\n"
    "screen              #Boundary Condition Check\n"
    "none\n"
    "none\n"
    "none\n"
    "none\n"
    "file   # pathways\n"
    "none\n"
    "includeAllJPi=1\n"
    "priorCovariance = prior.txt   # names\n"
    "format=3\n"
    "</config>\n";

void Synthetic() {
  const std::vector<std::string> lines = Lines(kSection);
  const ConfigSection c = ReadConfigSection(lines);
  check("closed", c.closed);
  check("eleven positional entries and three keys", c.entries.size() == 14);
  check("amatrix", c.Find("amatrix") && c.Find("amatrix")->value == "true" && c.Find("amatrix")->line == 1);
  check("output_dir", c.Find("output_dir")->value == "output/");
  check("check_pathways", c.Find("check_pathways")->value == "file" && c.Find("check_pathways")->line == 10);
  check("priorCovariance", c.Find("priorCovariance")->value == "prior.txt" && !c.Find("priorCovariance")->positional);
  check("unknown entry", c.Find("relativistic") == nullptr);

  std::vector<ConfigLineEdit> e;
  std::string why;
  check("output_dir set", ConfigEdits(lines, "output_dir", "out2/", ConfigChange::Set, e, why) && e.size() == 1 &&
                              e[0].text == "out2/               #Full Path to Output Directory");
  check("a longer output_dir keeps one blank before the comment",
        ConfigEdits(lines, "output_dir", "a/very/long/directory/", ConfigChange::Set, e, why) &&
            e[0].text == "a/very/long/directory/ #Full Path to Output Directory");
  check("a check set", ConfigEdits(lines, "check_data", "screen", ConfigChange::Set, e, why) && e[0].line == 6 &&
                           e[0].text == "screen");
  check("a key replaced in its line", ConfigEdits(lines, "priorCovariance", "p2.txt", ConfigChange::Set, e, why) &&
                                          e[0].text == "priorCovariance = p2.txt      # names");
  check("a key added before </config>", ConfigEdits(lines, "relativistic", "1", ConfigChange::Set, e, why) &&
                                            e.size() == 1 && e[0].insert && e[0].line == 15 && e[0].text == "relativistic=1");
  const std::vector<std::string> added = Apply(lines, e);
  check("the added key reads", ReadConfigSection(added).Find("relativistic")->value == "1");
  check("a key removed", ConfigEdits(lines, "includeAllJPi", "", ConfigChange::Remove, e, why) && e.size() == 1 &&
                             e[0].remove && e[0].line == 12);
  // Refusals, each with the reason.
  check("unknown name", !ConfigEdits(lines, "relativstic", "1", ConfigChange::Set, e, why) && why.find("not an entry") != std::string::npos);
  check("a removed key is no entry", !ConfigEdits(lines, "useAdaptiveGrid", "1", ConfigChange::Set, e, why) &&
                                         why.find("not an entry") != std::string::npos);
  check("bad bool", !ConfigEdits(lines, "relativistic", "yes", ConfigChange::Set, e, why));
  check("bad check word", !ConfigEdits(lines, "check_data", "print", ConfigChange::Set, e, why));
  check("'#' in a path", !ConfigEdits(lines, "output_dir", "out#1/", ConfigChange::Set, e, why));
  check("negative int", !ConfigEdits(lines, "maxOrbitalMomentum", "-1", ConfigChange::Set, e, why));
  check("format not edited", !ConfigEdits(lines, "format", "2", ConfigChange::Set, e, why));
  check("a positional entry is not removed", !ConfigEdits(lines, "output_dir", "", ConfigChange::Remove, e, why));
  // Two lines of one key: the last wins, the earlier one goes.
  std::vector<std::string> twice = lines;
  twice.insert(twice.begin() + 12, "relativistic=0");
  twice.insert(twice.begin() + 13, "relativistic=1");
  check("last line wins", ReadConfigSection(twice).Find("relativistic")->value == "1");
  check("two lines: the earlier removed, the last rewritten",
        ConfigEdits(twice, "relativistic", "0", ConfigChange::Set, e, why) && e.size() == 2 && e[0].remove &&
            e[0].line == 12 && e[1].line == 13 && e[1].text == "relativistic=0");
  // A short section is not closed.
  check("short section not closed", !ReadConfigSection(Lines("<config>\ntrue\nout/\n</config>\n")).closed);
}

}  // namespace

int main(int argc, char **argv) {
  Synthetic();
  long files = 0;
  for (int i = 1; i < argc; i++)
    for (const auto &e : fs::recursive_directory_iterator(argv[i])) {
      if (!e.is_regular_file() || e.path().extension() != ".azr" || e.path().string().find("/output/") != std::string::npos)
        continue;
      const Document d = Document::Parse(Read(e.path()));
      std::size_t open = 0;
      const std::vector<std::string> lines = ConfigLinesOf(d, open);
      if (lines.empty()) continue;
      files++;
      std::string text;
      for (const std::string &l : lines) text += l + "\n";
      const Reference ref = ReferenceRead(text);
      const ConfigSection c = ReadConfigSection(lines);
      const std::string where = e.path().string();
      check(where + ": closed as the reference", c.closed == ref.closed);
      for (const auto &kv : ref.values) {
        const ConfigEntry *x = c.Find(kv.first);
        check(where + ": " + kv.first + " as the reference", x && x->value == kv.second);
      }
      for (const ConfigEntry &x : c.entries) {
        if (x.name == "format") continue;
        std::vector<ConfigLineEdit> edits;
        std::string why;
        // A key of one line, or a positional entry: set to its own value, the same bytes.
        int count = 0;
        for (const ConfigEntry &y : c.entries) count += y.name == x.name;
        const bool known = x.positional || std::any_of(ConfigKeys().begin(), ConfigKeys().end(),
                                                       [&](const ConfigKey &k) { return x.name == k.name; });
        if (!known || count != 1 || (x.positional && x.value.empty())) continue;
        const bool ok = ConfigEdits(lines, x.name, x.value, ConfigChange::Set, edits, why);
        check(where + ": " + x.name + " written with its own value (" + why + ")", ok && Apply(lines, edits) == lines);
      }
    }
  std::printf("  %ld files\n", files);
  if (files == 0) check("no project file found", false);
  std::printf("%s: %ld checks, %d failure(s)\n", failures ? "FAIL" : "PASS", checks, failures);
  return failures ? 1 : 0;
}
