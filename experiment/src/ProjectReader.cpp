// The reading path of a format-3 project (ProjectReader.h): records -> Validate -> pairs of <levels> -> Build.
#include "azure2/experiment/ProjectReader.h"

#include <map>
#include <set>
#include <sstream>

#include "azure2/io/ConfigSection.h"
#include "azure2/io/Document.h"
#include "azure2/io/LevelRows.h"
#include "azure2/io/RecordReader.h"
#include "azure2/experiment/ThmSectionLines.h"

namespace azure2 {
namespace experiment {

bool Format3Project::Ok() const {
  for (const Diagnostic &d : diagnostics)
    if (d.severity == Severity::Error) return false;
  return true;
}

namespace {

// The <levels> section of a document read with io::ReadLevelSection; null when it is absent or not closed.
std::shared_ptr<const io::LevelSection> LevelSectionOf(const io::Document &doc) {
  const io::Document::Section *s = doc.Find("levels");
  if (!s || s->close == io::Document::Section::npos) return nullptr;
  std::ostringstream text;
  text << "<levels>\n";
  for (const std::string &line : doc.Body(*s)) text << line << "\n";
  text << "</levels>\n";
  std::istringstream in(text.str());
  return std::make_shared<const io::LevelSection>(io::ReadLevelSection(in));
}

// The names of the experiment[<name>] lines of <thm> (ReadThmSectionLine; the rest of the line is the THM block's
// grammar, read by the engine: src/ConfigThm.cpp).
std::set<std::string> ThmExperimentNames(const io::Document &doc) {
  std::set<std::string> names;
  const io::Document::Section *s = doc.Find("thm");
  if (!s) return names;
  for (const std::string &line : doc.Body(*s)) {
    const ThmSectionLine l = ReadThmSectionLine(line);
    if (l.experiment && l.text.find(']') != std::string::npos) names.insert(l.name);
  }
  return names;
}

// The old member and weight lists (segments= on an experiment line, weight[<k>]= / weightTest[<k>]=) do not exist in
// format 3 (project_format3.md section 3.9): the records say which data sets an experiment has and their weights.
// The engine would otherwise read such a list beside the members the records give.
void CheckThmOldLists(const io::Document &doc, std::vector<Diagnostic> &out) {
  const io::Document::Section *s = doc.Find("thm");
  if (!s) return;
  const std::vector<std::string> body = doc.Body(*s);
  for (std::size_t i = 0; i < body.size(); i++) {
    const ThmSectionLine l = ReadThmSectionLine(body[i]);
    std::string field;
    const std::string key = ThmOldListKey(l, field);
    if (key.empty()) continue;
    Diagnostic d;
    d.code = codes::kThmOldList;
    d.line = (int)(s->open + 2 + i);
    d.section = "thm";
    d.hasFound = true;
    d.found = key + "=";
    d.message = l.experiment ? "<thm> experiment[" + l.name + "]: segments= is the old member list; in format 3 a data set "
                                   "joins the experiment with its field '" + field + "'"
                             : "<thm> " + key + "= is the old weight list; in format 3 the data set or extrapolation names "
                                   "its weight table with its field '" + field + "'";
    out.push_back(d);
  }
}

// Every `thm <name>` of an active data set or extrapolation names an experiment of <thm> (`global` names none):
// otherwise the engine's THM set-up would refuse the experiment it does have for having no data sets, not the name that
// is wrong.  An inactive record is checked when it is switched on (allRecordsActive), as its experiment's membership is
// (Build adds active members only).
void CheckThmExperimentNames(const io::Document &doc, const ValidationResult &validated, bool allRecordsActive,
                             std::vector<Diagnostic> &out) {
  const std::set<std::string> names = ThmExperimentNames(doc);
  for (const ValidatedRecord &r : validated.records) {
    const Value *active = r.Get("active");
    if (!allRecordsActive && active && active->word == "no") continue;
    for (const ValidatedField &f : r.fields) {
      if (f.key != "thm" || f.value.word == "global" || names.count(f.value.word)) continue;
      Diagnostic d;
      d.code = codes::kThmExperimentUnknown;
      d.line = f.line;
      d.column = f.col;
      d.section = r.section;
      d.recordKind = r.kind;
      d.record = r.name;
      d.field = "thm";
      d.hasFound = true;
      d.found = f.value.word;
      d.allowed.assign(names.begin(), names.end());
      d.allowed.push_back("global");
      d.message = "<thm> has no experiment[" + f.value.word + "] line" +
                  (names.empty() ? std::string(" (it has no experiment)") : "");
      out.push_back(d);
    }
  }
}

// "use_brune=1 (this file)" or "(the default)": an option's value and where it came from, for a message.
std::string OptionWords(const char *name, bool on, bool written) {
  return std::string(name) + "=" + (on ? "1" : "0") + (written ? " (this file)" : " (the default)");
}

// The options of the calculation written in <config> (project_format3.md section 3.12): every value one the engine
// reads (E_CONFIG_OPTION_VALUE), and no two that cannot hold together, an option the file does not state taken at its
// default (E_CONFIG_OPTIONS_CONTRADICT).  A file without option keys gets no diagnostic here.
void CheckConfigOptions(const io::Document &doc, std::vector<Diagnostic> &out) {
  std::size_t open = 0;
  const io::ConfigSection config = io::ReadConfigSection(io::ConfigLinesOf(doc, open));
  struct Written {
    bool on;
    int line;
  };
  std::map<std::string, Written> written;  // the line that wins of each option the file states
  for (const io::ConfigEntry &e : config.entries) {
    if (e.positional || !io::ConfigOptionNamed(e.name)) continue;
    const int line = (int)(open + e.line) + 1;
    bool on = false;
    if (io::ConfigOptionValue(e.value, on)) {
      written[e.name] = {on, line};
      continue;
    }
    Diagnostic d;
    d.code = codes::kConfigOptionValue;
    d.line = line;
    d.section = "config";
    d.field = e.name;
    d.hasFound = true;
    d.found = e.value;
    d.allowed = {"0", "1", "true", "false"};
    d.message = e.name + "=" + e.value + ": an option of the calculation is 0, 1, true or false";
    out.push_back(d);
  }
  auto value = [&written](const char *name) {
    const auto w = written.find(name);
    return w != written.end() ? w->second.on : io::ConfigOptionNamed(name)->defaultOn;
  };
  for (const io::ConfigOptionContradiction &c : io::ConfigOptionContradictions()) {
    const auto a = written.find(c.a), b = written.find(c.b);
    const bool aWritten = a != written.end(), bWritten = b != written.end();
    if ((!aWritten && !bWritten) || value(c.a) != c.aValue || value(c.b) != c.bValue) continue;
    Diagnostic d;
    d.code = codes::kConfigOptionsContradict;
    d.line = aWritten ? a->second.line : b->second.line;
    d.section = "config";
    d.field = aWritten ? c.a : c.b;
    d.message = OptionWords(c.a, c.aValue, aWritten) + " with " + OptionWords(c.b, c.bValue, bWritten) +
                " cannot hold together: " + c.why;
    out.push_back(d);
  }
}

}  // namespace

std::vector<io::LevelRow> LevelRowsOf(const io::Document &doc, bool &ok) {
  const std::shared_ptr<const io::LevelSection> levels = LevelSectionOf(doc);
  ok = levels && !levels->malformed;
  return levels ? levels->rows : std::vector<io::LevelRow>();
}

bool RelativisticFlag(const io::Document &doc) {
  std::size_t open = 0;
  const io::ConfigSection config = io::ReadConfigSection(io::ConfigLinesOf(doc, open));
  const io::ConfigEntry *e = config.Find("relativistic");
  return config.closed && e && !e->positional && (e->value == "1" || e->value == "true");  // as Config::ReadConfigFile
}

Format3Project ReadFormat3Project(const io::Document &doc, const Format3ReadOptions &options) {
  Format3Project p;
  std::vector<SectionRecords> sections;
  for (const io::Document::Section &s : doc.Sections()) {
    if (s.depth != 0 || Schema::KindsOfSection(s.name).empty()) continue;
    const io::RecordReadResult r = io::ReadRecords(doc, s, Schema::KindsOfSection(s.name));
    for (const io::RecordDiagnostic &d : r.diagnostics) p.diagnostics.push_back(FromRecordDiagnostic(d, s.name));
    sections.push_back({s.name, r.set.records});
  }
  p.validated = Schema::Validate(sections);
  p.diagnostics.insert(p.diagnostics.end(), p.validated.diagnostics.begin(), p.validated.diagnostics.end());
  CheckThmExperimentNames(doc, p.validated, options.allRecordsActive, p.diagnostics);
  CheckThmOldLists(doc, p.diagnostics);
  // The old <potential> block is not read in format 3 (its default section would give a potential to pairs the file
  // does not name): the converter writes one potential record per pair.
  if (const io::Document::Section *old = doc.Find("potential")) {
    Diagnostic d;
    d.code = codes::kOldPotentialBlock;
    d.line = (int)old->open + 1;
    d.section = "potential";
    d.message = "<potential> is the old block of the hybrid potentials; format 3 states them as potential records in "
                "<potentials>: convert the original with AZURE2 --convert OLD NEW --to 3";
    p.diagnostics.push_back(d);
  }
  // useAdaptiveGrid chose a uniform grid of its own for a folded point (0) or the adaptive one (1): the uniform grid
  // is gone (it lost narrow resonances), so the key would say nothing; a file that still has it is refused, and the
  // line is deleted (the converter drops it from old files).
  {
    std::size_t open = 0;
    const io::ConfigSection config = io::ReadConfigSection(io::ConfigLinesOf(doc, open));
    for (const io::ConfigEntry &e : config.entries) {
      if (e.positional || e.name != "useAdaptiveGrid") continue;
      Diagnostic d;
      d.code = codes::kRemovedConfigKey;
      d.line = (int)(open + e.line) + 1;
      d.section = "config";
      d.field = e.name;
      d.message = "useAdaptiveGrid is no longer a key of <config>: delete the line.  Every folded point is on the "
                  "shared mesh of its nucleus, which is never coarser than the uniform grid of the same number of "
                  "points (the uniform grid, useAdaptiveGrid=0, lost narrow resonances and is gone); a value of 1 "
                  "said what is now the only grid.";
      p.diagnostics.push_back(d);
    }
  }
  CheckConfigOptions(doc, p.diagnostics);
  int buildChecks = 0;
  for (const ValidatedRecord &r : p.validated.records)
    if (r.quantity)
      for (const Refusal &x : r.quantity->refusals)
        if (x.checkedBy == "build") buildChecks++;
  p.levelSection = LevelSectionOf(doc);
  p.levelsRead = p.levelSection && !p.levelSection->malformed;
  if (p.levelSection) p.levels = LevelTablesOfRows(p.levelSection->rows);
  p.projectDir = options.projectDir;
  if (options.build && p.Ok() && p.levelsRead) {
    BuildOptions b;
    b.projectDir = options.projectDir;
    b.readDataFiles = options.readDataFiles;
    b.dataFiles = options.dataFiles;
    b.setUpChecks = options.setUpChecks;
    b.allRecordsActive = options.allRecordsActive;
    b.relativistic = RelativisticFlag(doc);
    p.built = Schema::Build(p.validated, p.levels, b);
    p.diagnostics.insert(p.diagnostics.end(), p.built.diagnostics.begin(), p.built.diagnostics.end());
    p.buildChecksSkipped = p.built.checksLeftToEngine;
  } else {
    p.buildChecksSkipped = buildChecks;
  }
  for (Diagnostic &d : p.diagnostics) d.file = options.fileName;
  return p;
}

Format3Project ReadFormat3Text(const std::string &text, const Format3ReadOptions &options) {
  return ReadFormat3Project(io::Document::Parse(text), options);
}

}  // namespace experiment
}  // namespace azure2
