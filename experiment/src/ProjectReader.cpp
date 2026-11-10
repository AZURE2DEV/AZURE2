// The reading path of a format-3 project (ProjectReader.h): records -> Validate -> pairs of <levels> -> Build.
#include "azure2/experiment/ProjectReader.h"

#include <sstream>

#include "azure2/io/Document.h"
#include "azure2/io/LevelRows.h"
#include "azure2/io/RecordReader.h"

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

}  // namespace

std::vector<io::LevelRow> LevelRowsOf(const io::Document &doc, bool &ok) {
  const std::shared_ptr<const io::LevelSection> levels = LevelSectionOf(doc);
  ok = levels && !levels->malformed;
  return levels ? levels->rows : std::vector<io::LevelRow>();
}

bool RelativisticFlag(const io::Document &doc) {
  const io::Document::Section *s = doc.Find("config");
  if (!s || s->close == io::Document::Section::npos) return false;
  bool relativistic = false;
  for (std::string entry : doc.Body(*s)) {
    const std::size_t hash = entry.find('#');
    if (hash != std::string::npos) entry = entry.substr(0, hash);
    const std::size_t eq = entry.find('=');
    if (eq == std::string::npos) continue;
    std::string key = entry.substr(0, eq), value = entry.substr(eq + 1);
    key.erase(0, key.find_first_not_of(" \t\r"));
    key.erase(key.find_last_not_of(" \t\r") + 1);
    value.erase(0, value.find_first_not_of(" \t\r"));
    value.erase(value.find_last_not_of(" \t\r") + 1);
    if (key == "relativistic") relativistic = value == "1" || value == "true";  // as Config::ReadConfigFile
  }
  return relativistic;
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
