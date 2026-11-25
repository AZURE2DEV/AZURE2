// The converter old -> format 3 (Convert.h): format 2 first (ConvertToFormat2, Format2.cpp), the old sections read
// (OldProject.cpp), the records built (Build.cpp), the document assembled here and validated with the table.
#include "azure2/convert/Convert.h"

#include <algorithm>
#include <sstream>

#include "Build.h"
#include "OldEra.h"
#include "OldProject.h"
#include "azure2/experiment/Build.h"
#include "azure2/experiment/Json.h"
#include "azure2/experiment/ProjectReader.h"
#include "azure2/experiment/Schema.h"
#include "azure2/io/ProjectFormat.h"
#include "azure2/io/Document.h"
#include "azure2/io/RecordReader.h"

namespace azure2 {
namespace convert {

namespace {

/// The old sections that format 3 replaces (section 8.4), in the order a refusal names them.
const std::vector<std::string> kOldSections = {"segmentsData", "segmentsTest", "targetInt", "abundances", "scatteringLengths"};

/// The sections of an old project the converter reads or removes (formats 1 and 2), besides the format-3 sections a
/// format-2 file may already hold (Schema::Sections, Schema::GrammarSections), and the sections each may hold.
const std::vector<std::pair<std::string, std::vector<std::string>>> kConvertedSections = {
    {"segmentsData", {}}, {"segmentsTest", {}}, {"targetInt", {}}, {"abundances", {}}, {"scatteringLengths", {}},
    {"potential", {}},    {"lastRun", {}},      {"mcmc", {"parameters"}}};

/// The refusal of a section the converter does not convert ("" when there is none): a section of a newer AZURE2 (its
/// dev branch), or one no AZURE2 the converter knows wrote.  Copying it through would give a project that runs as if
/// it were absent.
std::string UnconvertedSection(const io::Document &doc) {
  auto childrenOf = [](const std::string &name, bool &known) {
    known = true;
    for (const auto &c : kConvertedSections)
      if (c.first == name) return c.second;
    for (const experiment::GrammarSection &g : experiment::Schema::GrammarSections())
      if (g.name == name) return g.children;
    known = !experiment::Schema::KindsOfSection(name).empty();
    return std::vector<std::string>();
  };
  const std::vector<io::Document::Section> &all = doc.Sections();
  for (const io::Document::Section &s : all) {
    bool known = false;
    if (s.depth > 0) {
      const io::Document::Section &p = all[s.parent];
      if (p.depth > 0) continue;
      const std::vector<std::string> children = childrenOf(p.name, known);
      if (!known || std::find(children.begin(), children.end(), s.name) != children.end()) continue;
      return "<" + s.name + "> inside <" + p.name + "> (line " + std::to_string(s.open + 1) + ") is not a section the "
             "converter reads there; written through, nothing would read it and the converted project would compute "
             "as if it were absent.  Remove it from a copy of the project, or move its lines where they belong, and "
             "convert the copy";
    }
    childrenOf(s.name, known);
    if (known) continue;
    for (const auto &n : experiment::Schema::NewerSections())
      if (n.first == s.name)
        return "<" + s.name + "> (line " + std::to_string(s.open + 1) + ") is a section of a newer AZURE2 (its dev "
               "branch: " + n.second + ") that this build does not compute yet; converted without it, the project "
               "would give other numbers than the program that wrote it.  Keep running it with that program, or "
               "convert a copy without the section if those numbers may change";
    return "<" + s.name + "> (line " + std::to_string(s.open + 1) + ") is not a section of any AZURE2 project the "
           "converter knows; written through, nothing would read it and the converted project would compute as if it "
           "were absent.  Remove it from a copy of the project (or correct its name) and convert the copy";
  }
  return "";
}

/// Replaces section `name` (tags included) by `newName` holding `body`; nothing when it is absent.
bool ReplaceSection(io::Document &doc, const std::string &name, const std::string &newName, const std::vector<std::string> &body) {
  const io::Document::Section *s = doc.Find(name);
  if (!s || s->close == io::Document::Section::npos) return false;
  std::vector<std::string> lines;
  if (!newName.empty()) {
    lines.push_back("<" + newName + ">");
    lines.insert(lines.end(), body.begin(), body.end());
    lines.push_back("</" + newName + ">");
  }
  doc.ReplaceLines(s->open, s->close + 1, lines);
  return true;
}

/// Inserts a new section before the section `before` (or at the end of the document).
void InsertSection(io::Document &doc, const std::string &before, const std::string &name, const std::vector<std::string> &body) {
  std::vector<std::string> lines = {"<" + name + ">"};
  lines.insert(lines.end(), body.begin(), body.end());
  lines.push_back("</" + name + ">");
  const io::Document::Section *s = doc.Find(before);
  const std::size_t at = s ? s->open : doc.NumLines();
  doc.ReplaceLines(at, at, lines);
}

void ReplaceBody(io::Document &doc, const io::Document::Section &s, const std::vector<std::string> &body) {
  doc.ReplaceLines(s.open + 1, s.close, body);
}

std::string Trim(const std::string &s) {
  const std::size_t a = s.find_first_not_of(" \t");
  if (a == std::string::npos) return "";
  const std::size_t b = s.find_last_not_of(" \t\r");
  return s.substr(a, b - a + 1);
}

}  // namespace

std::string ReadAnswersJson(const std::string &text, std::map<std::string, std::string> &answers) {
  std::size_t i = 0;
  auto blanks = [&]() {
    while (i < text.size() && (text[i] == ' ' || text[i] == '\t' || text[i] == '\n' || text[i] == '\r')) i++;
  };
  auto string = [&](std::string &out) {
    if (i >= text.size() || text[i] != '"') return false;
    for (i++; i < text.size() && text[i] != '"'; i++) {
      if (text[i] == '\\') {
        if (++i >= text.size()) return false;
        const char c = text[i];
        out += c == 'n' ? '\n' : c == 't' ? '\t' : c;  // \" \\ \/ and the two control escapes a value can need
      } else {
        out += text[i];
      }
    }
    if (i >= text.size()) return false;
    i++;
    return true;
  };
  blanks();
  if (i >= text.size() || text[i] != '{') return "the answers are not a JSON object";
  i++;
  blanks();
  if (i < text.size() && text[i] == '}') return "";
  while (true) {
    blanks();
    std::string key, value;
    if (!string(key)) return "a key that is not a JSON string near byte " + std::to_string(i);
    blanks();
    if (i >= text.size() || text[i] != ':') return "no ':' after \"" + key + "\"";
    i++;
    blanks();
    if (i < text.size() && text[i] == '"') {
      if (!string(value)) return "the value of \"" + key + "\" is not a closed string";
    } else {
      const std::size_t a = i;
      while (i < text.size() && text[i] != ',' && text[i] != '}' && text[i] != ' ' && text[i] != '\n') i++;
      value = text.substr(a, i - a);
      if (value.empty() || value == "null" || value[0] == '[' || value[0] == '{')
        return "the value of \"" + key + "\" must be a string or a number";
    }
    answers[key] = value;
    blanks();
    if (i < text.size() && text[i] == ',') {
      i++;
      continue;
    }
    if (i < text.size() && text[i] == '}') return "";
    return "expected ',' or '}' after the value of \"" + key + "\"";
  }
}

std::string RenameParameterFile(const std::string &text, const std::vector<std::pair<std::string, std::string>> &names,
                                std::vector<std::pair<std::string, std::string>> &renamed) {
  const std::map<std::string, std::string> map(names.begin(), names.end());
  std::string out;
  std::size_t start = 0;
  while (start < text.size()) {
    std::size_t end = text.find('\n', start);
    end = end == std::string::npos ? text.size() : end + 1;
    const std::string line = text.substr(start, end - start);
    start = end;
    const std::size_t a = line.find_first_not_of(" \t");
    const std::size_t b = a == std::string::npos ? a : line.find_first_of(" \t\r\n", a);
    const auto it = a == std::string::npos ? map.end() : map.find(line.substr(a, b - a));
    if (it == map.end()) {
      out += line;
      continue;
    }
    renamed.push_back(*it);
    // the name ends where it ended (a right-aligned column), or one blank after a line start that had blanks
    const std::size_t width = b, size = it->second.size();
    const std::size_t pad = a == 0 ? 0 : (width > size + 1 ? width - size : 1);
    out += std::string(pad, ' ') + it->second + line.substr(b);
  }
  return out;
}

ParameterFileConversion ConvertParameterFile(const std::string &parameterText, const std::string &oldProjectText,
                                             const std::string &projectDir) {
  ParameterFileConversion out;
  if (io::FormatOf(io::Document::Parse(oldProjectText)) == 3) {
    out.refusal = "the project is format 3; give the old project the parameter file was written for (the converter "
                  "keeps it unchanged): the names come from its conversion";
    return out;
  }
  ConvertOptions o;
  o.projectDir = projectDir;
  o.allowQuestions = true;
  const ConvertResult r = ConvertToFormat3(oldProjectText, o);
  if (!r.ok) {
    out.refusal = "the project is not converted: " + r.refusal;
    return out;
  }
  out.text = RenameParameterFile(parameterText, r.parameterNames, out.renamed);
  return out;
}

namespace {

/// True when the old setup program's check box "use GSL Coulomb functions" was ticked in a format-1 project: the first
/// number of <lastRun> is the old program's option mask (AZURESetup::writeLastRun writes `paramMask` first, readLastRun
/// reads it with `>>`), and the box is bit 12 (Config::USE_GSL_COULOMB_FUNC = 1 << 12 = 4096) of it.
bool OldGslCoulombBoxTicked(const io::Document &doc) {
  if (io::FormatOf(doc) != 1) return false;  // ConvertToFormat2 removes <lastRun> from a format-1 file only
  const io::Document::Section *s = doc.Find("lastRun");
  if (!s || s->close == io::Document::Section::npos) return false;
  for (const std::string &line : doc.Body(*s)) {
    const std::size_t a = line.find_first_not_of(" \t\r");
    if (a == std::string::npos) continue;
    const std::size_t b = line.find_first_of(" \t\r", a);
    const std::string token = line.substr(a, b == std::string::npos ? std::string::npos : b - a);
    if (token.empty() || token.find_first_not_of("0123456789") != std::string::npos || token.size() > 10) return false;
    return (std::stoull(token) >> 12) & 1u;
  }
  return false;
}

}  // namespace

ConvertResult ConvertToFormat3(const std::string &givenText, const ConvertOptions &options) {
  ConvertResult res;
  // A UTF-8 byte-order mark (written by some Windows editors) before the first tag: no AZURE2 wrote one, and the
  // section reader would not see the first tag behind it.  Removed with a note; the lines keep their numbers.
  const bool bom = givenText.compare(0, 3, "\xEF\xBB\xBF") == 0;
  std::string oldText = bom ? givenText.substr(3) : givenText;
  std::vector<std::pair<std::string, std::string>> eraNotes;  // forms of earlier eras brought to a later one (OldEra.h)
  CloseOldLastRun(oldText, eraNotes);
  if (bom) {
    res.notes.push_back({"byte-order mark", "", "the file begins with a UTF-8 byte-order mark (a Windows editor's): removed"});
    res.constructs["byte-order mark"].converted++;
    res.constructs["byte-order mark"].notes++;
  }
  const io::Document given = io::Document::Parse(oldText);  // the lines as numbered in the file (for a refusal)
  io::Document doc = io::Document::Parse(oldText);
  if (io::FormatOf(doc) == 3) {
    res.refusal = "already format 3";
    return res;
  }
  for (const std::string &p : doc.Problems()) {
    res.refusal = "the project's sections cannot be read: " + p;
    return res;
  }
  if (io::FormatOf(doc) == 1) {
    res.refusal = FoldOldExternalCapture(doc, eraNotes);
    if (!res.refusal.empty()) return res;
  }
  for (const auto &n : eraNotes) {
    res.notes.push_back({n.first, "", n.second});
    res.constructs[n.first].converted++;
    res.constructs[n.first].notes++;
  }
  res.refusal = UnconvertedSection(doc);
  if (!res.refusal.empty()) return res;
  const bool gslCoulomb = OldGslCoulombBoxTicked(doc);  // before ConvertToFormat2 removes <lastRun>
  res.format2Changes = ConvertToFormat2(doc);

  OldProject old;
  std::string why = ReadOldProject(doc, old, options.projectDir);
  for (const auto &r : old.repairs) {
    res.notes.push_back({r.first, "", r.second});
    res.constructs[r.first].converted++;
    res.constructs[r.first].notes++;
  }
  if (!why.empty()) {
    res.refusal = why;
    // The refused line as numbered in the file: format 2's edits (ConvertToFormat2) leave these sections alone, so the
    // line keeps its place in its section.
    const io::Document::Section *now = old.refusedSection.empty() ? nullptr : doc.Find(old.refusedSection);
    const io::Document::Section *then = now ? given.Find(old.refusedSection) : nullptr;
    if (now && then && old.refusedLine > now->open) {
      LineProblem &p = res.lineProblem;
      p.section = old.refusedSection;
      const std::size_t at = then->open + (old.refusedLine - now->open);
      p.fileLine = (int)at + 1;
      p.text = given.Line(at);
      if (!p.text.empty() && p.text.back() == '\r') p.text.pop_back();
      p.problem = why;
      p.expected = ExpectedFields(p.section);
      p.readings = ProposeReadings(p.section, p.text);
      res.refusal += "\n" + p.Text();
    }
    return res;
  }
  Built built;
  why = BuildFormat3(doc, old, options, built, res);
  if (!why.empty()) {
    res.refusal = why;
    return res;
  }

  if (gslCoulomb) {
    const std::string construct = "<lastRun> use GSL Coulomb functions";
    res.notes.push_back(
        {construct, "",
         "the old program's check box \"use GSL Coulomb functions\" was ticked in <lastRun> (dropped with that section): "
         "it computed this project with GSL's Coulomb functions.  This program has its own Coulomb functions only "
         "(docs/upgrade/coulomb.md, sections 8.3 and 9.1: GSL's are wrong by up to 5 % in the penetrability inside the "
         "outer turning point), so results inside the Coulomb barrier can differ from the old ones by up to about 5 % in "
         "a penetrability"});
    res.constructs[construct].converted++;
    res.constructs[construct].notes++;
  }

  // The sections: each new one where its old one was (targets before data when there was no <targetInt>).
  const std::vector<std::string> data = io::WriteRecords(built.data), targets = io::WriteRecords(built.targets),
                                 extrapolations = io::WriteRecords(built.extrapolations),
                                 modelData = io::WriteRecords(built.modelData);
  if (!ReplaceSection(doc, "targetInt", "targets", targets) && !built.targets.records.empty())
    InsertSection(doc, "segmentsData", "targets", targets);
  if (!ReplaceSection(doc, "segmentsData", "data", data)) InsertSection(doc, "segmentsTest", "data", data);
  if (!ReplaceSection(doc, "segmentsTest", "extrapolations", extrapolations) && !built.extrapolations.records.empty())
    InsertSection(doc, "", "extrapolations", extrapolations);
  ReplaceSection(doc, "abundances", "", {});
  // One <modelData> section: a format-2 file may hold one already (format-3 datum records, copied as they are); the
  // converted <scatteringLengths> rows follow its records.
  const io::Document::Section *existing = doc.Find("modelData");
  if (existing && existing->close != io::Document::Section::npos && !built.modelData.records.empty()) {
    std::vector<std::string> body = doc.Body(*existing);
    if (!body.empty()) body.push_back("");
    body.insert(body.end(), modelData.begin(), modelData.end());
    ReplaceBody(doc, *existing, body);
    ReplaceSection(doc, "scatteringLengths", "", {});
  } else if (!ReplaceSection(doc, "scatteringLengths", built.modelData.records.empty() ? "" : "modelData", modelData) &&
             !built.modelData.records.empty()) {
    InsertSection(doc, "", "modelData", modelData);
  }
  // <potential> -> <potentials> where the block was (none when no pair gets a potential)
  ReplaceSection(doc, "potential", built.potentials.records.empty() ? "" : "potentials", io::WriteRecords(built.potentials));
  // dev's <channelRadii> lines -> `channel_radius` records, in the same section
  if (!built.channelRadii.records.empty()) ReplaceSection(doc, "channelRadii", "channelRadii", io::WriteRecords(built.channelRadii));
  // dev's <isospin> lines -> isospin_* records, in the same section
  if (!built.isospin.records.empty()) ReplaceSection(doc, "isospin", "isospin", io::WriteRecords(built.isospin));
  for (const auto &b : built.bodies) {
    if (b.first == "mcmc/parameters") {
      if (const io::Document::Section *m = doc.Find("mcmc"))
        if (const io::Document::Section *p = doc.FindChild(*m, "parameters")) ReplaceBody(doc, *p, b.second);
    } else if (const io::Document::Section *s = doc.Find(b.first)) {
      ReplaceBody(doc, *s, b.second);
    }
  }
  // <config>: key=value edits (format=3; removed and renamed keys)
  if (const io::Document::Section *c = doc.Find("config")) {
    std::vector<std::string> body;
    for (const std::string &line : doc.Body(*c)) {
      const std::size_t eq = line.find('=');
      const std::size_t hash = line.find('#');
      if (eq != std::string::npos && (hash == std::string::npos || eq < hash)) {
        const auto e = built.configEdits.find(Trim(line.substr(0, eq)));
        if (e != built.configEdits.end()) {
          if (!e->second.empty()) body.push_back(Trim(line.substr(0, eq)) + "=" + e->second);
          continue;
        }
      }
      body.push_back(line);
    }
    ReplaceBody(doc, *doc.Find("config"), body);
  }
  res.text = doc.Write();

  // Validate the output with the table (the converter never returns a file with errors).
  experiment::Format3ReadOptions read;
  read.build = false;  // the build needs the data files, which the conversion reads but the fixtures do not carry
  const experiment::Format3Project check = experiment::ReadFormat3Text(res.text, read);
  bool errors = false;
  for (const experiment::Diagnostic &d : check.diagnostics) {
    if (d.severity != experiment::Severity::Error) continue;
    if (options.allowQuestions && d.code == std::string(experiment::codes::kPlaceholder)) continue;
    errors = true;
  }
  res.diagnostics = check.diagnostics;
  if (errors) {
    res.refusal = "the converted project does not validate (" + std::to_string(check.diagnostics.size()) + " diagnostic(s))";
    return res;
  }
  // The final states of external capture, as the run's set-up checks them (experiment::ExternalCaptureFinalStates through
  // the engine's build; data files not required): a project the old program ran can be one this program refuses at
  // its first run.  Before 2026-10-26 (35793eb3) the engine created a final state that is not a level of <levels> as a
  // hidden level at the pair's E_x, with a free reduced width of 0.1 MeV^1/2 in each particle channel that no line of
  // the project held; that level is not written by guess here: the conversion stops with the old reason.
  {
    experiment::Format3ReadOptions build;
    build.projectDir = options.projectDir;
    build.dataFiles = experiment::DataFiles::NotUsed;
    const experiment::Format3Project built = experiment::ReadFormat3Text(res.text, build);
    for (const experiment::Diagnostic &d : built.diagnostics) {
      const std::string code = d.code;
      if (code != experiment::codes::kFinalStateNotLevel && code != experiment::codes::kFinalStateShared &&
          code != experiment::codes::kFinalStateUnbound)
        continue;
      res.diagnostics.push_back(d);
      res.refusal = d.message + " (" + code + ").  An AZURE2 before 2026-10-26 (35793eb3) ran such a project: it created "
                    "the final state as a hidden level at the pair's E_x, with a free reduced width of 0.1 MeV^1/2 in "
                    "each particle channel, which no line of the project held (its fit parameters moved the names of "
                    "every later level).  This program needs the final state as a level of <levels>: add a level of "
                    "that J^pi at that E_x (its widths fixed, or the values you mean) to a copy of the old project and "
                    "convert the copy, or switch the external capture of that pair off";
      return res;
    }
  }
  res.ok = true;
  return res;
}

std::string ConvertResult::ReportJson(const std::string &input, const std::string &output) const {
  experiment::JsonWriter j(2);
  j.BeginObject();
  j.Key("input").String(input);
  j.Key("output").String(output);
  j.Key("ok").Bool(ok);
  if (!refusal.empty()) j.Key("refusal").String(refusal);
  j.Key("format2_changes").StringArray(format2Changes);
  j.Key("constructs").BeginObject();
  for (const auto &c : constructs) {
    j.Key(c.first).BeginObject();
    j.Key("converted").Integer(c.second.converted);
    j.Key("notes").Integer(c.second.notes);
    j.Key("questions").Integer(c.second.questions);
    j.EndObject();
  }
  j.EndObject();
  j.Key("notes").BeginArray();
  for (const ConvertNote &n : notes) {
    j.BeginObject();
    j.Key("construct").String(n.construct);
    j.Key("record").String(n.record);
    j.Key("text").String(n.text);
    j.EndObject();
  }
  j.EndArray();
  j.Key("questions").BeginArray();
  for (const ConvertQuestion &q : questions) {
    j.BeginObject();
    j.Key("severity").String("ask");
    j.Key("construct").String(q.construct);
    j.Key("key").String(q.key);
    j.Key("message").String(q.text);
    j.Key("proposal").String(q.proposal);
    j.EndObject();
  }
  j.EndArray();
  j.Key("parameters").BeginObject();
  for (const auto &p : parameterNames) j.Key(p.first).String(p.second);
  j.EndObject();
  j.Key("side_files").BeginArray();
  for (const auto &f : sideFiles) j.String(f.first);
  j.EndArray();
  j.Key("diagnostics").BeginArray();
  j.EndArray();
  j.EndObject();
  // the diagnostics as their own JSON objects (Diagnostic::ToJson), inserted into the empty array above
  std::string s = j.Str();
  std::string list;
  for (std::size_t i = 0; i < diagnostics.size(); i++) list += (i ? ",\n    " : "\n    ") + diagnostics[i].ToJson();
  if (!list.empty()) {
    const std::string empty = "\"diagnostics\": []";
    const std::size_t at = s.rfind(empty);
    if (at != std::string::npos) s.replace(at, empty.size(), "\"diagnostics\": [" + list + "\n  ]");
  }
  return s + "\n";
}

}  // namespace convert
}  // namespace azure2
