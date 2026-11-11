// The converter old -> format 3 (Convert.h): format 2 first (ConvertToFormat2, Format2.cpp), the old sections read
// (OldProject.cpp), the records built (Build.cpp), the document assembled here and validated with the table.
#include "azure2/convert/Convert.h"

#include <sstream>

#include "Build.h"
#include "OldProject.h"
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

ConvertResult ConvertToFormat3(const std::string &oldText, const ConvertOptions &options) {
  ConvertResult res;
  io::Document doc = io::Document::Parse(oldText);
  if (io::FormatOf(doc) == 3) {
    res.refusal = "already format 3";
    return res;
  }
  for (const std::string &p : doc.Problems()) {
    res.refusal = "the project's sections cannot be read: " + p;
    return res;
  }
  res.format2Changes = ConvertToFormat2(doc);

  OldProject old;
  std::string why = ReadOldProject(doc, old);
  if (!why.empty()) {
    res.refusal = why;
    return res;
  }
  Built built;
  why = BuildFormat3(doc, old, options, built, res);
  if (!why.empty()) {
    res.refusal = why;
    return res;
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
