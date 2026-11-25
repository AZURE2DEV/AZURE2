// azure2::front::Project (Engine.h): the project file as io::Document, edited with the engine's grammars of <config>
// (io::ConfigEdits), <levels> (io::LevelRow, io::EditLevelRow) and the record sections (io::ReadRecords,
// io::WriteRecords), checked by ValidateProjectBytes (AZURE2 --validate).  docs/upgrade/qt_gui.md section 2.
#include "azure2/front/Engine.h"
#include "azure2/support/ClassicNumeric.h"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>

#include "ProjectCommands.h"
#include "azure2/experiment/ProjectPaths.h"
#include "azure2/experiment/Schema.h"
#include "azure2/io/ConfigSection.h"
#include "azure2/io/Document.h"
#include "azure2/io/LevelRowText.h"
#include "azure2/io/LevelRows.h"
#include "azure2/io/ProjectFormat.h"
#include "azure2/io/RecordReader.h"

namespace azure2 {
namespace front {

namespace {

namespace ex = azure2::experiment;

bool ReadBytes(const std::string &path, std::string &bytes) {
  std::ifstream in(path, std::ios::binary);
  if (!in) return false;
  std::ostringstream s;
  s << in.rdbuf();
  bytes = s.str();
  return true;
}

std::string Trimmed(const std::string &s) {
  const std::size_t a = s.find_first_not_of(" \t");
  if (a == std::string::npos) return std::string();
  const std::size_t b = s.find_last_not_of(" \t");
  return s.substr(a, b - a + 1);
}

/// The text of a field's value as it stands in its line (as pybind_project.cpp's ValueText): first token to the end
/// of the last, quotes included.
std::string ValueText(const io::Document &doc, const io::Field &f) {
  if (f.tokens.empty() || f.line < 1 || (std::size_t)f.line > doc.NumLines()) return std::string();
  const std::string &line = doc.Line(f.line - 1);
  const io::Token &last = f.tokens.back();
  const std::size_t from = f.tokens.front().col - 1, to = last.col - 1 + last.text.size() + (last.quoted ? 2 : 0);
  return to <= line.size() && from <= to ? line.substr(from, to - from) : std::string();
}

/// The tokens of a value text (the record grammar's lexer); false for a quote that is not closed.
bool TokensOf(const std::string &value, std::vector<io::Token> &tokens, std::string &why) {
  std::string comment;
  std::vector<io::RecordDiagnostic> problems;
  tokens = io::Tokenize(value, 0, comment, problems);
  if (!problems.empty()) {
    why = problems.front().message;
    return false;
  }
  if (!comment.empty()) {
    why = "a value holds no '#' (it starts a comment): " + value;
    return false;
  }
  return true;
}

bool IsToken(const std::string &s) {
  if (s.empty()) return false;
  for (char c : s)
    if (c == ' ' || c == '\t' || c == '#' || c == '"') return false;
  return true;
}

}  // namespace

std::string Record::Value(const std::string &key) const {
  for (const Field &f : fields)
    if (f.key == key) return f.value;
  return std::string();
}

struct Project::Impl {
  std::string path;
  io::Document doc;
  std::string saved;  ///< the bytes last read or written

  const io::Document::Section *Find(const std::string &name) const { return doc.Find(name); }

  /// The records of a section with the reader's own result.
  io::RecordReadResult Read(const std::string &section) const {
    const io::Document::Section *s = Find(section);
    if (!s || s->close == io::Document::Section::npos) return io::RecordReadResult();
    return io::ReadRecords(doc, *s, ex::Schema::KindsOfSection(section));
  }

  /// The line indices [first, last] of a record: its comment lines above the header to its last field line.
  static void Span(const io::Record &r, std::size_t &first, std::size_t &last) {
    int lastLine = r.line;
    for (const io::Field &f : r.fields) lastLine = std::max(lastLine, f.line);
    first = (std::size_t)(r.line - 1 - (int)r.comments.size());
    last = (std::size_t)(lastLine - 1);
  }

  /// Make a top-level section at the end of the file when it has none.
  const io::Document::Section *Ensure(const std::string &name) {
    if (const io::Document::Section *s = Find(name)) return s;
    std::size_t at = doc.NumLines();
    // a last line without its ending gets one before the new section (ReplaceLines gives each new line the ending of
    // the first line it replaces, the common one here)
    if (at > 0 && doc.Ending(at - 1).empty()) {
      const std::string last = doc.Line(at - 1);
      doc.ReplaceLines(at - 1, at, {last, "<" + name + ">", "</" + name + ">"});
    } else {
      doc.ReplaceLines(at, at, {"<" + name + ">", "</" + name + ">"});
    }
    return Find(name);
  }
};

Project::Project() : impl_(new Impl) {}
Project::~Project() = default;
Project::Project(Project &&) noexcept = default;
Project &Project::operator=(Project &&) noexcept = default;

bool Project::Open(const std::string &path, Project &out, std::string &why) {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  std::string bytes;
  if (!ReadBytes(path, bytes)) {
    why = "cannot read " + path;
    return false;
  }
  const std::string refusal = io::NotFormat3Message(bytes, path);
  if (!refusal.empty()) {
    why = refusal;
    return false;
  }
  out = FromText(bytes, path);
  return true;
}

Project Project::FromText(const std::string &text, const std::string &path) {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  Project p;
  p.impl_->path = path;
  p.impl_->doc = io::Document::Parse(text);
  p.impl_->saved = text;
  return p;
}

Project Project::Empty(const std::string &path) {
  std::string text = "<config>\n";
  for (const io::ConfigPositional &c : io::ConfigPositionals()) {
    const std::string name = c.name;
    if (name == "output_dir") text += "output/\n";
    else if (name == "checks_dir") text += "checks/\n";
    else {
      std::istringstream values(c.values);
      std::string first;
      values >> first;
      text += first + "\n";
    }
  }
  text += "format=3\n</config>\n<levels>\n</levels>\n<data>\n</data>\n<extrapolations>\n</extrapolations>\n"
          "<targets>\n</targets>\n";
  Project p = FromText(text, path);
  p.impl_->saved.clear();  // never written
  return p;
}

const std::string &Project::Path() const { return impl_->path; }
std::string Project::Directory() const { return ex::ProjectDirectoryOf(impl_->path); }
std::string Project::Text() const { return impl_->doc.Write(); }
bool Project::Modified() const { return Text() != impl_->saved; }

bool Project::Save(const std::string &path, std::string &why) {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  const std::string to = path.empty() ? impl_->path : path;
  const std::string text = Text();
  // a whole file or none: a temporary file beside it, renamed over it
  const std::string temporary = to + ".saving";
  {
    std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
    if (!out) {
      why = "cannot write " + temporary;
      return false;
    }
    out << text;
    if (!out.flush()) {
      why = "cannot write " + temporary;
      return false;
    }
  }
  std::remove(to.c_str());  // rename over an existing file (Windows)
  if (std::rename(temporary.c_str(), to.c_str()) != 0) {
    why = "cannot rename " + temporary + " to " + to;
    return false;
  }
  impl_->path = to;
  impl_->saved = text;
  return true;
}

std::vector<std::string> Project::Sections() const {
  std::vector<std::string> out;
  for (const io::Document::Section &s : impl_->doc.Sections())
    if (s.depth == 0) out.push_back(s.name);
  return out;
}

bool Project::HasSection(const std::string &name) const { return impl_->Find(name) != nullptr; }

std::vector<std::string> Project::SectionBody(const std::string &name) const {
  const io::Document::Section *s = impl_->Find(name);
  return s ? impl_->doc.Body(*s) : std::vector<std::string>();
}

void Project::SetSectionBody(const std::string &name, const std::vector<std::string> &lines) {
  const io::Document::Section *s = impl_->Find(name);
  if (!s && lines.empty()) return;
  if (!s) s = impl_->Ensure(name);
  if (s->close == io::Document::Section::npos) return;  // an unclosed section is reported by Validate
  if (impl_->doc.Body(*s) == lines) return;              // nothing changes: the bytes stay
  impl_->doc.ReplaceLines(s->open + 1, s->close, lines);
}

std::vector<ConfigEntry> Project::Config() const {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  std::size_t open = 0;
  const io::ConfigSection section = io::ReadConfigSection(io::ConfigLinesOf(impl_->doc, open));
  std::vector<ConfigEntry> out;
  for (const io::ConfigEntry &e : section.entries) out.push_back({e.name, e.value, e.positional});
  return out;
}

std::string Project::ConfigValue(const std::string &name, const std::string &absent) const {
  std::size_t open = 0;
  const io::ConfigSection section = io::ReadConfigSection(io::ConfigLinesOf(impl_->doc, open));
  const io::ConfigEntry *e = section.Find(name);
  return e ? e->value : absent;
}

bool Project::SetConfig(const std::string &name, const std::string &value, std::string &why, bool remove) {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  std::size_t open = 0;
  const std::vector<std::string> lines = io::ConfigLinesOf(impl_->doc, open);
  if (lines.empty()) {
    why = "the project has no <config> section";
    return false;
  }
  std::vector<io::ConfigLineEdit> edits;
  if (!io::ConfigEdits(lines, name, value, remove ? io::ConfigChange::Remove : io::ConfigChange::Set, edits, why))
    return false;
  // the edits name lines of `lines` (0 = the tag); applied from the last so that earlier indices hold
  std::stable_sort(edits.begin(), edits.end(),
                   [](const io::ConfigLineEdit &a, const io::ConfigLineEdit &b) { return a.line > b.line; });
  for (const io::ConfigLineEdit &e : edits) {
    const std::size_t at = open + (std::size_t)e.line;
    if (e.remove) impl_->doc.ReplaceLines(at, at + 1, {});
    else if (e.insert) impl_->doc.ReplaceLines(at, at, {e.text});
    else if (impl_->doc.Line(at) != e.text) impl_->doc.ReplaceLines(at, at + 1, {e.text});
  }
  return true;
}

std::vector<LevelRow> Project::Levels() const {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  std::vector<LevelRow> out;
  for (const std::string &line : SectionBody("levels")) {
    LevelRow row;
    io::LevelRowText text;
    std::string why;
    if (!Trimmed(line).empty() && io::SplitLevelRow(line, text, why))
      for (std::size_t c = 0; c < io::LevelColumnNames().size(); c++) row.tokens.push_back(text.Token((int)c));
    out.push_back(row);
  }
  return out;
}

bool Project::SetLevels(const std::vector<LevelRow> &rows, std::string &why) {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  const std::vector<std::string> old = SectionBody("levels");
  std::vector<std::string> lines;
  for (std::size_t i = 0; i < rows.size(); i++) {
    const LevelRow &row = rows[i];
    if (row.tokens.empty()) {
      lines.push_back(i < old.size() && Trimmed(old[i]).empty() ? old[i] : std::string());
      continue;
    }
    io::LevelRowText was;
    std::string notRow;
    std::string line;
    if (i < old.size() && io::SplitLevelRow(old[i], was, notRow)) {
      std::map<int, std::string> changed;
      for (std::size_t c = 0; c < row.tokens.size() && c < io::LevelColumnNames().size(); c++)
        if (row.tokens[c] != was.Token((int)c) && !(row.tokens[c].empty() && was.Token((int)c).empty()))
          changed[(int)c] = row.tokens[c];
      line = changed.empty() ? old[i] : io::EditLevelRow(was, changed);
    } else {
      // a new row: the 31 columns, then bindingE and gammaIsRWA together, then e1, each written when given
      const std::vector<std::string> &names = io::LevelColumnNames();
      for (std::size_t c = 0; c < row.tokens.size() && c < names.size(); c++) {
        if (row.tokens[c].empty()) continue;
        line += (line.empty() ? "" : " ") + row.tokens[c];
      }
    }
    // the engine's reader must take the row
    std::istringstream in(line);
    const io::LevelRow read(in);
    if (in.fail()) {
      why = "<levels> row " + std::to_string(i + 1) + " is not a row the engine reads: " + line;
      return false;
    }
    lines.push_back(line);
  }
  if (!HasSection("levels")) impl_->Ensure("levels");
  SetSectionBody("levels", lines);
  return true;
}

std::vector<Record> Project::Records(const std::string &section) const {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  std::vector<Record> out;
  const io::RecordReadResult read = impl_->Read(section);
  for (const io::Record &r : read.set.records) {
    Record x;
    x.kind = r.kind;
    x.name = r.name;
    x.line = r.line;
    for (const io::Field &f : r.fields) x.fields.push_back({f.key, ValueText(impl_->doc, f), f.trailingComment});
    out.push_back(x);
  }
  return out;
}

bool Project::SetRecord(const std::string &section, const Record &record, std::string &why) {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  const std::vector<std::string> kinds = ex::Schema::KindsOfSection(section);
  if (std::find(kinds.begin(), kinds.end(), record.kind) == kinds.end()) {
    why = "<" + section + "> holds no record of kind '" + record.kind + "'";
    return false;
  }
  if (!IsToken(record.name)) {
    why = "a record's name is one token without blanks, quotes or '#': '" + record.name + "'";
    return false;
  }
  io::Record made;
  made.kind = record.kind;
  made.name = record.name;
  for (const Field &f : record.fields) {
    if (!IsToken(f.key)) {
      why = record.name + ": a field's key is one token: '" + f.key + "'";
      return false;
    }
    io::Field field;
    field.key = f.key;
    if (!TokensOf(f.value, field.tokens, why)) {
      why = record.name + "." + f.key + ": " + why;
      return false;
    }
    field.trailingComment = f.comment;
    made.fields.push_back(field);
  }
  const io::RecordReadResult read = impl_->Read(section);
  const io::Record *old = nullptr;
  for (const io::Record &r : read.set.records)
    if (r.name == record.name) old = &r;
  if (old) {
    // the same keys in the same order: only the lines whose value changed are rewritten, from the value's column on
    bool sameKeys = old->kind == made.kind && old->fields.size() == made.fields.size();
    for (std::size_t i = 0; sameKeys && i < made.fields.size(); i++)
      sameKeys = old->fields[i].key == made.fields[i].key && !old->fields[i].tokens.empty() && !made.fields[i].tokens.empty();
    if (sameKeys) {
      for (std::size_t i = 0; i < made.fields.size(); i++) {
        const io::Field &o = old->fields[i];
        if (io::WriteTokens(o.tokens) == io::WriteTokens(made.fields[i].tokens) && o.trailingComment == made.fields[i].trailingComment)
          continue;
        const std::string &line = impl_->doc.Line(o.line - 1);
        std::string text = line.substr(0, o.tokens.front().col - 1) + io::WriteTokens(made.fields[i].tokens);
        if (!made.fields[i].trailingComment.empty()) text += "  " + made.fields[i].trailingComment;
        impl_->doc.ReplaceLines(o.line - 1, o.line, {text});
      }
      return true;
    }
    made.comments = old->comments;
    made.trailingComment = old->trailingComment;
    for (io::Field &f : made.fields)
      for (const io::Field &o : old->fields)
        if (o.key == f.key) {
          f.comments = o.comments;
          break;
        }
    std::size_t first = 0, last = 0;
    Impl::Span(*old, first, last);
    io::RecordSet one;
    one.records.push_back(made);
    impl_->doc.ReplaceLines(first, last + 1, io::WriteRecords(one));
    return true;
  }
  const io::Document::Section *s = impl_->Ensure(section);
  if (!s || s->close == io::Document::Section::npos) {
    why = "<" + section + "> is not closed";
    return false;
  }
  io::RecordSet one;
  one.records.push_back(made);
  std::vector<std::string> lines = io::WriteRecords(one);
  // a blank line between the last record (or the section's last line) and the new one
  if (s->close > s->open + 1 && !Trimmed(impl_->doc.Line(s->close - 1)).empty()) lines.insert(lines.begin(), "");
  impl_->doc.ReplaceLines(s->close, s->close, lines);
  return true;
}

bool Project::RenameRecord(const std::string &section, const std::string &from, const std::string &to, std::string &why) {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  if (!IsToken(to)) {
    why = "a record's name is one token without blanks, quotes or '#': '" + to + "'";
    return false;
  }
  const io::RecordReadResult read = impl_->Read(section);
  for (const io::Record &r : read.set.records)
    if (r.name == to) {
      why = "<" + section + "> has a record named " + to;
      return false;
    }
  for (const io::Record &r : read.set.records) {
    if (r.name != from) continue;
    const std::string &line = impl_->doc.Line(r.line - 1);
    const std::size_t at = (std::size_t)r.nameCol - 1;
    const std::size_t width = from.size() + (line.compare(at, 1, "\"") == 0 ? 2 : 0);
    impl_->doc.ReplaceLines(r.line - 1, r.line, {line.substr(0, at) + to + line.substr(at + width)});
    return true;
  }
  why = "<" + section + "> has no record named " + from;
  return false;
}

bool Project::RemoveRecord(const std::string &section, const std::string &name, std::string &why) {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  const io::RecordReadResult read = impl_->Read(section);
  for (const io::Record &r : read.set.records) {
    if (r.name != name) continue;
    std::size_t first = 0, last = 0;
    Impl::Span(r, first, last);
    // the blank line that separated it from the record before goes with it
    if (first > 0 && Trimmed(impl_->doc.Line(first - 1)).empty()) first--;
    impl_->doc.ReplaceLines(first, last + 1, {});
    return true;
  }
  why = "<" + section + "> has no record named " + name;
  return false;
}

bool Project::MoveRecord(const std::string &section, const std::string &name, int delta, std::string &why) {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  const io::RecordReadResult read = impl_->Read(section);
  const std::vector<io::Record> &records = read.set.records;
  std::size_t at = records.size();
  for (std::size_t i = 0; i < records.size(); i++)
    if (records[i].name == name) at = i;
  if (at == records.size()) {
    why = "<" + section + "> has no record named " + name;
    return false;
  }
  if ((delta < 0 && at == 0) || (delta > 0 && at + 1 >= records.size()) || delta == 0) {
    why = name + " is already the " + std::string(delta < 0 ? "first" : "last") + " record of <" + section + ">";
    return false;
  }
  const io::Record &upper = records[delta < 0 ? at - 1 : at], &lower = records[delta < 0 ? at : at + 1];
  std::size_t a1 = 0, a2 = 0, b1 = 0, b2 = 0;
  Impl::Span(upper, a1, a2);
  Impl::Span(lower, b1, b2);
  // the lower record's lines, the lines between the two (blank lines, comments of no record), the upper record's
  std::vector<std::string> lines;
  for (std::size_t i = b1; i <= b2; i++) lines.push_back(impl_->doc.Line(i));
  for (std::size_t i = a2 + 1; i < b1; i++) lines.push_back(impl_->doc.Line(i));
  for (std::size_t i = a1; i <= a2; i++) lines.push_back(impl_->doc.Line(i));
  impl_->doc.ReplaceLines(a1, b2 + 1, lines);
  return true;
}

std::string RecordOffMark() { return "#off "; }

namespace {

/// The runs of lines of a section body that carry the off mark, each with the line index (0-based, in the file) of its
/// first line, the mark taken away.
struct OffRun {
  std::size_t first = 0;
  std::vector<std::string> lines;
};

std::vector<OffRun> OffRuns(const io::Document &doc, const io::Document::Section &s) {
  std::vector<OffRun> runs;
  const std::string mark = RecordOffMark();
  for (std::size_t i = s.open + 1; i < s.close; i++) {
    const std::string &line = doc.Line(i);
    if (line.compare(0, mark.size(), mark) != 0) continue;
    if (runs.empty() || runs.back().first + runs.back().lines.size() != i) runs.push_back({i, {}});
    runs.back().lines.push_back(line.substr(mark.size()));
  }
  return runs;
}

}  // namespace

std::vector<Record> Project::RecordsOff(const std::string &section) const {
  std::vector<Record> out;
  const io::Document::Section *s = impl_->Find(section);
  if (!s || s->close == io::Document::Section::npos) return out;
  for (const OffRun &run : OffRuns(impl_->doc, *s)) {
    const io::RecordReadResult read = io::ReadRecords(run.lines, (int)run.first + 1, ex::Schema::KindsOfSection(section));
    for (const io::Record &r : read.set.records) {
      Record x;
      x.kind = r.kind;
      x.name = r.name;
      x.line = r.line;
      for (const io::Field &f : r.fields) {
        // the value as written: the line from the first token's column (the mark is not in these lines' columns)
        const std::string &line = run.lines[(std::size_t)f.line - run.first - 1];
        std::string value = f.tokens.empty() ? std::string() : line.substr((std::size_t)f.tokens.front().col - 1);
        if (!f.trailingComment.empty()) {
          const std::size_t at = value.rfind(f.trailingComment);
          if (at != std::string::npos) value = value.substr(0, at);
        }
        x.fields.push_back({f.key, Trimmed(value), f.trailingComment});
      }
      out.push_back(x);
    }
  }
  return out;
}

bool Project::SetRecordOn(const std::string &section, const std::string &name, bool on, std::string &why) {
  const io::Document::Section *s = impl_->Find(section);
  if (!s || s->close == io::Document::Section::npos) {
    why = "the project has no <" + section + "> section";
    return false;
  }
  const io::RecordReadResult read = impl_->Read(section);
  const io::Record *active = nullptr;
  for (const io::Record &r : read.set.records)
    if (r.name == name) active = &r;
  if (!on) {
    if (!active) {
      why = "<" + section + "> has no active record named " + name;
      return false;
    }
    std::size_t first = 0, last = 0;
    Impl::Span(*active, first, last);
    std::vector<std::string> lines;
    for (std::size_t i = first; i <= last; i++) lines.push_back(RecordOffMark() + impl_->doc.Line(i));
    impl_->doc.ReplaceLines(first, last + 1, lines);
    return true;
  }
  if (active) {
    why = "<" + section + "> has an active record named " + name + "; rename one of the two first";
    return false;
  }
  for (const OffRun &run : OffRuns(impl_->doc, *s)) {
    const io::RecordReadResult off = io::ReadRecords(run.lines, (int)run.first + 1, ex::Schema::KindsOfSection(section));
    for (const io::Record &r : off.set.records) {
      if (r.name != name) continue;
      std::size_t first = 0, last = 0;
      Impl::Span(r, first, last);
      std::vector<std::string> lines;
      for (std::size_t i = first; i <= last; i++) lines.push_back(impl_->doc.Line(i).substr(RecordOffMark().size()));
      impl_->doc.ReplaceLines(first, last + 1, lines);
      return true;
    }
  }
  why = "<" + section + "> has no record named " + name + " switched off";
  return false;
}

Validation Project::Validate() const {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  Validation out;
  const ProjectValidation v = ValidateProjectBytes(Text(), impl_->path);
  if (!v.notFormat3.empty()) {
    out.format3 = false;
    out.refusal = v.notFormat3;
    return out;
  }
  out.ok = v.ok;
  out.stoppingReport = v.stoppingText;
  for (const ex::Diagnostic &d : v.diagnostics) {
    Diagnostic x;
    x.severity = d.severity == ex::Severity::Error ? "error" : d.severity == ex::Severity::Warning ? "warning" : "ask";
    x.code = d.code;
    x.file = d.file;
    x.line = d.line;
    x.section = d.section;
    x.record = d.record;
    x.field = d.field;
    x.message = d.message;
    x.text = d.ToText();
    out.diagnostics.push_back(x);
  }
  return out;
}

}  // namespace front
}  // namespace azure2
