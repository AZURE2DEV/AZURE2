// The .azr file as a lossless document (Document.h; docs/upgrade/project.md section 2).
#include "azure2/io/Document.h"

#include <cctype>
#include <map>

namespace azure2 {
namespace io {

namespace {

std::string Trim(const std::string &s) {
  std::size_t a = 0, b = s.size();
  while (a < b && (s[a] == ' ' || s[a] == '\t')) a++;
  while (b > a && (s[b - 1] == ' ' || s[b - 1] == '\t')) b--;
  return s.substr(a, b - a);
}

}  // namespace

std::string TagOf(const std::string &line, bool &closing) {
  const std::string t = Trim(line);
  closing = false;
  if (t.size() < 3 || t.front() != '<' || t.back() != '>') return "";
  std::size_t a = 1;
  if (t[1] == '/') {
    closing = true;
    a = 2;
  }
  const std::string name = t.substr(a, t.size() - 1 - a);
  if (name.empty() || !std::isalpha(static_cast<unsigned char>(name[0]))) return "";
  for (char c : name)
    if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_') return "";
  return name;
}

Document Document::Parse(const std::string &bytes) {
  Document d;
  std::size_t i = 0;
  while (i < bytes.size()) {
    std::size_t j = i;
    while (j < bytes.size() && bytes[j] != '\n' && bytes[j] != '\r') j++;
    LineRec r;
    r.text = bytes.substr(i, j - i);
    if (j < bytes.size()) {
      if (bytes[j] == '\r' && j + 1 < bytes.size() && bytes[j + 1] == '\n') {
        r.ending = "\r\n";
        j += 2;
      } else {
        r.ending = std::string(1, bytes[j]);
        j += 1;
      }
    }
    d.lines_.push_back(r);
    i = j;
  }
  d.FindSections();
  return d;
}

std::string Document::Write() const {
  std::string out;
  for (const LineRec &r : lines_) {
    out += r.text;
    out += r.ending;
  }
  return out;
}

void Document::FindSections() {
  sections_.clear();
  problems_.clear();
  std::vector<std::size_t> open;  // indices into sections_
  for (std::size_t i = 0; i < lines_.size(); i++) {
    bool closing = false;
    const std::string name = TagOf(lines_[i].text, closing);
    if (name.empty()) continue;
    if (!closing) {
      Section s;
      s.name = name;
      s.open = i;
      s.close = Section::npos;
      s.depth = (int)open.size();
      s.parent = open.empty() ? Section::npos : open.back();
      sections_.push_back(s);
      open.push_back(sections_.size() - 1);
      continue;
    }
    // A closing tag: the innermost open section of that name.
    std::size_t k = open.size();
    while (k > 0 && sections_[open[k - 1]].name != name) k--;
    if (k == 0) {
      problems_.push_back("line " + std::to_string(i + 1) + ": </" + name + "> without <" + name + ">");
      continue;
    }
    for (std::size_t m = open.size(); m > k; m--)
      problems_.push_back("line " + std::to_string(i + 1) + ": </" + name + "> closes <" + name +
                          "> while <" + sections_[open[m - 1]].name + "> (line " +
                          std::to_string(sections_[open[m - 1]].open + 1) + ") is still open");
    sections_[open[k - 1]].close = i;
    open.resize(k - 1);
  }
  for (std::size_t k : open)
    problems_.push_back("line " + std::to_string(sections_[k].open + 1) + ": <" + sections_[k].name +
                        "> is not closed");
  // An unclosed section holds no sections: those opened after it belong to its parent, so one missing closing tag
  // does not hide the rest of the file (each reader then reports its own unterminated block).
  for (std::size_t k : open) {
    for (Section &s : sections_) {
      if (s.parent != k) continue;
      s.parent = sections_[k].parent;
    }
  }
  for (Section &s : sections_) {  // depths from the corrected parents (parents precede their children)
    s.depth = (s.parent == Section::npos) ? 0 : sections_[s.parent].depth + 1;
  }
  // A top-level section given twice: only the first is read.
  for (std::size_t a = 0; a < sections_.size(); a++)
    for (std::size_t b = a + 1; b < sections_.size(); b++)
      if (sections_[a].depth == 0 && sections_[b].depth == 0 && sections_[a].name == sections_[b].name) {
        problems_.push_back("line " + std::to_string(sections_[b].open + 1) + ": a second <" + sections_[b].name +
                            "> (the first is at line " + std::to_string(sections_[a].open + 1) +
                            "); only the first is read");
        break;
      }
}

const Document::Section *Document::Find(const std::string &name) const {
  for (const Section &s : sections_)
    if (s.depth == 0 && s.name == name) return &s;
  return nullptr;
}

const Document::Section *Document::FindChild(const Section &parent, const std::string &name) const {
  const std::size_t p = static_cast<std::size_t>(&parent - sections_.data());
  for (const Section &s : sections_)
    if (s.parent == p && s.name == name) return &s;
  return nullptr;
}

std::vector<std::string> Document::Body(const Section &s) const {
  std::vector<std::string> out;
  if (s.close == Section::npos) return out;
  for (std::size_t i = s.open + 1; i < s.close; i++) out.push_back(lines_[i].text);
  return out;
}

std::string Document::CommonEnding() const {
  std::map<std::string, int> n;
  for (const LineRec &r : lines_)
    if (!r.ending.empty()) n[r.ending]++;
  std::string best = "\n";
  int most = 0;
  for (const auto &kv : n)
    if (kv.second > most) {
      most = kv.second;
      best = kv.first;
    }
  return best;
}

void Document::ReplaceLines(std::size_t first, std::size_t last, const std::vector<std::string> &replacement) {
  if (first > last || last > lines_.size()) return;
  const std::string ending = (first < last && !lines_[first].ending.empty()) ? lines_[first].ending : CommonEnding();
  std::vector<LineRec> repl;
  for (const std::string &t : replacement) repl.push_back(LineRec{t, ending});
  // A replaced last line without an ending keeps having none.
  if (!repl.empty() && last == lines_.size() && last > first && lines_[last - 1].ending.empty()) repl.back().ending = "";
  lines_.erase(lines_.begin() + first, lines_.begin() + last);
  lines_.insert(lines_.begin() + first, repl.begin(), repl.end());
  FindSections();
}

void Document::RemoveSection(const Section &s) {
  if (s.close == Section::npos) return;
  const std::size_t first = s.open, last = s.close + 1;
  ReplaceLines(first, last, {});
}

}  // namespace io
}  // namespace azure2
