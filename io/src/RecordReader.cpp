// The generic reader and writer of format-3 record sections (RecordReader.h; docs/upgrade/project_format3.md
// sections 3.1-3.3).
#include "azure2/io/RecordReader.h"

#include <algorithm>
#include <cctype>

namespace azure2 {
namespace io {

namespace {

bool IsBlank(char c) { return c == ' ' || c == '\t'; }

bool IsKey(const std::string &s) {
  if (s.empty() || !std::isalpha(static_cast<unsigned char>(s[0]))) return false;
  for (char c : s)
    if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_') return false;
  return true;
}

std::string TrimRight(const std::string &s) {
  std::size_t b = s.size();
  while (b > 0 && IsBlank(s[b - 1])) b--;
  return s.substr(0, b);
}

bool NeedsQuotes(const std::string &s) {
  if (s.empty()) return true;
  for (char c : s)
    if (IsBlank(c) || c == '#' || c == '"') return true;
  return false;
}

std::string TokenText(const Token &t) {
  return (t.quoted || NeedsQuotes(t.text)) ? "\"" + t.text + "\"" : t.text;
}

std::string Joined(const std::vector<Token> &tokens, std::size_t from) {
  std::string s;
  for (std::size_t i = from; i < tokens.size(); i++) s += (i > from ? " " : "") + TokenText(tokens[i]);
  return s;
}

}  // namespace

std::vector<Token> Tokenize(const std::string &line, int lineNumber, std::string &comment,
                            std::vector<RecordDiagnostic> &problems) {
  std::vector<Token> tokens;
  comment.clear();
  std::size_t i = 0;
  const std::size_t n = line.size();
  while (i < n) {
    while (i < n && IsBlank(line[i])) i++;
    if (i >= n) break;
    if (line[i] == '#') {
      comment = TrimRight(line.substr(i));
      break;
    }
    Token t;
    t.col = static_cast<int>(i) + 1;
    if (line[i] == '"') {
      t.quoted = true;
      const std::size_t close = line.find('"', i + 1);
      if (close == std::string::npos) {
        t.text = line.substr(i + 1);
        problems.push_back({record_codes::kUnterminatedQuote, lineNumber, t.col, line.substr(i),
                            "the quoted token starting at column " + std::to_string(t.col) +
                                " has no closing '\"' on its line; a quoted token ends on the line it starts"});
        i = n;
      } else {
        t.text = line.substr(i + 1, close - i - 1);
        i = close + 1;
        if (i < n && !IsBlank(line[i]) && line[i] != '#') {
          std::size_t e = i;
          while (e < n && !IsBlank(line[e])) e++;
          problems.push_back({record_codes::kQuoteInToken, lineNumber, static_cast<int>(i) + 1,
                              line.substr(t.col - 1, e - (t.col - 1)),
                              "text follows the closing '\"' without a blank; quote the whole token or separate "
                              "the tokens"});
          t.text += line.substr(i, e - i);
          i = e;
        }
      }
    } else {
      std::size_t e = i;
      while (e < n && !IsBlank(line[e]) && line[e] != '#') e++;
      t.text = line.substr(i, e - i);
      const std::size_t q = t.text.find('"');
      if (q != std::string::npos)
        problems.push_back({record_codes::kQuoteInToken, lineNumber, t.col + static_cast<int>(q), t.text,
                            "'\"' inside the token '" + t.text + "'; a quote may only open a token"});
      i = e;
    }
    tokens.push_back(t);
  }
  return tokens;
}

RecordReadResult ReadRecords(const std::vector<std::string> &lines, int firstLine, const std::vector<std::string> &kinds) {
  RecordReadResult out;
  std::vector<std::string> pending;  // whole-line comments not yet attached
  Record *current = nullptr;
  std::string kindList;
  for (std::size_t k = 0; k < kinds.size(); k++) kindList += (k ? ", " : "") + kinds[k];

  for (std::size_t li = 0; li < lines.size(); li++) {
    const int lineNumber = firstLine + static_cast<int>(li);
    std::string comment;
    std::vector<Token> tokens = Tokenize(lines[li], lineNumber, comment, out.diagnostics);
    if (tokens.empty()) {
      if (!comment.empty()) pending.push_back(comment);
      continue;
    }
    const Token &first = tokens[0];
    const bool header = !first.quoted && std::find(kinds.begin(), kinds.end(), first.text) != kinds.end();
    if (header) {
      Record r;
      r.kind = first.text;
      r.line = lineNumber;
      r.col = first.col;
      r.comments.swap(pending);
      r.trailingComment = comment;
      if (tokens.size() >= 2) {
        r.name = tokens[1].text;
        r.nameCol = tokens[1].col;
      }
      if (tokens.size() != 2)
        out.diagnostics.push_back(
            {record_codes::kHeaderName, lineNumber, tokens.size() < 2 ? first.col : tokens[2].col,
             tokens.size() < 2 ? first.text : Joined(tokens, 2),
             tokens.size() < 2 ? "the header '" + first.text + "' has no name: a record starts with '" + first.text +
                                     " <name>'"
                               : "the header '" + first.text + " " + r.name + "' is followed by '" + Joined(tokens, 2) +
                                     "' on its line: a header is '<kind> <name>' and each field goes on its own line"});
      out.set.records.push_back(r);
      current = &out.set.records.back();
      continue;
    }
    if (!current) {
      out.diagnostics.push_back({record_codes::kUnknownKind, lineNumber, first.col, first.text,
                                 "'" + first.text + "' is not a record kind of this section (" + kindList +
                                     "); a record starts with a header '<kind> <name>' and its fields follow it"});
      pending.clear();
      continue;
    }
    if (first.quoted || !IsKey(first.text)) {
      out.diagnostics.push_back({record_codes::kKeySyntax, lineNumber, first.col, first.text,
                                 "'" + first.text + "' is not a field key: a key is a letter followed by letters, "
                                 "digits and '_'"});
      pending.clear();
      continue;
    }
    Field f;
    f.key = first.text;
    f.line = lineNumber;
    f.col = first.col;
    f.tokens.assign(tokens.begin() + 1, tokens.end());
    f.comments.swap(pending);
    f.trailingComment = comment;
    current->fields.push_back(f);
  }
  out.set.trailingComments.swap(pending);
  return out;
}

RecordReadResult ReadRecords(const Document &doc, const Document::Section &section, const std::vector<std::string> &kinds) {
  return ReadRecords(doc.Body(section), static_cast<int>(section.open) + 2, kinds);
}

std::vector<std::string> WriteRecords(const RecordSet &set) {
  std::vector<std::string> out;
  auto withComment = [](std::string s, const std::string &comment) {
    return comment.empty() ? s : s + "  " + comment;
  };
  for (std::size_t r = 0; r < set.records.size(); r++) {
    const Record &rec = set.records[r];
    if (r > 0) out.push_back("");
    for (const std::string &c : rec.comments) out.push_back(c);
    const std::string name = NeedsQuotes(rec.name) ? "\"" + rec.name + "\"" : rec.name;
    out.push_back(withComment(rec.name.empty() ? rec.kind : rec.kind + " " + name, rec.trailingComment));
    // The values of a record start in one column: column 20, or one blank after the record's longest key.
    std::size_t width = 15;
    for (const Field &f : rec.fields)
      if (!f.tokens.empty()) width = std::max(width, f.key.size() + 1);
    for (const Field &f : rec.fields) {
      for (const std::string &c : f.comments) out.push_back("    " + c);
      std::string line = "    " + f.key;
      if (!f.tokens.empty()) {
        line += std::string(width - f.key.size(), ' ');
        line += Joined(f.tokens, 0);
      }
      out.push_back(withComment(line, f.trailingComment));
    }
  }
  if (!set.trailingComments.empty()) {
    if (!out.empty()) out.push_back("");
    for (const std::string &c : set.trailingComments) out.push_back(c);
  }
  return out;
}

}  // namespace io
}  // namespace azure2
