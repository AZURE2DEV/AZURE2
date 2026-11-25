// Forms of an old project that only an era before the 2020 engine wrote (OldEra.h): docs/upgrade/project_format3.md
// section 8 and the history of the format in tests/convert/fixtures/history/README.md.
#include "OldEra.h"

#include "azure2/convert/LineRepair.h"

#include <cstdlib>
#include <map>
#include <sstream>
#include <vector>

namespace azure2 {
namespace convert {

namespace {

std::vector<std::string> Words(const std::string &line) {
  std::istringstream in(line);
  std::vector<std::string> out;
  for (std::string w; in >> w;) out.push_back(w);
  return out;
}

bool WholeInteger(const std::string &t, long &v) {
  if (t.empty()) return false;
  char *end = nullptr;
  v = std::strtol(t.c_str(), &end, 10);
  return *end == '\0';
}

bool Number(const std::string &t) {
  if (t.empty()) return false;
  char *end = nullptr;
  std::strtod(t.c_str(), &end);
  return end != t.c_str() && *end == '\0';
}

/// The line with its word number `index` (0-based) replaced by `text`, every other byte kept.
std::string ReplaceWord(const std::string &line, std::size_t index, const std::string &text) {
  std::size_t at = 0, n = 0;
  while (true) {
    const std::size_t a = line.find_first_not_of(" \t\r", at);
    if (a == std::string::npos) return line;
    const std::size_t b = line.find_first_of(" \t\r", a);
    if (n == index) return line.substr(0, a) + text + (b == std::string::npos ? "" : line.substr(b));
    if (b == std::string::npos) return line;
    at = b;
    n++;
  }
}

}  // namespace

std::string FoldOldExternalCapture(io::Document &doc, std::vector<std::pair<std::string, std::string>> &notes) {
  // The pair key -> external-capture multipolarity mask of the <externalCapture> lines, as the 2020 setup program
  // imported them (PairsTab::parseOldECSection, 6ecaee22): "isActive exitPairKey minJ maxJ multMask" a line, the mask
  // set on the exit pair; minJ and maxJ were not used by that import.
  std::map<long, long> mask;
  const io::Document::Section *ec = doc.Find("externalCapture");
  if (ec) {
    if (ec->close == io::Document::Section::npos) return "<externalCapture> is not closed";
    for (std::size_t i = ec->open + 1; i < ec->close; i++) {
      const std::string &line = doc.Line(i);
      const std::vector<std::string> w = Words(line);
      if (w.empty()) continue;
      long active = 0, pair = 0, m = 0;
      const std::string where = "<externalCapture> line " + std::to_string(i + 1) + " (\"" + line + "\")";
      if (w.size() != 5 || !WholeInteger(w[0], active) || !WholeInteger(w[1], pair) || !Number(w[2]) || !Number(w[3]) ||
          !WholeInteger(w[4], m))
        return where + " is not \"isActive exitPairKey minJ maxJ multMask\", the form the 2020 setup program imported";
      if (active != 1)
        // The 2020 import applied the mask of an inactive line too; the engine of the era that wrote the section is not
        // in the history, so whether an inactive line counted is not known: no guess.
        return where + " is inactive: the 2020 setup program's import applied its mask anyway, the engine that wrote "
               "the section is not known to have done so.  Delete the line (no external capture to pair " +
               std::to_string(pair) + ") or set its first field to 1 in a copy of the project, and convert the copy";
      if (mask.count(pair) && mask[pair] != m)
        return where + " gives pair " + std::to_string(pair) + " a second mask (" + std::to_string(mask[pair]) + " and " +
               std::to_string(m) + ")";
      mask[pair] = m;
      notes.push_back({"<externalCapture>", "line " + std::to_string(i + 1) + ": external capture to pair " +
                                                std::to_string(pair) + " with the multipolarity mask " + std::to_string(m) +
                                                " (the column ecMultMask of that pair's <levels> rows), as the 2020 setup "
                                                "program imported it; minJ " + w[2] + " and maxJ " + w[3] +
                                                " were not used by that import and are not converted"});
    }
  }
  // <levels>: rows of 30 columns (before ecMultMask existed) get it; a pair named by <externalCapture> takes its mask
  const io::Document::Section *lv = doc.Find("levels");
  int padded = 0, set = 0;
  if (lv && lv->close != io::Document::Section::npos) {
    std::vector<std::string> body;
    for (std::size_t i = lv->open + 1; i < lv->close; i++) {
      std::string line = doc.Line(i);
      const std::vector<std::string> w = Words(line);
      long ir = 0;
      bool glued = false;
      for (const std::string &x : w) glued = glued || GluedNumbers(x);
      if (w.size() >= 30 && !glued && WholeInteger(w[5], ir)) {
        const auto m = mask.find(ir);
        if (w.size() == 30) {
          // the trailing CR of a CRLF line stays at the end
          const bool cr = !line.empty() && line.back() == '\r';
          if (cr) line.pop_back();
          line += " " + std::to_string(m == mask.end() ? 0L : m->second);
          if (cr) line += '\r';
          padded++;
        } else if (m != mask.end() && w[30] != std::to_string(m->second)) {
          long old = 0;
          if (WholeInteger(w[30], old) && old != 0 && old != m->second)
            return "<levels> line " + std::to_string(i + 1) + " gives pair " + std::to_string(ir) + " the mask " + w[30] +
                   " and <externalCapture> gives it " + std::to_string(m->second);
          line = ReplaceWord(line, 30, std::to_string(m->second));
          set++;
        }
      }
      body.push_back(line);
    }
    if (padded || set) doc.ReplaceLines(lv->open + 1, lv->close, body);
  }
  if (padded)
    notes.push_back({"<levels> row", std::to_string(padded) + " <levels> row(s) of 30 columns (an AZURE2 before 2020): the "
                                     "31st column ecMultMask written, 0 unless <externalCapture> named the pair, as the "
                                     "2020 setup program read such a row (LevelsTab::readNuclearFile)"});
  if (ec) {
    const io::Document::Section *again = doc.Find("externalCapture");
    doc.ReplaceLines(again->open, again->close + 1, {});
  }
  return "";
}

std::string CloseOldLastRun(std::string &text, std::vector<std::pair<std::string, std::string>> &notes) {
  // <lastRun> holds the 2020 setup program's runtime choices (AZURESetup::writeLastRun), which no engine read; one that
  // is not closed at the end of the file changes nothing a run computed.  Any other section left open stays refused.
  std::istringstream in(text);
  bool open = false, closed = false;
  for (std::string line; std::getline(in, line);) {
    std::string t = line;
    t.erase(0, t.find_first_not_of(" \t"));
    t.erase(t.find_last_not_of(" \t\r") + 1);
    if (t == "<lastRun>") open = true;
    if (t == "</lastRun>") closed = true;
  }
  if (open && !closed) {
    const bool crlf = text.find("\r\n") != std::string::npos;
    if (!text.empty() && text.back() != '\n') text += crlf ? "\r\n" : "\n";
    text += crlf ? "</lastRun>\r\n" : "</lastRun>\n";
    notes.push_back({"<lastRun>", "<lastRun> is not closed at the end of the file: closed (the setup program's runtime "
                                  "choices; no engine read them, and the converter drops the section)"});
  }
  return "";
}

std::string TargetLineOfFirstEra(const std::string &line, std::string &note) {
  // <targetInt> of 2025-01-26 .. 2025-01-30 (17281967 .. 1c01a007): after the q coefficients "isConv nConv c0 ..
  // c(n-1)", sigma(E_lab) = c0 + c1 E + c2 E^2 + ... (TargetEffect::CalculateSigma of 17281967).  From 1c01a007 on the
  // writer put the equation of sigma between the flag and the count, always quoted (TargetIntTab::writeFile): a count
  // where the quoted equation stands is that first form.  The equation that gives the same sigma is written out.
  const std::vector<std::string> w = Words(line);
  if (w.size() < 11) return line;
  long nParams = 0, nQ = 0;
  if (!WholeInteger(w[8], nParams) || nParams < 0) return line;
  const std::size_t isQ = 9 + (std::size_t)nParams;
  if (w.size() < isQ + 2 || !WholeInteger(w[isQ + 1], nQ) || nQ < 0) return line;
  const std::size_t isConv = isQ + 2 + (std::size_t)nQ;
  if (w.size() < isConv + 2) return line;
  long flag = 0, n = 0;
  if (!WholeInteger(w[isConv], flag) || (flag != 0 && flag != 1)) return line;
  const std::string &next = w[isConv + 1];
  if (next.empty() || next[0] == '"' || !WholeInteger(next, n) || n < 0) return line;
  if (flag == 0 && n == 0) return line;  // "0 0": no coefficients either way, and the later grammar reads it so
  if (w.size() < isConv + 2 + (std::size_t)n) return line;
  for (long k = 0; k < n; k++)
    if (!Number(w[isConv + 2 + (std::size_t)k])) return line;
  // the rest after the coefficients: this era wrote nothing more
  if (w.size() != isConv + 2 + (std::size_t)n) return line;
  std::string eq;
  for (long k = 0; k < n; k++) {
    if (k) eq += "+";
    eq += "a" + std::to_string(k);
    if (k == 1) eq += "*x";
    if (k > 1) eq += "*x^" + std::to_string(k);
  }
  std::string out = line;
  // insert the equation after the flag (the word number isConv)
  std::size_t at = 0, idx = 0;
  while (idx <= isConv) {
    const std::size_t a = out.find_first_not_of(" \t", at);
    const std::size_t b = out.find_first_of(" \t\r", a);
    at = b;
    idx++;
  }
  out.insert(at, " \"" + eq + "\"");
  note = "the convolution coefficients of 2025-01-26 .. 01-30 (no equation token): sigma(E_lab) = " +
         (eq.empty() ? std::string("the constant sigma of the line") : eq) + " with a_k the coefficients, as that engine computed it";
  return out;
}

}  // namespace convert
}  // namespace azure2
