// Variants of old <segmentsData>, <segmentsTest> and <targetInt> lines (docs/upgrade/qt_gui.md section 19.5): those
// whose meaning is unambiguous are converted with a note saying what was read; the others are refused with the line's
// number in the file, the fields expected and the readings proposed, none applied by the converter.  Each proposed line
// is converted once to show that it reads, and the record it gives is compared with the values written out here.
//   convert_line_repair_test <scratch dir>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "azure2/convert/Convert.h"

using namespace azure2::convert;

namespace {

int failures = 0;
std::string scratch;

void check(const std::string &what, bool ok) {
  std::printf("%s  %s\n", ok ? "  ok  " : "  FAIL", what.c_str());
  if (!ok) failures++;
}

// The pairs of convert_constructs (p + 15N, a + 12C, ...): columns of a <levels> line.
const std::string kLevels =
    "<levels>\n"
    "1 -1 12.4 0 1 1 2 0 1 1 1 1.0 0.5 1 0.5 -1 0.0 1.007825 15.000109 1 7 12.1 12.1 0 0 0.0 0 4.85 0 0 0\n"
    "1 -1 12.4 0 1 2 0 2 1 1 0 1.0 0.0 1 0.0 1 0.0 4.002603 12.0 2 6 12.1 7.16 0 0 0.0 0 5.43 0 0 0\n"
    "</levels>\n";

std::string Project(const std::string &data, const std::string &test = "", const std::string &targetInt = "") {
  return "<config>\nformat=2\n</config>\n" + kLevels + "<segmentsData>\n" + data + "</segmentsData>\n<segmentsTest>\n" +
         test + "</segmentsTest>\n<targetInt>\n" + targetInt + "</targetInt>\n";
}

ConvertResult Convert(const std::string &text) {
  ConvertOptions o;
  o.projectDir = scratch;
  return ConvertToFormat3(text, o);
}

bool HasNote(const ConvertResult &r, const std::string &part) {
  for (const ConvertNote &n : r.notes)
    if (n.text.find(part) != std::string::npos) return true;
  return false;
}

/// The value of `key` in the record whose header is `header` ("" when absent).
std::string Field(const std::string &text, const std::string &header, const std::string &key) {
  std::istringstream in(text);
  std::string line;
  bool inside = false;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (!inside) {
      inside = line == header;
      continue;
    }
    if (line.empty() || line[0] == '<') break;
    std::istringstream w(line);
    std::string k;
    w >> k;
    if (k != key) continue;
    std::string rest;
    std::getline(w, rest);
    if (!rest.empty() && rest.back() == '\r') rest.pop_back();
    return rest.substr(rest.find_first_not_of(' '));
  }
  return "";
}

/// The 1-based line of `text` that equals `line` (0 when none).
int LineOf(const std::string &text, const std::string &line) {
  std::istringstream in(text);
  std::string l;
  for (int n = 1; std::getline(in, l); n++) {
    if (!l.empty() && l.back() == '\r') l.pop_back();
    if (l == line) return n;
  }
  return 0;
}

std::string Replace(std::string text, const std::string &from, const std::string &to) {
  const std::size_t at = text.find(from);
  if (at != std::string::npos) text.replace(at, from.size(), to);
  return text;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::printf("usage: convert_line_repair_test <scratch dir>\n");
    return 2;
  }
  scratch = argv[1];
  std::ofstream(scratch + "/my data.dat") << "0.3 1e-3 1e-4\n";
  const std::string good = "1 1 2 0.2 1.4 0 180 0 1 1 5 0.001 0.003 1 a01.dat\n";
  const ConvertResult base = Convert(Project(good));
  check("the plain line converts", base.ok);
  const std::string norm = Field(base.text, "dataset a01", "norm"), shift = Field(base.text, "dataset a01", "energy_shift");

  // ---- accepted, with a note
  {  // a UTF-8 byte-order mark before the first tag
    const ConvertResult r = Convert("\xEF\xBB\xBF" + Project(good));
    check("byte-order mark: converted", r.ok && HasNote(r, "byte-order mark"));
    check("byte-order mark: the same record", Field(r.text, "dataset a01", "norm") == norm);
  }
  {  // Windows line ends and tabs (the old reader's >> reads both)
    std::string crlf;
    for (char c : Project(good)) crlf += c == '\n' ? std::string("\r\n") : std::string(1, c);
    const ConvertResult r = Convert(crlf);
    check("CRLF line ends: converted, the same norm and shift",
          r.ok && Field(r.text, "dataset a01", "norm") == norm && Field(r.text, "dataset a01", "energy_shift") == shift);
    std::string tabs = good;
    for (char &c : tabs)
      if (c == ' ') c = '\t';
    const ConvertResult t = Convert(Project(tabs));
    check("tabs between the fields: converted, the same norm and shift",
          t.ok && Field(t.text, "dataset a01", "norm") == norm && Field(t.text, "dataset a01", "energy_shift") == shift);
  }
  {  // a comment line in each of the three sections
    const ConvertResult r = Convert(Project("# the data of A\n" + good, "# none\n", "# none either\n"));
    check("comment lines: converted, one note each", r.ok && HasNote(r, "a comment in <segmentsData>") &&
                                                         HasNote(r, "a comment in <segmentsTest>") &&
                                                         HasNote(r, "a comment in <targetInt>"));
    check("comment lines: the data line keeps its record", Field(r.text, "dataset a01", "norm") == norm);
  }
  {  // a data file whose name holds a blank, unquoted, and the file exists
    const ConvertResult r = Convert(Project("1 1 2 0.2 1.4 0 180 0 1 0 0 my data.dat\n"));
    check("file name with a blank, the file there: converted with a note",
          r.ok && HasNote(r, "the data file's name holds a blank"));
    check("file name with a blank: the record reads that file", r.text.find("my data.dat") != std::string::npos);
    // the 2020 and 2025-08 engines read the data file as the rest of the line: a name with blanks at the end of the
    // line is that name, file there or not
    const ConvertResult rest = Convert(Project("1 1 2 0.2 1.4 0 180 0 1 0 0 no such.dat\n"));
    check("file name with a blank to the end of the line, no such file: read as that name (the 2020 rule)",
          rest.ok && rest.text.find("\"no such.dat\"") != std::string::npos && HasNote(rest, "the rest of the line"));
    // followed by the later flags, with no such file: no era read it as one name: refused, the quoted name proposed
    const ConvertResult missing = Convert(Project("1 1 2 0.2 1.4 0 180 0 1 0 0 0 0 0 no such.dat 0 0\n"));
    check("file name with a blank, later flags after it, no such file: refused (no guess)", !missing.ok && missing.lineProblem.Found());
    bool quoted = false;
    for (const LineReading &g : missing.lineProblem.readings) quoted = quoted || g.line.find("\"no such.dat\"") != std::string::npos;
    check("file name with a blank, no such file: the quoted name proposed", quoted);
  }
  {  // a <targetInt> line: blanks inside its quoted segment list
    const std::string t = "1 \"1, 2\" 10 0 0 0 0 \"\" 0 1 5 1 0.943 0.834 0.688 0.517\n";
    const ConvertResult r = Convert(Project(good + "1 1 2 0.2 1.4 0 180 0 1 0 0 a02.dat\n", "", t));
    check("quoted segment list with a blank: converted with a note", r.ok && HasNote(r, "blanks inside a quoted token"));
    const ConvertResult same = Convert(Project(good + "1 1 2 0.2 1.4 0 180 0 1 0 0 a02.dat\n", "", Replace(t, "\"1, 2\"", "\"1,2\"")));
    check("quoted segment list with a blank: the same project as without it", same.ok && same.text == r.text);
  }
  {  // a <targetInt> line ending after its empty stopping equation: the parameter count can only be 0
    const ConvertResult r = Convert(Project(good, "", "1 \"1\" 10 0 0 0 0 \"\"\n"));
    const ConvertResult nine = Convert(Project(good, "", "1 \"1\" 10 0 0 0 0 \"\" 0\n"));
    check("target line without its parameter count after \"\": converted with a note", r.ok && HasNote(r, "parameter count"));
    check("target line without its parameter count: the same project as with 0", nine.ok && nine.text == r.text);
  }

  // ---- refused, with the line, the fields expected and the readings
  {  // no data file: refused without a crash (it ended the program with an uncaught exception before)
    const std::string line = "1 1 2 0.2 1.4 0 180 0 1 0 0 0 0 0";
    const std::string text = Project(good + line + "\n");
    const ConvertResult r = Convert(text);
    check("data line without a file: refused, 'names no data file'",
          !r.ok && r.refusal.find("names no data file") != std::string::npos);
    check("data line without a file: its line number in the file", r.lineProblem.fileLine == LineOf(text, line));
    check("data line without a file: no reading proposed", r.lineProblem.readings.empty());
  }
  {  // an energy shift and its error without the vary flag: two readings, neither applied
    const std::string line = "1 1 2 0.2 1.4 0 180 0 1 1 5 0.001 0.003 a01.dat";
    const std::string text = Project(line + "\n");
    const ConvertResult r = Convert(text);
    check("shift without its vary flag: refused", !r.ok);
    check("shift without its vary flag: the line, its number, the fields expected",
          r.lineProblem.section == "segmentsData" && r.lineProblem.fileLine == LineOf(text, line) &&
              r.lineProblem.text == line && r.lineProblem.expected.find("varyShift") != std::string::npos);
    check("shift without its vary flag: two readings", r.lineProblem.readings.size() == 2);
    check("shift without its vary flag: the refusal carries the diagnosis",
          r.refusal.find("line " + std::to_string(r.lineProblem.fileLine) + " of the file") != std::string::npos &&
              r.refusal.find("possible readings") != std::string::npos);
    if (r.lineProblem.readings.size() == 2) {
      // reading 1: the shift fixed -> the record of the full line with vary 0
      const ConvertResult a = Convert(Replace(text, line, r.lineProblem.readings[0].line));
      const ConvertResult want = Convert(Project("1 1 2 0.2 1.4 0 180 0 1 1 5 0.001 0.003 0 a01.dat\n"));
      check("reading 1 (shift fixed): converts to the line written with vary 0", a.ok && want.ok && a.text == want.text);
      // reading 2: the norm's error missing -> norm 1 fit, error 0, shift 5 with error 0.001 and vary 0.003 read as 0
      const ConvertResult b = Convert(Replace(text, line, r.lineProblem.readings[1].line));
      check("reading 2 (norm error 0): converts", b.ok && Field(b.text, "dataset a01", "norm").compare(0, 5, "1 fit") == 0);
    }
  }
  {  // a decimal comma
    const std::string line = "1 1 2 0.2 1.4 0 180 0 1 1 5 0,001 0.003 1 a01.dat";
    const ConvertResult r = Convert(Project(line + "\n"));
    check("decimal comma: refused with the reading 0,001 -> 0.001",
          !r.ok && !r.lineProblem.readings.empty() && r.lineProblem.readings[0].line.find("0.001 0.003 1 a01.dat") != std::string::npos);
    const ConvertResult a = Convert(Project(r.lineProblem.readings.empty() ? line : r.lineProblem.readings[0].line + "\n"));
    check("decimal comma: the reading converts to the plain line's record", a.ok && a.text == base.text);
  }
  {  // a test line without its data type
    const std::string line = "1 1 2 0.2 1.4 0.01 0 0 0";
    const std::string text = Project(good, line + "\n");
    const ConvertResult r = Convert(text);
    check("test line without isDiff: refused, its line number",
          !r.ok && r.lineProblem.section == "segmentsTest" && r.lineProblem.fileLine == LineOf(text, line));
    check("test line without isDiff: readings 0 and 1",
          r.lineProblem.readings.size() == 2 && r.lineProblem.readings[0].line == line + " 0" &&
              r.lineProblem.readings[1].line == line + " 1");
  }
  {  // the line number is the file's, also behind a byte-order mark and with CRLF line ends
    const std::string line = "1 1 2 0.2 1.4 0 180 0 1 0 0 0 0 0";
    std::string crlf;
    for (char c : Project(good + line + "\n")) crlf += c == '\n' ? std::string("\r\n") : std::string(1, c);
    const ConvertResult r = Convert("\xEF\xBB\xBF" + crlf);
    check("BOM and CRLF: the refused line's number and text (no CR)",
          !r.ok && r.lineProblem.fileLine == LineOf(crlf, line) && r.lineProblem.text == line);
  }

  std::printf("%s\n", failures ? "FAILED" : "all passed");
  return failures ? 1 : 0;
}
