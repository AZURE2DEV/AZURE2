// The generic record reader and writer (RecordReader.h; docs/upgrade/project_format3.md sections 3.1-3.3):
//  1. tokens, quotes, comments, blank lines, indentation, lines and columns on a synthetic section;
//  2. each reader code (E_UNTERMINATED_QUOTE, E_QUOTE_IN_TOKEN, E_UNKNOWN_KIND, E_HEADER_NAME, E_KEY_SYNTAX) triggered
//     alone, and all problems of a section returned at once;
//  3. round trips: Read(Write(Read(x))) == Read(x) (records, tokens, comments) and Write(Read(Write(r))) == Write(r)
//     byte for byte, on the synthetic section and on every record section of the format-3 fixtures given as
//     arguments (experiment/tests/data/*.azr).
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "azure2/io/Document.h"
#include "azure2/io/RecordReader.h"

using namespace azure2::io;

namespace {

int failures = 0;
int checks = 0;

void check(const std::string &what, bool ok) {
  checks++;
  if (!ok) {
    std::printf("  FAIL  %s\n", what.c_str());
    failures++;
  }
}

// Record kinds, written out here: io/ knows none (the table is experiment/'s).  The synthetic sections use all of
// them; a fixture section is read with the kinds of its section, as the engine reads it (within a section kinds and
// field keys are disjoint, across sections they are not: a data set has a field `target`).
const std::vector<std::string> kKinds = {"dataset", "normalization", "energy_scale", "extrapolation", "target", "beam",
                                         "detector", "datum"};
const std::vector<std::pair<std::string, std::vector<std::string>>> kSectionKinds = {
    {"data", {"dataset", "normalization", "energy_scale"}},
    {"targets", {"target", "beam", "detector"}},
    {"extrapolations", {"extrapolation"}},
    {"modelData", {"datum"}}};

std::vector<std::string> Lines(const std::string &text) {
  std::vector<std::string> lines;
  std::istringstream in(text);
  std::string l;
  while (std::getline(in, l)) lines.push_back(l);
  return lines;
}

// Records equal in everything but line and column numbers.
bool SameRecords(const RecordSet &a, const RecordSet &b) {
  if (a.records.size() != b.records.size() || a.trailingComments != b.trailingComments) return false;
  for (std::size_t r = 0; r < a.records.size(); r++) {
    const Record &x = a.records[r], &y = b.records[r];
    if (x.kind != y.kind || x.name != y.name || x.comments != y.comments || x.trailingComment != y.trailingComment ||
        x.fields.size() != y.fields.size())
      return false;
    for (std::size_t f = 0; f < x.fields.size(); f++) {
      const Field &p = x.fields[f], &q = y.fields[f];
      if (p.key != q.key || p.comments != q.comments || p.trailingComment != q.trailingComment ||
          p.tokens.size() != q.tokens.size())
        return false;
      for (std::size_t t = 0; t < p.tokens.size(); t++)
        if (p.tokens[t].text != q.tokens[t].text || p.tokens[t].quoted != q.tokens[t].quoted) return false;
    }
  }
  return true;
}

void RoundTrip(const std::string &what, const std::vector<std::string> &lines,
               const std::vector<std::string> &kinds = kKinds) {
  const RecordReadResult first = ReadRecords(lines, 1, kinds);
  const std::vector<std::string> written = WriteRecords(first.set);
  const RecordReadResult second = ReadRecords(written, 1, kinds);
  check(what + ": Read(Write(Read(x))) == Read(x)", SameRecords(first.set, second.set));
  check(what + ": the canonical text reads without problems", second.diagnostics.size() == first.diagnostics.size());
  check(what + ": Write(Read(Write(r))) == Write(r)", WriteRecords(second.set) == written);
}

// Exactly one diagnostic, with this code, at this line and column.
void OneCode(const std::string &text, const std::string &code, int line, int col) {
  const RecordReadResult r = ReadRecords(Lines(text), 1, kKinds);
  check(code + ": exactly one diagnostic", r.diagnostics.size() == 1);
  if (r.diagnostics.size() == 1) {
    check(code + ": its code", r.diagnostics[0].code == code);
    check(code + ": its line and column", r.diagnostics[0].line == line && r.diagnostics[0].col == col);
    check(code + ": a message", !r.diagnostics[0].message.empty());
  }
}

std::string Slurp(const std::string &path) {
  std::ifstream in(path, std::ios::binary);
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

}  // namespace

int main(int argc, char **argv) {
  const std::string section =
      "# about the first record\n"
      "dataset s01   # header comment\n"
      "    reaction       15N(p,a)12C\n"
      "\t# about the file\n"
      "  file  \"data/run 1.dat\"   # trailing\n"
      "\n"
      "    columns E_lab[MeV]   _ value error\n"
      "    source \"a # inside quotes\" more\n"
      "    flag\n"
      "normalization D\n"
      "value 1 fit prior 1 8%\n"
      "# at the end\n";
  {
    const RecordReadResult r = ReadRecords(Lines(section), 10, kKinds);
    check("synthetic: no problems", r.diagnostics.empty());
    check("synthetic: two records", r.set.records.size() == 2);
    if (r.set.records.size() == 2) {
      const Record &d = r.set.records[0];
      check("header: kind, name, line, columns", d.kind == "dataset" && d.name == "s01" && d.line == 11 && d.col == 1 &&
                                                     d.nameCol == 9);
      check("header: comment above and after", d.comments == std::vector<std::string>{"# about the first record"} &&
                                                   d.trailingComment == "# header comment");
      check("five fields", d.fields.size() == 5);
      if (d.fields.size() == 5) {
        check("field line and key column (indentation not significant)",
              d.fields[1].key == "file" && d.fields[1].line == 14 && d.fields[1].col == 3);
        check("quoted token with a blank", d.fields[1].tokens.size() == 1 && d.fields[1].tokens[0].text == "data/run 1.dat" &&
                                               d.fields[1].tokens[0].quoted && d.fields[1].tokens[0].col == 9);
        check("comment above a field and after its value",
              d.fields[1].comments == std::vector<std::string>{"# about the file"} &&
                  d.fields[1].trailingComment == "# trailing");
        check("blank line ignored, tokens and their columns",
              d.fields[2].tokens.size() == 4 && d.fields[2].tokens[0].text == "E_lab[MeV]" &&
                  d.fields[2].tokens[0].col == 13 && d.fields[2].tokens[1].text == "_" && d.fields[2].tokens[1].col == 26);
        check("'#' inside quotes is text", d.fields[3].tokens.size() == 2 && d.fields[3].tokens[0].text == "a # inside quotes" &&
                                               d.fields[3].trailingComment.empty());
        check("a key without value tokens", d.fields[4].key == "flag" && d.fields[4].tokens.empty());
      }
      check("second record: fields up to the end", r.set.records[1].fields.size() == 1 &&
                                                       r.set.records[1].fields[0].tokens.size() == 5);
      check("comments after the last field", r.set.trailingComments == std::vector<std::string>{"# at the end"});
    }
    RoundTrip("synthetic", Lines(section));
    const std::vector<std::string> canonical = WriteRecords(r.set);
    check("canonical form", canonical.size() > 3 && canonical[1] == "dataset s01  # header comment" &&
                                canonical[2] == "    reaction       15N(p,a)12C" &&
                                canonical[4] == "    file           \"data/run 1.dat\"  # trailing");
  }

  {
    // The values of one record start in one column: one blank after a key of 15 characters or more.
    const RecordReadResult r = ReadRecords(Lines("potential p1\n    pair p+12C\n    matching_radius 15 fm\n"
                                                 "beam b\n    spread 1 keV lab\n"),
                                           1, {"potential", "beam"});
    const std::vector<std::string> w = WriteRecords(r.set);
    check("a long key moves the record's value column",
          w.size() == 6 && w[1] == "    pair            p+12C" && w[2] == "    matching_radius 15 fm");
    check("a record without a long key keeps column 20", w.size() == 6 && w[5] == "    spread         1 keV lab");
  }

  // Each code alone.
  OneCode("dataset a\n    file \"x y\n", "E_UNTERMINATED_QUOTE", 2, 10);
  OneCode("dataset a\n    file \"x\"y\n", "E_QUOTE_IN_TOKEN", 2, 13);
  OneCode("dataset a\n    file x\"y\n", "E_QUOTE_IN_TOKEN", 2, 11);
  OneCode("datset a\n", "E_UNKNOWN_KIND", 1, 1);
  OneCode("dataset\n", "E_HEADER_NAME", 1, 1);
  OneCode("dataset a b\n", "E_HEADER_NAME", 1, 11);
  OneCode("dataset a\n    1file x\n", "E_KEY_SYNTAX", 2, 5);
  {
    // Every problem, not the first; the record is still read around them.
    const RecordReadResult r =
        ReadRecords(Lines("junk 1\ndataset a b\n    file \"x\n    -k 2\n    unit b\ntarget\n"), 1, kKinds);
    check("all problems returned", r.diagnostics.size() == 5);
    check("records read around problems",
          r.set.records.size() == 2 && r.set.records[0].name == "a" && r.set.records[0].fields.size() == 2);
  }
  {
    // A section of a Document: line numbers are the file's.
    const Document doc = Document::Parse("<config>\nformat=3\n</config>\n<data>\nnormalization D\n    value 1 fixed\n</data>\n");
    const RecordReadResult r = ReadRecords(doc, *doc.Find("data"), kKinds);
    check("document section: file line numbers",
          r.set.records.size() == 1 && r.set.records[0].line == 5 && r.set.records[0].fields[0].line == 6);
  }

  // The fixtures: every record section of every file given.
  int sections = 0;
  for (int a = 1; a < argc; a++) {
    const Document doc = Document::Parse(Slurp(argv[a]));
    for (const Document::Section &s : doc.Sections())
      for (const auto &sk : kSectionKinds) {
        if (s.name != sk.first) continue;
        RoundTrip(std::string(argv[a]) + " <" + s.name + ">", doc.Body(s), sk.second);
        check(std::string(argv[a]) + " <" + s.name + ">: no syntax problems", ReadRecords(doc, s, sk.second).diagnostics.empty());
        sections++;
      }
  }
  check("fixture sections found when files are given", argc == 1 || sections > 0);

  std::printf("%d checks, %d failures (%d fixture sections)\n", checks, failures, sections);
  return failures == 0 ? 0 : 1;
}
