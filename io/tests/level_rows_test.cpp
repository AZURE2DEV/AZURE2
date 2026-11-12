// io::LevelRow and io::ReadLevelSection (LevelRows.h), the grammar of the <levels> section: (1) synthetic rows and
// sections against hand values -- the optional columns told apart by their number, the halving of 2s and 2l, the
// failures, blank lines, the end tag; (2) every <levels> section of the project files under the directories given
// (tests/, examples/) against the engine's reader before it moved to io/ (NucLine and the loop of CNuc::Fill at
// a16f2a69, copied verbatim below): every field of every row bit for bit, and the same malformed / complete outcome.
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "azure2/io/Document.h"
#include "azure2/io/LevelRows.h"

namespace fs = std::filesystem;
using azure2::io::LevelRow;
using azure2::io::LevelSection;

namespace {
int failures = 0, checks = 0;
void check(const std::string &what, bool ok) {
  checks++;
  if (!ok) {
    std::printf("  FAIL  %s\n", what.c_str());
    failures++;
  }
}

// ---- the reader before the move (include/NucLine.h at a16f2a69, verbatim but for the class name, its comments
//      and the semicolons after the member functions) ----
class OldNucLine {
 public:
  OldNucLine(std::istream &stream) {
    stream >> levelJ_ >> levelPi_ >> levelE_ >> levelFix_ >> aa_ >> ir_ >> s_ >> l_ >> levelID_ >> isActive_ >> channelFix_ >> gamma_ >> j1_ >> pi1_ >> j2_ >> pi2_ >> e2_ >> m1_ >> m2_ >> z1_ >> z2_ >> entranceSepE_ >> sepE_ >> j3_ >> pi3_ >> e3_ >> pType_ >> chRad_ >> g1_ >> g2_ >> ecMultMask_;
    e1_ = 0.;
    bindingE_ = 0.0;
    gammaIsRWA_ = 0;
    if (!stream.fail() && !stream.eof()) {
      double v[3];
      int n = 0;
      stream >> std::ws;
      while (!stream.eof() && n < 4) {
        double x;
        if (!(stream >> x)) break;
        if (n < 3) v[n] = x;
        n++;
        if (!stream.eof()) stream >> std::ws;
      }
      if (!stream.fail()) {
        if (n == 1) {
          e1_ = v[0];
        } else if (n == 2 || n == 3) {
          bindingE_ = v[0];
          gammaIsRWA_ = (int)v[1];
          if (n == 3) e1_ = v[2];
        } else if (n > 3) {
          stream.setstate(std::ios::failbit);
        }
      }
    }
    s_ /= 2.;
    l_ /= 2;
  }
  double levelJ() const { return levelJ_; }
  int levelPi() const { return levelPi_; }
  double levelE() const { return levelE_; }
  int levelFix() const { return levelFix_; }
  int aa() const { return aa_; }
  int ir() const { return ir_; }
  double s() const { return s_; }
  int l() const { return l_; }
  int levelID() const { return levelID_; }
  int isActive() const { return isActive_; }
  int channelFix() const { return channelFix_; }
  double gamma() const { return gamma_; }
  double j1() const { return j1_; }
  int pi1() const { return pi1_; }
  double j2() const { return j2_; }
  int pi2() const { return pi2_; }
  double e2() const { return e2_; }
  double m1() const { return m1_; }
  double m2() const { return m2_; }
  int z1() const { return z1_; }
  int z2() const { return z2_; }
  double entranceSepE() const { return entranceSepE_; }
  double sepE() const { return sepE_; }
  int j3() const { return j3_; }
  int pi3() const { return pi3_; }
  double e3() const { return e3_; }
  int pType() const { return pType_; }
  double chRad() const { return chRad_; }
  double g1() const { return g1_; }
  double g2() const { return g2_; }
  double ecMultMask() const { return ecMultMask_; }
  double e1() const { return e1_; }
  double bindingE() const { return bindingE_; }
  int gammaIsRWA() const { return gammaIsRWA_; }

 private:
  double levelJ_;
  int levelPi_;
  double levelE_;
  int levelFix_;
  int aa_;
  int ir_;
  double s_;
  int l_;
  int levelID_;
  int isActive_;
  int channelFix_;
  double gamma_;
  double j1_;
  int pi1_;
  double j2_;
  int pi2_;
  double e2_;
  double m1_;
  double m2_;
  int z1_;
  int z2_;
  double entranceSepE_;
  double sepE_;
  int j3_;
  int pi3_;
  double e3_;
  int pType_;
  double chRad_;
  double g1_;
  double g2_;
  unsigned int ecMultMask_;
  double bindingE_;
  int gammaIsRWA_;
  double e1_;
};

/// The loop of CNuc::Fill at a16f2a69 over a <levels> section, keeping its rows and how it ended.
struct OldSection {
  std::vector<OldNucLine> rows;
  bool malformed = false, complete = false;
};
OldSection OldRead(std::istream &in) {
  OldSection out;
  std::string line = "";
  getline(in, line);  // <levels>
  line = "";
  while (!in.eof() && line != "</levels>") {
    getline(in, line);
    bool empty = true;
    for (unsigned int i = 0; i < line.size(); ++i)
      if (line[i] != ' ' && line[i] != '\t') {
        empty = false;
        break;
      }
    if (empty == true) continue;
    if (!in.eof() && line != "</levels>") {
      std::istringstream stm;
      stm.str(line);
      OldNucLine Line(stm);
      if (stm.rdstate() & (std::stringstream::failbit | std::stringstream::badbit)) {
        out.malformed = true;
        return out;
      }
      out.rows.push_back(Line);
    }
  }
  out.complete = (line == "</levels>");
  return out;
}

bool Same(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

/// Every field of the two rows, bit for bit.
bool SameRow(const OldNucLine &o, const LevelRow &n) {
  return Same(o.levelJ(), n.levelJ()) && o.levelPi() == n.levelPi() && Same(o.levelE(), n.levelE()) &&
         o.levelFix() == n.levelFix() && o.aa() == n.aa() && o.ir() == n.ir() && Same(o.s(), n.s()) && o.l() == n.l() &&
         o.levelID() == n.levelID() && o.isActive() == n.isActive() && o.channelFix() == n.channelFix() &&
         Same(o.gamma(), n.gamma()) && Same(o.j1(), n.j1()) && o.pi1() == n.pi1() && Same(o.j2(), n.j2()) &&
         o.pi2() == n.pi2() && Same(o.e2(), n.e2()) && Same(o.m1(), n.m1()) && Same(o.m2(), n.m2()) &&
         o.z1() == n.z1() && o.z2() == n.z2() && Same(o.entranceSepE(), n.entranceSepE()) && Same(o.sepE(), n.sepE()) &&
         o.j3() == n.j3() && o.pi3() == n.pi3() && Same(o.e3(), n.e3()) && o.pType() == n.pType() &&
         Same(o.chRad(), n.chRad()) && Same(o.g1(), n.g1()) && Same(o.g2(), n.g2()) &&
         Same(o.ecMultMask(), n.ecMultMask()) && Same(o.e1(), n.e1()) && Same(o.bindingE(), n.bindingE()) &&
         o.gammaIsRWA() == n.gammaIsRWA();
}

/// The <levels> section as Config::SectionStream gives it (trimmed tags), or "" without one.
std::string SectionText(const std::string &bytes) {
  azure2::io::Document doc = azure2::io::Document::Parse(bytes);
  const azure2::io::Document::Section *s = doc.Find("levels");
  if (!s) return "";
  std::string text;
  const std::size_t end = (s->close == azure2::io::Document::Section::npos) ? doc.NumLines() : s->close + 1;
  for (std::size_t i = s->open; i < end; i++) {
    bool closing = false;
    const std::string tag = azure2::io::TagOf(doc.Line(i), closing);
    text += (i == s->open || (i == s->close && !tag.empty())) ? (closing ? "</" : "<") + tag + ">" : doc.Line(i);
    text += "\n";
  }
  return text;
}

LevelSection Read(const std::string &text) {
  std::istringstream in(text);
  return azure2::io::ReadLevelSection(in);
}

const char *kRow = "1.5 -1 10.5 0 1 2 3 3 7 1 0 125.5 0.5 1 0 1 0.1 1.007825 12.0 1 6 0 1.94 0 1 0 0 4.5 5.586 0 3";
}  // namespace

int main(int argc, char **argv) {
  // (1) Rows.
  {
    std::istringstream in(kRow);
    LevelRow r(in);
    check("31 fields: no failbit", !in.fail());
    check("31 fields: values", r.levelJ() == 1.5 && r.levelPi() == -1 && r.levelE() == 10.5 && r.aa() == 1 &&
                                   r.ir() == 2 && r.levelID() == 7 && r.gamma() == 125.5 && r.m1() == 1.007825 &&
                                   r.z2() == 6 && r.sepE() == 1.94 && r.chRad() == 4.5 && r.g1() == 5.586 &&
                                   r.ecMultMask() == 3 && r.e2() == 0.1);
    check("2s = 3 -> s = 1.5; 2l = 3 -> l = 1 (integer division)", r.s() == 1.5 && r.l() == 1);
    check("31 fields: no optional values", r.e1() == 0.0 && r.bindingE() == 0.0 && r.gammaIsRWA() == 0);
  }
  {
    std::istringstream in(std::string(kRow) + " 0.25");
    LevelRow r(in);
    check("32 fields: e1", !in.fail() && r.e1() == 0.25 && r.bindingE() == 0.0 && r.gammaIsRWA() == 0);
  }
  {
    std::istringstream in(std::string(kRow) + " 2.224 1");
    LevelRow r(in);
    check("33 fields: binding energy and flag", !in.fail() && r.bindingE() == 2.224 && r.gammaIsRWA() == 1 && r.e1() == 0.0);
  }
  {
    std::istringstream in(std::string(kRow) + " 2.224 1 0.5  ");
    LevelRow r(in);
    check("34 fields (and trailing blanks): binding energy, flag, e1",
          !in.fail() && r.bindingE() == 2.224 && r.gammaIsRWA() == 1 && r.e1() == 0.5);
  }
  {
    std::istringstream in(std::string(kRow) + " 1 2 3 4");
    LevelRow r(in);
    check("35 fields: failbit", in.fail());
  }
  {
    std::istringstream in(std::string(kRow) + " x");
    LevelRow r(in);
    check("an optional field that is not a number: failbit", in.fail());
  }
  {
    std::istringstream in("1.5 -1 10.5 0 1 2 3 3 7 1 0 125.5");
    LevelRow r(in);
    check("12 fields: failbit", in.fail());
  }
  // (1) Sections.
  {
    const LevelSection s = Read(std::string("<levels>\n") + kRow + "\n  \t\n\n" + kRow + " 0.25\n</levels>\n");
    check("section: blank lines skipped, two rows, complete",
          s.rows.size() == 2 && !s.malformed && s.complete && s.rows[1].e1() == 0.25);
  }
  {
    const LevelSection s = Read(std::string("<levels>\n") + kRow + "\nnot a row\n" + kRow + "\n</levels>\n");
    check("section: a malformed line stops the reading, the rows before it kept",
          s.rows.size() == 1 && s.malformed && s.badLine == "not a row" && !s.complete);
  }
  {
    const LevelSection s = Read(std::string("<levels>\n") + kRow + "\n" + kRow);
    check("section: no </levels> -- incomplete, a last line without its newline not read",
          s.rows.size() == 1 && !s.malformed && !s.complete);
  }

  // (2) The project files: the new reader against the old one.
  long files = 0, rows = 0;
  for (int i = 1; i < argc; i++)
    for (const auto &e : fs::recursive_directory_iterator(argv[i])) {
      if (!e.is_regular_file() || e.path().extension() != ".azr" || e.path().string().find("/output/") != std::string::npos)
        continue;
      std::ifstream f(e.path(), std::ios::binary);
      std::ostringstream bytes;
      bytes << f.rdbuf();
      const std::string text = SectionText(bytes.str());
      if (text.empty()) continue;
      std::istringstream oldIn(text);
      const OldSection o = OldRead(oldIn);
      const LevelSection n = Read(text);
      bool same = o.rows.size() == n.rows.size() && o.malformed == n.malformed && o.complete == n.complete;
      for (std::size_t k = 0; same && k < o.rows.size(); k++) same = SameRow(o.rows[k], n.rows[k]);
      check(e.path().string() + ": rows as the reader before the move", same);
      files++;
      rows += (long)n.rows.size();
    }
  if (argc > 1) check("project files found", files > 100);
  std::printf("%s: %d checks (%ld project files, %ld rows), %d failure(s)\n", failures ? "FAIL" : "PASS", checks, files,
              rows, failures);
  return failures ? 1 : 0;
}
