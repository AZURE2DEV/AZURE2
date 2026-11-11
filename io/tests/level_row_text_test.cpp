// The text of a <levels> row (LevelRowText.h): synthetic rows (named tokens, a replaced token right-aligned in its
// width or widening the line, the optional columns appended in their groups, refusals), and every row of every
// project file under the given directories: split into its columns, written back unedited and with each column set
// to its own token, byte for byte; a replaced width reads back (LevelRow) with every other column unchanged.
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "azure2/io/Document.h"
#include "azure2/io/LevelRowText.h"
#include "azure2/io/LevelRows.h"

using namespace azure2::io;
namespace fs = std::filesystem;

namespace {

int failures = 0;
long checks = 0;

void check(const std::string &what, bool ok) {
  checks++;
  if (!ok) {
    if (failures < 30) std::printf("  FAIL  %s\n", what.c_str());
    failures++;
  }
}

std::string Read(const fs::path &p) {
  std::ifstream in(p, std::ios::binary);
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

LevelRow RowOf(const std::string &line) {
  std::istringstream in(line);
  return LevelRow(in);
}

bool ReadsAsRow(const std::string &line, LevelRow *out = nullptr) {
  std::istringstream in(line);
  LevelRow row(in);
  if (in.rdstate() & (std::ios::failbit | std::ios::badbit)) return false;
  if (out) *out = row;
  return true;
}

const std::string kRow =
    " 1.5   -1       3.5032    0    1    2    1    2    5    1    0               -0.49    1    1  0.5   -1            0 "
    "      0 13.0057    0    7      1.94351            0    0    0          0.0    10       0            0            0 "
    "      5";

void Synthetic() {
  LevelRowText t;
  std::string why;
  check("the 31-column row splits", SplitLevelRow(kRow, t, why));
  check("no optional columns", !t.hasThm && !t.hasE1);
  check("levelE token", t.Token(LevelColumnIndex("levelE")) == "3.5032");
  check("gamma token", t.Token(LevelColumnIndex("gamma")) == "-0.49");
  check("ecMultMask token", t.Token(LevelColumnIndex("ecMultMask")) == "5");
  check("an absent optional column reads empty", t.Token(LevelColumnIndex("e1")).empty());
  check("34 column names", LevelColumnNames().size() == 34 && LevelColumnIndex("nonsense") == -1);
  check("unedited: the same bytes", EditLevelRow(t, {}) == kRow);

  // A shorter token right-aligned in the old width; a longer one widens the line.
  const int g = LevelColumnIndex("gamma");
  const std::string shorter = EditLevelRow(t, {{g, "2"}});
  check("shorter: same length", shorter.size() == kRow.size());
  check("shorter: right-aligned", shorter.find("                   2    1    1") != std::string::npos);
  const std::string longer = EditLevelRow(t, {{g, "-0.123456789"}});
  check("longer: widened by the difference", longer.size() == kRow.size() + 7);
  LevelRow r = RowOf(kRow), before = RowOf(kRow);
  check("longer reads back", ReadsAsRow(longer, &r) && std::fabs(r.gamma() + 0.123456789) < 1e-15 &&
                                 r.levelE() == before.levelE() && r.ecMultMask() == before.ecMultMask());

  // The optional columns: e1 alone after one blank; the THM pair after two blanks; both.
  const std::string e1 = EditLevelRow(t, {{LevelColumnIndex("e1"), "0.5"}});
  check("e1 appended", e1 == kRow + " 0.5");
  check("e1 reads", ReadsAsRow(e1, &r) && r.e1() == 0.5 && r.bindingE() == 0.0);
  const std::string thm = EditLevelRow(t, {{LevelColumnIndex("bindingE"), "2.2"}});
  check("THM pair appended, flag 0", thm == kRow + "  2.2  0");
  check("THM pair reads", ReadsAsRow(thm, &r) && r.bindingE() == 2.2 && r.gammaIsRWA() == 0 && r.e1() == 0.0);
  LevelRowText tt;
  check("THM row splits", SplitLevelRow(thm, tt, why) && tt.hasThm && !tt.hasE1);
  const std::string both = EditLevelRow(tt, {{LevelColumnIndex("e1"), "0.3"}, {g, "7"}});
  check("e1 after the THM pair", both.substr(both.size() - 13) == "  2.2  0  0.3");
  check("both read", ReadsAsRow(both, &r) && r.bindingE() == 2.2 && r.e1() == 0.3 && r.gamma() == 7.0);
  check("in place on a THM row", EditLevelRow(tt, {{LevelColumnIndex("gammaIsRWA"), "1"}}) == kRow + "  2.2  1");

  // Refusals: too few and too many tokens.
  check("30 tokens refused", !SplitLevelRow(kRow.substr(0, kRow.rfind(' ')), t, why) && !why.empty());
  check("35 tokens refused", !SplitLevelRow(kRow + " 1 2 3 4", t, why));
}

}  // namespace

int main(int argc, char **argv) {
  Synthetic();
  long files = 0, rows = 0;
  for (int i = 1; i < argc; i++)
    for (const auto &e : fs::recursive_directory_iterator(argv[i])) {
      if (!e.is_regular_file() || e.path().extension() != ".azr" || e.path().string().find("/output/") != std::string::npos)
        continue;
      const Document d = Document::Parse(Read(e.path()));
      const Document::Section *s = d.Find("levels");
      if (!s || s->close == Document::Section::npos) continue;
      files++;
      for (std::size_t k = s->open + 1; k < s->close; k++) {
        const std::string &line = d.Line(k);
        if (line.find_first_not_of(" \t") == std::string::npos || !ReadsAsRow(line)) continue;
        rows++;
        const std::string where = e.path().string() + ":" + std::to_string(k + 1);
        LevelRowText t;
        std::string why;
        if (!SplitLevelRow(line, t, why)) {
          check(where + ": splits (" + why + ")", false);
          continue;
        }
        check(where + ": unedited, the same bytes", EditLevelRow(t, {}) == line);
        bool same = true;
        for (int c = 0; c < (int)LevelColumnNames().size(); c++)
          if (!t.Token(c).empty()) same = same && EditLevelRow(t, {{c, t.Token(c)}}) == line;
        check(where + ": each column written with its own token, the same bytes", same);
        LevelRow before = RowOf(line), after = before;
        const std::string edited = EditLevelRow(t, {{LevelColumnIndex("gamma"), "123.5"}});
        check(where + ": a replaced width reads back, the rest unchanged",
              ReadsAsRow(edited, &after) && after.gamma() == 123.5 && after.levelE() == before.levelE() &&
                  after.ir() == before.ir() && after.l() == before.l() && after.s() == before.s() &&
                  after.chRad() == before.chRad() && after.e1() == before.e1() && after.bindingE() == before.bindingE() &&
                  after.gammaIsRWA() == before.gammaIsRWA() && after.ecMultMask() == before.ecMultMask());
      }
    }
  std::printf("  %ld files, %ld rows\n", files, rows);
  if (files == 0) check("no project file found", false);
  std::printf("%s: %ld checks, %d failure(s)\n", failures ? "FAIL" : "PASS", checks, failures);
  return failures ? 1 : 0;
}
