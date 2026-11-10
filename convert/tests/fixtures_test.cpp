// The converter on real projects (docs/upgrade/project_format3.md section 8, evidence of step 3a):
//  1. the conversions of a fixed set of test projects equal the committed files of convert/tests/fixtures/ byte for
//     byte (and their side files), so a change of the converter's output is seen in review;
//  2. converting twice gives the same bytes; a format-3 input is refused ("already format 3"); every fixture
//     validates with zero errors;
//  3. T2's conversion has the records of the hand-written T2 of section 3.11 (experiment/tests/data/T2_format3.azr),
//     field by field (comments aside);
//  4. a project with questions is not converted without answers; with allowQuestions it is, with '?' placeholders.
//   convert_fixtures_test <repository root>
// With AZURE2_WRITE_CONVERT_FIXTURES=1 the fixtures are written instead of compared (after a deliberate change of
// the converter's output; the diff is then reviewed).
#include <algorithm>
#include <cstdio>
#include <iterator>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "azure2/convert/Convert.h"
#include "azure2/experiment/Schema.h"
#include "azure2/io/Document.h"
#include "azure2/io/RecordReader.h"

using namespace azure2::convert;

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

bool Slurp(const std::string &path, std::string &out) {
  std::ifstream in(path, std::ios::binary);
  if (!in) return false;
  std::ostringstream s;
  s << in.rdbuf();
  out = s.str();
  return true;
}

std::string Dir(const std::string &path) {
  const std::size_t slash = path.find_last_of('/');
  return slash == std::string::npos ? "." : path.substr(0, slash);
}

/// The records of the format-3 sections of a text, one line per field ("kind name: key tokens"), comments left out.
std::vector<std::string> Records(const std::string &text) {
  std::vector<std::string> out;
  const azure2::io::Document doc = azure2::io::Document::Parse(text);
  for (const azure2::io::Document::Section &s : doc.Sections()) {
    const std::vector<std::string> kinds = azure2::experiment::Schema::KindsOfSection(s.name);
    if (s.depth != 0 || kinds.empty()) continue;
    for (const azure2::io::Record &r : azure2::io::ReadRecords(doc, s, kinds).set.records)
      for (const azure2::io::Field &f : r.fields) {
        std::string line = r.kind + " " + r.name + ": " + f.key;
        for (const azure2::io::Token &t : f.tokens) line += " " + t.text;
        out.push_back(line);
      }
  }
  return out;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::printf("usage: convert_fixtures_test <repository root>\n");
    return 2;
  }
  const std::string root = argv[1];
  const std::string fixtures = root + "/convert/tests/fixtures/";
  const bool write = std::getenv("AZURE2_WRITE_CONVERT_FIXTURES") != nullptr;
  const std::vector<std::pair<std::string, std::string>> set = {
      {"tests/upgrade/T1/T1.azr", "T1"},
      {"tests/upgrade/T2/T2.azr", "T2"},
      {"tests/upgrade/T3/T3.azr", "T3"},
      {"tests/13N/13N.azr", "13N"},
      {"tests/natural_neon/natural_neon.azr", "natural_neon"},
      {"tests/18O_p_a_thm/18O_p_a_thm.azr", "18O_p_a_thm"},
      {"tests/beam_profile_kernel/beam_profile_kernel.azr", "beam_profile_kernel"},
      {"tests/norm_groups/norm_groups.azr", "norm_groups"},
      {"tests/shift_groups/shift_groups.azr", "shift_groups"},
      {"tests/pyazr/data/shift_fn_groups/shift_fn_groups.azr", "shift_fn_groups"},
      {"tests/target_density/fit_group.azr", "target_density_fit_group"},
  };
  std::string t2;
  for (const auto &item : set) {
    std::string old;
    check(item.second + ": read " + item.first, Slurp(root + "/" + item.first, old));
    ConvertOptions o;
    o.projectDir = Dir(root + "/" + item.first);
    const ConvertResult r = ConvertToFormat3(old, o);
    check(item.second + ": converted", r.ok);
    if (!r.ok) {
      std::printf("        %s\n", r.refusal.c_str());
      for (const auto &d : r.diagnostics) std::printf("        %s\n", d.ToText().c_str());
      continue;
    }
    if (item.second == "T2") t2 = r.text;
    // 1. the committed files
    std::vector<std::pair<std::string, std::string>> files = {{item.second + ".f3.azr", r.text}};
    for (const auto &side : r.sideFiles) files.push_back({item.second + "." + side.first, side.second});
    for (const auto &f : files) {
      if (write) {
        std::ofstream(fixtures + f.first, std::ios::binary) << f.second;
        continue;
      }
      std::string committed;
      check(item.second + ": " + f.first + " committed", Slurp(fixtures + f.first, committed));
      check(item.second + ": " + f.first + " byte for byte", committed == f.second);
    }
    // 2. idempotence, the format-3 refusal, validation
    check(item.second + ": converting twice gives the same bytes", ConvertToFormat3(old, o).text == r.text);
    const ConvertResult again = ConvertToFormat3(r.text, o);
    check(item.second + ": a format-3 input is refused", !again.ok && again.refusal == "already format 3");
    const Format3Check c = CheckFormat3Text(r.text, item.second);
    check(item.second + ": validates with zero errors", c.Ok());
  }
  // 3. T2 against the note's hand-written conversion
  {
    std::string hand;
    check("T2: read the hand-written T2_format3.azr", Slurp(root + "/experiment/tests/data/T2_format3.azr", hand));
    // the same fields of the same records, in any order of sections and records (fields keep their order per record)
    std::vector<std::string> a = Records(t2), b = Records(hand);
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    check("T2: the records of the hand-written conversion, field by field", a == b);
    std::vector<std::string> onlyA, onlyB;
    std::set_difference(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(onlyA));
    std::set_difference(b.begin(), b.end(), a.begin(), a.end(), std::back_inserter(onlyB));
    for (const std::string &x : onlyA) std::printf("        converter only: %s\n", x.c_str());
    for (const std::string &x : onlyB) std::printf("        hand only:      %s\n", x.c_str());
  }
  // 4. questions: none answered -> no result; allowQuestions -> '?' placeholders
  {
    std::string old;
    const std::string path = root + "/examples/c12c12_tumino2018/c12c12_tumino2018.azr";
    check("questions: read", Slurp(path, old));
    ConvertOptions o;
    o.projectDir = Dir(path);
    const ConvertResult r = ConvertToFormat3(old, o);
    check("questions: no result without answers", !r.ok && !r.questions.empty() && r.text.empty());
    o.allowQuestions = true;
    const ConvertResult q = ConvertToFormat3(old, o);
    check("questions: allowQuestions writes '?'", q.ok && q.text.find("    reaction       ?\n") != std::string::npos);
    bool onlyPlaceholders = !q.diagnostics.empty();
    for (const auto &d : q.diagnostics) onlyPlaceholders = onlyPlaceholders && d.code == std::string("E_PLACEHOLDER");
    check("questions: the placeholders are the only problems", onlyPlaceholders);
  }
  std::printf("convert_fixtures: %d checks, %d failures%s\n", checks, failures, write ? " (fixtures written)" : "");
  return failures ? 1 : 0;
}
