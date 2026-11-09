// Byte-identical round trips of real project files (docs/upgrade/project.md section 3, P1): every .azr under the
// given directories (outside output/), and variants made from each -- all line endings CRLF, all CR, no final
// newline, blank lines and trailing spaces and a comment added around every section, each top-level section removed
// in turn -- and, for every variant, Write(Parse(bytes)) == bytes; a removed section leaves every other byte in
// place (the bytes before and after it are those of the input).
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "azure2/io/Document.h"

using azure2::io::Document;
namespace fs = std::filesystem;

namespace {

int failures = 0;
long checks = 0;

void check(const std::string &what, bool ok) {
  checks++;
  if (!ok) {
    if (failures < 20) std::printf("  FAIL  %s\n", what.c_str());
    failures++;
  }
}

std::string Read(const fs::path &p) {
  std::ifstream in(p, std::ios::binary);
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

std::string ReplaceAll(std::string s, const std::string &a, const std::string &b) {
  std::string out;
  std::size_t i = 0;
  for (std::size_t j; (j = s.find(a, i)) != std::string::npos; i = j + a.size()) out += s.substr(i, j - i) + b;
  return out + s.substr(i);
}

}  // namespace

int main(int argc, char **argv) {
  std::vector<fs::path> files;
  for (int i = 1; i < argc; i++)
    for (const auto &e : fs::recursive_directory_iterator(argv[i]))
      if (e.is_regular_file() && e.path().extension() == ".azr" &&
          e.path().string().find("/output/") == std::string::npos)
        files.push_back(e.path());
  long variants = 0, sections = 0;
  for (const fs::path &f : files) {
    const std::string bytes = Read(f);
    const std::string name = f.string();
    std::vector<std::pair<std::string, std::string>> v = {
        {"as is", bytes},
        {"CRLF", ReplaceAll(ReplaceAll(bytes, "\r\n", "\n"), "\n", "\r\n")},
        {"CR", ReplaceAll(ReplaceAll(bytes, "\r\n", "\n"), "\n", "\r")},
        {"no final newline", (!bytes.empty() && bytes.back() == '\n') ? bytes.substr(0, bytes.size() - 1) : bytes},
    };
    // Blank lines, trailing spaces and a comment around every tag line.
    {
      Document d = Document::Parse(bytes);
      std::string out;
      for (std::size_t i = 0; i < d.NumLines(); i++) {
        bool closing = false;
        const bool tag = !azure2::io::TagOf(d.Line(i), closing).empty();
        if (tag && !closing) out += "\n# comment before a section\n";
        out += d.Line(i) + (tag ? "   \t" : "") + d.Ending(i);
        if (tag && closing) out += "\n\n";
      }
      v.push_back({"padded", out});
    }
    for (const auto &kv : v) {
      variants++;
      const Document d = Document::Parse(kv.second);
      check(name + " (" + kv.first + "): round trip", d.Write() == kv.second);
      if (kv.first == "as is") check(name + ": no structural problems", d.Problems().empty());
      if (kv.first == "padded") {
        const Document d0 = Document::Parse(bytes);
        check(name + " (padded): the same sections", d.Sections().size() == d0.Sections().size());
      }
    }
    // Each top-level section removed in turn: the bytes before and after it are untouched.
    const Document d0 = Document::Parse(bytes);
    for (std::size_t k = 0; k < d0.Sections().size(); k++) {
      const Document::Section s = d0.Sections()[k];
      if (s.depth != 0 || s.close == Document::Section::npos) continue;
      sections++;
      Document d = Document::Parse(bytes);
      d.RemoveSection(d.Sections()[k]);
      std::string before, after;
      for (std::size_t i = 0; i < s.open; i++) before += d0.Line(i) + d0.Ending(i);
      for (std::size_t i = s.close + 1; i < d0.NumLines(); i++) after += d0.Line(i) + d0.Ending(i);
      check(name + ": remove <" + s.name + ">", d.Write() == before + after);
    }
  }
  std::printf("  %zu files, %ld variants, %ld section removals\n", files.size(), variants, sections);
  if (files.empty()) {
    std::printf("  FAIL  no .azr found\n");
    failures++;
  }
  std::printf("%s: %ld checks, %d failure(s)\n", failures ? "FAIL" : "PASS", checks, failures);
  return failures ? 1 : 0;
}
