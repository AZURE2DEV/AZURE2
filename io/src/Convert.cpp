// Format-1 -> format-2 projects (Convert.h; docs/upgrade/project.md section 4).
#include "azure2/io/Convert.h"

#include <string>

namespace azure2 {
namespace io {

namespace {

std::string Trim(const std::string &s) {
  std::size_t a = 0, b = s.size();
  while (a < b && (s[a] == ' ' || s[a] == '\t')) a++;
  while (b > a && (s[b - 1] == ' ' || s[b - 1] == '\t')) b--;
  return s.substr(a, b - a);
}

// key=value of a line (comment after '#' dropped, both sides trimmed); false if the line has no '='.
bool KeyValue(const std::string &line, std::string &key, std::string &value) {
  std::string entry = line.substr(0, line.find('#'));
  const std::size_t eq = entry.find('=');
  if (eq == std::string::npos) return false;
  key = Trim(entry.substr(0, eq));
  value = Trim(entry.substr(eq + 1));
  return true;
}

}  // namespace

int FormatOf(const Document &doc) {
  const Document::Section *c = doc.Find("config");
  if (!c) return 1;
  for (const std::string &line : doc.Body(*c)) {
    std::string key, value;
    if (KeyValue(line, key, value) && key == "format") return value == "2" ? 2 : 1;
  }
  return 1;
}

std::vector<std::string> ConvertToFormat2(Document &doc) {
  std::vector<std::string> report;
  if (FormatOf(doc) == 2) return report;

  // <lastRun>
  if (const Document::Section *s = doc.Find("lastRun")) {
    if (s->close != Document::Section::npos) {
      doc.RemoveSection(*s);
      report.push_back("removed <lastRun> (the Qt GUI's runtime options; read by nothing)");
    }
  }

  // <mcmc>: keep its <parameters> only.
  if (const Document::Section *m = doc.Find("mcmc")) {
    if (m->close != Document::Section::npos) {
      const Document::Section *p = doc.FindChild(*m, "parameters");
      std::vector<std::string> keep;
      bool hasPriors = false;
      if (p && p->close != Document::Section::npos) {
        for (std::size_t i = p->open; i <= p->close; i++) keep.push_back(doc.Line(i));
        for (const std::string &line : doc.Body(*p))
          if (!Trim(line).empty()) hasPriors = true;
      }
      const std::size_t first = m->open, last = m->close + 1;
      if (hasPriors) {
        std::vector<std::string> lines;
        lines.push_back(doc.Line(first));
        lines.insert(lines.end(), keep.begin(), keep.end());
        lines.push_back(doc.Line(last - 1));
        doc.ReplaceLines(first, last, lines);
        report.push_back("removed <mcmc>'s run settings (the CLI asks for them); its <parameters> priors kept");
      } else {
        doc.ReplaceLines(first, last, {});
        report.push_back("removed <mcmc> (run settings only, no <parameters> priors)");
      }
    }
  }

  // <potential>: useAdaptiveGrid to <config>; the block only with a hybrid potential.
  std::string adaptiveGrid;
  if (const Document::Section *pot = doc.Find("potential")) {
    if (pot->close != Document::Section::npos) {
      bool hybrid = false;
      std::size_t gridLine = Document::Section::npos;
      for (std::size_t i = pot->open + 1; i < pot->close; i++) {
        std::string key, value;
        if (!KeyValue(doc.Line(i), key, value)) continue;
        if (key == "useHybridPotential" && value == "1") hybrid = true;  // any pair (pair= sections) with one
        if (key == "useAdaptiveGrid") {
          adaptiveGrid = value;
          gridLine = i;
        }
      }
      if (!hybrid) {
        doc.RemoveSection(*pot);
        report.push_back("removed <potential> (no hybrid potential)" +
                         std::string(adaptiveGrid.empty() ? "" : "; useAdaptiveGrid moved to <config>"));
      } else if (gridLine != Document::Section::npos) {
        doc.ReplaceLines(gridLine, gridLine + 1, {});
        report.push_back("moved useAdaptiveGrid from <potential> to <config>");
      }
    }
  }

  // <config>: the moved key and format=2, as its last lines.
  if (const Document::Section *c = doc.Find("config")) {
    if (c->close != Document::Section::npos) {
      std::vector<std::string> add;
      if (!adaptiveGrid.empty()) add.push_back("useAdaptiveGrid=" + adaptiveGrid);
      add.push_back("format=2");
      const std::size_t at = c->close;
      add.push_back(doc.Line(at));
      doc.ReplaceLines(at, at + 1, add);
      report.push_back("<config>: format=2");
    }
  }
  return report;
}

}  // namespace io
}  // namespace azure2
