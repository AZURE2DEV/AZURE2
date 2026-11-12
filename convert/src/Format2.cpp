// Format 1 -> format 2, the converter's first step (Convert.h; docs/upgrade/project.md section 4).
#include <string>

#include "azure2/convert/Convert.h"
#include "azure2/io/ProjectFormat.h"

namespace azure2 {
namespace convert {

std::vector<std::string> ConvertToFormat2(io::Document &doc) {
  std::vector<std::string> report;
  if (io::FormatOf(doc) != 1) return report;

  // <lastRun>
  if (const io::Document::Section *s = doc.Find("lastRun")) {
    if (s->close != io::Document::Section::npos) {
      doc.RemoveSection(*s);
      report.push_back("removed <lastRun> (the Qt GUI's runtime options; read by nothing)");
    }
  }

  // <mcmc>: keep its <parameters> only.
  if (const io::Document::Section *m = doc.Find("mcmc")) {
    if (m->close != io::Document::Section::npos) {
      const io::Document::Section *p = doc.FindChild(*m, "parameters");
      std::vector<std::string> keep;
      bool hasPriors = false;
      if (p && p->close != io::Document::Section::npos) {
        for (std::size_t i = p->open; i <= p->close; i++) keep.push_back(doc.Line(i));
        for (const std::string &line : doc.Body(*p))
          if (line.find_first_not_of(" \t") != std::string::npos) hasPriors = true;
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
  if (const io::Document::Section *pot = doc.Find("potential")) {
    if (pot->close != io::Document::Section::npos) {
      bool hybrid = false;
      std::size_t gridLine = io::Document::Section::npos;
      for (std::size_t i = pot->open + 1; i < pot->close; i++) {
        std::string key, value;
        if (!io::ConfigKeyValue(doc.Line(i), key, value)) continue;
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
      } else if (gridLine != io::Document::Section::npos) {
        doc.ReplaceLines(gridLine, gridLine + 1, {});
        report.push_back("moved useAdaptiveGrid from <potential> to <config>");
      }
    }
  }

  // <config>: the moved key and format=2, as its last lines.
  if (const io::Document::Section *c = doc.Find("config")) {
    if (c->close != io::Document::Section::npos) {
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

}  // namespace convert
}  // namespace azure2
