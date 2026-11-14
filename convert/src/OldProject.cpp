// The old sections of a project as the old engine read them (OldProject.h).  Every line goes through the one parser
// of its grammar (the headers of convert/include/azure2/convert/, which the engine also uses until step 3b); this
// file only walks the sections and refuses what the old engine refused.
#include "OldProject.h"

#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>

#include "azure2/io/LevelRows.h"  // the <levels> line (format 3 keeps that grammar); the converter reads only the pair columns
#include "azure2/convert/PriorCentreRow.h"
#include "azure2/convert/SpecsFromOldLines.h"
#include "azure2/convert/ThmLines.h"

namespace azure2 {
namespace convert {

namespace {

bool BlankLine(const std::string &line) {
  for (char c : line)
    if (c != ' ' && c != '\t') return false;
  return true;
}

std::string Trim(const std::string &s) {
  const std::size_t a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  const std::size_t b = s.find_last_not_of(" \t\r\n");
  return s.substr(a, b - a + 1);
}

}  // namespace

int OldPair::A1() const { return (int)std::lround(m1); }
int OldPair::A2() const { return (int)std::lround(m2); }

experiment::PairInfo OldPair::Info() const {
  experiment::PairInfo p;
  p.key = key;
  p.kind = pType;
  p.lightZ = z1;
  p.lightA = A1();
  p.heavyZ = z2;
  p.heavyA = A2();
  p.lightExcitation = e1;
  p.heavyExcitation = e2;
  p.lightMass = m1;
  p.heavyMass = m2;
  p.lightExcitationText = e1Text;
  p.heavyExcitationText = e2Text;
  return p;
}

std::string ShortestText(double value) {
  // The shortest decimal that reads back to the same double, in fixed or scientific form, whichever is shorter
  // (std::to_chars without a format: 90, 0.25, 5e+17).
  char buf[64];
  const std::to_chars_result r = std::to_chars(buf, buf + sizeof buf, value);
  return std::string(buf, r.ptr);
}

std::string NumberText(double value, const std::string &line) {
  // The longest such token, the first of equal length: "1.0" rather than a count "1" before it (any of them reads
  // back to the same double; the longest is the one written as a value).
  std::istringstream tokens(line);
  std::string t, best;
  while (tokens >> t) {
    char *end = nullptr;
    const double v = std::strtod(t.c_str(), &end);
    if (end != t.c_str() && *end == '\0' && std::memcmp(&v, &value, sizeof v) == 0 && t.size() > best.size()) best = t;
  }
  return best.empty() ? ShortestText(value) : best;
}

std::string SpinText(double j) {
  const long twice = std::lround(2 * j);
  if (twice % 2 == 0) return std::to_string(twice / 2);
  return std::to_string(twice) + "/2";
}

std::string ReadOldProject(const io::Document &doc, OldProject &out) {
  for (const char *name : {"levels", "segmentsData", "segmentsTest", "targetInt", "abundances", "scatteringLengths", "thm",
                           "parameterSettings", "config", "potential"})
    if (const io::Document::Section *s = doc.Find(name))
      if (s->close == io::Document::Section::npos) return std::string("<") + name + "> is not closed";
  // <levels>: the pairs, from the first line naming each (as CompoundSystem::Fill reads them), active or not.
  if (const io::Document::Section *s = doc.Find("levels")) {
    for (const std::string &line : doc.Body(*s)) {
      if (BlankLine(line)) continue;
      std::istringstream stm(line);
      io::LevelRow nl(stm);
      if (stm.rdstate() & (std::ios::failbit | std::ios::badbit)) return "a <levels> line the old engine cannot read: \"" + line + "\"";
      OldChannel ch;
      ch.pairKey = nl.ir();
      ch.J = nl.levelJ();
      ch.l = nl.l();
      ch.s = nl.s();
      ch.active = nl.isActive() == 1;
      out.channels.push_back(ch);
      if (out.pairs.count(nl.ir())) continue;
      OldPair p;
      p.key = nl.ir();
      p.pType = nl.pType();
      p.j1 = nl.j1();
      p.j2 = nl.j2();
      p.m1 = nl.m1();
      p.m2 = nl.m2();
      p.z1 = nl.z1();
      p.z2 = nl.z2();
      p.e1 = nl.e1();
      p.e2 = nl.e2();
      p.e1Text = NumberText(p.e1, line);
      p.e2Text = NumberText(p.e2, line);
      out.pairs[p.key] = p;
    }
  } else {
    return "the project has no <levels> section";
  }

  // <segmentsData>: WorkingCopy::Fill's refusals that do not depend on the level scheme.
  if (const io::Document::Section *s = doc.Find("segmentsData")) {
    int key = 0;
    for (std::size_t i = s->open + 1; i < s->close; i++) {
      const std::string &line = doc.Line(i);
      if (BlankLine(line)) continue;
      std::istringstream stm(line);
      auto parsed = std::make_shared<SegLine>(stm);
      key++;
      if (stm.rdstate() & (std::ios::failbit | std::ios::badbit))
        return "data segment #" + std::to_string(key) + " cannot be read" +
               (parsed->parseError().empty() ? "" : ": " + parsed->parseError());
      const int d = parsed->isDiff();
      if (d == 19) return "data segment #" + std::to_string(key) + " uses isDiff 19 (THM angular-distribution coefficients), which the old engine refused";
      if (d == 8 || d == 18) return "data segment #" + std::to_string(key) + " uses isDiff " + std::to_string(d) + " (polarization x cross section), which the old engine refused";
      if (parsed->isAdvanced() == 1 && parsed->operationType() != 0 && parsed->operationType() != 1)
        return "data segment #" + std::to_string(key) + ": its composite block has operation type " +
               std::to_string(parsed->operationType()) + ", which is neither 0 (sum) nor 1 (ratio)";
      OldDataLine dl;
      dl.key = key;
      dl.docLine = (int)i + 1;
      dl.text = line;
      dl.parsed = parsed;
      out.data.push_back(dl);
    }
  }

  // <segmentsTest>: WorkingCopy::MakePoints's refusals.
  if (const io::Document::Section *s = doc.Find("segmentsTest")) {
    int key = 0;
    for (std::size_t i = s->open + 1; i < s->close; i++) {
      const std::string &line = doc.Line(i);
      if (BlankLine(line)) continue;
      std::istringstream stm(line);
      auto parsed = std::make_shared<ExtrapLine>(stm);
      key++;
      if (stm.rdstate() & (std::ios::failbit | std::ios::badbit))
        return "test segment #" + std::to_string(key) + " cannot be read" +
               (parsed->parseError().empty() ? "" : ": " + parsed->parseError());
      const int d = parsed->isDiff();
      if (parsed->hasCovariance()) return "test segment #" + std::to_string(key) + " carries a covariance block, which the old engine refused";
      if (d == 9 || d == 19) return "test segment #" + std::to_string(key) + " uses isDiff " + std::to_string(d) + " (angular coefficients, a data-segment type), which the old engine refused";
      if (d == 8 || d == 18) return "test segment #" + std::to_string(key) + " uses isDiff " + std::to_string(d) + " (polarization x cross section), which the old engine refused";
      if (parsed->isAdvanced() == 1 && parsed->operationType() != 0 && parsed->operationType() != 1)
        return "test segment #" + std::to_string(key) + ": its composite block has operation type " +
               std::to_string(parsed->operationType()) + ", which is neither 0 (sum) nor 1 (ratio)";
      OldTestLine tl;
      tl.key = key;
      tl.docLine = (int)i + 1;
      tl.text = line;
      tl.parsed = parsed;
      out.tests.push_back(tl);
    }
  }

  // <targetInt>: WorkingCopy::ReadTargetEffectsFile's refusals of a line.
  if (const io::Document::Section *s = doc.Find("targetInt")) {
    int index = 0, active = 0;
    for (std::size_t i = s->open + 1; i < s->close; i++) {
      const std::string &line = doc.Line(i);
      if (BlankLine(line)) continue;
      std::istringstream stm(line);
      OldEffectLine e;
      ReadTargetIntLine(stm, e.line);
      index++;
      const std::string where = "<targetInt> line \"" + line + "\"";
      if (e.line.screeningError) return where + ": a screening block is 'screening Ue_eV vary sigma_eV'";
      if (e.line.dopplerError) return where + ": " + e.line.dopplerErrorText;
      if (e.line.densityFitError) return where + ": " + e.line.densityFitErrorText;
      if (stm.rdstate() & (std::ios::failbit | std::ios::badbit)) return where + " cannot be read";
      e.index = index;
      e.docLine = (int)i + 1;
      e.text = line;
      if (e.line.complete && e.line.isActive == 1) e.activeIndex = ++active;
      e.segments = TargetIntSegmentKeys(TargetIntUnquoted(e.line.segmentList));
      out.effects.push_back(e);
    }
  }

  // <abundances>
  if (const io::Document::Section *s = doc.Find("abundances")) {
    for (const std::string &line : doc.Body(*s)) {
      OldAbundance a;
      a.text = line;
      a.line = ReadAbundanceLine(line);
      if (a.line.kind == AbundanceLine::Blank) continue;
      if (a.line.kind == AbundanceLine::BadSum || a.line.kind == AbundanceLine::BadEntry)
        return "<abundances> line \"" + line + "\" cannot be read";
      out.abundances.push_back(a);
    }
  }

  // <potential>: the hybrid potentials, with the engine's parser (after the format-2 conversion it is there only when
  // a pair has the model on)
  if (const io::Document::Section *s = doc.Find("potential")) {
    std::string text = "<potential>\n";
    for (const std::string &line : doc.Body(*s)) text += line + "\n";
    text += "</potential>\n";
    std::istringstream in(text);
    out.potential = ReadPotentialBlock(in);
    out.hasPotential = true;
  }

  // <scatteringLengths>
  if (const io::Document::Section *s = doc.Find("scatteringLengths")) {
    int n = 0;
    for (const std::string &line : doc.Body(*s)) {
      n++;
      OldScatteringLength l;
      l.number = n;
      l.text = line;
      l.line = ReadScatteringLengthLine(line);
      if (l.line.kind == ScatteringLengthLine::Blank) continue;
      if (l.line.kind == ScatteringLengthLine::Malformed)
        return "<scatteringLengths> line " + std::to_string(n) + ": " + l.line.error;
      out.scatteringLengths.push_back(l);
    }
  }

  // <thm>: experiment[<name>] segments=..., weight[k]=, weightTest[k]=
  if (const io::Document::Section *s = doc.Find("thm")) {
    for (std::string line : doc.Body(*s)) {
      const std::size_t hash = line.find('#');
      if (hash != std::string::npos) line = line.substr(0, hash);
      const std::string trimmed = Trim(line);
      if (trimmed.empty()) continue;
      if (trimmed.compare(0, 11, "experiment[") == 0) {
        const std::size_t close = trimmed.find(']');
        if (close == std::string::npos) return "<thm> line \"" + trimmed + "\" cannot be read";
        const std::string name = trimmed.substr(11, close - 11);
        std::istringstream tokens(trimmed.substr(close + 1));
        std::string token;
        while (tokens >> token) {
          if (token.compare(0, 9, "segments=") != 0) continue;
          std::vector<int> keys;
          const std::string why = ParseThmSegmentList(token.substr(9), keys);
          if (!why.empty()) return "<thm> experiment[" + name + "]: " + why;
          // The member checks of the old engine (CheckThmExperiments, ThmNorm.cpp): a line in two experiments, a
          // number past the last data line, a line in use that is not THM.
          for (int k : keys) {
            const std::string where = "<thm> experiment[" + name + "]: segment " + std::to_string(k);
            const auto other = out.thm.experimentOfSegment.find(k);
            if (other != out.thm.experimentOfSegment.end() && other->second != name)
              return where + " is already in experiment[" + other->second + "], which the old engine refused";
            if (k > (int)out.data.size())
              return where + ": <segmentsData> has only " + std::to_string(out.data.size()) + " line(s), which the old engine refused";
            const SegLine &line = *out.data[k - 1].parsed;
            if (line.isActive() == 1 && line.isDiff() < 10)
              return where + " is not a THM segment (isDiff < 10), which the old engine refused";
            out.thm.experimentOfSegment[k] = name;
          }
        }
        continue;
      }
      const std::size_t eq = trimmed.find('=');
      if (eq == std::string::npos) continue;
      const std::string key = Trim(trimmed.substr(0, eq)), value = Trim(trimmed.substr(eq + 1));
      ThmWeightKey w;
      bool valid = false;
      if (ReadThmWeightKey(key, value, w, valid)) {
        if (!valid) return "<thm> line \"" + trimmed + "\" cannot be read";
        // The old engine's check of a data line's table (CheckThmWeights, WorkingCopy.cpp): a number past the last
        // line, or a line in use that is not THM, refuses every run with data.
        if (!w.test) {
          const std::string where = "<thm> weight[" + std::to_string(w.segment) + "]: ";
          if (w.segment > (int)out.data.size())
            return where + "<segmentsData> has only " + std::to_string(out.data.size()) + " line(s), which the old engine refused";
          const SegLine &line = *out.data[w.segment - 1].parsed;
          if (line.isActive() == 1 && line.isDiff() < 10)
            return where + "segment " + std::to_string(w.segment) + " is not a THM segment (isDiff < 10), which the old engine refused";
        }
        (w.test ? out.thm.weightOfTest : out.thm.weightOfData)[w.segment] = w.file;
      }
    }
  }

  // <parameterSettings>: the prior-centre rows.
  if (const io::Document::Section *s = doc.Find("parameterSettings")) {
    for (const std::string &line : doc.Body(*s)) {
      const std::string trimmed = Trim(line);
      if (trimmed.empty() || trimmed[0] == '#' || trimmed[0] == '<') continue;
      const PriorCentreRow row = ReadPriorCentreRow(trimmed);
      if (row.kind == PriorCentreRow::NotARow) continue;
      if (row.kind == PriorCentreRow::Malformed) return "<parameterSettings> row \"" + trimmed + "\" cannot be read";
      {
        // The engine's reading of the row (its whole name against the names of the <segmentsData> lines): a row it
        // refuses is refused here with its reason.
        experiment::SpecList<experiment::DatasetSpec> lines;
        lines.items.resize(out.data.size());
        std::istringstream one("<parameterSettings>\n" + trimmed + "\n");
        const std::string refusal = PriorCentresFromOldRows(&one, lines);
        if (!refusal.empty()) return refusal.substr(std::strlen("ERROR: "));
      }
      if (row.kind == PriorCentreRow::Norm) {
        if (row.centre <= 0.) return "<parameterSettings> row \"" + trimmed + "\": a normalization prior centre must be positive";
        out.normCentre[row.segment] = row.centre;
        out.normCentreText[row.segment] = row.centreText;
      } else {
        out.shiftCentre[row.segment] = row.centre;
        out.shiftCentreText[row.segment] = row.centreText;
      }
    }
  }

  // <config>: abundanceSumSigma= (Config::Read: key=value, '#' comments)
  if (const io::Document::Section *s = doc.Find("config")) {
    for (std::string line : doc.Body(*s)) {
      const std::size_t hash = line.find('#');
      if (hash != std::string::npos) line = line.substr(0, hash);
      const std::size_t eq = line.find('=');
      if (eq == std::string::npos) continue;
      const std::string key = Trim(line.substr(0, eq)), value = Trim(line.substr(eq + 1));
      if (key == "includeAllJPi") out.includeAllJPi = value == "1" || value == "true";
      if (key == "useAdaptiveGrid") out.adaptiveGridText = value;
      if (key != "abundanceSumSigma") continue;
      out.abundanceSumSigmaText = value;
      out.abundanceSumSigma = std::strtod(out.abundanceSumSigmaText.c_str(), nullptr);
    }
  }
  return "";
}

}  // namespace convert
}  // namespace azure2
