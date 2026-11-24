// The old sections of a project as the old engine read them (OldProject.h).  Every line goes through the one parser
// of its grammar (the headers of convert/include/azure2/convert/, which the engine also uses until step 3b); this
// file only walks the sections and refuses what the old engine refused.
#include "OldProject.h"

#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <map>
#include <stdexcept>
#include <sstream>
#include <vector>

#include "azure2/io/LevelRows.h"  // the <levels> line (format 3 keeps that grammar); the converter reads only the pair columns
#include "azure2/convert/LineRepair.h"
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

// The reason a <segmentsData> or <segmentsTest> line the old parser cannot read is one of another AZURE2 (its dev
// branch), "" otherwise: its type (the token at 0-based place `typeAt`) followed by a spin operator such as b:y,
// t:y, b:y*t:y or b:T20, where this program's grammar has the normalization.  dev's general spin observables are
// type 10 with the operator (read: SegLine.h, ExtrapLine.h; a spin_observable record); before 2026-11-11 dev wrote
// them as type 9, which is Legendre-coefficient data since, and that line is refused with the way out.
std::string NewerLineType(const std::string &line, std::size_t typeAt) {
  std::istringstream in(line);
  std::vector<std::string> tokens;
  for (std::string t; in >> t;) tokens.push_back(t);
  if (tokens.size() <= typeAt + 1 || tokens[typeAt + 1].find(':') == std::string::npos) return "";
  if (tokens[typeAt] == "10") return "";
  return " has type " + tokens[typeAt] + " followed by the spin operator '" + tokens[typeAt + 1] +
         "': a general spin observable (an analysing power, a spin correlation, a polarization or a tensor power) is "
         "type 10 followed by its operator (AZURE2's dev branch wrote type 9 before 2026-11-11).  Write 10 for the "
         "type on this line";
}

// The refusal of a line's end that no keyword block of the old grammar reads: the old engine stopped reading there
// without a word, a newer program may read it; converted, it would be lost.
std::string UnreadTail(const std::string &tail) {
  return " ends with '" + Trim(tail) + "', which no keyword block of this program's grammar reads (the old engine "
         "stopped reading the line there without a word; a newer AZURE2 may read it).  Correct the keyword, or delete "
         "it from a copy of the project if the program that wrote it ignored it too, and convert the copy";
}

/// A line whose first word starts with '#': a comment, which the old engine could not read (it refused the file).
bool CommentLine(const std::string &line) {
  const std::size_t a = line.find_first_not_of(" \t\r");
  return a != std::string::npos && line[a] == '#';
}

/// The line with the blanks inside its double-quoted tokens removed ("8, 9" -> "8,9"); the old GUI wrote the segment
/// list and the equations of a <targetInt> line without blanks, and the old reader read each as one token.
std::string QuotedBlanksRemoved(const std::string &line) {
  std::string out;
  bool inside = false;
  for (char c : line) {
    if (c == '"') inside = !inside;
    if (inside && (c == ' ' || c == '\t')) continue;
    out += c;
  }
  return out;
}

/// A data line whose data file's name holds blanks, unquoted: the line with the name quoted when the words from the
/// file's place to the next number name a file that exists under `dir` (and the first word alone does not), else "".
std::string QuotedDataFile(const std::string &line, const std::string &dir) {
  if (dir.empty()) return "";
  const std::vector<LineReading> readings = ProposeReadings("segmentsData", line);
  for (const LineReading &r : readings) {
    const std::size_t a = r.line.find('"');
    const std::size_t b = a == std::string::npos ? a : r.line.find('"', a + 1);
    if (b == std::string::npos || r.label.compare(0, 23, "the data file's name ho") != 0) continue;
    const std::string name = r.line.substr(a + 1, b - a - 1);
    std::error_code ec;
    const std::string first = name.substr(0, name.find(' '));
    if (std::filesystem::is_regular_file(std::filesystem::path(dir) / name, ec) &&
        !std::filesystem::exists(std::filesystem::path(dir) / first, ec))
      return r.line;
  }
  return "";
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

std::string ReadOldProject(const io::Document &doc, OldProject &out, const std::string &projectDir) {
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
      std::string line = doc.Line(i);
      if (BlankLine(line)) continue;
      if (CommentLine(line)) {
        out.repairs.push_back({"<segmentsData> line", "file line " + std::to_string(i + 1) + ", a comment in <segmentsData> (\"" +
                                                          Trim(line) + "\"), is skipped: the old engine could not read it"});
        continue;
      }
      out.refusedSection = "segmentsData";
      out.refusedLine = i;
      if (!OldLineReads("segmentsData", line)) {
        const std::string quoted = QuotedDataFile(line, projectDir);
        if (!quoted.empty()) {
          out.repairs.push_back({"<segmentsData> line", "file line " + std::to_string(i + 1) + ": the data file's name holds a "
                                 "blank and was not quoted; read as the file that exists: " + quoted.substr(quoted.find('"'))});
          line = quoted;
        }
      }
      std::istringstream stm(line);
      auto parsed = std::make_shared<SegLine>(stm);
      key++;
      // a newer program's type first: its operator token can leave the old parser's reading intact (a test line)
      const std::string newer = NewerLineType(line, 7);
      if (!newer.empty()) return "data segment #" + std::to_string(key) + newer;
      if (stm.rdstate() & (std::ios::failbit | std::ios::badbit)) {
        return "data segment #" + std::to_string(key) + " cannot be read" +
               (parsed->parseError().empty() ? "" : ": " + parsed->parseError());
      }
      if (!parsed->unreadTail().empty())
        return "data segment #" + std::to_string(key) + UnreadTail(parsed->unreadTail());
      if (parsed->dataFile().empty()) return "data segment #" + std::to_string(key) + " names no data file";
      if (DecimalCommaNumber(parsed->dataFile()))
        return "data segment #" + std::to_string(key) + ": '" + parsed->dataFile() + "' stands where the data file's name "
               "is: a number written with a decimal comma ends the numbers of the line there";
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
      if (CommentLine(line)) {
        out.repairs.push_back({"<segmentsTest> line", "file line " + std::to_string(i + 1) + ", a comment in <segmentsTest> (\"" +
                                                          Trim(line) + "\"), is skipped: the old engine could not read it"});
        continue;
      }
      out.refusedSection = "segmentsTest";
      out.refusedLine = i;
      std::istringstream stm(line);
      auto parsed = std::make_shared<ExtrapLine>(stm);
      key++;
      // a newer program's type first: its operator token can leave the old parser's reading intact (a test line)
      const std::string newer = NewerLineType(line, 9);
      if (!newer.empty()) return "test segment #" + std::to_string(key) + newer;
      if (stm.rdstate() & (std::ios::failbit | std::ios::badbit)) {
        return "test segment #" + std::to_string(key) + " cannot be read" +
               (parsed->parseError().empty() ? "" : ": " + parsed->parseError());
      }
      if (!parsed->unreadTail().empty())
        return "test segment #" + std::to_string(key) + UnreadTail(parsed->unreadTail());
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
      std::string line = doc.Line(i);
      if (BlankLine(line)) continue;
      if (CommentLine(line)) {
        out.repairs.push_back({"<targetInt> line", "file line " + std::to_string(i + 1) + ", a comment in <targetInt> (\"" + Trim(line) +
                                                       "\"), is skipped: the old engine could not read it"});
        continue;
      }
      out.refusedSection = "targetInt";
      out.refusedLine = i;
      const std::string unblanked = QuotedBlanksRemoved(line);
      if (unblanked != line && !OldLineReads("targetInt", line) && OldLineReads("targetInt", unblanked)) {
        out.repairs.push_back({"<targetInt> line", "file line " + std::to_string(i + 1) + ": blanks inside a quoted token removed "
                               "(the old reader read a quoted list or equation as one word): " + unblanked});
        line = unblanked;
      }
      {
        // the nine fields end with an empty stopping equation and no parameter count: the count can only be 0
        std::istringstream words(line);
        std::vector<std::string> w;
        for (std::string t; words >> t;) w.push_back(t);
        if (w.size() == 8 && w[7] == "\"\"" && !OldLineReads("targetInt", line) && OldLineReads("targetInt", line + " 0")) {
          out.repairs.push_back({"<targetInt> line", "file line " + std::to_string(i + 1) + " ends after an empty stopping-power "
                                 "equation: its parameter count, the ninth field, is 0"});
          line += " 0";
        }
      }
      std::istringstream stm(line);
      OldEffectLine e;
      ReadTargetIntLine(stm, e.line);
      index++;
      const std::string where = "<targetInt> line \"" + line + "\"";
      if (e.line.screeningError) return where + ": a screening block is 'screening Ue_eV vary sigma_eV'";
      if (e.line.dopplerError) return where + ": " + e.line.dopplerErrorText;
      if (e.line.densityFitError) return where + ": " + e.line.densityFitErrorText;
      if (!e.line.unknownKeyword.empty())
        return where + ": '" + e.line.unknownKeyword + "' is not a keyword of a <targetInt> line (beamprofile, "
               "sigmalab, screening, densityfit, densitygroup, doppler, udr); the old engine dropped it without a word "
               "at the end of a line, and a newer AZURE2 may read it.  Correct it, or delete it from a copy of the "
               "project if the program that wrote it ignored it too, and convert the copy";
      if (stm.rdstate() & (std::ios::failbit | std::ios::badbit)) return where + " cannot be read";
      e.index = index;
      e.docLine = (int)i + 1;
      e.text = line;
      if (e.line.complete && e.line.isActive == 1) e.activeIndex = ++active;
      e.segments = TargetIntSegmentKeys(TargetIntUnquoted(e.line.segmentList));
      out.effects.push_back(e);
    }
  }

  out.refusedSection.clear();  // the refusals below are not about one line of those three sections
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

  // <isospin> of dev (2026-10-17): "reference Z A key", "shift Z A key value vary centre sigma", "level id T",
  // "mix id key l s value vary sigma", "factor id key value", "pairT key t1 t2", "soft id sigmaE sigmaGamma", as
  // IsospinLevels::Read read them.  A block whose first line is an isospin_* record is format 3 already and is left
  // alone.
  if (const io::Document::Section *s = doc.Find("isospin")) {
    if (s->close == io::Document::Section::npos) return "<isospin> is not closed";
    const std::map<std::string, std::size_t> kFields = {{"reference", 3}, {"shift", 7}, {"level", 2}, {"mix", 7},
                                                        {"factor", 3},    {"pairT", 3}, {"soft", 3}};
    bool first = true;
    for (const std::string &raw : doc.Body(*s)) {
      const std::string line = raw.substr(0, raw.find('#'));
      if (BlankLine(line)) continue;
      std::istringstream in(line);
      OldProject::IsospinLine c;
      c.text = raw;
      in >> c.kind;
      if (first && c.kind.compare(0, 8, "isospin_") == 0) break;
      first = false;
      std::string field;
      while (in >> field) c.fields.push_back(field);
      const auto known = kFields.find(c.kind);
      if (known == kFields.end())
        return "<isospin> line \"" + raw + "\": unknown keyword \"" + c.kind + "\" (reference, shift, level, mix, factor, pairT, soft)";
      if (c.fields.size() != known->second)
        return "<isospin> line \"" + raw + "\": " + c.kind + " takes " + std::to_string(known->second) + " fields";
      for (const std::string &f : c.fields) {
        try {
          std::size_t used = 0;
          (void)std::stod(f, &used);
          if (used != f.size()) throw std::invalid_argument(f);
        } catch (...) {
          return "<isospin> line \"" + raw + "\": \"" + f + "\" is not a number";
        }
      }
      out.isospin.push_back(c);
    }
  }

  // <channelRadii> of dev (2026-11-10): "<pair key> <l> <radius> [s=<s>] [J=<J>]" a line, as CNuc read it.  A block
  // whose first line is a `channel_radius` record is format 3 already and is left alone.
  if (const io::Document::Section *s = doc.Find("channelRadii")) {
    if (s->close == io::Document::Section::npos) return "<channelRadii> is not closed";
    bool records = false, first = true;
    for (const std::string &raw : doc.Body(*s)) {
      const std::string line = raw.substr(0, raw.find('#'));
      if (BlankLine(line)) continue;
      std::istringstream stm(line);
      std::string word;
      stm >> word;
      if (first && word == "channel_radius") records = true;
      first = false;
      if (records) break;
      const std::string bad = "<channelRadii> line \"" + raw + "\" is not \"<pair key> <l> <radius> [s=<s>] [J=<J>]\"";
      OldProject::ChannelRadiusLine c;
      c.text = raw;
      std::istringstream in(line);
      double radius = 0.0;
      std::string radiusText;
      if (!(in >> c.key >> c.l >> radiusText)) return bad;
      try {
        std::size_t used = 0;
        radius = std::stod(radiusText, &used);
        if (used != radiusText.size()) return bad;
      } catch (...) {
        return bad;
      }
      if (c.l < 0 || !(radius > 0.0)) return "<channelRadii> line \"" + raw + "\": l must be >= 0 and the radius positive";
      c.radius = radiusText;
      std::string token;
      while (in >> token) {
        const bool isS = token.compare(0, 2, "s=") == 0, isJ = token.compare(0, 2, "J=") == 0;
        const std::string value = token.substr(2);
        bool ok = (isS || isJ) && !value.empty();
        // <n> or <n>/2
        const std::size_t slash = value.find('/');
        const std::string num = value.substr(0, slash);
        ok = ok && !num.empty() && num.find_first_not_of("0123456789") == std::string::npos &&
             (slash == std::string::npos || value.substr(slash) == "/2");
        if (!ok || (isS && !c.s.empty()) || (isJ && !c.J.empty())) return bad + ": \"" + token + "\" is not s=<s> or J=<J>";
        (isS ? c.s : c.J) = value;
      }
      out.channelRadii.push_back(c);
    }
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
