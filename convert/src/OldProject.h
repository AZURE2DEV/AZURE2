#ifndef AZURE2_CONVERT_OLDPROJECT_H
#define AZURE2_CONVERT_OLDPROJECT_H

// The old sections of a project, read with the old line parsers (convert/include/azure2/convert/*Line.h) and kept
// with the text of every line, so that numbers are written back as they were written (NumberText).  Internal to
// the converter.
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "azure2/convert/AbundanceLine.h"
#include "azure2/convert/ComponentList.h"
#include "azure2/convert/ExtrapLine.h"
#include "azure2/convert/PotentialBlock.h"
#include "azure2/convert/ScatteringLengthLine.h"
#include "azure2/convert/SegLine.h"
#include "azure2/convert/TargetIntLine.h"
#include "azure2/experiment/LevelTables.h"
#include "azure2/io/Document.h"

namespace azure2 {
namespace convert {

/// A particle pair of <levels> (columns 13-32 of the first line naming it).
struct OldPair {
  int key = 0;
  int pType = 0;       ///< 0 particles, 10 capture, 20 beta-delayed
  double j1 = 0, j2 = 0;
  double m1 = 0, m2 = 0;  ///< u
  int z1 = 0, z2 = 0;
  double e1 = 0, e2 = 0;  ///< MeV: light and heavy particle excitation
  std::string e1Text, e2Text;
  int A1() const;
  int A2() const;
  /// The pair as the format-3 tables hold it (the label's fields: particles, excitations as written).
  experiment::PairInfo Info() const;
};

/// A channel of a <levels> line: the partial wave (J, l, s) of a pair.
struct OldChannel {
  int pairKey = 0;
  double J = 0;
  int l = 0;
  double s = 0;
  bool active = false;
};

/// A <segmentsData> line.
struct OldDataLine {
  int key = 0;  ///< 1-based place among the non-blank lines (the engine's segment key)
  int docLine = 0;  ///< 1-based line of the file
  std::string text;
  std::shared_ptr<const SegLine> parsed;
};

/// A <segmentsTest> line.
struct OldTestLine {
  int key = 0;
  int docLine = 0;
  std::string text;
  std::shared_ptr<const ExtrapLine> parsed;
};

/// A <targetInt> line.
struct OldEffectLine {
  int index = 0;        ///< 1-based place among the non-blank lines
  int activeIndex = 0;  ///< 1-based place among the active lines (the n of density_<n>, screening_<n>), 0 if inactive
  int docLine = 0;
  std::string text;
  TargetIntLine line;
  std::vector<int> segments;  ///< the listed keys
};

struct OldAbundance {
  std::string text;
  AbundanceLine line;
};

struct OldScatteringLength {
  int number = 0;  ///< the line's number in the block (the engine's messages)
  std::string text;
  ScatteringLengthLine line;
};

/// The <thm> section's references to data lines.
struct OldThm {
  std::map<int, std::string> experimentOfSegment;  ///< data key -> experiment name
  std::map<int, std::string> weightOfData, weightOfTest;  ///< key -> weight file
};

struct OldProject {
  std::map<int, OldPair> pairs;
  std::vector<OldChannel> channels;
  bool includeAllJPi = false;  ///< <config> includeAllJPi=1: hard-sphere channels in every J^pi
  std::vector<OldDataLine> data;
  std::vector<OldTestLine> tests;
  std::vector<OldEffectLine> effects;
  std::vector<OldAbundance> abundances;
  std::vector<OldScatteringLength> scatteringLengths;
  OldThm thm;
  std::map<int, double> normCentre, shiftCentre;  ///< prior_centre rows by segment key
  std::map<int, std::string> normCentreText, shiftCentreText;
  /// A line of dev's <channelRadii> block: "<pair key> <l> <radius> [s=<s>] [J=<J>]", the numbers as written.
  struct ChannelRadiusLine {
    int key = 0, l = 0;
    std::string radius, s, J;  ///< s, J: "" when not given
    std::string text;
  };
  /// A line of dev's <isospin> block (2026-10-17): its keyword (reference, shift, level, mix, factor, pairT, soft) and
  /// its fields as written.
  struct IsospinLine {
    std::string kind;
    std::vector<std::string> fields;
    std::string text;
  };
  std::vector<IsospinLine> isospin;  ///< the lines of an old <isospin> block (none: no block, or records)
  std::vector<ChannelRadiusLine> channelRadii;  ///< the lines of an old <channelRadii> block (none: no block, or records)
  bool hasPotential = false;  ///< a <potential> block (the hybrid model on for some pair)
  PotentialBlock potential;   ///< as the engine reads it (its warnings included)
  double abundanceSumSigma = 0.0;  ///< <config> abundanceSumSigma= (the engine's global sum width)
  std::string abundanceSumSigmaText;
  /// <config> useAdaptiveGrid= (format 2 has it there; the last line wins): "" when absent, else its value.  Format 3
  /// has no such key: 0 chose a uniform sub-point grid of its own, which no longer exists.
  std::string adaptiveGridText;
  /// What the reader accepted that the old engine did not read, each unambiguous (docs/upgrade/qt_gui.md section
  /// 19.5): the construct and the note (a comment line skipped, blanks inside a quoted token, a file name with blanks
  /// quoted, an empty stopping equation's missing parameter count).
  std::vector<std::pair<std::string, std::string>> repairs;
  /// The line a refusal is about (refusals of a <segmentsData>, <segmentsTest> or <targetInt> line): its section and
  /// its index in the document read (doc.Line); "" when the refusal is about no such line.
  std::string refusedSection;
  std::size_t refusedLine = 0;
};

/// Reads the old sections of `doc` (format 2).  Returns "" or the reason the old engine refused the file.  `projectDir`:
/// where the data files are (a data file whose name holds blanks is recognised by its presence there).
std::string ReadOldProject(const io::Document &doc, OldProject &out, const std::string &projectDir = std::string());

/// The text of `value` as written in `line` when a whole blank-separated token of it parses to the same double
/// (bit for bit; the longest such token), else the shortest decimal that reads back to it.
std::string NumberText(double value, const std::string &line);
/// The shortest decimal that reads back to `value` (strtod of the result is `value` bit for bit).
std::string ShortestText(double value);
/// A spin as written in format 3: 0, 1, 3/2, ...
std::string SpinText(double j);

}  // namespace convert
}  // namespace azure2

#endif
