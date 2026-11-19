// azure2::front::ReadCurves (Engine.h; docs/upgrade/qt_gui.md section 5): the AZUREOut_* files of the last run as data,
// matched to the project's records by the engine's own build of the project.  The writer is
// WorkingCopy::WriteSegmentFiles (src/DataOutputFiles.cpp): one file per (entrance, exit) pair key, AZUREOut_aa=<a>_R=<r>
// (.out for the data sets, .extrap for the extrapolations, .acoeff for a Legendre extrapolation, _TOTAL_CAPTURE for an
// exit of every capture pair), one block per segment in segment order, each block's rows followed by two empty lines,
// and a .band sibling with the same blocks when the run wrote the uncertainty band.  The segments are the active data
// sets (or extrapolations) in record order, less a data set without points; chiSquared.out lists the keys of those
// that were used with their number of points, which this reader checks every block against.
#include "azure2/front/Engine.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>

#include "azure2/experiment/ProjectPaths.h"
#include "azure2/experiment/ProjectReader.h"

namespace azure2 {
namespace front {

namespace {

namespace ex = azure2::experiment;
const double kNaN = std::numeric_limits<double>::quiet_NaN();

using Block = std::vector<std::vector<double>>;

/// The blocks of an output file: rows of numbers, a block ended by an empty line (the writer's two empty lines make
/// one separator).  False when the file cannot be read.
bool ReadBlocks(const std::string &path, std::vector<Block> &blocks) {
  std::ifstream in(path);
  if (!in) return false;
  std::string line;
  Block current;
  while (std::getline(in, line)) {
    if (line.find_first_not_of(" \t\r") == std::string::npos) {
      if (!current.empty()) blocks.push_back(current);
      current.clear();
      continue;
    }
    if (line[line.find_first_not_of(" \t")] == '#') continue;  // a comment line (the sample decompositions)
    std::istringstream s(line);
    std::vector<double> row;
    std::string token;
    while (s >> token) {
      char *end = nullptr;
      const double v = std::strtod(token.c_str(), &end);
      row.push_back(end && *end == '\0' ? v : kNaN);
    }
    current.push_back(row);
  }
  if (!current.empty()) blocks.push_back(current);
  return true;
}

std::string FileName(int entrance, int exit, bool extrapolation, bool legendreExtrapolation) {
  std::string name = "AZUREOut_aa=" + std::to_string(entrance);
  name += exit == -1 ? std::string("_TOTAL_CAPTURE") : "_R=" + std::to_string(exit);
  if (!extrapolation) return name + ".out";
  if (exit != -1 && legendreExtrapolation) return name + ".acoeff";
  return name + ".extrap";
}

std::string ObservableWord(ex::Observable o) {
  switch (o) {
    case ex::Observable::CrossSection:
    case ex::Observable::CrossSectionE1:
    case ex::Observable::CrossSectionE2: return "cross_section";
    case ex::Observable::TotalCapture: return "total_capture";
    case ex::Observable::DifferentialLab:
    case ex::Observable::DifferentialCM: return "differential";
    case ex::Observable::AnalyzingPower: return "analyzing_power";
    case ex::Observable::PhaseShift: return "phase_shift";
    case ex::Observable::AngularCoefficients:
    case ex::Observable::AngularDistribution: return "legendre";
  }
  return "cross_section";
}

/// chiSquared.out: key -> N of the data sets the run used ("Segment#, Chi-Squared, N, Norm, ...").
std::map<int, int> UsedDataSets(const std::string &path, bool &present) {
  std::map<int, int> used;
  std::ifstream in(path);
  present = (bool)in;
  std::string line;
  while (std::getline(in, line)) {
    int key = 0, n = 0;
    double chi = 0;
    if (std::sscanf(line.c_str(), "%d,%lf,%d", &key, &chi, &n) == 3) used[key] = n;
  }
  return used;
}

struct Item {
  int key;
  bool active;
  int entrance, exit;
  ex::Observable observable;
};

}  // namespace

Curves ReadCurves(const Project &project) {
  Curves out;
  ex::Format3ReadOptions read;
  read.fileName = project.Path();
  read.projectDir = project.Directory();
  read.setUpChecks = false;
  const ex::Format3Project built = ex::ReadFormat3Text(project.Text(), read);
  if (!built.Ok()) {
    out.notes.push_back("the project has errors (Validate lists them): its curves cannot be matched to its records");
    return out;
  }
  out.outputDirectory = ex::PathToOpen(project.Directory(), project.ConfigValue("output_dir"));
  if (!out.outputDirectory.empty() && out.outputDirectory.back() != '/') out.outputDirectory += '/';
  {
    std::ifstream kin(ex::PathToOpen(project.Directory(), "kinflag.dat"));
    int flag = 0;
    if (kin >> flag) out.labFrame = flag != 0;
  }

  for (int pass = 0; pass < 2; pass++) {
    const bool extrapolation = pass == 1;
    const std::string section = extrapolation ? "extrapolations" : "data";
    const std::string kind = extrapolation ? "extrapolation" : "dataset";
    std::vector<Record> records;
    for (const Record &r : project.Records(section))
      if (r.kind == kind) records.push_back(r);
    std::vector<Item> items;
    if (extrapolation) {
      for (const auto &i : built.built.extrapolations.extrapolations.items)
        items.push_back({i.spec.key, i.spec.active, i.spec.reaction.entranceKey, i.spec.reaction.exitKey,
                         i.spec.reaction.observable});
    } else {
      for (const auto &i : built.built.data.datasets.items)
        items.push_back({i.spec.key, i.spec.active, i.spec.reaction.entranceKey, i.spec.reaction.exitKey,
                         i.spec.reaction.observable});
    }
    bool haveChi = false;
    const std::map<int, int> used = extrapolation ? std::map<int, int>()
                                                  : UsedDataSets(out.outputDirectory + "chiSquared.out", haveChi);
    std::map<std::string, std::vector<Block>> files, bands;
    std::map<std::string, int> nextBlock;
    std::vector<std::string> missing;
    for (const Item &item : items) {
      if (!item.active) continue;
      if (haveChi && !used.count(item.key)) continue;  // left out by the run (no points)
      const bool legendre = item.observable == ex::Observable::AngularDistribution ||
                            item.observable == ex::Observable::AngularCoefficients;
      const int exit = item.observable == ex::Observable::TotalCapture ? -1 : item.exit;
      const std::string file = FileName(item.entrance, exit, extrapolation,
                                        item.observable == ex::Observable::AngularDistribution);
      if (!files.count(file)) {
        std::vector<Block> blocks, bandBlocks;
        if (!ReadBlocks(out.outputDirectory + file, blocks)) missing.push_back(file);
        ReadBlocks(out.outputDirectory + file + ".band", bandBlocks);
        files[file] = blocks;
        bands[file] = bandBlocks;
      }
      const int block = nextBlock[file]++;
      Curve c;
      c.section = section;
      c.key = item.key;
      if (item.key >= 1 && item.key <= (int)records.size()) {
        const Record &r = records[item.key - 1];
        c.record = r.name;
        c.reaction = r.Value("reaction");
        c.quantity = r.Value("quantity");
        c.unit = r.Value("unit");
      }
      c.observable = ObservableWord(item.observable);
      c.sFactor = item.observable == ex::Observable::CrossSection || item.observable == ex::Observable::CrossSectionE1 ||
                  item.observable == ex::Observable::CrossSectionE2 || item.observable == ex::Observable::TotalCapture ||
                  item.observable == ex::Observable::DifferentialLab ||
                  item.observable == ex::Observable::DifferentialCM;
      c.positive = c.sFactor;
      c.file = file;
      c.block = block;
      const std::vector<Block> &blocks = files[file];
      if (block >= (int)blocks.size()) {
        if (!blocks.empty())
          out.notes.push_back(file + " has " + std::to_string(blocks.size()) + " block(s), fewer than the project's " +
                              section + " written to it: the file is older than the project (run again)");
        continue;
      }
      const Block &rows = blocks[block];
      if (haveChi && used.at(item.key) != (int)rows.size()) {
        out.notes.push_back(c.record + ": " + std::to_string(rows.size()) + " rows in " + file + ", " +
                            std::to_string(used.at(item.key)) + " points in chiSquared.out: not shown");
        continue;
      }
      if (legendre) {
        c.coefficients = rows;
      } else {
        const std::vector<Block> &bandBlocks = bands[file];
        const bool band = block < (int)bandBlocks.size() && bandBlocks[block].size() == rows.size();
        c.haveBand = band;
        for (std::size_t i = 0; i < rows.size(); i++) {
          const std::vector<double> &r = rows[i];
          auto at = [&r](std::size_t k) { return k < r.size() ? r[k] : kNaN; };
          CurvePoint p;
          p.energy = at(0);
          p.excitation = at(1);
          p.angle = at(2);
          p.model = at(3);
          p.modelS = at(4);
          p.data = extrapolation ? kNaN : at(5);
          p.dataError = extrapolation ? kNaN : at(6);
          p.dataS = extrapolation ? kNaN : at(7);
          p.dataSError = extrapolation ? kNaN : at(8);
          p.band = band && bandBlocks[block][i].size() > 6 ? bandBlocks[block][i][4] : kNaN;
          p.bandS = band && bandBlocks[block][i].size() > 6 ? bandBlocks[block][i][6] : kNaN;
          c.points.push_back(p);
        }
      }
      out.curves.push_back(c);
    }
    if (!missing.empty() && missing.size() == files.size()) {
      out.notes.push_back(std::string("no ") + (extrapolation ? "extrapolate" : "calculate or fit") +
                          " output in the output directory yet");
    } else {
      for (const std::string &f : missing) out.notes.push_back(f + ": not written by the last run");
    }
  }
  return out;
}

}  // namespace front
}  // namespace azure2
