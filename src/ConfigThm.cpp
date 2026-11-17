/*
 * The <thm> section of a project (Config::ReadThmBlock, ReadFormat3ThmSection) and the tables it and the records name:
 * the THM part of Config.
 */
#include "Config.h"
#include "ThmBlockEdits.h"
#include "azure2/experiment/ProjectPaths.h"
#include "azure2/experiment/ProjectReader.h"
#include "azure2/experiment/ThmSectionLines.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>

namespace {

// A file named in the section or by a record, relative to the project file unless absolute (the one rule:
// experiment::IsAbsolutePath, which knows a drive letter).
std::string ThmRelativePath(const std::string &projectFile, const std::string &value) {
  std::string dir;
  size_t slash = projectFile.find_last_of("/\\");
  if (slash != std::string::npos) dir = projectFile.substr(0, slash + 1);
  return (azure2::experiment::IsAbsolutePath(value) || dir.empty()) ? value : dir + value;
}

// A weight table named by a record (a data set's or an extrapolation's thm_weight), read when `readFiles`.
std::string ReadWeight(const std::string &projectFile, const std::string &file, const std::string &key, bool readFiles,
                       std::map<int, std::shared_ptr<const ThmWeightTable>> &tables, int recordKey) {
  std::shared_ptr<ThmWeightTable> table = std::make_shared<ThmWeightTable>();
  table->name = file;
  table->path = ThmRelativePath(projectFile, file);
  if (readFiles) {
    const std::string why = table->Read(table->path);
    if (!why.empty()) return "<thm> " + key + ": " + why;
  }
  tables[recordKey] = table;
  return "";
}

}  // namespace

std::string ReadFormat3ThmSection(const std::vector<std::string> &section,
                                  const azure2::experiment::ThmMembershipSpec &members, const std::string &projectFile,
                                  bool readFiles, ThmOptions &thm, std::vector<std::string> &idle,
                                  std::vector<std::string> &warnings) {
  bool any = false;
  for (const auto &x : members.experiments) any = any || !x.datasetKeys.empty();
  std::vector<std::string> lines = section;
  if (lines.empty()) {
    // Without the section, the records alone: an experiment no section defines is an internal guard below.
    if (!any && members.datasetWeights.empty() && members.extrapolationWeights.empty()) return "";
    lines = {"<thm>", "</thm>"};
  }
  std::string why = CheckFormat3ThmMembers(lines, members);
  if (!why.empty()) return why;
  idle = ThmExperimentsWithoutMembers(lines, members);
  for (std::size_t n = 1; n < lines.size(); n++) {
    std::string line = lines[n];
    size_t hash = line.find('#');
    if (hash != std::string::npos) line = line.substr(0, hash);
    size_t b = line.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) continue;
    size_t e = line.find_last_not_of(" \t\r\n");
    std::string trimmed = line.substr(b, e - b + 1);
    if (trimmed == "</thm>") {
      // The weight tables the records name, then the members of each experiment (DatasetSpec::key, ascending).
      for (const auto &w : members.datasetWeights) {
        why = ReadWeight(projectFile, w.second, "weight[" + std::to_string(w.first) + "]", readFiles, thm.weightBySegment, w.first);
        if (!why.empty()) return why;
      }
      for (const auto &w : members.extrapolationWeights) {
        why = ReadWeight(projectFile, w.second, "weightTest[" + std::to_string(w.first) + "]", readFiles,
                         thm.weightByTestSegment, w.first);
        if (!why.empty()) return why;
      }
      for (ThmExperiment &x : thm.experiments)
        for (const auto &m : members.experiments)
          if (m.name == x.name) x.segments = m.datasetKeys;
      why = CheckThmExperiments(thm.experiments);
      if (!why.empty()) return "<thm> " + why;
      // Coulomb consistency of the experiments with the global options (its warnings: the caller's).
      why = CheckThmCoulombConsistency(thm.experiments, thm.coulombIntegral, &warnings);
      if (!why.empty()) return "<thm> " + why;
      if (!readFiles) return "";
      // ps=table:<file>, relative to the project file as the weight tables.
      for (ThmExperiment &x : thm.experiments) {
        if (x.psKind != ThmExperiment::PS_TABLE) continue;
        why = ReadThmPsTable(ThmRelativePath(projectFile, x.psTable), x.psTableP, x.psTableW);
        if (!why.empty()) return "<thm> experiment[" + x.name + "]: ps: " + why;
        x.psMin = x.psTableP.front();
        x.psMax = x.psTableP.back();
      }
      // spectatorAngles=[cm:]table:<file>
      for (ThmExperiment &x : thm.experiments) {
        if (x.angleWindow != 2) continue;
        why = ReadThmAngleTable(ThmRelativePath(projectFile, x.angleTable), x.angleTableT, x.angleTableW);
        if (!why.empty()) return "<thm> experiment[" + x.name + "]: spectatorAngles: " + why;
        x.angleMin = x.angleTableT.front();
        x.angleMax = x.angleTableT.back();
      }
      // distortion=table:<file>, the weight-table format
      for (ThmExperiment &x : thm.experiments) {
        if (x.distortion != ThmExperiment::DIST_TABLE) continue;
        std::shared_ptr<ThmWeightTable> table = std::make_shared<ThmWeightTable>();
        table->name = x.distortionTable;
        table->path = ThmRelativePath(projectFile, x.distortionTable);
        why = table->Read(table->path);
        if (!why.empty()) return "<thm> experiment[" + x.name + "]: distortion: " + why;
        x.distortionWeights = table;
      }
      return "";
    }
    if (trimmed.compare(0, 11, "experiment[") == 0) {
      // experiment[<name>] key=value ...: several keys on one line; an experiment no active data set names is left out
      if (std::find(idle.begin(), idle.end(), azure2::experiment::ReadThmSectionLine(trimmed).name) != idle.end()) continue;
      why = ParseThmExperimentLine(trimmed, thm.experiments);
      if (!why.empty()) return "<thm> " + why;
      continue;
    }
    size_t eq = trimmed.find('=');
    std::string key = trimmed.substr(0, eq);
    std::string value = eq == std::string::npos ? std::string() : trimmed.substr(eq + 1);
    key.erase(key.find_last_not_of(" \t") + 1);
    value.erase(0, value.find_first_not_of(" \t"));
    // The global keys: ParseThmOptionLine.
    ThmOptionLine o;
    if (eq == std::string::npos || !ParseThmOptionLine(key, value, o)) return "<thm> line not understood: '" + trimmed + "'";
    switch (o.kind) {
      case ThmOptionLine::VERTEX:
        thm.vertex = o.word == "onshell" ? ThmOptions::ON_SHELL : o.word == "constant" ? ThmOptions::CONSTANT : ThmOptions::PER_LEVEL;
        break;
      case ThmOptionLine::KINEMATICS:
        thm.kinematics = o.word == "lacognata" ? ThmOptions::LA_COGNATA
                         : o.word == "triple"  ? ThmOptions::TRIPLE
                         : o.word == "kf3body" ? ThmOptions::KF_THREE_BODY
                                               : ThmOptions::LAMBDA32;
        break;
      case ThmOptionLine::ENTRANCE_L:
        thm.coherentL = o.word == "coherent";
        break;
      case ThmOptionLine::COULOMB_INTEGRAL:
        thm.coulombIntegral = o.flag;
        break;
      case ThmOptionLine::SPECTATOR_ENERGY:
        if (key == "spectatorEnergy")
          thm.spectatorEnergy = o.energy;
        else
          thm.spectatorEnergyByPair[o.pairKey] = o.energy;
        break;
    }
  }
  return "<thm> block is not terminated by </thm>.";
}

/*!
 * Reads the <thm> section of the project: `key=value` lines setting the variant of the THM (HOES) observable (see
 * ThmOptions) and `experiment[<name>] key=value ...` lines; the members of each experiment and the weight tables are
 * the data sets' and extrapolations' records (fields thm, thm_weight; ThmMembershipSpec), read into `thm` here
 * (ReadFormat3ThmSection).  A '#' starts a comment.
 *
 *   <thm>
 *   vertex=onshell
 *   kinematics=triple
 *   entranceL=incoherent
 *   coulombIntegral=1
 *   spectatorEnergy=0.4
 *   spectatorEnergy[1]=0.6
 *   experiment[LC09] beam=15N target=d spectator=n Ebeam=60
 *   </thm>
 *
 * An unknown key is an error rather than a silent default, so a misspelt option cannot change a fit without notice.
 */
int Config::ReadThmBlock() {
  std::vector<std::string> lines;
  if (std::unique_ptr<std::istream> section = SectionStream("thm"))
    for (std::string line; std::getline(*section, line);) lines.push_back(line);
  std::vector<std::string> idle, warnings;
  const std::string why = ReadFormat3ThmSection(lines, project->built.thm, configfile, true, thm, idle, warnings);
  for (const std::string &name : idle)
    outStream << "WARNING: <thm> experiment[" << name << "] has no data set in use; ignored." << std::endl;
  for (const std::string &w : warnings) outStream << "WARNING: <thm> " << w << std::endl;
  if (!why.empty()) {
    outStream << "ERROR: " << why << std::endl;
    return -1;
  }
  return 0;
}

/*!
 * Reads a THM weight table: two columns, E_cm of the THM entrance pair (MeV)
 * and w(E) > 0, one row per line, '#' starts a comment.  At least two rows,
 * strictly increasing in E.  Returns "" on success, else what is wrong.
 */
std::string ThmWeightTable::Read(const std::string &file) {
  e.clear();
  w.clear();
  lnw.clear();
  std::ifstream in(file.c_str());
  if (!in) return "cannot read the weight file '" + file + "'";
  std::string line;
  int lineNumber = 0;
  while (std::getline(in, line)) {
    lineNumber++;
    size_t hash = line.find('#');
    if (hash != std::string::npos) line = line.substr(0, hash);
    if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;
    std::istringstream ls(line);
    double energy, weight;
    std::string extra;
    std::ostringstream where;
    where << "'" << file << "' line " << lineNumber << ": ";
    if (!(ls >> energy >> weight) || (ls >> extra))
      return where.str() + "expected two numbers, E (MeV) and w";
    if (!std::isfinite(energy) || !std::isfinite(weight) || !(weight > 0.0))
      return where.str() + "the weight must be finite and > 0";
    if (!e.empty() && !(energy > e.back()))
      return where.str() + "the energies must be strictly increasing";
    e.push_back(energy);
    w.push_back(weight);
    lnw.push_back(std::log(weight));
  }
  if (e.size() < 2) return "'" + file + "' needs at least two rows (E w)";
  return "";
}

bool SameThmWeight(const ThmWeightTable *a, const ThmWeightTable *b) {
  if (a == b) return true;
  return a && b && a->e == b->e && a->w == b->w && a->lnw == b->lnw;
}

/*!
 * w(E), linear in E and in ln w between the rows of the table; the end value
 * beyond either end (*outside is then set).  Two rows of equal w give that w
 * exactly, so a constant table scales the model by exactly w.
 */
double ThmWeightTable::operator()(double energy, bool *outside) const {
  if (outside) *outside = !Covers(energy);
  if (energy <= e.front()) return w.front();
  if (energy >= e.back()) return w.back();
  size_t hi = std::upper_bound(e.begin(), e.end(), energy) - e.begin();
  size_t lo = hi - 1;
  double t = (energy - e[lo]) / (e[hi] - e[lo]);
  return w[lo] * std::exp(t * (lnw[hi] - lnw[lo]));
}
