/*
 * The <thm> block of the configuration (Config::ReadThmBlock) and the weight
 * tables it names: the THM part of Config, kept out of Config.cpp so that the
 * classic configuration stays as it is on dev.
 */
#include "Config.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <sstream>

/*!
 * Reads the optional <thm> block: `key=value` lines setting the variant of the
 * THM (HOES) observable, see Config::ThmOptions.  A '#' starts a comment.
 *
 *   <thm>
 *   vertex=onshell
 *   kinematics=triple
 *   entranceL=incoherent
 *   coulombIntegral=1
 *   spectatorEnergy=0.4
 *   spectatorEnergy[1]=0.6
 *   </thm>
 *
 * An unknown key is an error rather than a silent default, so a misspelt
 * option cannot change a fit without notice.
 */

int Config::ReadThmBlock() {
  std::ifstream in(configfile.c_str());
  if (!in) return -1;
  std::string line = "";
  while (!in.eof()) {
    getline(in, line);
    size_t b = line.find_first_not_of(" \t\r\n");
    if (b != std::string::npos && line.compare(b, 5, "<thm>") == 0) break;
    line = "";
  }
  if (line.find("<thm>") == std::string::npos) return 0;  // optional block, absent

  // A file named in the block, relative to the .azr.  Absolute is a leading
  // slash or backslash, or (Windows) a drive letter: "C:/..." or "C:\..."
  // handed to the native binary was otherwise prefixed with the project
  // directory.
  auto thmRelativePath = [this](const std::string &value) {
    std::string dir;
    size_t slash = configfile.find_last_of("/\\");
    if (slash != std::string::npos) dir = configfile.substr(0, slash + 1);
    bool absolute = !value.empty() && (value[0] == '/' || value[0] == '\\' ||
                                       (value.size() > 1 && value[1] == ':' && std::isalpha((unsigned char)value[0])));
    return (absolute || dir.empty()) ? value : dir + value;
  };
  auto flag = [](const std::string &v, bool &out) {
    if (v == "1" || v == "true" || v == "on") out = true;
    else if (v == "0" || v == "false" || v == "off") out = false;
    else return false;
    return true;
  };
  while (!in.eof()) {
    getline(in, line);
    size_t hash = line.find('#');
    if (hash != std::string::npos) line = line.substr(0, hash);
    size_t b = line.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) continue;
    size_t e = line.find_last_not_of(" \t\r\n");
    std::string trimmed = line.substr(b, e - b + 1);
    if (trimmed == "</thm>") {
      std::string why = CheckThmExperiments(thm.experiments);
      if (!why.empty()) {
        outStream << "ERROR: <thm> " << why << std::endl;
        return -1;
      }
      // Coulomb consistency of the experiments with the global options.
      std::vector<std::string> warnings;
      why = CheckThmCoulombConsistency(thm.experiments, thm.coulombIntegral, &warnings);
      if (!why.empty()) {
        outStream << "ERROR: <thm> " << why << std::endl;
        return -1;
      }
      for (const std::string &w : warnings) outStream << "WARNING: <thm> " << w << std::endl;
      // ps=table:<file>, relative to the .azr as weight[k]=.
      for (ThmExperiment &x : thm.experiments) {
        if (x.psKind != ThmExperiment::PS_TABLE) continue;
        why = ReadThmPsTable(thmRelativePath(x.psTable), x.psTableP, x.psTableW);
        if (!why.empty()) {
          outStream << "ERROR: <thm> experiment[" << x.name << "]: ps: " << why << std::endl;
          return -1;
        }
        x.psMin = x.psTableP.front();
        x.psMax = x.psTableP.back();
      }
      // spectatorAngles=[cm:]table:<file>, relative to the .azr as weight[k]=.
      for (ThmExperiment &x : thm.experiments) {
        if (x.angleWindow != 2) continue;
        why = ReadThmAngleTable(thmRelativePath(x.angleTable), x.angleTableT, x.angleTableW);
        if (!why.empty()) {
          outStream << "ERROR: <thm> experiment[" << x.name << "]: spectatorAngles: " << why << std::endl;
          return -1;
        }
        x.angleMin = x.angleTableT.front();
        x.angleMax = x.angleTableT.back();
      }
      // distortion=table:<file>, the weight[k] format, relative to the .azr.
      for (ThmExperiment &x : thm.experiments) {
        if (x.distortion != ThmExperiment::DIST_TABLE) continue;
        std::shared_ptr<ThmWeightTable> table = std::make_shared<ThmWeightTable>();
        table->name = x.distortionTable;
        table->path = thmRelativePath(x.distortionTable);
        why = table->Read(table->path);
        if (!why.empty()) {
          outStream << "ERROR: <thm> experiment[" << x.name << "]: distortion: " << why << std::endl;
          return -1;
        }
        x.distortionWeights = table;
      }
      return 0;
    }
    if (trimmed.compare(0, 11, "experiment[") == 0) {
      // experiment[<name>] key=value ...: several keys on one line.
      std::string why = ParseThmExperimentLine(trimmed, thm.experiments);
      if (!why.empty()) {
        outStream << "ERROR: <thm> " << why << std::endl;
        return -1;
      }
      continue;
    }
    size_t eq = trimmed.find('=');
    std::string key = trimmed.substr(0, eq);
    std::string value = eq == std::string::npos ? std::string() : trimmed.substr(eq + 1);
    key.erase(key.find_last_not_of(" \t") + 1);
    value.erase(0, value.find_first_not_of(" \t"));
    bool ok = eq != std::string::npos;
    if (!ok) {
    } else if (key == "vertex") {
      if (value == "onshell") thm.vertex = ThmOptions::ON_SHELL;
      else if (value == "constant") thm.vertex = ThmOptions::CONSTANT;
      else if (value == "perlevel" || value == "real") thm.vertex = ThmOptions::PER_LEVEL;
      else ok = false;
    } else if (key == "kinematics") {
      if (value == "lacognata") thm.kinematics = ThmOptions::LA_COGNATA;
      else if (value == "triple") thm.kinematics = ThmOptions::TRIPLE;
      else if (value == "kf3body") thm.kinematics = ThmOptions::KF_THREE_BODY;
      else if (value == "lambda32") thm.kinematics = ThmOptions::LAMBDA32;
      else ok = false;
    } else if (key == "entranceL") {
      if (value == "coherent") thm.coherentL = true;
      else if (value == "incoherent") thm.coherentL = false;
      else ok = false;
    } else if (key == "coulombIntegral") {
      ok = flag(value, thm.coulombIntegral);
    } else if (key.compare(0, 15, "spectatorEnergy") == 0) {
      std::istringstream vs(value);
      double x;
      std::string rest;
      ok = !!(vs >> x) && x >= 0.0 && !(vs >> rest);  // "0.4junk" is not a number
      if (ok && key == "spectatorEnergy") thm.spectatorEnergy = x;
      else if (ok && key.size() > 17 && key[15] == '[' && key.back() == ']') {
        std::istringstream ks(key.substr(16, key.size() - 17));
        int pairKey;
        ok = !!(ks >> pairKey);
        if (ok) thm.spectatorEnergyByPair[pairKey] = x;
      } else ok = false;
    } else if (key.compare(0, 6, "weight") == 0) {
      // weight[<k>]=<file> (k-th <segmentsData> line) or weightTest[<k>]=<file>
      // (k-th <segmentsTest> line).  Which segments exist, and whether they are
      // THM, is checked once the data are read (EData::Fill / MakePoints).
      bool test = key.compare(0, 11, "weightTest[") == 0;
      size_t open = test ? 10 : 6;
      int segKey = 0;
      ok = key.size() > open + 2 && key[open] == '[' && key.back() == ']' && !value.empty();
      if (ok) {
        std::istringstream ks(key.substr(open + 1, key.size() - open - 2));
        std::string rest;
        ok = !!(ks >> segKey) && segKey >= 1 && !(ks >> rest);
      }
      if (ok) {
        std::shared_ptr<ThmWeightTable> table = std::make_shared<ThmWeightTable>();
        table->name = value;
        table->path = thmRelativePath(value);
        std::string why = table->Read(table->path);
        if (!why.empty()) {
          outStream << "ERROR: <thm> " << key << ": " << why << std::endl;
          return -1;
        }
        (test ? thm.weightByTestSegment : thm.weightBySegment)[segKey] = table;
      }
    } else ok = false;
    if (!ok) {
      outStream << "ERROR: <thm> line not understood: '" << trimmed << "'" << std::endl;
      return -1;
    }
  }
  outStream << "ERROR: <thm> block is not terminated by </thm>." << std::endl;
  return -1;
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
