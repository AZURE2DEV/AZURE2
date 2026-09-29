#include "Config.h"
#include "NuclearPotentialManager.h"
#ifndef NO_STAT
#include <sys/stat.h>
#endif
#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

/*!
 * The constructor of the Config class sets defaults and the
 * stream reference for output.
 */

Config::Config(std::ostream &stream) :
  outStream(stream) {
  Reset();
}

/*!
 * This function resets Config structure.
 */

void Config::Reset() {
  chiVariance = 1.0;
  screenCheckMask = 0;
  fileCheckMask = 0;
  paramMask = 0;
  paramMask |= (USE_AMATRIX | USE_BRUNE_FORMALISM | IGNORE_ZERO_WIDTHS | TRANSFORM_PARAMETERS | CALCULATE_WITH_DATA | USE_LONGWAVELENGTH_APPROX);
  stopFlag = false;
  outputdir = "";
  checkdir = "";
  nloptAlgorithm = 0;       // Default to SBPLX
  useHybridMethod = false;  // Default to disabled
  useAdaptiveGrid = true;   // Default to adaptive grid
  thm = ThmOptions();
}

/*!
 * This funciton reads the configuration file and parses various options.
 */

int Config::ReadConfigFile() {
  std::string dummy;
  std::string temp;
  std::ifstream in(configfile.c_str());
  if (!in) return -1;
  std::string line = "";
  while (line != "<config>" && !in.eof()) getline(in, line);
  if (line != "<config>") return -1;
  in >> temp;
  getline(in, dummy);
  if (temp == "true")
    paramMask |= USE_AMATRIX;
  else
    paramMask &= ~USE_AMATRIX;
  getline(in, dummy);
  int poundSignPos = dummy.find_last_of('#');
  if (poundSignPos == std::string::npos)
    temp = dummy;
  else
    temp = dummy.substr(0, poundSignPos);
  int p2 = temp.find_last_not_of(" \n\t\r");
  if (p2 != std::string::npos) {
    int p1 = temp.find_first_not_of(" \n\t\r");
    if (p1 == std::string::npos) p1 = 0;
    outputdir = temp.substr(p1, (p2 - p1) + 1);
  } else
    outputdir = std::string();
  getline(in, dummy);
  poundSignPos = dummy.find_last_of('#');
  if (poundSignPos == std::string::npos)
    temp = dummy;
  else
    temp = dummy.substr(0, poundSignPos);
  p2 = temp.find_last_not_of(" \n\t\r");
  if (p2 != std::string::npos) {
    int p1 = temp.find_first_not_of(" \n\t\r");
    if (p1 == std::string::npos) p1 = 0;
    checkdir = temp.substr(p1, (p2 - p1) + 1);
  } else
    checkdir = std::string();
  in >> temp;
  getline(in, dummy);
  if (temp == "screen")
    screenCheckMask |= CHECK_COMPOUND_NUCLEUS;
  else if (temp == "file")
    fileCheckMask |= CHECK_COMPOUND_NUCLEUS;
  in >> temp;
  getline(in, dummy);
  if (temp == "screen")
    screenCheckMask |= CHECK_BOUNDARY_CONDITIONS;
  else if (temp == "file")
    fileCheckMask |= CHECK_BOUNDARY_CONDITIONS;
  in >> temp;
  getline(in, dummy);
  if (temp == "screen")
    screenCheckMask |= CHECK_DATA;
  else if (temp == "file")
    fileCheckMask |= CHECK_DATA;
  in >> temp;
  getline(in, dummy);
  if (temp == "screen")
    screenCheckMask |= CHECK_ENERGY_DEP;
  else if (temp == "file")
    fileCheckMask |= CHECK_ENERGY_DEP;
  in >> temp;
  getline(in, dummy);
  if (temp == "screen")
    screenCheckMask |= CHECK_LEGENDRE;
  else if (temp == "file")
    fileCheckMask |= CHECK_LEGENDRE;
  in >> temp;
  getline(in, dummy);
  if (temp == "screen")
    screenCheckMask |= CHECK_COUL_AMPLITUDES;
  else if (temp == "file")
    fileCheckMask |= CHECK_COUL_AMPLITUDES;
  in >> temp;
  getline(in, dummy);
  if (temp == "screen")
    screenCheckMask |= CHECK_PATHWAYS;
  else if (temp == "file")
    fileCheckMask |= CHECK_PATHWAYS;
  in >> temp;
  getline(in, dummy);
  if (temp == "screen")
    screenCheckMask |= CHECK_ANGULAR_DISTS;
  else if (temp == "file")
    fileCheckMask |= CHECK_ANGULAR_DISTS;
  line = "";
  while (line != "</config>" && !in.eof()) getline(in, line);
  if (line != "</config>") return -1;
  in.close();
  // A malformed optional block has already been reported: -2, so the caller
  // does not add a misleading "could not open" on top.
  if (this->ReadPotentialBlock() != 0) return -2;
  return this->ReadThmBlock() != 0 ? -2 : 0;
}

/*!
 * Reads the <potential> block of the configuration file and configures the
 * hybrid Coulomb method accordingly.
 *
 * The block is optional -- a file without one, or one with
 * useHybridPotential=0, leaves the defaults set by Reset() untouched, so this
 * is a no-op for every existing project.  The format is the one written by the
 * setup utility (gui/src/AZURESetup.cpp) and read back by
 * NuclearPotentialTab::readPotentialSettings:
 *
 *   <potential>
 *   useHybridPotential=1
 *   useAdaptiveGrid=1
 *   potentialType=0        # 0 = Woods-Saxon, 1 = Gaussian
 *   V0=80                  # depth, MeV
 *   R=3.6                  # radius, fm      (Woods-Saxon)
 *   a=0.6                  # diffuseness, fm (Woods-Saxon)
 *   r0=5.0                 # width, fm       (Gaussian)
 *   pair=2                 # everything below applies to pair 2 alone
 *   potentialType=1
 *   V0=90
 *   r0=4.0
 *   pair=5
 *   useHybridPotential=0   # pair 5 opts out
 *   </potential>
 *
 * A nuclear potential belongs to a particle pair -- it bends the radial wave
 * functions of that channel and no other -- so the keys are read per pair.
 * Everything before the first pair= is the default, which stands in for every
 * pair that does not name itself; a pair= line opens a section that starts
 * from that default and only has to state what differs.  A file with no pair=
 * line is therefore read exactly as it was when the model was global.
 *
 * useAdaptiveGrid is a property of the target-effect integration rather than
 * of any one channel, so it stays global wherever it appears.
 *
 * Parsing it here rather than in the setup utility is what makes the hybrid
 * model reachable from --no-gui and from pyazr: this function is on the path
 * both of them take.
 *
 * Returns 0 on success (including "no block present") and -1 if the block is
 * present but malformed.
 */

int Config::ReadPotentialBlock() {
  std::ifstream in(configfile.c_str());
  if (!in) return -1;

  std::string line = "";
  while (line != "<potential>" && !in.eof()) getline(in, line);
  if (line != "<potential>") return 0;  // optional block, absent

  NuclearPotentialManager &manager = NuclearPotentialManager::instance();
  manager.resetToDefault();

  NuclearPotentialSetting current = manager.getDefaultSetting();
  int currentPair = 0;  // 0 = the default section
  bool currentHasType = false, defaultHasType = false, closed = false;

  // Install whatever the section just read describes.  A section that asks for
  // the hybrid model without naming a potentialType, or names one that does not
  // exist, is refused rather than silently given the default shape.
  auto flush = [&]() {
    if (current.enabled && !currentHasType) {
      outStream << "WARNING: <potential> requests the hybrid method";
      if (currentPair) outStream << " for pair " << currentPair;
      outStream << " but gives no potentialType; it is disabled." << std::endl;
      current.enabled = false;
    }
    try {
      if (currentPair)
        manager.setSetting(currentPair, current);
      else
        manager.setDefaultSetting(current);
    } catch (const std::exception &e) {
      outStream << "WARNING: invalid potential parameters (" << e.what()
                << "); the hybrid method is disabled";
      if (currentPair) outStream << " for pair " << currentPair;
      outStream << "." << std::endl;
      current.enabled = false;
      if (currentPair)
        manager.clearPairSetting(currentPair);
      else
        manager.setDefaultEnabled(false);
    }
  };

  while (!in.eof()) {
    getline(in, line);
    size_t b = line.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) continue;
    size_t e = line.find_last_not_of(" \t\r\n");
    std::string trimmed = line.substr(b, e - b + 1);
    if (trimmed == "</potential>") {
      closed = true;
      break;
    }

    size_t eq = trimmed.find('=');
    if (eq == std::string::npos) continue;
    std::string key = trimmed.substr(0, eq);
    std::istringstream value(trimmed.substr(eq + 1));

    // A pair= line closes the section being read and opens the next one.  Every
    // section starts from the default, so a pair only has to state what differs.
    if (key == "pair") {
      int nextPair = 0;
      value >> nextPair;
      if (nextPair <= 0) {
        outStream << "WARNING: <potential> has pair=" << nextPair
                  << ", which is not a pair key; the section is ignored." << std::endl;
        continue;
      }
      flush();
      currentPair = nextPair;
      current = manager.getDefaultSetting();
      currentHasType = defaultHasType;  // inherited with the shape
      continue;
    }

    if (key == "useHybridPotential") {
      int v = 0;
      value >> v;
      current.enabled = (v == 1);
    } else if (key == "useAdaptiveGrid") {
      int v = 1;
      value >> v;
      useAdaptiveGrid = (v == 1);
    } else if (key == "potentialType") {
      int typeCode = 0;
      value >> typeCode;
      if (typeCode == 0) {
        current.type = "WoodsSaxon";
        currentHasType = true;
      } else if (typeCode == 1) {
        current.type = "Gaussian";
        currentHasType = true;
      } else {
        outStream << "WARNING: unknown potentialType " << typeCode
                  << " in <potential>";
        if (currentPair) outStream << " for pair " << currentPair;
        outStream << "; the hybrid method is disabled there." << std::endl;
        current.enabled = false;
        currentHasType = true;  // refused, not merely missing
      }
      if (!currentPair) defaultHasType = currentHasType;
    } else if (key == "V0")
      value >> current.V0;
    else if (key == "R")
      value >> current.R;
    else if (key == "a")
      value >> current.a;
    else if (key == "r0")
      value >> current.r0;
  }
  in.close();

  if (!closed) return -1;  // unterminated block
  flush();

  // The global switch stays the master: it gates the paramMask bit and the
  // GUI tab, while the manager decides which pairs the model applies to.
  useHybridMethod = manager.isAnyEnabled();
  return 0;
}

/*!
 * If stat() is enabled, this function checks for the output and checks
 * directories at runtime.
 */

#ifndef NO_STAT
int Config::CheckForInputFiles() {
  struct stat buffer;
  if (stat(outputdir.c_str(), &buffer) != 0) {
    outStream << "Could not find output directory: " << outputdir << ". Check that it exists." << std::endl;
    return -1;
  }
  if (stat(checkdir.c_str(), &buffer) != 0) {
    outStream << "Could not find checks directory: " << checkdir << ". Check that it exists." << std::endl;
    return -1;
  }
  return 0;
}
#endif

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
      ok = !!(vs >> x) && x >= 0.0;
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
