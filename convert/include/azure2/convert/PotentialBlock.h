#ifndef AZURE2_CONVERT_POTENTIALBLOCK_H
#define AZURE2_CONVERT_POTENTIALBLOCK_H

#include <istream>
#include <sstream>
#include <string>
#include <vector>

#include "azure2/experiment/HybridPotentialSettings.h"

namespace azure2 {
namespace convert {

/*!
 * The old <potential> block (formats 1 and 2), read by the engine (Config::ReadPotentialBlock) and by the converter to
 * format 3 with this one parser:
 *
 *   <potential>
 *   useHybridPotential=1
 *   useAdaptiveGrid=1      # format 1; format 2 has it in <config>
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
 * Everything before the first pair= is the default, which stands in for every pair that does not name itself; a pair=
 * line opens a section that starts from that default and states only what differs.  A section that asks for the
 * hybrid model without naming a potentialType, or names one that does not exist, is disabled with a warning rather
 * than given the default shape; a section whose parameters do not make a potential (HybridPotentialSetting::Potential
 * throws) is not installed: a pair falls back to the default, the default is off.  The block has no Coulomb radius,
 * matching radius or step: the settings keep 3.5 fm, 15 fm and 0.0275 fm.
 */
struct PotentialBlock {
  experiment::HybridPotentialSettings settings;
  bool hasAdaptiveGrid = false;  ///< a useAdaptiveGrid= line was read (the last one wins)
  bool adaptiveGrid = true;      ///< its value (== 1)
  bool closed = false;           ///< </potential> was read; without it the last section is not installed
  /// The warnings, in the order the lines that cause them were read, each a whole "WARNING: ..." line.
  std::vector<std::string> warnings;
};

/// Reads the block from `in`, which holds it from its <potential> line (that first line is skipped).
inline PotentialBlock ReadPotentialBlock(std::istream &in) {
  PotentialBlock out;
  experiment::HybridPotentialSettings &potentials = out.settings;
  std::string line = "";
  getline(in, line);  // <potential>

  experiment::HybridPotentialSetting current = potentials.defaults;
  int currentPair = 0;  // 0 = the default section
  bool currentHasType = false, defaultHasType = false;
  auto forPair = [&]() { return currentPair ? " for pair " + std::to_string(currentPair) : std::string(); };

  // Install whatever the section just read describes.
  auto flush = [&]() {
    if (current.enabled && !currentHasType) {
      out.warnings.push_back("WARNING: <potential> requests the hybrid method" + forPair() +
                             " but gives no potentialType; it is disabled.");
      current.enabled = false;
    }
    try {
      current.Potential();
      if (currentPair)
        potentials.pairs[currentPair] = current;
      else
        potentials.defaults = current;
    } catch (const std::exception &e) {
      out.warnings.push_back(std::string("WARNING: invalid potential parameters (") +
                             (currentPair ? "Failed to set potential for pair: " : "Failed to set default potential: ") +
                             e.what() + "); the hybrid method is disabled" + forPair() + ".");
      if (currentPair)
        potentials.pairs.erase(currentPair);
      else
        potentials.defaults.enabled = false;
    }
  };

  while (!in.eof()) {
    getline(in, line);
    size_t b = line.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) continue;
    size_t e = line.find_last_not_of(" \t\r\n");
    std::string trimmed = line.substr(b, e - b + 1);
    if (trimmed == "</potential>") {
      out.closed = true;
      break;
    }

    size_t eq = trimmed.find('=');
    if (eq == std::string::npos) continue;
    std::string key = trimmed.substr(0, eq);
    std::istringstream value(trimmed.substr(eq + 1));

    // A pair= line closes the section being read and opens the next one.  Every section starts from the default, so
    // a pair only has to state what differs.
    if (key == "pair") {
      int nextPair = 0;
      value >> nextPair;
      if (nextPair <= 0) {
        out.warnings.push_back("WARNING: <potential> has pair=" + std::to_string(nextPair) +
                               ", which is not a pair key; the section is ignored.");
        continue;
      }
      flush();
      currentPair = nextPair;
      current = potentials.defaults;
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
      out.hasAdaptiveGrid = true;
      out.adaptiveGrid = (v == 1);
    } else if (key == "potentialType") {
      int typeCode = 0;
      value >> typeCode;
      if (typeCode == 0) {
        current.shape = physics::FormShape::WoodsSaxon;
        currentHasType = true;
      } else if (typeCode == 1) {
        current.shape = physics::FormShape::Gaussian;
        currentHasType = true;
      } else {
        out.warnings.push_back("WARNING: unknown potentialType " + std::to_string(typeCode) + " in <potential>" + forPair() +
                               "; the hybrid method is disabled there.");
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

  if (out.closed) flush();
  return out;
}

}  // namespace convert
}  // namespace azure2

#endif
