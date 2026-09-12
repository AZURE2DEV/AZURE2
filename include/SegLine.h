#ifndef SEGLINE_H
#define SEGLINE_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

/// A class to read and store a line from the data segments input file.

/*!
 * The SegLine class reads and stores a line from the data segments input file.
 */

/*!
 * Reads a data-file name that may be wrapped in double quotes.
 *
 * The name sits in the middle of the segment line: optional fields describing
 * composite and UPOS segments follow it, so the rest of the line cannot simply
 * be taken as the name, and a bare ">>" stops at the first space -- which
 * silently truncates any path containing one to its first word.  A quoted name
 * is therefore read as a single unit.  Unquoted names are returned exactly as
 * before, so every existing .azr file keeps parsing identically; only files
 * written with a space in the path use the quoted form.
 */
inline bool ReadPathToken(std::istream &stream, std::string &path) {
  std::string token;
  if (!(stream >> token)) return false;
  if (token.empty() || token[0] != '"') {
    path = token;
    return true;
  }
  path = token.substr(1);
  if (!path.empty() && path[path.size() - 1] == '"') {
    path.erase(path.size() - 1);  // the whole name fitted in one token
    return true;
  }
  // ">>" stopped at the space inside the name; take the remainder up to the
  // closing quote.  It still carries that space, so it is appended as read.
  std::string rest;
  if (std::getline(stream, rest, '"')) path += rest;
  return true;
}

/*!
 * Removes a matching pair of surrounding double quotes, if present.
 */
inline std::string StripPathQuotes(const std::string &path) {
  if (path.size() >= 2 && path[0] == '"' && path[path.size() - 1] == '"')
    return path.substr(1, path.size() - 2);
  return path;
}

class SegLine {
 public:
  /*!
   * Constructor fill the SegLine object from an input stream.
   */
  SegLine(std::istream &stream) {
    stream >> isActive_ >> entranceKey_ >> exitKey_ >> minE_ >> maxE_ >> minA_ >> maxA_ >> isDiff_;
    if (isDiff_ == 2) stream >> phaseJ_ >> phaseL_;
    stream >> dataNorm_ >> varyNorm_ >> dataNormError_;
    std::string dummyString;
    getline(stream, dummyString);

    // Initialize advanced segment parameters to defaults
    isAdvanced_ = 0;
    operationType_ = 0;
    componentsList_ = "";
    // Initialize UPOS parameters to defaults
    isUPOS_ = 0;
    secondaryDecayL_ = 0;
    Ic_ = 0.0;
    delta_ = 0.0;

    // Try to parse energy shift from the remaining line, default to 0.0 for backward compatibility
    std::istringstream restStream(dummyString);
    double tempEnergyShift;
    if (restStream >> tempEnergyShift) {
      energyShift_ = tempEnergyShift;
      // Try to read energyShiftError and varyEnergyShift
      if (restStream >> energyShiftError_ >> varyEnergyShift_) {
        // Successfully read all three values
      } else {
        // Only energy shift was available (old format)
        energyShiftError_ = 0.0;
        varyEnergyShift_ = 0;
      }
      // Get remaining string after energy shift parameters
      std::string remainingString;
      getline(restStream, remainingString);

      // Parse the file path and advanced segment data
      // Format: "filepath isAdvanced operationType numComponents entrance1 exit1 ..."
      std::istringstream advancedStream(remainingString);
      std::string tempDataFile;
      if (ReadPathToken(advancedStream, tempDataFile)) {
        dataFile_ = tempDataFile;

        // Try to parse advanced segment parameters
        int tempIsAdvanced;
        if (advancedStream >> tempIsAdvanced && tempIsAdvanced == 1) {
          isAdvanced_ = 1;
          if (advancedStream >> operationType_) {
            int numComponents;
            if (advancedStream >> numComponents) {
              // A marker value of -1 indicates the newer format that adds a
              // per-component "scaling" token. In that case the next integer
              // is the real component count and each component carries
              // (entrance, exit, angle, scaling) instead of (entrance, exit, angle).
              bool hasScalingToken = false;
              if (numComponents == -1) {
                hasScalingToken = true;
                if (!(advancedStream >> numComponents)) numComponents = 0;
              }
              std::vector<std::string> components;
              for (int i = 0; i < numComponents; i++) {
                int entrance, exit;
                double angle;
                if (advancedStream >> entrance >> exit >> angle) {
                  double scaling = 1.0;
                  if (hasScalingToken) {
                    if (!(advancedStream >> scaling)) scaling = 1.0;
                  }
                  std::string componentStr = "Entrance: " + std::to_string(entrance) + ", Exit: " + std::to_string(exit);
                  if (angle > -900.0) {
                    componentStr += ", Angle: " + std::to_string(angle);
                  }
                  if (hasScalingToken && scaling > -900.0) {
                    componentStr += ", Scaling: " + std::to_string(scaling);
                  }
                  components.push_back(componentStr);
                } else {
                  // Fallback to old format (entrance, exit only) for backward compatibility
                  advancedStream.clear();
                  if (advancedStream >> entrance >> exit) {
                    components.push_back("Entrance: " + std::to_string(entrance) + ", Exit: " + std::to_string(exit));
                  }
                }
              }
              if (!components.empty()) {
                componentsList_ = components[0];
                for (size_t i = 1; i < components.size(); i++) {
                  componentsList_ += ";" + components[i];
                }
              }
            }
          }
        }
        // Try to parse optional UPOS parameters after file path and advanced data
        int tempIsUPOS;
        if (advancedStream >> tempIsUPOS) {
          isUPOS_ = tempIsUPOS;
          if (advancedStream >> secondaryDecayL_ >> Ic_ >> delta_) {
            // successfully read all UPOS fields
          } else {
            secondaryDecayL_ = 0;
            Ic_ = 0.0;
            delta_ = 0.0;
          }
        }
      } else {
        // If no file path found, use trimmed string
        int p2 = remainingString.find_last_not_of(" \n\t\r");
        if (p2 != std::string::npos) {
          int p1 = remainingString.find_first_not_of(" \n\t\r");
          if (p1 == std::string::npos) p1 = 0;
          dataFile_ = StripPathQuotes(remainingString.substr(p1, (p2 - p1) + 1));
        } else
          dataFile_ = std::string();
      }
    } else {
      // No energy shift found, treat entire dummyString as dataFile (backward compatibility)
      energyShift_ = 0.0;
      energyShiftError_ = 0.0;
      varyEnergyShift_ = 0;
      int p2 = dummyString.find_last_not_of(" \n\t\r");
      if (p2 != std::string::npos) {
        int p1 = dummyString.find_first_not_of(" \n\t\r");
        if (p1 == std::string::npos) p1 = 0;
        dataFile_ = StripPathQuotes(dummyString.substr(p1, (p2 - p1) + 1));
      } else
        dataFile_ = std::string();
    }
  };
  /*!
   * Returns non-zero if the line is to be included in the calculation.
   */
  int isActive() const { return isActive_; };
  /*!
   * Returns the particle pair key corresponding to the
   * entrance channel for the data segment.
   */
  int entranceKey() const { return entranceKey_; };
  /*!
   * Returns the particle pair key corresponding to the
   * exit channel for the data segment.
   */
  int exitKey() const { return exitKey_; };
  /*!
   * Returns the minimum energy to be included in the segment from the data.
   */
  double minE() const { return minE_; };
  /*!
   * Returns the maximum energy to be included in the segment from the data.
   */
  double maxE() const { return maxE_; };
  /*!
   * Returns the minimum angle to be included in the segment from the data.
   */
  double minA() const { return minA_; };
  /*!
   * Returns the maximum angle to be included in the segment from the data.
   */
  double maxA() const { return maxA_; };
  /*!
   * Return 0 if the segment is angle-integrated cross section, 1 for
   * differential cross section, and 2 for phase shift.
   */
  int isDiff() const { return isDiff_; };
  /*!
   * Returns the path of the data file for the segment.
   */
  std::string dataFile() const { return dataFile_; };
  /*!
   * Returns the data normalization for the segment.
   */
  double dataNorm() const { return dataNorm_; };
  /*!
   * Returns the data normalization error for the segment.
   */
  double dataNormError() const { return dataNormError_; };
  /*!
   * Returns non-zero of the normalization is to be fit.
   */
  int varyNorm() const { return varyNorm_; };
  /*!
   * Returns the spin value for the segment if the segment contains phase shift.
   */
  double phaseJ() const { return phaseJ_; };
  /*!
   * Returns the orbital angular momentum value for the segment
   * if the segment contains phase shift.
   */
  int phaseL() const { return phaseL_; };
  /*!
   * Returns the energy shift value for the segment.
   */
  double energyShift() const { return energyShift_; };
  /*!
   * Returns the energy shift error for the segment.
   */
  double energyShiftError() const { return energyShiftError_; };
  /*!
   * Returns non-zero if the energy shift is to be fit.
   */
  int varyEnergyShift() const { return varyEnergyShift_; };
  /*!
   * Returns non-zero if this is an advanced segment (sum/ratio).
   */
  int isAdvanced() const { return isAdvanced_; };
  /*!
   * Returns the operation type (0 for sum, 1 for ratio).
   */
  int operationType() const { return operationType_; };
  /*!
   * Returns the semicolon-separated list of components.
   */
  std::string componentsList() const { return componentsList_; };
  /*!
   * Returns flag for if this segment is an unobserved primary transition (1 = is, 0 = not).
   */
  int isUPOS() const { return isUPOS_; };
  /*!
   * Returns the angular momentum of the decay for unobserved primary, observed secondary.
   */
  int secondaryDecayL() const { return secondaryDecayL_; };
  /*!
   * Returns the spin of the final state of the decay for unobserved primary, observed secondary.
   */
  double Ic() const { return Ic_; };
  /*!
   * Returns the multipole mixing ratio of the decay for unobserved primary, observed secondary.
   */
  double delta() const { return delta_; };

 private:
  int isActive_;
  int entranceKey_;
  int exitKey_;
  double minE_;
  double maxE_;
  double minA_;
  double maxA_;
  int isDiff_;
  std::string dataFile_;
  double dataNorm_;
  double dataNormError_;
  int varyNorm_;
  double phaseJ_;
  int phaseL_;
  double energyShift_;
  double energyShiftError_;
  int varyEnergyShift_;
  int isAdvanced_;
  int operationType_;
  std::string componentsList_;
  int isUPOS_;
  int secondaryDecayL_;
  double Ic_;
  double delta_;
};

#endif
