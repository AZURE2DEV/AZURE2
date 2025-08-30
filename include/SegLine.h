#ifndef SEGLINE_H
#define SEGLINE_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

///A class to read and store a line from the data segments input file.

/*!
 * The SegLine class reads and stores a line from the data segments input file.
 */

class SegLine {
 public:
  /*!
   * Constructor fill the SegLine object from an input stream.
   */
  SegLine(std::istream &stream) {
    stream >> isActive_ >>  entranceKey_ >> exitKey_ >> minE_ >> maxE_ >> minA_ >> maxA_ >> isDiff_;
    if(isDiff_==2) stream >> phaseJ_ >> phaseL_;
    stream >> dataNorm_ >> varyNorm_ >> dataNormError_;
    std::string dummyString;
    getline(stream,dummyString);
    
    // Initialize advanced segment parameters to defaults
    isAdvanced_ = 0;
    operationType_ = 0;
    componentsList_ = "";
    
    // Try to parse energy shift from the remaining line, default to 0.0 for backward compatibility
    std::istringstream restStream(dummyString);
    double tempEnergyShift;
    if(restStream >> tempEnergyShift) {
      energyShift_ = tempEnergyShift;
      // Try to read energyShiftError and varyEnergyShift
      if(restStream >> energyShiftError_ >> varyEnergyShift_) {
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
      if(advancedStream >> tempDataFile) {
        dataFile_ = tempDataFile;
        
        // Try to parse advanced segment parameters
        int tempIsAdvanced;
        if(advancedStream >> tempIsAdvanced && tempIsAdvanced == 1) {
          isAdvanced_ = 1;
          if(advancedStream >> operationType_) {
            int numComponents;
            if(advancedStream >> numComponents) {
              std::vector<std::string> components;
              for(int i = 0; i < numComponents; i++) {
                int entrance, exit;
                if(advancedStream >> entrance >> exit) {
                  components.push_back("Entrance: " + std::to_string(entrance) + ", Exit: " + std::to_string(exit));
                }
              }
              if(!components.empty()) {
                componentsList_ = components[0];
                for(size_t i = 1; i < components.size(); i++) {
                  componentsList_ += ";" + components[i];
                }
              }
            }
          }
        }
      } else {
        // If no file path found, use trimmed string
        int p2 = remainingString.find_last_not_of(" \n\t\r");
        if (p2 != std::string::npos) {  
          int p1 = remainingString.find_first_not_of(" \n\t\r");
          if (p1 == std::string::npos) p1 = 0;
          dataFile_=remainingString.substr(p1,(p2-p1)+1);
        } else dataFile_=std::string();  
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
        dataFile_=dummyString.substr(p1,(p2-p1)+1);
      } else dataFile_=std::string();  
    }
  };
  /*!
   * Returns non-zero if the line is to be included in the calculation.
   */
  int isActive() const {return isActive_;};
  /*!
   * Returns the particle pair key corresponding to the 
   * entrance channel for the data segment.
   */
  int entranceKey() const {return entranceKey_;};
  /*!
   * Returns the particle pair key corresponding to the 
   * exit channel for the data segment.
   */
  int exitKey() const {return exitKey_;};
  /*!
   * Returns the minimum energy to be included in the segment from the data.
   */
  double minE() const {return minE_;};
  /*!
   * Returns the maximum energy to be included in the segment from the data.
   */
  double maxE() const {return maxE_;};
  /*!
   * Returns the minimum angle to be included in the segment from the data.
   */
  double minA() const {return minA_;};
  /*!
   * Returns the maximum angle to be included in the segment from the data.
   */
  double maxA() const {return maxA_;};
  /*!
   * Return 0 if the segment is angle-integrated cross section, 1 for 
   * differential cross section, and 2 for phase shift.
   */
  int isDiff() const {return isDiff_;};
  /*!
   * Returns the path of the data file for the segment.
   */
  std::string dataFile() const {return dataFile_;};
  /*!
   * Returns the data normalization for the segment.
   */
  double dataNorm() const {return dataNorm_;};
  /*!
   * Returns the data normalization error for the segment.
   */
  double dataNormError() const {return dataNormError_;};
  /*! 
   * Returns non-zero of the normalization is to be fit.
   */
  int varyNorm() const {return varyNorm_;};
  /*!
   * Returns the spin value for the segment if the segment contains phase shift.
   */
  double phaseJ() const {return phaseJ_;};
  /*!
   * Returns the orbital angular momentum value for the segment
   * if the segment contains phase shift.
   */
  int phaseL() const {return phaseL_;};
  /*!
   * Returns the energy shift value for the segment.
   */
  double energyShift() const {return energyShift_;};
  /*!
   * Returns the energy shift error for the segment.
   */
  double energyShiftError() const {return energyShiftError_;};
  /*!
   * Returns non-zero if the energy shift is to be fit.
   */
  int varyEnergyShift() const {return varyEnergyShift_;};
  /*!
   * Returns non-zero if this is an advanced segment (sum/ratio).
   */
  int isAdvanced() const {return isAdvanced_;};
  /*!
   * Returns the operation type (0 for sum, 1 for ratio).
   */
  int operationType() const {return operationType_;};
  /*!
   * Returns the semicolon-separated list of components.
   */
  std::string componentsList() const {return componentsList_;};
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
};

#endif
