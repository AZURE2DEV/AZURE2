#ifndef EXTRAPLINE_H
#define EXTRAPLINE_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

///A class to read and store a line from the extrapolation input file.

/*!
 * The ExtrapLine class reads and stores a line from the extrapolation input file.
 */

class ExtrapLine {
 public:
  /*!
   * Constructor fill the ExtrapLine object from an input stream.
   */
  ExtrapLine(std::istream &stream) {
    stream >> isActive_ >> entranceKey_ >> exitKey_ >> minE_
	   >> maxE_ >> eStep_ >> minA_ >> maxA_ >> aStep_ >> isDiff_;
    if(isDiff_==2) stream >> phaseJ_ >> phaseL_;
    else if(isDiff_==3) stream >> maxAngDistOrder_;

    // Initialize advanced segment parameters to defaults
    isAdvanced_ = 0;
    operationType_ = 0;
    componentsList_ = "";

    // Try to read advanced segment data from the remaining line
    std::string dummyString;
    getline(stream, dummyString);

    if(!dummyString.empty()) {
      std::istringstream advancedStream(dummyString);
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
   * Returns the minimum energy to be generated.
   */
  double minE() const {return minE_;};
  /*!
   * Returns the maximum energy to be generated.
   */
  double maxE() const {return maxE_;};
  /*!
   * Returns the minimum angle to be generated.
   */
  double minA() const {return minA_;};
  /*!
   * Returns the maximum angle to be generated.
   */
  double maxA() const {return maxA_;};
  /*!
   * Returns the size energy step between generated points.
   */
  double eStep() const {return eStep_;};
  /*!
   * Returns the size angle step between generated points.
   */
  double aStep() const {return aStep_;};
  /*!
   * Return 0 if the segment is angle-integrated cross section, 1 for 
   * differential cross section, and 2 for phase shift.
   */
  int isDiff() const {return isDiff_;};
  /*!
   * Returns the spin value for the segment if the segment is to contain 
   * phase shift.
   */
  double phaseJ() const {return phaseJ_;};
  /*!
   * Returns the orbital angular momentum value for the segment
   * if the segment is to contain phase shift.
   */
  int phaseL() const {return phaseL_;};
  /*!
   * Returns the maximum polynomial order if segment is
   * angular distribution.
   */
  int maxAngDistOrder() const {return maxAngDistOrder_;};
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
  double eStep_;
  double aStep_;
  int isDiff_;
  double phaseJ_;
  int phaseL_;
  int maxAngDistOrder_;
  int isAdvanced_;
  int operationType_;
  std::string componentsList_;
};

#endif
