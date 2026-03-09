#ifndef SEGLINE_H
#define SEGLINE_H

#include <iostream>
#include <string>

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
    stream >> dataNorm_ >> varyNorm_ >> dataNormError_ >> isUPOS_ >> secondaryDecayL_ >> Ic_>>delta_;
//    std::cout<<isUPOS_<<","<<secondaryDecayL_<<","<<Ic_<<","<<delta_<<std::endl;
    std::string dummyString;
    getline(stream,dummyString);
    int p2 = dummyString.find_last_not_of(" \n\t\r");
    if (p2 != std::string::npos) {  
      int p1 = dummyString.find_first_not_of(" \n\t\r");
      if (p1 == std::string::npos) p1 = 0;
      dataFile_=dummyString.substr(p1,(p2-p1)+1);
    } else dataFile_=std::string();  
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
   * Return 0 if the segment is angle-integrated cross section, 1 for lab
   * differential cross section, 2 for phase shift, 3 for angle integrated
   * total capture, and 4 for C.M. differential cross section for data.
   * for extrapolation segments return 0 for angle-integrated cross section,
   * 1 for lab differential cross section, 2 for phase shifts, 3 for 
   * angular distrabution coefficients, 4 for angle-integrated total
   * capture, and 5 for C.M. differential cross sections.
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
   * Returns flag for if this segment is an unobserved primary transition (1 = is, 0 = not).
   */
  int isUPOS() const {return isUPOS_;}; //isUPOS_
  /*!
   * Returns the angular momentum of the decay for unobserved particle
   * observed secondary decay
   */
  int secondaryDecayL() const {return secondaryDecayL_;}; //secondaryDecayL_
  /*!
   * Returns the spin of the final state of the decay for unobserved particle
   * observed secondary decay
   */
  double Ic() const {return Ic_;}; //Ic_
  /*!
   * Returns the multipole mixing ratio of the decay for unobserved particle
   * observed secondary decay
   */
  double delta() const {return delta_;}; //Ic_

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
  int isUPOS_;
  int secondaryDecayL_;
  double Ic_;
  double delta_;
};

#endif
