#ifndef EXTRAPLINE_H
#define EXTRAPLINE_H

#include <iostream>
#include <string>

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
    stream >> isUPOS_ >> secondaryDecayL_ >> Ic_ >> delta_;

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
   * Returns flag for if the segment is for an unobserved primary,
   * observed secondary.
   */
  int isUPOS() const {return isUPOS_;}; //1 if it is, 0 if not
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
   * Returns the multipolarity mixing ratio for unobserved particle
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
  double eStep_;
  double aStep_;
  int isDiff_;
  double phaseJ_;
  int phaseL_;
  int maxAngDistOrder_;
  int isUPOS_;
  int secondaryDecayL_;
  double Ic_;
  double delta_;
};

#endif
