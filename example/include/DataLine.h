#ifndef DATALINE_H
#define DATALINE_H

#include <fstream>

///A class to read and store a line from a data file.

/*!
 * The DataLine class reads and stores a formatted line from a data file. 
 */

class DataLine {
 public:
  /*!
   * Constructor fills the DataLine object from an input stream.
   */
  DataLine(std::ifstream &stream) {
    stream >> energy_ >> angle_ >> crossSection_ >> error_;
  };
  /*!
   * Returns the angle for the read in data point.
   */
  double angle() const {return angle_;};
  /*!
   * Returns the energy for the read in data point.
   */
  double energy() const {return energy_;};
  /*!
   * Returns the cross section for the read in data point.
   */
  double crossSection() const {return crossSection_;};
  /*!
   * Returns the cross section error for the read in data point.
   */
  double error() const {return error_;};
 private:
  double angle_;
  double energy_;
  double crossSection_;
  double error_;
};

#endif
