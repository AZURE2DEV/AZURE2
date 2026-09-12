#ifndef DATALINE_H
#define DATALINE_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

/// A class to read and store a line from a data file.

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
    // Optional numeric columns after the four required ones (for example the
    // per-point energy window of a beam-profile experimental effect).  The
    // rest of the line is consumed here, which for an ordinary four-column
    // file only swallows the newline the next read would have skipped anyway.
    if (stream.good()) {
      std::string rest;
      std::getline(stream, rest);
      std::istringstream rs(rest);
      double value;
      while (rs >> value) extra_.push_back(value);
    }
  };
  /*!
   * Number of optional extra columns read after the four required ones.
   */
  int numExtra() const { return static_cast<int>(extra_.size()); };
  /*!
   * Value of optional extra column \p i (0-based).
   */
  double extra(int i) const { return extra_[i]; };
  /*!
   * Returns the angle for the read in data point.
   */
  double angle() const { return angle_; };
  /*!
   * Returns the energy for the read in data point.
   */
  double energy() const { return energy_; };
  /*!
   * Returns the cross section for the read in data point.
   */
  double crossSection() const { return crossSection_; };
  /*!
   * Returns the cross section error for the read in data point.
   */
  double error() const { return error_; };

 private:
  double angle_;
  double energy_;
  double crossSection_;
  double error_;
  std::vector<double> extra_;
};

#endif
