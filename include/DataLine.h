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
    // One physical line per DataLine, parsed from its own string stream.
    //
    // The old form read four numbers straight off the file stream inside a
    // "while (!in.eof())" loop.  A line that does not parse -- a '#' header,
    // a blank line, a stray word -- sets failbit without consuming anything,
    // eof is then never reached, and AZURE2 spins at 100% CPU forever with no
    // message.  Reading a whole line first means the stream always advances,
    // so the caller's loop always terminates, and the file's contents cannot
    // leave the stream in a failed state.
    //
    // Blank lines and lines whose first non-blank character is '#' are
    // comments and are skipped.  Anything else must carry at least the four
    // required numeric columns; if it does not, valid() is false and raw()
    // holds the text so the caller can say which line was wrong.
    //
    // A line ends at LF, CRLF or a bare CR.  Files written with classic-Mac
    // line endings (CR only) are otherwise read as a single line: its first
    // four numbers become one point and every later number is taken as an
    // optional extra column, so all but one point of the segment are lost
    // without a message.
    std::string raw;
    while (ReadLine(stream, raw)) {
      linesConsumed_++;
      std::string::size_type k = raw.find_first_not_of(" \t");
      if (k == std::string::npos || raw[k] == '#') continue;   // blank or comment
      std::istringstream ls(raw);
      if (ls >> energy_ >> angle_ >> crossSection_ >> error_) {
        // Optional numeric columns after the four required ones (for example
        // the per-point energy window of a beam-profile experimental effect).
        double value;
        while (ls >> value) extra_.push_back(value);
        valid_ = true;
      } else {
        raw_ = raw;          // report it; valid_ stays false
      }
      return;
    }
    atEnd_ = true;           // no data line left in the file
  };
  /*!
   * True when the constructor found no further line to read: the file is
   * exhausted.  The caller's read loop should stop here.
   */
  bool atEnd() const { return atEnd_; };
  /*!
   * Reads one line terminated by LF, CRLF or CR; false at end of stream with
   * nothing read.
   */
  static bool ReadLine(std::istream &in, std::string &out) {
    out.clear();
    char c;
    bool any = false;
    while (in.get(c)) {
      any = true;
      if (c == '\n') return true;
      if (c == '\r') {
        if (in.peek() == '\n') in.get(c);
        return true;
      }
      out.push_back(c);
    }
    return any;
  };
  /*!
   * True when a data line was read and its four required columns parsed.
   * False, with atEnd() also false, means a non-comment line that could not
   * be parsed; raw() then holds its text.
   */
  bool valid() const { return valid_; };
  /*!
   * The text of an unparseable line, for the caller's error message.
   */
  const std::string &raw() const { return raw_; };
  /*!
   * Physical lines consumed from the stream, comments and blanks included,
   * so the caller can keep a running line number for messages.
   */
  int linesConsumed() const { return linesConsumed_; };
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
  bool atEnd_ = false;
  bool valid_ = false;
  std::string raw_;
  int linesConsumed_ = 0;
};

#endif
