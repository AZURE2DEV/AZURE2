#ifndef EXTRAPLINE_H
#define EXTRAPLINE_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "ComponentList.h"
#include "SegLine.h"  // ReadThicknessBlock

/// A class to read and store a line from the extrapolation input file.

/*!
 * The ExtrapLine class reads and stores a line from the extrapolation input file.
 */

class ExtrapLine {
 public:
  /*!
   * Constructor fill the ExtrapLine object from an input stream.
   */
  ExtrapLine(std::istream &stream) {
    stream >> isActive_ >> entranceKey_ >> exitKey_ >> minE_ >> maxE_ >> eStep_ >> minA_ >> maxA_ >> aStep_ >> isDiff_;
    if (isDiff_ == 2)
      stream >> phaseJ_ >> phaseL_;
    else if (isDiff_ == 3)
      stream >> maxAngDistOrder_;
    else if (isDiff_ == 10)
      stream >> spinObservable_;  // a spin observable (SpinObservable.h)

    // Initialize advanced segment parameters to defaults
    isAdvanced_ = 0;
    operationType_ = 0;
    componentsList_ = "";
    parseError_ = "";
    parseWarning_ = "";
    natural_ = "0";
    total_ = 0;
    transmission_ = 0.0;
    yield_ = 0.0;
    phaseS_ = -1.0;  // phase-shift channel spin: < 0 = not given ("spin <s>" keyword)
    // Initialize UPOS parameters to defaults
    isUPOS_ = 0;
    secondaryDecayL_ = 0;
    Ic_ = 0.0;
    delta_ = 0.0;

    // Try to read advanced segment data from the remaining line. The fields
    // above are the whole line for a plain segment, so the stream is already at
    // eof here; getline would then set failbit and the caller reads that as a
    // malformed line, silently rejecting a perfectly valid segment. Only read
    // when something is actually left.
    std::string dummyString;
    if (!stream.eof()) getline(stream, dummyString);

    if (!dummyString.empty()) {
      std::istringstream advancedStream(dummyString);
      int tempIsAdvanced;
      if (advancedStream >> tempIsAdvanced && tempIsAdvanced == 1) {
        isAdvanced_ = 1;
        if (advancedStream >> operationType_) {
          const std::string componentError = ReadComponentList(advancedStream, componentsList_, parseWarning_);
          if (!componentError.empty()) {
            parseError_ = componentError;
            stream.setstate(std::ios::failbit);
          }
        }
      }
      // Optional keyword blocks at the end of the line (SegLine.h): a plain
      // line has no advanced block, so the integer read above failed on the
      // keyword itself and left the stream positioned on it.
      advancedStream.clear();
      std::string keyword;
      while (advancedStream >> keyword) {
        if (keyword == "sample" || keyword == "natural") {
          if (!(advancedStream >> natural_)) {
            natural_ = "0";
            break;
          }
        } else if (keyword == "spin") {
          if (!(advancedStream >> phaseS_) || phaseS_ < 0.0) {
            stream.setstate(std::ios::failbit);
            break;
          }
        } else if (keyword == "total") {
          if (!(advancedStream >> total_) || (total_ != 0 && total_ != 1)) {
            stream.setstate(std::ios::failbit);
            break;
          }
        } else if (keyword == "transmission" || keyword == "yield") {
          // The fit form "n vary sigma" is read and its vary / sigma ignored: a test segment has no parameter and
          // uses n (EData::MakePoints notes it).
          double &n = (keyword == "transmission") ? transmission_ : yield_;
          bool fit = false;
          int vary = 0;
          double sigma = 0.0;
          parseError_ = ReadThicknessBlock(advancedStream, keyword, n, fit, vary, sigma);
          if (!parseError_.empty()) {
            stream.setstate(std::ios::failbit);
            break;
          }
          if (fit) thicknessFitIgnored_ = true;
        } else if (keyword == "covariance") {
          hasCovariance_ = true;  // a data-segment block only: EData::MakePoints refuses the line
          break;
        } else if (keyword == "sqrtshift") {
          // data-segment blocks: a test segment has no shift, norm or shift function.  They used to stop the
          // keyword loop here, losing a "sample" or "total" written after them (audit 2026-10-25, item 29);
          // now they are skipped with their values and the loop goes on.
          std::string a, b, c;
          if (!(advancedStream >> a >> b >> c)) break;
        } else if (keyword == "shiftgroup" || keyword == "shiftfngroup" || keyword == "normgroup" || keyword == "absolute" ||
                   keyword == "thicknessgroup") {
          std::string a;
          if (!(advancedStream >> a)) break;
          if (keyword == "thicknessgroup") thicknessFitIgnored_ = true;
        } else if (keyword == "shiftfn") {
          std::string expression;
          int numParams = 0;
          if (!(advancedStream >> expression)) break;
          if (expression.size() > 1 && expression[0] == '"' && expression[expression.size() - 1] != '"') {
            std::string rest;
            if (std::getline(advancedStream, rest, '"')) expression += rest;
          }
          if (!(advancedStream >> numParams)) break;
          bool okParams = true;
          for (int j = 0; j < 3 * numParams && okParams; j++) {
            std::string v;
            if (!(advancedStream >> v)) okParams = false;
          }
          if (!okParams) break;
        } else {
          break;
        }
      }
    }
  };
  /*!
   * Returns non-zero if the line is to be included in the calculation.
   */
  int isActive() const { return isActive_; };
  /*!
   * Returns true if the line carries a "covariance" block, which only data segments take.
   */
  bool hasCovariance() const { return hasCovariance_; };
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
   * Returns the minimum energy to be generated.
   */
  double minE() const { return minE_; };
  /*!
   * Returns the maximum energy to be generated.
   */
  double maxE() const { return maxE_; };
  /*!
   * Returns the minimum angle to be generated.
   */
  double minA() const { return minA_; };
  /*!
   * Returns the maximum angle to be generated.
   */
  double maxA() const { return maxA_; };
  /*!
   * Returns the size energy step between generated points.
   */
  double eStep() const { return eStep_; };
  /*!
   * Returns the size angle step between generated points.
   */
  double aStep() const { return aStep_; };
  /*!
   * Return 0 if the segment is angle-integrated cross section, 1 for
   * differential cross section, and 2 for phase shift.
   */
  int isDiff() const { return isDiff_; };
  /*!
   * Returns the spin value for the segment if the segment is to contain
   * phase shift.
   */
  double phaseJ() const { return phaseJ_; };
  /*!
   * Returns the orbital angular momentum value for the segment
   * if the segment is to contain phase shift.
   */
  int phaseL() const { return phaseL_; };
  //! The observable token of an isDiff 10 segment (empty otherwise).
  std::string spinObservable() const { return spinObservable_; };
  /*!
   * The channel spin of a phase-shift test segment ("spin <s>" keyword), or -1 when not given.
   */
  double phaseS() const { return phaseS_; };
  /*!
   * Returns the maximum polynomial order if segment is
   * angular distribution.
   */
  int maxAngDistOrder() const { return maxAngDistOrder_; };
  /*!
   * Returns non-zero if this is an advanced segment (sum/ratio).
   */
  int isAdvanced() const { return isAdvanced_; };
  /*!
   * Returns the operation type (0 for sum, 1 for ratio).
   */
  int operationType() const { return operationType_; };
  /*!
   * Returns the sample token ("sample <t>", alias "natural <t>"): "0" (none), "1" (the default sample) or a sample name.
   */
  const std::string &natural() const { return natural_; };
  /// 1 for "total 1" (SegLine::total).
  int total() const { return total_; };
  /// "transmission <n>" / "yield <n>" (SegLine.h): 0 when absent.
  double transmission() const { return transmission_; };
  double yieldThickness() const { return yield_; };
  /// True when the line carries the fit form "n vary sigma" or a thicknessgroup block: a test segment ignores both
  /// and uses its own n.
  bool thicknessFitIgnored() const { return thicknessFitIgnored_; };
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
   * Returns the spin of the final state for unobserved primary, observed secondary.
   */
  double Ic() const { return Ic_; };
  /*!
   * Returns the multipole mixing ratio for unobserved primary, observed secondary.
   */
  double delta() const { return delta_; };

  /// What could not be read when the line failed (empty when it did not, or for an older failure without a message).
  std::string parseError() const { return parseError_; };
  /// A note about the line that did not fail it (printed by EData with the segment number).
  std::string parseWarning() const { return parseWarning_; };

 private:
  std::string parseError_;
  std::string parseWarning_;
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
  std::string spinObservable_;
  double phaseS_;
  int maxAngDistOrder_;
  int isAdvanced_;
  int operationType_;
  std::string natural_;
  int total_;
  double transmission_;
  double yield_;
  bool thicknessFitIgnored_ = false;
  bool hasCovariance_ = false;
  std::string componentsList_;
  int isUPOS_;
  int secondaryDecayL_;
  double Ic_;
  double delta_;
};

#endif
