#ifndef AZURE2_CONVERT_TARGETINTLINE_H
#define AZURE2_CONVERT_TARGETINTLINE_H

#include <cctype>
#include <cstdlib>
#include <istream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace azure2 {
namespace convert {

/*!
 * One line of the old <targetInt> section as it is written, before any meaning is given to it: the grammar of
 * TargetEffect's constructor (src/TargetEffect.cpp), moved here unchanged so that the engine (until the format-3
 * reader replaces the section, step 3b) and the converter to format 3 read the line with one parser.
 *
 *   active "segments" points isConvolution sigma isTargetIntegration density "stopping" nParams p_1 .. p_n
 *     isQ nQ Q_1 .. Q_n  isConvCoefficients "expression" nConv c_1 .. c_n
 *     [isStraggling [coefficient [widthMultiplier [pointsPerWidth]]]]  ["lo-hi,...|auto[:tol[:margin]]" [blend [autoTol]]]
 *     { beamprofile n {xi omega alpha weight} tpc truncation db | sigmalab | screening Ue vary sigma
 *     | densityfit vary sigma[%] [max] | densitygroup label | doppler T [vary sigma [Tmax]]
 *     | udr "file" flightPath burst channel centred }
 *
 * Fields keep the values and spellings as read; whether they make sense (a densityfit on a thick target, a doppler
 * fit form) is the reader's to decide.  The stream's state after the call is what the engine checks (failbit: a
 * malformed line), with the eof-only failure of an exactly consumed line cleared as before.
 */
struct TargetIntLine {
  struct ProfileComponent {
    double xi = 0, omega = 0, alpha = 0, weight = 0;
  };

  /// false when the line ends inside the first nine fields (the engine then sets nothing else)
  bool complete = false;
  int isActive = 0;
  std::string segmentList;  ///< as read, quotes included
  int numIntegrationPoints = 0;
  int isConvolution = 0;
  double sigma = 0.0;
  int isTargetIntegration = 0;
  double density = 0.0;
  std::string stoppingPowerEq;  ///< as read, quotes included
  int numParameters = 0;
  std::vector<double> parameters;
  int isQCoefficients = 0;
  int numQCoefficients = 0;
  std::vector<double> qCoefficients;
  int isConvCoefficients = 0;
  int numConvCoefficients = 0;
  std::vector<double> convCoefficients;
  std::string convolutionEq;  ///< as read, quotes included

  // the optional straggling and adaptive-grid fields
  bool hasStragglingFlag = false;
  int isStraggling = 0;
  bool hasStragglingCoefficient = false;
  double stragglingCoefficient = 0.0;
  bool hasResonanceWidthMultiplier = false;
  double resonanceWidthMultiplier = 0.0;
  bool hasPointsPerWidth = false;
  double pointsPerWidth = 0.0;

  // the optional energy windows
  bool hasRangesToken = false;
  std::string rangesToken;  ///< as read, quotes removed
  std::vector<std::pair<double, double>> ranges;
  bool autoWindows = false;
  bool hasAutoWindowTolerance = false;
  double autoWindowTolerance = 0.0;
  bool hasAutoWindowMargin = false;
  double autoWindowMargin = 0.0;
  double transitionWidth = 0.0;  ///< 0 when absent
  double autoTolerance = 0.0;    ///< 0 when absent

  // the keyword blocks
  std::vector<ProfileComponent> beamProfile;
  double beamTpcSigma = 0.0;
  double beamTruncation = 0.0;
  int beamDetailedBalance = 0;
  bool isBeamProfile = false;
  bool sigmaLab = false;
  bool hasScreening = false;
  bool screeningError = false;
  double screeningNominal = 0.0;
  int screeningVary = 0;
  double screeningSigma = 0.0;
  bool densityFitRead = false;
  int densityVary = 0;
  std::string densitySigmaToken;
  double densityMaxRead = 0.0;
  bool densityFitError = false;
  std::string densityFitErrorText;
  std::string densityGroup;
  bool isDoppler = false;
  bool dopplerError = false;
  std::string dopplerErrorText;
  double dopplerTemperature = 0.0;
  bool dopplerFitRead = false;
  std::string dopplerVaryToken, dopplerSigmaToken, dopplerMaxToken;
  bool isUdr = false;
  std::string udrFile;
  double udrFlightPath = 0.0;
  double udrBurstFwhm = 0.0;
  double udrChannelWidth = 0.0;
  int udrCentred = 1;
};

/// Strips every '"' of a token (the engine's reading of the quoted list and expressions).
inline std::string TargetIntUnquoted(std::string text) {
  std::size_t found = 0;
  while ((found = text.find('\"')) != std::string::npos) text.erase(found, 1);
  return text;
}

/// The segment keys of the line's list ("3,6,7-14", quotes removed): 3 6 7 8 ... 14, in the order written.
inline std::vector<int> TargetIntSegmentKeys(const std::string &list) {
  std::vector<int> keys;
  std::size_t i = 0;
  int lastSegNum = 0;
  bool inclusive = false;
  while (i < list.length()) {
    if (list[i] >= '0' && list[i] <= '9') {
      std::string tempString;
      while (list[i] != ',' && list[i] != '-' && i < list.length()) {
        tempString += list[i];
        i++;
      }
      std::istringstream stm;
      stm.str(tempString);
      int tempSegNum;
      stm >> tempSegNum;
      if (inclusive == true)
        for (int j = lastSegNum + 1; j <= tempSegNum; j++) keys.push_back(j);
      else
        keys.push_back(tempSegNum);
      lastSegNum = tempSegNum;
    }
    if (list[i] == '-')
      inclusive = true;
    else
      inclusive = false;
    i++;
  }
  return keys;
}

/// Reads one <targetInt> line from `stream` into `out` (see TargetIntLine).
inline void ReadTargetIntLine(std::istream &stream, TargetIntLine &out) {
  stream >> out.isActive >> out.segmentList >> out.numIntegrationPoints >> out.isConvolution >> out.sigma >>
      out.isTargetIntegration >> out.density >> out.stoppingPowerEq >> out.numParameters;
  if (stream.eof()) return;
  out.complete = true;
  for (int i = 0; i < out.numParameters; i++) {
    double tempParameter = 0.0;
    stream >> tempParameter;
    out.parameters.push_back(tempParameter);
  }

  stream >> out.isQCoefficients >> out.numQCoefficients;
  for (int i = 0; i < out.numQCoefficients && stream.good(); i++) {
    double tempQCoefficient = 0.0;
    stream >> tempQCoefficient;
    out.qCoefficients.push_back(tempQCoefficient);
  }

  stream >> out.isConvCoefficients >> out.convolutionEq >> out.numConvCoefficients;
  for (int i = 0; i < out.numConvCoefficients && stream.good(); i++) {
    double tempConvCoefficient = 0.0;
    stream >> tempConvCoefficient;
    out.convCoefficients.push_back(tempConvCoefficient);
  }

  // Straggling flag and coefficient, then the adaptive-grid fields, each optional (older lines stop earlier).
  // Skip whitespace and check for a digit (not a '<', which would be </targetInt>).
  if (stream.good()) stream >> std::ws;
  if (stream.good() && stream.peek() != '<' && std::isdigit(stream.peek())) {
    stream >> out.isStraggling;
    out.hasStragglingFlag = true;
    if (stream.good()) stream >> std::ws;
    if (stream.good() && stream.peek() != '<' && (std::isdigit(stream.peek()) || stream.peek() == '.' || stream.peek() == '-')) {
      stream >> out.stragglingCoefficient;
      out.hasStragglingCoefficient = true;
      if (stream.good()) stream >> std::ws;
      if (stream.good() && stream.peek() != '<' && (std::isdigit(stream.peek()) || stream.peek() == '.' || stream.peek() == '-')) {
        double rwm;
        stream >> rwm;
        if (stream) {
          out.resonanceWidthMultiplier = rwm;
          out.hasResonanceWidthMultiplier = true;
          if (stream.good()) stream >> std::ws;
          if (stream.good() && stream.peek() != '<' && (std::isdigit(stream.peek()) || stream.peek() == '.' || stream.peek() == '-')) {
            double ppw;
            stream >> ppw;
            if (stream) {
              out.pointsPerWidth = ppw;
              out.hasPointsPerWidth = true;
            }
          }
        }
      }
    }
  }

  // Optional restriction of the effect to lab-energy windows: a quoted token "lo1-hi1,lo2-hi2" (empty "" means
  // unrestricted), then optionally the blend width at each edge (MeV) and the automatic-application tolerance.
  // Every probe is guarded on good(): probing past the end of an exactly consumed line is not an error.
  if (stream.good()) {
    stream >> std::ws;
    if (stream.good() && stream.peek() == '"') {
      std::string rangesToken;
      stream >> rangesToken;
      out.hasRangesToken = true;
      rangesToken = TargetIntUnquoted(rangesToken);
      out.rangesToken = rangesToken;
      std::istringstream rs(rangesToken);
      std::string item;
      while (std::getline(rs, item, ',')) {
        if (item.compare(0, 4, "auto") == 0) {
          // "auto", "auto:<tol>", "auto:<tol>:<margin>": windows from the level scheme
          out.autoWindows = true;
          std::istringstream as(item.substr(4));
          char colon = 0;
          double v = 0.0;
          if (as >> colon && colon == ':' && as >> v && v > 0.) {
            out.autoWindowTolerance = v;
            out.hasAutoWindowTolerance = true;
          }
          colon = 0;
          if (as >> colon && colon == ':' && as >> v && v > 0.) {
            out.autoWindowMargin = v;
            out.hasAutoWindowMargin = true;
          }
          continue;
        }
        std::size_t dash = item.find('-', 1);  // energies are positive; skip a leading sign
        if (dash == std::string::npos) continue;
        double lo = std::atof(item.substr(0, dash).c_str());
        double hi = std::atof(item.substr(dash + 1).c_str());
        if (hi > lo) out.ranges.push_back(std::pair<double, double>(lo, hi));
      }
      if (stream.good()) {
        stream >> std::ws;
        if (stream.good() && stream.peek() != '<' && (std::isdigit(stream.peek()) || stream.peek() == '.')) {
          double tw;
          stream >> tw;
          if (!stream.fail() && tw > 0.) out.transitionWidth = tw;
          if (stream.good()) {
            stream >> std::ws;
            if (stream.good() && stream.peek() != '<' && (std::isdigit(stream.peek()) || stream.peek() == '.')) {
              double tol;
              stream >> tol;
              if (!stream.fail() && tol > 0.) out.autoTolerance = tol;
            }
          }
        }
      }
    }
  }

  // Optional keyword blocks, introduced by a word so that older readers (which only ever probe for digits or a
  // quote) stop in front of them.
  while (stream.good()) {
    stream >> std::ws;
    if (!(stream.good() && std::isalpha(stream.peek()))) break;
    std::string keyword;
    stream >> keyword;
    if (keyword == "beamprofile") {
      int numComponents = 0;
      stream >> numComponents;
      for (int i = 0; i < numComponents && stream.good(); i++) {
        TargetIntLine::ProfileComponent c;
        stream >> c.xi >> c.omega >> c.alpha >> c.weight;
        out.beamProfile.push_back(c);
      }
      int dbFlag = 0;
      stream >> out.beamTpcSigma >> out.beamTruncation >> dbFlag;
      out.beamDetailedBalance = dbFlag;
      out.isBeamProfile = !out.beamProfile.empty() && !stream.fail();
      if (out.beamProfile.empty()) stream.setstate(std::ios_base::failbit);
    } else if (keyword == "sigmalab") {
      out.sigmaLab = true;
    } else if (keyword == "screening") {
      int vary = 0;
      stream >> out.screeningNominal >> vary >> out.screeningSigma;
      if (stream.fail() || out.screeningSigma < 0.0) {
        // Flagged as well as failed: at the end of the line the fail bit is cleared below.
        out.screeningError = true;
        stream.setstate(std::ios_base::failbit);
      } else {
        out.hasScreening = true;
        out.screeningVary = vary;
      }
    } else if (keyword == "densityfit") {
      int vary = 0;
      stream >> vary >> out.densitySigmaToken;
      if (stream.fail() || out.densitySigmaToken.empty()) {
        out.densityFitError = true;
        out.densityFitErrorText =
            "a densityfit block is 'densityfit vary sigma' (sigma in atoms/cm^2 or a percentage, \"2.3%\")";
        stream.setstate(std::ios_base::failbit);
      } else {
        out.densityFitRead = true;
        out.densityVary = vary;
        // optional: the largest density the grids must span (atoms/cm^2)
        if (stream.good()) stream >> std::ws;
        if (stream.good() && (std::isdigit(stream.peek()) || stream.peek() == '.')) {
          stream >> out.densityMaxRead;
          if (stream.fail()) {
            out.densityFitError = true;
            out.densityFitErrorText = "the optional densityfit maximum is not a number";
          }
        }
      }
    } else if (keyword == "densitygroup") {
      stream >> out.densityGroup;
      if (stream.fail() || out.densityGroup.empty()) {
        out.densityFitError = true;
        out.densityFitErrorText = "a densitygroup block is 'densitygroup <label>'";
        stream.setstate(std::ios_base::failbit);
      }
    } else if (keyword == "doppler") {
      stream >> out.dopplerTemperature;
      if (stream.fail() || !(out.dopplerTemperature > 0.0)) {
        out.dopplerError = true;  // flagged as well as failed: at the end of the line the fail bit is cleared below
        out.dopplerErrorText = "a doppler block is 'doppler T' or 'doppler T vary sigma [Tmax]' with the effective"
                               " temperature T > 0 in kelvin";
        stream.setstate(std::ios_base::failbit);
      } else {
        out.isDoppler = true;
        // optional fit form: vary sigma [Tmax], numbers (a following keyword starts with a letter)
        auto numberToken = [&](std::string &tok) {
          if (stream.good()) stream >> std::ws;
          if (!stream.good()) return false;
          const int c = stream.peek();
          if (!(std::isdigit(c) || c == '.' || c == '-' || c == '+')) return false;
          stream >> tok;
          return !stream.fail();
        };
        if (numberToken(out.dopplerVaryToken)) {
          out.dopplerFitRead = true;
          if (!numberToken(out.dopplerSigmaToken)) {
            out.dopplerError = true;
            out.dopplerErrorText = "a doppler block is 'doppler T' or 'doppler T vary sigma [Tmax]': the fit form needs"
                                   " both vary (0 or 1) and sigma (the prior width in K, 0 = none)";
            stream.setstate(std::ios_base::failbit);
          } else {
            numberToken(out.dopplerMaxToken);
          }
        }
      }
    } else if (keyword == "udr") {
      stream >> std::ws;
      std::string file;
      if (stream.peek() == '"') {
        stream.get();
        std::getline(stream, file, '"');
      } else {
        stream >> file;
      }
      int centred = 1;
      stream >> out.udrFlightPath >> out.udrBurstFwhm >> out.udrChannelWidth >> centred;
      out.udrFile = file;
      out.udrCentred = centred;
      out.isUdr = !out.udrFile.empty() && out.udrFlightPath > 0.0 && !stream.fail();
      if (!out.isUdr) stream.setstate(std::ios_base::failbit);
    } else {
      stream.setstate(std::ios_base::failbit);
    }
  }
  // A line consumed exactly to its end during the optional probes is not a parse error; genuinely malformed fields
  // fail without reaching eof.
  if (stream.fail() && stream.eof()) stream.clear(std::ios_base::eofbit);
}

}  // namespace convert
}  // namespace azure2

#endif
