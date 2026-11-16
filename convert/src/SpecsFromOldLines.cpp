// TEMPORARY (step 3d deletes it): the old sections as set-up specs (SpecsFromOldLines.h).  The loops are the old
// engine's reading loops (WorkingCopy::Fill, WorkingCopy::MakePoints, WorkingCopy::ReadTargetEffectsFile, CNuc::ReadAbundances,
// WorkingCopy::ReadPriorCentres) with the set-up taken out: what they printed while reading a line is printed here, in the
// same words, and everything else is left to the engine, which sets up the specs.
#include "azure2/convert/SpecsFromOldLines.h"

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <vector>

#include "azure2/convert/AbundanceLine.h"
#include "azure2/convert/ExtrapLine.h"
#include "azure2/convert/SegLine.h"
#include "azure2/convert/TargetIntLine.h"

namespace azure2 {
namespace convert {

using experiment::ComponentSpec;
using experiment::DatasetSpec;
using experiment::ExtrapolationSpec;
using experiment::Observable;
using experiment::ParameterNames;
using experiment::ReactionSpec;
using experiment::SpecList;
using experiment::TargetEffectSpec;

namespace {

bool Blank(const std::string &line) {
  for (char c : line)
    if (c != ' ' && c != '\t') return false;
  return true;
}

// isDiff of a <segmentsData> line without the THM offset of 10.  A code the old engine gave no meaning computed the
// angle-integrated cross section (no flag of Segment set): it is kept so.
Observable DataObservable(int code) {
  switch (code) {
    case 1: return Observable::DifferentialLab;
    case 2: return Observable::PhaseShift;
    case 3: return Observable::TotalCapture;
    case 4: return Observable::DifferentialCM;
    case 5: return Observable::CrossSectionE1;
    case 6: return Observable::CrossSectionE2;
    case 7: return Observable::AnalyzingPower;
    case 9: return Observable::AngularCoefficients;
    default: return Observable::CrossSection;
  }
}

// isDiff of a <segmentsTest> line without the THM offset (its codes differ from the data lines' from 3 on).
Observable ExtrapolationObservable(int code) {
  switch (code) {
    case 1: return Observable::DifferentialLab;
    case 2: return Observable::PhaseShift;
    case 3: return Observable::AngularDistribution;
    case 4: return Observable::TotalCapture;
    case 5: return Observable::DifferentialCM;
    case 7: return Observable::AnalyzingPower;
    default: return Observable::CrossSection;
  }
}

// The fields a data line and a test line share.  Returns "" or why the line cannot be set up.
template <class Line>
std::string ReadReaction(const Line &line, Observable observable, bool thm, ReactionSpec &r) {
  r.entranceKey = line.entranceKey();
  r.exitKey = line.exitKey();
  r.observable = observable;
  r.thm = thm;
  r.lineWindow_keV = line.lineWindow();
  if (observable == Observable::PhaseShift) {
    r.phaseJ = line.phaseJ();
    r.phaseL = line.phaseL();
    r.phaseS = line.phaseS();
  }
  if (line.isAdvanced() == 1) {
    if (line.operationType() != 0 && line.operationType() != 1)
      return "its composite block has operation type " + std::to_string(line.operationType()) +
             ", which is neither 0 (sum) nor 1 (ratio)";
    r.combination = line.operationType() == 1 ? ReactionSpec::Combination::Ratio : ReactionSpec::Combination::Sum;
  }
  for (const ComponentEntry &e : line.components()) {
    ComponentSpec c;
    c.entranceKey = e.entrance;
    c.exitKey = e.exit;
    c.hasAngle = e.hasAngle;
    c.angle = e.hasAngle ? e.angle : 0.0;
    c.scaling = e.hasScaling ? e.scaling : 1.0;
    r.components.push_back(c);
  }
  r.sample = line.natural();
  r.total = line.total() == 1;
  r.transmission = line.transmission();
  r.yieldThickness = line.yieldThickness();
  return std::string();
}

// The old names of a segment's fit parameters: segment_<key>_norm, ... (the names param files, <parameterSettings>
// rows and prior_centre rows of an old project use).
ParameterNames OldParameterNames(int key, int shiftFunctionTerms) {
  const std::string stem = "segment_" + std::to_string(key) + "_";
  ParameterNames n;
  n.norm = stem + "norm";
  n.energyShift = stem + "energy_shift";
  n.energyShiftSqrt = stem + "energy_shift_sqrt";
  for (int j = 0; j < shiftFunctionTerms; j++) n.shiftFunction.push_back(stem + "energy_shift_fn" + std::to_string(j));
  n.thickness = stem + "thickness";
  return n;
}

DatasetSpec DatasetOf(const SegLine &s, int key, std::string &why) {
  DatasetSpec d;
  d.key = key;
  d.active = s.isActive() == 1;
  const bool thm = s.isDiff() >= 10;
  why = ReadReaction(s, DataObservable(thm ? s.isDiff() - 10 : s.isDiff()), thm, d.reaction);
  d.minE = s.minE();
  d.maxE = s.maxE();
  d.minA = s.minA();
  d.maxA = s.maxA();
  d.dataFile = s.dataFile();
  d.absolute = s.absolute() == 1;
  if (d.absolute) d.reaction.legendre = azure2::experiment::LegendreExpansion::Absolute;  // "absolute 1": A_k (b/sr)
  d.norm = s.dataNorm();
  d.normError = s.dataNormError();
  d.varyNorm = s.varyNorm() == 1;
  d.normGroup = s.normGroup();
  d.energyShift = s.energyShift();
  d.energyShiftError = s.energyShiftError();
  d.varyEnergyShift = s.varyEnergyShift() == 1;
  d.hasSqrtShift = s.hasSqrtShift();
  d.energyShiftSqrt = s.energyShiftSqrt();
  d.energyShiftSqrtError = s.energyShiftSqrtError();
  d.varyEnergyShiftSqrt = s.varyEnergyShiftSqrt() == 1;
  d.shiftFunction = s.shiftFn();
  d.shiftFunctionValues = s.shiftFnValues();
  d.names = OldParameterNames(key, (int)s.shiftFnValues().size());
  d.shiftFunctionErrors = s.shiftFnErrors();
  for (int v : s.shiftFnVary()) d.shiftFunctionVary.push_back(v == 1);
  d.shiftGroup = s.shiftGroup();
  d.shiftFunctionGroup = s.shiftFnGroup();
  d.covarianceFile = s.covarianceFile();
  d.thicknessFit = s.thicknessFit();
  d.thicknessVary = s.thicknessVary() == 1;
  d.thicknessSigma = s.thicknessSigma();
  d.thicknessGroup = s.thicknessGroup();
  d.upos = s.isUPOS() == 1;
  d.secondaryDecayL = s.secondaryDecayL();
  d.secondaryFinalSpin = s.Ic();
  d.secondaryMixingRatio = s.delta();
  return d;
}

ExtrapolationSpec ExtrapolationOf(const ExtrapLine &s, int key, std::string &why) {
  ExtrapolationSpec x;
  x.key = key;
  x.names = OldParameterNames(key, 0);
  x.active = s.isActive() == 1;
  const bool thm = s.isDiff() >= 10;
  why = ReadReaction(s, ExtrapolationObservable(thm ? s.isDiff() - 10 : s.isDiff()), thm, x.reaction);
  x.minE = s.minE();
  x.maxE = s.maxE();
  x.minA = s.minA();
  x.maxA = s.maxA();
  // The old loops: energy += eStep from Emin while energy <= Emax, inside angle += aStep from Amin while angle <= Amax;
  // a step of 0 gives one value.
  for (double e = s.minE(); e <= s.maxE(); e += s.eStep()) {
    x.energies.push_back(e);
    if (s.eStep() == 0.0) break;
  }
  for (double a = s.minA(); a <= s.maxA(); a += s.aStep()) {
    x.angles.push_back(a);
    if (s.aStep() == 0.0) break;
  }
  if (x.reaction.observable == Observable::AngularDistribution) x.maxAngDistOrder = s.maxAngDistOrder();
  return x;
}

}  // namespace

SpecList<DatasetSpec> DatasetSpecsFromOldLines(std::istream &in, bool relativistic) {
  SpecList<DatasetSpec> out;
  std::string line = "";
  getline(in, line);  // <segmentsData>
  line = "";
  int numTotalSegments = 0;
  while (!in.eof() && line != "</segmentsData>") {
    getline(in, line);
    if (Blank(line)) continue;
    if (in.eof() || line == "</segmentsData>") continue;
    std::istringstream stm(line);
    SegLine segment(stm);
    std::ostringstream message;
    auto stop = [&]() {
      out.failed = true;
      out.error = message.str();
      return out;
    };
    if (stm.rdstate() & (std::stringstream::failbit | std::stringstream::badbit)) {
      if (!segment.parseError().empty())
        message << "ERROR: Data segment #" << numTotalSegments + 1 << ": " << segment.parseError() << ".";
      return stop();
    }
    numTotalSegments++;
    std::string warning;
    if (!segment.parseWarning().empty())
      warning = "Data segment #" + std::to_string(numTotalSegments) + ": " + segment.parseWarning() + ".";
    if (!warning.empty()) message << "WARNING: " << warning << std::endl;
    if (segment.isDiff() == 19) {
      // isDiff 9 + the THM offset: the angular-distribution coefficients are not a THM observable.
      message << "ERROR: Data segment #" << numTotalSegments
              << " uses isDiff 19 (THM angular-distribution coefficients), which is not supported.";
      return stop();
    }
    if (segment.isDiff() >= 10 && relativistic) {
      // The THM machinery (vertex, line shape, distortion, the HOES kinematics) is non-relativistic.
      message << "ERROR: Data segment #" << numTotalSegments
              << " is a THM segment; THM segments are not supported with relativistic=1.";
      return stop();
    }
    if (segment.isDiff() >= 10 && !segment.covarianceFile().empty()) {
      // A THM segment's chi^2 is its experiment's (profiled norm, ThmNorm.cpp), not a sum a covariance could amend.
      message << "ERROR: Data segment #" << numTotalSegments
              << " is a THM segment with a covariance block; a data covariance is not supported on"
              << " THM segments.";
      return stop();
    }
    if (segment.isDiff() == 8 || segment.isDiff() == 18) {
      // Polarization x Cross Section was removed: the outgoing polarization of A(a,b)B is the analyzing power of
      // the inverse reaction B(b,a)A.
      message << "ERROR: Data segment #" << numTotalSegments << " uses isDiff " << segment.isDiff() << " (Polarization x"
              << " Cross Section), which is no longer supported." << std::endl
              << "       Divide the data by the differential cross section and enter the"
              << " polarization as an Analyzing Power (isDiff 7) segment on the inverse"
              << " channel (entrance = the measured exit pair, exit = the measured entrance pair).";
      return stop();
    }
    std::string why;
    SpecList<DatasetSpec>::Item item;
    item.spec = DatasetOf(segment, numTotalSegments, why);
    if (!why.empty()) {
      message << "ERROR: Data segment #" << numTotalSegments << ": " << why << ".";
      return stop();
    }
    item.warning = warning;
    out.items.push_back(item);
  }
  if (line != "</segmentsData>") out.failed = true;
  return out;
}

SpecList<ExtrapolationSpec> ExtrapolationSpecsFromOldLines(std::istream &in, bool relativistic) {
  SpecList<ExtrapolationSpec> out;
  std::string line = "";
  getline(in, line);  // <segmentsTest>
  line = "";
  int numTotalSegments = 0;
  while (!in.eof() && line != "</segmentsTest>") {
    getline(in, line);
    if (Blank(line)) continue;
    if (in.eof() || line == "</segmentsTest>") continue;
    std::istringstream stm(line);
    ExtrapLine segment(stm);
    std::ostringstream message;
    auto stop = [&]() {
      out.failed = true;
      out.error = message.str();
      return out;
    };
    if (stm.rdstate() & (std::stringstream::failbit | std::stringstream::badbit)) {
      if (!segment.parseError().empty())
        message << "ERROR: Test segment #" << numTotalSegments + 1 << ": " << segment.parseError() << ".";
      return stop();
    }
    numTotalSegments++;
    std::string warning;
    if (!segment.parseWarning().empty())
      warning = "Test segment #" + std::to_string(numTotalSegments) + ": " + segment.parseWarning() + ".";
    if (!warning.empty()) message << "WARNING: " << warning << std::endl;
    if (segment.hasCovariance()) {
      message << "ERROR: Test segment #" << numTotalSegments
              << " carries a covariance block; a data covariance belongs on data segments only.";
      return stop();
    }
    if (segment.isDiff() >= 10 && relativistic) {
      message << "ERROR: Test segment #" << numTotalSegments
              << " is a THM segment; THM segments are not supported with relativistic=1.";
      return stop();
    }
    if (segment.isDiff() == 9 || segment.isDiff() == 19) {
      message << "ERROR: Test segment #" << numTotalSegments << " uses isDiff " << segment.isDiff()
              << " (angular coefficients), a data-segment type; a test segment"
              << " computes the coefficients with isDiff 3 and a maximum order.";
      return stop();
    }
    if (segment.isDiff() == 8 || segment.isDiff() == 18) {
      message << "ERROR: Test segment #" << numTotalSegments << " uses isDiff " << segment.isDiff() << " (Polarization x"
              << " Cross Section), which is no longer supported; use an Analyzing Power (isDiff 7)"
              << " segment on the inverse channel instead.";
      return stop();
    }
    std::string why;
    SpecList<ExtrapolationSpec>::Item item;
    item.spec = ExtrapolationOf(segment, numTotalSegments, why);
    if (!why.empty()) {
      message << "ERROR: Test segment #" << numTotalSegments << ": " << why << ".";
      return stop();
    }
    item.warning = warning;
    if (segment.isActive() == 1 && segment.thicknessFitIgnored())
      item.note = "Note: test segment #" + std::to_string(numTotalSegments) +
                  " carries the fit form of its transmission / yield block or a thicknessgroup; a test segment "
                  "has no fit parameter and uses its own sample thickness n.";
    out.items.push_back(item);
  }
  if (line != "</segmentsTest>") out.failed = true;
  return out;
}

TargetEffectSpec TargetEffectSpecFromOldLine(std::istream &stream, bool &readFailed) {
  TargetIntLine line;
  ReadTargetIntLine(stream, line);
  readFailed = (stream.rdstate() & (std::stringstream::failbit | std::stringstream::badbit)) != 0;
  TargetEffectSpec t;
  // A line that ends after the parameter count was read as an inactive effect with every default (no message).
  if (!line.complete) return t;
  t.active = line.isActive == 1;
  t.datasetKeys = TargetIntSegmentKeys(TargetIntUnquoted(line.segmentList));
  t.numIntegrationPoints = line.numIntegrationPoints;
  t.convolution = line.isConvolution == 1;
  t.sigma = line.sigma;
  t.targetIntegration = line.isTargetIntegration == 1;
  t.density = line.density;
  t.stoppingPowerEquation = TargetIntUnquoted(line.stoppingPowerEq);
  t.stoppingPowerParameters = line.parameters;
  t.qCoefficients = line.isQCoefficients == 1;
  t.qCoefficientValues = line.qCoefficients;
  t.convolutionCoefficients = line.isConvCoefficients == 1;
  t.convolutionEquation = TargetIntUnquoted(line.convolutionEq);
  t.convolutionCoefficientValues = line.convCoefficients;
  // The old reader stopped filling the coefficients at the end of the line, and the old engine then read
  // numConvCoefficients of them, past the end of a short list (undefined); the missing ones are 0 here.
  t.convolutionCoefficientValues.resize(line.numConvCoefficients > 0 ? line.numConvCoefficients : 0, 0.0);
  t.hasStragglingFlag = line.hasStragglingFlag;
  t.straggling = line.isStraggling == 1;
  t.hasStragglingCoefficient = line.hasStragglingCoefficient;
  t.stragglingCoefficient = line.stragglingCoefficient;
  t.hasResonanceWidthMultiplier = line.hasResonanceWidthMultiplier;
  t.resonanceWidthMultiplier = line.resonanceWidthMultiplier;
  t.hasPointsPerWidth = line.hasPointsPerWidth;
  t.pointsPerWidth = line.pointsPerWidth;
  t.beamProfile = line.isBeamProfile;
  for (const TargetIntLine::ProfileComponent &p : line.beamProfile) {
    TargetEffectSpec::ProfileComponent c;
    c.xi = p.xi;
    c.omega = p.omega;
    c.alpha = p.alpha;
    c.weight = p.weight;
    t.beamProfileComponents.push_back(c);
  }
  t.beamTpcSigma = line.beamTpcSigma;
  t.beamTruncation = line.beamTruncation;
  t.beamDetailedBalance = line.beamDetailedBalance == 1;
  t.sigmaLab = line.sigmaLab;
  t.screening = line.hasScreening;
  t.screeningMalformed = line.screeningError;
  t.screeningUe = line.screeningNominal;
  t.screeningVary = line.screeningVary != 0;
  t.screeningSigma = line.screeningSigma;
  t.densityFit = line.densityFitRead;
  t.densityVary = line.densityVary != 0;
  t.densitySigmaText = line.densitySigmaToken;
  t.densityMax = line.densityMaxRead;
  if (line.densityFitError) t.densityReadError = line.densityFitErrorText;
  t.densityGroup = line.densityGroup;
  t.doppler = line.isDoppler;
  t.dopplerTemperature = line.dopplerTemperature;
  t.dopplerFit = line.dopplerFitRead;
  t.dopplerVaryText = line.dopplerVaryToken;
  t.dopplerSigmaText = line.dopplerSigmaToken;
  t.dopplerMaxText = line.dopplerMaxToken;
  if (line.dopplerError) t.dopplerReadError = line.dopplerErrorText;
  t.udr = line.isUdr;
  t.udrFile = line.udrFile;
  t.udrFlightPath = line.udrFlightPath;
  t.udrBurstFwhm = line.udrBurstFwhm;
  t.udrChannelWidth = line.udrChannelWidth;
  t.udrCentred = line.udrCentred != 0;
  return t;
}

SpecList<TargetEffectSpec> TargetEffectSpecsFromOldLines(std::istream &in) {
  SpecList<TargetEffectSpec> out;
  std::string line = "";
  getline(in, line);  // <targetInt>
  line = "";
  while (line != "</targetInt>" && !in.eof()) {
    getline(in, line);
    if (Blank(line)) continue;
    if (line == "</targetInt>" || in.eof()) continue;
    std::istringstream stm(line);
    SpecList<TargetEffectSpec>::Item item;
    item.spec = TargetEffectSpecFromOldLine(stm, item.readFailed);
    item.spec.source = line;
    out.items.push_back(item);
  }
  if (line != "</targetInt>") out.failed = true;
  return out;
}

experiment::AbundanceSpecs AbundanceSpecsFromOldLines(std::istream *section) {
  experiment::AbundanceSpecs out;
  if (!section) return out;
  std::istream &in = *section;
  std::string line = "";
  getline(in, line);  // <abundances>
  out.abundances.present = true;
  line = "";
  while (!in.eof() && line != "</abundances>") {
    getline(in, line);
    if (line == "</abundances>") break;
    const AbundanceLine entry = ReadAbundanceLine(line);
    if (entry.kind == AbundanceLine::Blank) continue;
    if (entry.kind == AbundanceLine::BadSum) {
      out.failed = true;
      out.error = "ERROR: <abundances> line \"" + line + "\" is not \"[<sample>:] sum <sigma>\".";
      return out;
    }
    if (entry.kind == AbundanceLine::BadEntry) {
      out.failed = true;
      out.error = "ERROR: <abundances> line \"" + line +
                  "\" is not \"[<sample>:] <Z> <A> <abundance> <vary> <sigma>\" or \"[<sample>:] sum <sigma>\".";
      return out;
    }
    experiment::AbundanceRowSpec row;
    row.sample = entry.sample;
    if (entry.kind == AbundanceLine::Sum) {
      row.isSum = true;
      row.sumSigma = entry.sumSigma;
    } else {
      row.z = entry.z;
      row.a = entry.a;
      row.value = entry.value;
      row.vary = entry.vary == 1;
      row.sigma = entry.sigma;
    }
    out.abundances.rows.push_back(row);
  }
  if (line != "</abundances>") out.failed = true;
  return out;
}

std::string OldSegmentParameterName(int key, bool isNorm) {
  const ParameterNames n = OldParameterNames(key, 0);
  return isNorm ? n.norm : n.energyShift;
}

std::string PriorCentresFromOldRows(std::istream *section, SpecList<DatasetSpec> &datasets) {
  if (!section) return std::string();
  const int numSegmentLines = (int)datasets.items.size();
  std::istream &in = *section;
  std::string line;
  bool inside = false;
  while (getline(in, line)) {
    const size_t a = line.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) continue;
    const size_t b = line.find_last_not_of(" \t\r\n");
    const std::string trimmed = line.substr(a, b - a + 1);
    if (trimmed == "<parameterSettings>") {
      inside = true;
      continue;
    }
    if (trimmed[0] == '<') {
      inside = false;
      continue;
    }
    if (!inside || trimmed[0] == '#') continue;
    std::istringstream stm(trimmed);
    std::vector<std::string> tok;
    std::string t;
    while (stm >> t) tok.push_back(t);
    if (tok.size() != 3 || tok[1] != "prior_centre") continue;

    // The row names its parameter by the whole name (the old engine matched it against its parameter names).
    const std::string &name = tok[0];
    int key = 0;
    bool isNorm = false;
    for (int k = 1; k <= numSegmentLines && key == 0; k++) {
      if (name == OldSegmentParameterName(k, true)) {
        key = k;
        isNorm = true;
      } else if (name == OldSegmentParameterName(k, false)) {
        key = k;
      }
    }
    if (key == 0 && numSegmentLines < 1)
      return "ERROR: <parameterSettings> row \"" + trimmed + "\": a prior_centre row moves the prior of a data "
             "segment, and the project has no <segmentsData> line.";
    if (key == 0)
      return "ERROR: <parameterSettings> row \"" + trimmed + "\" names no parameter a prior_centre row can move: it "
             "takes " + OldSegmentParameterName(1, true) + " ... " + OldSegmentParameterName(numSegmentLines, true) +
             " or " + OldSegmentParameterName(1, false) + " ... " + OldSegmentParameterName(numSegmentLines, false) +
             " (a <segmentsData> line's key, counting inactive lines).  Delete the row, or rename it after the data "
             "line it means (segment_<N>_norm, N the line's place in <segmentsData>).";
    char *end = NULL;
    const double centre = strtod(tok[2].c_str(), &end);
    if (end == tok[2].c_str() || *end != '\0' || !std::isfinite(centre))
      return "ERROR: <parameterSettings> row \"" + trimmed +
             "\" is not \"segment_N_norm prior_centre value\" or \"segment_N_energy_shift prior_centre value\".";
    if (isNorm && centre <= 0.)
      return "ERROR: <parameterSettings> row \"" + trimmed + "\": a normalization prior centre must be positive.";
    DatasetSpec &d = datasets.items[key - 1].spec;  // keys are 1, 2, ... in item order
    if (isNorm) {
      d.normCentre.written = true;
      d.normCentre.value = centre;
      d.normCentre.source = trimmed;
    } else {
      d.shiftCentre.written = true;
      d.shiftCentre.value = centre;
      d.shiftCentre.source = trimmed;
    }
  }
  return std::string();
}

}  // namespace convert
}  // namespace azure2
