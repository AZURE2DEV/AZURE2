// The run's output files (WorkingCopy::WriteOutputFiles): AZUREOut_* with their .band siblings and the sample
// decompositions, chiSquared.out, normalizations.out and shifts.out; thm_experiments.out through
// WorkingCopy::WriteThmExperiments.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "AZUREOutput.h"
#include "Config.h"
#include "CovarianceBand.h"
#include "WorkingCopy.h"

/*!
 * The decomposition of a sample segment, one file per segment since
 * AZUREOut files are named by pair keys only: the energy, the total the fit
 * compares with the data, and the abundance-weighted term of every compound
 * nucleus (azure2::EvaluateModelPoint).  Same frame as the
 * AZUREOut file (kinflag.dat), one line per point.
 */
namespace {
/// The cross-section and S-factor columns of an AZUREOut row from a model or data value v: sigma (b) and S = sigma E
/// exp(2 pi eta + g E) (MeV b; the point's factor times exp(g E), g = 0 but for the modified S factor, whose rows give
/// S* in this column).  For a data set or extrapolation of the S factor, v = c S (c the scale's constant,
/// Segment::ModelScaleFactor), so sigma = (v / c) / (E exp(2 pi eta + g E)) and S = v / c; otherwise v is the cross
/// section (times the scale's constant of an extrapolation in another unit, as before) and S = v E exp(2 pi eta).
struct SigmaAndS {
  bool modelIsS;
  double scale;
  double shape_perMeV;
  double Conversion(const Unit &point) const {
    const double c = point.GetSFactorConversion();
    return shape_perMeV == 0.0 ? c : c * std::exp(shape_perMeV * point.GetCMEnergy());
  }
  double Sigma(double v, const Unit &point) const { return modelIsS ? v / scale / Conversion(point) : v; }
  double S(double v, const Unit &point) const { return modelIsS ? v / scale : v * Conversion(point); }
};
}  // namespace

static void WriteNaturalDecomposition(WorkingCopy &data, const Config &configure, const azure2::Results &results, bool isExtrap,
                                      int kinflag, const std::vector<std::string> &nucleusLabels, const Segment *segment) {
  std::string name = configure.outputdir + "AZUREOut_sample_seg=" + std::to_string(segment->GetSegmentKey()) +
                     (isExtrap ? ".extrap" : ".out");
  std::ofstream out(configure.OpenPath(name).c_str());
  if (!out) return;
  out << "# Segment " << segment->GetSegmentKey() << ": abundance-weighted sum over compound nuclei (sample "
      << segment->GetNaturalToken() << ")" << std::endl;
  if (segment->IsTransmission())
    out << "# transmission of n = " << segment->GetSampleThickness(data.State(*segment)) << " atoms/b: the total column is exp(-n s), the"
        << " nucleus columns the abundance-weighted cross sections s_k (b), s = sum_k s_k" << std::endl;
  out << "# E_" << (kinflag ? "lab" : "cm") << " [MeV]   total";
  for (size_t k = 0; k < nucleusLabels.size(); k++) out << "   " << nucleusLabels[k];
  out << std::endl;
  out.precision(10);
  for (int i = 1; i <= segment->NumPoints(); i++) {
    Unit *point = segment->GetPoint(data, i);
    double kin = (kinflag != 0) ? Segment::Jacobian(data.State(*segment), 0, i - 1) : 1.0;
    if (kin == 0.0) kin = 1.0;
    out << std::setw(18) << std::scientific << (kinflag != 0 ? point->GetLabEnergy() : point->GetCMEnergy())
        << std::setw(18) << std::scientific << point->GetModelValue(results) / kin;
    const std::vector<double> &terms = point->GetNaturalTerms(results);
    for (size_t k = 0; k < nucleusLabels.size(); k++)
      out << std::setw(18) << std::scientific << (k < terms.size() ? terms[k] / kin : 0.0);
    out << std::endl;
  }
}

/*!
 * Writes the output files for the calculation.  The output files are all in center of mass frame, and contain
 * columns for energy, angle, calculated cross section, calculated s-factor, experimental cross section and error
 * and experimental s-factor and error.
 */

// One ".band" line: energy, excitation, angle, xs, d(xs), S-factor, d(S-factor).
static void WriteBandLine(std::ostream &o, double energy, double excitation,
                          double angle, double xs, double dxs, double conv) {
  o << std::setw(18) << std::scientific << energy
    << std::setw(18) << std::scientific << excitation
    << std::setw(18) << std::scientific << angle
    << std::setw(18) << std::scientific << xs
    << std::setw(18) << std::scientific << dxs
    << std::setw(18) << std::scientific << xs * conv
    << std::setw(18) << std::scientific << dxs * conv
    << std::endl;
}

void WorkingCopy::WriteOutputFiles(const Config &configure, azure2::Results &results, bool isFit, const BandData *band) {
  // use kinflag.dat to control output option
  std::ifstream kinoption;
  kinoption.open(configure.OpenPath("kinflag.dat"));
  int kinflag = 0;  // 0 cm output, 1 lab output
  if (kinoption) {
    kinoption >> kinflag;
    kinoption.close();
  }
  if (kinflag != 0) configure.outStream << "Using alternate output format..." << std::endl;

  // THM experiments and profiled THM norms: their profile from the model values in `results` (a snapshot written
  // during a fit has not been profiled), so that the norm and the background written are those of these models.
  if (configure.paramMask & Config::CALCULATE_WITH_DATA) {
    for (int g = 0; g < Data().NumThmGroups(); g++)
      if (!set_->thmGroups_[g].trivial) ProfileThmGroup(g, results);
    for (int i = 1; i <= Data().NumSegments(); i++)
      if (Data().ThmGroupOf(i) < 0 && Data().GetSegment(i)->IsProfiledNorm())
        results.objective.segmentChi2[i] = Data().GetSegment(i)->ProfileNormChiSquared(*this, results, results.objective.norm[i]);
  }

  WriteSegmentFiles(configure, results, kinflag, band);
  if (!isFit && (configure.paramMask & Config::CALCULATE_WITH_DATA)) WriteChiSquaredFile(configure, results);
  if (!isFit && (configure.paramMask & Config::CALCULATE_WITH_DATA) && Data().NumThmGroups() > 0)
    WriteThmExperiments(configure, results);
  WriteNormalizationsFile(configure, results);
  WriteShiftsFile(configure, results);
}

void WorkingCopy::WriteSegmentFiles(const Config &configure, const azure2::Results &results, int kinflag,
                              const BandData *band) {
  // The files are opened under the output directory as resolved (Config::OpenPath), a prefix ending in '/'.
  std::string outputPrefix = configure.OpenPath(configure.outputdir);
  if (!outputPrefix.empty() && outputPrefix.back() != '/') outputPrefix += '/';
  AZUREOutput output(outputPrefix);
  if (!(configure.paramMask & Config::CALCULATE_WITH_DATA)) output.SetExtrap();

  // When a covariance is available, write a sibling ".band" file per output file,
  // with the same block/point structure so the GUI can pair them.
  bool writeBand = band && (configure.paramMask & Config::CALCULATE_COVARIANCE_BAND) && !band->grad.empty();
  std::map<std::string, std::ofstream> bandFiles;
  auto bandStreamFor = [&](int aa, int ir) -> std::ofstream * {
    std::string name = configure.outputdir + "AZUREOut_aa=" + std::to_string(aa);
    if (ir == -1)
      name += "_TOTAL_CAPTURE";
    else
      name += "_R=" + std::to_string(ir);
    name += output.IsExtrap() ? ".extrap.band" : ".out.band";
    std::ofstream &f = bandFiles[name];
    if (!f.is_open()) {
      f.open(configure.OpenPath(name).c_str());
      f.precision(10);
    }
    return f.is_open() ? &f : nullptr;
  };
  // Look up a point's parameter-sensitivity row, or nullptr if absent.
  auto gradFor = [&](Unit *p) -> const std::vector<double> * {
    if (!band) return nullptr;
    std::map<Unit *, std::vector<double>>::const_iterator it = band->grad.find(p);
    return (it == band->grad.end()) ? nullptr : &it->second;
  };

  SegmentIterator firstSumIterator = Data().GetSegments().end();
  for (SegmentIterator segment = Data().GetSegments().begin();
       segment < Data().GetSegments().end(); segment++) {
    if (segment->IsTotalCapture()) {
      firstSumIterator = segment;
      // With component segments, no need to skip - total capture is handled within the segment
    }
    int aa = segment->GetEntranceKey();
    int ir = segment->GetExitKey();
    std::filebuf *buf;
    if (firstSumIterator != Data().GetSegments().end())
      buf = output(aa, -1);
    else {
      if (segment->IsAngularDist() &&
          !(configure.paramMask & Config::CALCULATE_WITH_DATA))
        buf = output(aa, ir, true);
      else
        buf = output(aa, ir);
    }
    std::ostream out(buf);
    // A THM experiment with a background: the fitted curve is the model plus b(E).
    const int segmentIndex = (int)(segment - Data().GetSegments().begin()) + 1;
    const bool thmBackground = !output.IsExtrap() && Data().ThmGroupOf(segmentIndex) >= 0 &&
                               set_->thmGroups_[Data().ThmGroupOf(segmentIndex)].terms > 0;
    SegmentIterator thisSegment = segment;
    if (firstSumIterator != Data().GetSegments().end()) thisSegment = firstSumIterator;
    const SigmaAndS columns{segment->ModelIsSFactor(), segment->ModelScaleFactor(), segment->ModelSFactorShape()};

    // Band stream for this segment (skip angular-coefficient extrap files, which
    // have no scalar cross section; other segments get a block for alignment).
    bool bandThisSegment = writeBand && !(segment->IsAngularDist() && output.IsExtrap());
    std::ofstream *bandOut = bandThisSegment ? bandStreamFor(aa, (firstSumIterator != Data().GetSegments().end()) ? -1 : ir) : nullptr;
    // Point gradient, summed across total-capture components like the value is.
    auto pointBandGrad = [&](UnitIterator point) -> std::vector<double> {
      std::vector<double> g;
      const std::vector<double> *g0 = gradFor(&*point);
      if (g0) g = *g0;
      if (firstSumIterator != Data().GetSegments().end()) {
        int pointIndex = point - segment->GetPoints(*this).begin() + 1;
        for (SegmentIterator it = firstSumIterator; it < segment; it++) {
          const std::vector<double> *gi = gradFor(it->GetPoint(*this, pointIndex));
          if (!gi) continue;
          if (g.empty())
            g = *gi;
          else
            for (size_t k = 0; k < g.size() && k < gi->size(); k++) g[k] += (*gi)[k];
        }
      }
      return g;
    };

    if (kinflag == 0) {
      for (UnitIterator point = segment->GetPoints(*this).begin(); point < segment->GetPoints(*this).end(); point++) {
        out.precision(10);
        if (segment->IsAngularDist()) {
          out << std::setw(18) << std::scientific << point->GetCMEnergy();
          for (double c : segment->LegendreCoefficients(*point, results)) out << std::setw(18) << c;
          out << std::endl;
          if (bandOut) WriteBandLine(*bandOut, point->GetCMEnergy(), point->GetExcitationEnergy(),
                                     point->GetCMAngle(), 0., 0., 0.);
        } else {
          double fitCrossSection = point->GetModelValue(results);
          if (thmBackground) fitCrossSection += ThmBackgroundAt(segmentIndex, point->GetCMEnergy(), results);
          if (firstSumIterator != Data().GetSegments().end()) {
            int pointIndex = point - segment->GetPoints(*this).begin() + 1;
            for (SegmentIterator it = firstSumIterator; it < segment; it++)
              fitCrossSection += it->GetPoint(*this, pointIndex)->GetModelValue(results);
          }
          out << std::setw(18) << std::scientific << point->GetCMEnergy()
              << std::setw(18) << std::scientific << point->GetExcitationEnergy()
              << std::setw(18) << std::scientific << point->GetCMAngle()
              << std::setw(18) << std::scientific << columns.Sigma(fitCrossSection, *point)
              << std::setw(18) << std::scientific << columns.S(fitCrossSection, *point);
          if (!output.IsExtrap()) {
            double dataNorm = NormUsed((int)(thisSegment - Data().GetSegments().begin()) + 1, results);
            const int row = (int)(point - segment->GetPoints(*this).begin());
            out << std::setw(18) << std::scientific << columns.Sigma(segment->CMValue(State(*segment), row) * dataNorm, *point)
                << std::setw(18) << std::scientific << columns.Sigma(segment->CMError(State(*segment), row) * dataNorm, *point)
                << std::setw(18) << std::scientific << columns.S(segment->CMValue(State(*segment), row) * dataNorm, *point)
                << std::setw(18) << std::scientific << columns.S(segment->CMError(State(*segment), row) * dataNorm, *point)
                << std::endl;
          } else
            out << std::endl;
          if (bandOut) {
            std::vector<double> g = pointBandGrad(point);
            double dxs = g.empty() ? 0. : band->dXS(g);
            WriteBandLine(*bandOut, point->GetCMEnergy(), point->GetExcitationEnergy(),
                          point->GetCMAngle(), fitCrossSection, dxs, columns.Conversion(*point));
          }
        }
      }
    }
    if (kinflag == 1) {
      for (UnitIterator point = segment->GetPoints(*this).begin(); point < segment->GetPoints(*this).end(); point++) {
        out.precision(10);
        if (segment->IsAngularDist()) {
          out << std::setw(18) << std::scientific << point->GetCMEnergy();
          for (double c : segment->LegendreCoefficients(*point, results)) out << std::setw(18) << c;
          out << std::endl;
          if (bandOut) WriteBandLine(*bandOut, point->GetLabEnergy(), point->GetExcitationEnergy(),
                                     point->GetLabAngle(), 0., 0., 0.);
        } else {
          double fitCrossSection = (point->GetModelValue(results) +
                                    (thmBackground ? ThmBackgroundAt(segmentIndex, point->GetCMEnergy(), results) : 0.0)) /
                                   Segment::Jacobian(State(*segment), 0, (int)(point - segment->GetPoints(*this).begin()));
          if (firstSumIterator != Data().GetSegments().end()) {
            int pointIndex = point - segment->GetPoints(*this).begin() + 1;
            for (SegmentIterator it = firstSumIterator; it < segment; it++)
              fitCrossSection += it->GetPoint(*this, pointIndex)->GetModelValue(results);
          }
          out << std::setw(18) << std::scientific << point->GetLabEnergy()
              << std::setw(18) << std::scientific << point->GetExcitationEnergy()
              << std::setw(18) << std::scientific << point->GetLabAngle()
              << std::setw(18) << std::scientific << columns.Sigma(fitCrossSection, *point)
              << std::setw(18) << std::scientific << columns.S(fitCrossSection, *point);
          if (!output.IsExtrap()) {
            double dataNorm = NormUsed((int)(thisSegment - Data().GetSegments().begin()) + 1, results);
            const int row = (int)(point - segment->GetPoints(*this).begin());
            out << std::setw(18) << std::scientific << columns.Sigma(segment->LabValue(row) * dataNorm, *point)
                << std::setw(18) << std::scientific << columns.Sigma(segment->LabError(row) * dataNorm, *point)
                << std::setw(18) << std::scientific << columns.S(segment->LabValue(row) * dataNorm, *point)
                << std::setw(18) << std::scientific << columns.S(segment->LabError(row) * dataNorm, *point)
                << std::endl;
          } else
            out << std::endl;
          if (bandOut) {
            std::vector<double> g = pointBandGrad(point);
            double kin = Segment::Jacobian(State(*segment), 0, (int)(point - segment->GetPoints(*this).begin()));
            double dxs = g.empty() ? 0. : band->dXS(g) / kin;
            WriteBandLine(*bandOut, point->GetLabEnergy(), point->GetExcitationEnergy(),
                          point->GetLabAngle(), fitCrossSection, dxs, columns.Conversion(*point));
          }
        }
      }
    }
    out << std::endl
        << std::endl;
    out.flush();
    if (segment->IsNatural())
      WriteNaturalDecomposition(*this, configure, results, output.IsExtrap(), kinflag, set_->nucleusLabels_, &*segment);
    if (bandOut) {
      *bandOut << std::endl
               << std::endl;
      bandOut->flush();
    }
    firstSumIterator = Data().GetSegments().end();
  }
}

void WorkingCopy::WriteChiSquaredFile(const Config &configure, const azure2::Results &results) {
  std::ofstream chiOut;
  std::string chiOutFile = configure.outputdir + "chiSquared.out";
  chiOut.open(configure.OpenPath(chiOutFile).c_str());
  chiOut << "Segment#," << " Chi-Squared, " << " N, " << " Norm, " << " Norm-Chi-Squared " << std::endl;

  double totalChiSquared = 0.;
  double totalNormChiSquared = 0.;
  double totalN = 0.;
  for (SegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++) {
    SegmentIterator thisSegment = segment;
    // The normalization penalty as the fit sees it (AZURECalc::operator()):
    // the deviation from the nominal normalization in units of its uncertainty,
    // which the segment stores as a percentage.  This used to add a literal 1
    // per segment, so the reported total was just the number of segments.
    // A fixed norm pays no prior (docs/upgrade/objective.md D5), nor does a profiled THM norm, which is not a
    // parameter: the fit adds no prior for it (D4).
    double normChiSquared = 0.;
    double normError = thisSegment->GetNominalNorm() / 100. * thisSegment->GetNormError();
    if (normError != 0. && thisSegment->IsVaryNorm() && !thisSegment->IsProfiledNorm() &&
        !thisSegment->IsNormPriorReplaced()) normChiSquared =
                             pow((State(*thisSegment).norm - thisSegment->GetNominalNorm()) / normError, 2.0);
    const int thisIndex = (int)(thisSegment - Data().GetSegments().begin()) + 1;
    totalChiSquared += results.objective.segmentChi2[thisIndex];
    totalNormChiSquared += normChiSquared;
    totalN += thisSegment->NumPoints();
    chiOut << thisSegment->GetSegmentKey()
           << ","
           << results.objective.segmentChi2[thisIndex]
           << ","
           << thisSegment->NumPoints()
           << ","
           << NormUsed(thisIndex, results)
           << ","
           << normChiSquared
           << std::endl;
  }

  chiOut << "Total-Chi-Squared: "
         << totalChiSquared
         << " Total-Norm-Chi-Squared: "
         << totalNormChiSquared
         << " Total-N: "
         << totalN;
  // Park formalism: the J > 0 wall, when a level sits beyond it.
  if (configure.paramMask & Config::USE_PARK_FORMALISM)
    chiOut << " Total-Park-Chi-Squared: " << results.objective.park;
  // The energy-shift penalties the fit minimised (Segment::AddShiftPenalty), when
  // any shift is a fit parameter.  Without this the file reported less than the
  // minimiser's own total by exactly this amount.
  bool anyShiftVaried = false;
  double totalShiftChiSquared = 0.;
  for (SegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++) {
    if (segment->IsVaryEnergyShift() || segment->IsVaryEnergyShiftSqrt()) anyShiftVaried = true;
    if (segment->OwnsShiftFn())
      for (int q = 0; q < segment->NumShiftFnParams(); q++)
        if (segment->IsVaryShiftFn(q)) anyShiftVaried = true;
    segment->AddShiftPenalty(State(*segment), totalShiftChiSquared);
  }
  if (anyShiftVaried) chiOut << " Total-Shift-Chi-Squared: " << totalShiftChiSquared;
  // Isotopic abundances: their priors and the sum-to-one term.
  if (results.objective.abundance != 0.0)
    chiOut << " Total-Abundance-Chi-Squared: " << results.objective.abundance;
  // Electron-screening priors, when a screening potential has one.
  {
    bool anyPrior = false;
    for (int k = 1; k <= Data().NumScreeningParams(); k++)
      if (Data().GetScreeningTargetEffect(k)->GetScreeningSigma() > 0.0) anyPrior = true;
    if (anyPrior) chiOut << " Total-Screening-Chi-Squared: " << this->ScreeningPenalty();
  }
  // Fitted target densities: their priors, and the densities themselves (which param.sav also holds).
  if (Data().NumDensityParams() > 0) {
    bool anyPrior = false;
    for (int k = 1; k <= Data().NumDensityParams(); k++)
      if (Data().GetDensityTargetEffect(k)->GetDensitySigma() > 0.0) anyPrior = true;
    if (anyPrior) chiOut << " Total-Density-Chi-Squared: " << this->DensityPenalty();
  }
  // Apparatus parameters (stopping scales, a kernel's width): one line per row with a fitted prior.
  WriteApparatusChiSquared(*this, chiOut);
  // Fitted Doppler temperatures: their priors.
  if (Data().NumDopplerParams() > 0) {
    bool anyPrior = false;
    for (int k = 1; k <= Data().NumDopplerParams(); k++)
      if (Data().GetDopplerTargetEffect(k)->GetDopplerSigma() > 0.0) anyPrior = true;
    if (anyPrior) chiOut << " Total-Doppler-Chi-Squared: " << this->DopplerPenalty();
  }
  // Fitted sample thicknesses: their priors.
  if (Data().NumThicknessParams() > 0) {
    bool anyPrior = false;
    for (int k = 1; k <= Data().NumThicknessParams(); k++) {
      const Segment *seg = &set_->segments_[Data().GetThicknessSegment(k) - 1];
      if (seg->GetThicknessSigma() > 0.0 && !seg->IsThicknessGroupFollower()) anyPrior = true;
    }
    if (anyPrior) chiOut << " Total-Thickness-Chi-Squared: " << this->ThicknessPenalty();
  }
  // Nuisance priors (<parameterSettings> useAsNuisance), part of the objective in every path (objective.md D6).
  if (results.objective.nuisance != 0.0) chiOut << " Total-Nuisance-Chi-Squared: " << results.objective.nuisance;
  // The model data (T9): the <modelData> data (modelData.out).
  if (results.objective.hasModelData) chiOut << " Total-Model-Data-Chi-Squared: " << results.objective.modelData;
  // Isospin levels: the priors on the internal Coulomb shifts and isospin-forbidden amplitudes, and the soft ties.
  if (results.objective.hasIsospin) chiOut << " Total-Isospin-Chi-Squared: " << results.objective.isospin;
  // Prior parameter covariance (<config> priorCovariance=): (p - mu)^T C^-1 (p - mu), whenever the key is set.
  // Written with 15 significant digits (a new field: the older ones keep their 6).
  if (!configure.priorCovariance.empty()) {
    const std::streamsize prec = chiOut.precision(15);
    chiOut << " Total-Prior-Covariance-Chi-Squared: " << results.objective.priorCovariance;
    chiOut.precision(prec);
  }
  chiOut << std::endl
         << std::endl;
  chiOut.flush();
  chiOut.close();
  if (Data().HasCovariance()) this->WriteCovarianceSummary(configure, results);
}

void WorkingCopy::WriteNormalizationsFile(const Config &configure, const azure2::Results &results) {
  bool isVaryNorm = false;
  for (SegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++)
    if (segment->IsVaryNorm()) isVaryNorm = true;
  if (isVaryNorm) {
    std::string outputfile = configure.outputdir + "normalizations.out";
    std::ofstream out(configure.OpenPath(outputfile).c_str());
    if (out) {
      out.precision(4);
      out << std::scientific;
      out << "segment_key_#, file_name, low_angle_bound, high_angle_bound, aframe, low_energy_bound, high_energy_bound, eframe, norm, shift" << std::endl;
      // One complete row per segment whose normalization is varied (a norm-group member keeps its own
      // flag).  Before 2026-10-27 only the first four fields were conditional, so every other segment
      // left a row of six fields without key, file or angles.
      for (SegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++) {
        if (!segment->IsVaryNorm()) continue;
        out << segment->GetSegmentKey() << ","
            << segment->GetDataFile() << ","
            << segment->GetMinAngle() << ","
            << segment->GetMaxAngle() << ",";
        if (segment->IsCMDifferential() == true) {
          out << "CM,";
        } else
          out << "Lab,";
        out << segment->GetMinEnergy() << ","
            << segment->GetMaxEnergy() << ",";
        if (segment->IsCMDifferential() == true) {
          out << "CM,";
        } else
          out << "Lab,";
        out << NormUsed((int)(segment - Data().GetSegments().begin()) + 1, results) << ","
            << State(*segment).energyShift << std::endl;
      }
      out.flush();
      out.close();
    } else
      configure.outStream << "Could not write normalization file." << std::endl;
  }
}

void WorkingCopy::WriteShiftsFile(const Config &configure, const azure2::Results &results) {
  // Check if any energy shifts (constant or sqrt(E) term) are varied
  bool isVaryEnergyShift = false;
  for (SegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++) {
    if (segment->IsVaryEnergyShift() || segment->IsVaryEnergyShiftSqrt() || segment->IsShiftGroupFollower() ||
        segment->HasShiftFn()) {
      isVaryEnergyShift = true;
      break;
    }
  }

  // Write energy shifts file
  if (isVaryEnergyShift) {
    std::string outputfile = configure.outputdir + "shifts.out";
    std::ofstream out(configure.OpenPath(outputfile).c_str());
    if (out) {
      out.precision(4);
      out << std::scientific;
      // One complete row per segment with a free shift term; sqrt_shift is the
      // coefficient b of E' = E + shift + b*sqrt(E/MeV), in MeV^1/2.
      out << "segment_key_#, file_name, low_angle_bound, high_angle_bound, aframe, low_energy_bound, high_energy_bound, eframe, norm, shift, sqrt_shift" << std::endl;
      for (SegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++) {
        if (!(segment->IsVaryEnergyShift() || segment->IsVaryEnergyShiftSqrt() || segment->IsShiftGroupFollower() ||
              segment->HasShiftFn()))
          continue;
        out << segment->GetSegmentKey() << ","
            << segment->GetDataFile() << ","
            << segment->GetMinAngle() << ","
            << segment->GetMaxAngle() << ",";
        if (segment->IsCMDifferential() == true) {
          out << "CM,";
        } else
          out << "Lab,";
        out << segment->GetMinEnergy() << ","
            << segment->GetMaxEnergy() << ",";
        if (segment->IsCMDifferential() == true) {
          out << "CM,";
        } else
          out << "Lab,";
        out << NormUsed((int)(segment - Data().GetSegments().begin()) + 1, results) << ","
            << State(*segment).energyShift << ","
            << State(*segment).energyShiftSqrt << std::endl;
      }
      // Shift functions: one line per function, after the per-segment rows.  The
      // lines start with '#' so that readers of the rows above skip them.
      bool header = false;
      double totalShiftChi = 0.;
      for (SegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++)
        segment->AddShiftPenalty(State(*segment), totalShiftChi);
      for (int pos = 0; pos < (int)set_->segments_.size(); pos++) {
        const Segment &owner = set_->segments_[pos];
        if (!owner.OwnsShiftFn() || !owner.HasShiftFn()) continue;
        if (!header) {
          out << "# shift functions: dE(x) = shift + f(x), x = original lab energy of a point (MeV), dE in MeV" << std::endl
              << "# segment_key_#, member_keys, f(x), shift, a0 a1 ..., x_low, x_high, dE(x_low), dE(x_high), penalty"
              << std::endl;
          header = true;
        }
        std::ostringstream keys;
        double lo = std::numeric_limits<double>::infinity(), hi = -lo;
        for (int m = 0; m < (int)set_->segments_.size(); m++) {
          if (set_->segments_[m].GetShiftFnSourcePos() != pos) continue;
          keys << (m == pos ? "" : " ") << set_->segments_[m].GetSegmentKey();
          for (int pid = 1; pid <= set_->segments_[m].NumPoints(); pid++) {
            Unit *pt = set_->segments_[m].GetPoint(*this, pid);
            const double x = pt->GetOriginalEnergy() > 0 ? pt->GetOriginalEnergy() : pt->GetLabEnergy();
            lo = std::min(lo, x);
            hi = std::max(hi, x);
          }
        }
        // The function's own penalty: the owner's shift penalties minus those of its
        // constant and sqrt(E) terms.
        double pen = 0.0, penOwn = 0.0;
        owner.AddShiftPenalty(State(owner), pen);
        if (owner.IsVaryEnergyShift() && owner.GetEnergyShiftError() != 0.0)
          penOwn = pow((State(owner).energyShift - owner.GetNominalEnergyShift()) / owner.GetEnergyShiftError(), 2.0);  // as AddShiftPenalty
        out << "# " << owner.GetSegmentKey() << "," << keys.str() << "," << owner.GetShiftFnExpression() << ","
            << State(owner).energyShift << ",";
        for (int q = 0; q < owner.NumShiftFnParams(); q++) out << (q ? " " : "") << State(owner).shiftFnValues[q];
        out << "," << lo << "," << hi << "," << owner.TotalEnergyShift(State(owner), lo) << "," << owner.TotalEnergyShift(State(owner), hi) << ","
            << pen - penOwn << std::endl;
      }
      if (header) out << "# Total-Shift-Chi-Squared: " << totalShiftChi << std::endl;
      out.flush();
      out.close();
    } else
      configure.outStream << "Could not write energy shifts file." << std::endl;
  }
}
