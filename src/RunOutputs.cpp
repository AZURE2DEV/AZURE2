// The output stages of a data or extrapolation run (ProgramRun.h): the final evaluation and its chi^2 lines, the
// band, the output files, scattering lengths and <modelData>.
#include "ProgramRun.h"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "CompoundSystem.h"
#include "Config.h"
#include "ModelData.h"
#include "InferenceTypes.h"
#include "Profile.h"
#include "ScatteringLength.h"
#include "Segment.h"
#include "WorkingCopy.h"

namespace {

// The free parameters as the band sees them.
ParamIndexMap FreeParameterMap(const ProgramRun &run) {
  const int nMn = run.params.GetMinuitParams().Params().size();
  std::vector<bool> fixedMask(nMn);
  for (int fi = 0; fi < nMn; fi++) fixedMask[fi] = run.params.GetMinuitParams().Parameter(fi).IsFixed();
  return BuildParamIndexMap(run.compound.get(), run.data, fixedMask);
}

}  // namespace

void FinalEvaluation(ProgramRun &run) {
  run.data->SetFit(false);
  run.data->SetErrorAnalysis(false);
  double chiSquared = 0.0;
  {
    AZUREProfile::ScopedReport profileReport(run.configure.profile, (run.configure.paramMask & Config::PERFORM_FIT) ? "final evaluation" : "calculation");
    chiSquared = (*run.objective)(run.params.GetMinuitParams().Params());
  }
  if (run.configure.paramMask & Config::CALCULATE_WITH_DATA) {
    run.configure.outStream << std::endl
                          << std::endl;
    for (SegmentIterator segment = run.data->Data().GetSegments().begin();
         segment < run.data->Data().GetSegments().end(); segment++) {
      run.configure.outStream << "Segment #"
                            << segment->GetSegmentKey()
                            << " Chi-Squared/N: "
                            << run.results.objective.segmentChi2[segment - run.data->Data().GetSegments().begin() + 1] / segment->NumPoints()
                            << std::endl;
      // With component segments, no need to skip - total capture is handled within the segment
    }
    // Twelve significant digits: a unit change of a prior must be readable
    // from this line against a data chi-squared of thousands.
    std::ostringstream total;
    total << std::setprecision(12) << chiSquared;
    run.configure.outStream << "Total Chi-Squared: "
                          << total.str() << std::endl
                          << std::endl;
  }
}

void WriteRunOutputs(ProgramRun &run) {
  // Parameter covariance for the cross-section band, saved for reuse by later
  // extrapolation runs.  Any fit (MIGRAD or LM) yields it; MINOS is not needed.
  const std::string covPath = run.configure.outputdir + "covariance.dat";
  // Write Output Files
  run.configure.outStream << "Writing output files..." << std::endl;
  // Band covariance: this run's fit if available, else one saved by a prior fit.
  BandData bandData;
  const BandData *bandPtr = nullptr;
  if (run.configure.paramMask & Config::CALCULATE_COVARIANCE_BAND) {
    // Optionally inflate to reduced chi^2 = 1 (Birge ratio), then persist.
    if (run.haveBandCov) {
      if ((run.configure.paramMask & Config::SCALE_COVARIANCE_BY_CHI2) &&
          (run.configure.paramMask & Config::CALCULATE_WITH_DATA)) {
        double chi2data = 0.;
        int nData = 0;
        for (int s = 1; s <= run.data->Data().NumSegments(); s++) {
          const Segment *seg = run.data->Data().GetSegment(s);
          if (!seg) continue;
          chi2data += run.results.objective.segmentChi2[s];
          nData += seg->NumPoints();
        }
        // Degrees of freedom use the total number of free fit parameters
        // (including norms/shifts), even though the covariance matrix itself
        // holds only the R-matrix sub-block (run.bandCov.size()).
        int nFree = 0;
        const int nMnAll = run.params.GetMinuitParams().Params().size();
        for (int i = 0; i < nMnAll; i++)
          if (!run.params.GetMinuitParams().Parameter(i).IsFixed()) nFree++;
        int ndf = nData - nFree;
        double redChi2 = (ndf > 0) ? chi2data / ndf : 1.0;
        double factor = (redChi2 > 1.0) ? redChi2 : 1.0;  // inflate only
        run.configure.outStream << "Calculating covariance scaling: chi-Squared = " << chi2data
                              << ", N = " << nData << ", free run.params = " << nFree
                              << ", nu = " << ndf << ", reduced chi-squared = " << redChi2 << "." << std::endl;
        if (factor > 1.0) {
          const int nMat = run.bandCov.size();
          for (int a = 0; a < nMat; a++)
            for (int b = 0; b < nMat; b++) run.bandCov.M[a][b] *= factor;
          run.bandCov.scale = factor;
          run.configure.outStream << "Covariance scaled by " << factor << "." << std::endl;
        } else {
          run.configure.outStream << "Reduced chi-squared <= 1; covariance left unscaled."
                                << std::endl;
        }
      }
      // Name the columns ("#columns" header, which LoadBandCovariance and numpy.loadtxt skip): the free level
      // energies and reduced-width amplitudes in parameter order, as param.sav holds them (pyazr.prior_from_fit).
      {
        const ROOT::Minuit2::MnUserParameters &mp = run.params.GetMinuitParams();
        std::vector<std::string> names;
        for (unsigned int i = 0; i < mp.Params().size(); i++) {
          const ParamKind kind = run.params.Kind((int)i);
          if (!mp.Parameter(i).IsFixed() && IsRMatrixParameter(kind))
            names.push_back(mp.GetName(i));
        }
        if ((int)names.size() == run.bandCov.size()) run.bandCov.names = names;
      }
      if (SaveBandCovariance(run.configure.OpenPath(covPath), run.bandCov))
        run.configure.outStream << "Saved parameter covariance for cross-section bands." << std::endl;
    }
    BandCovariance cov;
    bool fitRun = run.configure.paramMask & Config::PERFORM_FIT;  // fit or MINOS
    if (run.haveBandCov)
      cov = run.bandCov;
    else if (fitRun) {
      // The fit ran but produced no usable covariance (already reported above);
      // don't fall back to a stale saved file.
      if (BandPackedColumns(FreeParameterMap(run)).empty())
        run.configure.outStream << "Uncertainty band skipped: the fit has no free level parameter. The band "
                                 "propagates the covariance of level energies, widths and the THM coherent "
                                 "background; the parameters of targets and data sets are not in it."
                              << std::endl;
      else
        run.configure.outStream << "Uncertainty band skipped: this fit produced no usable covariance."
                              << std::endl;
    } else if (!LoadBandCovariance(run.configure.OpenPath(covPath), cov)) {
      // Extrapolation / no-fit run: the band needs a covariance from a prior fit.
      run.configure.outStream << "No saved parameter covariance found (" << covPath
                            << "); writing the normal output without an uncertainty band. "
                               "To get a band here, first run a fit (or MINOS) with the "
                               "\"Uncertainty band\" option enabled, then rerun this extrapolation."
                            << std::endl;
    }
    if (!cov.empty()) {
      if (BuildBandData(run.compound.get(), run.state, run.inputs, run.data, run.results, run.configure, cov, bandData, &run.fitFixedMask)) {
        bandPtr = &bandData;
        run.configure.outStream << "Writing cross-section uncertainty bands..." << std::endl;
      } else {
        run.configure.outStream << "Cross-section bands unavailable: the analytic model "
                                 "sensitivities are not supported for this configuration."
                              << std::endl;
      }
      // BuildBandData computes the per-point sensitivities on the real
      // compound/data objects (its gradient map is keyed by the real Unit*),
      // and that forward+adjoint solve overwrites each point's stored best-fit
      // cross section.  Re-run the best-fit forward pass so the plotted and
      // written segments are the best fit, not the last gradient probe; the
      // band itself lives in bandData and is unaffected.
      (*run.objective)(run.params.GetMinuitParams().Params());
    }
  }
  run.data->WriteOutputFiles(run.configure, run.results, false, bandPtr);
  // s-wave scattering lengths of every particle pair at the parameters just written: the table goes to
  // scatteringLengths.out only, the run log (GUI execution screen) names the file.
  try {
    std::string report = ScatteringLength::Report(run.compound.get(), run.state, run.configure);
    std::ofstream slOut(run.configure.OpenPath(run.configure.outputdir + "scatteringLengths.out").c_str());
    slOut << report;
    run.configure.outStream << "Scattering lengths written to " << run.configure.outputdir << "scatteringLengths.out"
                          << std::endl;
  } catch (const std::exception &e) {
    run.configure.outStream << "Scattering lengths not calculated: " << e.what() << std::endl;
  }
  // The <modelData> data at the same parameters (ModelData.h): modelData.out and the run log, only with the section.
  const std::string modelData = ModelData::Report(run.compound.get(), run.state, run.inputs, run.configure,
                                                    run.compound->GetModelData());
  if (!modelData.empty()) {
    std::ofstream mdOut(run.configure.OpenPath(run.configure.outputdir + "modelData.out").c_str());
    mdOut << modelData;
    run.configure.outStream << modelData;
  }
}
