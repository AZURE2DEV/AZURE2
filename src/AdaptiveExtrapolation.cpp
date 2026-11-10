// The extrapolation run's adaptive energy grids (AdaptiveExtrapolation.h).
#include "AdaptiveExtrapolation.h"

#include <cmath>
#include <filesystem>
#include <random>
#include <sstream>
#include <system_error>

#include "AZURECalc.h"
#include "AZUREParams.h"
#include "CNuc.h"
#include "Config.h"
#include "EData.h"
#include "ExperimentAdapter.h"
#include "Kinematics.h"
#include "Results.h"
#include "azure2/experiment/AdaptiveGrid.h"
#include "azure2/experiment/Frames.h"
#include "azure2/experiment/ProjectReader.h"
#include "azure2/physics/LevelParameters.h"

namespace {

using azure2::experiment::ExtrapolationSpec;
using azure2::experiment::ExtrapolationSpecs;

// The positions of the run's data parameters (EData::FillMnParams sets them from the length of the vector before):
// the compound nucleus's parameters, then the data set's, as the run laid them out.
void LayOutParameters(CNuc *compound, EData *data, const Config &configure) {
  AZUREParams layout;
  compound->FillMnParams(layout, &configure);
  data->FillMnParams(layout);
}

// The specs with every message cleared: the run's MakePoints printed them once.
ExtrapolationSpecs Quiet(ExtrapolationSpecs specs) {
  for (auto &item : specs.extrapolations.items) {
    item.warning.clear();
    item.note.clear();
  }
  for (auto &item : specs.targetEffects.items) {
    item.warning.clear();
    item.note.clear();
  }
  return specs;
}

// A scratch directory for the external-capture integrals of the sweeps, removed with the object.
class ScratchDirectory {
 public:
  ScratchDirectory() {
    std::error_code ec;
    const std::filesystem::path base = std::filesystem::temp_directory_path(ec);
    std::random_device random;
    for (int attempt = 0; attempt < 100 && !ec; attempt++) {
      const std::filesystem::path p = base / ("azure2_adaptive_" + std::to_string(random()));
      if (std::filesystem::create_directory(p, ec) && !ec) {
        path_ = p;
        return;
      }
    }
  }
  ~ScratchDirectory() {
    std::error_code ec;
    if (!path_.empty()) std::filesystem::remove_all(path_, ec);
  }
  bool Ok() const { return !path_.empty(); }
  std::string Dir() const { return path_.string() + "/"; }

 private:
  std::filesystem::path path_;
};

// The model values of one extrapolation at a batch of lab energies, one component per angle, by the run's point pass.
class SweepModel {
 public:
  SweepModel(const Config &sweepConfig, CNuc *compound, const vector_r &params, const ExtrapolationSpecs &specs,
             std::size_t index)
      : config_(sweepConfig), compound_(compound), params_(params), specs_(specs), index_(index) {}

  std::vector<std::vector<double>> operator()(const std::vector<double> &energies) const {
    ExtrapolationSpecs s = specs_;
    for (std::size_t i = 0; i < s.extrapolations.items.size(); i++) {
      ExtrapolationSpec &x = s.extrapolations.items[i].spec;
      x.adaptive.on = false;
      if (i == index_) {
        x.energies = energies;
      } else if (x.active && !x.energies.empty()) {
        x.energies.resize(1);
        x.angles.resize(1);
      }
    }
    const ExtrapolationSpec &mine = s.extrapolations.items[index_].spec;
    EData data;
    if (data.MakePoints(config_, compound_, s) == -1) throw std::runtime_error("the points of a sweep could not be made");
    LayOutParameters(compound_, &data, config_);
    if (data.Initialize(compound_, config_) == -1) throw std::runtime_error("the points of a sweep could not be initialized");
    data.SetFit(false);
    data.SetErrorAnalysis(false);
    azure2::Results results;
    AZURECalc calculate(&data, results, compound_, config_);
    calculate(params_);
    ESegment *segment = nullptr;
    for (int k = 1; k <= data.NumSegments(); k++)
      if (data.GetSegment(k)->GetSegmentKey() == mine.key) segment = data.GetSegment(k);
    const std::size_t nE = energies.size(), nA = mine.angles.size();
    if (!segment || (std::size_t)segment->NumPoints() != nE * nA) throw std::runtime_error("a sweep lost its points");
    std::vector<std::vector<double>> values(nE, std::vector<double>(nA));
    for (std::size_t a = 0; a < nA; a++)
      for (std::size_t e = 0; e < nE; e++)
        values[e][a] = segment->GetPoint((int)(a * nE + e) + 1)->GetModelValue(results);
    return values;
  }

 private:
  const Config &config_;
  CNuc *compound_;
  const vector_r &params_;
  const ExtrapolationSpecs &specs_;
  std::size_t index_;
};

// E_lambda + k Gamma_lambda (k = 0, +-1/2, +-1, +-2, +-4) of the levels of the entrance pair's compound nucleus, as lab
// energies in [lo, hi].
std::vector<double> Seeds(CNuc *compound, int entranceKey, double lo, double hi, bool observed) {
  std::vector<double> seeds;
  if (!compound->IsPairKey(entranceKey)) return seeds;
  const PPair *pair = compound->GetPair(compound->GetPairNumFromKey(entranceKey));
  if (pair->GetPType() != 0) return seeds;
  const azure2::experiment::Partition partition = PartitionOf(pair);
  const double threshold = pair->GetSepE() + pair->GetExE();
  for (int j = 1; j <= compound->NumJGroups(); j++) {
    JGroup *group = compound->GetJGroup(j);
    if (group->GetNucleus() != pair->GetNucleus() || !group->IsInRMatrix()) continue;
    for (int la = 1; la <= group->NumLevels(); la++) {
      const ALevel *level = group->GetLevel(la);
      if (!level->IsInRMatrix()) continue;
      double width = 0.0;
      if (observed)
        for (int ch = 1; ch <= group->NumChannels(); ch++)
          if (static_cast<azure2::physics::ObservedKind>(level->GetObservedKind(ch)) == azure2::physics::ObservedKind::Width)
            width += std::fabs(level->GetBigGamma(ch));
      const double eCM = (observed ? level->GetTransformE() : level->GetFitE()) - threshold;
      for (double k : {0.0, -0.5, 0.5, -1.0, 1.0, -2.0, 2.0, -4.0, 4.0}) {
        if (k != 0.0 && !(width > 0.0)) continue;
        const double c = eCM + k * width;
        if (!(c > 0.0)) continue;
        const double lab = azure2::experiment::BeamLabFromCM(partition, c, Kinematics::Relativistic());
        if (lab >= lo && lab <= hi) seeds.push_back(lab);
      }
    }
  }
  return seeds;
}

}  // namespace

bool HasAdaptiveGrids(const Config &configure) {
  if (!configure.format3) return false;
  for (const auto &item : configure.format3->built.extrapolations.extrapolations.items)
    if (item.spec.active && item.spec.adaptive.on) return true;
  return false;
}

bool SessionRefusesAdaptiveGrids(const Config &configure) {
  if (!HasAdaptiveGrids(configure)) return false;
  configure.outStream << "ERROR: an extrapolation with an adaptive grid: the engine's extrapolation run (menu 3) makes "
                         "that grid at its parameters; a pyazr session takes a static grid (range, list, file), and "
                         "pyazr.refine_grid refines one with a Python model."
                      << std::endl;
  return true;
}

std::unique_ptr<EData> MakeRefinedExtrapolations(const Config &configure, CNuc *compound, const vector_r &params) {
  ExtrapolationSpecs specs = Quiet(configure.format3->built.extrapolations);
  ScratchDirectory scratch;
  if (!scratch.Ok()) {
    configure.outStream << "ERROR: adaptive grids: no scratch directory for the external-capture integrals of the "
                           "sweeps could be made in the system's temporary directory."
                        << std::endl;
    return nullptr;
  }
  // The sweeps' configuration: integrals computed afresh in the scratch directory, no check files, no messages.
  Config sweepConfig(configure);
  sweepConfig.outputdir = scratch.Dir();
  sweepConfig.paramMask &= ~Config::USE_PREVIOUS_INTEGRALS;
  sweepConfig.paramMask |= Config::USE_API;
  sweepConfig.fileCheckMask = 0;
  sweepConfig.screenCheckMask = 0;
  const bool observed = compound->ObservedWidthsAt(configure, params);
  for (std::size_t i = 0; i < specs.extrapolations.items.size(); i++) {
    ExtrapolationSpec &x = specs.extrapolations.items[i].spec;
    if (!x.active || !x.adaptive.on) continue;
    std::vector<double> start = x.energies;
    const std::size_t seeded = start.size();
    for (double e : Seeds(compound, x.reaction.entranceKey, x.minE, x.maxE, observed)) start.push_back(e);
    azure2::experiment::AdaptiveGridResult grid;
    try {
      grid = azure2::experiment::RefineGrid(start, x.adaptive.tolerance, x.adaptive.maxPoints,
                                            SweepModel(sweepConfig, compound, params, specs, i));
    } catch (const std::exception &e) {
      configure.outStream << "ERROR: extrapolation #" << x.key << ": its adaptive grid could not be refined: " << e.what()
                          << "." << std::endl;
      return nullptr;
    }
    configure.outStream << "Adaptive grid of extrapolation #" << x.key << ": " << grid.x.size() << " energies ("
                        << seeded << " of the start grid, " << start.size() - seeded << " level seeds, " << grid.sweeps
                        << " sweeps, " << grid.evaluations << " evaluations)." << std::endl;
    if (grid.budgetReached) {
      std::ostringstream at;
      at.precision(8);
      at << grid.largestDeviationAt;
      configure.outStream << "WARNING: extrapolation #" << x.key << ": the adaptive grid reached max_points = "
                          << x.adaptive.maxPoints << " with intervals left to split; the largest deviation from "
                          << "linear interpolation left is " << grid.largestDeviation << " (relative; tolerance "
                          << x.adaptive.tolerance << ") at E_lab = " << at.str() << " MeV." << std::endl;
    }
    x.energies = grid.x;
    x.adaptive.on = false;
  }
  std::unique_ptr<EData> data(new EData);
  if (data->MakePoints(configure, compound, specs) == -1) return nullptr;
  LayOutParameters(compound, data.get(), configure);
  return data;
}
