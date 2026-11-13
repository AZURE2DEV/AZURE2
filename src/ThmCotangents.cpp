// THM data in the scalar gradient pass: the model cotangents of profiled THM data (include/ThmCotangents.h).
#include "ThmCotangents.h"

#include <cstddef>
#include <vector>

#include "AZUREGrad.h"
#include "Config.h"
#include "ModelPass.h"
#include "Results.h"
#include "Segment.h"
#include "THMMatrixFunc.h"
#include "ThmProfile.h"
#include "Unit.h"
#include "WorkingCopy.h"

namespace azure2 {

bool ThmThroughEnergyTasks(const WorkingCopy &data, const Config &config) {
  if (ThmFiniteDifferencesRequested(config) || (config.paramMask & Config::USE_RMC_FORMALISM)) return false;
  if (data.Data().NumThmCoherentParams() > 0) return false;
  for (int s = 1; s <= data.Data().NumSegments(); s++) {
    const Segment *segment = data.Data().GetSegment(s);
    if (segment && segment->IsTHM() && segment->HasComponents()) return false;
  }
  return true;
}

std::unique_ptr<ThmAdjointContext> ThmGradientContext(CNuc *compound, const WorkingCopy &data, const Config &config,
                                                      const vector_matrix_r *shiftDeriv) {
  bool any = false, lineshape = false;
  for (int s = 1; s <= data.Data().NumSegments(); s++) {
    const Segment *segment = data.Data().GetSegment(s);
    if (!segment || !segment->IsTHM()) continue;
    any = true;
    lineshape = lineshape || segment->GetThmLineshape();
  }
  if (!any || !ThmThroughEnergyTasks(data, config)) return nullptr;
  return std::unique_ptr<ThmAdjointContext>(
      new ThmAdjointContext(BuildThmAdjointContext(compound, config, shiftDeriv, lineshape)));
}

void ThmModelCotangents(WorkingCopy &data, Results &results, const std::vector<GradPoint> &points,
                        const std::vector<double> &models, std::vector<double> &fitBar) {
  // The points of each THM data set (indices into `points`, in data order); their model values in the record, which
  // the experiments' profiles read.
  std::vector<std::vector<std::size_t>> bySegment(data.Data().NumSegments() + 1);
  for (std::size_t k = 0; k < points.size(); k++)
    if (points[k].segment->IsTHM()) {
      bySegment[(std::size_t)points[k].segNum].push_back(k);
      points[k].point->SetModelValue(results, models[k]);
    }
  for (int g = 0; g < data.Data().NumThmGroups(); g++) data.ProfileThmGroup(g, results);
  const double one = 1.0;
  for (int s = 1; s <= data.Data().NumSegments(); s++) {
    const std::vector<std::size_t> &mine = bySegment[(std::size_t)s];
    if (mine.empty()) continue;
    const Segment &segment = *data.Data().GetSegment(s);
    // the model in the datum's frame, f m, and f
    std::vector<double> m, f, d, e;
    for (std::size_t k : mine) {
      AppendInDataFrame(data, segment, points[k].pointIdx, models[k], &one, 1, m, f);
      d.push_back(segment.LabValue(points[k].pointIdx));
      e.push_back(segment.LabError(points[k].pointIdx));
    }
    const int group = data.Data().ThmGroupOf(s);
    double scale = 1.0;
    if (group >= 0) {
      scale = results.objective.thmProfiles[(std::size_t)group].s;
    } else if (segment.IsProfiledNorm()) {
      double norm = 1.0;  // the objective's profile of the models set above (Segment::ProfileNormChiSquared)
      segment.ProfileNormChiSquared(data, results, norm);
      scale = norm != 0.0 ? 1.0 / norm : 0.0;
    } else {
      const double norm = data.State(segment).norm;
      scale = norm != 0.0 ? 1.0 / norm : 0.0;
    }
    for (std::size_t i = 0; i < mine.size(); i++) {
      if (e[i] == 0.0) {
        fitBar[mine[i]] = 0.0;
        continue;
      }
      const std::size_t k = mine[i];
      const double r = group >= 0 ? results.objective.thmProfiles[(std::size_t)group].Residual(
                                        m[i], d[i], e[i], points[k].point->GetCMEnergy())
                                  : (scale * m[i] - d[i]) / e[i];
      fitBar[k] = 2.0 * r * scale * f[i] / e[i];
    }
  }
}

}  // namespace azure2
