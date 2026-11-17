// AZURECalc, the diagnostic dump of a failed evaluation (AZURE_DUMP_NONFINITE, reference/command_line.rst): what
// operator() calls when the switch is set and a chi-squared is non-finite or above the threshold.  A unit of its
// own because nothing of a normal run passes through it; it reads the last gradient of the caller's FitterStates
// and writes nonfinite_dump.txt under a lock of its own.
#include "AZURECalc.h"
#include "RunDiagnostics.h"
#include "azure2/support/Environment.h"
#include "ModelPass.h"
#include "Objective.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>
#include "Config.h"
#include "CompoundSystem.h"
#include "WorkingCopy.h"
#include "Segment.h"
#include "Unit.h"
#include "AZUREParams.h"

/*!
 * AZURE_DUMP_NONFINITE=1 (optional threshold: AZURE_DUMP_NONFINITE=<chi2>, default 1e10): whenever a
 * chi-squared evaluation is non-finite or above the threshold, append to nonfinite_dump.txt (in the
 * project's directory, Config::baseDirectory) the parameters that moved most from the last gradient point, the segments with
 * the largest or non-finite contributions, and their offending points (with any non-finite
 * target-effect sub-points).  At most 20 dumps per run.
 */
void AZURECalc::DumpNonFinite(const vector_r &p, double chi, const CompoundSystem *lc, WorkingCopy *ld, const azure2::Results &results) const {
  double threshold = azure2::support::ProcessEnvironment().Number(azure2::support::Switch::kDumpNonFinite, 0.0);
  if (!(threshold > 1.0)) threshold = 1.0e10;
  if (std::isfinite(chi) && chi <= threshold) return;
  const int d = ++configure().diagnostics->nonFiniteDumps;
  if (d > 20) return;
  AZUREParams names;
  compound()->FillMnParams(names, state_, inputs_, &configure());
  data()->FillMnParams(names);
  const auto &mp = names.GetMinuitParams();
  std::ostringstream out;
  out << "==== dump " << d << ": chi2 " << chi << " (threshold " << threshold << ")\n";
  std::vector<double> base, g2 = states_.gradient.g2, gl;
  {
    std::lock_guard<std::mutex> lk(states_.gradient.mutex);
    base = states_.gradient.at;
    gl = states_.gradient.value;
  }
  out << "-- parameters, sorted by |p - p_lastgradient| / max(|p_lastgradient|, 1e-3)"
      << (base.size() == p.size() ? "" : " (no gradient point yet: |p| only)") << "\n";
  std::vector<std::pair<double, int>> moved;
  for (size_t i = 0; i < p.size(); i++) {
    if ((int)i < (int)mp.Params().size() && mp.Parameter(i).IsFixed()) continue;
    const double b = base.size() == p.size() ? base[i] : 0.0;
    moved.push_back({std::fabs(p[i] - b) / std::max(std::fabs(b), 1.0e-3), (int)i});
  }
  std::sort(moved.rbegin(), moved.rend());
  for (size_t k = 0; k < moved.size() && k < 40; k++) {
    const int i = moved[k].second;
    const double b = base.size() == p.size() ? base[i] : 0.0;
    out << "  idx " << std::setw(4) << i << "  " << std::setw(34) << std::left
        << ((int)i < (int)mp.Params().size() ? mp.GetName(i) : std::string("?")) << std::right
        << "  now " << std::setw(14) << p[i] << "  lastgrad " << std::setw(14) << b
        << "  rel.move " << moved[k].first;
    if (gl.size() == p.size()) out << "  grad " << gl[i];
    if (g2.size() == p.size()) out << "  G2 " << g2[i] << "  -grad/G2 " << (g2[i] != 0.0 && gl.size() == p.size() ? -gl[i] / g2[i] : 0.0);
    out << "\n";
  }
  out << "-- segments with non-finite or largest chi2\n";
  std::vector<std::pair<double, int>> segs;
  for (int i = 1; i <= ld->Data().NumSegments(); i++) {
    const Segment *seg = ld->Data().GetSegment(i);
    if (!seg) continue;
    double c = 0.0;
    for (int pid = 1; pid <= seg->NumPoints(); pid++) {
      Unit *pt = seg->GetPoint(*ld, pid);
      if (pt) c += azure2::PointResidual(*ld, *seg, pid - 1, results, ld->State(*seg).norm).Chi2();
    }
    segs.push_back({std::isfinite(c) ? c : std::numeric_limits<double>::max(), i});
  }
  std::sort(segs.rbegin(), segs.rend());
  for (size_t k = 0; k < segs.size() && k < 5; k++) {
    const Segment *seg = ld->Data().GetSegment(segs[k].second);
    out << "  segment " << segs[k].second << " key " << seg->GetSegmentKey() << " (" << seg->GetEntranceKey() << "->"
        << seg->GetExitKey() << ") chi2 " << (segs[k].first == std::numeric_limits<double>::max() ? std::string("NON-FINITE") : std::to_string(segs[k].first))
        << "  norm " << ld->State(*seg).norm << " shift " << ld->State(*seg).energyShift << " sqrtshift " << ld->State(*seg).energyShiftSqrt << "\n";
    int shown = 0;
    for (int pid = 1; pid <= seg->NumPoints() && shown < 6; pid++) {
      Unit *pt = seg->GetPoint(*ld, pid);
      if (!pt) continue;
      const double fit = pt->GetModelValue(results);
      const double pull = azure2::PointResidual(*ld, *seg, pid - 1, results, ld->State(*seg).norm).Pull();
      if (std::isfinite(fit) && std::fabs(pull) < 1.0e4) continue;
      shown++;
      out << "    point " << pid << " Elab " << pt->GetLabEnergy() << " Ecm " << pt->GetCMEnergy() << " angle "
          << pt->GetCMAngle() << " fit " << fit << " data " << seg->CMValue(ld->State(*seg), pid - 1) << " pull " << pull;
      int nBad = 0;
      double firstBadE = 0.0, firstBadS = 0.0;
      for (int j = 1; j <= pt->NumSubPoints(); j++) {
        const double sv = pt->GetSubPoint(*ld, j)->ValueFor(*pt, j, results);
        if (!std::isfinite(sv) || std::fabs(sv) > 1.0e10) {
          if (!nBad) { firstBadE = pt->GetSubPoint(*ld, j)->GetCMEnergy(); firstBadS = sv; }
          nBad++;
        }
      }
      if (pt->NumSubPoints() > 0)
        out << "  sub-points " << pt->NumSubPoints() << " bad " << nBad
            << (nBad ? " first at Ecm " + std::to_string(firstBadE) + " value " + std::to_string(firstBadS) : std::string());
      out << "\n";
    }
  }
  static std::mutex fileMutex;
  std::lock_guard<std::mutex> lk(fileMutex);
  std::ofstream f(configure().OpenPath("nonfinite_dump.txt"), std::ios::app);
  f << out.str();
  std::cerr << "[nonfinite] chi2 " << chi << ": dump " << d << " written to nonfinite_dump.txt" << std::endl;
}
