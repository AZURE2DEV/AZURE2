#include "AZURECalc.h"
#include "Profile.h"
#include "ScatteringLength.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <string>
#include <map>
#include <set>
#include <fstream>
#include <mutex>
#include <sstream>
#include "ParameterLabel.h"
#include "Config.h"
#include "CNuc.h"
#include "EData.h"
#include "ESegment.h"
#include "EPoint.h"
#include "ParameterLimitsManager.h"
#include "AZUREParams.h"
#include "AZUREGrad.h"
#include "CovarianceBand.h"
#include "GSLException.h"
#include <iostream>
#include <iomanip>
#include <thread>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

#ifdef _OPENMP
#include <omp.h>

// setenv/unsetenv are POSIX and missing from the Windows (MinGW) runtime.
// There, _putenv("NAME=value") sets a variable for getenv and "NAME=" removes it.
namespace {
void setEnvVar(const char *name, const char *value) {
#ifdef _WIN32
  static thread_local std::string entry;
  entry = std::string(name) + "=" + value;
  _putenv(entry.c_str());
#else
  setenv(name, value, 1);
#endif
}
void unsetEnvVar(const char *name) {
#ifdef _WIN32
  static thread_local std::string entry;
  entry = std::string(name) + "=";
  _putenv(entry.c_str());
#else
  unsetenv(name);
#endif
}
}  // namespace
#endif

namespace {

// Value returned to Minuit2 for a non-finite or absurdly large chi-squared.
constexpr double kChi2Ceiling = 1.0e100;

bool InParallelRegion() {
#ifdef _OPENMP
  return omp_in_parallel();
#else
  return false;
#endif
}

// Central-difference step for parameter x0.  The usual scale*(|x0| + 1) is
// absolute for |x0| < 1, which is right for widths and energies that pass
// through zero but not for a normalization: a yield normalized in absolute
// units (11B+a CASPAR thick-target yields, n ~ 2e-5) moved by 1e-6, i.e. 5%,
// and its finite-difference gradient came out 14x too large.  A norm is
// positive and scale-free, so its step is relative.
double FdStep(ParamKind kind, double x0, double scale) {
  if (kind == ParamKind::Norm && x0 != 0.0) return scale * std::fabs(x0);
  return scale * (std::fabs(x0) + 1.0);
}

// AZURE_GRAD_CHECK on the Jacobian path: the Jacobian-path gradient, handed to
// the finite-difference check that the adjoint gradient runs, so both are
// compared with the same finite differences before any exit.
std::vector<double> g_jacGradForCheck;
bool g_deferCheckExit = false;

// chi^2 contribution of a single segment: its data term plus its own
// energy-shift nuisance penalty, evaluated at the segment's current shift.
// Mirrors exactly the per-segment math in AZURECalc::Chi2Value (data term +
// lines guarded by IsVaryEnergyShift), so a finite difference of this over the
// segment's shift equals a finite difference of the whole chi^2 -- *provided*
// no other segment's contribution depends on this shift.  That holds when the
// segment is self-contained: no components, and no point mapping onto segments
// the shift does not move (see the eligibility guard in Gradient()).  All other segments then produce
// identical +h/-h contributions that cancel in the central difference.
double SegmentLocalChi2(ESegment *segment, CNuc *lc, const Config &config,
                        EData *ld) {
  double chi = 0.0;
  const double norm = segment->GetNorm();
  for (int pid = 0; pid < segment->NumPoints(); pid++) {
    double th = segment->CalculateTheoreticalCrossSection(pid, lc, config, ld);
    EPoint *pt = segment->GetPoint(pid + 1);
    if (!pt) continue;
    pt->SetFitCrossSection(th);
    double r = th - pt->GetCMCrossSection() * norm;
    double err = pt->GetCMCrossSectionError() * norm;
    if (err != 0.0) chi += (r * r) / (err * err);
  }
  segment->AddShiftPenalty(chi);
  return chi;
}

// The shift penalties as (full parameter index, nominal, sigma), for the code
// that needs a penalty's derivatives rather than its value: d/dx = 2 (x - nominal)
// / sigma^2, curvature 2 / sigma^2.  Exactly the terms ESegment::AddShiftPenalty
// adds, so the gradient, G2 and the least-squares priors stay consistent with
// chi^2.
template <class F>
void ForEachShiftPenalty(EData *data, const ParamIndexMap &pmap, F f) {
  for (int s = 1; s <= data->NumSegments(); s++) {
    ESegment *seg = data->GetSegment(s);
    if (!seg) continue;
    int idx;
    if (seg->IsVaryEnergyShift() && seg->GetEnergyShiftError() != 0.0 && (idx = pmap.EnergyShiftIndex(s)) >= 0)
      f(idx, seg->GetNominalEnergyShift(), seg->GetEnergyShiftError());
    if (seg->IsVaryEnergyShiftSqrt() && seg->GetEnergyShiftSqrtError() != 0.0 &&
        (idx = pmap.EnergyShiftSqrtIndex(s)) >= 0)
      f(idx, seg->GetNominalEnergyShiftSqrt(), seg->GetEnergyShiftSqrtError());
    if (seg->OwnsShiftFn() && seg->HasShiftFn())
      for (int j = 0; j < seg->NumShiftFnParams(); j++)
        if (seg->IsVaryShiftFn(j) && seg->GetShiftFnErrors()[j] != 0.0 && (idx = pmap.EnergyShiftFnIndex(s, j)) >= 0)
          f(idx, seg->GetNominalShiftFnValues()[j], seg->GetShiftFnErrors()[j]);
  }
}

// The electron-screening priors ((U_e - nominal)/sigma)^2, U_e in eV: f(full index, nominal, sigma) for each
// screening parameter with sigma > 0 (EData::ScreeningPenalty is their sum).
template <class F>
void ForEachScreeningPrior(EData *data, const ParamIndexMap &pmap, F f) {
  for (int k = 1; k <= data->NumScreeningParams(); k++) {
    const TargetEffect *te = data->GetTargetEffect(data->GetScreeningEffect(k));
    if (te->GetScreeningSigma() <= 0.0) continue;
    const int idx = pmap.ScreeningIndex(k);
    if (idx >= 0) f(idx, te->GetScreeningNominal(), te->GetScreeningSigma());
  }
}

// The segments (1-based) whose points one energy-shift parameter moves, and how to
// set it on them: its own segment, plus every follower of its shift group for a
// group's constant shift and for each shift-function parameter.
struct ShiftParamTarget {
  ParamKind kind;
  int fnTerm;
  std::vector<int> segs;
  double Get(EData *ld) const {
    ESegment *seg = ld->GetSegment(segs.front());
    if (kind == ParamKind::EnergyShiftSqrt) return seg->GetEnergyShiftSqrt();
    if (kind == ParamKind::EnergyShiftFn) return seg->GetShiftFnValues()[fnTerm];
    return seg->GetEnergyShift();
  }
  void Set(EData *ld, double v) const {
    for (int s : segs) {
      ESegment *seg = ld->GetSegment(s);
      if (kind == ParamKind::EnergyShiftSqrt) seg->SetEnergyShiftSqrt(v);
      else if (kind == ParamKind::EnergyShiftFn) seg->SetShiftFnValue(fnTerm, v);
      else seg->SetEnergyShift(v);
    }
  }
  // Every member self-contained (no components, not a total-capture sum), so the
  // parameter moves these segments' rows and nothing else.
  bool Local(EData *ld) const {
    for (int s : segs) {
      ESegment *seg = ld->GetSegment(s);
      if (seg->HasComponents() || seg->IsTotalCapture()) return false;
    }
    return true;
  }
};

ShiftParamTarget ShiftTargets(EData *ld, const ParamDesc &d) {
  ShiftParamTarget t{d.kind, d.fnTerm, {d.segment}};
  const int pos = d.segment - 1;
  for (int s = 1; s <= ld->NumSegments(); s++) {
    if (s == d.segment) continue;
    ESegment *seg = ld->GetSegment(s);
    if ((d.kind == ParamKind::EnergyShift && seg->GetShiftLeaderPos() == pos) ||
        (d.kind == ParamKind::EnergyShiftFn && seg->GetShiftFnSourcePos() == pos))
      t.segs.push_back(s);
  }
  return t;
}

// For each segment (1-based), the other segments its points share a calculation
// with through point mapping across segments (EData::MapData): those whose
// points are mapped onto its points, and those its points are mapped onto.
std::vector<std::vector<int>> CrossSegmentPartners(EData *ld) {
  std::vector<std::vector<int>> partners(ld->NumSegments() + 1);
  for (int s = 1; s <= ld->NumSegments(); s++) {
    ESegment *seg = ld->GetSegment(s);
    if (!seg) continue;
    for (int pid = 1; pid <= seg->NumPoints(); pid++) {
      EPoint *pt = seg->GetPoint(pid);
      if (!pt || !pt->IsMapped()) continue;
      const int h = pt->GetMap().segment;
      if (h == s || h < 1 || h > ld->NumSegments()) continue;
      partners[s].push_back(h);
      partners[h].push_back(s);
    }
  }
  for (std::vector<int> &v : partners) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
  }
  return partners;
}

// True when no point of `segs` shares a calculation with a point outside them.
// A shift parameter that moves exactly `segs` can then be differenced on their
// rows alone: EData::MapData maps across segments only where the shifts agree
// as functions of the parameters, so mapped partners move together, and the
// host -- first in data order -- is calculated first (`segs` is ascending).
bool MappingConfined(const std::vector<std::vector<int>> &partners, const std::vector<int> &segs) {
  for (int s : segs)
    for (int q : partners[s])
      if (std::find(segs.begin(), segs.end(), q) == segs.end()) return false;
  return true;
}

}  // namespace

/*!
 * AZURE_DUMP_NONFINITE=1 (optional threshold: AZURE_DUMP_NONFINITE=<chi2>, default 1e10): whenever a
 * chi-squared evaluation is non-finite or above the threshold, append to nonfinite_dump.txt (in the
 * working directory) the parameters that moved most from the last gradient point, the segments with
 * the largest or non-finite contributions, and their offending points (with any non-finite
 * target-effect sub-points).  At most 20 dumps per run.
 */
void AZURECalc::DumpNonFinite(const vector_r &p, double chi, CNuc *lc, EData *ld) const {
  double threshold = std::atof(std::getenv("AZURE_DUMP_NONFINITE"));
  if (!(threshold > 1.0)) threshold = 1.0e10;
  if (std::isfinite(chi) && chi <= threshold) return;
  static std::atomic<int> nDump{0};
  const int d = ++nDump;
  if (d > 20) return;
  AZUREParams names;
  compound()->FillMnParams(names.GetMinuitParams(), &configure());
  data()->FillMnParams(names.GetMinuitParams());
  const auto &mp = names.GetMinuitParams();
  std::ostringstream out;
  out << "==== dump " << d << ": chi2 " << chi << " (threshold " << threshold << ")\n";
  std::vector<double> base, g2 = g2Cache_, gl;
  {
    std::lock_guard<std::mutex> lk(lastGradMutex_);
    base = lastGradP_;
    gl = lastGrad_;
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
  for (int i = 1; i <= ld->NumSegments(); i++) {
    ESegment *seg = ld->GetSegment(i);
    if (!seg) continue;
    double c = 0.0;
    for (int pid = 1; pid <= seg->NumPoints(); pid++) {
      EPoint *pt = seg->GetPoint(pid);
      if (!pt) continue;
      const double r = pt->GetFitCrossSection() - pt->GetCMCrossSection() * seg->GetNorm();
      const double e = pt->GetCMCrossSectionError() * seg->GetNorm();
      if (e != 0.0) c += r * r / (e * e);
    }
    segs.push_back({std::isfinite(c) ? c : std::numeric_limits<double>::max(), i});
  }
  std::sort(segs.rbegin(), segs.rend());
  for (size_t k = 0; k < segs.size() && k < 5; k++) {
    ESegment *seg = ld->GetSegment(segs[k].second);
    out << "  segment " << segs[k].second << " key " << seg->GetSegmentKey() << " (" << seg->GetEntranceKey() << "->"
        << seg->GetExitKey() << ") chi2 " << (segs[k].first == std::numeric_limits<double>::max() ? std::string("NON-FINITE") : std::to_string(segs[k].first))
        << "  norm " << seg->GetNorm() << " shift " << seg->GetEnergyShift() << " sqrtshift " << seg->GetEnergyShiftSqrt() << "\n";
    int shown = 0;
    for (int pid = 1; pid <= seg->NumPoints() && shown < 6; pid++) {
      EPoint *pt = seg->GetPoint(pid);
      if (!pt) continue;
      const double fit = pt->GetFitCrossSection();
      const double e = pt->GetCMCrossSectionError() * seg->GetNorm();
      const double pull = (e != 0.0) ? (fit - pt->GetCMCrossSection() * seg->GetNorm()) / e : 0.0;
      if (std::isfinite(fit) && std::fabs(pull) < 1.0e4) continue;
      shown++;
      out << "    point " << pid << " Elab " << pt->GetLabEnergy() << " Ecm " << pt->GetCMEnergy() << " angle "
          << pt->GetCMAngle() << " fit " << fit << " data " << pt->GetCMCrossSection() << " pull " << pull;
      int nBad = 0;
      double firstBadE = 0.0, firstBadS = 0.0;
      for (int j = 1; j <= pt->NumSubPoints(); j++) {
        const double sv = pt->GetSubPoint(j)->GetFitCrossSection();
        if (!std::isfinite(sv) || std::fabs(sv) > 1.0e10) {
          if (!nBad) { firstBadE = pt->GetSubPoint(j)->GetCMEnergy(); firstBadS = sv; }
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
  std::ofstream f("nonfinite_dump.txt", std::ios::app);
  f << out.str();
  std::cerr << "[nonfinite] chi2 " << chi << ": dump " << d << " written to nonfinite_dump.txt" << std::endl;
}

double AZURECalc::operator()(const vector_r &p) const {
  AZUREProfile::ScopedReport profileReport("evaluation");
  AZUREProfile::Timer profileTimer(AZUREProfile::kEvaluation);
  int thisIteration = data()->Iterations();
  data()->Iterate();
  bool isFit = data()->IsFit();

  CNuc *localCompound = NULL;
  EData *localData = NULL;
  if (isFit) {
    // New multithreading with object pools
    if (!pools_initialized_) {
      InitializePools();
    }

    // Get objects from pool
    localCompound = GetPooledCNuc();
    localData = GetPooledEData();

    // Old multithreading
    // localCompound = compound()->Clone();
    // localData = data()->Clone();
  } else {
    localCompound = compound();
    localData = data();
  }

  // Fill Compound Nucleus From Minuit Parameters
  const auto te0 = std::chrono::steady_clock::now();
  localCompound->FillCompoundFromParams(p);
  localData->FillNormsFromParams(p);
  localData->FillEnergyShiftsFromParams(p, localData, localCompound, &configure());
  // A chi-squared alone: no Park second derivatives (CNuc::CalcShiftFunctions).
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure(), false);
  const auto te1 = std::chrono::steady_clock::now();

  // Calculate all points in one parallel pass when not already inside a parallel
  // region (MIGRAD with the analytic gradient, line searches, plain calculations).
  const bool pre = PrecalculatePoints(localCompound, localData);
  if (std::getenv("AZURE_GRAD_TIMING")) {
    static std::atomic<int> evals{0};
    const int c = ++evals;
    if (c <= 20 || c % 50 == 0) {
      auto ms = [](auto a, auto b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
      std::cerr << "[eval-timing] evaluation " << c << ": parameters + energy shifts " << ms(te0, te1)
                << " ms, points " << ms(te1, std::chrono::steady_clock::now()) << " ms" << std::endl;
    }
  }

  // Process segments with components - use new integrated calculation method
  double chiSquared = 0.0;
  for (int i = 1; i <= localData->NumSegments(); i++) {
    ESegment *segment = localData->GetSegment(i);
    if (segment) {
      // Recalculate points using the new combined calculation method
      for (int pointIdx = 0; pointIdx < segment->NumPoints(); pointIdx++) {
        double theoreticalValue = segment->CalculateTheoreticalCrossSection(pointIdx, localCompound, configure(), localData, pre);
        EPoint *point = segment->GetPoint(pointIdx + 1);
        if (point) {
          point->SetFitCrossSection(theoreticalValue);
        }
      }

      // Recalculate chi-squared for this segment with components
      double segmentChiSquared = 0.0;
      for (int pointIdx = 0; pointIdx < segment->NumPoints(); pointIdx++) {
        EPoint *point = segment->GetPoint(pointIdx + 1);
        if (point) {
          double residual = point->GetFitCrossSection() - point->GetCMCrossSection() * segment->GetNorm();
          double error = point->GetCMCrossSectionError() * segment->GetNorm();
          if (error != 0.0) {
            segmentChiSquared += (residual * residual) / (error * error);
          }
        }
      }

      // Add normalization chi-squared contribution
      double dataNorm = segment->GetNorm();
      double dataNormNominal = segment->GetNominalNorm();
      double dataNormError = dataNormNominal / 100. * segment->GetNormError();
      if (dataNormError != 0.) {
        chiSquared += pow((dataNorm - dataNormNominal) / dataNormError, 2.0);
      }

      // Add the energy-shift penalties
      segment->AddShiftPenalty(chiSquared);

      segment->SetSegmentChiSquared(segmentChiSquared);
      chiSquared += segmentChiSquared;
    }
  }

  // Add nuisance parameter chi-squared contributions
  if (limitsManager_) {
    chiSquared += CalculateNuisanceChiSquared(p);
  }
  // Park formalism: keep every level's overlap J positive.
  if (configure().paramMask & Config::USE_PARK_FORMALISM) {
    double parkPenalty = localCompound->ParkNormPenalty();
    chiSquared += parkPenalty;
    localData->SetParkPenalty(parkPenalty);
  }
  // Electron-screening potentials: their priors.
  chiSquared += localData->ScreeningPenalty();
  // Isotopic abundances: their priors and the sum-to-one term.
  if (localCompound->HasAbundanceBlock()) {
    double abundancePenalty = localCompound->AbundancePenalty(configure().abundanceSumSigma);
    chiSquared += abundancePenalty;
    localData->SetAbundancePenalty(abundancePenalty);
  }
  // Measured scattering lengths (<scatteringLengths>), evaluated at each pair's threshold.
  if (localCompound->HasScatteringLengthConstraints()) {
    double slChi2 = ScatteringLength::Chi2(localCompound, configure(), localCompound->ScatteringLengthConstraints());
    chiSquared += slChi2;
    localData->SetScatteringLengthChi2(slChi2);
  }

  if (std::getenv("AZURE_DUMP_NONFINITE")) DumpNonFinite(p, chiSquared, localCompound, localData);

  if (!localData->IsErrorAnalysis() && thisIteration != 0) {
    if (thisIteration % 10 == 0) configure().outStream
        << "\r\tIteration: " << std::setw(6) << thisIteration
        << " Chi-Squared: " << chiSquared;
    configure().outStream.flush();

    if (thisIteration % kOutputInterval == 0) WriteIterationOutput(p);
  }
  if (isFit) {
    // New multithreading with object pools
    ReturnPooledCNuc(localCompound);
    ReturnPooledEData(localData);

    // Old multithreading
    // delete localCompound;
    // delete localData;
  }

  // A non-finite (or absurdly large) chi-squared must reach Minuit2 as a large
  // FINITE number.  MnLineSearch fits a parabola through the function values;
  // +inf (which this used to return for NaN) turns that parabola into NaN, the
  // line search returns garbage and MIGRAD ends with no status message -- the
  // 11B+a pol_iso --use-gradient fit stopped that way after 30 evaluations
  // (2026-10-11).  A finite ceiling makes the step simply look very bad, so the
  // line search backtracks.  The ceiling stays far above any real chi-squared
  // and far below overflow in the parabola arithmetic.
  if (!std::isfinite(chiSquared) || chiSquared > kChi2Ceiling) chiSquared = kChi2Ceiling;

  if (configure().stopFlag && isFit)
    return 0.;
  else
    return chiSquared;
}

void AZURECalc::WriteIterationOutput(const vector_r &p) const {
  // Work on private copies: this runs in the middle of a fit and must not
  // disturb the objects the minimizer is stepping.  The copies are made once
  // and reused (every snapshot refills them from p below, exactly as the pooled
  // objects are refilled in operator()), instead of being cloned and deleted at
  // every snapshot -- see output_data_ in AZURECalc.h.  The lock also keeps two
  // threads from writing the snapshot files at the same time.
  std::lock_guard<std::mutex> outputLock(output_mutex_);
  if (!output_compound_) output_compound_.reset(compound()->Clone());
  if (!output_data_) output_data_.reset(data()->Clone());
  CNuc *lc = output_compound_.get();
  EData *ld = output_data_.get();

  try {
    lc->FillCompoundFromParams(p);
    ld->FillNormsFromParams(p);
    ld->FillEnergyShiftsFromParams(p, ld, lc, &configure());
    if (configure().paramMask & Config::USE_BRUNE_FORMALISM) lc->CalcShiftFunctions(configure());

    for (int i = 1; i <= ld->NumSegments(); i++) {
      ESegment *segment = ld->GetSegment(i);
      if (!segment) continue;
      for (int pid = 0; pid < segment->NumPoints(); pid++) {
        EPoint *pt = segment->GetPoint(pid + 1);
        if (pt) pt->SetFitCrossSection(
            segment->CalculateTheoreticalCrossSection(pid, lc, configure(), ld));
      }
    }

    AZUREParams params;
    lc->FillMnParams(params.GetMinuitParams(), &configure());
    ld->FillMnParams(params.GetMinuitParams());
    // FillMnParams reads the levels' input values (GetE/GetGamma), not the
    // point being evaluated; copy p in so param.fit matches this snapshot.
    ROOT::Minuit2::MnUserParameters &mp = params.GetMinuitParams();
    if (mp.Params().size() == p.size())
      for (unsigned int i = 0; i < p.size(); i++) mp.SetValue(i, p[i]);
    ld->SyncShiftGroupParams(mp);
    ld->SyncNormGroupParams(mp);
    WriteParameters(params, configure());
    ld->WriteOutputFiles(configure(), true);
    lc->TransformOut(configure());
    lc->PrintTransformParams(configure());
  } catch (...) {
    // An intermediate snapshot is a convenience, never a reason to abort a fit.
    // A failure may leave the copies half-filled: drop them so the next
    // snapshot starts from fresh clones.
    output_compound_.reset();
    output_data_.reset();
  }
}

bool AZURECalc::PrecalculatePoints(CNuc *lc, EData *ld) const {
#ifdef _OPENMP
  if (omp_in_parallel()) return false;
#else
  return false;
#endif
  if (std::getenv("AZURE_SERIAL_CHI2")) return false;
  std::vector<EPoint *> work;
  std::set<ESegment *> seen;
  auto add = [&](ESegment *seg) {
    if (!seg || !seen.insert(seg).second) return;
    for (int p = 1; p <= seg->NumPoints(); p++) {
      EPoint *pt = seg->GetPoint(p);
      if (pt && !pt->IsMapped()) work.push_back(pt);
    }
  };
  for (int i = 1; i <= ld->NumSegments(); i++) {
    ESegment *seg = ld->GetSegment(i);
    if (!seg) continue;
    add(seg);
    for (ESegment *comp : seg->GetComponentSegments()) add(comp);
  }
  // A host also calculates every point mapped onto it (EData::MapData), so
  // its work grows with that number: one Fowler-type angular distribution at
  // a single energy is one item.  Largest items first, one at a time, so such
  // an item does not end up last on one thread.  Each item is calculated the
  // same way whatever the order.
  std::vector<double> cost(work.size());
  for (size_t k = 0; k < work.size(); k++)
    cost[k] = (1.0 + work[k]->NumSubPoints()) * (1.0 + work[k]->NumLocalMappedPoints());
  std::vector<size_t> order(work.size());
  for (size_t k = 0; k < order.size(); k++) order[k] = k;
  std::stable_sort(order.begin(), order.end(), [&cost](size_t a, size_t b) { return cost[a] > cost[b]; });
  std::vector<EPoint *> sorted(work.size());
  for (size_t k = 0; k < order.size(); k++) sorted[k] = work[order[k]];
  work.swap(sorted);
  const Config &cfg = configure();
  AZUREProfile::Timer profileTimer(AZUREProfile::kPoints);
#pragma omp parallel for schedule(dynamic, 1)
  for (int k = 0; k < (int)work.size(); k++) {
    try {
      work[k]->Calculate(lc, cfg);
    } catch (...) {
      work[k]->SetFitCrossSection(0.0);  // as the per-point path, which returns 0 on failure
    }
  }
  return true;
}

double AZURECalc::Chi2On(CNuc *lc, EData *ld, const vector_r &p, std::vector<double> *segChis) const {
  lc->FillCompoundFromParams(p);
  ld->FillNormsFromParams(p);
  ld->FillEnergyShiftsFromParams(p, ld, lc, &configure());
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) lc->CalcShiftFunctions(configure(), false);
  const bool pre = PrecalculatePoints(lc, ld);

  double chiSquared = 0.0;
  if (segChis) segChis->assign(ld->NumSegments() + 1, 0.0);
  for (int i = 1; i <= ld->NumSegments(); i++) {
    ESegment *segment = ld->GetSegment(i);
    if (!segment) continue;
    double segChi = 0.0;
    for (int pid = 0; pid < segment->NumPoints(); pid++) {
      double th = segment->CalculateTheoreticalCrossSection(pid, lc, configure(), ld, pre);
      EPoint *pt = segment->GetPoint(pid + 1);
      if (pt) {
        pt->SetFitCrossSection(th);
        double r = th - pt->GetCMCrossSection() * segment->GetNorm();
        double err = pt->GetCMCrossSectionError() * segment->GetNorm();
        if (err != 0.0) segChi += (r * r) / (err * err);
      }
    }
    if (segChis) (*segChis)[i] = segChi;
    double dataNorm = segment->GetNorm();
    double nom = segment->GetNominalNorm();
    double nerr = nom / 100.0 * segment->GetNormError();
    if (nerr != 0.0) chiSquared += pow((dataNorm - nom) / nerr, 2.0);
    segment->AddShiftPenalty(chiSquared);
    chiSquared += segChi;
  }
  if (limitsManager_) chiSquared += CalculateNuisanceChiSquared(p);
  // Park: the J > 0 wall, as operator() adds it.  Without it the finite-difference gradient columns (and the
  // AZURE_GRAD_CHECK reference) missed the wall's gradient wherever a level is past it.
  if (configure().paramMask & Config::USE_PARK_FORMALISM) chiSquared += lc->ParkNormPenalty();
  // The working copy lc holds the abundances of p (FillCompoundFromParams above); compound() holds whatever the
  // last full evaluation left there, so using it here froze the abundance prior in the LM/GSL cost and in the
  // finite-difference gradient columns (found 2026-10-17 by the isospin-levels session).
  if (lc->HasAbundanceBlock()) chiSquared += lc->AbundancePenalty(configure().abundanceSumSigma);
  if (lc->HasScatteringLengthConstraints())
    chiSquared += ScatteringLength::Chi2(lc, configure(), lc->ScatteringLengthConstraints());
  chiSquared += ld->ScreeningPenalty();
  return chiSquared;
}

double AZURECalc::Chi2Value(const vector_r &p) const {
  // Side-effect-free chi-squared (mirrors the chi-squared math of operator()).
  // Pooled working copies keep their applied energy shifts and mapping between
  // calls; a fresh Clone() re-shifts every shifted segment and rebuilds the
  // mapping each time.
  CNuc *lc = GetPooledCNuc();
  EData *ld = GetPooledEData();
  std::vector<double> segPooled;
  const bool poolCheck = std::getenv("AZURE_POOL_CHECK") != nullptr;
  const double chi = Chi2On(lc, ld, p, poolCheck ? &segPooled : nullptr);
  ReturnPooledCNuc(lc);
  ReturnPooledEData(ld);
  if (poolCheck) {
    static std::atomic<int> pc{0};
    if (++pc <= 5) {
      CNuc *fc = compound()->Clone();
      EData *fd = data()->Clone();
      std::vector<double> segFresh;
      setEnvVar("AZURE_SERIAL_CHI2", "1");  // reference: per-point (serial) path on a fresh copy
      const double chiFresh = Chi2On(fc, fd, p, &segFresh);
      unsetEnvVar("AZURE_SERIAL_CHI2");
      delete fc;
      delete fd;
      std::cerr << "[pool-check] pooled+parallel " << chi << " fresh+serial " << chiFresh << std::endl;
      std::vector<std::pair<double, int>> d;
      for (size_t i = 1; i < segFresh.size() && i < segPooled.size(); i++)
        d.push_back({std::fabs(segPooled[i] - segFresh[i]), (int)i});
      std::sort(d.rbegin(), d.rend());
      for (int k = 0; k < 6 && k < (int)d.size() && d[k].first > 1e-6; k++)
        std::cerr << "[pool-check]   segment " << d[k].second << ": pooled " << segPooled[d[k].second] << " fresh "
                  << segFresh[d[k].second] << std::endl;
    }
  }
  return chi;
}

std::vector<double> AZURECalc::Gradient(const std::vector<double> &p) const {
  AZUREProfile::ScopedReport profileReport("gradient");
  AZUREProfile::Timer profileTimer(AZUREProfile::kGradient);
  {
    std::lock_guard<std::mutex> lk(lastGradMutex_);
    lastGradP_ = p;
  }
  std::vector<double> grad(p.size(), 0.0);
  const auto tg0 = std::chrono::steady_clock::now();
  // With G2 enabled, one residual Jacobian gives both the gradient (2 J^T r +
  // penalties) and Minuit's G2, which it requests at the same point next.
  if (g2Enabled_ && !std::getenv("AZURE_GRAD_ADJOINT") && GradientFromJacobian(p, grad)) {
    if (std::getenv("AZURE_GRAD_TIMING")) {
      static std::atomic<int> jc{0};
      int c = ++jc;
      if (c <= 5 || c % 50 == 0)
        std::cerr << "[grad-timing] call " << c << " (Jacobian path, gradient + G2): "
                  << std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - tg0).count() << " ms" << std::endl;
    }
    // AZURE_GRAD_CHECK: on the first call also take the adjoint gradient (which
    // itself runs the sampled finite-difference check) and compare all components.
    static std::atomic<bool> jChecked{false};
    if (std::getenv("AZURE_GRAD_CHECK") && !jChecked.exchange(true)) {
      setEnvVar("AZURE_GRAD_ADJOINT", "1");
      g_jacGradForCheck = grad;
      g_deferCheckExit = true;
      std::vector<double> ga = Gradient(p);
      g_deferCheckExit = false;
      g_jacGradForCheck.clear();
      unsetEnvVar("AZURE_GRAD_ADJOINT");
      AZUREParams fpc;
      compound()->FillMnParams(fpc.GetMinuitParams(), &configure());
      data()->FillMnParams(fpc.GetMinuitParams());
      const int nMnc = fpc.GetMinuitParams().Params().size();
      double worst = 0.0; int worstIdx = -1, nFree = 0;
      for (size_t i = 0; i < grad.size() && i < ga.size(); i++) {
        if ((int)i < nMnc && fpc.GetMinuitParams().Parameter(i).IsFixed()) continue;  // Minuit ignores these
        nFree++;
        if (grad[i] == 0.0 && ga[i] == 0.0) continue;
        const double rel = std::fabs(grad[i] - ga[i]) / (std::fabs(grad[i]) + std::fabs(ga[i]) + 1.0);
        if (rel > worst) { worst = rel; worstIdx = (int)i; }
      }
      std::cerr << "[grad-check] Jacobian-path vs adjoint gradient over " << nFree << " free parameters: worst relative difference " << worst << " at idx "
                << worstIdx << (worstIdx >= 0 ? " (" + std::to_string(grad[worstIdx]) + " vs " + std::to_string(ga[worstIdx]) + ")" : "")
                << std::endl;
      if (std::string(std::getenv("AZURE_GRAD_CHECK")) == "exit") std::exit(0);
    }
    {
      std::lock_guard<std::mutex> lk(lastGradMutex_);
      lastGrad_ = grad;
    }
    return grad;
  }
  const bool brune = (configure().paramMask & Config::USE_BRUNE_FORMALISM);

  // --- Analytic energy / reduced-width block via the shared adjoint engine. ---
  // Pooled working copies keep their applied energy shifts and mapping between calls; a fresh
  // Clone() re-shifted every shifted segment and rebuilt the mapping each time (~1 s here).
  CNuc *lc = GetPooledCNuc();
  EData *ld = GetPooledEData();
  lc->FillCompoundFromParams(p);
  ld->FillNormsFromParams(p);
  ld->FillEnergyShiftsFromParams(p, ld, lc, &configure());

  vector_matrix_r shiftDeriv;
  const vector_matrix_r *sdp = nullptr;
  if (brune) {
    lc->CalcShiftFunctions(configure());
    shiftDeriv = BuildShiftDerivTable(lc, configure());
    sdp = &shiftDeriv;
  }

  ParamIndexMap pmap = BuildParamIndexMap(lc, ld, std::vector<bool>());
  GradAccum accum;
  accum.Init(lc);

  // chi^2 = sum (fit - data*n)^2 / (err*n)^2 + penalties.  The data-term model
  // cotangent is d(chi2)/d(model) = 2 r / err^2.  The same per-point model also
  // gives the analytic normalization gradient, accumulated per segment here so
  // it needs no extra forward pass:
  //   d(chi2)/d(n_s) = sum_pts [ -2 r data/(cmErr^2 n^2) - 2 r^2/(cmErr^2 n^3) ].
  std::vector<double> normData(ld->NumSegments() + 1, 0.0);
  FitBarFn fitBarFn = [&](ESegment *seg, int i, int pid, double model) -> double {
    EPoint *pt = seg->GetPoint(pid + 1);
    if (!pt) return 0.0;
    double norm = seg->GetNorm();
    double dataval = pt->GetCMCrossSection();
    double cmErr = pt->GetCMCrossSectionError();
    double r = model - dataval * norm;
    double err = cmErr * norm;
    if (err == 0.0) return 0.0;
    // A norm-group follower's points move its leader's norm (o), the one parameter.
    const int o = ld->NormOwner(i);
    if (norm != 0.0 && o >= 1 && o < (int)normData.size() && ld->GetSegment(o)->IsVaryNorm()) {
      double e2 = cmErr * cmErr;
      // fitBarFn runs inside the parallel point loop of AccumulateEGammaGradient,
      // and all points of a segment share the same normData[o], so guard the
      // accumulation.
      double dNorm = -2.0 * r * dataval / (e2 * norm * norm) - 2.0 * r * r / (e2 * norm * norm * norm);
#pragma omp atomic
      normData[o] += dNorm;
    }
    return 2.0 * r / (err * err);
  };

  bool eg = AccumulateEGammaGradient(lc, ld, configure(), pmap, sdp, fitBarFn, accum);
  const auto tg1 = std::chrono::steady_clock::now();
  if (!eg && std::getenv("AZURE_GRAD_DEBUG")) {
    std::cerr << "[grad] analytic energy/gamma adjoint bailed -> full finite "
                 "differences (no speed-up). An unsupported segment/config is "
                 "present (RMC, or a not-yet-handled case)."
              << std::endl;
  }
  // Measured scattering lengths: exact through the level-matrix adjoint at each threshold; where that
  // is unsupported (RMC, levels left out of the R-matrix) every parameter falls back to finite differences.
  if (eg && lc->HasScatteringLengthConstraints() &&
      !ScatteringLength::Chi2Gradient(lc, configure(), lc->ScatteringLengthConstraints(), accum, sdp))
    eg = false;
  if (eg) {
    AddParkPenaltyGradient(lc, configure(), accum); // Park J > 0 wall
    AddAbundancePenaltyGradient(lc, configure(), accum);  // abundance priors, sum-to-one
    accum.Scatter(pmap, grad);                      // energies + reduced widths, abundances, screening
    AddNuisanceGradient(p, grad);                   // nuisance penalty
    ForEachScreeningPrior(ld, pmap, [&](int idx, double nominal, double e) {  // screening priors
      if (idx < (int)grad.size()) grad[idx] += 2.0 * (p[idx] - nominal) / (e * e);
    });
    for (int s = 1; s <= ld->NumSegments(); s++) {  // normalizations
      ESegment *seg = ld->GetSegment(s);
      if (!seg || !seg->IsVaryNorm() || seg->IsNormGroupFollower()) continue;  // a follower's norm is fixed
      int idx = pmap.NormIndex(s);
      if (idx < 0 || idx >= (int)grad.size()) continue;
      double g = normData[s];
      double n0 = seg->GetNominalNorm();
      double nerr = n0 / 100.0 * seg->GetNormError();
      if (nerr != 0.0) g += 2.0 * (seg->GetNorm() - n0) / (nerr * nerr);
      grad[idx] = g;
    }
  }

  // Fixed-parameter mask: Minuit ignores the gradient of fixed parameters, so
  // do not waste finite differences on them (large fits fix most energy shifts).
  AZUREParams fp;
  compound()->FillMnParams(fp.GetMinuitParams(), &configure());
  data()->FillMnParams(fp.GetMinuitParams());
  const int nMn = fp.GetMinuitParams().Params().size();

  // Segment-local energy-shift gradient is exact only when a shift's effect is
  // confined to its own segments' data term.  Point mapping across segments
  // (EData::MapData) keeps it so as long as the mapped partners of those
  // segments are among them (MappingConfined): partners move together, and
  // groupChi re-energises them all before it calculates them, host first.
  // Nuisance penalties handled by the limits manager can couple a shift to the
  // global chi^2, so disable the local path when one is active.
  const std::vector<std::vector<int>> partners = CrossSegmentPartners(ld);
  const bool localShiftOk = !HasNuisanceParameters();

  // --- Finite differences only for what is left: non-fixed energy shifts, and
  //     (if the analytic path bailed) the energy/gamma and norm blocks too. ---
  for (int idx = 0; idx < (int)p.size() && idx < pmap.NumFull(); idx++) {
    if (idx < nMn && fp.GetMinuitParams().Parameter(idx).IsFixed()) continue;
    ParamKind kind = pmap.Desc(idx).kind;
    if (eg && (kind == ParamKind::LevelEnergy || kind == ParamKind::Gamma || kind == ParamKind::Norm ||
               kind == ParamKind::Screening)) continue;  // analytic
    double x0 = p[idx];
    double h = FdStep(kind, x0, 1.0e-6);

    // Fast path: a free energy shift on self-contained segments only moves those
    // segments' points (its own, plus its shift group's followers for a group's
    // constant shift or a shift-function parameter), so finite-difference just
    // their chi^2 on the already-filled clone instead of re-solving the entire
    // dataset twice.
    if ((kind == ParamKind::EnergyShift || kind == ParamKind::EnergyShiftSqrt || kind == ParamKind::EnergyShiftFn) &&
        localShiftOk && pmap.Desc(idx).segment >= 1) {
      const int s = pmap.Desc(idx).segment;
      const ShiftParamTarget target = ShiftTargets(ld, pmap.Desc(idx));
      if (target.Local(ld) && MappingConfined(partners, target.segs)) {
        auto groupChi = [&](double v) {
          target.Set(ld, v);
          // Every segment is moved before any is calculated: a host calculates
          // the observables of points mapped onto it from other segments, with
          // their energy-dependent values (EPoint::Calculate).
          for (int m : target.segs) ld->GetSegment(m)->UpdatePointEnergiesWithShift(lc, &configure());
          double c = 0.0;
          for (int m : target.segs) c += SegmentLocalChi2(ld->GetSegment(m), lc, configure(), ld);
          return c;
        };
        const double d0 = target.Get(ld);
        double chiP = groupChi(x0 + h);
        double chiM = groupChi(x0 - h);
        target.Set(ld, d0);  // restore base state on the clone
        for (int m : target.segs) ld->GetSegment(m)->UpdatePointEnergiesWithShift(lc, &configure());
        grad[idx] = (chiP - chiM) / (2.0 * h);

        // Optional self-check: compare against the exact full-dataset finite
        // difference for this shift.  Enable with AZURE_SHIFT_GRAD_CHECK=1 to
        // validate the fast path on real data (one gradient evaluation), then
        // disable it for the actual fit.
        if (std::getenv("AZURE_SHIFT_GRAD_CHECK")) {
          vector_r pp = p;
          pp[idx] = x0 + h;
          vector_r pm = p;
          pm[idx] = x0 - h;
          double gFull = (Chi2Value(pp) - Chi2Value(pm)) / (2.0 * h);
          double denom = std::max(1.0, std::fabs(gFull));
          std::cerr << "[shift-grad-check] seg " << s << " param " << idx
                    << ": local=" << grad[idx] << " full=" << gFull
                    << " reldiff=" << std::fabs(grad[idx] - gFull) / denom
                    << std::endl;
        }
        continue;
      }
    }

    vector_r pp = p;
    pp[idx] = x0 + h;
    vector_r pm = p;
    pm[idx] = x0 - h;
    grad[idx] = (Chi2Value(pp) - Chi2Value(pm)) / (2.0 * h);
  }

  const auto tg2 = std::chrono::steady_clock::now();
  static std::atomic<int> gradCalls{0};
  if (std::getenv("AZURE_GRAD_CHECK") || std::getenv("AZURE_GRAD_TIMING")) {
    int c = ++gradCalls;
    if (c <= 5 || c % 50 == 0) {
      auto ms = [](auto a, auto b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
      int nfd = 0;
      for (int idx = 0; idx < (int)p.size() && idx < pmap.NumFull(); idx++) {
        if (idx < nMn && fp.GetMinuitParams().Parameter(idx).IsFixed()) continue;
        ParamKind kind = pmap.Desc(idx).kind;
        if (!(eg && (kind == ParamKind::LevelEnergy || kind == ParamKind::Gamma || kind == ParamKind::Norm ||
                     kind == ParamKind::Screening))) nfd++;
      }
      std::cerr << "[grad-timing] call " << c << ": setup+adjoint " << ms(tg0, tg1) << " ms (analytic " << (eg ? "ok" : "BAILED")
                << "), finite-difference part " << ms(tg1, tg2) << " ms for " << nfd << " parameters" << std::endl;
    }
  }

  // Optional full self-check, AZURE_GRAD_CHECK=1 (or =exit to stop the run after
  // it): on the first call, compare every non-fixed gradient component with a
  // central finite difference of the side-effect-free chi-squared and list the
  // worst disagreements.  Minuit itself never checks (CheckGradient() is false).
  static std::atomic<bool> gradChecked{false};
  if (const char *gc = std::getenv("AZURE_GRAD_CHECK")) {
    if (!gradChecked.exchange(true)) {
      static const char *kindName[] = {"energy", "gamma", "norm", "shift", "shift_sqrt", "shift_fn", "abundance", "screening"};
      struct Row { double score, g, fd; int idx; ParamDesc d; double gj, scoreJ; };
      std::vector<Row> rows;
      constexpr int kNumKinds = sizeof(kindName) / sizeof(kindName[0]);  // one per ParamKind
      int nBad[kNumKinds] = {}, nAll[kNumKinds] = {}, nBadJ[kNumKinds] = {};
      const bool haveJ = g_jacGradForCheck.size() == grad.size();
      const auto tc0 = std::chrono::steady_clock::now();
      const double f0 = Chi2Value(p);
      std::cerr << "[grad-check] one Chi2Value: " << std::chrono::duration<double>(std::chrono::steady_clock::now() - tc0).count() << " s" << std::endl;
      std::vector<int> freeIdx;
      for (int idx = 0; idx < (int)p.size() && idx < pmap.NumFull(); idx++)
        if (!(idx < nMn && fp.GetMinuitParams().Parameter(idx).IsFixed())) freeIdx.push_back(idx);
      const char *cn = std::getenv("AZURE_GRAD_CHECK_N");
      const int nCheck = cn ? std::atoi(cn) : 16;
      const int stride = std::max(1, (int)freeIdx.size() / std::max(1, nCheck));
      std::vector<int> checkIdx;
      for (int q = 0; q < (int)freeIdx.size(); q += stride) checkIdx.push_back(freeIdx[q]);
      // Targeted additions: AZURE_GRAD_CHECK_IDX = "i,j,a-b" (parameter
      // indices) and AZURE_GRAD_CHECK_SEGS = "s,a-b" (1-based segment numbers
      // in data order: their norm and energy-shift parameters).  Fixed
      // parameters are skipped.
      auto parseList = [](const char *spec) {
        std::vector<int> out;
        std::stringstream ss(spec);
        std::string item;
        while (std::getline(ss, item, ',')) {
          if (item.empty()) continue;
          size_t dash = item.find('-', 1);
          int lo = std::atoi(item.substr(0, dash).c_str());
          int hi = (dash == std::string::npos) ? lo : std::atoi(item.substr(dash + 1).c_str());
          for (int v = lo; v <= hi; v++) out.push_back(v);
        }
        return out;
      };
      std::vector<int> extra;
      if (const char *ci = std::getenv("AZURE_GRAD_CHECK_IDX"))
        for (int v : parseList(ci)) extra.push_back(v);
      if (const char *cs = std::getenv("AZURE_GRAD_CHECK_SEGS"))
        for (int sNum : parseList(cs)) {
          if (sNum < 1 || sNum > ld->NumSegments()) continue;
          extra.push_back(pmap.NormIndex(sNum));
          extra.push_back(pmap.EnergyShiftIndex(sNum));
          extra.push_back(pmap.EnergyShiftSqrtIndex(sNum));
          for (int j = 0; pmap.EnergyShiftFnIndex(sNum, j) >= 0; j++) extra.push_back(pmap.EnergyShiftFnIndex(sNum, j));
        }
      for (int v : extra)
        if (std::find(freeIdx.begin(), freeIdx.end(), v) != freeIdx.end() &&
            std::find(checkIdx.begin(), checkIdx.end(), v) == checkIdx.end())
          checkIdx.push_back(v);
      for (int idx : checkIdx) {
        // AZURE_GRAD_CHECK_H scales the finite-difference step (default 1e-6
        // relative): a component that moves with it is truncation error of the
        // FD (e.g. near a narrow level), not of the analytic value.
        static const double hScale = std::getenv("AZURE_GRAD_CHECK_H") ? std::atof(std::getenv("AZURE_GRAD_CHECK_H")) : 1.0e-6;
        double x0 = p[idx], h = FdStep(pmap.Desc(idx).kind, x0, hScale);
        vector_r pp = p, pm = p;
        pp[idx] = x0 + h;
        pm[idx] = x0 - h;
        double fd = (Chi2Value(pp) - Chi2Value(pm)) / (2.0 * h);
        auto scoreOf = [&](double g) {
          return std::fabs(g - fd) / (std::fabs(fd) + std::fabs(g) + 1.0e-6 * std::fabs(f0) + 1.0);
        };
        double score = scoreOf(grad[idx]);
        double gj = haveJ ? g_jacGradForCheck[idx] : 0.0, scoreJ = haveJ ? scoreOf(gj) : 0.0;
        ParamDesc d = pmap.Desc(idx);
        int k = (int)d.kind;
        nAll[k]++;
        if (score > 1.0e-2) nBad[k]++;
        if (scoreJ > 1.0e-2) nBadJ[k]++;
        rows.push_back({score, grad[idx], fd, idx, d, gj, scoreJ});
      }
      std::sort(rows.begin(), rows.end(),
                [](const Row &a, const Row &b) { return std::max(a.score, a.scoreJ) > std::max(b.score, b.scoreJ); });
      std::cerr << "[grad-check] chi2 " << f0 << ", analytic block " << (eg ? "used" : "BAILED (all finite differences)")
                << "; components with relative disagreement > 1e-2:";
      for (int k = 0; k < kNumKinds; k++)
        if (nAll[k]) std::cerr << " " << kindName[k] << " " << nBad[k] << "/" << nAll[k];
      if (haveJ) {
        std::cerr << "; Jacobian path:";
        for (int k = 0; k < kNumKinds; k++)
          if (nAll[k]) std::cerr << " " << kindName[k] << " " << nBadJ[k] << "/" << nAll[k];
      }
      std::cerr << std::endl;
      // Full precision, so two runs (e.g. AZURE_NO_GROUPED_GRADIENT) can be
      // compared value by value; AZURE_GRAD_CHECK_ROWS lists more than 60.
      const std::streamsize prec = std::cerr.precision(15);
      const int maxRows = std::getenv("AZURE_GRAD_CHECK_ROWS") ? std::atoi(std::getenv("AZURE_GRAD_CHECK_ROWS")) : 60;
      for (int r = 0; r < (int)rows.size() && r < maxRows; r++) {
        const Row &w = rows[r];
        std::cerr << "[grad-check] idx " << w.idx << " " << kindName[(int)w.d.kind] << " J" << w.d.jGroup << " L" << w.d.level
                  << " ch" << w.d.channel << " seg" << w.d.segment << "  analytic " << w.g << "  fd " << w.fd
                  << "  rel " << w.score;
        if (haveJ) std::cerr << "  jacobian " << w.gj << "  relJ " << w.scoreJ;
        std::cerr << std::endl;
      }
      std::cerr.precision(prec);
      if (std::string(gc) == "exit" && !g_deferCheckExit) std::exit(0);
    }
  }

  ReturnPooledCNuc(lc);
  ReturnPooledEData(ld);

  return grad;
}

bool AZURECalc::GradientFromJacobian(const std::vector<double> &p, std::vector<double> &grad) const {
  vector_r r, jac;
  std::vector<int> packedToFull;
  if (!ResidualJacobian(p, r, jac, packedToFull) || r.empty()) return false;
  const size_t nCols = packedToFull.size(), nRes = r.size();
  std::vector<double> g2(p.size(), 0.0);
  grad.assign(p.size(), 0.0);
  for (size_t i = 0; i < nRes; i++)
    for (size_t a = 0; a < nCols; a++) {
      const double d = jac[i * nCols + a];
      const int f = packedToFull[a];
      if (f < 0 || f >= (int)p.size()) continue;
      grad[f] += 2.0 * r[i] * d;
      g2[f] += 2.0 * d * d;
    }
  // Penalties (x - x0)^2 / sigma^2: gradient 2 (x - x0) / sigma^2, curvature 2 / sigma^2.
  ParamIndexMap pmap = BuildParamIndexMap(compound(), data(), std::vector<bool>());
  for (int s = 1; s <= data()->NumSegments(); s++) {
    ESegment *seg = data()->GetSegment(s);
    if (!seg) continue;
    int idx;
    if (seg->IsVaryNorm() && (idx = pmap.NormIndex(s)) >= 0 && idx < (int)p.size()) {
      const double n0 = seg->GetNominalNorm(), nerr = n0 / 100.0 * seg->GetNormError();
      if (nerr != 0.0) { grad[idx] += 2.0 * (p[idx] - n0) / (nerr * nerr); g2[idx] += 2.0 / (nerr * nerr); }
    }
  }
  // Park: the J > 0 wall, sum over the levels with J < 0 of (J/s)^2, is not a
  // residual.  Its gradient (zero while every J > 0) as the adjoint path adds it.
  if (configure().paramMask & Config::USE_PARK_FORMALISM) {
    CNuc *lc = GetPooledCNuc();
    lc->FillCompoundFromParams(p);
    lc->CalcShiftFunctions(configure());
    if (lc->ParkNormPenalty() > 0.0) {
      GradAccum accum;
      accum.Init(lc);
      AddParkPenaltyGradient(lc, configure(), accum);
      accum.Scatter(pmap, grad);
    }
    ReturnPooledCNuc(lc);
  }
  ForEachShiftPenalty(data(), pmap, [&](int idx, double nominal, double e) {
    if (idx >= (int)p.size()) return;
    grad[idx] += 2.0 * (p[idx] - nominal) / (e * e);
    g2[idx] += 2.0 / (e * e);
  });
  ForEachScreeningPrior(data(), pmap, [&](int idx, double nominal, double e) {
    if (idx >= (int)p.size()) return;
    grad[idx] += 2.0 * (p[idx] - nominal) / (e * e);
    g2[idx] += 2.0 / (e * e);
  });
  // Isotopic abundances: prior and sum-to-one gradient and (diagonal) curvature.  The sample sums are taken at
  // p: compound() holds the abundances of the last full evaluation, so its AbundanceSum was stale here.
  std::vector<double> sumAtP;
  if (compound()->HasAbundanceBlock()) {
    sumAtP.assign(compound()->NumSamples(), 0.0);
    for (int smp = 1; smp <= compound()->NumSamples(); smp++) sumAtP[smp - 1] = compound()->AbundanceSum(smp);
    for (int k = 1; k <= compound()->NumAbundanceParams(); k++) {
      const int idx = pmap.AbundanceIndex(k);
      if (idx < 0 || idx >= (int)p.size()) continue;
      int sample, nucleus;
      compound()->AbundanceOf(k, sample, nucleus);
      sumAtP[sample - 1] += p[idx] - compound()->GetSample(sample)->entries[nucleus - 1].abundance;
    }
  }
  if (compound()->HasAbundanceBlock()) {
    for (int k = 1; k <= compound()->NumAbundanceParams(); k++) {
      const int idx = pmap.AbundanceIndex(k);
      if (idx < 0 || idx >= (int)p.size()) continue;
      int sample, nucleus;
      compound()->AbundanceOf(k, sample, nucleus);
      const SampleAbundance &e = compound()->GetSample(sample)->entries[nucleus - 1];
      if (e.sigma > 0.0) {
        grad[idx] += 2.0 * (p[idx] - e.nominal) / (e.sigma * e.sigma);
        g2[idx] += 2.0 / (e.sigma * e.sigma);
      }
      const double s = compound()->SampleSumSigma(sample, configure().abundanceSumSigma);
      if (s > 0.0) {
        grad[idx] += 2.0 * (sumAtP[sample - 1] - 1.0) / (s * s);
        g2[idx] += 2.0 / (s * s);
      }
    }
  }
  g2CacheP_ = p;
  g2Cache_ = g2;
  return true;
}

std::vector<double> AZURECalc::G2(const std::vector<double> &p) const {
  if (!g2Cache_.empty() && g2CacheP_ == p) return g2Cache_;
  const auto t0 = std::chrono::steady_clock::now();
  std::vector<double> g2(p.size(), 0.0);
  vector_r r, jac;
  std::vector<int> packedToFull;
  if (!ResidualJacobian(p, r, jac, packedToFull) || r.empty()) return std::vector<double>();
  const size_t nCols = packedToFull.size(), nRes = r.size();
  for (size_t i = 0; i < nRes; i++)
    for (size_t a = 0; a < nCols; a++) {
      const double d = jac[i * nCols + a];
      const int f = packedToFull[a];
      if (f >= 0 && f < (int)g2.size()) g2[f] += 2.0 * d * d;
    }
  // Penalty curvatures, (x - x0)^2 / sigma^2 -> 2 / sigma^2.
  ParamIndexMap pmap = BuildParamIndexMap(compound(), data(), std::vector<bool>());
  for (int s = 1; s <= data()->NumSegments(); s++) {
    ESegment *seg = data()->GetSegment(s);
    if (!seg) continue;
    int idx;
    if (seg->IsVaryNorm() && (idx = pmap.NormIndex(s)) >= 0 && idx < (int)g2.size()) {
      const double nerr = seg->GetNominalNorm() / 100.0 * seg->GetNormError();
      if (nerr != 0.0) g2[idx] += 2.0 / (nerr * nerr);
    }
  }
  ForEachShiftPenalty(data(), pmap, [&](int idx, double, double e) {
    if (idx < (int)g2.size()) g2[idx] += 2.0 / (e * e);
  });
  ForEachScreeningPrior(data(), pmap, [&](int idx, double, double e) {
    if (idx < (int)g2.size()) g2[idx] += 2.0 / (e * e);
  });
  if (compound()->HasAbundanceBlock()) {
    for (int k = 1; k <= compound()->NumAbundanceParams(); k++) {
      const int idx = pmap.AbundanceIndex(k);
      if (idx < 0 || idx >= (int)g2.size()) continue;
      int sample, nucleus;
      compound()->AbundanceOf(k, sample, nucleus);
      const SampleAbundance &e = compound()->GetSample(sample)->entries[nucleus - 1];
      if (e.sigma > 0.0) g2[idx] += 2.0 / (e.sigma * e.sigma);
      const double s = compound()->SampleSumSigma(sample, configure().abundanceSumSigma);
      if (s > 0.0) g2[idx] += 2.0 / (s * s);
    }
  }
  static std::atomic<int> g2Calls{0};
  if (std::getenv("AZURE_GRAD_TIMING")) {
    int c = ++g2Calls;
    if (c <= 5 || c % 50 == 0)
      std::cerr << "[g2-timing] call " << c << ": " << std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count()
                << " ms" << std::endl;
  }
  return g2;
}

bool AZURECalc::HasNuisanceParameters() const {
  // The CLI always hands AZURECalc a ParameterLimitsManager, so "a manager
  // exists" is not the question; whether it declares any nuisance parameter is.
  if (!limitsManager_) return false;
  AZUREParams fp;
  compound()->FillMnParams(fp.GetMinuitParams(), &configure());
  data()->FillMnParams(fp.GetMinuitParams());
  // Index the parameters exactly as CalculateNuisanceChiSquared and
  // AddNuisanceGradient do, and ignore the ones they skip: a nuisance flag on a
  // norm or shift adds nothing to chi^2 or the gradient (those penalties come
  // from the segment's own norm/shift error, which G2 already includes), so it
  // must not switch off the analytic G2 or the segment-local shift derivative.
  // The GUI writes nuisance = 1 on penalized norms, so without this filter a
  // model with norm uncertainties lost both (11B+a, 2026-10-10).
  int nf = 0;
  for (const auto &mp : fp.GetMinuitParams().Parameters()) {
    const std::string &name = mp.GetName();
    if (mp.IsFixed() && name.find("segment") == std::string::npos) continue;
    const bool inert = name.find("norm") != std::string::npos || name.find("shift") != std::string::npos ||
                       name.find("abundance") != std::string::npos || name.find("screening") != std::string::npos;
    if (!inert && limitsManager_->IsNuisanceParameterByIndex(nf)) {
      if (std::getenv("AZURE_GRAD_DEBUG"))
        std::cerr << "[grad] analytic G2 off: nuisance prior on '" << name << "'" << std::endl;
      return true;
    }
    nf++;
  }
  return false;
}

bool AZURECalc::EnableG2(const std::vector<double> &p) {
  g2Enabled_ = false;
  if (HasNuisanceParameters()) return false;  // nuisance penalties are not in the Jacobian
  std::vector<double> g2 = G2(p);
  g2Enabled_ = !g2.empty();
  return g2Enabled_;
}

bool AZURECalc::ResidualJacobian(const vector_r &full, vector_r &residuals,
                                 vector_r &jac, std::vector<int> &packedToFull) const {
  // Pooled working copies keep their applied energy shifts and mapping between calls; a fresh
  // Clone() re-shifted every shifted segment and rebuilt the mapping each time (~1 s here).
  CNuc *lc = GetPooledCNuc();
  EData *ld = GetPooledEData();
  lc->FillCompoundFromParams(full);
  ld->FillNormsFromParams(full);
  ld->FillEnergyShiftsFromParams(full, ld, lc, &configure());

  vector_matrix_r shiftDeriv;
  const vector_matrix_r *sdp = nullptr;
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) {
    lc->CalcShiftFunctions(configure());
    shiftDeriv = BuildShiftDerivTable(lc, configure());
    sdp = &shiftDeriv;
  }

  // Fixed-parameter mask, in the Minuit parameter order.
  AZUREParams tp;
  compound()->FillMnParams(tp.GetMinuitParams(), &configure());
  data()->FillMnParams(tp.GetMinuitParams());
  int nMn = tp.GetMinuitParams().Params().size();
  std::vector<bool> fixed(nMn);
  for (int i = 0; i < nMn; i++) fixed[i] = tp.GetMinuitParams().Parameter(i).IsFixed();

  ParamIndexMap pmap = BuildParamIndexMap(lc, ld, fixed);
  int nCols = 0;
  const auto tj0 = std::chrono::steady_clock::now();
  bool ok = ComputeResidualJacobian(lc, ld, configure(), pmap, sdp, residuals, jac, nCols);
  const auto tj1 = std::chrono::steady_clock::now();
  int nShiftLocal = 0, nShiftFull = 0;

  // Energy-shift columns by central differences, as AZUREAPI::
  // CalculateResidualJacobianRWA does for pyazr.  The adjoint does not cover a
  // shift (it moves the energy of every quantity of the segment), and left as
  // zero the Levenberg-Marquardt and GSL drivers flag the shifts as
  // data-insensitive and never move them.  Differencing the residual vector
  // of the same function keeps the rows aligned.
  if (ok) {
    const size_t nRes = residuals.size();
    vector_r rPlus, rMinus, jTmp;
    int nc = 0;

    // A shift of a self-contained segment (no components, not a total-capture
    // sum, no point mapping onto other segments) moves only that segment's rows, so
    // its column is differenced on those rows alone; every other row of the
    // full re-evaluation below would difference to exactly zero.  The full
    // path re-solved the WHOLE Jacobian, adjoints included, twice per shift:
    // with 55 free shifts on 11B+a one gradient took over an hour.  Same
    // criterion as the segment-local shift gradient in Gradient().
    bool rowLayoutUnknown = false;
    const std::vector<std::vector<int>> partners = CrossSegmentPartners(ld);
    std::vector<int> rowStart(ld->NumSegments() + 2, 0);
    {
      int row = 0;
      for (int s = 1; s <= ld->NumSegments(); s++) {
        rowStart[s] = row;
        ESegment *seg = ld->GetSegment(s);
        if (!seg) continue;
        for (int pid = 0; pid < seg->NumPoints(); pid++) {
          EPoint *pt = seg->GetPoint(pid + 1);
          if (!pt) continue;
          row++;
        }
      }
      rowStart[ld->NumSegments() + 1] = row;
      if ((size_t)row != nRes) rowLayoutUnknown = true;  // row layout not understood: use the full path
    }
    // Residuals of one segment's rows, in row order, exactly as
    // ComputeResidualJacobian forms them.
    auto segmentResiduals = [&](ESegment *seg, vector_r &out, bool parallel) {
      const int nPts = seg->NumPoints();
      const double norm = seg->GetNorm();
      std::vector<double> byPoint(nPts, 0.0);
      std::vector<char> present(nPts, 0);
      auto one = [&](int pid) {
        EPoint *pt = seg->GetPoint(pid + 1);
        if (!pt) return;
        present[pid] = 1;
        const double model = seg->CalculateTheoreticalCrossSection(pid, lc, configure(), ld);
        const double denom = pt->GetCMCrossSectionError() * norm;
        byPoint[pid] = (denom != 0.0) ? (model - pt->GetCMCrossSection() * norm) / denom : 0.0;
      };
      // No parallel construct inside an enclosing parallel region (see
      // CoulFuncCache::localCache).  In parallel each point is calculated on
      // its own, without point mapping (ScopedPointMappingOff, as in the
      // AZUREGrad.cpp loops): a mapped point would otherwise read a value its
      // host may not have calculated yet.  Serially, index order puts every
      // host first and the sharing is kept.
      if (parallel && !InParallelRegion()) {
#pragma omp parallel for schedule(dynamic)
        for (int pid = 0; pid < nPts; pid++) {
          ScopedPointMappingOff noPointMapping;
          one(pid);
        }
      } else {
        for (int pid = 0; pid < nPts; pid++) one(pid);
      }
      out.clear();
      for (int pid = 0; pid < nPts; pid++)
        if (present[pid]) out.push_back(byPoint[pid]);
    };

    // Local shift columns: a parameter whose segments (its own, plus its shift
    // group's followers for a group's constant shift or a shift-function
    // parameter) are all self-contained is differenced on their rows alone.
    struct LocalShift { int f, packed; ShiftParamTarget target; };
    std::vector<LocalShift> localCols;
    std::vector<int> fullShifts;
    for (int f = 0; f < pmap.NumFull() && f < (int)full.size(); f++) {
      if (fixed[f]) continue;
      const ParamKind kind = pmap.Desc(f).kind;
      if (kind != ParamKind::EnergyShift && kind != ParamKind::EnergyShiftSqrt && kind != ParamKind::EnergyShiftFn) continue;
      const int packed = pmap.FullToPacked(f);
      if (packed < 0 || packed >= nCols) continue;
      const int sNum = pmap.Desc(f).segment;
      if (!rowLayoutUnknown && sNum >= 1 && sNum <= ld->NumSegments()) {
        const ShiftParamTarget target = ShiftTargets(ld, pmap.Desc(f));
        // Mapping across segments keeps the column local while it stays among
        // the parameter's segments (MappingConfined): they are all moved, then
        // all evaluated, below.
        if (target.Local(ld) && MappingConfined(partners, target.segs)) {
          localCols.push_back({f, packed, target});
          continue;
        }
      }
      fullShifts.push_back(f);
    }
    // Segments in sequence, each one's points in parallel: re-energising a
    // point (UpdatePointEnergyWithShift) and evaluating it touch only that
    // point.  (Parallel over segments instead left the largest segment running
    // alone at the end.)
    auto updateSegment = [&](ESegment *seg) {
      if (!InParallelRegion()) {
#pragma omp parallel for schedule(dynamic, 1)
        for (int pid = 0; pid < seg->NumPoints(); pid++) seg->UpdatePointGroupEnergyWithShift(pid, lc, &configure());
      } else {
        seg->UpdatePointEnergiesWithShift(lc, &configure());
      }
    };
    for (const LocalShift &ls : localCols) {
      const std::vector<int> &segs = ls.target.segs;
      std::vector<vector_r> rP(segs.size()), rM(segs.size());
      const double x0 = full[ls.f];
      const double h = 1.0e-6 * (std::fabs(x0) + 1.0);
      const double d0 = ls.target.Get(ld);
      ls.target.Set(ld, x0 + h);
      for (size_t m = 0; m < segs.size(); m++) updateSegment(ld->GetSegment(segs[m]));
      for (size_t m = 0; m < segs.size(); m++) segmentResiduals(ld->GetSegment(segs[m]), rP[m], true);
      ls.target.Set(ld, x0 - h);
      for (size_t m = 0; m < segs.size(); m++) updateSegment(ld->GetSegment(segs[m]));
      for (size_t m = 0; m < segs.size(); m++) segmentResiduals(ld->GetSegment(segs[m]), rM[m], true);
      ls.target.Set(ld, d0);
      for (size_t m = 0; m < segs.size(); m++) updateSegment(ld->GetSegment(segs[m]));
      for (size_t m = 0; m < segs.size(); m++) {
        const int sNum = segs[m];
        const size_t r0 = (size_t)rowStart[sNum], nSegRows = (size_t)(rowStart[sNum + 1] - rowStart[sNum]);
        if (rP[m].size() == nSegRows && rM[m].size() == nSegRows)
          for (size_t r = 0; r < nSegRows; r++)
            jac[(r0 + r) * (size_t)nCols + (size_t)ls.packed] = (rP[m][r] - rM[m][r]) / (2.0 * h);
      }
    }
    nShiftLocal = (int)localCols.size();

    for (int f : fullShifts) {
      const int packed = pmap.FullToPacked(f);
      const double x0 = full[f];
      const double h = 1.0e-6 * (std::fabs(x0) + 1.0);
      nShiftFull++;
      auto residualsAt = [&](double value, vector_r &out) -> bool {
        vector_r fk = full;
        fk[f] = value;
        lc->FillCompoundFromParams(fk);
        ld->FillNormsFromParams(fk);
        ld->FillEnergyShiftsFromParams(fk, ld, lc, &configure());
        out.clear();
        jTmp.clear();
        return ComputeResidualJacobian(lc, ld, configure(), pmap, sdp, out, jTmp, nc);
      };
      const bool okp = residualsAt(x0 + h, rPlus);
      const bool okm = residualsAt(x0 - h, rMinus);
      if (okp && okm && rPlus.size() == nRes && rMinus.size() == nRes)
        for (size_t r = 0; r < nRes; r++)
          jac[r * (size_t)nCols + (size_t)packed] = (rPlus[r] - rMinus[r]) / (2.0 * h);
    }
    lc->FillCompoundFromParams(full);
    ld->FillNormsFromParams(full);
    ld->FillEnergyShiftsFromParams(full, ld, lc, &configure());
  }
  if (std::getenv("AZURE_GRAD_TIMING")) {
    static std::atomic<int> rjCalls{0};
    const int c = ++rjCalls;
    if (c <= 5 || c % 50 == 0) {
      auto ms = [](auto a, auto b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
      std::cerr << "[grad-timing] residual Jacobian " << c << ": adjoint rows " << ms(tj0, tj1) << " ms, energy-shift columns "
                << ms(tj1, std::chrono::steady_clock::now()) << " ms (" << nShiftLocal << " segment-local, " << nShiftFull
                << " full)" << std::endl;
    }
  }

  // Measured scattering lengths: one row per constrained part (Re, Im) after the data rows, in the
  // order of ScatteringLength::Residuals, with exact rows from the level-matrix adjoint at each
  // threshold.  They depend on no energy shift (zero columns there).
  if (ok && lc->HasScatteringLengthConstraints()) {
    if (configure().paramMask & Config::USE_BRUNE_FORMALISM) lc->CalcShiftFunctions(configure());
    const std::vector<ScatteringLength::Constraint> &cons = lc->ScatteringLengthConstraints();
    vector_r rs;
    ScatteringLength::Residuals(lc, configure(), cons, rs);
    for (size_t q = 0; q < rs.size() && ok; q++) {
      GradAccum acc;
      acc.Init(lc);
      if (!ScatteringLength::Chi2Gradient(lc, configure(), cons, acc, sdp, (int)q)) {
        ok = false;
        break;
      }
      vector_r rowFull(pmap.NumFull(), 0.0);
      acc.Scatter(pmap, rowFull);
      residuals.push_back(rs[q]);
      const size_t base = jac.size();
      jac.resize(base + (size_t)nCols, 0.0);
      for (int f = 0; f < pmap.NumFull(); f++) {
        const int a = pmap.FullToPacked(f);
        if (a >= 0 && a < nCols) jac[base + (size_t)a] = rowFull[f];
      }
    }
  }

  packedToFull.resize(pmap.NumPacked());
  for (int a = 0; a < pmap.NumPacked(); a++) packedToFull[a] = pmap.PackedToFull(a);

  ReturnPooledCNuc(lc);
  ReturnPooledEData(ld);
  return ok;
}

double AZURECalc::RunLevenbergMarquardt(AZUREParams &params, int maxIter,
                                        BandCovariance *bandCovOut) const {
  ROOT::Minuit2::MnUserParameters &mp = params.GetMinuitParams();
  const int nFull = mp.Params().size();
  vector_r full(nFull);
  for (int i = 0; i < nFull; i++) full[i] = mp.Value(i);

  // First Jacobian fixes the column <-> parameter mapping; also collects the
  // free-parameter start values, projection limits, and Gaussian priors.
  vector_r r, J, x;
  std::vector<int> p2f;
  std::vector<double> lo, hi, pen_nom, pen_inv2;
  const int nFree = PrepareFreeParams(full, mp, r, J, p2f, x, lo, hi, pen_nom, pen_inv2);
  if (nFree < 0) return -1.0;  // analytic Jacobian unsupported -> caller falls back
  if (nFree == 0) return Chi2Value(full);

  double cost = Chi2Value(full);
  double lambda = 1.0e-3;

  std::vector<double> H(nFree * nFree), g(nFree), Hl(nFree * nFree), dx(nFree);

  for (int iter = 0; iter < maxIter; iter++) {
    if (iter > 0 && !ResidualJacobian(full, r, J, p2f)) break;
    // MIGRAD reaches WriteIterationOutput through operator(); this path
    // evaluates through the side-effect-free Chi2Value and so must drive the
    // snapshots itself.
    if (iter > 0 && iter % kIterationOutputInterval == 0) WriteIterationOutput(full);
    const int nRes = (int)r.size();

    // Gauss-Newton normal equations: H = J^T J (+ penalty diag), g = J^T r (+ penalty grad).
    std::fill(H.begin(), H.end(), 0.0);
    std::fill(g.begin(), g.end(), 0.0);
    for (int i = 0; i < nRes; i++) {
      const double *Ji = &J[(size_t)i * nFree];
      double ri = r[i];
      for (int a = 0; a < nFree; a++) {
        g[a] += Ji[a] * ri;
        double Jia = Ji[a];
        double *Ha = &H[(size_t)a * nFree];
        for (int b = a; b < nFree; b++) Ha[b] += Jia * Ji[b];
      }
    }
    for (int a = 0; a < nFree; a++) {
      for (int b = 0; b < a; b++) H[(size_t)a * nFree + b] = H[(size_t)b * nFree + a];  // symmetrize
      if (pen_inv2[a] != 0.0) {
        H[(size_t)a * nFree + a] += pen_inv2[a];
        g[a] += (x[a] - pen_nom[a]) * pen_inv2[a];
      }
    }

    // Convergence on the (projected) gradient.
    double gmax = 0.0;
    for (int a = 0; a < nFree; a++) gmax = std::max(gmax, std::fabs(g[a]));
    if (gmax < 1.0e-8) break;

    // Inner loop: damp until a step decreases the cost.
    bool accepted = false;
    for (int tries = 0; tries < 40 && !accepted; tries++) {
      // (H + lambda*diag(H)) dx = -g
      Hl = H;
      for (int a = 0; a < nFree; a++) Hl[(size_t)a * nFree + a] += lambda * H[(size_t)a * nFree + a];

      gsl_matrix *A = gsl_matrix_alloc(nFree, nFree);
      for (int a = 0; a < nFree; a++)
        for (int b = 0; b < nFree; b++) gsl_matrix_set(A, a, b, Hl[(size_t)a * nFree + b]);
      gsl_vector *bvec = gsl_vector_alloc(nFree);
      for (int a = 0; a < nFree; a++) gsl_vector_set(bvec, a, -g[a]);
      gsl_vector *sol = gsl_vector_alloc(nFree);

      gsl_set_error_handler_off();
      int status = gsl_linalg_cholesky_decomp1(A);
      if (status == 0) status = gsl_linalg_cholesky_solve(A, bvec, sol);

      bool solved = (status == 0);
      if (solved)
        for (int a = 0; a < nFree; a++) dx[a] = gsl_vector_get(sol, a);
      gsl_matrix_free(A);
      gsl_vector_free(bvec);
      gsl_vector_free(sol);

      if (!solved) {
        lambda *= 4.0;
        if (lambda > 1.0e14) {
          accepted = false;
          break;
        }
        continue;
      }

      vector_r xnew(nFree), fullNew = full;
      for (int a = 0; a < nFree; a++) {
        xnew[a] = std::min(std::max(x[a] + dx[a], lo[a]), hi[a]);  // clamp to limits
        fullNew[p2f[a]] = xnew[a];
      }
      double costNew = Chi2Value(fullNew);

      if (costNew < cost && std::isfinite(costNew)) {
        double rel = (cost - costNew) / std::max(std::fabs(cost), 1.0);
        full = fullNew;
        x = xnew;
        cost = costNew;
        lambda = std::max(lambda * 0.3, 1.0e-12);
        accepted = true;
        configure().outStream << "\r\tLM iteration: " << std::setw(4) << iter + 1
                              << "  Chi-Squared: " << cost << "        ";
        configure().outStream.flush();
        if (rel < 1.0e-9) iter = maxIter;  // converged
      } else {
        lambda *= 4.0;
        if (lambda > 1.0e14) break;
      }
    }
    if (!accepted) break;
  }
  configure().outStream << std::endl;

  // Write best-fit values back, and parameter errors from (J^T J + penalties)^{-1}.
  for (int i = 0; i < nFull; i++) mp.SetValue(i, full[i]);
  FinalizeLeastSquaresCovariance(full, p2f, pen_inv2, mp, bandCovOut);

  return cost;
}

int AZURECalc::PrepareFreeParams(const vector_r &full,
                                 const ROOT::Minuit2::MnUserParameters &mp,
                                 vector_r &r, vector_r &J, std::vector<int> &p2f,
                                 vector_r &x, std::vector<double> &lo,
                                 std::vector<double> &hi, std::vector<double> &pen_nom,
                                 std::vector<double> &pen_inv2) const {
  // First Jacobian fixes the column <-> parameter mapping.
  if (!ResidualJacobian(full, r, J, p2f)) return -1;  // unsupported -> caller falls back
  const int nFree = (int)p2f.size();

  // Free-parameter values, limits, and quadratic-penalty (prior) terms.
  x.assign(nFree, 0.0);
  lo.assign(nFree, -std::numeric_limits<double>::infinity());
  hi.assign(nFree, std::numeric_limits<double>::infinity());
  pen_nom.assign(nFree, 0.0);
  pen_inv2.assign(nFree, 0.0);
  // Shift priors by full parameter index (p2f holds full indices of the same map).
  std::map<int, std::pair<double, double>> shiftPrior;  // full index -> (nominal, sigma)
  ForEachShiftPenalty(data(), BuildParamIndexMap(compound(), data(), std::vector<bool>()),
                      [&](int idx, double nominal, double sig) { shiftPrior[idx] = {nominal, sig}; });
  for (int a = 0; a < nFree; a++) {
    int f = p2f[a];
    x[a] = full[f];
    const auto &par = mp.Parameter(f);
    if (par.HasLowerLimit()) lo[a] = par.LowerLimit();
    if (par.HasUpperLimit()) hi[a] = par.UpperLimit();
    // Gaussian penalties contribute (1/sigma^2) to the diagonal of J^T J and
    // (x-nominal)/sigma^2 to the gradient.  Norm / shift / nuisance penalties
    // are identified the same way operator()/CalculateNuisanceChiSquared do.
    std::string name = par.GetName();
    if (name.find("norm") != std::string::npos) {
      for (int s = 1; s <= data()->NumSegments(); s++) {
        ESegment *seg = data()->GetSegment(s);
        if (!seg || !seg->IsVaryNorm()) continue;
        char vn[64];
        snprintf(vn, sizeof(vn), "segment_%d_norm", seg->GetSegmentKey());
        if (name == vn) {
          double n0 = seg->GetNominalNorm();
          double sig = n0 / 100.0 * seg->GetNormError();
          if (sig != 0.0) {
            pen_nom[a] = n0;
            pen_inv2[a] = 1.0 / (sig * sig);
          }
          break;
        }
      }
    } else if (name.find("shift") != std::string::npos) {
      auto it = shiftPrior.find(f);
      if (it != shiftPrior.end()) {
        pen_nom[a] = it->second.first;
        pen_inv2[a] = 1.0 / (it->second.second * it->second.second);
      }
    } else if (name.find("screening_") == 0) {
      // Electron-screening potential: its prior is a diagonal penalty row.
      for (int k = 1; k <= data()->NumScreeningParams(); k++) {
        if (name != "screening_" + std::to_string(data()->GetScreeningEffect(k))) continue;
        const TargetEffect *te = data()->GetTargetEffect(data()->GetScreeningEffect(k));
        if (te->GetScreeningSigma() > 0.0) {
          pen_nom[a] = te->GetScreeningNominal();
          pen_inv2[a] = 1.0 / (te->GetScreeningSigma() * te->GetScreeningSigma());
        }
      }
    } else if (name.find("abundance_") == 0) {
      // Isotopic abundance: its <abundances> prior is a diagonal penalty row.
      // The sum-to-one term couples the abundances and is not a per-parameter
      // penalty: the least-squares drivers leave it out (warned once), MIGRAD
      // carries it exactly (AZURECalc::operator()).
      bool anySum = false;
      for (int k = 1; k <= compound()->NumAbundanceParams(); k++) {
        int sample, nucleus;
        compound()->AbundanceOf(k, sample, nucleus);
        if (compound()->SampleSumSigma(sample, configure().abundanceSumSigma) > 0.0) anySum = true;
        if (name != compound()->AbundanceParamName(sample, nucleus)) continue;
        const SampleAbundance &e = compound()->GetSample(sample)->entries[nucleus - 1];
        if (e.sigma > 0.0) {
          pen_nom[a] = e.nominal;
          pen_inv2[a] = 1.0 / (e.sigma * e.sigma);
        }
      }
      static std::atomic<bool> warnedSum{false};
      if (anySum && !warnedSum.exchange(true))
        configure().outStream << "WARNING: a sum-to-one width is set for the abundances, but the least-squares fitter applies only "
                              << "the per-nucleus abundance priors; the sum-to-one penalty is carried by MIGRAD only." << std::endl;
    } else if (limitsManager_ && limitsManager_->IsNuisanceParameterByIndex(a)) {
      double sig = limitsManager_->GetConvertedErrorByIndex(a);
      if (sig > 0.0) {
        pen_nom[a] = limitsManager_->GetConvertedNominalValueByIndex(a);
        pen_inv2[a] = 1.0 / (sig * sig);
      }
    }
  }
  return nFree;
}

void AZURECalc::FinalizeLeastSquaresCovariance(const vector_r &full,
                                               const std::vector<int> &p2f,
                                               const std::vector<double> &pen_inv2,
                                               ROOT::Minuit2::MnUserParameters &mp,
                                               BandCovariance *bandCovOut) const {
  const int nFree = (int)p2f.size();
  vector_r r, J;
  std::vector<int> p2fLocal;
  if (!ResidualJacobian(full, r, J, p2fLocal)) return;

  std::vector<double> H(nFree * nFree, 0.0);
  const int nRes = (int)r.size();
  for (int i = 0; i < nRes; i++) {
    const double *Ji = &J[(size_t)i * nFree];
    for (int a = 0; a < nFree; a++) {
      double *Ha = &H[(size_t)a * nFree];
      for (int b = a; b < nFree; b++) Ha[b] += Ji[a] * Ji[b];
    }
  }
  // Pure J^T J diagonal = each parameter's data sensitivity (column-norm^2 in
  // J), captured before priors are folded into H.  Near zero => unconstrained.
  std::vector<double> jtjDiag(nFree);
  for (int a = 0; a < nFree; a++) {
    for (int b = 0; b < a; b++) H[(size_t)a * nFree + b] = H[(size_t)b * nFree + a];
    jtjDiag[a] = H[(size_t)a * nFree + a];
    if (pen_inv2[a] != 0.0) H[(size_t)a * nFree + a] += pen_inv2[a];
  }

  // The Minuit name alone ("width_1_2") does not say which level or channel is
  // meant, so pair it with the level's J^pi and energy and the channel's pair,
  // L and S.
  auto describe = [&](int a) {
    const int minuitIndex = p2f[a];
    std::string text = mp.GetName(minuitIndex);
    const std::string label = AZURELabel::Parameter(compound(), data(), minuitIndex);
    if (!label.empty()) text += "  =  " + label;
    return text;
  };

  // Name the free parameters the data barely constrain (near-zero J column).
  auto reportInsensitive = [&]() {
    double mx = 0.0;
    for (int a = 0; a < nFree; a++) mx = std::max(mx, jtjDiag[a]);
    int nInsensitive = 0;
    for (int a = 0; a < nFree; a++) {
      if (mx <= 0.0 || jtjDiag[a] <= 1.e-10 * mx) {
        if (!nInsensitive)
          configure().outStream << "  Data-insensitive parameter(s) -- fix, constrain, or add "
                                   "a prior to one of these:"
                                << std::endl;
        configure().outStream << "    " << describe(a) << std::endl;
        nInsensitive++;
      }
    }
    if (!nInsensitive)
      configure().outStream << "No single parameter is insensitive; the flat direction is a "
                               "degeneracy between two or more correlated parameters."
                            << std::endl;
  };

  // Covariance = (J^T J + priors)^{-1}.  If it is singular (unconstrained or
  // degenerate parameter), retry with a small ridge so the flat direction gets
  // a large finite variance instead of losing the errors and band; warn.
  double maxDiag = 0.0;
  for (int a = 0; a < nFree; a++) maxDiag = std::max(maxDiag, H[(size_t)a * nFree + a]);
  const double ridges[4] = {0.0, maxDiag * 1.e-9, maxDiag * 1.e-6, maxDiag * 1.e-3};
  gsl_matrix *A = gsl_matrix_alloc(nFree, nFree);
  gsl_set_error_handler_off();
  int used = -1;
  for (int t = 0; t < 4; t++) {
    for (int a = 0; a < nFree; a++)
      for (int b = 0; b < nFree; b++)
        gsl_matrix_set(A, a, b, H[(size_t)a * nFree + b] + (a == b ? ridges[t] : 0.0));
    if (gsl_linalg_cholesky_decomp1(A) == 0 && gsl_linalg_cholesky_invert(A) == 0) {
      used = t;
      break;
    }
  }
  if (used < 0) {
    configure().outStream << "Warning: the fit curvature matrix (J^T J) is singular; parameter errors "
                             "and the uncertainty band are unavailable."
                          << std::endl;
    reportInsensitive();
  } else {
    if (used > 0) {
      configure().outStream << "Warning: the fit curvature matrix was near-singular; regularized it "
                               "with a small ridge to obtain the covariance. Uncertainties along the "
                               "poorly-constrained direction(s) are large and approximate."
                            << std::endl;
      reportInsensitive();
      // Strongly correlated pairs (|rho| near 1) are the degenerate combinations.
      struct CorrPair {
        double c;
        int a, b;
      };
      std::vector<CorrPair> pairs;
      for (int a = 0; a < nFree; a++) {
        double va = gsl_matrix_get(A, a, a);
        if (va <= 0.0) continue;
        for (int b = a + 1; b < nFree; b++) {
          double vb = gsl_matrix_get(A, b, b);
          if (vb <= 0.0) continue;
          double c = gsl_matrix_get(A, a, b) / std::sqrt(va * vb);
          if (std::fabs(c) > 0.95) pairs.push_back({c, a, b});
        }
      }
      std::sort(pairs.begin(), pairs.end(),
                [](const CorrPair &x, const CorrPair &y) { return std::fabs(x.c) > std::fabs(y.c); });
      for (size_t k = 0; k < pairs.size() && k < 5; k++)
        configure().outStream << "  Correlated (rho=" << pairs[k].c << "):" << std::endl
                              << "    " << describe(pairs[k].a) << std::endl
                              << "    " << describe(pairs[k].b) << std::endl;
    }
    for (int a = 0; a < nFree; a++) {
      double v = gsl_matrix_get(A, a, a);
      if (v > 0.0) mp.SetError(p2f[a], std::sqrt(v));
    }
    // Name free parameters whose error dwarfs their value (they dominate the band).
    {
      int n = 0;
      for (int a = 0; a < nFree; a++) {
        double v = gsl_matrix_get(A, a, a);
        if (v <= 0.0) continue;
        double val = mp.Value(p2f[a]);
        if (std::sqrt(v) > 10.0 * std::max(std::fabs(val), 1.e-30)) {
          if (!n)
            configure().outStream << "Note: weakly-determined parameter(s) (error > 1000% of value) "
                                     "dominate the uncertainty band. Constrain or fix these to "
                                     "tighten it:"
                                  << std::endl;
          configure().outStream << "    " << describe(a) << std::endl;
          n++;
        }
      }
    }
    // Export the full covariance for cross-section bands.  Columns are the
    // packed free parameters, in the same order as p2f, so we tag them with
    // their parameter identities via a matching index map.
    if (bandCovOut) {
      AZUREParams tp;
      compound()->FillMnParams(tp.GetMinuitParams(), &configure());
      data()->FillMnParams(tp.GetMinuitParams());
      int nMn = tp.GetMinuitParams().Params().size();
      std::vector<bool> fixed(nMn);
      for (int i = 0; i < nMn; i++) fixed[i] = tp.GetMinuitParams().Parameter(i).IsFixed();
      ParamIndexMap pmap = BuildParamIndexMap(compound(), data(), fixed);
      if (pmap.NumPacked() == nFree) {
        // Keep only the R-matrix sub-block: the band is insensitive to norms
        // and energy shifts, so they are dropped from the saved covariance.
        const std::vector<int> rc = RMatrixPackedColumns(pmap);
        const int m = (int)rc.size();
        bandCovOut->cols.resize(m);
        for (int a = 0; a < m; a++) bandCovOut->cols[a] = pmap.Desc(pmap.PackedToFull(rc[a]));
        bandCovOut->M.assign(m, std::vector<double>(m, 0.0));
        for (int a = 0; a < m; a++)
          for (int b = 0; b < m; b++) bandCovOut->M[a][b] = gsl_matrix_get(A, rc[a], rc[b]);
      }
    }
  }
  gsl_matrix_free(A);
}

bool AZURECalc::ResidualsOnly(const vector_r &full, vector_r &residuals) const {
  // Forward-only counterpart of ResidualJacobian (no adjoint): same row order
  // and residual r = (model - data*n)/(cmErr*n) as ComputeResidualJacobian.
  // Pooled working copies keep their applied energy shifts and mapping between calls; a fresh
  // Clone() re-shifted every shifted segment and rebuilt the mapping each time (~1 s here).
  CNuc *lc = GetPooledCNuc();
  EData *ld = GetPooledEData();
  lc->FillCompoundFromParams(full);
  ld->FillNormsFromParams(full);
  ld->FillEnergyShiftsFromParams(full, ld, lc, &configure());
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) lc->CalcShiftFunctions(configure(), false);

  residuals.clear();
  for (int i = 1; i <= ld->NumSegments(); i++) {
    ESegment *segment = ld->GetSegment(i);
    if (!segment) continue;
    const double norm = segment->GetNorm();
    for (int pid = 0; pid < segment->NumPoints(); pid++) {
      EPoint *pt = segment->GetPoint(pid + 1);
      if (!pt) continue;
      double th = segment->CalculateTheoreticalCrossSection(pid, lc, configure(), ld);
      double denom = pt->GetCMCrossSectionError() * norm;
      residuals.push_back(denom != 0.0 ? (th - pt->GetCMCrossSection() * norm) / denom : 0.0);
    }
  }
  // Measured scattering lengths, as the last rows of ResidualJacobian.
  if (lc->HasScatteringLengthConstraints())
    ScatteringLength::Residuals(lc, configure(), lc->ScatteringLengthConstraints(), residuals);
  ReturnPooledCNuc(lc);
  ReturnPooledEData(ld);
  return true;
}

namespace {

// State shared with the GSL least-squares callbacks; owned by RunGSLNonlinear.
struct GSLNlsData {
  const AZURECalc *self;
  int nFree;
  int nData;                    // number of data residual rows
  const std::vector<int> *p2f;  // packed free param -> full-vector index
  const std::vector<double> *lo;
  const std::vector<double> *hi;
  const std::vector<double> *pen_nom;
  const std::vector<double> *pen_inv2;
  const std::vector<int> *penRows;  // free-param indices that carry a prior
  vector_r baseFull;                // full vector holding the fixed params
};

// Project x into [lo,hi] and scatter it into the full parameter vector.
void gslExpandFull(const GSLNlsData *d, const gsl_vector *x, vector_r &full) {
  full = d->baseFull;
  for (int a = 0; a < d->nFree; a++) {
    double v = gsl_vector_get(x, a);
    v = std::min(std::max(v, (*d->lo)[a]), (*d->hi)[a]);
    full[(*d->p2f)[a]] = v;
  }
}

// Residual vector: data rows plus Gaussian-prior rows sqrt(1/sigma^2)*(x-nom).
int gslResidualF(const gsl_vector *x, void *params, gsl_vector *f) {
  const GSLNlsData *d = static_cast<const GSLNlsData *>(params);
  vector_r full;
  gslExpandFull(d, x, full);
  vector_r res;
  if (!d->self->ResidualsOnly(full, res) || (int)res.size() != d->nData) return GSL_EBADFUNC;
  for (int i = 0; i < d->nData; i++) gsl_vector_set(f, i, res[i]);
  for (size_t k = 0; k < d->penRows->size(); k++) {
    int a = (*d->penRows)[k];
    double s = std::sqrt((*d->pen_inv2)[a]);
    gsl_vector_set(f, d->nData + (int)k, s * (full[(*d->p2f)[a]] - (*d->pen_nom)[a]));
  }
  return GSL_SUCCESS;
}

// Residual Jacobian: analytic data-row block plus the constant prior rows.
int gslResidualDf(const gsl_vector *x, void *params, gsl_matrix *Jm) {
  const GSLNlsData *d = static_cast<const GSLNlsData *>(params);
  vector_r full;
  gslExpandFull(d, x, full);
  vector_r r, Jflat;
  std::vector<int> p2f;
  if (!d->self->ResidualJacobian(full, r, Jflat, p2f)) return GSL_EBADFUNC;
  if ((int)p2f.size() != d->nFree || (int)r.size() != d->nData) return GSL_EBADFUNC;
  const int nFree = d->nFree;
  for (int i = 0; i < d->nData; i++)
    for (int a = 0; a < nFree; a++)
      gsl_matrix_set(Jm, i, a, Jflat[(size_t)i * nFree + a]);
  for (size_t k = 0; k < d->penRows->size(); k++) {
    int row = d->nData + (int)k;
    for (int a = 0; a < nFree; a++) gsl_matrix_set(Jm, row, a, 0.0);
    int a = (*d->penRows)[k];
    gsl_matrix_set(Jm, row, a, std::sqrt((*d->pen_inv2)[a]));
  }
  return GSL_SUCCESS;
}

// Per-iteration progress line (chi^2 = ||residual||^2, priors included).
void gslProgress(const size_t iter, void *params,
                 const gsl_multifit_nlinear_workspace *w) {
  const GSLNlsData *d = static_cast<const GSLNlsData *>(params);
  gsl_vector *f = gsl_multifit_nlinear_residual(w);
  double chi2 = 0.0;
  gsl_blas_ddot(f, f, &chi2);
  d->self->configure().outStream << "\r\tGSL-LM iteration: " << std::setw(4) << iter
                                 << "  Chi-Squared: " << chi2 << "        ";
  d->self->configure().outStream.flush();

  // Intermediate output, so a long GSL fit also leaves usable files behind
  // before it converges rather than only at the end.
  if (iter > 0 && iter % AZURECalc::kIterationOutputInterval == 0) {
    vector_r full = d->baseFull;
    const gsl_vector *x = gsl_multifit_nlinear_position(w);
    for (int a = 0; a < d->nFree; a++) full[(*d->p2f)[a]] = gsl_vector_get(x, a);
    d->self->WriteIterationOutput(full);
  }
}

}  // namespace

double AZURECalc::RunGSLNonlinear(AZUREParams &params, int maxIter,
                                  BandCovariance *bandCovOut) const {
  ROOT::Minuit2::MnUserParameters &mp = params.GetMinuitParams();
  const int nFull = mp.Params().size();
  vector_r full(nFull);
  for (int i = 0; i < nFull; i++) full[i] = mp.Value(i);

  vector_r r0, J0, x;
  std::vector<int> p2f;
  std::vector<double> lo, hi, pen_nom, pen_inv2;
  const int nFree = PrepareFreeParams(full, mp, r0, J0, p2f, x, lo, hi, pen_nom, pen_inv2);
  if (nFree < 0) return -1.0;  // analytic Jacobian unsupported -> caller falls back
  if (nFree == 0) return Chi2Value(full);
  const int nData = (int)r0.size();

  // Gaussian priors on norms / shifts / nuisance params become extra residual rows.
  std::vector<int> penRows;
  for (int a = 0; a < nFree; a++)
    if (pen_inv2[a] > 0.0) penRows.push_back(a);
  const int nRes = nData + (int)penRows.size();

  GSLNlsData d;
  d.self = this;
  d.nFree = nFree;
  d.nData = nData;
  d.p2f = &p2f;
  d.lo = &lo;
  d.hi = &hi;
  d.pen_nom = &pen_nom;
  d.pen_inv2 = &pen_inv2;
  d.penRows = &penRows;
  d.baseFull = full;

  gsl_multifit_nlinear_fdf fdf;
  fdf.f = gslResidualF;
  fdf.df = gslResidualDf;
  fdf.fvv = nullptr;  // second directional derivative via finite differences of f
  fdf.n = nRes;
  fdf.p = nFree;
  fdf.params = &d;

  // Geodesic-accelerated LM in a trust region; the default solver factorizes J
  // directly (QR) rather than forming J^T J.
  gsl_multifit_nlinear_parameters gp = gsl_multifit_nlinear_default_parameters();
  gp.trs = gsl_multifit_nlinear_trs_lmaccel;

  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_workspace *w = gsl_multifit_nlinear_alloc(T, &gp, nRes, nFree);
  if (!w) return -1.0;

  gsl_vector *x0 = gsl_vector_alloc(nFree);
  for (int a = 0; a < nFree; a++) gsl_vector_set(x0, a, x[a]);

  gsl_set_error_handler_off();
  gsl_multifit_nlinear_init(x0, &fdf, w);

  const double xtol = 1.0e-9, gtol = 1.0e-9, ftol = 1.0e-9;
  int info = 0;
  int status = gsl_multifit_nlinear_driver(maxIter, xtol, gtol, ftol,
                                           gslProgress, &d, &info, w);
  configure().outStream << std::endl;

  std::string reason = "reached the maximum number of iterations";
  if (status == GSL_SUCCESS)
    reason = (info == 1) ? "converged (small step)" : "converged (small gradient)";
  else if (status != GSL_EMAXITER)
    reason = std::string("stopped: ") + gsl_strerror(status);
  configure().outStream << "GSL trust-region solver " << reason << " after "
                        << gsl_multifit_nlinear_niter(w) << " iterations." << std::endl;

  // Read back the (projected) best-fit parameters.
  gsl_vector *xf = gsl_multifit_nlinear_position(w);
  for (int a = 0; a < nFree; a++) {
    double v = gsl_vector_get(xf, a);
    v = std::min(std::max(v, lo[a]), hi[a]);
    full[p2f[a]] = v;
  }
  for (int i = 0; i < nFull; i++) mp.SetValue(i, full[i]);

  gsl_vector_free(x0);
  gsl_multifit_nlinear_free(w);

  double cost = Chi2Value(full);

  // Errors and band covariance from (J^T J + priors)^{-1}, as the LM path does.
  FinalizeLeastSquaresCovariance(full, p2f, pen_inv2, mp, bandCovOut);

  return cost;
}

void AZURECalc::AddNuisanceGradient(const vector_r &p, std::vector<double> &grad) const {
  if (!limitsManager_) return;

  AZUREParams tempParams;
  compound()->FillMnParams(tempParams.GetMinuitParams(), &configure());
  data()->FillMnParams(tempParams.GetMinuitParams());

  std::vector<int> nonFixedToActualIndex;
  for (int i = 0; i < tempParams.GetMinuitParams().Params().size(); i++) {
    if (!tempParams.GetMinuitParams().Parameter(i).IsFixed() ||
        tempParams.GetMinuitParams().Parameter(i).GetName().find("segment") != std::string::npos) {
      nonFixedToActualIndex.push_back(i);
    }
  }

  for (int nf = 0; nf < (int)nonFixedToActualIndex.size() && nf < (int)p.size(); nf++) {
    int actualIndex = nonFixedToActualIndex[nf];
    std::string paramName = tempParams.GetMinuitParams().Parameter(actualIndex).GetName();
    // Norms, shifts and abundances carry their own priors (segment errors, the <abundances> sigma).
    if (paramName.find("norm") != std::string::npos || paramName.find("shift") != std::string::npos ||
        paramName.find("abundance") != std::string::npos || paramName.find("screening") != std::string::npos) continue;
    if (!limitsManager_->IsNuisanceParameterByIndex(nf)) continue;
    double nominalValue = limitsManager_->GetConvertedNominalValueByIndex(nf);
    double paramError = limitsManager_->GetConvertedErrorByIndex(nf);
    if (paramError > 0.0) {
      // d/d p[nf] of ((p[nf]-nominal)/error)^2.
      grad[nf] += 2.0 * (p[nf] - nominalValue) / (paramError * paramError);
    }
  }
}

double AZURECalc::CalculateNuisanceChiSquared(const vector_r &p) const {
  double nuisanceChiSquared = 0.0;

  // Create temporary AZUREParams to get parameter names
  AZUREParams tempParams;
  compound()->FillMnParams(tempParams.GetMinuitParams(), &configure());
  data()->FillMnParams(tempParams.GetMinuitParams());

  // Build mapping from non-fixed parameter index to actual parameter index
  std::vector<int> nonFixedToActualIndex;
  for (int i = 0; i < tempParams.GetMinuitParams().Params().size(); i++) {
    if (!tempParams.GetMinuitParams().Parameter(i).IsFixed() || tempParams.GetMinuitParams().Parameter(i).GetName().find("segment") != std::string::npos) {
      nonFixedToActualIndex.push_back(i);
    }
  }

  // Check each non-fixed parameter to see if it's marked as nuisance
  for (int nonFixedIndex = 0; nonFixedIndex < nonFixedToActualIndex.size() && nonFixedIndex < p.size(); nonFixedIndex++) {
    int actualIndex = nonFixedToActualIndex[nonFixedIndex];
    std::string paramName = tempParams.GetMinuitParams().Parameter(actualIndex).GetName();

    // If norm, shift or abundance in param name, skip: those carry their own priors
    if (paramName.find("norm") != std::string::npos || paramName.find("shift") != std::string::npos ||
        paramName.find("abundance") != std::string::npos || paramName.find("screening") != std::string::npos) {
      continue;
    }

    // First check if this parameter is marked as nuisance (fast check)
    if (!limitsManager_->IsNuisanceParameterByIndex(nonFixedIndex)) {
      continue;  // Skip if not a nuisance parameter
    }

    // Only do expensive conversions if parameter is marked as nuisance
    double nominalValue = limitsManager_->GetConvertedNominalValueByIndex(nonFixedIndex);
    double paramError = limitsManager_->GetConvertedErrorByIndex(nonFixedIndex);

    // If we got valid values (non-zero error means this parameter has valid nuisance settings)
    if (paramError > 0.0) {
      double paramValue = p[nonFixedIndex];
      double deviation = (paramValue - nominalValue) / paramError;
      nuisanceChiSquared += deviation * deviation;
    }
  }

  return nuisanceChiSquared;
}

/*!
 * Initialize object pools with pre-allocated CNuc and EData objects
 */
void AZURECalc::InitializePools() const {
  std::lock_guard<std::mutex> lock(pool_mutex_);
  if (pools_initialized_) return;

  // Calculate pool size based on OpenMP threads (fixes interaction with OpenMP)
  int pool_size = 4;  // default minimum
#ifdef _OPENMP
  pool_size = std::max(4, omp_get_max_threads());
#else
  pool_size = std::max(4, static_cast<int>(std::thread::hardware_concurrency()));
#endif

  // Pre-allocate CNuc objects by cloning once
  for (int i = 0; i < pool_size; ++i) {
    cnuc_pool_.push(std::unique_ptr<CNuc>(compound()->Clone()));
  }

  // Pre-allocate EData objects by cloning once
  for (int i = 0; i < pool_size; ++i) {
    edata_pool_.push(std::unique_ptr<EData>(data()->Clone()));
  }

  pools_initialized_ = true;
}

/*!
 * Get a CNuc object from the pool, creating new if pool is empty
 */
CNuc *AZURECalc::GetPooledCNuc() const {
  std::lock_guard<std::mutex> lock(pool_mutex_);

  if (!cnuc_pool_.empty()) {
    auto obj = std::move(cnuc_pool_.top());
    cnuc_pool_.pop();
    return obj.release();
  }

  // Fallback: create new if pool is empty (shouldn't happen often)
  return compound_->Clone();
}

/*!
 * Get an EData object from the pool, creating new if pool is empty
 */
EData *AZURECalc::GetPooledEData() const {
  std::lock_guard<std::mutex> lock(pool_mutex_);

  if (!edata_pool_.empty()) {
    auto obj = std::move(edata_pool_.top());
    edata_pool_.pop();
    return obj.release();
  }

  // Fallback: create new if pool is empty (shouldn't happen often)
  return data_->Clone();
}

/*!
 * Return a CNuc object to the pool for reuse
 */
void AZURECalc::ReturnPooledCNuc(CNuc *obj) const {
  if (!obj) return;

  std::lock_guard<std::mutex> lock(pool_mutex_);
  cnuc_pool_.push(std::unique_ptr<CNuc>(obj));
}

/*!
 * Return an EData object to the pool for reuse
 */
void AZURECalc::ReturnPooledEData(EData *obj) const {
  if (!obj) return;

  std::lock_guard<std::mutex> lock(pool_mutex_);
  edata_pool_.push(std::unique_ptr<EData>(obj));
}

/*!
 * Write parameters to file
 */
void AZURECalc::WriteParameters(AZUREParams &params, const Config &configure) const {
  char filename[256];
  snprintf(filename, sizeof(filename), "%sparam.fit", configure.outputdir.c_str());
  std::ofstream out;
  out.open(filename);
  if (out) {
    out.precision(7);
    for (int i = 0; i < params.GetMinuitParams().Params().size(); i++) {
      out << std::setw(20) << params.GetMinuitParams().GetName(i)
          << std::scientific << std::setw(20) << params.GetMinuitParams().Value(i)
          << std::scientific << std::setw(20) << params.GetMinuitParams().Error(i) << std::endl;
    }
    out.flush();
    out.close();
  } else
    configure.outStream << "Could not save param.fit file." << std::endl;
}