// AZUREGrad, THM (HOES) rows: the rows of the residual Jacobian and of the chi^2 gradient of the THM data
// (ComputeTHMRows), with the HOES adjoint of one point and its fold (ThmPointModelGradient).  The rest of the
// gradient (accumulator, parameter-index map, adjoint orchestration) is AZUREGrad.cpp.
#include "AZUREGrad.h"
#include "Objective.h"
#include "RunDiagnostics.h"
#include "azure2/support/Environment.h"
#include "GradEnergyTasks.h"
#include "ModelPass.h"
#include "Profile.h"
#include <atomic>
#include <chrono>
#include <iostream>
#include <cstdlib>
#include "CompoundSystem.h"
#include "CalculationPlan.h"
#include "WorkingCopy.h"
#include "JGroup.h"
#include "ALevel.h"
#include "AChannel.h"
#include "PPair.h"
#include "Segment.h"
#include "Unit.h"
#include "TargetEffect.h"
#include "AMatrixFunc.h"
#include "Config.h"
#include "ChannelFunc.h"
#include "CoulFunc.h"
#include "ShftFunc.h"
#include "THMMatrixFunc.h"
#include "PointFold.h"
#include "azure2/support/FloatingPointPrecision.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <string>
#include <unordered_map>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

// ---------------------------------------------------------------------------
// THM (HOES) model derivative of one point -- a point or a mapped point, with
// its folding if it has one -- by the HOES adjoint (THMMatrixFunc::
// THMPointAdjoint): d m/d(E, gamma) into `accum`, d m/d(coherent value k)
// into coherentBar[k].  A folded point is linear in its sub-point cross
// sections; the weights come from the combiner as in GradTargetEffectAdjoint.
// A folded analyzing power (a ratio of two folds): not covered (false).
bool ThmPointModelGradient(WorkingCopy &data, Unit *point, azure2::Results &results, const CompoundSystem *compound, const LevelState &state, const Config &config,
                           const ThmAdjointContext &ctx, THMMatrixFunc &f,
                           GradAccum &accum, std::vector<double> &coherentBar) {
  f.Reset(compound, state, config);
  auto adjoint = [&](Site *pt, double weight) -> bool {
    f.ClearMatrices();
    f.FillMatrices(pt);
    f.InvertMatrices();
    return f.THMPointAdjoint(pt, weight, accum, ctx, &coherentBar);
  };
  const TargetEffect *te = point->IsTargetEffect() ? data.Data().GetTargetEffect(point->GetTargetEffectNum()) : nullptr;
  if (!(te && te->IsSubPointEffect())) return adjoint(point, 1.0);
  if (point->IsAnalyzingPower()) return false;
  const int nsp = point->NumSubPoints();
  if (nsp <= 0) return true;
  std::vector<double> sigma(nsp), c(nsp, 0.0);
  for (int i = 1; i <= nsp; i++) sigma[i - 1] = point->GetSubPoint(data, i)->GetFitCrossSection(results);
  std::vector<double> v(sigma);  // the slots' records are read, never written (SlotGroup.h)
  for (int i = 1; i <= nsp; i++) {
    const double s0 = sigma[i - 1], h = 1.0e-3 * (std::fabs(s0) + 1.0);
    v[i - 1] = s0 + h;
    const double yp = FoldValues(*point, data, config, v);
    v[i - 1] = s0 - h;
    const double ym = FoldValues(*point, data, config, v);
    v[i - 1] = s0;
    c[i - 1] = (yp - ym) / (2.0 * h);
  }
  for (int i = 1; i <= nsp; i++)
    if (c[i - 1] != 0.0 && !adjoint(point->GetSubPoint(data, i), c[i - 1])) return false;
  return true;
}

}  // namespace

// THM (HOES) rows: see the declaration in AZUREGrad.h.  Shared by the API
// (residual Jacobian, chi2 gradient, model gradients) and the CLI band.
std::vector<THMRows> ComputeTHMRows(const CompoundSystem *lc, LevelState &state, const LevelInputs &inputs, WorkingCopy *ld, azure2::Results &results, const Config &config,
                                    const vector_r &full, const ParamIndexMap &pmap) {
  std::vector<THMRows> rows;
  const bool brune = (config.paramMask & Config::USE_BRUNE_FORMALISM);
  const int nCols = pmap.NumPacked();

  // Row bookkeeping mirrors ComputeResidualJacobian / EvaluateFilledChi2: one
  // row per non-null point, segments in order.
  int row = 0;
  for (int s = 1; s <= ld->Data().NumSegments(); s++) {
    const Segment *seg = ld->Data().GetSegment(s);
    if (!seg) continue;
    if (seg->IsTHM()) {
      THMRows tr;
      tr.segment = s;
      tr.firstRow = row;
      rows.push_back(tr);
    }
    for (int pid = 0; pid < seg->NumPoints(); pid++)
      if (seg->GetPoint(*ld, pid + 1)) row++;
  }
  if (rows.empty()) return rows;

  // Models of every THM point at the currently filled compound.  The points
  // are independent (as in every other forward pass, the compound is only
  // read), so they are evaluated in parallel: this runs 2 x (free E and gamma
  // columns) + 1 times per Jacobian, and serially it left all but one thread
  // idle -- 280 s for the 230 columns of the 12C+12C THM example against a
  // 0.67 s forward pass.
  std::vector<std::pair<int, int>> order;  // (row block t, point index pid)
  for (size_t t = 0; t < rows.size(); t++) {
    const Segment *seg = ld->Data().GetSegment(rows[t].segment);
    for (int pid = 0; pid < seg->NumPoints(); pid++)
      if (seg->GetPoint(*ld, pid + 1)) order.push_back(std::make_pair((int)t, pid));
  }
  auto models = [&](std::vector<vector_r> &out) {
    out.assign(rows.size(), vector_r());
    std::vector<size_t> slot(order.size());
    for (size_t k = 0; k < order.size(); k++) {
      slot[k] = out[order[k].first].size();
      out[order[k].first].push_back(0.0);
    }
    const long n = (long)order.size();
#pragma omp parallel for schedule(dynamic)
    for (long k = 0; k < n; k++) { const azure2::support::ExtendedPrecisionScope fpPrecision;
      const Segment *seg = ld->Data().GetSegment(rows[order[k].first].segment);
      out[order[k].first][slot[k]] =
          azure2::EvaluateModelPoint(*ld, *seg, order[k].second, lc, state, config, results, azure2::LeafValues::Calculate);
    }
  };
  auto fillCompound = [&](const vector_r &p) {
    lc->FillCompoundFromParams(p, state);
    if (brune) lc->CalcShiftFunctions(config, state);
  };

  // Norms and shifts are the caller's; the THM model depends on the R-matrix
  // parameters (and the shift already applied to the point energies).
  ld->FillNormsFromParams(full);
  ld->FillEnergyShiftsFromParams(full, results, lc, &state, &inputs, &config);
  fillCompound(full);
  std::vector<vector_r> m0, mp, mm;
  models(m0);
  for (size_t t = 0; t < rows.size(); t++) {
    rows[t].m = m0[t];
    rows[t].Jm.assign(m0[t].size() * (size_t)nCols, 0.0);
    rows[t].J.assign(m0[t].size() * (size_t)nCols, 0.0);
  }

  // d m / d p, E and gamma columns -- and the THM coherent backgrounds
  // (cbackground=), which only THM points see -- by the HOES adjoint: one
  // reverse sweep per point (ThmPointModelGradient).  A configuration it does
  // not cover, or AZURE2_THM_FD=1, takes the central differences below.
  bool analytic = !ThmFiniteDifferencesRequested(config);
  for (size_t t = 0; t < rows.size() && analytic; t++)
    if (ld->Data().GetSegment(rows[t].segment)->HasComponents()) analytic = false;
  if (analytic) {
    vector_matrix_r shiftDeriv;
    const vector_matrix_r *sdp = nullptr;
    if (brune) {
      shiftDeriv = BuildShiftDerivTable(lc, state, config);
      sdp = &shiftDeriv;
    }
    bool lineshape = false;
    for (size_t t = 0; t < rows.size(); t++)
      if (ld->Data().GetSegment(rows[t].segment)->GetThmLineshape()) lineshape = true;
    const ThmAdjointContext ctx = BuildThmAdjointContext(lc, state, ld->ThmCoherentValues(), config, sdp, lineshape);
    const int nCoherent = ld->Data().NumThmCoherentParams();
    std::vector<size_t> slot(order.size());
    {
      std::vector<size_t> next(rows.size(), 0);
      for (size_t k = 0; k < order.size(); k++) slot[k] = next[order[k].first]++;
    }
    bool bail = false;
    const long n = (long)order.size();
#pragma omp parallel
    { const azure2::support::ExtendedPrecisionScope fpPrecision;
      GradAccum accum;
      accum.Init(lc, state);
      vector_r fullRow(pmap.NumFull(), 0.0);
      std::vector<double> coherentBar(nCoherent, 0.0);
      THMMatrixFunc thm(lc, state, config);  // this thread's solver for the region (re-pointed per point)
#pragma omp for schedule(dynamic)
      for (long k = 0; k < n; k++) {
        if (bail) continue;
        const Segment *seg = ld->Data().GetSegment(rows[order[k].first].segment);
        Unit *pt = seg->GetPoint(*ld, order[k].second + 1);
        accum.Zero();
        std::fill(coherentBar.begin(), coherentBar.end(), 0.0);
        if (!ThmPointModelGradient(*ld, pt, results, lc, state, config, ctx, thm, accum, coherentBar)) {
#pragma omp critical
          bail = true;
          continue;
        }
        std::fill(fullRow.begin(), fullRow.end(), 0.0);
        accum.Scatter(pmap, fullRow);
        double *row = &rows[order[k].first].Jm[slot[k] * nCols];
        for (int c = 0; c < nCols; c++) {
          const int f = pmap.PackedToFull(c);
          const ParamDesc &desc = pmap.Desc(f);
          if (IsRMatrixParameter(desc.kind))
            row[c] = fullRow[f];
          else if (desc.kind == ParamKind::ThmCoherent && desc.channel >= 1 && desc.channel <= nCoherent)
            row[c] = coherentBar[desc.channel - 1];
        }
      }
    }
    if (bail) {
      analytic = false;
      for (size_t t = 0; t < rows.size(); t++) std::fill(rows[t].Jm.begin(), rows[t].Jm.end(), 0.0);
      if (azure2::support::ProcessEnvironment().Set(azure2::support::Switch::kGradDebug))
        std::cerr << "[grad] THM adjoint: a point outside its path (a target effect with an automatic "
                     "tolerance, RMC, ...) -> central differences for the THM rows" << std::endl;
    }
  }

  // Central differences (AZURE2_THM_FD=1, or the adjoint declined).
  for (int c = 0; c < nCols && !analytic; c++) {
    const int f = pmap.PackedToFull(c);
    const ParamKind kind = pmap.Desc(f).kind;
    if (kind == ParamKind::ThmCoherent) {
      const double x0 = full[f];
      const double h = 1.0e-6 * (std::fabs(x0) + 1.0);
      vector_r pp = full;
      pp[f] = x0 + h;
      ld->FillThmCoherentFromParams(pp);
      models(mp);
      pp[f] = x0 - h;
      ld->FillThmCoherentFromParams(pp);
      models(mm);
      ld->FillThmCoherentFromParams(full);
      for (size_t t = 0; t < rows.size(); t++)
        for (size_t i = 0; i < m0[t].size(); i++)
          rows[t].Jm[i * nCols + c] = (mp[t][i] - mm[t][i]) / (2.0 * h);
      continue;
    }
    if (!IsRMatrixParameter(kind)) continue;
    const double x0 = full[f];
    const double h = 1.0e-6 * (std::fabs(x0) + 1.0);
    vector_r pp = full;
    pp[f] = x0 + h;
    fillCompound(pp);
    models(mp);
    pp[f] = x0 - h;
    fillCompound(pp);
    models(mm);
    for (size_t t = 0; t < rows.size(); t++)
      for (size_t i = 0; i < m0[t].size(); i++)
        rows[t].Jm[i * nCols + c] = (mp[t][i] - mm[t][i]) / (2.0 * h);
  }
  fillCompound(full);

  // Residuals and the chain rule through the profiled scale.
  for (size_t t = 0; t < rows.size(); t++) {
    THMRows &tr = rows[t];
    const Segment *seg = ld->Data().GetSegment(tr.segment);
    std::vector<double> d, e, ml, Jml;
    for (int pid = 0; pid < seg->NumPoints(); pid++) {
      Unit *pt = seg->GetPoint(*ld, pid + 1);
      if (!pt) continue;
      d.push_back(seg->LabValue(pid));
      e.push_back(seg->LabError(pid));
      pt->SetModelValue(results, tr.m[d.size() - 1]);  // the residuals below take the model in the datum's frame
      azure2::AppendInDataFrame(*ld, *seg, pid, tr.m[d.size() - 1], tr.Jm.data() + (d.size() - 1) * nCols, nCols, ml, Jml);
    }
    const size_t n = d.size();
    if (ld->Data().ThmGroupOf(tr.segment) >= 0) continue;  // a THM experiment: below
    double Smm = 0.0, Smd = 0.0;
    for (size_t i = 0; i < n; i++) {
      if (e[i] == 0.0) continue;
      double w = 1.0 / (e[i] * e[i]);
      Smm += ml[i] * ml[i] * w;
      Smd += ml[i] * d[i] * w;
    }
    // s multiplies the model (s = 1/n, n the norm on the data).
    bool varyScale = false;
    double s;
    if (seg->IsProfiledNorm()) {
      varyScale = (Smd > 0.0 && Smm > 0.0);
      s = varyScale ? Smd / Smm : 1.0;
    } else {
      double norm = ld->State(*seg).norm;
      s = (norm != 0.0) ? 1.0 / norm : 0.0;
    }
    tr.r.assign(n, 0.0);
    for (size_t i = 0; i < n; i++)
      if (e[i] != 0.0) tr.r[i] = (s * ml[i] - d[i]) / e[i];

    vector_r ds(nCols, 0.0);
    if (varyScale) {
      for (int c = 0; c < nCols; c++) {
        double sdJ = 0.0, smJ = 0.0;
        for (size_t i = 0; i < n; i++) {
          if (e[i] == 0.0) continue;
          double w = 1.0 / (e[i] * e[i]);
          sdJ += d[i] * Jml[i * nCols + c] * w;
          smJ += ml[i] * Jml[i * nCols + c] * w;
        }
        ds[c] = (sdJ - 2.0 * s * smJ) / Smm;
      }
    }
    for (size_t i = 0; i < n; i++) {
      if (e[i] == 0.0) continue;
      for (int c = 0; c < nCols; c++)
        tr.J[i * nCols + c] = (s * Jml[i * nCols + c] + ml[i] * ds[c]) / e[i];
    }
    tr.s = s;
    tr.profiled = varyScale;
    tr.ds = ds;
  }

  // THM experiments: one linear profile (shared norm, background) over the
  // points of all their segments, differentiated through the profile.
  for (int g = 0; g < ld->Data().NumThmGroups(); g++) {
    std::vector<size_t> members;  // indices into rows
    for (size_t t = 0; t < rows.size(); t++)
      if (ld->Data().ThmGroupOf(rows[t].segment) == g) members.push_back(t);
    if (members.empty()) continue;
    ld->ProfileThmGroup(g, results);  // norms and profile from the models set above
    const ThmProfile &prof = results.objective.thmProfiles[g];
    std::vector<double> m, d, e, energy, Jm;
    for (size_t t : members) {
      const Segment *seg = ld->Data().GetSegment(rows[t].segment);
      size_t i = 0;
      for (int pid = 0; pid < seg->NumPoints(); pid++) {
        Unit *pt = seg->GetPoint(*ld, pid + 1);
        if (!pt) continue;
        d.push_back(seg->LabValue(pid));
        e.push_back(seg->LabError(pid));
        energy.push_back(pt->GetCMEnergy());
        azure2::AppendInDataFrame(*ld, *seg, pid, rows[t].m[i], rows[t].Jm.data() + i * nCols, nCols, m, Jm);
        i++;
      }
    }
    std::vector<double> J, G;
    ThmProfileDerivative(prof, m, d, e, energy, Jm, nCols, J, G);
    size_t k = 0;
    for (size_t t : members) {
      THMRows &tr = rows[t];
      const size_t n = tr.m.size();
      tr.r.assign(n, 0.0);
      for (size_t i = 0; i < n; i++) tr.r[i] = prof.Residual(m[k + i], d[k + i], e[k + i], energy[k + i]);
      tr.J.assign(J.begin() + k * nCols, J.begin() + (k + n) * nCols);
      tr.G.assign(G.begin() + k * nCols, G.begin() + (k + n) * nCols);
      tr.s = prof.s;
      tr.profiled = false;  // the whole derivative is in J and G
      tr.ds.assign(nCols, 0.0);
      k += n;
    }
  }
  return rows;
}
