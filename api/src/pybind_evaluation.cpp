// The objective and the evaluation of a session for pyazr's Objective and Evaluation (pyazr/objective.py,
// pyazr/evaluate.py; docs/upgrade/pyazr3.md section 7): one evaluation read as the output writer reads it, the
// objective's least-squares rows, its log-probability, and the table of the numbers a point carries beside its model
// value.  A new point output is one row of PointOutputs(); pyazr.kinds()["data_outputs"] lists the rows and
// Evaluation.data[i].extra holds those that apply to a data set.
#include <algorithm>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <limits>
#include <string>
#include <vector>

#include "AZUREAPI.h"
#include "AZURECalcMCMC.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "ModelPass.h"
#include "PPair.h"
#include "Results.h"
#include "Segment.h"
#include "Unit.h"
#include "WorkingCopy.h"
#include "pybind_engine_call.h"

namespace py = pybind11;

void AddModelData(py::dict &terms, AZUREAPI &api, const vector_r &free);  // pybind_model_data.cpp
std::vector<double> NuclearAmplitudeCoefficients(const Unit &point, const azure2::Results &results);  // pybind_elastic_amplitudes.cpp

namespace {

py::array_t<double> Array(const std::vector<double> &v) {
  return py::array_t<double>(static_cast<py::ssize_t>(v.size()), v.data());
}

// What a point output's reader sees: the data set as set up, the segment, the point and its row, the evaluation's
// results buffer.
struct PointView {
  const WorkingCopy &data;
  const Segment &segment;
  const Unit &point;
  int row;
  const azure2::Results &results;
};

// One number a point carries beside its model value, or a list of them (ragged: a count that varies by point).
struct PointOutput {
  const char *name;
  const char *unit;  // "model": the data set's engine unit
  const char *when;  // the condition on the data set, in words
  const char *doc;
  bool ragged;
  bool (*applies)(const Segment &, const CompoundSystem &);
  std::vector<double> (*read)(const PointView &);
};

bool Always(const Segment &, const CompoundSystem &) { return true; }
// A capture reaction: the exit pair is a photon pair (type 10), or the total capture of the entrance's nucleus.
bool Capture(const Segment &s, const CompoundSystem &nuc) {
  if (s.GetExitKey() == -1) return true;
  const PPair *exit = nuc.GetPairFromKey(s.GetExitKey());
  return exit && exit->GetPType() == 10;
}

// Elastic scattering of a charged pair at an angle, in a sample of one compound nucleus.
bool ChargedElastic(const Segment &s, const CompoundSystem &nuc) {
  if (!s.IsDifferential() || s.IsNatural() || s.GetEntranceKey() != s.GetExitKey()) return false;
  const PPair *pair = nuc.GetPairFromKey(s.GetEntranceKey());
  return pair && pair->GetPType() == 0 && pair->GetZ(1) > 0 && pair->GetZ(2) > 0;
}

const std::vector<PointOutput> &PointOutputs() {
  static const std::vector<PointOutput> rows = {
      {"energy_lab", "MeV", "every data set", "the beam's lab energy of the point (after the data set's energy shift)",
       false, Always, [](const PointView &v) { return std::vector<double>{v.point.GetLabEnergy()}; }},
      {"angle_lab", "deg", "a differential data set", "the lab angle of the point",
       false, [](const Segment &s, const CompoundSystem &) { return s.IsDifferential(); },
       [](const PointView &v) { return std::vector<double>{v.point.GetLabAngle()}; }},
      {"sfactor_conversion", "MeV", "every data set",
       "the point's factor E exp(2 pi eta) from a cross section in b to the S factor in MeV b (0 where it has none)",
       false, Always, [](const PointView &v) { return std::vector<double>{v.point.GetSFactorConversion()}; }},
      {"model_e1", "model", "a capture reaction", "the E1 part of the point's capture cross section (c.m.)",
       false, Capture, [](const PointView &v) { return std::vector<double>{v.point.GetFitE1CrossSection(v.results)}; }},
      {"model_e2", "model", "a capture reaction", "the E2 part of the point's capture cross section (c.m.)",
       false, Capture, [](const PointView &v) { return std::vector<double>{v.point.GetFitE2CrossSection(v.results)}; }},
      {"legendre", "model", "an angular distribution (Legendre coefficients)",
       "the c.m. Legendre coefficients A_L of the point, L = 0, 1, ... (a list per point)",
       true, [](const Segment &s, const CompoundSystem &) { return s.IsAngularDist(); },
       [](const PointView &v) { return v.segment.LegendreCoefficients(v.point, v.results); }},
      {"nuclear_amplitudes", "b/sr", "charged-particle elastic scattering at an angle, a point without a target-effect fold",
       "the coefficients b_l and a_l (re, im) of ENDF-6's nuclear amplitude expansion (LAW=5 LTP=1, ENDF-102 eqs. "
       "6.13-6.14) at the point's energy, as the LIST record orders them (a list per point; empty for a folded point)",
       true, [](const Segment &s, const CompoundSystem &nuc) { return ChargedElastic(s, nuc); },
       [](const PointView &v) { return NuclearAmplitudeCoefficients(v.point, v.results); }},
      {"nucleus_terms", "model", "a sample of several compound nuclei",
       "the abundance-weighted term of each compound nucleus of the sample (a list per point, the nuclei's order)",
       true, [](const Segment &s, const CompoundSystem &) { return s.IsNatural(); },
       [](const PointView &v) { return v.point.GetNaturalTerms(v.results); }},
  };
  return rows;
}

// The arrays of the reported segment at position `first` (0-based; the first of its run of equal keys, as
// AZUREAPI::ReportedSegments) and of the run's rows of chiSquared.out.
py::dict DataSetDict(AZUREAPI &api, const AZUREAPI::ObjectiveBreakdown &b, int first, int last) {
  WorkingCopy &wc = *api.data();
  const CompoundSystem &nuc = *api.compound();
  const azure2::Results &res = api.results();
  const Segment &seg = wc.Data().GetSegments()[first];
  const int index = first + 1;  // the engine's segment number
  double chi2 = 0.0, normChi2 = 0.0, shiftChi2 = 0.0, shiftConstant = 0.0, shiftSqrt = 0.0, shiftFunction = 0.0;
  for (int i = first + 1; i <= last + 1; i++) {
    chi2 += res.objective.segmentChi2[i];
    if (i < (int)b.segments.norm.size()) normChi2 += b.segments.norm[i];
    if (i < (int)b.segments.shift.size()) {
      shiftChi2 += b.segments.shift[i] + b.segments.shiftSqrt[i] + b.segments.shiftFn[i];
      shiftConstant += b.segments.shift[i];
      shiftSqrt += b.segments.shiftSqrt[i];
      shiftFunction += b.segments.shiftFn[i];
    }
  }
  // The data set's rows of the objective's least-squares rows (Objective.residuals): one per point of every segment,
  // in segment order, the data rows first.
  const std::vector<Segment> &all = wc.Data().GetSegments();
  int firstRow = 0, rowCount = 0;
  for (int i = 0; i < first; i++) firstRow += (int)all[i].GetPoints(wc).size();
  for (int i = first; i <= last; i++) rowCount += (int)all[i].GetPoints(wc).size();
  const double norm = wc.NormUsed(index, res);
  std::vector<double> energies, angles, excitation, model, values, errors, residuals;
  const int thmGroup = wc.Data().ThmGroupOf(index);
  if (thmGroup >= 0) wc.ThmGroupResiduals(thmGroup, index, res, residuals);  // the shared profile's rows
  const UnitRange<Unit> points = seg.GetPoints(wc);
  for (int k = 0; k < (int)points.size(); k++) {
    const Unit &p = points[k];
    energies.push_back(p.GetCMEnergy());
    angles.push_back(p.GetCMAngle());
    excitation.push_back(p.GetExcitationEnergy());
    model.push_back(p.GetModelValue(res));
    values.push_back(seg.CMValue(wc.State(seg), k));
    errors.push_back(seg.CMError(wc.State(seg), k));
    if (thmGroup < 0) residuals.push_back(azure2::PointResidual(wc, seg, k, res, norm).Pull());
  }
  py::dict extra;
  for (const PointOutput &o : PointOutputs()) {
    if (!o.applies(seg, nuc)) continue;
    if (o.ragged) {
      py::list perPoint;
      for (int k = 0; k < (int)points.size(); k++) perPoint.append(Array(o.read(PointView{wc, seg, points[k], k, res})));
      extra[o.name] = perPoint;
    } else {
      std::vector<double> column;
      for (int k = 0; k < (int)points.size(); k++) column.push_back(o.read(PointView{wc, seg, points[k], k, res})[0]);
      extra[o.name] = Array(column);
    }
  }
  py::dict d;
  d["segment"] = index;
  d["segment_key"] = seg.GetSegmentKey();
  d["chi2"] = chi2;
  d["norm_chi2"] = normChi2;
  d["shift_chi2"] = shiftChi2;
  d["shift_constant_chi2"] = shiftConstant;
  d["shift_sqrt_chi2"] = shiftSqrt;
  d["shift_function_chi2"] = shiftFunction;
  d["first_row"] = firstRow;
  d["rows"] = rowCount;
  d["norm"] = norm;
  d["energies"] = Array(energies);
  d["angles"] = Array(angles);
  d["excitation_energies"] = Array(excitation);
  d["model"] = Array(model);
  d["values"] = Array(values);
  d["errors"] = Array(errors);
  d["residuals"] = Array(residuals);
  d["extra"] = extra;
  return d;
}

}  // namespace

py::list PointOutputFields() {
  py::list out;
  for (const PointOutput &o : PointOutputs()) {
    py::dict d;
    d["name"] = o.name;
    d["unit"] = o.unit;
    d["when"] = o.when;
    d["doc"] = o.doc;
    d["per_point"] = o.ragged ? "list" : "number";
    out.append(d);
  }
  return out;
}

py::dict EvaluationDict(AZUREAPI &api, const vector_r &free) {
  AZUREAPI::ObjectiveBreakdown b;
  {
    EngineCallScope engine;
    b = api.EvaluateObjectiveRWA(free);
  }
  const WorkingCopy &wc = *api.data();
  const azure2::Results &res = api.results();
  py::dict d;
  d["objective"] = b.total;
  d["term_names"] = api.ObjectiveTermNames();
  d["term_values"] = Array(b.terms);
  d["excluded"] = b.excluded;
  d["log_prob"] = b.excluded ? -std::numeric_limits<double>::infinity() : -0.5 * b.total;
  // The data chi^2 as chiSquared.out's Total-Chi-Squared sums it: the stored row of every segment, in order.
  double chi2 = 0.0;
  for (int i = 1; i <= wc.Data().NumSegments(); i++) chi2 += res.objective.segmentChi2[i];
  d["chi2"] = chi2;
  py::list data;
  const std::vector<Segment> &segments = wc.Data().GetSegments();
  for (int i = 0; i < (int)segments.size();) {
    int last = i;
    while (last + 1 < (int)segments.size() && segments[last + 1].GetSegmentKey() == segments[i].GetSegmentKey()) last++;
    data.append(DataSetDict(api, b, i, last));
    i = last + 1;
  }
  d["data"] = data;
  AddModelData(d, api, free);  // the <modelData> rows: level quantities, no data point is evaluated again
  return d;
}

py::dict ObjectiveRowsDict(AZUREAPI &api, const vector_r &free, bool withJacobian,
                           const std::vector<int> &skip) {
  vector_r rows, jac;
  bool ok;
  {
    EngineCallScope engine;
    ok = api.ObjectiveRowsRWA(free, rows, withJacobian ? &jac : nullptr, skip);
  }
  py::dict d;
  d["ok"] = ok;
  d["rows"] = Array(rows);
  if (withJacobian && ok) {
    py::array_t<double> j({(py::ssize_t)rows.size(), (py::ssize_t)free.size()});
    if (!jac.empty()) std::copy(jac.begin(), jac.end(), j.mutable_data());
    d["jacobian"] = j;
  }
  return d;
}

double LogProbability(AZUREAPI &api, const vector_r &free) {
  EngineCallScope engine;
  return api.LogProbabilityRWA(free);
}

// The engine's ensemble sampler set up on the session (as RunMCMCSampling sets it up at the start x0, fit basis) and
// asked for ln P at x: (ln P, its constant LogPriorNormalisation()).  The tie of pyazr's log_prob to the sampler in
// tests/pyazr/objective_object_test.py.
py::tuple SamplerLogProbability(AZUREAPI &api, const vector_r &x0, const vector_r &x) {
  double lp, constant;
  {
    EngineCallScope engine;
    SamplerStates walkers;
    AZURECalcMCMC sampler(api.data(), api.compound(), api.state(), api.inputs(), walkers, api.configure());
    sampler.PrepareRun(x0, true);
    lp = sampler.LogProbability(x, true);
    constant = sampler.LogPriorNormalisation();
  }
  return py::make_tuple(lp, constant);
}
