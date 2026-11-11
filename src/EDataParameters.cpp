// The data-side fit parameters: where each block (normalizations, energy shifts, sqrt(E) and shift-function
// coefficients, sample thicknesses, densities, Doppler temperatures) starts in the parameter vector, their names
// and start values (FillMnParams), the shift, shift-function, normalization and thickness groups, the density
// parameters of the target effects, and the prior centres given by a project's specs.
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "AZUREParams.h"
#include "Config.h"
#include "EData.h"
#include "GSLException.h"
#include "Minuit2/MnUserParameters.h"

/*!
 * Moves the nominal normalization / energy shift of the segments whose specs give a prior centre (after the groups
 * are resolved, so that a group's follower keeps its own centre).  A centre for a data set that has no segment (not
 * active, or skipped at load) is reported and ignored.  A refusal of the reader (specs.priorCentreError) is reported
 * after the centres before it are applied, and fails the run.  Segment keys are unique in the data set, so at most one
 * segment matches a spec.
 */

int EData::ApplyPriorCentres(const Config &configure, const azure2::experiment::DataSpecs &specs) {
  for (const auto &item : specs.datasets.items) {
    const azure2::experiment::DatasetSpec &d = item.spec;
    if (!d.hasNormCentre && !d.hasShiftCentre) continue;
    ESegment *segment = nullptr;
    for (int s = 1; s <= Draft().NumSegments(); s++)
      if (Draft().GetSegment(s)->GetSegmentKey() == d.key) segment = Draft().GetSegment(s);
    if (d.hasNormCentre) {
      if (segment)
        segment->SetNominalNorm(d.normCentre);
      else
        configure.outStream << "WARNING: <parameterSettings> row \"" << d.normCentreSource
                            << "\" names no active data segment; ignored." << std::endl;
    }
    if (d.hasShiftCentre) {
      if (segment)
        segment->SetNominalEnergyShift(d.shiftCentre);
      else
        configure.outStream << "WARNING: <parameterSettings> row \"" << d.shiftCentreSource
                            << "\" names no active data segment; ignored." << std::endl;
    }
  }
  if (!specs.priorCentreError.empty()) {
    configure.outStream << specs.priorCentreError << std::endl;
    return -1;
  }
  return 0;
}

/*!
 * Links the segments of each shift group ("shiftgroup <label>" on the segment
 * line) to the group's first segment in file order: that segment keeps its own
 * constant-shift parameter (value, uncertainty, vary flag, penalty), and every
 * other member takes that parameter's value instead of its own, which is fixed
 * in the fit.  A note per group lists the members.
 */

int EData::ResolveShiftGroups(const Config &configure) {
  std::vector<std::pair<std::string, int>> leaders;  // label, leader position
  for (int pos = 0; pos < (int)Draft().segments_.size(); pos++) {
    ESegment &seg = Draft().segments_[pos];
    const std::string label = seg.GetShiftGroup();
    if (label.empty()) continue;
    auto it = std::find_if(leaders.begin(), leaders.end(),
                           [&label](const std::pair<std::string, int> &e) { return e.first == label; });
    if (it == leaders.end()) {
      leaders.push_back({label, pos});
      continue;
    }
    ESegment &leader = Draft().segments_[it->second];
    seg.SetShiftGroupFollower(it->second, State(leader).energyShift);
    State(seg).energyShift = State(leader).energyShift;
  }
  for (const auto &g : leaders) {
    ESegment &leader = Draft().segments_[g.second];
    std::ostringstream members;
    int n = 0;
    for (ESegment &seg : Draft().segments_)
      if (seg.GetShiftLeaderPos() == g.second) {
        members << " " << seg.GetSegmentKey();
        n++;
      }
    if (n == 0) {
      configure.Note("Note: shift group '" + g.first + "' has only segment #" +
                     std::to_string(leader.GetSegmentKey()) + "; nothing is shared.");
      continue;
    }
    configure.Note("Note: shift group '" + g.first + "': segments #" + std::to_string(leader.GetSegmentKey()) + " +" +
                   members.str() + " share one constant energy shift (segment #" +
                   std::to_string(leader.GetSegmentKey()) + "'s parameter, " +
                   (leader.IsVaryEnergyShift() ? "free" : "fixed") + ").");
  }
  return ResolveShiftFunctions(configure);
}

namespace {

// Indices j of the parameters "aj" an Equation expression refers to.  The Equation
// class reads "a" followed by digits as a parameter (asin, acos, atan are not).
std::set<int> ShiftFnParamRefs(const std::string &expression) {
  std::set<int> refs;
  for (size_t k = 0; k + 1 < expression.size(); k++) {
    if (expression[k] != 'a' || !std::isdigit((unsigned char)expression[k + 1])) continue;
    size_t e = k + 1;
    while (e < expression.size() && std::isdigit((unsigned char)expression[e])) e++;
    refs.insert(std::atoi(expression.substr(k + 1, e - k - 1).c_str()));
    k = e - 1;
  }
  return refs;
}

}  // namespace

/*!
 * Shift functions ("shiftfn" block, design note 10-13-26): the segment declaring one
 * owns its parameters; it and every follower of its shift group apply
 * dE(x) = a + f(x; a0..), x the point's original lab energy.  Rejected: a follower
 * declaring a function (the group's calibration belongs to its first segment), a
 * function together with a sqrt(E) term on any member (two descriptions of one
 * calibration), and an expression referring to an undeclared parameter or leaving a
 * declared one unused.  A shift-function group ("shiftfngroup <label>") also applies
 * the function of its first segment to every member, but each member keeps its own
 * constant shift; such a group may not overlap a shift group.  For each function the
 * free shift terms (each free constant a of its members, and every varied aj) are
 * checked for degeneracy over the members' data energies; a term with no effect there
 * is an error.  One log line per function.
 */

int EData::ResolveShiftFunctions(const Config &configure) {
  // Shift-function groups: the first segment with a label (file order) declares the
  // function; the others apply it and keep their own constant shifts.
  std::vector<int> fnGroupLeader(Draft().segments_.size(), -1);
  {
    std::vector<std::pair<std::string, int>> fnLeaders;  // label, leader position
    for (int pos = 0; pos < (int)Draft().segments_.size(); pos++) {
      ESegment &seg = Draft().segments_[pos];
      const std::string label = seg.GetShiftFnGroup();
      if (label.empty()) continue;
      if (!seg.GetShiftGroup().empty()) {
        configure.outStream << "ERROR: Data segment #" << seg.GetSegmentKey() << " is in shift group '"
                            << seg.GetShiftGroup() << "' and in shift-function group '" << label
                            << "': use one of them (shiftgroup shares the constant shift as well)." << std::endl;
        return -1;
      }
      auto it = std::find_if(fnLeaders.begin(), fnLeaders.end(),
                             [&label](const std::pair<std::string, int> &e) { return e.first == label; });
      if (it == fnLeaders.end()) {
        fnLeaders.push_back({label, pos});
        if (!seg.OwnsShiftFn()) {
          configure.outStream << "ERROR: Data segment #" << seg.GetSegmentKey() << " is the first segment of "
                              << "shift-function group '" << label << "' but declares no shiftfn: put the "
                              << "group's shiftfn block on it." << std::endl;
          return -1;
        }
        continue;
      }
      if (seg.OwnsShiftFn()) {
        configure.outStream << "ERROR: Data segment #" << seg.GetSegmentKey() << " declares a shiftfn, but it is not "
                            << "the first segment of shift-function group '" << label << "' (segment #"
                            << Draft().segments_[it->second].GetSegmentKey() << " is): one function per group." << std::endl;
        return -1;
      }
      fnGroupLeader[pos] = it->second;
    }
  }
  int fnStart = 0;
  for (int pos = 0; pos < (int)Draft().segments_.size(); pos++) {
    ESegment &owner = Draft().segments_[pos];
    if (!owner.OwnsShiftFn()) continue;
    const int key = owner.GetSegmentKey();
    const std::string expression = owner.GetDeclaredShiftFn();
    const int n = (int)owner.GetNominalShiftFnValues().size();
    if (owner.IsShiftGroupFollower()) {
      configure.outStream << "ERROR: Data segment #" << key << " declares a shiftfn, but it is not the first segment of "
                          << "shift group '" << owner.GetShiftGroup() << "' (segment #"
                          << Draft().segments_[owner.GetShiftLeaderPos()].GetSegmentKey()
                          << " is): put the block on that segment; it applies to the whole group." << std::endl;
      return -1;
    }
    std::vector<int> members{pos};
    for (int f = 0; f < (int)Draft().segments_.size(); f++)
      if (Draft().segments_[f].GetShiftLeaderPos() == pos || fnGroupLeader[f] == pos) members.push_back(f);
    for (int m : members)
      if (Draft().segments_[m].HasSqrtShiftBlock()) {
        configure.outStream << "ERROR: Data segment #" << Draft().segments_[m].GetSegmentKey() << " has a sqrtshift block and "
                            << "shift function \"" << expression << "\" (segment #" << key
                            << "): two descriptions of one energy calibration.  Write the sqrt(E) term into the function, "
                            << "e.g. a" << n << "*sqrt(x), and drop the sqrtshift block." << std::endl;
        return -1;
      }
    const std::set<int> refs = ShiftFnParamRefs(expression);
    for (int j : refs)
      if (j >= n) {
        configure.outStream << "ERROR: Data segment #" << key << ": shift function \"" << expression << "\" uses a" << j
                            << " but declares " << n << " parameter" << (n == 1 ? "" : "s") << " (a0..a" << n - 1 << ")."
                            << std::endl;
        return -1;
      }
    for (int j = 0; j < n; j++)
      if (!refs.count(j)) {
        configure.outStream << "ERROR: Data segment #" << key << ": shift function \"" << expression << "\" declares a"
                            << j << " but never uses it." << std::endl;
        return -1;
      }
    for (int m : members)
      Draft().segments_[m].InstallShiftFn(State(Draft().segments_[m]), owner, pos, configure);
    owner.SetShiftFnParamStart(fnStart);
    fnStart += n;

    // Degeneracy check: rows = the members' distinct original lab energies (per member,
    // since members of a shift-function group have their own constants); columns =
    // d(dE)/d(term) for every free constant (one per constant owner) and every varied aj.
    std::vector<std::pair<int, double>> rows;  // constant owner position, energy
    for (int m : members) {
      const int constOwner = Draft().segments_[m].GetShiftLeaderPos() >= 0 ? Draft().segments_[m].GetShiftLeaderPos() : m;
      std::vector<double> e;
      for (int pid = 1; pid <= Draft().segments_[m].NumPoints(); pid++) {
        EPoint *pt = Draft().segments_[m].GetPoint(*this, pid);
        e.push_back(pt->GetOriginalEnergy() > 0 ? pt->GetOriginalEnergy() : pt->GetLabEnergy());
      }
      std::sort(e.begin(), e.end());
      e.erase(std::unique(e.begin(), e.end()), e.end());
      for (double x : e) rows.push_back({constOwner, x});
    }
    std::vector<double> energies;
    for (const auto &r : rows) energies.push_back(r.second);
    std::sort(energies.begin(), energies.end());
    energies.erase(std::unique(energies.begin(), energies.end()), energies.end());
    std::vector<std::string> names;
    std::vector<std::vector<double>> cols;
    std::vector<int> constOwners;
    for (const auto &r : rows)
      if (std::find(constOwners.begin(), constOwners.end(), r.first) == constOwners.end()) constOwners.push_back(r.first);
    for (int c : constOwners) {
      if (!Draft().segments_[c].IsVaryEnergyShift()) continue;
      names.push_back(constOwners.size() == 1 ? std::string("a") : "a#" + std::to_string(Draft().segments_[c].GetSegmentKey()));
      std::vector<double> col(rows.size());
      for (size_t r = 0; r < rows.size(); r++) col[r] = (rows[r].first == c) ? 1.0 : 0.0;
      cols.push_back(col);
    }
    SegmentState probe = State(owner);  // a scratch copy: the owner's own values stay untouched
    for (int j = 0; j < n; j++) {
      if (!owner.IsVaryShiftFn(j)) continue;
      const double v0 = probe.shiftFnValues[j], h = 1.0e-6 * (std::fabs(v0) + 1.0);
      std::vector<double> c(rows.size());
      for (size_t r = 0; r < rows.size(); r++) {
        owner.SetShiftFnValue(probe, j, v0 + h);
        const double fp = owner.ShiftFnTerm(probe, rows[r].second);
        owner.SetShiftFnValue(probe, j, v0 - h);
        c[r] = (fp - owner.ShiftFnTerm(probe, rows[r].second)) / (2.0 * h);
      }
      owner.SetShiftFnValue(probe, j, v0);
      double norm = 0.0;
      for (double v : c) norm += v * v;
      if (!(norm > 0.0) || !std::isfinite(norm)) {
        configure.outStream << "ERROR: Data segment #" << key << ": in shift function \"" << expression << "\", a" << j
                            << " has no effect at this group's energies (" << energies.front() << "-" << energies.back()
                            << " MeV); check the expression." << std::endl;
        return -1;
      }
      names.push_back("a" + std::to_string(j));
      cols.push_back(c);
    }
    std::ostringstream note;
    note << "Note: shift function on segment #" << key;
    if (members.size() > 1) {
      const bool fnGroup = !owner.GetShiftFnGroup().empty();
      note << " (" << (fnGroup ? "shift-function group '" + owner.GetShiftFnGroup() : "group '" + owner.GetShiftGroup())
           << "', segments #" << key;
      for (size_t m = 1; m < members.size(); m++) note << " " << Draft().segments_[members[m]].GetSegmentKey();
      note << (fnGroup ? "; each keeps its own constant shift)" : ")");
    }
    note << ": dE = a + " << expression << ", x = lab energy (MeV);";
    for (int j = 0; j < n; j++) {
      const double err = owner.GetShiftFnErrors()[j];
      note << " a" << j << " " << owner.GetNominalShiftFnValues()[j]
           << (owner.IsVaryShiftFn(j) ? (err != 0.0 ? " free (penalty width " + std::to_string(err) + ")" : " free, unpenalised")
                                      : " fixed")
           << (j + 1 < n ? "," : ".");
    }
    const int m = (int)rows.size(), k = (int)cols.size();
    if (k >= 2) {
      // Condition number of the column-normalised matrix (scale-free).
      double cond = std::numeric_limits<double>::infinity();
      if (m >= k) {
        gsl_matrix *A = gsl_matrix_alloc(m, k);
        gsl_matrix *V = gsl_matrix_alloc(k, k);
        gsl_vector *S = gsl_vector_alloc(k), *work = gsl_vector_alloc(k);
        for (int c = 0; c < k; c++) {
          double norm = 0.0;
          for (int r = 0; r < m; r++) norm += cols[c][r] * cols[c][r];
          norm = std::sqrt(norm);
          for (int r = 0; r < m; r++) gsl_matrix_set(A, r, c, cols[c][r] / norm);
        }
        int svd;
        {
          GslQuiet quiet;  // scoped: a process-wide handler swap would outlive this check
          svd = gsl_linalg_SV_decomp(A, V, S, work);
        }
        if (svd == 0 && gsl_vector_get(S, k - 1) > 0.0)
          cond = gsl_vector_get(S, 0) / gsl_vector_get(S, k - 1);
        gsl_matrix_free(A);
        gsl_matrix_free(V);
        gsl_vector_free(S);
        gsl_vector_free(work);
      }
      note << " Condition number of the free terms (";
      for (int c = 0; c < k; c++) note << names[c] << (c + 1 < k ? ", " : "");
      note << ") over " << m << " member energ" << (m == 1 ? "y" : "ies") << " " << energies.front() << "-" << energies.back()
           << " MeV: " << cond << ".";
      configure.Note(note.str());
      if (!(cond < 1.0e8))
        configure.outStream << "WARNING: shift function on segment #" << key << ": its free terms cannot be told apart "
                            << "over this group's energies (condition number " << cond << "); fix one of them or "
                            << "constrain it with a penalty." << std::endl;
    } else {
      configure.Note(note.str());
    }
  }
  return 0;
}

/*!
 * Writes each shift-group follower's parameter as the value it actually uses,
 * its group leader's, so that param.sav / param.fit read back unchanged.
 */

void EData::SyncShiftGroupParams(ROOT::Minuit2::MnUserParameters &p) const {
  const int off = set_->energyShiftParamOffset_;
  if (off < 0) return;
  for (int pos = 0; pos < (int)set_->segments_.size(); pos++) {
    const int lead = set_->segments_[pos].GetShiftLeaderPos();
    if (lead < 0 || off + pos >= (int)p.Params().size() || off + lead >= (int)p.Params().size()) continue;
    p.SetValue(off + pos, p.Value(off + lead));
  }
}

/*!
 * Links the segments of each normalization group ("normgroup <label>" on the
 * segment line) to the group's first segment in file order: that segment keeps
 * its own norm (value, uncertainty, vary flag, penalty), and every other member
 * takes its value instead.  A follower whose line varies its norm keeps that
 * parameter, fixed, so the parameter layout does not depend on the groups.  A
 * note per group lists the members.
 */

void EData::ResolveNormGroups(const Config &configure) {
  std::vector<std::pair<std::string, int>> leaders;  // label, leader position
  for (int pos = 0; pos < (int)Draft().segments_.size(); pos++) {
    ESegment &seg = Draft().segments_[pos];
    const std::string label = seg.GetNormGroup();
    if (label.empty()) continue;
    auto it = std::find_if(leaders.begin(), leaders.end(),
                           [&label](const std::pair<std::string, int> &e) { return e.first == label; });
    if (it == leaders.end()) {
      leaders.push_back({label, pos});
      continue;
    }
    ESegment &leader = Draft().segments_[it->second];
    seg.SetNormGroupFollower(it->second, State(leader).norm, leader.GetNominalNorm());
    State(seg).norm = State(leader).norm;
  }
  for (const auto &g : leaders) {
    ESegment &leader = Draft().segments_[g.second];
    std::ostringstream members;
    int n = 0;
    for (ESegment &seg : Draft().segments_)
      if (seg.GetNormLeaderPos() == g.second) {
        members << " " << seg.GetSegmentKey();
        n++;
      }
    if (n == 0) {
      configure.Note("Note: normalization group '" + g.first + "' has only segment #" +
                     std::to_string(leader.GetSegmentKey()) + "; nothing is shared.");
      continue;
    }
    configure.Note("Note: normalization group '" + g.first + "': segments #" + std::to_string(leader.GetSegmentKey()) +
                   " +" + members.str() + " share one normalization (segment #" +
                   std::to_string(leader.GetSegmentKey()) + "'s, " + (leader.IsVaryNorm() ? "free" : "fixed") + ").");
  }
}

/*!
 * Writes each norm-group follower's (fixed) norm parameter as the value it
 * actually uses, its group leader's, so that param.sav / param.fit read back
 * unchanged.
 */

void EData::SyncNormGroupParams(ROOT::Minuit2::MnUserParameters &p) const {
  const int off = Data().GetNormParamOffset();
  if (off < 0) return;
  std::vector<int> index(set_->segments_.size(), -1);  // norm parameter of each segment, -1 = none
  int i = off;
  for (int pos = 0; pos < (int)set_->segments_.size(); pos++)
    if (set_->segments_[pos].IsVaryNorm()) index[pos] = i++;
  for (int pos = 0; pos < (int)set_->segments_.size(); pos++) {
    const int lead = set_->segments_[pos].GetNormLeaderPos();
    if (lead < 0 || index[pos] < 0 || index[pos] >= (int)p.Params().size()) continue;
    const double value = (index[lead] >= 0 && index[lead] < (int)p.Params().size()) ? p.Value(index[lead])
                                                                                     : State(Data().segments_[lead]).norm;
    p.SetValue(index[pos], value);
  }
  // Thickness-group followers likewise: their (fixed) segment_<key>_thickness is written as the leader's value.
  if (set_->thicknessParamOffset_ >= 0) {
    for (int k = 1; k <= Data().NumThicknessParams(); k++) {
      const ESegment &seg = set_->segments_[set_->thicknessSegments_[k - 1] - 1];
      const int idx = set_->thicknessParamOffset_ + k - 1;
      if (!seg.IsThicknessGroupFollower() || idx >= (int)p.Params().size()) continue;
      const int lk = Data().ThicknessParamOfSegment(set_->thicknessSegments_[k - 1]);
      const int lidx = (lk > 0) ? set_->thicknessParamOffset_ + lk - 1 : -1;
      p.SetValue(idx, (lidx >= 0 && lidx < (int)p.Params().size()) ? p.Value(lidx)
                                                                     : Data().segments_[seg.GetThicknessLeaderPos()].GetSampleThickness(State(Data().segments_[seg.GetThicknessLeaderPos()])));
    }
  }
}

/*!
 * Fills the Minuit parameter array from initial values in the EData object.
 */

namespace {

// The name of a target effect's parameter: from its target in format 3 (<target>.areal_density), else the old
// <kind>_<effect number>.
std::string EffectParameterName(const TargetEffect &te, const char *oldKind, const char *field, int effect) {
  return te.GetTargetName().empty() ? std::string(oldKind) + "_" + std::to_string(effect) : te.GetTargetName() + field;
}

}  // namespace

void EData::FillMnParams(AZUREParams &params) {
  ROOT::Minuit2::MnUserParameters &p = params.GetMinuitParams();
  // 1-based position of a segment in the list (the records' and ParamIndexMap's numbering).
  auto SegmentPosition = [this](ESegmentIterator it) { return (int)(it - Data().GetSegments().begin()) + 1; };
  // The layout is the data set's; laid out again by a later call (pyazr) only if it moved (ReviseDataSet).
  struct Layout {
    int norm, shift, shiftSqrt, shiftFn, screening, density, doppler, thickness, thmCoherent;
  } layout;
  layout.norm = (int)p.Params().size();
  for (ESegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++) {
    if (segment->IsVaryNorm() && !segment->IsProfiledNorm()) {
      const std::string &varname = segment->GetParameterNames().norm;
      const double norm = State(*segment).norm;
      params.Add(varname, norm, norm * 0.05,
                 ParameterRecord::OfSegment(ParamKind::Norm, SegmentPosition(segment), segment->GetSegmentKey()));
      p.SetLowerLimit(varname, 0.0);
      if (segment->IsNormGroupFollower()) p.Fix(varname);  // the group's norm is its leader's parameter

      // Parameter settings will be applied by ParameterLimitsManager during fit
    }
    // With component segments, no need to skip - total capture is handled within the segment
  }

  // Add energy shift parameters FOR ALL SEGMENTS (like normalization)
  layout.shift = (int)p.Params().size();
  for (ESegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++) {
    // Always add energy shift parameter, regardless of IsVaryEnergyShift()
    const std::string &varname = segment->GetParameterNames().energyShift;
    double stepSize = (segment->GetEnergyShiftError() > 0.0) ? segment->GetEnergyShiftError() * 0.01 : 0.0005;
    params.Add(varname, State(*segment).energyShift, stepSize,
               ParameterRecord::OfSegment(ParamKind::EnergyShift, SegmentPosition(segment), segment->GetSegmentKey()));

    if (!segment->IsVaryEnergyShift()) {
      p.Fix(varname);  // Fix parameter if not varying
    }
  }

  // sqrt(E) energy-shift coefficients, again one per segment so the layout
  // (and every parameter index downstream) depends only on the segment count.
  // A segment without a "sqrtshift" block carries a fixed 0 here, so a model
  // written before the term existed fits exactly as before.
  layout.shiftSqrt = (int)p.Params().size();
  for (ESegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++) {
    const std::string &varname = segment->GetParameterNames().energyShiftSqrt;
    double stepSize = (segment->GetEnergyShiftSqrtError() > 0.0) ? segment->GetEnergyShiftSqrtError() * 0.01 : 0.0005;
    params.Add(varname, State(*segment).energyShiftSqrt, stepSize,
               ParameterRecord::OfSegment(ParamKind::EnergyShiftSqrt, SegmentPosition(segment), segment->GetSegmentKey()));
    if (!segment->IsVaryEnergyShiftSqrt()) p.Fix(varname);
  }

  // Shift-function parameters: only segments that declare a "shiftfn" block have
  // any, so a model without one gets no new parameter at all.
  layout.shiftFn = (int)p.Params().size();
  for (ESegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++) {
    if (!segment->OwnsShiftFn() || !segment->HasShiftFn()) continue;
    for (int j = 0; j < segment->NumShiftFnParams(); j++) {
      const std::string &varname = segment->GetParameterNames().shiftFunction.at(j);
      const double err = segment->GetShiftFnErrors()[j];
      params.Add(varname, State(*segment).shiftFnValues[j], (err > 0.0) ? err * 0.01 : 0.0005,
                 ParameterRecord::ShiftFnTerm(SegmentPosition(segment), segment->GetSegmentKey(), j));
      if (!segment->IsVaryShiftFn(j)) p.Fix(varname);
    }
  }

  // Electron-screening potentials (eV): one per target effect with a "screening" block, none otherwise.
  layout.screening = (int)p.Params().size();
  for (int k = 1; k <= Data().NumScreeningParams(); k++) {
    const TargetEffect *te = Data().GetScreeningTargetEffect(k);
    const std::string varname = EffectParameterName(*te, "screening", ".screening", Data().GetScreeningEffect(k));
    const double ue = EffectValues(Data().ScreeningCopies(k).front()).screeningUe;
    params.Add(varname, ue, std::max(10.0, 0.05 * std::fabs(ue)),
               ParameterRecord::OfEffect(ParamKind::Screening, Data().GetScreeningEffect(k)));
    p.SetLowerLimit(varname, 0.0);
    if (!te->IsScreeningVaried()) p.Fix(varname);
  }

  // Fitted target densities (atoms/cm^2): one per densityfit line or density group, none otherwise.  Bounded by
  // the density the sub-point grids span (TargetEffect::GetLayoutDensity) and a hundredth of the line's.
  layout.density = (int)p.Params().size();
  for (int k = 1; k <= Data().NumDensityParams(); k++) {
    const TargetEffect *te = Data().GetDensityTargetEffect(k);
    const std::string varname = EffectParameterName(*te, "density", ".areal_density", Data().GetDensityEffect(k));
    const double n = EffectValues(Data().GetDensityLeadCopy(k)).density;
    params.Add(varname, n, (te->GetDensitySigma() > 0.0) ? 0.1 * te->GetDensitySigma() : 0.01 * n,
               ParameterRecord::OfEffect(ParamKind::Density, Data().GetDensityEffect(k)));
    p.SetLimits(varname, 0.01 * te->GetDensityNominal(), te->GetLayoutDensity());
    if (!te->IsDensityVaried()) p.Fix(varname);
  }

  // Fitted Doppler temperatures (K): one per line with the fit form 'doppler T vary sigma [Tmax]', none otherwise.
  // Bounded by the temperature the grids span (TargetEffect::GetDopplerLayoutTemperature) and T0/4.
  layout.doppler = (int)p.Params().size();
  for (int k = 1; k <= Data().NumDopplerParams(); k++) {
    const TargetEffect *te = Data().GetDopplerTargetEffect(k);
    const std::string varname = EffectParameterName(*te, "doppler", ".doppler", Data().GetDopplerEffect(k));
    const double t = EffectValues(Data().DopplerCopies(k).front()).dopplerTemperature;
    params.Add(varname, t, (te->GetDopplerSigma() > 0.0) ? 0.1 * te->GetDopplerSigma() : 0.01 * t,
               ParameterRecord::OfEffect(ParamKind::Doppler, Data().GetDopplerEffect(k)));
    p.SetLimits(varname, te->GetDopplerMinimum(), te->GetDopplerLayoutTemperature());
    if (!te->IsDopplerVaried()) p.Fix(varname);
  }

  // Fitted sample thicknesses (atoms/barn): one per data segment whose transmission / yield block has the fit form
  // 'n vary sigma', none otherwise.  A thickness-group follower's is fixed (its value is its leader's).
  layout.thickness = (int)p.Params().size();
  for (int k = 1; k <= Data().NumThicknessParams(); k++) {
    const ESegment *seg = Data().GetSegment(Data().GetThicknessSegment(k));
    const std::string &varname = seg->GetParameterNames().thickness;
    const double n = seg->GetSampleThickness(State(*seg));
    params.Add(varname, n, (seg->GetThicknessSigma() > 0.0) ? 0.1 * seg->GetThicknessSigma() : 0.01 * n,
               ParameterRecord::OfSegment(ParamKind::Thickness, Data().GetThicknessSegment(k), seg->GetSegmentKey()));
    p.SetLowerLimit(varname, 0.0);
    if (!seg->IsThicknessVaried()) p.Fix(varname);
  }

  // THM coherent backgrounds (cbackground=): Re/Im of each amplitude, last.
  layout.thmCoherent = Data().thmCoherentParams_.empty() ? -1 : (int)p.Params().size();
  for (int k = 0; k < Data().NumThmCoherentParams(); k++) {
    const ThmCoherentParam &c = Data().GetThmCoherentParam(k);
    const double value = thmCoherentValues_[k];
    params.Add(c.name, value, value != 0.0 ? 0.1 * std::fabs(value) : 0.1,
               ParameterRecord::ThmCoherent(c.experiment, c.jGroup, c.entrance, c.exit, c.term));
    if (c.fixed) p.Fix(c.name);
  }
  const DataSet &d = Data();
  if (d.normParamOffset_ != layout.norm || d.energyShiftParamOffset_ != layout.shift ||
      d.energyShiftSqrtParamOffset_ != layout.shiftSqrt || d.energyShiftFnParamOffset_ != layout.shiftFn ||
      d.screeningParamOffset_ != layout.screening || d.densityParamOffset_ != layout.density ||
      d.dopplerParamOffset_ != layout.doppler || d.thicknessParamOffset_ != layout.thickness ||
      d.thmCoherentParamOffset_ != layout.thmCoherent)
    ReviseDataSet([&layout](DataSet &next) {
      next.normParamOffset_ = layout.norm;
      next.energyShiftParamOffset_ = layout.shift;
      next.energyShiftSqrtParamOffset_ = layout.shiftSqrt;
      next.energyShiftFnParamOffset_ = layout.shiftFn;
      next.screeningParamOffset_ = layout.screening;
      next.densityParamOffset_ = layout.density;
      next.dopplerParamOffset_ = layout.doppler;
      next.thicknessParamOffset_ = layout.thickness;
      next.thmCoherentParamOffset_ = layout.thmCoherent;
    });
}

int EData::ResolveThicknessGroups(const Config &configure) {
  std::vector<std::pair<std::string, int>> leaders;  // label, leader position (0-based)
  for (int pos = 0; pos < (int)Draft().segments_.size(); pos++) {
    ESegment &seg = Draft().segments_[pos];
    const std::string label = seg.GetThicknessGroup();
    if (label.empty()) continue;
    if (!seg.IsTransformed()) {
      configure.outStream << "ERROR: segment #" << seg.GetSegmentKey() << ": thicknessgroup '" << label
                          << "' needs the segment's own transmission or yield block (the group shares the sample "
                          << "thickness of transmission / yield segments)." << std::endl;
      return -1;
    }
    auto it = std::find_if(leaders.begin(), leaders.end(),
                           [&label](const std::pair<std::string, int> &e) { return e.first == label; });
    if (it == leaders.end()) {
      leaders.push_back({label, pos});
      continue;
    }
    const ESegment &leader = Draft().segments_[it->second];
    seg.SetThicknessGroupFollower(it->second, leader.GetSampleThickness(State(leader)));
    State(seg).thickness = leader.GetSampleThickness(State(leader));
  }
  for (const auto &g : leaders) {
    ESegment &leader = Draft().segments_[g.second];
    std::ostringstream members;
    int n = 0;
    for (ESegment &seg : Draft().segments_)
      if (seg.GetThicknessLeaderPos() == g.second) {
        members << " " << seg.GetSegmentKey();
        n++;
      }
    if (n == 0) {
      configure.Note("Note: thickness group '" + g.first + "' has only segment #" +
                     std::to_string(leader.GetSegmentKey()) + "; nothing is shared.");
      continue;
    }
    configure.Note("Note: thickness group '" + g.first + "': segments #" + std::to_string(leader.GetSegmentKey()) +
                   " +" + members.str() + " share one sample thickness (segment #" +
                   std::to_string(leader.GetSegmentKey()) + "'s, " +
                   (leader.HasThicknessFit() ? "parameter " + leader.GetParameterNames().thickness
                                             : std::string("fixed")) +
                   ").");
  }
  Draft().thicknessSegments_.clear();
  for (int pos = 0; pos < (int)Draft().segments_.size(); pos++)
    if (Draft().segments_[pos].HasThicknessFit()) Draft().thicknessSegments_.push_back(pos + 1);
  return 0;
}

/*!
 * Groups the densityfit / densitygroup lines into parameters.  SetUpTargetEffects gives every segment a line
 * lists its own copy of the effect (and every component of such a segment one more), so the effects are visited
 * by line (TargetEffect::GetLineNumber): the first copy of a line decides, the other copies follow it.  A
 * line that lists no existing segment has no copy and so no parameter (dev makes one that moves nothing).
 */
int EData::SetUpDensityParams(const Config &configure) {
  Draft().densityEffects_.clear();
  Draft().densityLeadCopies_.clear();
  Draft().densityParamOfEffect_.assign(Draft().NumTargetEffects(), 0);
  std::map<int, int> firstCopyOfLine;       // line number -> its first copy (target effect number)
  std::map<std::string, int> groupLeader;   // label -> first copy of the group's first line
  for (int n = 1; n <= Draft().NumTargetEffects(); n++) {
    TargetEffect *te = Draft().GetTargetEffect(n);
    const int line = te->GetLineNumber();
    auto seen = firstCopyOfLine.find(line);
    if (seen != firstCopyOfLine.end()) {
      // Another copy of a line already grouped: the same parameter, and the group's values if it joined one.
      const TargetEffect *lineFirst = Draft().GetTargetEffect(seen->second);
      if (lineFirst->HasDensityFit()) te->JoinDensityGroup(*lineFirst);
      Draft().densityParamOfEffect_[n - 1] = Draft().densityParamOfEffect_[seen->second - 1];
      continue;
    }
    firstCopyOfLine[line] = n;
    const std::string &group = te->GetDensityGroup();
    if (!group.empty()) {
      auto it = groupLeader.find(group);
      if (it != groupLeader.end()) {
        const TargetEffect *lead = Draft().GetTargetEffect(it->second);
        const int leadLine = lead->GetLineNumber();
        if (std::fabs(te->GetDensity() - lead->GetDensity()) > 1.0e-9 * lead->GetDensity()) {
          std::ostringstream note;
          note << "Note: <targetInt> effect " << line << " (density group " << group << ") has density "
               << te->GetDensity() << "; the group's density is its first line's, " << lead->GetDensity()
               << " (effect " << leadLine << ").";
          configure.Note(note.str());
        }
        if (te->HasDensityFit() && (te->IsDensityVaried() != lead->IsDensityVaried() ||
                                    std::fabs(te->GetDensitySigma() - lead->GetDensitySigma()) > 1.0e-9 * lead->GetDensity()))
          configure.Note("Note: <targetInt> effect " + std::to_string(line) + " (density group " + group +
                         "): its densityfit block is ignored, the group's is its first line's (effect " +
                         std::to_string(leadLine) + ").");
        te->JoinDensityGroup(*lead);
        Draft().densityParamOfEffect_[n - 1] = Draft().densityParamOfEffect_[it->second - 1];
        continue;
      }
      if (!te->HasDensityFit()) {
        configure.outStream << "ERROR: density group " << group << ": its first <targetInt> line (effect " << line
                            << ") has no 'densityfit vary sigma' block; the group's parameter and prior are read from it."
                            << std::endl;
        return -1;
      }
      groupLeader[group] = n;
    }
    if (!te->HasDensityFit()) continue;
    Draft().densityEffects_.push_back(line);
    Draft().densityLeadCopies_.push_back(n);
    Draft().densityParamOfEffect_[n - 1] = (int)Draft().densityEffects_.size();
  }
  return 0;
}
