// The stopping of the targets from the stopping/ library at set-up (its refusals and its source in the log), and the
// stopping-scale nuisance s (docs/upgrade/energy_grid.md section 6.2, decision D6): one parameter per (material,
// projectile), its hand-over to the effects' sub-points and points, and its prior.
#include <algorithm>
#include <atomic>
#include <cmath>
#include <map>
#include <set>
#include <sstream>

#include "AZUREGrad.h"
#include "Config.h"
#include "WorkingCopy.h"
#include "azure2/experiment/StoppingChecks.h"

void WorkingCopy::SetUpStoppingScaleParams(const Config &configure) {
  Draft().stoppingScaleLeads_.clear();
  Draft().stoppingScaleParamOfEffect_.assign(Draft().NumTargetEffects(), 0);
  std::map<std::string, int> byKey;
  for (int n = 1; n <= Draft().NumTargetEffects(); n++) {
    const TargetEffect *te = Draft().GetTargetEffect(n);
    if (!te->IsTargetIntegration() || !te->HasStoppingScale()) continue;
    auto it = byKey.find(te->GetStoppingScaleKey());
    if (it == byKey.end()) {
      Draft().stoppingScaleLeads_.push_back(n);
      it = byKey.emplace(te->GetStoppingScaleKey(), (int)Draft().stoppingScaleLeads_.size()).first;
      if (te->IsStoppingScaleVaried() && te->IsStoppingScaleSigmaDefault()) {
        std::ostringstream note;
        note << "Note: stopping scale " << te->GetStoppingScaleName() << " is fitted with the prior width "
             << te->GetStoppingScaleSigma() << ", SRIM's stated accuracy for the projectile class (stopping.md section 7); "
             << "write 'scale <s> fit prior <s> <width>' to give another.";
        configure.Note(note.str());
      }
    }
    Draft().stoppingScaleParamOfEffect_[n - 1] = it->second;
  }
}

void WorkingCopy::FillStoppingScalesFromParams(const vector_r &p, const Config &configure) {
  if (Data().NumStoppingScaleParams() == 0 || set_->stoppingScaleParamOffset_ < 0) return;
  std::vector<char> changed(Data().NumTargetEffects(), 0);
  bool any = false;
  for (int k = 1; k <= Data().NumStoppingScaleParams(); k++) {
    const int idx = set_->stoppingScaleParamOffset_ + k - 1;
    if (idx < 0 || idx >= (int)p.size()) return;  // a vector without the stopping-scale block
    const TargetEffect *lead = Data().GetStoppingScaleTargetEffect(k);
    double v = p[idx];
    // The sub-point grids span the layout scale (as the layout density): beyond it the depth integral would run off.
    const double vMax = lead->GetLayoutStoppingScale(), vMin = lead->GetStoppingScaleMinimum();
    if (v > vMax || !(v >= vMin)) {
      static std::atomic<bool> warned{false};
      if (!warned.exchange(true))
        configure.outStream << "WARNING: " << lead->GetStoppingScaleName() << " = " << v << " is outside [" << vMin << ", "
                            << vMax << "], the range its sub-point grids span; the calculation uses the nearest end."
                            << " (Printed once.)" << std::endl;
      v = (v > vMax) ? vMax : vMin;
    }
    for (int n = 1; n <= Data().NumTargetEffects(); n++) {
      if (Data().StoppingScaleParamOfEffect(n) != k) continue;
      TargetEffectValues &values = EffectValues(n);
      if (values.stoppingScale == v) continue;
      values.stoppingScale = v;
      changed[n - 1] = 1;
      any = true;
    }
  }
  if (!any) return;
  // Where s x epsilon enters: the stopping of every sub-point (read by the target kernel and the layer's E_out), and
  // the thickness and the sigma(E_bar) depth of every point, each recomputed from its energy as the set-up made it.
  auto update = [&](const TermRef &seg) {
    if (!seg.segment->IsTargetEffect()) return;
    const int n = seg.GetTargetEffectNum();
    if (n < 1 || n > (int)changed.size() || !changed[n - 1]) return;
    const TargetEffect *te = Data().GetTargetEffect(n);
    const double s = EffectValues(n).stoppingScale, density = EffectValues(n).density;
    for (UnitIterator point = seg.GetPoints(*this).begin(); point < seg.GetPoints(*this).end(); point++) {
      if (point->GetTargetEffectNum() != n) continue;
      const double f = point->GetCMFactor();
      for (int j = 1; j <= point->NumSubPoints(); j++) {
        Unit *sub = point->GetSubPoint(*this, j);
        sub->SetStoppingPower(f * te->StoppingAt(sub->GetCMEnergy() / f, configure, s));
      }
      if (!(point->GetTargetThickness() > 0.0) || te->IsThickTarget()) continue;
      point->SetTargetThickness(f * te->TargetThickness(point->GetTargetLabEnergy(), configure, density, s));
      point->SetTargetFallbackDepth(f * te->MeanEnergyLoss(point->GetTargetLabEnergy(), configure, density, s));
    }
  };
  for (const Segment &seg : set_->segments_) update({&seg, 0});
  for (const LeafRef &leaf : set_->leafRefs_) update({&set_->segments_[leaf.first], leaf.second});
}

double WorkingCopy::StoppingScalePenalty() const {
  double chi2 = 0.0;
  for (int k = 1; k <= Data().NumStoppingScaleParams(); k++)
    chi2 += Data().GetStoppingScaleTargetEffect(k)->StoppingScalePenalty(EffectValues(Data().GetStoppingScaleLeadCopy(k)).stoppingScale);
  return chi2;
}

bool WorkingCopy::CheckTargetStoppings(const Config &configure) {
  std::set<std::string> logged;
  auto refuse = [&](const TargetEffect &te, double eLab) {
    bool dataMissing = false;
    const std::string why = te.GetStopping().Check(eLab, dataMissing);
    if (why.empty()) return false;
    const azure2::experiment::Diagnostic d =
        dataMissing ? azure2::experiment::StoppingDataMissing(te.GetTargetName(), why)
                    : azure2::experiment::StoppingOutOfRange(te.GetTargetName(), te.GetStopping().GetProjectile().label, eLab, why);
    configure.outStream << "ERROR: " << d.ToText() << std::endl;
    return true;
  };
  for (const Segment &segment : Data().GetSegments()) {
    if (!segment.IsTargetEffect()) continue;
    for (int t = 0; t < segment.NumTerms(); t++) {
      const TargetEffect *te = Data().GetTargetEffect(segment.TermTargetEffectNum(t));
      if (!te || !te->IsTargetIntegration() || !te->GetStopping().IsLibrary()) continue;
      for (const Unit &point : segment.TermPoints(*this, t))
        if (refuse(*te, point.GetLabEnergy())) return false;
      // the source of the stopping, once per target and projectile (SRIM is cited in it)
      const std::string label = te->GetTargetName() + "|" + te->GetStopping().GetProjectile().label;
      if (logged.insert(label).second)
        configure.outStream << "Stopping of target " << te->GetTargetName() << " for " << te->GetStopping().GetProjectile().label
                            << ":\n" << te->GetStopping().Provenance() << std::flush;
    }
  }
  return true;
}

bool WorkingCopy::CheckSubPointStoppings(const Config &configure) {
  for (const Segment &segment : Data().GetSegments()) {
    if (!segment.IsTargetEffect()) continue;
    for (int t = 0; t < segment.NumTerms(); t++) {
      const TargetEffect *te = Data().GetTargetEffect(segment.TermTargetEffectNum(t));
      if (!te || !te->IsTargetIntegration() || !te->GetStopping().IsLibrary()) continue;
      for (const Unit &point : segment.TermPoints(*this, t)) {
        if (!point.IsTargetEffect()) continue;
        const double f = point.GetCMFactor();
        for (int j = 1; j <= point.NumSubPoints(); j++) {
          const Unit *sub = point.GetSubPoint(*this, j);
          if (std::isfinite(sub->GetStoppingPower())) continue;
          bool dataMissing = false;
          const double eLab = sub->GetCMEnergy() / f;
          const std::string why = te->GetStopping().Check(eLab, dataMissing);
          configure.outStream << "ERROR: "
                              << azure2::experiment::StoppingOutOfRange(te->GetTargetName(), te->GetStopping().GetProjectile().label,
                                                                        eLab, why.empty() ? "not finite" : why)
                                     .ToText()
                              << std::endl;
          return false;
        }
      }
    }
  }
  return true;
}

void AppendStoppingScaleParams(ParamIndexMap &map, WorkingCopy *data, int &i) {
  if (!data) return;
  for (int k = 1; k <= data->Data().NumStoppingScaleParams(); k++, i++) {
    map.stoppingScaleIndex_[k] = i;
    ParamDesc d{ParamKind::StoppingScale, -1, -1, -1, -1, false};
    d.stoppingScale = k;
    map.desc_.push_back(d);
  }
}
