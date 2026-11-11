// The stopping of the targets from the stopping/ library at set-up (its refusals and its source in the log), and the
// stopping scale s (docs/upgrade/energy_grid.md section 6.2, decision D6): its row of the apparatus parameters
// (ApparatusParameters.h), one parameter per (material, projectile), and what follows a change of s.
#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>

#include "ApparatusParameters.h"
#include "Config.h"
#include "WorkingCopy.h"
#include "azure2/experiment/StoppingChecks.h"

namespace {

// One declaration per target-integration effect with a scale, keyed by (material, projectile) (BuildEffects.cpp
// ScaleOf): stopping_scale.<material>.<projectile>.
void DeclareStoppingScale(const TargetEffect &te, std::vector<ApparatusDeclaration> &out) {
  if (!te.IsTargetIntegration() || !te.HasStoppingScale()) return;
  ApparatusDeclaration d;
  d.key = te.GetStoppingScaleKey();
  d.owner = {{"quantity", "stopping_scale"}, {"material", te.GetStoppingScaleMaterial()}, {"projectile", te.GetStoppingScaleProjectile()}};
  d.unit = "1";
  d.nominal = te.GetStoppingScaleNominal();
  d.vary = te.IsStoppingScaleVaried();
  d.sigma = te.GetStoppingScaleSigma();
  d.layout = te.GetLayoutStoppingScale();
  if (d.vary && te.IsStoppingScaleSigmaDefault()) {
    std::ostringstream note;
    note << "Note: stopping scale " << te.GetStoppingScaleName() << " is fitted with the prior width " << d.sigma
         << ", SRIM's stated accuracy for the projectile class (stopping.md section 7); "
         << "write 'scale <s> fit prior <s> <width>' to give another.";
    d.note = note.str();
  }
  out.push_back(d);
}

bool WriteStoppingScale(const TargetEffect &, int, double v, TargetEffectValues &values) {
  if (values.stoppingScale == v) return false;
  values.stoppingScale = v;
  return true;
}

double ReadStoppingScale(const TargetEffect &, int, const TargetEffectValues &values) { return values.stoppingScale; }

void FollowStoppingScale(WorkingCopy &data, const std::vector<char> &changed, const Config &configure) {
  data.FollowStoppingScales(changed, configure);
}

}  // namespace

ApparatusRow StoppingScaleRow() {
  return {"stopping_scale", "stopping scale", "Stopping-Scale", "epsilon -> s epsilon for that material and projectile",
          DeclareStoppingScale, WriteStoppingScale, ReadStoppingScale, FollowStoppingScale};
}

// Where s x epsilon enters: the stopping of every sub-point (read by the target kernel and the layer's E_out), and the
// thickness and the sigma(E_bar) depth of every point, each recomputed from its energy as the set-up made it.
void WorkingCopy::FollowStoppingScales(const std::vector<char> &changed, const Config &configure) {
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
