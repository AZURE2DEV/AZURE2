#include "CNucExternalCapture.h"

#include <atomic>
#include <cmath>
#include <iostream>
#include <memory>

#include "AChannel.h"
#include "GSLException.h"
#include "ALevel.h"
#include "CNuc.h"
#include "CNucLevelParameters.h"
#include "CNucSystem.h"
#include "Config.h"
#include "ECIntegral.h"
#include "JGroup.h"
#include "PPair.h"
#include "azure2/physics/ExternalCapture.h"
#include "azure2/physics/LevelParameters.h"

namespace phys = azure2::physics;

namespace {
/*!
 * The formal -> Brune iteration inside a formal-basis evaluation (level_parameters.md section 9): 1e-12 MeV or 1e-12
 * relative, the last iterate at the rounding floor, NaN instead of formal values.
 */
phys::FormalToBruneOptions EvaluationIteration(bool wantTransform) {
  phys::FormalToBruneOptions options;
  options.energy_tolerance_MeV = 1e-12;
  options.relative_tolerance = 1e-12;
  options.evaluation = true;
  options.want_transform = wantTransform;
  return options;
}

/// The first warning of an evaluation-time iteration, once per process (a fit evaluates thousands of times).
void WarnOnce(const phys::LevelSet &set) {
  static std::atomic<bool> told{false};
  if (set.warnings.empty() || told.exchange(true)) return;
  std::cerr << "WARNING (formal basis, evaluation): " << set.warnings.front()
            << " (further such warnings are not shown)" << std::endl;
}

/// At setup a level that cannot be transformed stops the run, as TransformIn does.
void ThrowIfFailed(const phys::LevelSet &set) {
  for (std::size_t i = 0; i < set.status.size(); i++)
    if (set.status[i] == phys::LevelSet::Status::EigenFailed || set.status[i] == phys::LevelSet::Status::NotConverged)
      throw GSLException("External-capture final states (formal basis): " +
                         (i < set.warnings.size() ? set.warnings[i] : std::string("a level could not be transformed")));
}
}  // namespace

double SqrtNFFactor(CNuc &compound, JGroup *group, const ALevel *level, double n0, const std::vector<double> &gamma) {
  std::vector<phys::BoundChannel> channels;
  for (int ch = 1; ch <= group->NumChannels() && ch <= level->NumNFIntegrals(); ch++) {
    PPair *pair = compound.GetPair(group->GetChannel(ch)->GetPairNum());
    channels.push_back({pair->GetChRad(), pair->GetRedMass(), gamma[ch - 1], level->GetNFIntegral(ch)});
  }
  return 1.0 / std::sqrt(phys::bound_state_normalization(n0, channels));
}

complex ExternalReducedWidth(CNuc &compound, JGroup *group, ALevel *level, AChannel *channel, bool initial,
                             const Config &configure) {
  const char type = channel->GetRadType();
  if (type != 'E' && !(type == 'M' && channel->GetL() == 1)) return complex(0.0, 0.0);
  JGroup *finalGroup = nullptr;
  ALevel *finalLevel = nullptr;
  for (int j = 1; j <= compound.NumJGroups() && !finalLevel; j++)
    for (int la = 1; la <= compound.GetJGroup(j)->NumLevels() && !finalLevel; la++)
      if (compound.GetJGroup(j)->GetLevel(la)->IsECLevel() &&
          channel->GetPairNum() == compound.GetJGroup(j)->GetLevel(la)->GetECPairNum()) {
        finalGroup = compound.GetJGroup(j);
        finalLevel = finalGroup->GetLevel(la);
      }
  if (!finalLevel) return complex(0.0, 0.0);
  const int L = channel->GetL();
  const unsigned char mask = finalLevel->GetECMultMask();
  if (!((type == 'E' && L == 1 && (mask & isE1)) || (type == 'M' && L == 1 && (mask & isM1)) ||
        (type == 'E' && L == 2 && (mask & isE2))))
    return complex(0.0, 0.0);
  const double levelEnergy = initial ? level->GetE() : level->GetTransformE();
  const double finalEnergy = initial ? finalLevel->GetE() : finalLevel->GetTransformE();
  std::vector<phys::CaptureChannel> in;
  std::vector<AChannel *> inChannels;
  for (int ch = 1; ch <= group->NumChannels(); ch++) {
    AChannel *c = group->GetChannel(ch);
    if (c->GetRadType() != 'P') continue;
    in.push_back({c->GetPairNum(), c->GetL(), c->GetS(), initial ? level->GetGamma(ch) : level->GetTransformGamma(ch)});
    inChannels.push_back(c);
  }
  std::vector<phys::FinalChannel> out;
  std::vector<AChannel *> outChannels;
  for (int chp = 1; chp <= finalGroup->NumChannels(); chp++) {
    AChannel *c = finalGroup->GetChannel(chp);
    if (c->GetRadType() != 'P') continue;
    const double g = initial ? finalLevel->GetGamma(chp) : finalLevel->GetTransformGamma(chp);
    // s_fc = gamma~_fc N_f^(-1/2) K_c (Eq. 2), in the engine's order.
    out.push_back({c->GetPairNum(), c->GetL(), c->GetS(),
                   g * finalLevel->GetSqrtNFFactor() * finalLevel->GetECConversionFactor(chp)});
    outChannels.push_back(c);
  }
  const phys::RadialIntegral integral = [&](std::size_t c, std::size_t cf) {
    ECIntegral theECIntegral(compound.GetPair(outChannels[cf]->GetPairNum()), configure);
    return theECIntegral(inChannels[c]->GetL(), outChannels[cf]->GetL(), inChannels[c]->GetS(),
                         outChannels[cf]->GetS(), group->GetJ(), finalGroup->GetJ(), L, type, levelEnergy,
                         finalEnergy, true);
  };
  return phys::external_reduced_width(type == 'M' ? phys::Multipole::M : phys::Multipole::E, L, finalGroup->GetJ(),
                                      in, out, integral);
}

// ------------------------------------------------------------------ formal basis: Brune quantities per evaluation
// (CNuc members: the final states' Brune amplitudes of step 2 and the --ignore-externals projectors of step 4.)

/*!
 * The R-matrix blocks of the formal basis that hold an external-capture final state with other levels
 * (cnuc_split.md section 8): in them the final state's formal amplitude is a mixture of the block's poles, and its
 * Brune amplitude, the residue of its own pole, is what external capture needs.  Setup state, shared by the clones.
 */
struct CNuc::ECFinalBlocks {
  CNucSystem structure;
  std::vector<std::size_t> blocks;
  bool gslCoulomb = false, rmc = false;
};

void CNuc::SetUpECFinalBlocks(const Config &configure) {
  ecFinalBlocks_.reset();
  // The final state's energy in external capture (the photon energy, the bound-state wave): its Brune energy, the
  // level energy now except in the formal basis (below).
  for (int j = 1; j <= this->NumJGroups(); j++)
    for (int la = 1; la <= this->GetJGroup(j)->NumLevels(); la++)
      if (this->GetJGroup(j)->GetLevel(la)->IsECLevel())
        this->GetJGroup(j)->GetLevel(la)->SetECEnergy(this->GetJGroup(j)->GetLevel(la)->GetE());
  if ((configure.paramMask & Config::USE_BRUNE_FORMALISM) || !(configure.paramMask & Config::USE_EXTERNAL_CAPTURE))
    return;
  auto record = std::make_shared<ECFinalBlocks>(ECFinalBlocks{
      BuildSystemFromCNuc(*this, configure, CNucSystemContent::Structure), {},
      !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC), !!(configure.paramMask & Config::USE_RMC_FORMALISM)});
  for (std::size_t b = 0; b < record->structure.system.num_blocks(); b++) {
    const std::vector<int> &levels = record->structure.map.levelOfBlockLevel[b];
    // The lowest level of a block (B_c = S_c at its energy) has formal parameters equal to its Brune ones: the
    // iteration of brune_from_formal moves B_c by S_c(E) - B_c = 0 (level_parameters.md section 3), as long as its
    // energy stays there.  So does the only level of a block.
    JGroup *group = this->GetJGroup(record->structure.map.groupOfBlock[b]);
    std::size_t lowest = 0;
    for (std::size_t lb = 1; lb < levels.size(); lb++)
      if (group->GetLevel(levels[lb])->GetE() < group->GetLevel(levels[lowest])->GetE()) lowest = lb;
    // A lowest-level final state whose energy a fit may move leaves the energy where B_c was taken, and its formal
    // amplitude then differs from its Brune one at second order: it is recomputed too (external_capture.md section 4).
    for (std::size_t lb = 0; lb < levels.size(); lb++)
      if ((lb != lowest || !group->GetLevel(levels[lb])->EnergyFixed()) && group->GetLevel(levels[lb])->IsECLevel()) {
        record->blocks.push_back(b);
        break;
      }
  }
  if (record->blocks.empty()) return;
  ecFinalBlocks_ = record;
  // The Brune energies of their final states, from the formal parameters (EvaluationIteration).
  const phys::ChannelFunctions cf = ChannelFunctionsFromCNuc(*this, record->structure.map, record->gslCoulomb, record->rmc);
  const phys::FormalToBruneOptions options = EvaluationIteration(false);
  for (std::size_t b : record->blocks) {
    JGroup *group = this->GetJGroup(record->structure.map.groupOfBlock[b]);
    const std::vector<int> &levels = record->structure.map.levelOfBlockLevel[b];
    const std::size_t nc = group->NumChannels();
    std::vector<double> E, gamma, B(nc);
    for (int la : levels) {
      E.push_back(group->GetLevel(la)->GetE());
      for (std::size_t k = 0; k < nc; k++) gamma.push_back(group->GetLevel(la)->GetGamma((int)k + 1));
    }
    for (std::size_t k = 0; k < nc; k++) B[k] = group->GetChannel((int)k + 1)->GetBoundaryCondition();
    const phys::LevelSet brune = phys::brune_from_formal(record->structure.system, b, E, gamma, B, cf, options);
    ThrowIfFailed(brune);
    WarnOnce(brune);
    for (std::size_t lb = 0; lb < levels.size(); lb++)
      if (group->GetLevel(levels[lb])->IsECLevel()) group->GetLevel(levels[lb])->SetECEnergy(brune.E[lb]);
  }
}

/*!
 * Formal basis: the Brune amplitudes of the final states in ecFinalBlocks_ at the current fit parameters (physics/
 * brune_from_formal), as their external-capture amplitude, with N_f^1/2 formed from them.  The iteration is
 * EvaluationIteration's, so that the capture agrees with the Brune basis to rounding (the reports use 1e-6 MeV).
 */
void CNuc::ECFinalStatesFromFormal() {
  const ECFinalBlocks &r = *ecFinalBlocks_;
  const phys::ChannelFunctions cf = ChannelFunctionsFromCNuc(*this, r.structure.map, r.gslCoulomb, r.rmc);
  const phys::FormalToBruneOptions options = EvaluationIteration(false);
  for (std::size_t b : r.blocks) {
    JGroup *group = this->GetJGroup(r.structure.map.groupOfBlock[b]);
    const std::vector<int> &levels = r.structure.map.levelOfBlockLevel[b];
    const std::size_t nc = group->NumChannels();
    std::vector<double> E, gamma, B(nc);
    for (int la : levels) {
      E.push_back(group->GetLevel(la)->GetFitE());
      for (std::size_t k = 0; k < nc; k++) gamma.push_back(group->GetLevel(la)->GetFitGamma((int)k + 1));
    }
    for (std::size_t k = 0; k < nc; k++) B[k] = group->GetChannel((int)k + 1)->GetBoundaryCondition();
    const phys::LevelSet brune = phys::brune_from_formal(r.structure.system, b, E, gamma, B, cf, options);
    WarnOnce(brune);
    for (std::size_t lb = 0; lb < levels.size(); lb++) {
      ALevel *level = group->GetLevel(levels[lb]);
      if (!level->IsECLevel()) continue;
      const std::vector<double> g(brune.gamma.begin() + lb * nc, brune.gamma.begin() + (lb + 1) * nc);
      for (std::size_t k = 0; k < nc; k++) level->SetECGamma((int)k + 1, g[k]);
      level->SetSqrtNFFactor(SqrtNFFactor(*this, group, level, 1.0, g));
    }
  }
}

/*!
 * Formal basis with --ignore-externals: every (block, photon channel) through which channel capture runs, for a block
 * of more than one level (external_capture.md section 3).  Setup state, shared by the clones; the projectors
 * themselves are recomputed with the fit parameters (FormalCaptureMasksFromParams).
 */
struct CNuc::CaptureMaskSetup {
  CNucSystem structure;
  struct Slot {
    std::size_t block;
    int group, channel;
  };
  std::vector<Slot> slots;
  bool gslCoulomb = false, rmc = false;
};

void CNuc::SetUpFormalCaptureMasks(const Config &configure) {
  captureMaskSetup_.reset();
  captureMasks_.clear();
  if ((configure.paramMask & Config::USE_BRUNE_FORMALISM) || !(configure.paramMask & Config::IGNORE_ZERO_WIDTHS) ||
      !(configure.paramMask & Config::TRANSFORM_PARAMETERS))
    return;
  auto setup = std::make_shared<CaptureMaskSetup>(CaptureMaskSetup{
      BuildSystemFromCNuc(*this, configure, CNucSystemContent::Structure), {},
      !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC), !!(configure.paramMask & Config::USE_RMC_FORMALISM)});
  const CNucSystemMap &map = setup->structure.map;
  for (int p = 1; p <= this->NumPairs(); p++) {
    const azure2::EntrancePathways &P = this->Pathways().Entrance(p);
    for (const azure2::ExternalRow &ec : P.external) {
      if (!ec.channel_capture || !ec.photon_channel) continue;
      const int group = P.internal[ec.cc_row].group, channel = ec.photon_channel;
      for (std::size_t b = 0; b < map.groupOfBlock.size(); b++) {
        if (map.groupOfBlock[b] != group || map.levelOfBlockLevel[b].size() < 2) continue;
        bool known = false;
        for (const CaptureMaskSetup::Slot &s : setup->slots) known = known || (s.group == group && s.channel == channel);
        if (!known) setup->slots.push_back({b, group, channel});
      }
    }
  }
  if (setup->slots.empty()) return;
  captureMaskSetup_ = setup;
  for (const CaptureMaskSetup::Slot &s : setup->slots) {
    CaptureMask mask;
    mask.group = s.group;
    mask.channel = s.channel;
    mask.levels = map.levelOfBlockLevel[s.block];
    captureMasks_.push_back(mask);
  }
}

/*!
 * The projector Pi = T^-1 P T of each slot at the current fit parameters (physics/ capture_mask, Eq. 5), T from
 * brune_from_formal (EvaluationIteration) as for the final states; empty when no level of the block is masked.
 */
void CNuc::FormalCaptureMasksFromParams() {
  const CaptureMaskSetup &r = *captureMaskSetup_;
  const phys::ChannelFunctions cf = ChannelFunctionsFromCNuc(*this, r.structure.map, r.gslCoulomb, r.rmc);
  const phys::FormalToBruneOptions options = EvaluationIteration(true);
  for (std::size_t i = 0; i < r.slots.size(); i++) {
    const CaptureMaskSetup::Slot &s = r.slots[i];
    CaptureMask &mask = captureMasks_[i];
    JGroup *group = this->GetJGroup(s.group);
    const std::size_t n = mask.levels.size(), nc = group->NumChannels();
    std::vector<unsigned char> masked(n, 0);
    bool any = false;
    for (std::size_t lb = 0; lb < n; lb++) {
      masked[lb] = group->GetLevel(mask.levels[lb])->ZeroInputWidth(s.channel) ? 1 : 0;
      any = any || masked[lb];
    }
    mask.projector.clear();
    if (!any) continue;
    // A masked level whose particle amplitudes are all zero is decoupled: no shift-matrix element links it to another
    // level, so T is the identity on its row and column and Pi = P exactly (external_capture.md section 7).
    bool coupled = false;
    for (std::size_t lb = 0; lb < n && !coupled; lb++)
      if (masked[lb]) coupled = MaskedLevelCoupled(group, mask.levels[lb]);
    if (!coupled) continue;
    std::vector<double> E, gamma, B(nc);
    for (int la : mask.levels) {
      E.push_back(group->GetLevel(la)->GetFitE());
      for (std::size_t k = 0; k < nc; k++) gamma.push_back(group->GetLevel(la)->GetFitGamma((int)k + 1));
    }
    for (std::size_t k = 0; k < nc; k++) B[k] = group->GetChannel((int)k + 1)->GetBoundaryCondition();
    const phys::LevelSet brune = phys::brune_from_formal(r.structure.system, s.block, E, gamma, B, cf, options);
    WarnOnce(brune);
    mask.projector = phys::capture_mask(brune.transform, n, masked);
  }
}

/*!
 * Does the adjoint miss a dependence of a formal-basis capture mask on the fit parameters?  It does when a slot has a
 * projector (T depends on the parameters), or when a masked level is decoupled now but one of its particle amplitudes
 * is free: P is exact at the point, but dPi = [P, dT] is not zero there (external_capture.md section 7).
 */
bool CNuc::FormalCaptureMasksOutsideAdjoint(const std::function<bool(int, int, int)> &freeGamma) {
  for (const CaptureMask &mask : captureMasks_) {
    if (!mask.projector.empty()) return true;
    JGroup *group = this->GetJGroup(mask.group);
    for (int la : mask.levels) {
      if (!group->GetLevel(la)->ZeroInputWidth(mask.channel)) continue;
      for (int ch = 1; ch <= group->NumChannels(); ch++)
        if (group->GetChannel(ch)->GetRadType() == 'P' && freeGamma(mask.group, la, ch)) return true;
    }
  }
  return false;
}

/// Has level `la` of `group` a particle amplitude other than zero (current fit values)?
bool CNuc::MaskedLevelCoupled(JGroup *group, int la) const {
  for (int ch = 1; ch <= group->NumChannels(); ch++)
    if (group->GetChannel(ch)->GetRadType() == 'P' && group->GetLevel(la)->GetFitGamma(ch) != 0.0) return true;
  return false;
}

const CNuc::CaptureMask *CNuc::FormalCaptureMask(int group, int channel) const {
  for (const CaptureMask &m : captureMasks_)
    if (m.group == group && m.channel == channel) return &m;
  return nullptr;
}

void FinalStateOfAnotherPair(CNuc &compound, const Config &configure, int otherPairKey, int level, int pairNum) {
  PPair *pair = compound.GetPair(pairNum);
  configure.outStream << "ERROR: external capture into pairs " << otherPairKey
                      << " and " << pair->GetPairKey() << ": both have the final state J^pi = " << pair->GetJ(2)
                      << (pair->GetPi(2) > 0 ? "+" : "-") << " at E_x = " << pair->GetHeavyExE() << " MeV (level " << level
                      << " of that J^pi group).  A photon pair has its own final state: remove the duplicated pair, or"
                      << " give it the energy of another level." << std::endl;
}
