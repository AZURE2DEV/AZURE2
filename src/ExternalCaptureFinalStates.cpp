// The external-capture set-up of the compound system (CompoundSystem.h): the final states of the photon pairs
// (ParseExternalCapture), the formal-basis final-state blocks and capture masks, and their values at a fit vector.
#include "ExternalCaptureAmplitudes.h"

#include <map>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include <cmath>
#include <iostream>
#include <memory>

#include "AChannel.h"
#include "GSLException.h"
#include "ALevel.h"
#include "CompoundSystem.h"
#include "EngineChannelFunctions.h"
#include "SystemOfCompound.h"
#include "Config.h"
#include "CoulFunc.h"
#include "ECIntegral.h"
#include "Kinematics.h"
#include "JGroup.h"
#include "PPair.h"
#include "NFIntegral.h"
#include "ProjectSetup.h"
#include "WhitFunc.h"
#include "azure2/experiment/FinalStates.h"
#include "azure2/experiment/ProjectReader.h"
#include "azure2/io/LevelRows.h"
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

/// The first warning of an evaluation-time iteration, once per run (a fit evaluates thousands of times).
void WarnOnce(const phys::LevelSet &set, const std::shared_ptr<Config::NoteLog> &printedOnce) {
  if (set.warnings.empty() || !Config::FirstTimeIn(printedOnce, "formal basis, evaluation")) return;
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

namespace {
/*!
 * The two factors an external-capture final state needs in one particle channel (pair, l) of radius `radius`: the conversion from its
 * reduced width amplitude to the normalization of its exterior wave function, and the N_f integral,
 * int_a^inf u^2 / u(a)^2 dr.
 *
 * Bound (E_f = levelEnergy - threshold < 0): exterior C W(r), so the conversion is sqrt(2 mu a)/hbar / W(a)
 * (gamma -> ANC) and the integral that of W^2 / W(a)^2.
 *
 * Unbound (E_f > 0), a narrow resonance: ECIntegral integrates the outgoing wave normalized at a, u_f = O(r)/O(a), so
 * the conversion is sqrt(2 mu a)/hbar alone, and the N_f integral is continued as (hbar^2 / 2 mu a) dS/dE at E_f --
 * for a bound state that is exactly int W^2/W(a)^2 (checked to 1e-8 for p+12C l = 1 at three energies on dev), and S
 * and dS/dE are smooth through the threshold, so N_f joins the bound case continuously.
 */
void ECFinalStateFactors(const PPair *pair, double radius, int l, int twoJ, double levelEnergy, const Config &configure,
                         double &nfIntegral, double &conversion) {
  const double finalEnergy = levelEnergy - pair->GetSepE() - pair->GetExE();
  const double toAmplitude = sqrt(2.0 * pair->GetRedMass() * radius * uconv / pow(hbarc, 2.0));
  if (finalEnergy > 0.0) {
    CoulFunc coul(pair);
    nfIntegral = pow(hbarc, 2.0) / (2.0 * pair->GetRedMass() * uconv * radius) * coul.PEShift_dE(l, twoJ, radius, finalEnergy);
    // With a potential that also acts inside the radius the surface amplitude carries |(f - b)/(f - L)| =
    // sqrt(P'/P) (PEShift_dE above is already S').
    conversion = toAmplitude;
    if (pair->GetHybridPotential(l) && pair->GetHybridPotential(l)->interior_boundary) {
      const double primed = coul.Penetrability(l, twoJ, radius, finalEnergy);
      const double outer = coul.OuterPenetrability(l, twoJ, radius, finalEnergy);
      if (outer > 0.0) conversion *= sqrt(primed / outer);
    }
  } else {
    NFIntegral newNFIntegral(pair, radius);
    nfIntegral = newNFIntegral(l, twoJ, levelEnergy);
    WhitFunc newWhitFunc(pair);
    conversion = toAmplitude / newWhitFunc(l, twoJ, radius, fabs(finalEnergy));
    // A potential that also acts inside the radius (interior_potential.md section 8): the surface amplitude is
    // gamma |1/D| (the sign stays the ANC's, as in the level conversion), and the norm integral becomes
    // (hbar^2 / 2 mu a) dS'/dE = (N_f + the potential wave's interior part) / D^2.
    const BoundInteriorFactor interior = BoundInterior(*pair, l, twoJ, radius, fabs(finalEnergy));
    conversion *= fabs(interior.dInverse);
    nfIntegral = (nfIntegral + interior.interiorNorm) * interior.dInverse * interior.dInverse;
  }
}

/*!
 * An external-capture final state above the threshold of `particlePair`, a pair with a channel in its J^pi group: a
 * narrow resonance (ECIntegral::IntegrateUnbound).  Said once per capture pair and particle pair (`noted`); refused
 * (false) where the outgoing wave cannot be evaluated: eta beyond the open-channel limit, i.e. within about a keV of
 * the threshold of a charged pair.
 */
bool NoteUnboundFinalState(const Config &configure, const PPair *exitPair, const PPair *particlePair,
                           std::set<std::pair<int, int>> &noted) {
  const double finalEnergy = exitPair->GetHeavyExE() - particlePair->GetSepE() - particlePair->GetExE();
  if (!(finalEnergy > 0.0)) return true;
  const double eta = Kinematics::Channel(particlePair).Eta(finalEnergy);
  if (eta > 100.0) {
    configure.outStream << "ERROR: external capture into pair " << exitPair->GetPairKey() << ": its final state E_x = "
                        << exitPair->GetHeavyExE() << " MeV lies only " << finalEnergy * 1e3 << " keV above the threshold of pair "
                        << particlePair->GetPairKey() << " (eta = " << eta << " > 100), too close for its outgoing wave."
                        << "  Treat it as bound (E_x just below the threshold) or move it." << std::endl;
    return false;
  }
  if (noted.insert(std::make_pair(exitPair->GetPairKey(), particlePair->GetPairKey())).second)
    configure.outStream << "NOTE: external capture into pair " << exitPair->GetPairKey() << ": its final state E_x = "
                        << exitPair->GetHeavyExE() << " MeV is UNBOUND, " << finalEnergy * 1e3 << " keV above the threshold of pair "
                        << particlePair->GetPairKey() << ".  It is treated as a narrow resonance: outgoing-wave exterior,"
                        << " strength from the level's observed width in that pair.  Valid while that width is small next to the photon"
                        << " energy and next to " << finalEnergy * 1e3 << " keV." << std::endl;
  return true;
}
}  // namespace

/*!
 * Marks the final state of each photon pair with external capture (`ecPairs`: pair number -> multipole mask) and
 * stores its N_f integrals and Whittaker conversion factors; the final states are decided by
 * experiment::ExternalCaptureFinalStates on the <levels> `rows`, as the validation decides them.
 */

int CompoundSystem::ParseExternalCapture(const Config &configure, std::map<int, int> &ecPairs,
                               const std::vector<azure2::io::LevelRow> &rows) {
  std::vector<int> ecKeys;
  for (const auto &ec : ecPairs) ecKeys.push_back(this->Pair(ec.first)->GetPairKey());
  const azure2::experiment::LevelTables tables = azure2::experiment::LevelTablesOfRows(rows);
  const std::vector<azure2::experiment::FinalStateOfPair> finals = azure2::experiment::ExternalCaptureFinalStates(tables, ecKeys);
  std::set<int> usedExitKeys = ProjectExitPairKeys(configure);  // the exit pairs of the active data sets (their specs)
  std::set<std::pair<int, int>> unboundNoted;  // (capture pair, particle pair) an unbound final state was noted for
  // A total-capture rate (exit key -1) exits into every capture pair of its entrance nucleus.
  if ((configure.paramMask & Config::CALCULATE_REACTION_RATE) && configure.rateParams.exitPair == -1)
    for (int key : CapturePairKeys(configure.rateParams.entrancePair)) usedExitKeys.insert(key);
  for (std::map<int, int>::iterator ec = ecPairs.begin(); ec != ecPairs.end(); ec++) {
    PPair *exitPair = this->Pair(ec->first);
    if (exitPair->GetPType() != 10) {
      configure.outStream << "Final state is not a capture pair." << std::endl;
      continue;
    }
    // The final state must be a level of the project (user guide, "Gamma-Ray Calculations"; it carries the ANC and the
    // gamma widths): refused for a pair an active data set (or the rate) exits into, else warned about and switched off,
    // so no level is ever made up.  tests/ec_final_state_listed.
    double jValue = exitPair->GetJ(2);
    int parity = exitPair->GetPi(2);
    JGroup newJGroup(jValue, parity);
    newJGroup.SetNucleus(exitPair->GetNucleus(), exitPair->GetNucleusZ(), exitPair->GetNucleusA());
    int jGroupNum = this->IsJGroup(newJGroup);
    const azure2::experiment::FinalStateOfPair &finalState = *azure2::experiment::FinalStateOf(finals, exitPair->GetPairKey());
    if (finalState.levels.empty()) {
      const bool used = usedExitKeys.count(exitPair->GetPairKey()) > 0;
      configure.outStream << (used ? "ERROR" : "WARNING") << ": external capture into pair " << exitPair->GetPairKey()
                          << " needs its final state, J^pi = " << jValue << (parity > 0 ? "+" : "-") << " at E_x = "
                          << exitPair->GetHeavyExE() << " MeV, as a level of the project (with its ANC and gamma widths); none is listed"
                          << (jGroupNum ? " within 1 keV in that J^pi group." : " and that J^pi group does not exist.")
                          << (used ? "" : "  No segment uses this pair: external capture into it is switched off for this run.")
                          << std::endl;
      if (used) return -1;
      continue;
    }
    // Two levels of one J^pi within 1 keV of each other (different level IDs): the final state is ambiguous.
    if (finalState.levels.size() > 1) {
      configure.outStream << "ERROR: external capture into pair " << exitPair->GetPairKey() << ": " << finalState.levels.size()
                          << " levels of J^pi = " << jValue << (parity > 0 ? "+" : "-") << " lie within 1 keV of its final-state energy E_x = "
                          << exitPair->GetHeavyExE() << " MeV, so the final state is ambiguous.  Move the extra level(s) by more than 1 keV"
                          << " or give the pair the energy of the one it decays to." << std::endl;
      return -1;
    }
    const azure2::experiment::LevelInfo &level = tables.levels[finalState.levels[0]];
    if (level.group != jGroupNum) throw std::logic_error("ParseExternalCapture: the level tables number J^pi groups unlike Fill");
    const int levelNum = level.level;
    if (finalState.sharedWithPairKey != 0) {
      FinalStateOfAnotherPair(*this, configure, finalState.sharedWithPairKey, levelNum, ec->first);
      return -1;
    }
    this->Group(jGroupNum)->GetLevel(levelNum)->SetECParams(ec->first, ec->second);
    for (int ch = 1; ch <= this->Group(jGroupNum)->NumChannels(); ch++) {
      PPair *theFinalPair = this->Pair(this->Group(jGroupNum)->GetChannel(ch)->GetPairNum());
      const double finalRadius = this->Group(jGroupNum)->GetChannel(ch)->GetRadius();
      double nfIntegralValue = 0.;
      double ecConvert = 0.;
      if (theFinalPair->GetPType() == 0) {
        if (!NoteUnboundFinalState(configure, exitPair, theFinalPair, unboundNoted)) return -1;
        ECFinalStateFactors(theFinalPair, finalRadius, this->Group(jGroupNum)->GetChannel(ch)->GetL(), theFinalPair->SpinOrbitTwoJ(jValue),
                            exitPair->GetHeavyExE(), configure, nfIntegralValue, ecConvert);
      }
      this->Group(jGroupNum)->GetLevel(levelNum)->AddNFIntegral(nfIntegralValue);
      this->Group(jGroupNum)->GetLevel(levelNum)->AddECConversionFactor(ecConvert);
    }
  }
  return 0;
}


// ------------------------------------------------------------------ formal basis: Brune quantities per evaluation
// (CompoundSystem members: the final states' Brune amplitudes of step 2 and the --ignore-externals projectors of step 4.)

/*!
 * The R-matrix blocks of the formal basis that hold an external-capture final state with other levels
 * (cnuc_split.md section 8): in them the final state's formal amplitude is a mixture of the block's poles, and its
 * Brune amplitude, the residue of its own pole, is what external capture needs.  Setup state, shared by the clones.
 */
struct CompoundSystem::ECFinalBlocks {
  std::shared_ptr<const SystemOfCompound> structure;
  std::vector<std::size_t> blocks;
  bool rmc = false;
  std::shared_ptr<Config::NoteLog> printedOnce;  ///< the run's once-only warnings (Config::printedOnce)
};

void CompoundSystem::SetUpECFinalBlocks(const Config &configure, const LevelInputs &inputs) {
  ecFinalBlocks_.reset();
  // The final state's energy in external capture (the photon energy, the bound-state wave): its Brune energy, the
  // level energy now except in the formal basis (below).
  for (int j = 1; j <= this->NumJGroups(); j++)
    for (int la = 1; la <= this->Group(j)->NumLevels(); la++)
      if (this->Group(j)->GetLevel(la)->IsECLevel())
        this->Group(j)->GetLevel(la)->SetECEnergy(inputs.E(j, la));
  if ((configure.paramMask & Config::USE_BRUNE_FORMALISM) || !(configure.paramMask & Config::USE_EXTERNAL_CAPTURE))
    return;
  auto record = std::make_shared<ECFinalBlocks>(ECFinalBlocks{
      this->SystemStructureHandle(), {},
      !!(configure.paramMask & Config::USE_RMC_FORMALISM),
      configure.printedOnce});
  for (std::size_t b = 0; b < record->structure->system.num_blocks(); b++) {
    const std::vector<int> &levels = record->structure->map.levelOfBlockLevel[b];
    // The lowest level of a block (B_c = S_c at its energy) has formal parameters equal to its Brune ones: the
    // iteration of brune_from_formal moves B_c by S_c(E) - B_c = 0 (level_parameters.md section 3), as long as its
    // energy stays there.  So does the only level of a block.
    const int j = record->structure->map.groupOfBlock[b];
    const JGroup *group = this->Group(j);
    std::size_t lowest = 0;
    for (std::size_t lb = 1; lb < levels.size(); lb++)
      if (inputs.E(j, levels[lb]) < inputs.E(j, levels[lowest])) lowest = lb;
    // A lowest-level final state whose energy a fit may move leaves the energy where B_c was taken, and its formal
    // amplitude then differs from its Brune one at second order: it is recomputed too (external_capture.md section 4).
    // A group that holds an isospin level has B_c = -l_c (IsospinTies.h): none of its levels sits where B_c was taken.
    const bool atBoundary = !isospin_->Groups().count(j);
    for (std::size_t lb = 0; lb < levels.size(); lb++)
      if ((lb != lowest || !atBoundary || !group->GetLevel(levels[lb])->EnergyFixed()) && group->GetLevel(levels[lb])->IsECLevel()) {
        record->blocks.push_back(b);
        break;
      }
  }
  if (record->blocks.empty()) return;
  ecFinalBlocks_ = record;
  // The Brune energies of their final states, from the formal parameters (EvaluationIteration).
  const phys::ChannelFunctions cf = EngineChannelFunctions(*this, record->structure->map, record->rmc);
  const phys::FormalToBruneOptions options = EvaluationIteration(false);
  for (std::size_t b : record->blocks) {
    const int j = record->structure->map.groupOfBlock[b];
    JGroup *group = this->Group(j);
    const std::vector<int> &levels = record->structure->map.levelOfBlockLevel[b];
    const std::size_t nc = group->NumChannels();
    std::vector<double> E, gamma, B(nc);
    for (int la : levels) {
      E.push_back(inputs.E(j, la));
      for (std::size_t k = 0; k < nc; k++) gamma.push_back(inputs.Gamma(j, la, (int)k + 1));
    }
    for (std::size_t k = 0; k < nc; k++) B[k] = group->GetChannel((int)k + 1)->GetBoundaryCondition();
    const phys::LevelSet brune = phys::brune_from_formal(record->structure->system, b, E, gamma, B, cf, options);
    ThrowIfFailed(brune);
    WarnOnce(brune, record->printedOnce);
    for (std::size_t lb = 0; lb < levels.size(); lb++)
      if (group->GetLevel(levels[lb])->IsECLevel()) group->GetLevel(levels[lb])->SetECEnergy(brune.E[lb]);
  }
}

/*!
 * Formal basis: the Brune amplitudes of the final states in ecFinalBlocks_ at the current fit parameters (physics/
 * brune_from_formal), as their external-capture amplitude, with N_f^1/2 formed from them.  The iteration is
 * EvaluationIteration's, so that the capture agrees with the Brune basis to rounding (the reports use 1e-6 MeV).
 */
void CompoundSystem::ECFinalStatesFromFormal(LevelState &state) const {
  const ECFinalBlocks &r = *ecFinalBlocks_;
  const phys::ChannelFunctions cf = EngineChannelFunctions(*this, r.structure->map, r.rmc);
  const phys::FormalToBruneOptions options = EvaluationIteration(false);
  for (std::size_t b : r.blocks) {
    const int j = r.structure->map.groupOfBlock[b];
    const JGroup *group = this->GetJGroup(j);
    const std::vector<int> &levels = r.structure->map.levelOfBlockLevel[b];
    const std::size_t nc = group->NumChannels();
    std::vector<double> E, gamma, B(nc);
    for (int la : levels) {
      E.push_back(state.E(j, la));
      for (std::size_t k = 0; k < nc; k++) gamma.push_back(state.Gamma(j, la, (int)k + 1));
    }
    for (std::size_t k = 0; k < nc; k++) B[k] = group->GetChannel((int)k + 1)->GetBoundaryCondition();
    const phys::LevelSet brune = phys::brune_from_formal(r.structure->system, b, E, gamma, B, cf, options);
    WarnOnce(brune, r.printedOnce);
    for (std::size_t lb = 0; lb < levels.size(); lb++) {
      const ALevel *level = group->GetLevel(levels[lb]);
      if (!level->IsECLevel()) continue;
      const std::vector<double> g(brune.gamma.begin() + lb * nc, brune.gamma.begin() + (lb + 1) * nc);
      for (std::size_t k = 0; k < nc; k++) state.SetECGamma(j, levels[lb], (int)k + 1, g[k]);
      state.SetSqrtNF(j, levels[lb], SqrtNFFactor(*this, group, level, 1.0, g));
    }
  }
}

/*!
 * Formal basis with --ignore-externals: every (block, photon channel) through which channel capture runs, for a block
 * of more than one level (external_capture.md section 3).  Setup state, shared by the clones; the projectors
 * themselves are recomputed with the fit parameters (FormalCaptureMasksFromParams).
 */
struct CompoundSystem::CaptureMaskSetup {
  std::shared_ptr<const SystemOfCompound> structure;
  std::vector<CaptureMaskSlot> slots;
  bool rmc = false;
  std::shared_ptr<Config::NoteLog> printedOnce;  ///< the run's once-only warnings (Config::printedOnce)
};

void CompoundSystem::SetUpFormalCaptureMasks(const Config &configure, LevelState &state) {
  captureMaskSetup_.reset();
  state.ShapeCaptureProjectors(0);
  if ((configure.paramMask & Config::USE_BRUNE_FORMALISM) || !(configure.paramMask & Config::IGNORE_ZERO_WIDTHS) ||
      !(configure.paramMask & Config::TRANSFORM_PARAMETERS))
    return;
  auto setup = std::make_shared<CaptureMaskSetup>(CaptureMaskSetup{
      this->SystemStructureHandle(), {},
      !!(configure.paramMask & Config::USE_RMC_FORMALISM),
      configure.printedOnce});
  const SystemOfCompoundMap &map = setup->structure->map;
  for (int p = 1; p <= this->NumPairs(); p++) {
    const azure2::EntrancePathways &P = this->Pathways().Entrance(p);
    for (const azure2::ExternalRow &ec : P.external) {
      if (!ec.channel_capture || !ec.photon_channel) continue;
      const int group = P.internal[ec.cc_row].group, channel = ec.photon_channel;
      for (std::size_t b = 0; b < map.groupOfBlock.size(); b++) {
        if (map.groupOfBlock[b] != group || map.levelOfBlockLevel[b].size() < 2) continue;
        bool known = false;
        for (const CaptureMaskSlot &s : setup->slots) known = known || (s.group == group && s.channel == channel);
        if (!known) setup->slots.push_back({b, group, channel, map.levelOfBlockLevel[b]});
      }
    }
  }
  if (setup->slots.empty()) return;
  captureMaskSetup_ = setup;
  state.ShapeCaptureProjectors(setup->slots.size());
}

/*!
 * The projector Pi = T^-1 P T of each slot at the current fit parameters (physics/ capture_mask, Eq. 5), T from
 * brune_from_formal (EvaluationIteration) as for the final states; empty when no level of the block is masked.
 */
void CompoundSystem::FormalCaptureMasksFromParams(LevelState &state) const {
  const CaptureMaskSetup &r = *captureMaskSetup_;
  const phys::ChannelFunctions cf = EngineChannelFunctions(*this, r.structure->map, r.rmc);
  const phys::FormalToBruneOptions options = EvaluationIteration(true);
  for (std::size_t i = 0; i < r.slots.size(); i++) {
    const CaptureMaskSlot &mask = r.slots[i];
    const JGroup *group = this->GetJGroup(mask.group);
    const std::size_t n = mask.levels.size(), nc = group->NumChannels();
    std::vector<unsigned char> masked(n, 0);
    bool any = false;
    for (std::size_t lb = 0; lb < n; lb++) {
      masked[lb] = state.ZeroInputWidth(mask.group, mask.levels[lb], mask.channel) ? 1 : 0;
      any = any || masked[lb];
    }
    state.SetCaptureProjector((int)i, std::vector<double>());
    if (!any) continue;
    // A masked level whose particle amplitudes are all zero is decoupled: no shift-matrix element links it to another
    // level, so T is the identity on its row and column and Pi = P exactly (external_capture.md section 7).
    bool coupled = false;
    for (std::size_t lb = 0; lb < n && !coupled; lb++)
      if (masked[lb]) coupled = MaskedLevelCoupled(state, mask.group, mask.levels[lb]);
    if (!coupled) continue;
    std::vector<double> E, gamma, B(nc);
    for (int la : mask.levels) {
      E.push_back(state.E(mask.group, la));
      for (std::size_t k = 0; k < nc; k++) gamma.push_back(state.Gamma(mask.group, la, (int)k + 1));
    }
    for (std::size_t k = 0; k < nc; k++) B[k] = group->GetChannel((int)k + 1)->GetBoundaryCondition();
    const phys::LevelSet brune = phys::brune_from_formal(r.structure->system, mask.block, E, gamma, B, cf, options);
    WarnOnce(brune, r.printedOnce);
    state.SetCaptureProjector((int)i, phys::capture_mask(brune.transform, n, masked));
  }
}

/*!
 * Does the adjoint miss a dependence of a formal-basis capture mask on the fit parameters?  It does when a slot has a
 * projector (T depends on the parameters), or when a masked level is decoupled now but one of its particle amplitudes
 * is free: P is exact at the point, but dPi = [P, dT] is not zero there (external_capture.md section 7).
 */
bool CompoundSystem::FormalCaptureMasksOutsideAdjoint(const LevelState &state,
                                             const std::function<bool(int, int, int)> &freeGamma) const {
  if (!captureMaskSetup_) return false;
  for (std::size_t i = 0; i < captureMaskSetup_->slots.size(); i++) {
    const CaptureMaskSlot &mask = captureMaskSetup_->slots[i];
    if (!state.CaptureProjector((int)i).empty()) return true;
    const JGroup *group = this->GetJGroup(mask.group);
    for (int la : mask.levels) {
      if (!state.ZeroInputWidth(mask.group, la, mask.channel)) continue;
      for (int ch = 1; ch <= group->NumChannels(); ch++)
        if (group->GetChannel(ch)->GetRadType() == 'P' && freeGamma(mask.group, la, ch)) return true;
    }
  }
  return false;
}

/// Has level `la` of J-group `j` a particle amplitude other than zero (current fit values)?
bool CompoundSystem::MaskedLevelCoupled(const LevelState &state, int j, int la) const {
  const JGroup *group = this->GetJGroup(j);
  for (int ch = 1; ch <= group->NumChannels(); ch++)
    if (group->GetChannel(ch)->GetRadType() == 'P' && state.Gamma(j, la, ch) != 0.0) return true;
  return false;
}

int CompoundSystem::FormalCaptureMaskSlot(int group, int channel) const {
  RefuseUnfinished("CompoundSystem::FormalCaptureMaskSlot");
  if (!captureMaskSetup_) return -1;
  for (std::size_t i = 0; i < captureMaskSetup_->slots.size(); i++)
    if (captureMaskSetup_->slots[i].group == group && captureMaskSetup_->slots[i].channel == channel) return (int)i;
  return -1;
}

const CompoundSystem::CaptureMaskSlot &CompoundSystem::FormalCaptureMaskOf(int slot) const {
  RefuseUnfinished("CompoundSystem::FormalCaptureMaskOf");
  return captureMaskSetup_->slots[slot];
}
