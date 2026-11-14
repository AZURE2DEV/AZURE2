// CompoundSystem: the observed quantities of the current fit values (a table), and one fit value set with what derives
// from it, for the model data (ModelData.cpp), which compare observed widths and energies with literature values
// inside a fit.
#include <exception>
#include <vector>

#include "ALevel.h"
#include "CompoundSystem.h"
#include "ExternalCaptureAmplitudes.h"
#include "Config.h"
#include "JGroup.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <initializer_list>
#include <sstream>
#include <algorithm>
#include <set>
#include <map>
#include <tuple>
#include "RunScratch.h"
#include "RunDiagnostics.h"
#include "azure2/physics/Constants.h"
#include <mutex>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "ParameterLabel.h"
#include "ChannelFunc.h"
#include "EngineChannelFunctions.h"
#include "SystemOfCompound.h"
#include "CompoundSystemChecks.h"
#include "LevelReport.h"
#include "LevelInputMessages.h"
#include "CoulFunc.h"
#include "ECIntegral.h"
#include "GSLException.h"
#include "ProjectSetup.h"
#include "azure2/experiment/FinalStates.h"
#include "azure2/experiment/ProjectReader.h"
#include "azure2/io/LevelRows.h"
#include "Minuit2/MnUserParameters.h"
#include "AZUREParams.h"
#include "NFIntegral.h"
#include "ShftFunc.h"
#include "azure2/io/LevelGrouping.h"
#include "azure2/physics/ObservableCoupling.h"
#include "azure2/physics/Pathways.h"

ObservedLevelTable CompoundSystem::EmptyObservedTable() const {
  std::vector<int> numLevels, numChannels;
  for (int j = 1; j <= this->NumJGroups(); j++) {
    numLevels.push_back(this->GetJGroup(j)->NumLevels());
    numChannels.push_back(this->GetJGroup(j)->NumChannels());
  }
  return ObservedLevelTable(numLevels, numChannels);
}

std::optional<ObservedLevelTable> CompoundSystem::ObservedOfFitValues(const Config &configure, LevelState &state,
                                                          const LevelInputs &inputs) const {
  // In Park's basis TransformToBrune recomputes the shift functions, which set the bound-state normalizations
  // N_f^-1/2 with Park's J; a calculation reads the ones its own evaluation set, so they are put back.
  std::vector<double> sqrtNF;
  for (int j = 1; j <= this->NumJGroups(); j++)
    for (int la = 1; la <= this->GetJGroup(j)->NumLevels(); la++)
      sqrtNF.push_back(state.SqrtNF(j, la));
  std::optional<ObservedLevelTable> table = this->EmptyObservedTable();
  try {
    this->TransformToBrune(configure, Messages::kQuiet, state, inputs, *table);
    this->ObservedWidths(configure, state, *table);
  } catch (const std::exception &) {
    table.reset();
  }
  std::size_t k = 0;
  for (int j = 1; j <= this->NumJGroups(); j++)
    for (int la = 1; la <= this->GetJGroup(j)->NumLevels(); la++)
      state.SetSqrtNF(j, la, sqrtNF[k++]);
  return table;
}

void CompoundSystem::SetFitLevelParameter(LevelState &state, const LevelParameterLayout::Slot &slot, double value) const {
  const int j = slot.group, la = slot.level;
  const JGroup *group = this->GetJGroup(j);
  const ALevel *level = group->GetLevel(la);
  if (slot.kind == LevelParameterLayout::Slot::Kind::kEnergy) {
    state.SetE(j, la, value);
  } else {
    state.SetGamma(j, la, slot.channel, value);
    state.SetECGamma(j, la, slot.channel, value);
  }
  // What FillCompoundFromParams derives from the level's values.
  std::vector<double> g;
  for (int c = 1; c <= group->NumChannels(); c++) g.push_back(state.Gamma(j, la, c));
  state.SetSqrtNF(j, la, SqrtNFFactor(*this, group, level, 1.0, g));
  if (ecFinalBlocks_) this->ECFinalStatesFromFormal(state);
  if (captureMaskSetup_) this->FormalCaptureMasksFromParams(state);
}


namespace phys = azure2::physics;

namespace {

//! 2j as an integer, for the physics/ coupling coefficients (docs/upgrade/cnuc_split.md section 9).
int Twice(double j) { return (int)std::lround(2.0 * j); }

}  // namespace

/*!
 * Performs the initial parameter transformations from physical to formal parameters: the observed widths, ANCs and
 * amplitudes of the levels become Brune's (or Park's) reduced width amplitudes, the external-capture amplitude is
 * subtracted from the photon channels, and, with the Brune parametrization off, the levels become formal ones.  The
 * arithmetic is physics/ level_parameters' (docs/upgrade/level_parameters.md sections 2-4), with the engine's channel
 * functions (EngineChannelFunctions); the formal parameters take B_c at the lowest level of a group, and each formal
 * level the orientation of its Brune level (HANDOFF decision 10, cnuc_split.md section 7).  This function words the
 * warnings and keeps the bound-state normalization of external capture.
 */

bool CompoundSystem::TransformIn(const Config &configure, LevelState &state, LevelInputs &inputs) const {
  inputs.StartTransformation();
  // Park (2021): the fit parameter is the observed reduced width amplitude, Gamma = 2 P gamma^2 (Eqs. 31-33), so the
  // factor that turns an observed width into Brune's reduced width is not applied; its J (Eq. 28) must stay positive.
  const bool park = (configure.paramMask & Config::USE_PARK_FORMALISM);
  const bool formal = !(configure.paramMask & Config::USE_BRUNE_FORMALISM);
  try {
    const SystemOfCompound &structure = this->SystemStructure();
    const phys::System &system = structure.system;
    const SystemOfCompoundMap &map = structure.map;
    const phys::ChannelFunctions cf = EngineChannelFunctions(*this, map, configure);

    // 1. Observed widths, ANCs and amplitudes -> Brune or Park amplitudes (Eqs. 2.1-2.3).
    for (std::size_t b = 0; b < system.num_blocks(); b++) {
      const int j = map.groupOfBlock[b];
      const JGroup *theJGroup = this->GetJGroup(j);
      const std::vector<int> &levels = map.levelOfBlockLevel[b];
      const std::size_t nc = theJGroup->NumChannels();
      const phys::ProjectLevels in = ProjectLevelsOfInputs(*this, inputs, map, b);
      const phys::BasisAmplitudes amp = phys::amplitudes_from_observed(
          system, b, phys::internal_channels(system, b), in, park ? phys::LevelBasis::Park : phys::LevelBasis::Brune, cf);
      for (std::size_t lb = 0; lb < levels.size(); lb++) {
        const int la = levels[lb];
        const ALevel *theLevel = theJGroup->GetLevel(la);
        for (std::size_t k = 0; k < nc; k++) {
          const int ch = (int)k + 1;
          const double input = in.value[lb * nc + k];
          // The --ignore-externals mask is a property of the input (external_capture.md section 3).
          state.SetZeroInputWidth(j, la, ch, input == 0.0);
          WarnInputColumn(*this, inputs, configure, j, la, ch, amp.notes[lb * nc + k], input);
        }
        WarnObservedOutOfReach(*this, inputs, configure, j, la, amp, in, lb, nc);
        const double parkNorm = amp.park_norm[lb];
        // The bound-state normalization of an external-capture final state is defined with Brune's reduced widths
        // gamma_obs^2 / J, so with Park's amplitudes it becomes 1 / sqrt(J + sum_c A_c gamma_obs^2) (see
        // CalcShiftFunctions).
        if (park) state.SetParkNorm(j, la, parkNorm);
        const std::vector<double> g(amp.gamma.begin() + lb * nc, amp.gamma.begin() + (lb + 1) * nc);
        for (std::size_t k = 0; k < nc; k++) inputs.SetGamma(j, la, (int)k + 1, g[k]);
        state.SetSqrtNF(j, la, SqrtNFFactor(*this, theJGroup, theLevel, park ? parkNorm : 1.0, g));
      }
    }

    // 2. External capture subtracted from the photon amplitudes (Eq. 3.2) of every group, while every level still
    // holds its Brune (Park) amplitudes: the final state enters external capture with its own pole's amplitude in
    // every basis (cnuc_split.md section 8).  The subtraction changes photon amplitudes only, which no other
    // level's external capture reads.
    for (std::size_t b = 0; b < system.num_blocks(); b++) {
      const int j = map.groupOfBlock[b];
      const JGroup *theJGroup = this->GetJGroup(j);
      const std::vector<int> &levels = map.levelOfBlockLevel[b];
      const std::size_t nc = theJGroup->NumChannels();
      std::vector<double> E, gamma;
      for (int la : levels) {
        E.push_back(inputs.E(j, la));
        for (std::size_t k = 0; k < nc; k++) gamma.push_back(inputs.Gamma(j, la, (int)k + 1));
      }
      const std::vector<phys::ExternalWidth> ext = ExternalWidthsOfInputs(*this, state, inputs, map, b, configure);
      std::vector<std::string> warnings;
      std::vector<std::size_t> exceeded;
      phys::subtract_external_capture(system, b, gamma, ext, warnings, &exceeded);
      for (std::size_t i : exceeded) {
        WarnExternalExceedsTotal(*this, inputs, configure, j, levels[i / nc], (int)(i % nc) + 1);
      }
      for (std::size_t lb = 0; lb < levels.size(); lb++)
        for (std::size_t k = 0; k < nc; k++) inputs.SetGamma(j, levels[lb], (int)k + 1, gamma[lb * nc + k]);
    }
    // 3. Without the Brune parametrization, the formal parameters (Eq. 4.1), group by group.
    for (std::size_t b = 0; b < system.num_blocks() && formal; b++) {
      const int j = map.groupOfBlock[b];
      const JGroup *theJGroup = this->GetJGroup(j);
      const std::vector<int> &levels = map.levelOfBlockLevel[b];
      const std::size_t nc = theJGroup->NumChannels();
      std::vector<double> E, gamma;
      for (int la : levels) {
        E.push_back(inputs.E(j, la));
        for (std::size_t k = 0; k < nc; k++) gamma.push_back(inputs.Gamma(j, la, (int)k + 1));
      }
      // A hard-sphere group (includeAllJPi) has channels but no levels: nothing to transform.
      if (!E.empty()) {
        std::vector<double> B(nc);
        for (std::size_t k = 0; k < nc; k++) B[k] = theJGroup->GetChannel((int)k + 1)->GetBoundaryCondition();
        phys::LevelSet f;
        try {
          // Each level keeps its position and receives the eigenvalue of the rank of its Brune energy, so a level
          // listed first (a background pole, say) keeps its own energy, widths and fix flags.  TransformOut uses the
          // same rank mapping in the other direction.
          f = phys::formal_from_brune(system, b, E, gamma, B, cf);
        } catch (const std::domain_error &) {
          // N a = E M a has a real solution only for a positive-definite M (Brune, PRC 66, 044611 (2002)); M is
          // 1 - sum_c gamma gamma' dS_c/dE, so widths too large for the channel radii make it indefinite.
          std::ostringstream msg;
          msg << "No standard R-matrix (constant boundary conditions) corresponds to the levels of J^pi = "
              << theJGroup->GetJ() << (theJGroup->GetPi() > 0 ? "+" : "-")
              << ": the observed widths are too large for the channel radii (the matrix 1 - sum gamma gamma' dS/dE "
              << "of the transformation is not positive definite).  Use the Brune parametrization (the default; "
              << "drop --no-brune) or smaller widths for these levels:";
          for (int la : levels)
            msg << std::endl
                << "    " << AZURELabel::Level(theJGroup, inputs.E(j, la), j, la);
          inputs.SetTransformError(msg.str());
          return false;
        }
        E.swap(f.E);
        gamma.swap(f.gamma);
        for (std::size_t lb = 0; lb < levels.size(); lb++) inputs.SetE(j, levels[lb], E[lb]);
      }
      for (std::size_t lb = 0; lb < levels.size(); lb++)
        for (std::size_t k = 0; k < nc; k++) inputs.SetGamma(j, levels[lb], (int)k + 1, gamma[lb * nc + k]);
    }
  } catch (const std::exception &e) {
    inputs.SetTransformError(std::string("Parameter transformation: ") + e.what());
    return false;
  }
  return true;
}

/*!
 * Performs the final parameter transformations from formal to physical parameters.
 */

ObservedLevelTable CompoundSystem::TransformOut(const Config &configure, LevelState &state, const LevelInputs &inputs) const {
  ObservedLevelTable table = this->EmptyObservedTable();
  this->TransformToBrune(configure, Messages::kPrint, state, inputs, table);
  this->ObservedWidths(configure, state, table);
  return table;
}

/*!
 * The Brune (Park) values of the fit parameters, E and gamma, with N_f^-1/2 formed from them, in `table`: in the
 * formal basis by physics/ brune_from_formal, otherwise a copy (N_f^-1/2 the fit's).  kPrint: print the warnings
 * for levels that do not convert.
 */

void CompoundSystem::TransformToBrune(const Config &configure, Messages messages, LevelState &state, const LevelInputs &inputs,
                            ObservedLevelTable &table) const {
  // Park (2021): the fit parameters already are the observed reduced width
  // amplitudes, so Gamma = 2 P gamma^2 with no 1 + sum gamma^2 dS/dE factor.
  // Refresh J and the bound-state normalization for the current parameters.
  const bool park = (configure.paramMask & Config::USE_PARK_FORMALISM);
  if (park) {
    this->CalcShiftFunctions(configure, state);
    for (int j = 1; j <= this->NumJGroups(); j++)
      for (int la = 1; la <= this->GetJGroup(j)->NumLevels(); la++) {
        const ALevel *theLevel = this->GetJGroup(j)->GetLevel(la);
        if (messages == Messages::kPrint && theLevel->IsInRMatrix() && state.ParkNorm(j, la) <= 0.0) WarnParkNormNotPositive(*this, state, inputs, configure, j, la);
      }
  }
  // Every level first as a copy of its fit values; the formal basis then converts the levels of the R-matrix.
  for (int j = 1; j <= this->NumJGroups(); j++)
    for (int la = 1; la <= this->GetJGroup(j)->NumLevels(); la++) {
      const ALevel *theLevel = this->GetJGroup(j)->GetLevel(la);
      ObservedLevel &out = table.Level(j, la);
      out.E = state.E(j, la);
      out.iterations = 0;
      out.sqrtNF = state.SqrtNF(j, la);
      for (int ch = 1; ch <= this->GetJGroup(j)->NumChannels(); ch++) out.gamma[ch - 1] = state.Gamma(j, la, ch);
    }
  if (!(configure.paramMask & Config::USE_BRUNE_FORMALISM)) {
    // Formal -> Brune by the fixed-point iteration of physics/ brune_from_formal (level_parameters.md section 4), on
    // the fit values.  A level not in the R-matrix keeps its values.
    const phys::FormalToBruneOptions options;  // the engine's 1e-6 MeV and 1000 iterations
    const SystemOfCompound &structure = this->SystemStructure();
    const phys::ChannelFunctions cf = EngineChannelFunctions(*this, structure.map, configure);
    for (std::size_t b = 0; b < structure.system.num_blocks(); b++) {
      const int j = structure.map.groupOfBlock[b];
      const JGroup *theJGroup = this->GetJGroup(j);
      const std::vector<int> &levels = structure.map.levelOfBlockLevel[b];
      const std::size_t nc = theJGroup->NumChannels();
      std::vector<double> E, gamma, B(nc);
      for (int la : levels) {
        E.push_back(state.E(j, la));
        for (std::size_t k = 0; k < nc; k++) gamma.push_back(state.Gamma(j, la, (int)k + 1));
      }
      for (std::size_t k = 0; k < nc; k++) B[k] = theJGroup->GetChannel((int)k + 1)->GetBoundaryCondition();
      const phys::LevelSet converted = phys::brune_from_formal(structure.system, b, E, gamma, B, cf, options);
      for (std::size_t lb = 0; lb < levels.size(); lb++) {
        const int la = levels[lb];
        if (messages == Messages::kPrint) WarnFormalToBrune(*this, inputs, configure, j, la, converted, lb, nc, options);
        ObservedLevel &out = table.Level(j, la);
        out.E = converted.E[lb];
        out.iterations = converted.iterations[lb];
        const std::vector<double> g(converted.gamma.begin() + lb * nc, converted.gamma.begin() + (lb + 1) * nc);
        for (std::size_t k = 0; k < nc; k++) out.gamma[k] = g[k];
        out.sqrtNF = SqrtNFFactor(*this, theJGroup, theJGroup->GetLevel(la), 1.0, g);
      }
    }
  }
}

/*!
 * The observed quantities of every level (cnuc_split.md section 8) from its Brune (Park) values in `table`
 * (TransformToBrune) and, with external capture, the external part of its photon amplitudes: physics/
 * observed_levels.  Sets the external amplitude, W and the reported value with its kind in `table`.
 */

void CompoundSystem::ObservedWidths(const Config &configure, const LevelState &state, ObservedLevelTable &table) const {
  const bool park = (configure.paramMask & Config::USE_PARK_FORMALISM);
  const bool external = (configure.paramMask & Config::USE_EXTERNAL_CAPTURE);
  const bool ignoreZero = (configure.paramMask & Config::IGNORE_ZERO_WIDTHS);
  const SystemOfCompound &all = this->SystemAllGroups();
  const phys::ChannelFunctions cf = EngineChannelFunctions(*this, all.map, configure);
  for (std::size_t b = 0; b < all.system.num_blocks(); b++) {
    const int j = all.map.groupOfBlock[b];
    const JGroup *theJGroup = this->GetJGroup(j);
    const std::vector<int> &levels = all.map.levelOfBlockLevel[b];
    const std::size_t nc = theJGroup->NumChannels();
    std::vector<double> E, gamma;
    std::vector<unsigned char> rwa;
    std::vector<phys::ExternalWidth> ext(levels.size() * nc);
    for (std::size_t lb = 0; lb < levels.size(); lb++) {
      const ALevel *theLevel = theJGroup->GetLevel(levels[lb]);
      ObservedLevel &out = table.Level(j, levels[lb]);
      E.push_back(out.E);
      for (std::size_t k = 0; k < nc; k++) {
        const int ch = (int)k + 1;
        const AChannel *theChannel = theJGroup->GetChannel(ch);
        gamma.push_back(out.gamma[k]);
        rwa.push_back(theLevel->GammaIsRWA(ch) ? 1 : 0);
        complex externalWidth(0.0, 0.0);
        if ((theChannel->GetRadType() == 'M' || theChannel->GetRadType() == 'E') && theLevel->IsInRMatrix() &&
            external && !(ignoreZero && state.ZeroInputWidth(j, levels[lb], ch))) {
          externalWidth = ExternalReducedWidthOfTable(*this, table, j, levels[lb], theChannel, configure);
          ext[lb * nc + k].present = true;
          ext[lb * nc + k].value = externalWidth;
        }
        out.external[k] = externalWidth;
      }
    }
    const phys::ObservedLevels obs =
        phys::observed_levels(all.system, b, phys::internal_channels(all.system, b), E, gamma, ext, rwa,
                              park ? phys::LevelBasis::Park : phys::LevelBasis::Brune, cf);
    for (std::size_t lb = 0; lb < levels.size(); lb++)
      for (std::size_t k = 0; k < nc; k++) {
        ObservedLevel &out = table.Level(j, levels[lb]);
        out.W[k] = obs.W[lb * nc + k];
        out.value[k] = obs.value[lb * nc + k];
        out.kind[k] = (unsigned char)obs.kind[lb * nc + k];
      }
  }
}

/*!
 * The observed widths of a parameter vector, as parameters.out reports them (TransformOut without its messages):
 * particle and photon widths with the level's normalization and, with external capture, the external part of the
 * photon amplitude (physics/ observed_levels).  The fit state is set back to the vector afterwards, since
 * TransformToBrune in Park's basis recomputes the shift functions and the bound-state normalizations.
 */

std::optional<ObservedLevelTable> CompoundSystem::ObservedWidthsAt(const Config &configure, const vector_r &params, LevelState &state,
                                                       const LevelInputs &inputs) const {
  this->FillCompoundFromParams(params, state);
  ObservedLevelTable table = this->EmptyObservedTable();
  try {
    this->TransformToBrune(configure, Messages::kQuiet, state, inputs, table);
    this->ObservedWidths(configure, state, table);
  } catch (const std::exception &) {
    this->FillCompoundFromParams(params, state);
    return std::nullopt;
  }
  this->FillCompoundFromParams(params, state);
  return table;
}
