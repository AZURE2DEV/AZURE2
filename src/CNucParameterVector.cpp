// CNuc: the parameter vector of the level block, read and written along the layout (LevelParameterLayout): the
// vector seeded from the levels' input values, a physical vector taken as input values, the observed vector of a
// table, and the change of amplitude basis between Brune and Park.
#include <cmath>
#include <string>
#include <vector>
#include "CNuc.h"
#include "AZUREParams.h"
#include "CNucExternalCapture.h"
#include "CNucSystem.h"
#include "Config.h"
#include "Minuit2/MnUserParameters.h"
#include "CNucLevelParameters.h"

namespace phys = azure2::physics;

/*!
 * Fills the CNuc object from the parameter array.
 */

void CNuc::FillCompoundFromParamsPhysical(const vector_r &observed, LevelState &state, LevelInputs &inputs) const {
  const vector_r &p = gammaDecays_->Empty() ? observed : gammaDecays_->WithPartialWidths(observed);
  if (!gammaDecays_->Empty()) inputs.SetGammaDecayNumbers(gammaDecays_->NumbersAt(observed));
  for (const LevelParameterLayout::Level &slots : this->ParameterLayout().Levels()) {
    const int j = slots.group, la = slots.level;
    const ALevel *level = this->GetJGroup(j)->GetLevel(la);
    inputs.SetE(j, la, p[slots.energy]);
    state.SetInputEnergy(j, la, p[slots.energy]);
    std::vector<double> g;
    for (int ch = 1; ch <= (int)slots.amplitude.size(); ch++) {
      inputs.SetGamma(j, la, ch, p[slots.amplitude[ch - 1]]);
      g.push_back(p[slots.amplitude[ch - 1]]);
    }
    state.SetSqrtNF(j, la, SqrtNFFactor(*this, this->GetJGroup(j), level, 1.0, g));
  }
}

/*!
 * Fills the Minuit parameter array from initial values in the CNuc object.
 */

void CNuc::FillMnParams(AZUREParams &params, const LevelState &state, const LevelInputs &inputs,
                        const Config *config) const {
  ROOT::Minuit2::MnUserParameters &p = params.GetMinuitParams();
  for (const LevelParameterLayout::Slot &slot : this->ParameterLayout().Slots()) {
    const int j = slot.group, la = slot.level, ch = slot.channel;
    const ALevel *level = this->GetJGroup(j)->GetLevel(la);
    const std::string varname = LevelParameterLayout::Name(slot);
    if (slot.kind == LevelParameterLayout::Slot::Kind::kEnergy) {
      const double energy = inputs.E(j, la);
      params.Add(varname, energy, 0.1 * energy, ParameterRecord::LevelEnergy(j, la, slot.levelNumber));
      bool isUnbound = false;
      for (int ir = 1; ir <= this->NumPairs(); ir++) {
        const PPair *pair = this->GetPair(ir);
        if (pair->GetNucleus() != this->GetJGroup(j)->GetNucleus()) continue;   // thresholds of this nucleus only
        if (pair->GetPType() == 0 &&
            energy > (pair->GetSepE() + pair->GetExE())) isUnbound = true;
      }
      // Jakub's Fix: it can create problems with the fit with pyazr and brick
      if (!isUnbound) p.Fix(varname);
      if (level->EnergyFixed() && !p.Parameter(p.Index(varname)).IsFixed()) p.Fix(varname);
      continue;
    }
    // Parameter settings will be applied by ParameterLimitsManager during fit
    if (gammaDecays_->AddParameter(params, inputs.GammaDecayNumbers(), j, la, ch, slot.levelNumber))
      continue;  // a gamma_decay record's number
    const double gamma = inputs.Gamma(j, la, ch);
    params.Add(varname, gamma, 0.1 * gamma, ParameterRecord::Width(j, la, ch, slot.levelNumber));
    // A width given as 0 is not fitted, in every basis: the input decides, not the transformed amplitude, which
    // in the formal basis is 0 or rounding by the order of the sums (external_capture.md section 9).  Without a
    // transformation the input is the parameter.
    if (inputs.TransformedIn() ? state.ZeroInputWidth(j, la, ch) : gamma == 0.0) p.Fix(varname);
    if (level->ChannelFixed(ch) && !p.Parameter(p.Index(varname)).IsFixed()) p.Fix(varname);
    // Wigner-limit bounds: GetWignerLimit is gamma^2_W = 3 hbar^2/(2 mu a^2) (MeV), the limit on the reduced width
    // squared; the fit parameter is the amplitude (MeV^1/2), so the bound is its square root (gamma^2_W itself
    // would let theta^2 = gamma^2/gamma^2_W run past 1 by a factor gamma^2_W).
    if (config && (config->paramMask & Config::USE_WIGNER_LIMITS)) {
      const AChannel *channel = this->GetJGroup(j)->GetChannel(ch);
      double wignerLimit = channel->GetWignerLimit();
      if (channel->GetRadType() == 'P' && wignerLimit > 0.0 && std::isfinite(wignerLimit)) {
        double gammaLimit = sqrt(wignerLimit);
        p.SetLimits(varname, -gammaLimit, gammaLimit);
      }
    }
  }
  // The isotopic abundances follow, written by WorkingCopy::FillMnParams (AbundanceSamples::AddParameters).
}


vector_r CNuc::GetTransformParams(const LevelState &state, const ObservedLevelTable &table) const {
  // The observed values TransformOut computed (physics/ observed_levels, cnuc_split.md section 8): Gamma (eV),
  // C (fm^-1/2), mu, Q, B and the Brune amplitude of a gammaIsRWA column, signed, as a <levels> line holds them.
  const LevelParameterLayout &layout = this->ParameterLayout();
  vector_r params(layout.Size());
  for (int i = 0; i < layout.Size(); i++) {
    const LevelParameterLayout::Slot &slot = layout.Slots()[i];
    const ObservedLevel &level = table.Level(slot.group, slot.level);
    params[i] =
        slot.kind == LevelParameterLayout::Slot::Kind::kEnergy
            ? level.E
            : gammaDecays_->ValueOr(*this, state, slot.group, slot.level, slot.channel, level.value[slot.channel - 1]);
  }
  // The level block only: the caller appends the abundance block of the fit's vector (AbundanceSamples::AppendTo).
  return params;
}

/*!
 * Converts a parameter set read from a file written in the other alternative
 * parametrization.  Brune -> Park multiplies every amplitude of a level by
 * sqrt(J), J = 1/(1 + sum_c gamma~^2 dS_c/dE); Park -> Brune divides by
 * sqrt(J), J = 1 - sum_c gamma^2 dS_c/dE.  Photon and beta amplitudes rescale
 * with the level like the particle ones.  Levels outside the R matrix
 * (external-capture final states) are the same in both and are left alone.
 * The positions are the layout's (ParameterLayout).
 */

bool CNuc::ConvertAmplitudeBasis(ROOT::Minuit2::MnUserParameters &p, int from, int to, const Config &configure) {
  if (from == to) return true;
  const bool toPark = (to == 2);
  if (!((from == 1 && to == 2) || (from == 2 && to == 1))) return false;
  // physics/ park_from_brune, brune_from_park on each R-matrix block (level_parameters.md Eq. 2.3); the errors
  // scale with the amplitudes, by J^1/2 or J^-1/2.
  const CNucSystem &structure = this->SystemStructure();
  const phys::ChannelFunctions cf = ChannelFunctionsFromCNuc(*this, structure.map, configure);
  for (std::size_t b = 0; b < structure.system.num_blocks(); b++) {
    const int j = structure.map.groupOfBlock[b];
    const std::vector<int> &levels = structure.map.levelOfBlockLevel[b];
    const int nc = this->Group(j)->NumChannels();
    std::vector<double> E, gamma, J;
    for (int la : levels) {
      const LevelParameterLayout::Level &slots = this->ParameterLayout().Of(j, la);
      E.push_back(p.Value(slots.energy));
      for (int ch = 1; ch <= nc; ch++) gamma.push_back(p.Value(slots.amplitude[ch - 1]));
    }
    std::vector<double> out;
    try {
      out = toPark ? phys::park_from_brune(structure.system, b, E, gamma, cf, &J)
                   : phys::brune_from_park(structure.system, b, E, gamma, cf, &J);
    } catch (const std::domain_error &error) {
      configure.outStream << "**ERROR: " << error.what() << std::endl
                          << "  No Brune (or standard) parameters correspond to these widths." << std::endl;
      return false;
    }
    for (std::size_t lb = 0; lb < levels.size(); lb++) {
      const LevelParameterLayout::Level &slots = this->ParameterLayout().Of(j, levels[lb]);
      const double factor = toPark ? std::sqrt(J[lb]) : 1.0 / std::sqrt(J[lb]);
      for (int ch = 1; ch <= nc; ch++) {
        const int a = slots.amplitude[ch - 1];
        p.SetValue(a, out[lb * nc + ch - 1]);
        p.SetError(a, p.Error(a) * factor);
      }
    }
  }
  return true;
}
