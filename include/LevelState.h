#ifndef LEVEL_STATE_H
#define LEVEL_STATE_H

#include <cstdint>
#include <vector>

#include "GammaDecayFill.h"
#include "azure2/physics/LevelPairShift.h"

/*!
 * The levels at one parameter vector and what follows from it: the values an evaluation writes and every
 * calculation reads.  Per level the fit energy and reduced-width amplitudes, the external-capture amplitudes and
 * N_f^1/2 with which it enters external capture as a final state, the shift functions at its energy (with the energy
 * they were computed at), and in Park's basis J = 1 - sum_c gamma_c^2 dS_c/dE with the first and second derivatives.
 * Beside the levels: what the fill left for each gamma_decay record (GammaDecayFill) and, in the formal basis with
 * --ignore-externals, the projector of each capture mask (CompoundSystem::FormalCaptureMaskSlot).
 *
 * Indexed like the compound's levels: J-group j and level la, both 1-based, and the group's channels ch, 1-based.  The
 * per-channel values of a group lie level after level, (la - 1) n_c + (ch - 1), as in physics/ JpiBlock.
 *
 * Owned by whoever evaluates (the program's run, the session, a fit's pooled evaluation, a walker, the rate, ...: a
 * value, copied from its caller's); CompoundSystem's fills (FillCompoundFromParams with GammaDecaySet::ToEngine,
 * CalcShiftFunctions, TransformIn, SetFitLevelParameter) write the one they are given, readers get it const.  Its
 * shape (levels, channels, capture masks) is set by the compound's set-up (CompoundSystem::Fill, Initialize).
 */
class LevelState {
 public:
  LevelState() = default;
  /// numLevels[j - 1] levels of numChannels[j - 1] channels in J-group j, each with the values of a level not yet
  /// filled: E and amplitudes 0, N_f^1/2 and J 1, no shift functions computed.
  LevelState(const std::vector<int> &numLevels, const std::vector<int> &numChannels);

  // -- the fit values ---------------------------------------------------------------------------------------------
  /// Fit level energy, MeV (E_x).
  double E(int j, int la) const { return groups_[j - 1].E[la - 1]; }
  /// Fit reduced-width amplitude, MeV^1/2.
  double Gamma(int j, int la, int ch) const { return groups_[j - 1].gamma[At(j, la, ch)]; }

  // -- external capture -------------------------------------------------------------------------------------------
  /// The amplitude with which the level enters external capture as a final state (MeV^1/2): its Brune (Park)
  /// amplitude in every basis (cnuc_split.md section 8); the fit amplitude except in the formal basis.
  double ECGamma(int j, int la, int ch) const { return groups_[j - 1].ecGamma[At(j, la, ch)]; }
  /// The N_f^1/2 normalization of the level's bound-state wave function.
  double SqrtNF(int j, int la) const { return groups_[j - 1].sqrtNF[la - 1]; }
  /// Was the channel's input width exactly zero (a zero in the <levels> line, a channel the level does not list, or a
  /// zero in the physical vector handed to CompoundSystem::TransformIn)?  With IGNORE_ZERO_WIDTHS (--ignore-externals, the
  /// default) such a radiative channel has no external part: neither in the observed width (TransformIn,
  /// TransformOut) nor in the resonant channel capture through the level (AMatrixFunc).  A property of the input,
  /// written where the input is read (CompoundSystem::Fill, CompoundSystem::TransformIn), not a test on the fit amplitude.
  bool ZeroInputWidth(int j, int la, int ch) const { return groups_[j - 1].zero[At(j, la, ch)] != 0; }

  // -- shift functions at the level energy ------------------------------------------------------------------------
  /// S_c(E_lambda) of a channel.
  double Shift(int j, int la, int ch) const { return groups_[j - 1].shift[At(j, la, ch)]; }
  /// dS_c/dE and d^2S_c/dE^2 at E_lambda, MeV^-1 and MeV^-2: in Park's basis and for a level coincident with another
  /// of its J-group (azure2::physics::levels_coincide; the limits of the divided differences need them), 0 otherwise;
  /// the second derivative once it was asked for (SecondDerivativesValidAt).
  double ShiftDerivative(int j, int la, int ch) const { return groups_[j - 1].d1[At(j, la, ch)]; }
  double ShiftSecondDerivative(int j, int la, int ch) const { return groups_[j - 1].d2[At(j, la, ch)]; }
  /// Park (2021) Eq. (28): J = 1 - sum_c gamma_c^2 (dS_c/dE) at E_lambda; 1 outside Park's basis.
  double ParkNorm(int j, int la) const { return groups_[j - 1].parkNorm[la - 1]; }
  /// Were the shift functions computed at exactly this energy?
  bool ShiftsValidAt(int j, int la, double energy) const { return groups_[j - 1].shiftEnergy[la - 1] == energy; }
  /// Were the first derivatives computed at exactly this energy?
  bool FirstDerivativesValidAt(int j, int la, double energy) const {
    return groups_[j - 1].d1Energy[la - 1] == energy;
  }
  /// Were the second derivatives computed at exactly this energy?
  bool SecondDerivativesValidAt(int j, int la, double energy) const {
    return groups_[j - 1].d2Energy[la - 1] == energy;
  }

  /*!
   * A number that changes whenever the energy, an amplitude, a shift function, J or a second derivative of the level
   * is set, unique over all levels of the process: equal versions mean those values are equal (ChannelSpace keys its
   * memos on the state, the level and the version).
   */
  std::uint64_t Version(int j, int la) const { return groups_[j - 1].version[la - 1]; }

  // -- beside the levels -------------------------------------------------------------------------------------------
  /// Has a fill stored what it computed for gamma_decay record r (GammaDecaySet::ToEngine)?
  bool HasGammaDecayFill(int r) const { return r >= 0 && r < (int)gammaDecayFills_.size(); }
  /// What the last fill left for record r: its values (slot order) and the parts of its photon amplitudes.
  const GammaDecayFill &GammaDecayFillOf(int r) const { return gammaDecayFills_[r]; }
  /// The projector Pi = T^-1 P T of capture-mask slot `slot` (CompoundSystem::FormalCaptureMaskOf), row-major over the slot's
  /// levels; empty: Pi = P (no masked level, or every masked level decoupled).
  const std::vector<double> &CaptureProjector(int slot) const { return captureProjectors_[slot]; }

  // -- writing (CompoundSystem's fills) -------------------------------------------------------------------------------------
  void SetE(int j, int la, double energy);
  void SetGamma(int j, int la, int ch, double gamma);
  void SetECGamma(int j, int la, int ch, double gamma) { groups_[j - 1].ecGamma[At(j, la, ch)] = gamma; }
  void SetSqrtNF(int j, int la, double value) { groups_[j - 1].sqrtNF[la - 1] = value; }
  void SetZeroInputWidth(int j, int la, int ch, bool zero) { groups_[j - 1].zero[At(j, la, ch)] = zero ? 1 : 0; }
  void SetShift(int j, int la, int ch, double value);
  void SetShiftDerivative(int j, int la, int ch, double value);
  void SetShiftSecondDerivative(int j, int la, int ch, double value);
  void SetParkNorm(int j, int la, double value);
  void SetShiftsEnergy(int j, int la, double energy) { groups_[j - 1].shiftEnergy[la - 1] = energy; }
  void SetFirstDerivativesEnergy(int j, int la, double energy) { groups_[j - 1].d1Energy[la - 1] = energy; }
  void SetSecondDerivativesEnergy(int j, int la, double energy) { groups_[j - 1].d2Energy[la - 1] = energy; }
  void SetGammaDecayFill(int r, GammaDecayFill fill);
  /// `slots` capture-mask slots, every projector empty (CompoundSystem's set-up).
  void ShapeCaptureProjectors(std::size_t slots) { captureProjectors_.assign(slots, std::vector<double>()); }
  void SetCaptureProjector(int slot, std::vector<double> projector) { captureProjectors_[slot] = std::move(projector); }

 private:
  struct Group {
    int numChannels = 0;
    std::vector<double> E, sqrtNF, parkNorm, shiftEnergy, d1Energy, d2Energy;  ///< per level
    std::vector<double> gamma, ecGamma, shift, d1, d2;               ///< per level and channel
    std::vector<unsigned char> zero;                                  ///< per level and channel
    std::vector<std::uint64_t> version;                               ///< per level
  };
  std::size_t At(int j, int la, int ch) const {
    return (std::size_t)(la - 1) * groups_[j - 1].numChannels + (std::size_t)(ch - 1);
  }
  void Touch(int j, int la);
  std::vector<Group> groups_;
  std::vector<GammaDecayFill> gammaDecayFills_;          ///< per gamma_decay record, once a fill has stored it
  std::vector<std::vector<double>> captureProjectors_;  ///< per capture-mask slot
};

/// The shift function of channel ch between levels la and lb of J-group j at this state, as the level matrices hold
/// it (azure2::physics::LevelPairShift: E, S, S', S'' of both levels).
inline azure2::physics::LevelPairShift<double> PairShift(const LevelState &state, int j, int la, int lb, int ch) {
  auto at = [&](int l) {
    return azure2::physics::ShiftAtLevel<double>{state.E(j, l), state.Shift(j, l, ch), state.ShiftDerivative(j, l, ch),
                                                 state.ShiftSecondDerivative(j, l, ch)};
  };
  return azure2::physics::level_pair_shift(at(la), at(lb));
}

#endif
