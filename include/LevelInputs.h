#ifndef LEVEL_INPUTS_H
#define LEVEL_INPUTS_H

#include <string>
#include <vector>

/*!
 * The levels' starting values in the engine's basis: what a transformation into the fit's basis reads and writes.
 * Per level the input energy (MeV, E_x) and per channel the input width or amplitude, with the flag that the input
 * width was exactly zero (the --ignore-externals mask, external_capture.md section 3); the input numbers of the
 * gamma_decay records; whether the transformation ran and why it failed.
 *
 * Until the transformation has run the widths are the file's: partial widths in eV, ANCs, amplitudes where the
 * <levels> line says so (a gamma_decay record's level holds the partial widths of its numbers).  After it they are
 * the reduced-width amplitudes of the fit's basis (MeV^1/2), and in the formal basis the energies are the formal ones.
 *
 * Indexed like the compound's levels: J-group j and level la, both 1-based, and the group's channels ch, 1-based; the
 * per-channel values of a group lie level after level, (la - 1) n_c + (ch - 1), as in LevelState.
 *
 * Held by the compound nucleus (CNuc::Inputs); made by CNuc::Fill from the <levels> rows (ALevel's file values) and
 * the gamma_decay records, rewritten by CNuc::TransformIn and CNuc::FillCompoundFromParamsPhysical.  FillMnParams and
 * the reports read it.
 */
class LevelInputs {
 public:
  LevelInputs() = default;
  /// numLevels[j - 1] levels of numChannels[j - 1] channels in J-group j, every value 0, no flag set.
  LevelInputs(const std::vector<int> &numLevels, const std::vector<int> &numChannels);

  /// Input level energy, MeV (E_x).
  double E(int j, int la) const { return groups_[j - 1].E[la - 1]; }
  /// Input width or amplitude of a channel (see the class comment for its unit).
  double Gamma(int j, int la, int ch) const { return groups_[j - 1].gamma[At(j, la, ch)]; }
  /// Was the channel's input width exactly zero (a zero in the <levels> line, a channel the level does not list, or a
  /// zero in the physical vector handed to TransformIn)?  With IGNORE_ZERO_WIDTHS (--ignore-externals, the default)
  /// such a radiative channel has no external part: neither in the observed width (TransformIn / TransformOut) nor in
  /// the resonant channel capture through the level (AMatrixFunc).  A property of the input, fixed when the width is
  /// read, not a test on the working amplitude.
  bool ZeroInputWidth(int j, int la, int ch) const { return groups_[j - 1].zero[At(j, la, ch)] != 0; }

  /// The input numbers of the gamma_decay records (per record, slot order; Gamma_gamma in eV), which FillMnParams puts
  /// into the parameter vector as it puts a level's input widths: the file's, or the last physical vector's.
  const std::vector<std::vector<double>> &GammaDecayNumbers() const { return gammaDecayNumbers_; }

  /// Has the transformation into the fit's basis run on these values?  Until it has, the widths are the file's.
  bool TransformedIn() const { return transformedIn_; }
  /// Why the last transformation failed (empty if it did not).
  const std::string &TransformError() const { return transformError_; }

  // -- writing (CNuc::Fill, TransformIn, FillCompoundFromParamsPhysical) --------------------------------------------
  void SetE(int j, int la, double energy) { groups_[j - 1].E[la - 1] = energy; }
  void SetGamma(int j, int la, int ch, double gamma) { groups_[j - 1].gamma[At(j, la, ch)] = gamma; }
  void SetZeroInputWidth(int j, int la, int ch, bool zero) { groups_[j - 1].zero[At(j, la, ch)] = zero ? 1 : 0; }
  void SetGammaDecayNumbers(std::vector<std::vector<double>> numbers) { gammaDecayNumbers_ = std::move(numbers); }
  /// The transformation starts (TransformedIn true, no error).
  void StartTransformation() {
    transformedIn_ = true;
    transformError_.clear();
  }
  void SetTransformError(std::string why) { transformError_ = std::move(why); }

 private:
  struct Group {
    int numChannels = 0;
    std::vector<double> E;              ///< per level
    std::vector<double> gamma;          ///< per level and channel
    std::vector<unsigned char> zero;    ///< per level and channel
  };
  std::size_t At(int j, int la, int ch) const {
    return (std::size_t)(la - 1) * groups_[j - 1].numChannels + (std::size_t)(ch - 1);
  }
  std::vector<Group> groups_;
  std::vector<std::vector<double>> gammaDecayNumbers_;
  bool transformedIn_ = false;
  std::string transformError_;
};

#endif
