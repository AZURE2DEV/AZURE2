#ifndef ALEVEL_H
#define ALEVEL_H

#include "Constants.h"
#include <vector>


namespace azure2 {
namespace io {
class LevelRow;
}
}  // namespace azure2

/// An AZURE level object.

/*!
 * An R-matrix level  represents a specific eigenstate of the compound nucleus.
 *
 * A level carries its set-up data only: IDs, fix flags, the external-capture final-state data, and the values of its
 * <levels> lines (FileE, FileGamma, FileZeroWidth), which CNuc::Fill reads to group the levels and to make the
 * owner's LevelInputs.  Nothing of a level is written after set-up.  The input values in the engine's basis (what
 * CNuc::TransformIn makes of the file's) are a LevelInputs of the owner (the run, the session, ...); the values a
 * calculation reads at a parameter vector (the fit energy and amplitudes, the shift functions, Park's J, the
 * external-capture amplitudes and N_f^1/2) are a LevelState of the owner.
 *
 * The observed values of the reports are not the level's: CNuc::TransformOut returns them as an
 * ObservedLevelTable.
 *
 * Channel numbers are 1-based and index the channels of the owning JGroup, so
 * the same index means the same channel for every level in the group.
 */

class ALevel {
 public:
  /// Build from a line of the nuclear input file.
  explicit ALevel(const azure2::io::LevelRow &);
  /// Build with only an energy; widths are added afterwards (PadChannels, SetGamma).
  ALevel(double);
  /// Is this level part of the A-/R-Matrix calculation?
  /// False for a level that exists only as an external-capture final state,
  /// such as a subthreshold bound state.
  bool IsInRMatrix() const;
  /// Is the level energy held fixed during fitting?
  bool EnergyFixed() const;
  /// Is this channel's reduced width amplitude held fixed during fitting?
  bool ChannelFixed(int) const;
  /// Does this channel's <levels> gamma field hold a reduced-width amplitude
  /// directly, rather than a physical width? (The 33-field THM line's flag.)
  bool GammaIsRWA(int) const;
  /// Was this channel's width in the <levels> lines exactly zero (or the level lists no line for it)?  The start of
  /// LevelInputs::ZeroInputWidth, which TransformIn sets again from the values it transforms.
  bool FileZeroWidth(int) const;
  /// Is the level a final state for external capture?
  bool IsECLevel() const;
  int NumNFIntegrals() const;
  /// 1-based pair number this level is the bound state of; 0 if it is not one.
  int GetECPairNum() const;
  /// Bit mask of the external-capture multipolarities feeding this level.
  unsigned char GetECMultMask() const;

  // -- the <levels> values --------------------------------------------------
  /// Level energy as the <levels> line gives it (MeV, E_x).
  double FileE() const;
  /// The channel's width as the <levels> lines give it (partial width in eV, ANC, or an amplitude: GammaIsRWA); 0
  /// for a channel the level does not list.
  double FileGamma(int) const;

  // -- external capture ------------------------------------------------------
  /// Channel integral in the denominator of \f$N_f^{1/2}\f$,
  /// \f$\int_a^\infty [W_c(kr)/W_c(ka_c)]^2\f$.
  double GetNFIntegral(int) const;
  /// The energy (MeV, E_x) with which it enters external capture: its Brune energy at setup (cnuc_split.md section 8).
  double GetECEnergy() const;
  /// Factor converting this channel's reduced width amplitude to an ANC.
  double GetECConversionFactor(int) const;

  /// Put the width of an input file line at its J-group channel number (1-based), padding any channel
  /// before it with a zero, fixed width (CNuc::Fill): a level may list its channels in another order than its
  /// group, or leave one out.
  void SetGammaFromLine(int channelNum, const azure2::io::LevelRow &nucLine);
  /// Pad the width vectors with zero, fixed widths up to n channels.
  void PadChannels(int n);
  void AddNFIntegral(double);
  void SetECEnergy(double);
  void AddECConversionFactor(double);
  /// Mark the level as an external-capture final state of a pair, with its
  /// multipolarity mask.
  void SetECParams(int, unsigned char);
  /// The level IDs (column 9 of <levels>) of the lines read into this level, in file order without repeats: the
  /// names a <modelData> datum uses for it (docs/upgrade/project_format3.md F13).
  const std::vector<int> &IDs() const { return ids_; }
  void AddID(int id);

 private:
  std::vector<int> ids_;
  bool isinrmatrix_;
  bool energyfixed_;
  bool isECLevel_;
  int ecPairNum_;
  unsigned char ecMultMask_;
  double level_e_;
  std::vector<bool> channelfixed_;
  /// Per-channel input convention: true = gamma column was a reduced width
  /// amplitude (MeV^(1/2)), false = physical partial width/ANC (legacy default).
  std::vector<bool> gamma_is_rwa_;
  /// See FileZeroWidth.
  std::vector<bool> zero_input_;
  vector_r gammas_;
  double ec_energy_ = 0.0;
  vector_r nf_integrals_;
  vector_r ec_conv_factors_;
};

#endif
