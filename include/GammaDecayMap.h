#ifndef GAMMADECAYMAP_H
#define GAMMADECAYMAP_H

#include <complex>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "CNucSystem.h"
#include "Constants.h"
#include "azure2/experiment/GammaDecay.h"

class AZUREParams;
class CNuc;
class Config;
class ParamIndexMap;
struct GradAccum;
namespace azure2 {
namespace io {
class LevelRow;
}
}  // namespace azure2

/*!
 * The radiative decay of a level written as Gamma_gamma, b_f and delta_f (a gamma_decay record,
 * docs/upgrade/gamma_parametrisation.md) between the fit's parameter vector and the level parameters the engine has:
 * the record's numbers sit in the vector at the positions of the level's photon amplitudes (one number per photon
 * channel), and every fill (CNuc::FillCompoundFromParams) replaces them by the amplitudes of the partial widths
 * they give (Eq. 2.3), converted by the engine's own conversion (physics/ amplitudes_from_observed) with the level's
 * current energy and particle amplitudes; a channel with external capture has its gamma^ext at the current parameters
 * subtracted as TransformIn does (section 2.7).  The engine's physics is not involved: it receives amplitudes as
 * before.
 *
 * Set-up state, read-only after CNuc::Fill and shared by the copies of the compound (CNuc::Clone).  What a fill
 * computes is the compound's (CNuc::GammaDecayFillOf), which the gradient's pull-back, the priors and the report read.
 */

/// What a fill leaves in the compound for one record (CNuc::GammaDecayFillOf): its numbers and, per slot (photon
/// channel), the parts of the engine's amplitude gamma_int = sign(g) (g^2 - (Im gamma^ext)^2)^1/2 - Re gamma^ext.
struct GammaDecayFill {
  std::vector<double> q;                         ///< the record's numbers, slot order (Gamma_gamma in eV)
  std::vector<double> total;                     ///< g_k, the amplitude of the record's observed width (MeV^1/2)
  std::vector<std::complex<double>> external;    ///< gamma^ext_k (MeV^1/2); 0 without external part
  std::vector<unsigned char> floor;              ///< 1: g_k^2 < (Im gamma^ext_k)^2, gamma_int = -Re gamma^ext_k
  /// d gamma^ext_k / d gamma_lambda c per channel of the level's group, and d gamma^ext_k / d gamma~_f d per channel
  /// of the final state's group (N_f's dependence included); empty without external part.
  std::vector<std::vector<std::complex<double>>> dLevel, dFinal;
};

/// One record placed in the compound.
struct GammaDecayLevel {
  azure2::experiment::GammaDecaySpec spec;
  enum class Part { TotalWidth, Branching, Mixing };
  struct Slot {
    Part part = Part::TotalWidth;
    int final = -1;  ///< index into spec.finals (Branching, Mixing)
  };
  std::vector<int> photonChannels;  ///< the group's photon channels (1-based), ascending: the record's slots
  std::vector<Slot> slots;          ///< what the number in each slot is, in the same order
  int energyIndex = 0;              ///< full-vector index of the level's energy (its channels follow)
  int levelNumber = 0;              ///< the level counted over all groups (names of the vector, from 1)
  std::size_t block = 0;            ///< the level's block of the structure System
  std::vector<int> channelOfBlockChannel;  ///< the group channel of each block channel
  /// The final state of a slot's external part (ExternalCaptureFinalState; group 0: none), the full-vector index of
  /// its energy and the number of its channels, which follow the energy in the vector.
  struct External {
    int group = 0, level = 0, energyIndex = -1, channels = 0;
  };
  std::vector<External> external;  ///< per slot
  bool HasExternal() const {
    for (const External &e : external)
      if (e.group) return true;
    return false;
  }
};

class GammaDecaySet {
 public:
  /*!
   * The <gammaDecays> section of the project (as ModelData::Read: records validated by the format-3 table, built
   * against <levels>), placed in the compound, with the partial widths of every record written as the level's input
   * widths (eV) for TransformIn.  Refused with a message (nullptr): the reader's diagnostics, the formal basis or
   * --no-transform, a photon channel that is a static moment, a level outside the R-matrix.  An empty set without the
   * section.  The set keeps `configure` (the run's, which outlives the compound) for the external-capture integrals.
   */
  static std::shared_ptr<const GammaDecaySet> Read(CNuc &compound, const Config &configure,
                                                   const std::vector<azure2::io::LevelRow> &rows);

  bool Empty() const { return levels_.empty(); }
  const std::vector<GammaDecayLevel> &Levels() const { return levels_; }
  /// The record of level (j, la), or -1; with the slot of channel ch in it (or -1 when ch is no photon channel).
  int RecordOf(int j, int la) const;
  int SlotOf(int record, int ch) const;
  bool OwnsSlot(int j, int la, int ch) const;

  /// CNuc::FillMnParams: the record's parameter at photon channel ch of level (j, la), if (j, la, ch) is a record
  /// slot (true), else nothing (false); its value from `inputs` (the compound's input numbers, InputNumbers' layout).
  bool AddParameter(AZUREParams &params, const std::vector<std::vector<double>> &inputs, int j, int la, int ch,
                    int levelNumber) const;
  /// The numbers of every record as the file gives them (per record, slot order; Gamma_gamma in eV): the compound's
  /// input numbers until FillCompoundFromParamsPhysical sets others.
  std::vector<std::vector<double>> FileNumbers() const;
  /// The numbers of every record at its slots of a full vector (the fit's or the observed one: the record's slots
  /// hold the same numbers in both, gamma_parametrisation.md section 2.8).
  std::vector<std::vector<double>> NumbersAt(const vector_r &p) const;
  /// The vector of observed parameters (CNuc::GetTransformParams' form) with every record's numbers replaced by the
  /// partial widths (eV) they give (PartialWidths): the levels' observed widths, which CNuc::
  /// FillCompoundFromParamsPhysical writes into the levels.  The one place a record slot of that vector is read.
  vector_r WithPartialWidths(const vector_r &observed) const;
  /// CNuc::FillCompoundFromParams: the full vector p with the record slots replaced by the engine's amplitudes; what
  /// the fill computed stored in the compound (GammaDecayFill).
  vector_r ToEngine(CNuc &compound, const vector_r &p) const;
  /// CNuc::GetTransformParams: the record's value of the last fill at a record slot, else `otherwise`.
  double ValueOr(const CNuc &compound, int j, int la, int ch, double otherwise) const;
  /// Partial widths (eV, signed, one per photon channel in slot order) of a record from its values q (slot order);
  /// with dGdq, also dGamma_k/dq_i (row-major k * m + i); excess: sum of the explicit branchings minus 1 (> 0: the
  /// remainder is 0 and the others are scaled, section 2.6).
  static std::vector<double> PartialWidths(const GammaDecayLevel &level, const std::vector<double> &q,
                                           std::vector<double> *dGdq = nullptr, double *excess = nullptr);

  /*!
   * The gradient's pull-back (Eq. 2.5): the cotangents of the record levels' photon amplitudes in `accum` (of the
   * compound filled at the point) turned into those of the record's numbers, and the level's energy and particle
   * amplitudes given their share; added to gradFull at the indices of `map`.  GradAccum::Scatter calls it and leaves
   * the record slots' own cotangents out.
   */
  void PullBack(CNuc &compound, const GradAccum &accum, const ParamIndexMap &map, vector_r &gradFull) const;
  /// The full-vector indices of everything the photon amplitudes of record r depend on: the level's energy and
  /// channels (the record's numbers, its particle amplitudes) and the channels of the final states of its external
  /// parts.  The gradient asks for the photon amplitudes' cotangents when one of them is free.
  std::vector<int> InputIndices(const ParamIndexMap &map, int r) const;

  /// The priors of the records and the remainder's wall as residuals (r, and dr/dq as (full index, value) pairs).
  struct PenaltyRow {
    double r = 0.0;
    std::vector<std::pair<int, double>> d;
  };
  std::vector<PenaltyRow> PenaltyRows(const CNuc &compound) const;
  /// Their chi^2, sum r^2.
  double Penalty(const CNuc &compound) const;
  /// d chi^2 / d q, as (full index, value) pairs (GradAccum::direct).
  void AddPenaltyGradient(const CNuc &compound, std::vector<std::pair<int, double>> &grad) const;

  /// The block of parameters.out: per record Gamma_gamma, every b_f, delta_f and the partial widths; a note for a
  /// branching or remainder at 0 and for a sum of branchings above 1.
  std::string Report(const CNuc &compound) const;

 private:
  struct Point;
  /// The map of record level L at energy E (MeV), its numbers q, the level's amplitudes `particle` (per group channel)
  /// and the final states' energies and amplitudes per slot (Point in GammaDecayMap.cpp); with slopes, the
  /// derivatives of gamma^ext.
  Point Map(CNuc &compound, const GammaDecayLevel &L, double E, const std::vector<double> &q,
            const std::vector<double> &particle, const std::vector<std::pair<double, std::vector<double>>> &finals,
            bool slopes) const;

  bool park_ = false, gslCoulomb_ = false, rmc_ = false;
  const Config *configure_ = nullptr;  ///< the run's (external-capture integrals: its per-thread scratch)
  std::shared_ptr<const CNucSystem> structure_;  ///< the structure System (CNucSystemContent::Structure)
  std::vector<GammaDecayLevel> levels_;
};

#endif
