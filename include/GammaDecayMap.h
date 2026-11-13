#ifndef GAMMADECAYMAP_H
#define GAMMADECAYMAP_H

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
 * current energy and particle amplitudes.  The engine's physics is not involved: it receives amplitudes as before.
 *
 * Set-up state, read-only after CNuc::Fill and shared by the copies of the compound (CNuc::Clone).  The record's
 * values of the last fill are the compound's (CNuc::GammaDecayValues), which the gradient's pull-back and the priors
 * read.
 */

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
};

class GammaDecaySet {
 public:
  /*!
   * The <gammaDecays> section of the project (as ModelData::Read: records validated by the format-3 table, built
   * against <levels>), placed in the compound, with the partial widths of every record written as the level's input
   * widths (eV) for TransformIn.  Refused with a message (nullptr): the reader's diagnostics, the formal basis or
   * --no-transform, a record channel with external capture in the run, a photon channel that is a static moment, a
   * level outside the R-matrix.  An empty set without the section.
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
  /// slot (true), else nothing (false).
  bool AddParameter(AZUREParams &params, int j, int la, int ch, int levelNumber) const;
  /// CNuc::FillCompoundFromParams: the full vector p with the record slots replaced by the engine's amplitudes; the
  /// record's values stored in the compound.
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
  bool park_ = false, gslCoulomb_ = false, rmc_ = false;
  std::shared_ptr<const CNucSystem> structure_;  ///< the structure System (CNucSystemContent::Structure)
  std::vector<GammaDecayLevel> levels_;
};

#endif
