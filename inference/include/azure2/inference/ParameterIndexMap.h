#ifndef AZURE2_INFERENCE_PARAMETERINDEXMAP_H
#define AZURE2_INFERENCE_PARAMETERINDEXMAP_H

#include <map>
#include <utility>
#include <vector>

#include "azure2/inference/ParameterKind.h"

namespace azure2 {
namespace inference {

/// Description of one entry in the full (unfiltered) flat parameter vector.
struct ParamDesc {
  ParamKind kind;
  int jGroup;   ///< 1-based JGroup index   (LevelEnergy/Gamma only, else -1)
  int level;    ///< 1-based level index    (LevelEnergy/Gamma only, else -1)
  int channel;  ///< 1-based channel index  (Gamma; ThmCoherent: its 1-based index in that block; else -1)
  int segment;  ///< 1-based segment index  (Norm/EnergyShift/EnergyShiftSqrt only, else -1)
  bool fixed;   ///< whether this parameter is fixed (excluded from packed vector)
  int nucleus = -1;  ///< 1-based compound nucleus (Abundance only, else -1)
  int sample = -1;   ///< 1-based sample of the <abundances> block (Abundance only, else -1)
  int fnTerm = -1;  ///< j of a shift-function parameter aj (EnergyShiftFn only, else -1)
  int screening = -1;  ///< 1-based screening parameter, Screening only, else -1
  int density = -1;    ///< 1-based density parameter, Density only, else -1
  int doppler = -1;    ///< 1-based Doppler parameter, Doppler only, else -1
  int thickness = -1;  ///< 1-based thickness parameter, Thickness only, else -1
  int apparatus = -1;  ///< 1-based apparatus parameter, Apparatus only, else -1
};

/*!
 * \brief Bidirectional map between what a parameter belongs to and its place in the flat parameter vector.
 *
 * The flat vector layout is, in order (the order of the parameter file):
 *   - for each JGroup j, for each level la:  E_{j,la}, then gamma_{j,la,ch} for ch=1..NumChannels
 *   - the isotopic abundances
 *   - the norm parameters
 *   - the energy-shift parameters
 *   - the sqrt(E) energy-shift coefficients
 *   - the shift-function parameters a0.. of each segment declaring a "shiftfn" block, in segment order
 *   - the screening potentials, the fitted densities, the Doppler temperatures, the sample thicknesses
 *   - the THM coherent-background parameters (cbackground=; none without it)
 *   - the apparatus parameters
 *
 * "Full" indices run over every entry above.  The "packed" index is the position within the non-fixed subset
 * (the vector a sampler or an optimizer passes).
 *
 * The map does not know the model: its maker appends the entries in the vector's order, one call per parameter
 * (the engine's BuildParamIndexMap walks the compound system and the data set and calls the appends below), then
 * calls Pack with the fixed flags.  An append returns the full index it gave.
 */
class ParamIndexMap {
 public:
  /// Full index of the norm parameter for a 1-based segment, or -1 if that
  /// segment does not carry a (varying) norm parameter.
  int NormIndex(int segment) const { return Find(normIndex_, segment); }
  /// Full index of the energy-shift parameter for a 1-based segment, or -1.
  int EnergyShiftIndex(int segment) const { return Find(shiftIndex_, segment); }
  /// Full index of the sqrt(E) energy-shift coefficient for a 1-based segment, or -1.
  int EnergyShiftSqrtIndex(int segment) const { return Find(sqrtShiftIndex_, segment); }
  /// Full index of the abundance parameter with this 1-based ordinal, or -1.
  int AbundanceIndex(int ordinal) const { return Find(abundanceIndex_, ordinal); }
  /// Full index of parameter aj of the shift function declared on a 1-based segment, or -1.
  int EnergyShiftFnIndex(int segment, int j) const {
    auto it = fnShiftIndex_.find({segment, j});
    return (it == fnShiftIndex_.end()) ? -1 : it->second;
  }
  /// Full index of the 1-based screening parameter k, or -1.
  int ScreeningIndex(int k) const { return Find(screeningIndex_, k); }
  /// Full index of the 1-based density parameter k, or -1.
  int DensityIndex(int k) const { return Find(densityIndex_, k); }
  /// Full index of the 1-based Doppler-temperature parameter k, or -1.
  int DopplerIndex(int k) const { return Find(dopplerIndex_, k); }
  /// Full index of the 1-based apparatus parameter k, or -1.
  int ApparatusIndex(int k) const { return Find(apparatusIndex_, k); }
  /// Full index of the 1-based sample-thickness parameter k, or -1.
  int ThicknessIndex(int k) const { return Find(thicknessIndex_, k); }

  /// Number of entries in the full (unfiltered) parameter vector.
  int NumFull() const { return (int)desc_.size(); }
  /// Number of non-fixed (packed) parameters.
  int NumPacked() const { return numPacked_; }

  /// Description of full-index i.
  const ParamDesc &Desc(int fullIndex) const { return desc_[fullIndex]; }

  /// Packed index for a full index, or -1 if that parameter is fixed.
  int FullToPacked(int fullIndex) const { return fullToPacked_[fullIndex]; }
  /// Full index for a packed index.
  int PackedToFull(int packedIndex) const { return packedToFull_[packedIndex]; }

  /// Full index of E_{j,la} / gamma_{j,la,ch} (1-based) without a map lookup, or -1.
  int EIndexFlat(int jGroup, int level) const { return eFlat_[jGroup - 1][level - 1]; }
  int GammaIndexFlat(int jGroup, int level, int channel) const { return gammaFlat_[jGroup - 1][level - 1][channel - 1]; }

  // ---- Making the map: in the vector's order.
  /// The sizes of the level block: numLevels[j - 1] levels of numChannels[j - 1] channels in J-group j.  First.
  void SizeLevelBlock(const std::vector<int> &numLevels, const std::vector<int> &numChannels);
  /// The energy of level la of J-group j.
  int AppendLevelEnergy(int jGroup, int level);
  /// The amplitude of channel ch of level la of J-group j; `radiative`: the slot holds a number of a gamma_decay
  /// record (kind RadiativeDecay) instead of the amplitude itself.
  int AppendAmplitude(int jGroup, int level, int channel, bool radiative);
  /// The abundance with this 1-based ordinal, of nucleus `nucleus` in sample `sample`.
  int AppendAbundance(int ordinal, int sample, int nucleus);
  int AppendNorm(int segment);
  int AppendEnergyShift(int segment);
  int AppendEnergyShiftSqrt(int segment);
  int AppendEnergyShiftFn(int segment, int j);
  int AppendScreening(int k);
  int AppendDensity(int k);
  int AppendDoppler(int k);
  /// Thickness parameter k, of segment `segment`.
  int AppendThickness(int k, int segment);
  /// The k-th (1-based) parameter of the THM coherent-background block (no lookup: the block is read in order).
  int AppendThmCoherent(int k);
  int AppendApparatus(int k);
  /// Last: the fixed flag of every entry (by full index; a vector shorter than the map leaves the rest free) and
  /// the packed <-> full lookups.
  void Pack(const std::vector<bool> &fixed);

 private:
  static int Find(const std::map<int, int> &index, int key) {
    auto it = index.find(key);
    return (it == index.end()) ? -1 : it->second;
  }
  int Append(const ParamDesc &d) {
    desc_.push_back(d);
    return (int)desc_.size() - 1;
  }

  std::vector<std::vector<int>> eFlat_;                  ///< [j-1][la-1] -> full idx
  std::vector<std::vector<std::vector<int>>> gammaFlat_; ///< [j-1][la-1][ch-1] -> full idx
  std::vector<ParamDesc> desc_;                          ///< per full index
  std::map<int, int> normIndex_;                         ///< segment -> full idx
  std::map<int, int> shiftIndex_;                        ///< segment -> full idx
  std::map<int, int> sqrtShiftIndex_;                    ///< segment -> full idx
  std::map<int, int> abundanceIndex_;                    ///< abundance ordinal -> full idx
  std::map<std::pair<int, int>, int> fnShiftIndex_;      ///< (segment, j) -> full idx
  std::map<int, int> screeningIndex_;                    ///< screening parameter k -> full idx
  std::map<int, int> densityIndex_;                      ///< density parameter k -> full idx
  std::map<int, int> dopplerIndex_;                      ///< Doppler parameter k -> full idx
  std::map<int, int> thicknessIndex_;                    ///< thickness parameter k -> full idx
  std::map<int, int> apparatusIndex_;                    ///< apparatus parameter k -> full idx
  std::vector<int> fullToPacked_;
  std::vector<int> packedToFull_;
  int numPacked_ = 0;
};

}  // namespace inference
}  // namespace azure2

#endif  // AZURE2_INFERENCE_PARAMETERINDEXMAP_H
