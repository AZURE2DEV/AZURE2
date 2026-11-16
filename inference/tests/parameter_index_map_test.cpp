// inference_parameter_index_map: the index map of the full parameter vector (azure2/inference/ParameterIndexMap.h),
// made by hand for a model with a parameter of every kind, in the order of the parameter file: the level block, the
// abundances, the normalizations, the three kinds of energy shift, screening, density, Doppler, thickness, the THM
// coherent background, the apparatus parameters.  Every lookup, the description of every entry and the packed <->
// full maps with two entries fixed are written out; an append that registers the wrong lookup or a Pack that
// counts a fixed entry fails here.
#include "azure2/inference/ParameterIndexMap.h"

#include <vector>

#include "check.h"

using namespace azure2::inference;

int main() {
  ParamIndexMap map;
  // Two J groups: one level of two channels (its second amplitude a gamma_decay number), one level of one channel.
  map.SizeLevelBlock({1, 1}, {2, 1});
  check("before any append: no index", map.NumFull() == 0 && map.EIndexFlat(1, 1) == -1 && map.GammaIndexFlat(1, 1, 2) == -1 &&
                                           map.NormIndex(1) == -1 && map.EnergyShiftFnIndex(1, 0) == -1);
  bool returned = map.AppendLevelEnergy(1, 1) == 0;
  returned = returned && map.AppendAmplitude(1, 1, 1, false) == 1;
  returned = returned && map.AppendAmplitude(1, 1, 2, true) == 2;
  returned = returned && map.AppendLevelEnergy(2, 1) == 3;
  returned = returned && map.AppendAmplitude(2, 1, 1, false) == 4;
  returned = returned && map.AppendAbundance(1, 1, 2) == 5;  // ordinal 1: nucleus 2 of sample 1
  returned = returned && map.AppendNorm(2) == 6;             // segment 1 has no free norm
  returned = returned && map.AppendEnergyShift(1) == 7;
  returned = returned && map.AppendEnergyShift(2) == 8;
  returned = returned && map.AppendEnergyShiftSqrt(1) == 9;
  returned = returned && map.AppendEnergyShiftSqrt(2) == 10;
  returned = returned && map.AppendEnergyShiftFn(2, 0) == 11;
  returned = returned && map.AppendEnergyShiftFn(2, 1) == 12;
  returned = returned && map.AppendScreening(1) == 13;
  returned = returned && map.AppendDensity(1) == 14;
  returned = returned && map.AppendDoppler(1) == 15;
  returned = returned && map.AppendThickness(1, 2) == 16;
  returned = returned && map.AppendThmCoherent(1) == 17;
  returned = returned && map.AppendThmCoherent(2) == 18;
  returned = returned && map.AppendApparatus(1) == 19;
  check("each append returns the next full index", returned);
  check("twenty entries", map.NumFull() == 20);

  check("level block", map.EIndexFlat(1, 1) == 0 && map.GammaIndexFlat(1, 1, 1) == 1 && map.GammaIndexFlat(1, 1, 2) == 2 &&
                           map.EIndexFlat(2, 1) == 3 && map.GammaIndexFlat(2, 1, 1) == 4);
  check("abundance by ordinal", map.AbundanceIndex(1) == 5 && map.AbundanceIndex(2) == -1);
  check("norm by segment", map.NormIndex(2) == 6 && map.NormIndex(1) == -1);
  check("shifts by segment", map.EnergyShiftIndex(1) == 7 && map.EnergyShiftIndex(2) == 8 && map.EnergyShiftIndex(3) == -1 &&
                                 map.EnergyShiftSqrtIndex(1) == 9 && map.EnergyShiftSqrtIndex(2) == 10);
  check("shift function by (segment, j)", map.EnergyShiftFnIndex(2, 0) == 11 && map.EnergyShiftFnIndex(2, 1) == 12 &&
                                              map.EnergyShiftFnIndex(1, 0) == -1 && map.EnergyShiftFnIndex(2, 2) == -1);
  check("target effects by their number", map.ScreeningIndex(1) == 13 && map.DensityIndex(1) == 14 && map.DopplerIndex(1) == 15 &&
                                              map.ThicknessIndex(1) == 16 && map.ScreeningIndex(2) == -1 &&
                                              map.DensityIndex(0) == -1 && map.DopplerIndex(2) == -1 && map.ThicknessIndex(2) == -1);
  check("apparatus by its number", map.ApparatusIndex(1) == 19 && map.ApparatusIndex(2) == -1);

  const ParamKind kinds[] = {ParamKind::LevelEnergy, ParamKind::Gamma, ParamKind::RadiativeDecay, ParamKind::LevelEnergy,
                             ParamKind::Gamma, ParamKind::Abundance, ParamKind::Norm, ParamKind::EnergyShift,
                             ParamKind::EnergyShift, ParamKind::EnergyShiftSqrt, ParamKind::EnergyShiftSqrt,
                             ParamKind::EnergyShiftFn, ParamKind::EnergyShiftFn, ParamKind::Screening, ParamKind::Density,
                             ParamKind::Doppler, ParamKind::Thickness, ParamKind::ThmCoherent, ParamKind::ThmCoherent,
                             ParamKind::Apparatus};
  bool kindsOk = true;
  for (int i = 0; i < map.NumFull(); i++) kindsOk = kindsOk && map.Desc(i).kind == kinds[i];
  check("the kind of every entry", kindsOk);
  check("description: level entries", map.Desc(2).jGroup == 1 && map.Desc(2).level == 1 && map.Desc(2).channel == 2 &&
                                          map.Desc(2).segment == -1 && map.Desc(3).jGroup == 2 && map.Desc(3).channel == -1);
  check("description: abundance", map.Desc(5).sample == 1 && map.Desc(5).nucleus == 2 && map.Desc(5).jGroup == -1);
  check("description: segment entries", map.Desc(6).segment == 2 && map.Desc(8).segment == 2 && map.Desc(12).segment == 2 &&
                                            map.Desc(12).fnTerm == 1 && map.Desc(11).fnTerm == 0 && map.Desc(6).fnTerm == -1);
  check("description: target effects", map.Desc(13).screening == 1 && map.Desc(14).density == 1 && map.Desc(15).doppler == 1 &&
                                           map.Desc(16).thickness == 1 && map.Desc(16).segment == 2 && map.Desc(13).density == -1);
  check("description: cbkg and apparatus", map.Desc(17).channel == 1 && map.Desc(18).channel == 2 && map.Desc(19).apparatus == 1);

  // Packing: entries 3 (an energy) and 7 (a shift) fixed; a mask shorter than the vector leaves the rest free.
  std::vector<bool> fixed(8, false);
  fixed[3] = fixed[7] = true;
  map.Pack(fixed);
  check("eighteen free", map.NumPacked() == 18);
  bool packed = true;
  for (int f = 0, a = 0; f < map.NumFull(); f++) {
    const bool isFixed = f == 3 || f == 7;
    packed = packed && map.Desc(f).fixed == isFixed && map.FullToPacked(f) == (isFixed ? -1 : a);
    if (!isFixed) packed = packed && map.PackedToFull(a++) == f;
  }
  check("packed <-> full, fixed entries left out", packed);
  // Packing again with another mask replaces the first.
  map.Pack({});
  check("no mask: every entry free", map.NumPacked() == 20 && map.FullToPacked(3) == 3 && !map.Desc(7).fixed);
  return result();
}
