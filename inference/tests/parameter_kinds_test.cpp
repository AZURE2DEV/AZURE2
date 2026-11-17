// inference_parameter_kinds: the table of parameter kinds and the record of one parameter
// (azure2/inference/ParameterKind.h).  A row moved out of the order of the codes, a key used twice or a factory that
// leaves an owner field unset fails here.
#include "azure2/inference/ParameterKind.h"

#include <set>
#include <string>

#include "check.h"

using namespace azure2::inference;

int main() {
  const std::vector<ParamKindRow> &table = ParamKindTable();
  check("fourteen kinds", table.size() == 14);
  bool ordered = true, filled = true;
  std::set<std::string> keys;
  for (std::size_t i = 0; i < table.size(); i++) {
    ordered = ordered && (int)table[i].kind == (int)i;
    filled = filled && table[i].key && *table[i].key && table[i].unit && *table[i].unit && table[i].owner &&
             *table[i].owner && table[i].doc && *table[i].doc;
    keys.insert(table[i].key);
  }
  check("row i is the kind of code i", ordered);
  check("every row has key, unit, owner fields and doc", filled);
  check("keys are unique", keys.size() == table.size());
  check("KindRow finds the row", std::string(KindRow(ParamKind::Norm).key) == "norm" &&
                                     std::string(KindRow(ParamKind::RadiativeDecay).key) == "radiative" &&
                                     std::string(KindRow(ParamKind::LevelEnergy).unit) == "MeV" &&
                                     std::string(KindRow(ParamKind::Gamma).unit) == "MeV^1/2");
  // The codes are part of the parameter records pyazr reads: they never change.
  check("the codes", (int)ParamKind::LevelEnergy == 0 && (int)ParamKind::Gamma == 1 && (int)ParamKind::Norm == 2 &&
                         (int)ParamKind::EnergyShift == 3 && (int)ParamKind::EnergyShiftSqrt == 4 &&
                         (int)ParamKind::EnergyShiftFn == 5 && (int)ParamKind::Abundance == 6 &&
                         (int)ParamKind::Screening == 7 && (int)ParamKind::Density == 8 && (int)ParamKind::Doppler == 9 &&
                         (int)ParamKind::Thickness == 10 && (int)ParamKind::ThmCoherent == 11 &&
                         (int)ParamKind::Apparatus == 12 && (int)ParamKind::RadiativeDecay == 13);
  // A level parameter and a normalization are freed by a parameter file's nonzero value; the others by the project.
  check("own priors: every kind but energy, width and cbkg",
        !KindRow(ParamKind::LevelEnergy).ownPrior && !KindRow(ParamKind::Gamma).ownPrior &&
            !KindRow(ParamKind::ThmCoherent).ownPrior && KindRow(ParamKind::Norm).ownPrior &&
            KindRow(ParamKind::Density).ownPrior && KindRow(ParamKind::RadiativeDecay).ownPrior);
  check("segment fields: norm, the three shifts, thickness",
        KindRow(ParamKind::Norm).segmentField && KindRow(ParamKind::EnergyShift).segmentField &&
            KindRow(ParamKind::EnergyShiftSqrt).segmentField && KindRow(ParamKind::EnergyShiftFn).segmentField &&
            KindRow(ParamKind::Thickness).segmentField && !KindRow(ParamKind::Abundance).segmentField &&
            !KindRow(ParamKind::Screening).segmentField && !KindRow(ParamKind::Gamma).segmentField);
  check("level slots: energy, width, radiative", IsLevelSlot(ParamKind::LevelEnergy) && IsLevelSlot(ParamKind::Gamma) &&
                                                     IsLevelSlot(ParamKind::RadiativeDecay) &&
                                                     !IsLevelSlot(ParamKind::Norm) && !IsLevelSlot(ParamKind::ThmCoherent));

  // The factories: the fields of the kind are set, the others stay -1 or empty.
  const ParameterRecord e = ParameterRecord::LevelEnergy(2, 3, 7);
  check("LevelEnergy", e.kind == ParamKind::LevelEnergy && e.jGroup == 2 && e.level == 3 && e.levelNumber == 7 &&
                           e.channel == -1 && e.segment == -1);
  const ParameterRecord w = ParameterRecord::Width(2, 3, 4, 7);
  check("Width", w.kind == ParamKind::Gamma && w.jGroup == 2 && w.level == 3 && w.channel == 4 && w.levelNumber == 7);
  const ParameterRecord r = ParameterRecord::Radiative(1, 2, 5, 3, "decay1", "branching", "16O(0+)");
  check("Radiative", r.kind == ParamKind::RadiativeDecay && r.channel == 5 && r.decayRecord == "decay1" &&
                         r.decayPart == "branching" && r.finalState == "16O(0+)" && std::string(r.Row().key) == "radiative");
  const ParameterRecord n = ParameterRecord::OfSegment(ParamKind::Norm, 4, 12);
  check("OfSegment", n.kind == ParamKind::Norm && n.segment == 4 && n.segmentKey == 12 && n.jGroup == -1);
  const ParameterRecord f = ParameterRecord::ShiftFnTerm(4, 12, 2);
  check("ShiftFnTerm", f.kind == ParamKind::EnergyShiftFn && f.segment == 4 && f.fnTerm == 2);
  const ParameterRecord t = ParameterRecord::OfEffect(ParamKind::Screening, 3);
  check("OfEffect", t.kind == ParamKind::Screening && t.effect == 3 && t.segment == -1);
  const ParameterRecord a = ParameterRecord::Apparatus(2, "stopping_scale", "1", {{"material", "TiN"}, {"projectile", "1H"}});
  check("Apparatus", a.kind == ParamKind::Apparatus && a.effect == 2 && a.apparatus == "stopping_scale" &&
                         a.apparatusUnit == "1" && a.owner.size() == 2 && a.owner[1].second == "1H");
  const ParameterRecord b = ParameterRecord::Abundance(1, "natural", 2, 10, 22);
  check("Abundance", b.kind == ParamKind::Abundance && b.sample == 1 && b.sampleName == "natural" && b.nucleus == 2 &&
                         b.nucleusZ == 10 && b.nucleusA == 22);
  const ParameterRecord c = ParameterRecord::ThmCoherent("exp1", 1, 2, 3, 1);
  check("ThmCoherent", c.kind == ParamKind::ThmCoherent && c.thmExperiment == "exp1" && c.jGroup == 1 && c.channel == 2 &&
                           c.exitChannel == 3 && c.thmTerm == 1);
  return result();
}
