// Schema::Build, the units of a data set (BuildContext.h): the unit of the values, the columns of its data file and
// the factor on the model value (project_format3.md sections 3.3, 3.4.1, 7.2; data_kinds.md sections 4.1 and 4.6).
//
// A data set's numbers are converted to the engine's units where the points are made (SetupSpecs.h ColumnConversion,
// experiment/DataRows.h): the data are compared in MeV, deg and the quantity's engine unit.  A quantity that is a
// scale of an observable the engine computes (an S factor, a yield per unit charge) is a factor on the model value
// (SetupSpecs.h ValueScale; the tree kind Scaled), never a conversion of the data, so that a fitted energy shift moves
// the factor.  An extrapolation has no data: its model is given in the unit it writes (a constant factor).
#include <cmath>
#include <limits>
#include <sstream>

#include "azure2/experiment/BuildContext.h"
#include "azure2/experiment/Frames.h"
#include "azure2/physics/Constants.h"

namespace azure2 {
namespace experiment {
namespace build {

namespace {

std::string Number(double x) {
  std::ostringstream s;
  s.precision(10);
  s << x;
  return s.str();
}

// A yield (data_kinds.md section 4.6; docs/upgrade/energy_grid.md section 5): the engine computes reactions per
// incident particle, 1e-24 int_{E_out}^{E_0} <sigma>/epsilon dE, on an infinite target (E_out = 1 keV, eq. 3) or on a
// layer (int_{E_out}^{E_0} dE/epsilon = N; the target effect's layerYield, BuildEffects.cpp); per unit charge the model
// is that times the particles per unit charge, 1e-6 C / (q e) for 1/uC; counts are reactions times the efficiency.
// On a layer, energy windows are refused: outside them the integration falls back to a cross section at one energy.
bool YieldScale(Context &c, const ValidatedRecord &r, const Unit &unit, ReactionSpec &spec) {
  bool ok = true;
  const Value *t = r.Get("target");
  const ValidatedRecord *target = t ? c.Record(t->word, "target") : nullptr;
  const Value *kind = target ? target->Get("kind") : nullptr;
  if (kind && kind->word == "layer" && r.Has("fold_windows")) {
    c.Error(codes::kFieldNotApplicable, r, "fold_windows",
            "fold_windows on a yield on target " + target->name + " (kind layer): outside a window the engine takes the "
            "cross section at the depth-averaged energy, which is not a yield; drop fold_windows (the whole layer is "
            "integrated at every point)");
    ok = false;
  }
  const std::string dimension = unit.dimension;
  if (dimension != "per_charge" && dimension != "counts_per_charge") return ok;  // per particle, per sr, arbitrary
  const Value *q = r.Get("charge_state");
  if (!q) return false;  // the validator refused it (E_FIELD_REQUIRED)
  const int chargeState = static_cast<int>(q->integer);
  spec.scale.factor = ParticlesPerUnitCharge(unit, chargeState);
  spec.scale.text = "x " + Number(spec.scale.factor) + " beam particles per unit of charge of " + std::string(unit.token) +
                    " (charge state " + std::to_string(chargeState) + ", e = " + Number(physics::kElementaryCharge_C) + " C)";
  if (dimension == "counts_per_charge") {
    const Value *e = r.Get("efficiency");
    if (!e) return false;  // the validator refused it
    const Nuisance &n = e->nuisance;
    if (n.mode != NuisanceMode::Fixed || n.hasPrior) {
      c.Error(codes::kNotBuilt, r, "efficiency",
              "a fitted efficiency: it multiplies the model as the normalization divides the data, so the two cannot be "
              "told apart; fix the efficiency and give its uncertainty as the prior of the norm");
      return false;
    }
    spec.scale.factor *= n.value.value;
    spec.scale.text += " x efficiency " + n.value.text;
  }
  return ok;
}

}  // namespace

bool ValueUnit(Context &c, const ValidatedRecord &r, ReactionSpec &spec, ColumnConversion *columns) {
  const Value *u = r.Get("unit");
  if (!u || !u->unit || !r.quantity) return true;
  const Unit &unit = *u->unit;
  if (r.quantity->name == "yield") return YieldScale(c, r, unit, spec);
  // Every engine unit of the table is the base of its dimension (Units.cpp), and a quantity's other units (1 or % of
  // a ratio of reactions, arbitrary) are bases or scales of theirs: the factor to the base is the factor to the engine.
  if (unit.factor == 1.0) return true;
  if (columns) {
    if (r.Has("covariance")) {
      c.Error(codes::kNotBuilt, r, "covariance",
              std::string("values in ") + unit.token + " with a covariance file: the file's absolute entries are compared "
              "with the values in the engine unit (" + r.quantity->engineUnit + "), and data sets that share it share them; "
              "write the values in " + r.quantity->engineUnit);
      return false;
    }
    columns->value = unit.factor;
    return true;
  }
  // an extrapolation: the model in the written unit
  spec.scale.factor = Rescale(1.0, 1.0 / unit.factor);
  spec.scale.text = "x " + Number(spec.scale.factor) + ": the values in " + unit.token + " (engine unit " + r.quantity->engineUnit + ")";
  return true;
}

bool DataColumns(Context &c, const ValidatedRecord &r, const Columns &cols, const std::string &secondRole, DatasetSpec &d) {
  const bool beta = r.quantity && r.quantity->name == "beta_delayed_spectrum";
  bool ok = !cols.fromFile && cols.roles.size() >= 4;
  const ColumnRole *energy = ok ? &cols.roles[0] : nullptr;
  if (ok) {
    const std::string first = energy->role;
    ok = beta ? first == "E_decay" : (first == "E_lab" || first == "E_cm");
  }
  for (std::size_t i = 1; ok && i < cols.roles.size(); i++) {
    const ColumnRole &role = cols.roles[i];
    const std::string expected = i == 1 ? secondRole : i == 2 ? "value" : i == 3 ? "error" : i == 4 ? "window_lo" : i == 5 ? "window_hi" : "_";
    if (role.role != expected && !(i >= 4 && role.role == "_") && !(i == 1 && role.role == "_" && expected == "_")) ok = false;
  }
  if (!ok) {
    c.Error(codes::kNotBuilt, r, "columns",
            std::string("this engine reads the columns ") + (beta ? "E_decay" : "E_lab|E_cm") + " " + secondRole +
                " value error[ window_lo window_hi], in this order ('from_file' and other orders: the format-3 reading "
                "of data files, step 3d)");
    return false;
  }
  ColumnConversion &conv = d.columns;
  conv.energy = energy->unit ? energy->unit->factor : 1.0;
  if (energy->role == "E_cm") {
    if (cols.Find("window_lo")) {
      c.Error(codes::kNotBuilt, r, "columns",
              "an E_cm column with window_lo, window_hi: the windows of a beam-profile kernel are reconstructed lab "
              "energies; give the energies as E_lab");
      ok = false;
    }
    const PairInfo *p = c.Pair(d.reaction.entranceKey);
    if (p) {
      conv.energyCM = true;
      conv.lightMass = p->lightMass;
      conv.heavyMass = p->heavyMass;
      conv.relativistic = c.opt_.relativistic;
    }
  }
  const ColumnRole *angle = cols.AngleRole();
  if (angle && angle->unit) conv.angle = angle->unit->factor;
  return ok;
}

void EnergyWindow(const ValidatedRecord &r, DatasetSpec &d) {
  const Value *w = r.Get("energy_window");
  if (!w) {
    d.minE = 0.0;
    d.maxE = std::numeric_limits<double>::max();
    return;
  }
  d.minE = ToBase(w->lo.value, w->unit);
  d.maxE = ToBase(w->hi.value, w->unit);
  const ColumnConversion &conv = d.columns;
  if (!conv.energyCM) return;
  Partition entrance;
  entrance.m1 = conv.lightMass;
  entrance.m2 = conv.heavyMass;
  d.minE = BeamLabFromCM(entrance, d.minE, conv.relativistic);
  d.maxE = BeamLabFromCM(entrance, d.maxE, conv.relativistic);
}

}  // namespace build
}  // namespace experiment
}  // namespace azure2
