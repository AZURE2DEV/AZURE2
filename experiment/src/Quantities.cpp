// The quantity table of format 3 (Quantities.h; docs/upgrade/project_format3.md section 7.2).  One function per
// quantity fills its row member by member; Quantities() lists them in schema order.  The validator reads nothing else
// about a quantity.  Rows follow the note's tables; where the note leaves a choice open the row says so in a comment.
//
// A new quantity is one more function here, named after it, and its entry in Quantities(), for example:
//
//   Quantity SFactor() {
//     Quantity q;
//     q.name = "sfactor";
//     q.title = "Astrophysical S factor";
//     q.doc = "S(E) = E sigma(E) exp(2 pi eta) with the engine's eta and atomic masses (data_kinds.md section 4.1).";
//     q.recordKinds = kDataKinds;
//     q.engineUnit = "MeV*b";
//     q.units = {AnyUnitOf("sfactor")};
//     q.energyRoles = {"E_cm", "E_lab"};
//     q.reaction = {R::kChargedEntrance, "charged entrance"};
//     q.fields = {Allows("multipole", kGammaExit)};
//     q.tree = "Scale_SF(E) o Sigma(r)";
//     q.exfor = {",SIG,,SFC", "PAR,SIG,,SFC", "EP/PAR,SIG,,SFC"};
//     q.refusals = {RefusedAtValidation("E_REACTION_NOT_ALLOWED", "neutral",
//                                       "the S factor is defined for charged particles (Gamow factor)")};
//     return q;
//   }
//
// Members a row does not set keep the defaults of Quantities.h.
#include "azure2/experiment/Quantities.h"

#include <utility>

namespace azure2 {
namespace experiment {

namespace {

using R = ReactionRule;

/// Every unit of a dimension (Units.cpp), allowed when `when` holds.
UnitRule AnyUnitOf(const std::string &dimension, Condition when = {}) { return {"", dimension, when}; }
/// One unit token, allowed when `when` holds.
UnitRule OneUnit(const std::string &unit, Condition when = {}) { return {unit, "", when}; }
/// A field the quantity allows (when `when` holds); a field that no QuantityField names is refused for the quantity.
QuantityField Allows(const std::string &field, Condition when = {}) { return {field, false, when}; }
/// A field the quantity requires (when `when` holds).
QuantityField Requires(const std::string &field, Condition when = {}) { return {field, true, when}; }

/// Norm modes allowed when `when` holds, with or without a prior.
NormRule NormModesWithPrior(Condition when, std::vector<std::string> modes, const std::string &message) {
  return {std::move(when), std::move(modes), true, message};
}
NormRule NormModesWithoutPrior(Condition when, std::vector<std::string> modes, const std::string &message) {
  return {std::move(when), std::move(modes), false, message};
}

/// A refusal and where it is checked: by Schema::Validate, by step 3's Build (needs <levels> or the data file), or by
/// the EXFOR importer.  `rule` is a rule name of ReactionRule or the condition in words.
Refusal RefusedAtValidation(const std::string &code, const std::string &rule, const std::string &message) {
  return {code, rule, message, "validate"};
}
Refusal RefusedAtBuild(const std::string &code, const std::string &rule, const std::string &message) {
  return {code, rule, message, "build"};
}
/// Checked by the engine's set-up, which evaluates what the table cannot (a stopping formula on a sub-point grid).
Refusal RefusedAtSetUp(const std::string &code, const std::string &rule, const std::string &message) {
  return {code, rule, message, "setup"};
}
Refusal RefusedByImporter(const std::string &code, const std::string &rule, const std::string &message) {
  return {code, rule, message, "importer"};
}

const Condition kNoThm = Condition::Absent("thm");
const Condition kGammaExit = Condition::Of(Condition::GammaExit);
const Condition kRatio = Condition::Of(Condition::RatioReaction);
const Condition kAngle = Condition::Of(Condition::AngleAxis);
const Condition kNoAngle = Condition::Of(Condition::AngleAxis, true);

const std::vector<std::string> kDataKinds = {"dataset", "extrapolation"};
const std::vector<std::string> kDatum = {"datum"};

// ---------------------------------------------------------------- data sets and extrapolations

Quantity Sigma() {
  Quantity q;
  q.name = "sigma";
  q.aliases = {"cross_section"};
  q.title = "Cross section";
  q.doc = "Angle-integrated cross section sigma(E) of the reaction; on a thin layer target the target-averaged sigma "
          "(data_kinds.md eq. 2).";
  q.recordKinds = kDataKinds;
  q.engineUnit = "b";
  q.units = {AnyUnitOf("cross_section"), OneUnit("1", kRatio), OneUnit("arbitrary")};
  q.energyRoles = {"E_lab", "E_cm", "E_decay"};
  q.windowRoles = true;
  q.reaction = {R::kSum | R::kRatio | R::kNoChargedElastic, "sums and ratios; no charged-particle elastic"};
  q.fields = {Allows("target", kNoThm), Allows("beam"), Allows("detector"), Allows("multipole", kGammaExit),
              Allows("secondary_gamma"), Allows("thm")};
  q.targetKinds = {"layer", "none"};
  q.tree = "Sigma(r) [amplitude filter]; Sum/Scale/Ratio from the reaction; Weighted if the target has isotopes; folds "
           "inside the leaf (M4)";
  q.exfor = {",SIG", "PAR,SIG", ",SIG,,TTA"};
  q.refusals = {
      RefusedAtValidation("E_REACTION_NOT_ALLOWED", "charged_elastic",
                          "the angle-integrated elastic cross section of charged particles is infinite: give "
                          "dsigma_domega or ratio_to_rutherford (LEXFOR p. S.2); the integral of the nuclear amplitude alone is "
                          "quantity sigma_nuclear"),
      RefusedAtValidation("E_FIELD_NOT_APPLICABLE", "multipole with a particle exit",
                          "a multipole filter selects photon multipoles"),
      RefusedAtValidation("E_TARGET_KIND", "target kind infinite",
                          "a thick (infinite) target gives a yield: use quantity yield")};
  return q;
}

Quantity SigmaNuclear() {
  Quantity q;
  q.name = "sigma_nuclear";
  q.title = "Nuclear elastic integral";
  q.doc = "Angle-integrated elastic cross section of the nuclear amplitude alone, pi/k^2 sum g |T|^2; for a charged pair "
          "it leaves out the Coulomb amplitude and its interference and is not a measurable cross section: for tests, "
          "diagnostics and comparison with phase-shift formulas.";
  q.recordKinds = kDataKinds;
  q.engineUnit = "b";
  q.units = {AnyUnitOf("cross_section"), OneUnit("1", kRatio)};
  q.energyRoles = {"E_lab", "E_cm"};
  // A neutral pair's elastic sigma is the same integral and keeps the one name 'sigma'.  Sums and ratios of such terms
  // as for sigma (a natural target's isotopes, tests/natural_target).
  q.reaction = {R::kSum | R::kRatio | R::kElastic | R::kChargedEntrance, "elastic, charged entrance; sums and ratios"};
  q.fields = {Allows("target"), Allows("beam"), Allows("detector")};
  q.targetKinds = {"layer", "none"};
  q.tree = "SigmaNuclear(r) = pi/k^2 sum g |T|^2 (GenMatrixFunc, angle-integrated, entrance = exit)";
  q.refusals = {RefusedAtValidation("E_REACTION_NOT_ALLOWED", "inelastic", "the nuclear elastic integral is an elastic quantity"),
                RefusedAtValidation("E_REACTION_NOT_ALLOWED", "neutral",
                                    "for a neutral pair the elastic cross section is quantity sigma")};
  q.datasetWarning = "sigma_nuclear leaves out the Coulomb amplitude and its interference: it is not a measurable cross "
                     "section (for tests, diagnostics and comparison with phase-shift formulas)";
  return q;
}

Quantity SigmaTotal() {
  Quantity q;
  q.name = "sigma_total";
  q.title = "Total cross section";
  q.doc = "Total cross section of a neutral projectile, from the optical theorem; on a target with isotopes the "
          "abundance-weighted sum of its nuclei's totals, one elastic term per nucleus (20Ne(n,n)20Ne + 22Ne(n,n)22Ne).";
  q.recordKinds = kDataKinds;
  q.engineUnit = "b";
  q.units = {AnyUnitOf("cross_section")};
  q.energyRoles = {"E_lab", "E_cm"};
  // A sum only on a target with isotopes (the set-up refuses it otherwise: BuildEffects.cpp NaturalTotals).
  q.reaction = {R::kSum | R::kElastic | R::kNeutralEntrance, "elastic, neutral entrance; on a target with isotopes a sum, one term per nucleus"};
  // Apparatus allowed as for sigma (the note's row lists only the reaction rule).
  q.fields = {Allows("target"), Allows("beam"), Allows("detector")};
  q.targetKinds = {"layer", "none"};
  q.tree = "SigmaTotal(a -> a) by the optical theorem (GM 103-119)";
  q.exfor = {",SIG (SF3 TOT)"};
  q.refusals = {RefusedAtValidation("E_REACTION_NOT_ALLOWED", "charged_entrance",
                                    "the total cross section of charged particles is infinite")};
  return q;
}

Quantity DsigmaDomega() {
  Quantity q;
  q.name = "dsigma_domega";
  q.aliases = {"dsdo"};
  q.title = "Differential cross section";
  q.doc = "dsigma/dOmega at the angle of the angle column; the value's frame is the angle's (theta_cm: c.m. value).";
  q.recordKinds = kDataKinds;
  q.engineUnit = "b/sr";
  // '1' with a ratio: section 3.5's example 12C(p,p)12C / 12C(p,p)12C @ 90 deg cm (correction to section 7.2).
  q.units = {AnyUnitOf("diff_cross_section"), OneUnit("1", kRatio)};
  q.energyRoles = {"E_lab", "E_cm"};
  q.angleRoles = {"theta_lab", "theta_cm"};
  q.angle = Quantity::RequiredAngle;
  q.windowRoles = true;  // a beam-profile kernel on differential data (tests/beam_profile_kernel)
  q.reaction = {R::kSum | R::kRatio | R::kNoTotal, "sums and ratios; no [total] exit"};
  q.fields = {Allows("target", kNoThm), Allows("beam"), Allows("detector"), Allows("secondary_gamma"), Allows("thm")};
  q.targetKinds = {"layer", "none"};
  q.tree = "dsdo_cm(r, theta); lab data converted at load";
  q.exfor = {",DA", "PAR,DA"};
  q.refusals = {RefusedAtValidation("E_REACTION_NOT_ALLOWED", "total",
                                    "[total] sums angle-integrated photon cross sections: use quantity sigma"),
                RefusedByImporter("E_MODE_NOT_ALLOWED", ",DA,,REL without norm fit",
                                  "relative data need a fitted normalization")};
  return q;
}

Quantity SFactor() {
  Quantity q;
  q.name = "sfactor";
  q.title = "Astrophysical S factor";
  q.doc = "S(E) = E sigma(E) exp(2 pi eta) with the engine's eta and atomic masses (data_kinds.md section 4.1).";
  q.recordKinds = kDataKinds;
  q.engineUnit = "MeV*b";
  q.units = {AnyUnitOf("sfactor")};
  q.energyRoles = {"E_cm", "E_lab"};
  q.reaction = {R::kChargedEntrance, "charged entrance"};
  q.fields = {Allows("multipole", kGammaExit)};
  q.tree = "Scale_SF(E) o Sigma(r)";
  q.exfor = {",SIG,,SFC", "PAR,SIG,,SFC", "EP/PAR,SIG,,SFC"};
  q.refusals = {RefusedAtValidation("E_REACTION_NOT_ALLOWED", "neutral",
                                    "the S factor is defined for charged particles (Gamow factor)")};
  return q;
}

Quantity ModifiedSFactor() {
  Quantity q;
  q.name = "modified_sfactor";
  q.title = "Modified S factor of heavy-ion fusion";
  q.doc = "S*(E) = E sigma(E) exp(2 pi eta + g E), g the field 'shape' (Patterson et al. 1969; data_kinds.md section 4.1).";
  q.recordKinds = kDataKinds;
  q.engineUnit = "MeV*b";
  q.units = {AnyUnitOf("sfactor")};
  q.energyRoles = {"E_cm", "E_lab"};
  q.reaction = {R::kChargedEntrance, "charged entrance"};
  q.fields = {Requires("shape"), Allows("multipole", kGammaExit)};
  q.tree = "Scale_SF*(E; g) o Sigma(r)";
  // no EXFOR code: S* is not an EXFOR quantity (its data are SIG, which the importer maps to sigma)
  q.refusals = {RefusedAtValidation("E_REACTION_NOT_ALLOWED", "neutral",
                                    "the modified S factor is defined for charged particles (Gamow factor)"),
                RefusedAtValidation("E_FIELD_REQUIRED", "no shape", "the constant g of exp(g E) is given, never guessed")};
  return q;
}

Quantity RatioToReference() {
  Quantity q;
  q.name = "ratio_to_reference";
  q.title = "Ratio to a reference cross section";
  q.doc = "sigma(E) of the reaction divided by a tabulated reference sigma_ref(E) (a monitor reaction or a standard); "
          "the table's energies are in the data set's energy role, frame and unit (data_kinds.md section 4.9).";
  q.recordKinds = kDataKinds;
  q.engineUnit = "1";
  q.units = {OneUnit("1")};
  q.energyRoles = {"E_lab", "E_cm"};
  q.reaction = {R::kNoChargedElastic, "no charged-particle elastic"};
  q.fields = {Requires("reference"), Allows("target"), Allows("beam"), Allows("detector")};
  q.targetKinds = {"layer", "none"};
  q.tree = "Ratio(Fold[sigma(r)], Fold[sigma_ref(table)])";
  q.exfor = {};  // EXFOR states the monitor in MONITOR/MONIT-REF with the data as SIG: the importer does not map it
  q.refusals = {RefusedAtValidation("E_FIELD_REQUIRED", "no reference", "the reference table is named, never assumed"),
                RefusedAtBuild("E_REFERENCE_TABLE", "table",
                               "the table cannot be read, is not increasing in energy, has a value <= 0, or does not "
                               "cover an energy of the data set or extrapolation")};
  return q;
}

Quantity RatioToRutherford() {
  Quantity q;
  q.name = "ratio_to_rutherford";
  q.aliases = {"rth"};
  q.title = "Ratio to Rutherford";
  q.doc = "dsigma/dOmega(c.m.) divided by the Rutherford cross section at the same c.m. angle.";
  q.recordKinds = kDataKinds;
  q.engineUnit = "1";
  q.units = {OneUnit("1")};
  q.energyRoles = {"E_lab", "E_cm"};
  q.angleRoles = {"theta_lab", "theta_cm"};
  q.angle = Quantity::RequiredAngle;
  q.reaction = {R::kElastic | R::kChargedEntrance | R::kNotIdentical, "elastic, charged, not identical"};
  q.fields = {Allows("target"), Allows("beam"), Allows("detector")};
  q.targetKinds = {"layer", "none"};
  q.tree = "Ratio(Fold[dsdo_cm(r)], Fold[rutherford(r)])";
  q.exfor = {",DA,,RTH", "PAR,DA,,RTH"};
  q.refusals = {
      RefusedAtValidation("E_REACTION_NOT_ALLOWED", "identical",
                          "has identical particles; the Coulomb reference is Mott-symmetrised: use quantity "
                          "'ratio_to_mott' (K2)"),
      RefusedAtValidation("E_REACTION_NOT_ALLOWED", "neutral", "has no Coulomb scattering: give dsigma_domega"),
      RefusedAtValidation("E_REACTION_NOT_ALLOWED", "inelastic",
                          "the ratio to Rutherford is not defined: give dsigma_domega")};
  return q;
}

Quantity RatioToMott() {
  Quantity q;
  q.name = "ratio_to_mott";
  q.title = "Ratio to Mott";
  q.doc = "dsigma/dOmega(c.m.) of identical particles divided by the Mott cross section.";
  q.recordKinds = kDataKinds;
  q.engineUnit = "1";
  q.units = {OneUnit("1")};
  q.energyRoles = {"E_lab", "E_cm"};
  q.angleRoles = {"theta_lab", "theta_cm"};
  q.angle = Quantity::RequiredAngle;
  q.reaction = {R::kElastic | R::kChargedEntrance | R::kIdentical, "elastic, identical particles"};
  q.fields = {Allows("target"), Allows("beam"), Allows("detector")};
  q.targetKinds = {"layer", "none"};
  q.tree = "Ratio(Fold[dsdo_cm(r)], Fold[mott(r)])";
  q.exfor = {",DA,,MOT"};
  q.refusals = {RefusedAtValidation("E_REACTION_NOT_ALLOWED", "non_identical",
                                    "has no identical particles: use quantity 'ratio_to_rutherford'"),
                RefusedAtValidation("E_REACTION_NOT_ALLOWED", "inelastic",
                                    "the ratio to Mott is not defined: give dsigma_domega")};
  return q;
}

Quantity Legendre() {
  Quantity q;
  q.name = "legendre";
  q.title = "Legendre coefficients";
  q.doc = "c.m. Legendre coefficients of the angular distribution, in the expansion named by 'expansion' "
          "(data_kinds.md section 5.2; 2l2 in b/sr as LEG/2L2 is in EXFOR).";
  q.recordKinds = kDataKinds;
  q.engineUnit = "1";
  q.units = {OneUnit("1", Condition::Is("expansion", {"a_over_a0", "rsl"})),
             AnyUnitOf("diff_cross_section", Condition::Is("expansion", {"absolute", "2l2"})),
             AnyUnitOf("cross_section", Condition::Is("expansion", {"4pi"}))};
  q.energyRoles = {"E_lab", "E_cm"};
  q.orderRole = true;
  // The apparatus as for dsigma_domega: the engine folds coefficient data with a target, beam or detector
  // (tests/angular_coefficients, conv.azr), and the converter writes what it folded.
  q.fields = {Requires("expansion"), Requires("orders"), Allows("target"), Allows("beam"), Allows("detector")};
  q.targetKinds = {"layer", "none"};
  q.norm = {NormModesWithPrior(Condition::Is("expansion", {"a_over_a0", "rsl"}), {"fixed"},
                               "normalised coefficients carry no normalization: norm is fixed (data_formats.rst, "
                               "Legendre data)")};
  q.tree = "linear maps of the c.m. coefficient leaves A_L";
  q.exfor = {",DA,,LEG", "LEG/RS", "LEG/RSL", "LEG/2L2", "LEG/4PI"};
  q.refusals = {RefusedByImporter("E_ROLE_NOT_ALLOWED", "lab-frame coefficients",
                                  "Legendre coefficients are c.m. quantities"),
                RefusedByImporter("E_NOT_IN_LIST", "RSD, RS0, 1K2, LEG/RTH", "not supported (K3)")};
  return q;
}

Quantity AnalyzingPower() {
  Quantity q;
  q.name = "analyzing_power";
  q.aliases = {"ay"};
  q.title = "Analyzing power";
  q.doc = "A_y of a spin-1/2 projectile; folded as A_y sigma over sigma.";
  q.recordKinds = kDataKinds;
  q.engineUnit = "1";
  q.units = {OneUnit("1")};
  q.energyRoles = {"E_lab", "E_cm"};
  q.angleRoles = {"theta_lab", "theta_cm"};
  q.angle = Quantity::RequiredAngle;
  q.reaction = {R::kSpinHalfProjectile, "spin-1/2 projectile"};
  q.fields = {Allows("target"), Allows("beam"), Allows("detector")};
  q.targetKinds = {"layer", "none"};
  q.tree = "A_y(r, theta), folded as A_y sigma over sigma";
  q.exfor = {",POL/DA,,ANA"};
  q.refusals = {RefusedAtValidation("E_REACTION_NOT_ALLOWED", "spin",
                                    "A_y is the vector analysing power of a spin-1/2 projectile (n, p, t, h); "
                                    "tensor analysing powers are not supported")};
  return q;
}

Quantity PhaseShift() {
  Quantity q;
  q.name = "phase_shift";
  q.title = "Phase shift";
  q.doc = "The elastic phase shift delta_{J l s} (data_formats.rst, phase shifts).";
  q.recordKinds = kDataKinds;
  q.engineUnit = "deg";
  q.units = {AnyUnitOf("angle")};
  q.energyRoles = {"E_lab", "E_cm"};
  q.reaction = {R::kElastic, "elastic"};
  q.fields = {Requires("J"), Requires("l"), Requires("s")};
  q.tree = "delta_{J l s}(r)";
  q.refusals = {RefusedAtValidation("E_REACTION_NOT_ALLOWED", "inelastic", "a phase shift is an elastic quantity"),
                RefusedAtBuild("E_FIELD_NOT_APPLICABLE", "(J, l, s) not allowed by the pair",
                               "checked against <levels>")};
  return q;
}

Quantity Yield() {
  Quantity q;
  q.name = "yield";
  q.title = "Yield";
  q.doc = "Reactions (or detected counts) per incident particle or per unit of beam charge, on a target layer, "
          "1e-24 int_{E_out}^{E_0} <sigma>/epsilon dE with int_{E_out}^{E_0} dE/epsilon = N, or on an infinite target "
          "(E_out = 1 keV); epsilon the lab stopping per active nucleus (data_kinds.md eqs. 3, 6; energy_grid.md section 5).";
  q.recordKinds = kDataKinds;
  q.engineUnit = "1/particle";
  q.units = {AnyUnitOf("yield", kNoAngle), AnyUnitOf("yield_per_sr", kAngle), AnyUnitOf("per_charge"),
             AnyUnitOf("counts_per_charge"), OneUnit("arbitrary")};
  q.energyRoles = {"E_lab"};
  q.angleRoles = {"theta_lab"};
  q.angle = Quantity::OptionalAngle;
  q.fields = {Requires("target"), Allows("beam"), Allows("detector"),
              Requires("charge_state", Condition::Dimension({"per_charge", "counts_per_charge"})),
              Requires("efficiency", Condition::Dimension({"counts_per_charge"}))};
  q.targetKinds = {"layer", "infinite"};
  q.tree = "Fold_target o sigma (or dsdo): 1e-24 integral of sigma/epsilon from E_out to E_0 (a layer: E_out where "
           "the integral of 1/epsilon reaches N; infinite: 1 keV); x 1e-6 C/(q e) per uC; x efficiency for counts";
  q.exfor = {",MLT,,TT", ",PY,,TT", "TT/CH", "PAR,MLT,G,TT/CH"};
  q.refusals = {RefusedAtValidation("E_TARGET_KIND", "target kind none or sample",
                                    "a yield needs a target with a thickness (layer) or a thick (infinite) target"),
                RefusedAtValidation("E_FIELD_REQUIRED", "per-charge unit without charge_state",
                                    "the number of beam particles per uC depends on the charge state q"),
                RefusedAtSetUp("E_LAYER_STOPPING", "a stopping not positive and finite in the layer",
                               "give a stopping defined at every energy the beam has in the layer"),
                RefusedAtSetUp("E_LAYER_EXIT_BELOW_MODEL", "the beam stops in the layer",
                               "use a target of kind infinite, or the areal density the beam crosses"),
                RefusedByImporter("E_NOT_IN_LIST", "TTY activation yields", "not supported (K9)")};
  return q;
}

Quantity Transmission() {
  Quantity q;
  q.name = "transmission";
  q.title = "Transmission";
  q.doc = "Neutron transmission T = exp(-n sigma) through a sample of thickness n (atoms/b).";
  q.recordKinds = kDataKinds;
  q.engineUnit = "1";
  q.units = {OneUnit("1")};
  q.energyRoles = {"E_lab"};
  q.reaction = {R::kSum | R::kNeutralEntrance, "neutral entrance; sums"};
  q.fields = {Requires("target"), Allows("beam"), Allows("detector")};
  q.targetKinds = {"sample"};
  q.tree = "Transform_T(n; Sum) (data_model.md section 3)";
  q.exfor = {",TRN"};
  q.refusals = {RefusedAtValidation("E_REACTION_NOT_ALLOWED", "charged_entrance",
                                    "transmission is a neutron quantity")};
  return q;
}

Quantity NeutronYield() {
  Quantity q;
  q.name = "neutron_yield";
  q.title = "Neutron-sample yield";
  q.doc = "The self-shielded yield of a neutron sample (SAMMY's yield).";
  q.recordKinds = kDataKinds;
  q.engineUnit = "1";
  q.units = {OneUnit("1")};
  q.energyRoles = {"E_lab"};
  q.reaction = {R::kSum | R::kNeutralEntrance, "neutral entrance; sums"};
  q.fields = {Requires("target"), Allows("beam"), Allows("detector")};
  q.targetKinds = {"sample"};
  q.tree = "Transform_Y(n; ...)";
  q.refusals = {RefusedAtValidation("E_REACTION_NOT_ALLOWED", "charged_entrance",
                                    "the sample yield is a neutron quantity")};
  return q;
}

Quantity BetaDelayedSpectrum() {
  Quantity q;
  q.name = "beta_delayed_spectrum";
  q.title = "beta-delayed particle spectrum";
  q.doc = "The spectrum of particles emitted after beta decay, in arbitrary units.";
  q.recordKinds = kDataKinds;
  q.engineUnit = "arbitrary";
  q.units = {OneUnit("arbitrary")};
  q.energyRoles = {"E_decay"};
  q.reaction = {0, "beta entrance pair (step 3)"};
  q.fields = {Allows("detector")};
  // The feeding amplitudes of the levels carry the spectrum's scale, so its norm may be fixed (tests/
  // beta_delayed_alpha fixes it); this rule decides before the general one for arbitrary units.
  q.norm = {NormModesWithPrior(Condition(), {"fixed", "fit"},
                               "the spectrum's scale is carried by the beta-feeding amplitudes")};
  q.tree = "the leaves with EnergyIn = DecayLab";
  q.refusals = {RefusedAtBuild("E_REACTION_NOT_ALLOWED", "non-beta entrance", "needs a beta entrance pair")};
  return q;
}

// ---------------------------------------------------------------- model data (no beam energy)

Quantity ResonanceStrength() {
  Quantity q;
  q.name = "resonance_strength";
  q.title = "Resonance strength";
  q.doc = "omega gamma of a level from observed widths at E_r (data_kinds.md section 4.7, eq. 8).";
  q.recordKinds = kDatum;
  q.engineUnit = "eV";
  q.units = {OneUnit("eV"), OneUnit("meV"), OneUnit("ueV"), OneUnit("keV")};
  q.reaction = {R::kNoTotal, "one reaction, no [total] exit"};
  q.fields = {Requires("level"), Requires("J"), Requires("reaction")};
  q.tree = "omega Gamma_a Gamma_b / Gamma at E_r";
  q.exfor = {",WID/STR", "PAR,WID/STR"};
  q.refusals = {RefusedAtBuild("E_OUT_OF_RANGE", "E_r < 0", "a subthreshold level has no strength: use anc")};
  return q;
}

Quantity ResonanceEnergy() {
  Quantity q;
  q.name = "resonance_energy";
  q.title = "Resonance energy";
  q.doc = "The observed level energy in the stated frame.";
  q.recordKinds = kDatum;
  q.engineUnit = "MeV";
  q.units = {OneUnit("MeV"), OneUnit("keV"), OneUnit("eV")};
  q.fields = {Requires("level"), Requires("J"), Requires("frame"),
              Requires("pair", Condition::Is("frame", {"E_cm", "E_lab"}))};
  q.tree = "E_lambda (observed) in the frame";
  q.exfor = {",EN"};
  return q;
}

Quantity TotalWidth() {
  Quantity q;
  q.name = "total_width";
  q.title = "Total width";
  q.doc = "The observed total width at E_r.";
  q.recordKinds = kDatum;
  q.engineUnit = "eV";
  q.units = {OneUnit("eV"), OneUnit("keV"), OneUnit("MeV")};
  q.fields = {Requires("level"), Requires("J")};
  q.tree = "Gamma (observed) at E_r";
  q.exfor = {",WID"};
  return q;
}

Quantity PartialWidth() {
  Quantity q;
  q.name = "partial_width";
  q.title = "Partial width";
  q.doc = "The observed partial width of a pair at E_r.";
  q.recordKinds = kDatum;
  q.engineUnit = "eV";
  q.units = {OneUnit("eV"), OneUnit("keV"), OneUnit("MeV")};
  q.fields = {Requires("level"), Requires("J"), Requires("pair")};
  q.tree = "Gamma_c (observed) at E_r";
  q.exfor = {"PAR,WID"};
  q.refusals = {RefusedByImporter("E_NOT_IN_LIST", ",WID,,G and ,WID,,2G", "unless the importer removes g")};
  return q;
}

Quantity BranchingRatio() {
  Quantity q;
  q.name = "branching_ratio";
  q.title = "Branching ratio";
  q.doc = "Gamma_b / Gamma of a level.";
  q.recordKinds = kDatum;
  q.engineUnit = "1";
  q.units = {OneUnit("1"), OneUnit("%")};
  q.reaction = {R::kNoTotal, "one reaction, no [total] exit"};
  q.fields = {Requires("level"), Requires("J"), Requires("reaction")};
  q.tree = "Gamma_b / Gamma";
  return q;
}

Quantity Anc() {
  Quantity q;
  q.name = "anc";
  q.title = "Asymptotic normalization coefficient";
  q.doc = "C of a bound level from its reduced width through the Whittaker function at a_c.";
  q.recordKinds = kDatum;
  q.engineUnit = "fm^-1/2";
  q.units = {OneUnit("fm^-1/2")};
  q.fields = {Requires("level"), Requires("J"), Requires("pair"), Requires("l"), Requires("s")};
  q.tree = "C from the reduced width through W(a_c)";
  q.refusals = {RefusedAtBuild("E_FIELD_NOT_APPLICABLE", "open channel", "an ANC belongs to a closed channel")};
  return q;
}

Quantity ScatteringLength() {
  Quantity q;
  q.name = "scattering_length";
  q.title = "Scattering length";
  q.doc = "Re and Im of the s-wave scattering length of a particle pair, Coulomb-modified for a charged pair "
          "(include/ScatteringLength.h).";
  q.recordKinds = kDatum;
  q.engineUnit = "fm";
  q.units = {OneUnit("fm")};
  q.fields = {Requires("pair"), Requires("s"), Allows("imag"), Requires("imag_error", Condition::Present("imag"))};
  q.tree = "a_0 (Re, Im)";
  q.exfor = {",AMP", "COH,AMP"};
  return q;
}

}  // namespace

const std::vector<NormRule> &GeneralNormRules() {
  static const std::vector<NormRule> rules = {
      NormModesWithoutPrior(Condition::Dimension({"arbitrary"}), {"fit", "profiled"},
                            "the value is in arbitrary units: its normalization is free (fit without a prior) or "
                            "profiled"),
  };
  return rules;
}

const std::vector<Quantity> &Quantities() {
  static const std::vector<Quantity> rows = {
      // data sets and extrapolations
      Sigma(), SigmaNuclear(), SigmaTotal(), DsigmaDomega(), SFactor(), ModifiedSFactor(), RatioToReference(), RatioToRutherford(), RatioToMott(), Legendre(),
      AnalyzingPower(), PhaseShift(), Yield(), Transmission(), NeutronYield(), BetaDelayedSpectrum(),
      // model data
      ResonanceStrength(), ResonanceEnergy(), TotalWidth(), PartialWidth(), BranchingRatio(), Anc(),
      ScatteringLength(),
  };
  return rows;
}

const Quantity *FindQuantity(const std::string &nameOrAlias) {
  for (const Quantity &q : Quantities()) {
    if (q.name == nameOrAlias) return &q;
    for (const std::string &a : q.aliases)
      if (a == nameOrAlias) return &q;
  }
  return nullptr;
}

std::vector<std::string> QuantityNames(const std::string &recordKind) {
  std::vector<std::string> out;
  for (const Quantity &q : Quantities())
    for (const std::string &k : q.recordKinds)
      if (k == recordKind) out.push_back(q.name);
  return out;
}

}  // namespace experiment
}  // namespace azure2
