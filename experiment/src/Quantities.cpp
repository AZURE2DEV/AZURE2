// The quantity table of format 3 (Quantities.h; docs/upgrade/project_format3.md section 7.2).  One row per quantity;
// the validator reads nothing else about a quantity.  Rows follow the note's tables; where the note leaves a choice
// open the row says so in a comment.
#include "azure2/experiment/Quantities.h"

namespace azure2 {
namespace experiment {

namespace {

using R = ReactionRule;

UnitRule Dim(const std::string &dimension, Condition when = {}) { return {"", dimension, when}; }
UnitRule Tok(const std::string &unit, Condition when = {}) { return {unit, "", when}; }
QuantityField Allows(const std::string &field, Condition when = {}) { return {field, false, when}; }
QuantityField Requires(const std::string &field, Condition when = {}) { return {field, true, when}; }

const Condition kNoThm = Condition::Absent("thm");
const Condition kGammaExit = Condition::Of(Condition::GammaExit);
const Condition kRatio = Condition::Of(Condition::RatioReaction);
const Condition kAngle = Condition::Of(Condition::AngleAxis);
const Condition kNoAngle = Condition::Of(Condition::AngleAxis, true);

const std::vector<std::string> kDataKinds = {"dataset", "extrapolation"};
const std::vector<std::string> kDatum = {"datum"};

}  // namespace

const std::vector<NormRule> &GeneralNormRules() {
  static const std::vector<NormRule> rules = {
      {Condition::Dimension({"arbitrary"}), {"fit", "profiled"}, false,
       "the value is in arbitrary units: its normalization is free (fit without a prior) or profiled"},
  };
  return rules;
}

const std::vector<Quantity> &Quantities() {
  static const std::vector<Quantity> rows = {
      // ---------------------------------------------------------------- data sets and extrapolations
      {"sigma", {"cross_section"}, "Cross section",
       "Angle-integrated cross section sigma(E) of the reaction; on a thin layer target the target-averaged sigma "
       "(data_kinds.md eq. 2).",
       kDataKinds, "b",
       {Dim("cross_section"), Tok("1", kRatio), Tok("arbitrary")},
       {"E_lab", "E_cm", "E_decay"}, {}, Quantity::NoAngle, false, true,
       {R::kSum | R::kRatio | R::kNoChargedElastic, "sums and ratios; no charged-particle elastic"},
       {Allows("target", kNoThm), Allows("beam"), Allows("detector"), Allows("multipole", kGammaExit),
        Allows("secondary_gamma"), Allows("thm")},
       {"layer", "none"}, {},
       "Sigma(r) [amplitude filter]; Sum/Scale/Ratio from the reaction; Weighted if the target has isotopes; folds "
       "inside the leaf (M4)",
       {",SIG", "PAR,SIG", ",SIG,,TTA"},
       {{"E_REACTION_NOT_ALLOWED", "charged_elastic",
         "the angle-integrated elastic cross section of charged particles is infinite: give dsigma_domega or "
         "ratio_to_rutherford (LEXFOR p. S.2)", "validate"},
        {"E_FIELD_NOT_APPLICABLE", "multipole with a particle exit", "a multipole filter selects photon multipoles",
         "validate"},
        {"E_TARGET_KIND", "target kind infinite", "a thick (infinite) target gives a yield: use quantity yield",
         "validate"}}},
      {"sigma_total", {}, "Total cross section",
       "Total cross section of a neutral projectile, from the optical theorem.", kDataKinds, "b",
       {Dim("cross_section")}, {"E_lab", "E_cm"}, {}, Quantity::NoAngle, false, false,
       {R::kElastic | R::kNeutralEntrance, "elastic, neutral entrance"},
       // Apparatus allowed as for sigma (the note's row lists only the reaction rule).
       {Allows("target"), Allows("beam"), Allows("detector")},
       {"layer", "none"}, {}, "SigmaTotal(a -> a) by the optical theorem (GM 103-119)", {",SIG (SF3 TOT)"},
       {{"E_REACTION_NOT_ALLOWED", "charged_entrance", "the total cross section of charged particles is infinite",
         "validate"}}},
      {"dsigma_domega", {"dsdo"}, "Differential cross section",
       "dsigma/dOmega at the angle of the angle column; the value's frame is the angle's (theta_cm: c.m. value).",
       kDataKinds, "b/sr",
       // '1' with a ratio: section 3.5's example 12C(p,p)12C / 12C(p,p)12C @ 90 deg cm (correction to section 7.2).
       {Dim("diff_cross_section"), Tok("1", kRatio)},
       {"E_lab", "E_cm"}, {"theta_lab", "theta_cm"}, Quantity::RequiredAngle, false, false,
       {R::kSum | R::kRatio | R::kNoTotal, "sums and ratios; no [total] exit"},
       {Allows("target", kNoThm), Allows("beam"), Allows("detector"), Allows("secondary_gamma"), Allows("thm")},
       {"layer", "none"}, {}, "dsdo_cm(r, theta); lab data converted at load", {",DA", "PAR,DA"},
       {{"E_REACTION_NOT_ALLOWED", "total", "[total] sums angle-integrated photon cross sections: use quantity sigma",
         "validate"},
        {"E_MODE_NOT_ALLOWED", ",DA,,REL without norm fit", "relative data need a fitted normalization", "importer"}}},
      {"sfactor", {}, "Astrophysical S factor",
       "S(E) = E sigma(E) exp(2 pi eta) with the engine's eta and atomic masses (data_kinds.md section 4.1).",
       kDataKinds, "MeV*b", {Dim("sfactor")}, {"E_cm", "E_lab"}, {}, Quantity::NoAngle, false, false,
       {R::kChargedEntrance, "charged entrance"},
       {Allows("multipole", kGammaExit)}, {}, {}, "Scale_SF(E) o Sigma(r)",
       {",SIG,,SFC", "PAR,SIG,,SFC", "EP/PAR,SIG,,SFC"},
       {{"E_REACTION_NOT_ALLOWED", "neutral", "the S factor is defined for charged particles (Gamow factor)",
         "validate"}}},
      {"ratio_to_rutherford", {"rth"}, "Ratio to Rutherford",
       "dsigma/dOmega(c.m.) divided by the Rutherford cross section at the same c.m. angle.", kDataKinds, "1",
       {Tok("1")}, {"E_lab", "E_cm"}, {"theta_lab", "theta_cm"}, Quantity::RequiredAngle, false, false,
       {R::kElastic | R::kChargedEntrance | R::kNotIdentical, "elastic, charged, not identical"},
       {Allows("target"), Allows("beam"), Allows("detector")}, {"layer", "none"}, {},
       "Ratio(Fold[dsdo_cm(r)], Fold[rutherford(r)])", {",DA,,RTH", "PAR,DA,,RTH"},
       {{"E_REACTION_NOT_ALLOWED", "identical",
         "has identical particles; the Coulomb reference is Mott-symmetrised: use quantity 'ratio_to_mott' (K2)",
         "validate"},
        {"E_REACTION_NOT_ALLOWED", "neutral", "has no Coulomb scattering: give dsigma_domega", "validate"},
        {"E_REACTION_NOT_ALLOWED", "inelastic", "the ratio to Rutherford is not defined: give dsigma_domega",
         "validate"}}},
      {"ratio_to_mott", {}, "Ratio to Mott",
       "dsigma/dOmega(c.m.) of identical particles divided by the Mott cross section.", kDataKinds, "1", {Tok("1")},
       {"E_lab", "E_cm"}, {"theta_lab", "theta_cm"}, Quantity::RequiredAngle, false, false,
       {R::kElastic | R::kChargedEntrance | R::kIdentical, "elastic, identical particles"},
       {Allows("target"), Allows("beam"), Allows("detector")}, {"layer", "none"}, {},
       "Ratio(Fold[dsdo_cm(r)], Fold[mott(r)])", {",DA,,MOT"},
       {{"E_REACTION_NOT_ALLOWED", "non_identical", "has no identical particles: use quantity 'ratio_to_rutherford'",
         "validate"},
        {"E_REACTION_NOT_ALLOWED", "inelastic", "the ratio to Mott is not defined: give dsigma_domega", "validate"}}},
      {"legendre", {}, "Legendre coefficients",
       "c.m. Legendre coefficients of the angular distribution, in the expansion named by 'expansion' "
       "(data_kinds.md section 5.2; 2l2 in b/sr as LEG/2L2 is in EXFOR).",
       kDataKinds, "1",
       {Tok("1", Condition::Is("expansion", {"a_over_a0", "rsl"})),
        Dim("diff_cross_section", Condition::Is("expansion", {"absolute", "2l2"})),
        Dim("cross_section", Condition::Is("expansion", {"4pi"}))},
       {"E_lab", "E_cm"}, {}, Quantity::NoAngle, true, false,
       {0, ""},
       {Requires("expansion"), Requires("orders")}, {},
       {{Condition::Is("expansion", {"a_over_a0", "rsl"}), {"fixed"}, true,
         "normalised coefficients carry no normalization: norm is fixed (data_formats.rst, Legendre data)"}},
       "linear maps of the c.m. coefficient leaves A_L", {",DA,,LEG", "LEG/RS", "LEG/RSL", "LEG/2L2", "LEG/4PI"},
       {{"E_ROLE_NOT_ALLOWED", "lab-frame coefficients", "Legendre coefficients are c.m. quantities", "importer"},
        {"E_NOT_IN_LIST", "RSD, RS0, 1K2, LEG/RTH", "not supported (K3)", "importer"}}},
      {"analyzing_power", {"ay"}, "Analyzing power",
       "A_y of a spin-1/2 projectile; folded as A_y sigma over sigma.", kDataKinds, "1", {Tok("1")},
       {"E_lab", "E_cm"}, {"theta_lab", "theta_cm"}, Quantity::RequiredAngle, false, false,
       {R::kSpinHalfProjectile, "spin-1/2 projectile"},
       {Allows("target"), Allows("beam"), Allows("detector")}, {"layer", "none"}, {},
       "A_y(r, theta), folded as A_y sigma over sigma", {",POL/DA,,ANA"},
       {{"E_REACTION_NOT_ALLOWED", "spin", "A_y is the vector analysing power of a spin-1/2 projectile (n, p, t, h); "
         "tensor analysing powers are not supported", "validate"}}},
      {"phase_shift", {}, "Phase shift", "The elastic phase shift delta_{J l s} (data_formats.rst, phase shifts).",
       kDataKinds, "deg", {Dim("angle")}, {"E_lab", "E_cm"}, {}, Quantity::NoAngle, false, false,
       {R::kElastic, "elastic"},
       {Requires("J"), Requires("l"), Requires("s")}, {}, {}, "delta_{J l s}(r)", {},
       {{"E_REACTION_NOT_ALLOWED", "inelastic", "a phase shift is an elastic quantity", "validate"},
        {"E_FIELD_NOT_APPLICABLE", "(J, l, s) not allowed by the pair", "checked against <levels>", "build"}}},
      {"yield", {}, "Yield",
       "Reactions (or detected counts) per incident particle or per unit of beam charge, on a target layer or an "
       "infinite target (data_kinds.md eqs. 3, 6).",
       kDataKinds, "1/particle",
       {Dim("yield", kNoAngle), Dim("yield_per_sr", kAngle), Dim("per_charge"), Dim("counts_per_charge"),
        Tok("arbitrary")},
       {"E_lab"}, {"theta_lab"}, Quantity::OptionalAngle, false, false,
       {0, ""},
       {Requires("target"), Allows("beam"), Allows("detector"),
        Requires("charge_state", Condition::Dimension({"per_charge", "counts_per_charge"})),
        Requires("efficiency", Condition::Dimension({"counts_per_charge"}))},
       {"layer", "infinite"}, {},
       "Fold_target o sigma (or dsdo) x N 1e-24 for a layer, integral sigma/epsilon for an infinite target; "
       "x 1e-6 C/(q e) per uC; x efficiency for counts",
       {",MLT,,TT", ",PY,,TT", "TT/CH", "PAR,MLT,G,TT/CH"},
       {{"E_TARGET_KIND", "target kind none or sample", "a yield needs a target with a thickness (layer) or a thick "
         "(infinite) target", "validate"},
        {"E_FIELD_REQUIRED", "per-charge unit without charge_state", "the number of beam particles per uC depends on "
         "the charge state q", "validate"},
        {"E_NOT_IN_LIST", "TTY activation yields", "not supported (K9)", "importer"}}},
      {"transmission", {}, "Transmission",
       "Neutron transmission T = exp(-n sigma) through a sample of thickness n (atoms/b).", kDataKinds, "1",
       {Tok("1")}, {"E_lab"}, {}, Quantity::NoAngle, false, false,
       {R::kSum | R::kNeutralEntrance, "neutral entrance; sums"},
       {Requires("target"), Allows("beam"), Allows("detector")}, {"sample"}, {}, "Transform_T(n; Sum) (data_model.md "
       "section 3)", {",TRN"},
       {{"E_REACTION_NOT_ALLOWED", "charged_entrance", "transmission is a neutron quantity", "validate"}}},
      {"neutron_yield", {}, "Neutron-sample yield",
       "The self-shielded yield of a neutron sample (SAMMY's yield).", kDataKinds, "1", {Tok("1")}, {"E_lab"}, {},
       Quantity::NoAngle, false, false,
       {R::kSum | R::kNeutralEntrance, "neutral entrance; sums"},
       {Requires("target"), Allows("beam"), Allows("detector")}, {"sample"}, {}, "Transform_Y(n; ...)", {},
       {{"E_REACTION_NOT_ALLOWED", "charged_entrance", "the sample yield is a neutron quantity", "validate"}}},
      {"beta_delayed_spectrum", {}, "beta-delayed particle spectrum",
       "The spectrum of particles emitted after beta decay, in arbitrary units.", kDataKinds, "arbitrary",
       {Tok("arbitrary")}, {"E_decay"}, {}, Quantity::NoAngle, false, false,
       {0, "beta entrance pair (step 3)"},
       {Allows("detector")}, {}, {}, "the leaves with EnergyIn = DecayLab", {},
       {{"E_REACTION_NOT_ALLOWED", "non-beta entrance", "needs a beta entrance pair", "build"}}},
      // ---------------------------------------------------------------- model data (no beam energy)
      {"resonance_strength", {}, "Resonance strength",
       "omega gamma of a level from observed widths at E_r (data_kinds.md section 4.7, eq. 8).", kDatum, "eV",
       {Tok("eV"), Tok("meV"), Tok("ueV"), Tok("keV")}, {}, {}, Quantity::NoAngle, false, false,
       {0, "one reaction"},
       {Requires("level"), Requires("J"), Requires("reaction")}, {}, {}, "omega Gamma_a Gamma_b / Gamma at E_r",
       {",WID/STR", "PAR,WID/STR"},
       {{"E_OUT_OF_RANGE", "E_r < 0", "a subthreshold level has no strength: use anc", "build"}}},
      {"resonance_energy", {}, "Resonance energy", "The observed level energy in the stated frame.", kDatum, "MeV",
       {Tok("MeV"), Tok("keV"), Tok("eV")}, {}, {}, Quantity::NoAngle, false, false, {},
       {Requires("level"), Requires("J"), Requires("frame"), Requires("pair", Condition::Is("frame", {"E_cm", "E_lab"}))},
       {}, {}, "E_lambda (observed) in the frame", {",EN"}, {}},
      {"total_width", {}, "Total width", "The observed total width at E_r.", kDatum, "eV",
       {Tok("eV"), Tok("keV"), Tok("MeV")}, {}, {}, Quantity::NoAngle, false, false, {},
       {Requires("level"), Requires("J")}, {}, {}, "Gamma (observed) at E_r", {",WID"}, {}},
      {"partial_width", {}, "Partial width", "The observed partial width of a pair at E_r.", kDatum, "eV",
       {Tok("eV"), Tok("keV"), Tok("MeV")}, {}, {}, Quantity::NoAngle, false, false, {},
       {Requires("level"), Requires("J"), Requires("pair")}, {}, {}, "Gamma_c (observed) at E_r", {"PAR,WID"},
       {{"E_NOT_IN_LIST", ",WID,,G and ,WID,,2G", "unless the importer removes g", "importer"}}},
      {"branching_ratio", {}, "Branching ratio", "Gamma_b / Gamma of a level.", kDatum, "1", {Tok("1"), Tok("%")}, {},
       {}, Quantity::NoAngle, false, false, {},
       {Requires("level"), Requires("J"), Requires("reaction")}, {}, {}, "Gamma_b / Gamma", {}, {}},
      {"anc", {}, "Asymptotic normalization coefficient",
       "C of a bound level from its reduced width through the Whittaker function at a_c.", kDatum, "fm^-1/2",
       {Tok("fm^-1/2")}, {}, {}, Quantity::NoAngle, false, false, {},
       {Requires("level"), Requires("J"), Requires("pair"), Requires("l"), Requires("s")}, {}, {},
       "C from the reduced width through W(a_c)", {},
       {{"E_FIELD_NOT_APPLICABLE", "open channel", "an ANC belongs to a closed channel", "build"}}},
      {"scattering_length", {}, "Scattering length",
       "Re and Im of the scattering length of a neutral pair (include/ScatteringLength.h).", kDatum, "fm",
       {Tok("fm")}, {}, {}, Quantity::NoAngle, false, false,
       {R::kNeutralPair, "neutral pair"},
       {Requires("pair"), Requires("s"), Allows("imag"), Requires("imag_error", Condition::Present("imag"))}, {}, {},
       "a_0 (Re, Im)", {",AMP", "COH,AMP"},
       {{"E_REACTION_NOT_ALLOWED", "charged_pair", "the scattering length is defined for a neutral pair", "validate"}}},
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
