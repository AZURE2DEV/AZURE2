// Materials (docs/upgrade/stopping.md sections 2, 6, 7): formulas, Bragg additivity and its uncertainty, user tables
// and formulas overriding the evaluation with their provenance, unit conversions of a user table.
#include <cmath>
#include <stdexcept>
#include <string>

#include "azure2/stopping/Electronic.h"
#include "azure2/stopping/Nuclear.h"
#include "azure2/stopping/Stopping.h"
#include "check.h"

using namespace azure2::stopping;

int main(int, char **) {
  const auto ta2o5 = ParseFormula("Ta2O5");
  check::True("Ta2O5: two constituents, Ta first", ta2o5.size() == 2 && ta2o5[0].z == 73 && ta2o5[1].z == 8);
  check::Rel("Ta2O5: Ta fraction 2/7", ta2o5[0].atom_fraction, 2.0 / 7.0, 1e-15);
  const auto tialn = ParseFormula("Ti0.5Al0.5N");
  check::Rel("Ti0.5Al0.5N: N fraction 1/2", tialn[2].atom_fraction, 0.5, 1e-15);
  check::Rel("CH3OH sums repeated H: H fraction 4/6", ParseFormula("CH3OH")[1].atom_fraction, 4.0 / 6.0, 1e-15);
  check::Throws<std::invalid_argument>("unknown symbol refused", [] { ParseFormula("Xx2"); });
  check::Throws<std::invalid_argument>("lower-case start refused", [] { ParseFormula("ta2o5"); });
  check::Throws<std::invalid_argument>("element above Z = 92 refused", [] { ParseFormula("Pu"); });
  check::Throws<std::invalid_argument>("zero count refused", [] { ParseFormula("Ta0O5"); });

  const Material m = MaterialFromFormula("Ta2O5", Phase::kSolid, 8.2);
  check::Rel("mean atomic weight of Ta2O5", MeanAtomicWeight(m), (2 * 180.95 + 5 * 15.999) / 7, 1e-12);
  check::Rel("atom fraction of O", AtomFraction(m, 8), 5.0 / 7.0, 1e-15);
  for (const Projectile &p : {Proton(), Alpha()}) {
    for (double e : {0.02, 0.5, 2.0, 9.0}) {
      const StoppingValue v = Stopping(m, p, e);
      const auto ta = ElectronicStopping(p, 73, e);
      const auto o = ElectronicStopping(p, 8, e);
      const double nta = NuclearStopping(p, 73, Element(73).atomic_weight_u, e);
      const double no = NuclearStopping(p, 8, Element(8).atomic_weight_u, e);
      const std::string tag = p.label + " at " + std::to_string(e) + " MeV";
      check::Rel("Bragg: Ta2O5 electronic = (2 Ta + 5 O)/7, " + tag, v.electronic, (2 * ta.value + 5 * o.value) / 7, 1e-14);
      check::Rel("Bragg: Ta2O5 nuclear = (2 Ta + 5 O)/7, " + tag, v.nuclear, (2 * nta + 5 * no) / 7, 1e-14);
      const double unc = std::hypot((2 * ta.value * ta.rel_uncertainty + 5 * o.value * o.rel_uncertainty) / 7,
                                    kNuclearRelUncertainty * v.nuclear) / v.total;
      check::Rel("uncertainty: electronic parts correlated, nuclear in quadrature, " + tag, v.rel_uncertainty, unc, 1e-14);
    }
  }
  check::Rel("per active atom = per atom / fraction: Ta2O5 per Ta = Ta + 2.5 O (electronic, 1 MeV p)",
             Stopping(m, Proton(), 1.0).electronic / AtomFraction(m, 73),
             ElectronicStopping(Proton(), 73, 1.0).value + 2.5 * ElectronicStopping(Proton(), 8, 1.0).value, 1e-14);
  check::Rel("compound correction factor is 1 in this branch", CompoundCorrectionFactor(m, Proton(), 1.0), 1.0, 0.0);

  // A user table replaces the evaluation for its projectile only; nodes exact, log-log between them.
  Material t = MaterialFromFormula("TiN", Phase::kSolid, 5.21);
  t.user_stopping.push_back(ReadStoppingTable("# measured\nE[keV] eps[eV/1e15atoms/cm2]\n100 40\n400 20\n1600 10\n", t,
                                              Proton(), 0.03, "test table, TiN, protons"));
  check::Rel("user table at a node", Stopping(t, Proton(), 0.4).total, 20.0, 1e-14);
  check::Rel("user table log-log between nodes (E^-1/2: 200 keV -> 20 sqrt 2)", Stopping(t, Proton(), 0.2).total,
             20.0 * std::sqrt(2.0), 1e-12);
  check::True("user table flagged and its uncertainty used",
              Stopping(t, Proton(), 0.2).from_user && Stopping(t, Proton(), 0.2).rel_uncertainty == 0.03);
  check::True("user table does not apply to alphas", !Stopping(t, Alpha(), 0.2).from_user);
  Material h3 = MaterialFromFormula("TiN", Phase::kSolid, 5.21);
  h3.user_stopping.push_back(ReadStoppingTable("E[MeV] eps[eV/1e15atoms/cm2]\n0.1 80\n2 40\n", h3, Helium3(), 0.03, "3He table"));
  check::True("a 3He table does not apply to tritons (same mass, other charge)", !Stopping(h3, Triton(), 0.5).from_user &&
                                                                                   Stopping(h3, Helium3(), 0.5).from_user);
  check::Throws<std::out_of_range>("user table not extrapolated", [&] { Stopping(t, Proton(), 2.0); });
  check::True("provenance names the user table", StoppingProvenance(t, Proton()).find("test table, TiN") != std::string::npos);
  const std::string prov = StoppingProvenance(m, Alpha());
  check::True("provenance names the table row and the phase of O", prov.find("electronic_coefficients.tsv:") != std::string::npos &&
                                                                     prov.find("as gas") != std::string::npos);

  // A user formula.
  Material f = MaterialFromFormula("C", Phase::kSolid, 2.0);
  UserStopping formula;
  formula.projectile_z = 2;
  formula.projectile_mass_u = Alpha().mass_u;
  formula.formula = [](double e) { return 50.0 / std::sqrt(e); };
  formula.rel_uncertainty = 0.05;
  formula.provenance = "test formula";
  f.user_stopping.push_back(formula);
  check::Rel("user formula", Stopping(f, Alpha(), 4.0).total, 25.0, 1e-15);

  // Unit conversions of a table: MeV cm^2/g with the mean atomic weight, keV/um with the density.
  const double a_si = Element(14).atomic_weight_u;
  Material si = MaterialFromFormula("Si", Phase::kSolid, 2.33);
  const auto mass = ReadStoppingTable("E[MeV] S[MeV*cm2/g]\n1 100\n2 80\n", si, Proton(), 0.02, "mass stopping");
  check::Rel("MeV cm^2/g -> eV/(1e15 atoms/cm^2): x A_r / N_A / 1e-21", mass.stopping[0], 100.0 * a_si / 6.02214076e23 / 1e-21, 1e-14);
  const auto lin = ReadStoppingTable("E[MeV] S[keV/um]\n1 10\n2 8\n", si, Proton(), 0.02, "linear stopping");
  const double n_si = 2.33 * 6.02214076e23 / a_si;  // atoms/cm^3
  check::Rel("keV/um -> eV/(1e15 atoms/cm^2): 10 keV/um = 1e8 eV/cm / n", lin.stopping[0], 1e8 / n_si * 1e15, 1e-14);
  Material nodensity = MaterialFromFormula("Si", Phase::kSolid, 0.0);
  check::Throws<std::invalid_argument>("keV/um without a density refused",
                                       [&] { ReadStoppingTable("E[MeV] S[keV/um]\n1 10\n2 8\n", nodensity, Proton(), 0.02, "x"); });
  check::Throws<std::invalid_argument>("unknown unit refused",
                                       [&] { ReadStoppingTable("E[MeV] S[eV/A]\n1 10\n2 8\n", si, Proton(), 0.02, "x"); });
  check::Throws<std::invalid_argument>("descending energies refused",
                                       [&] { ReadStoppingTable("E[MeV] S[MeV*cm2/g]\n2 10\n1 8\n", si, Proton(), 0.02, "x"); });
  return check::Done();
}
