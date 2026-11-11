// Materials (docs/upgrade/stopping.md sections 2, 6, 7): formulas, SRIM's default phase, the Bragg sum with phase
// factor and compound correction against LunaSRIM (TiN, Ta2O5, TiH2, LiF, CaF2, H2O, CH4 and CO2 as gases, C = 1,
// 0.93, 1.08), the uncertainty, user tables and formulas overriding the evaluation, the provenance text.
#include <cmath>
#include <stdexcept>
#include <string>

#include "azure2/stopping/Electronic.h"
#include "azure2/stopping/Nuclear.h"
#include "azure2/stopping/Stopping.h"
#include "check.h"

using namespace azure2::stopping;

int main(int argc, char **argv) {
  const std::string ref = argc > 1 ? argv[1] : "srim_reference.csv";
  const auto ta2o5 = ParseFormula("Ta2O5");
  check::True("Ta2O5: two constituents, Ta first", ta2o5.size() == 2 && ta2o5[0].z == 73 && ta2o5[1].z == 8);
  check::Rel("Ta2O5: Ta fraction 2/7", ta2o5[0].atom_fraction, 2.0 / 7.0, 1e-15);
  check::Rel("Ti0.5Al0.5N: N fraction 1/2", ParseFormula("Ti0.5Al0.5N")[2].atom_fraction, 0.5, 1e-15);
  check::Rel("CH3OH sums repeated H: H fraction 4/6", ParseFormula("CH3OH")[1].atom_fraction, 4.0 / 6.0, 1e-15);
  check::Throws<std::invalid_argument>("unknown symbol refused", [] { ParseFormula("Xx2"); });
  check::Throws<std::invalid_argument>("lower-case start refused", [] { ParseFormula("ta2o5"); });
  check::Throws<std::invalid_argument>("element above Z = 92 refused", [] { ParseFormula("Pu"); });
  check::Throws<std::invalid_argument>("zero count refused", [] { ParseFormula("Ta0O5"); });

  // SRIM's default phase: gas for one natural-gas element, solid otherwise (a gaseous compound must be said).
  check::True("default phase: O gas, Ne gas, Ta solid, H2O solid, CH4 solid",
              MaterialFromFormula("O").phase == Phase::kGas && MaterialFromFormula("Ne").phase == Phase::kGas &&
                  MaterialFromFormula("Ta").phase == Phase::kSolid && MaterialFromFormula("H2O").phase == Phase::kSolid &&
                  MaterialFromFormula("CH4").phase == Phase::kSolid);
  check::True("a new material has no compound correction (C = 1) and no density",
              MaterialFromFormula("TiN").compound_correction == 1.0 && MaterialFromFormula("TiN").mass_density_g_cm3 == 0.0);

  const Material m = MaterialFromFormula("Ta2O5");
  check::Rel("mean atomic weight of Ta2O5 with SRIM's weights (180.948, 15.999)", MeanAtomicWeight(m),
             (2 * ElementData(73).atomic_weight_u + 5 * ElementData(8).atomic_weight_u) / 7, 1e-15);
  check::Rel("atom fraction of O", AtomFraction(m, 8), 5.0 / 7.0, 1e-15);

  // Compounds against LunaSRIM: Bragg's rule with the phase factor and the compound correction, applied once.
  check::Largest largest;
  int n = 0;
  for (const auto &r : check::Rows(ref, "compound")) {
    Material c = MaterialFromFormula(r[0]);
    c.phase = r[1] == "1" ? Phase::kGas : Phase::kSolid;
    c.compound_correction = std::stod(r[2]);
    const Projectile p = Ion(std::stoi(r[3]), std::stod(r[4]), "ion");
    const StoppingValue v = Stopping(c, p, std::stod(r[5]) * 1e-3);
    const std::string tag = r[0] + " (phase " + r[1] + ", C = " + r[2] + ") Z1=" + r[3] + " E=" + r[5] + " keV";
    check::RelQuiet("S_e " + tag, v.electronic, std::stod(r[6]), 1e-10);
    check::RelQuiet("S_n " + tag, v.nuclear, std::stod(r[7]), 1e-10);
    largest.Track("S_e " + r[0] + (r[2] == "1.0" ? "" : " C=" + r[2]), v.electronic, std::stod(r[6]));
    largest.Track("S_n " + r[0], v.nuclear, std::stod(r[7]));
    ++n;
  }
  check::True("compounds: " + std::to_string(n) + " values (S_e and S_n) against LunaSRIM to 1e-10", n > 900);
  largest.Summary("compounds against LunaSRIM");

  for (const Projectile &p : {Proton(), Alpha()}) {
    for (double e : {0.02, 0.5, 2.0, 9.0}) {
      const StoppingValue v = Stopping(m, p, e);
      const double ta = ElectronicStopping(p, 73, e, Phase::kSolid, 1.0).value;
      const double o = ElectronicStopping(p, 8, e, Phase::kSolid, 1.0).value;
      const std::string tag = p.label + " at " + std::to_string(e) + " MeV";
      check::Rel("Bragg: Ta2O5 electronic = (2 Ta + 5 O)/7 (O as a solid: phase factor 0.92), " + tag, v.electronic, (2 * ta + 5 * o) / 7, 1e-14);
      const double unc = std::hypot(ElectronicRelUncertainty(p.z) * v.electronic, kNuclearRelUncertainty * v.nuclear) / v.total;
      check::Rel("uncertainty: electronic parts correlated, nuclear in quadrature, " + tag, v.rel_uncertainty, unc, 1e-14);
    }
  }
  Material bad = m;
  bad.compound_correction = 0.0;
  check::Throws<std::invalid_argument>("a compound correction of 0 refused", [&] { Stopping(bad, Proton(), 1.0); });

  // A user table replaces the evaluation for its projectile only; nodes exact, log-log between them.
  Material t = MaterialFromFormula("TiN");
  t.mass_density_g_cm3 = 5.21;
  t.user_stopping.push_back(ReadStoppingTable("# measured\nE[keV] eps[eV/1e15atoms/cm2]\n100 40\n400 20\n1600 10\n", t,
                                              Proton(), 0.03, "test table, TiN, protons"));
  check::Rel("user table at a node", Stopping(t, Proton(), 0.4).total, 20.0, 1e-14);
  check::Rel("user table log-log between nodes (E^-1/2: 200 keV -> 20 sqrt 2)", Stopping(t, Proton(), 0.2).total,
             20.0 * std::sqrt(2.0), 1e-12);
  check::True("user table flagged and its uncertainty used",
              Stopping(t, Proton(), 0.2).from_user && Stopping(t, Proton(), 0.2).rel_uncertainty == 0.03);
  check::True("user table does not apply to alphas", !Stopping(t, Alpha(), 0.2).from_user);
  Material h3 = MaterialFromFormula("TiN");
  h3.user_stopping.push_back(ReadStoppingTable("E[MeV] eps[eV/1e15atoms/cm2]\n0.1 80\n2 40\n", h3, Helium3(), 0.03, "3He table"));
  check::True("a 3He table does not apply to tritons (same mass, other charge)", !Stopping(h3, Triton(), 0.5).from_user &&
                                                                                   Stopping(h3, Helium3(), 0.5).from_user);
  check::Throws<std::out_of_range>("user table not extrapolated", [&] { Stopping(t, Proton(), 2.0); });
  check::True("provenance names the user table", StoppingProvenance(t, Proton()).find("test table, TiN") != std::string::npos);
  Material co2 = MaterialFromFormula("CO2");
  co2.phase = Phase::kGas;
  co2.compound_correction = 0.93;
  const std::string prov = StoppingProvenance(co2, Alpha());
  check::True("provenance: SRIM-2013 with its citation, LunaSRIM, the data source",
              prov.find("SRIM-2013") != std::string::npos && prov.find("268 (2010) 1818") != std::string::npos &&
                  prov.find("LunaSRIM") != std::string::npos && prov.find("SCOEF03.DAT") != std::string::npos);
  check::True("provenance: phase and compound correction used",
              prov.find("phase gas") != std::string::npos && prov.find("C = 0.93") != std::string::npos);
  std::printf("%s", prov.c_str());

  // A user formula.
  Material f = MaterialFromFormula("C");
  UserStopping formula;
  formula.projectile_z = 2;
  formula.projectile_mass_u = Alpha().mass_u;
  formula.formula = [](double e) { return 50.0 / std::sqrt(e); };
  formula.rel_uncertainty = 0.05;
  formula.provenance = "test formula";
  f.user_stopping.push_back(formula);
  check::Rel("user formula", Stopping(f, Alpha(), 4.0).total, 25.0, 1e-15);

  // Unit conversions of a table: MeV cm^2/g with the mean atomic weight, keV/um with the density.
  const double a_si = ElementData(14).atomic_weight_u;
  Material si = MaterialFromFormula("Si");
  si.mass_density_g_cm3 = 2.33;
  const auto mass = ReadStoppingTable("E[MeV] S[MeV*cm2/g]\n1 100\n2 80\n", si, Proton(), 0.02, "mass stopping");
  check::Rel("MeV cm^2/g -> eV/(1e15 atoms/cm^2): x A_r / N_A / 1e-21", mass.stopping[0], 100.0 * a_si / 6.02214076e23 / 1e-21, 1e-14);
  const auto lin = ReadStoppingTable("E[MeV] S[keV/um]\n1 10\n2 8\n", si, Proton(), 0.02, "linear stopping");
  const double n_si = 2.33 * 6.02214076e23 / a_si;  // atoms/cm^3
  check::Rel("keV/um -> eV/(1e15 atoms/cm^2): 10 keV/um = 1e8 eV/cm / n", lin.stopping[0], 1e8 / n_si * 1e15, 1e-14);
  const Material nodensity = MaterialFromFormula("Si");
  check::Throws<std::invalid_argument>("keV/um without a density refused",
                                       [&] { ReadStoppingTable("E[MeV] S[keV/um]\n1 10\n2 8\n", nodensity, Proton(), 0.02, "x"); });
  check::Throws<std::invalid_argument>("unknown unit refused",
                                       [&] { ReadStoppingTable("E[MeV] S[eV/A]\n1 10\n2 8\n", si, Proton(), 0.02, "x"); });
  check::Throws<std::invalid_argument>("descending energies refused",
                                       [&] { ReadStoppingTable("E[MeV] S[MeV*cm2/g]\n2 10\n1 8\n", si, Proton(), 0.02, "x"); });
  return check::Done();
}
