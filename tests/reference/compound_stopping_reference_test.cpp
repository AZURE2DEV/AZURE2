// Reference check for the compound stopping cross section (Bragg's rule).
//
// Target integration uses the stopping cross section and the areal density
// only through their product, so the two must count the same entity. The
// checks below are identities against the single-element values:
//
//   per active atom A of A_xB_y :  eps = eps_A + (y/x) eps_B
//   per average atom (no active) :  eps = (x eps_A + y eps_B)/(x+y)
//
// and the generated AZURE2 equation must carry the same weights as the
// numeric function. Needs erya/data/SRIM2013.xml, so the working directory
// is external/ (see CMakeLists.txt).

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "SRIMUtilities.h"

static int failures = 0;

static void check(const char *what, double got, double want, double rel = 1e-12) {
  double denom = std::fabs(want) > 0 ? std::fabs(want) : 1.0;
  bool ok = std::fabs(got - want) / denom <= rel;
  std::printf("%-58s %s  (%.10e vs %.10e)\n", what, ok ? "PASS" : "FAIL", got, want);
  if (!ok) failures++;
}

static void check_contains(const char *what, const std::string &s, const std::string &needle) {
  bool ok = s.find(needle) != std::string::npos;
  std::printf("%-58s %s\n", what, ok ? "PASS" : "FAIL");
  if (!ok) failures++;
}

int main() {
  SRIMUtilities srim;
  if (!srim.isDataLoaded()) {
    std::printf("SRIM2013.xml not found -- run from external/\n");
    return 2;
  }
  const int Si = 14, O = 8, Ta = 73;
  SRIMElementData dSi = srim.readSRIMDataForElement(Si);
  SRIMElementData dO = srim.readSRIMDataForElement(O);
  SRIMElementData dTa = srim.readSRIMDataForElement(Ta);
  if (!dSi.isValid || !dO.isValid || !dTa.isValid) {
    std::printf("missing element data\n");
    return 2;
  }

  std::vector<CompoundElement> sio2 = srim.parseCompoundFormula("SiO2");
  std::vector<CompoundElement> ta2o5 = srim.parseCompoundFormula("Ta2O5");
  check("SiO2 parses to two elements", sio2.size(), 2);
  check("Ta2O5 parses to two elements", ta2o5.size(), 2);

  for (double E : {300.0, 1000.0, 3000.0}) {
    double eSi = srim.calculateZieglerStoppingPower(E, dSi);
    double eO = srim.calculateZieglerStoppingPower(E, dO);
    double eTa = srim.calculateZieglerStoppingPower(E, dTa);
    char label[96];

    std::snprintf(label, sizeof label, "SiO2 per Si atom = eps_Si + 2 eps_O  (%g keV)", E);
    check(label, srim.calculateCompoundStoppingPower(E, sio2, Si), eSi + 2.0 * eO);
    std::snprintf(label, sizeof label, "SiO2 per O atom  = eps_O + eps_Si/2  (%g keV)", E);
    check(label, srim.calculateCompoundStoppingPower(E, sio2, O), eO + 0.5 * eSi);
    std::snprintf(label, sizeof label, "SiO2 average atom = (eps_Si + 2 eps_O)/3  (%g keV)", E);
    check(label, srim.calculateCompoundStoppingPower(E, sio2), (eSi + 2.0 * eO) / 3.0);
    std::snprintf(label, sizeof label, "SiO2 with a Z not in the formula falls back to average (%g keV)", E);
    check(label, srim.calculateCompoundStoppingPower(E, sio2, 6), (eSi + 2.0 * eO) / 3.0);
    std::snprintf(label, sizeof label, "Ta2O5 per Ta atom = eps_Ta + 2.5 eps_O  (%g keV)", E);
    check(label, srim.calculateCompoundStoppingPower(E, ta2o5, Ta), eTa + 2.5 * eO);
  }

  // The equation the GUI writes into the .azr must carry the same weights.
  std::vector<double> p;
  std::string eq = srim.generateCompoundAZUREEquation(sio2, p, Si);
  check("SiO2/Si equation has 16 parameters", p.size(), 16);
  check_contains("SiO2/Si equation weights Si by 1", eq, "(1*((a0*");
  check_contains("SiO2/Si equation weights O by 2", eq, "+ 2*((a8*");
  eq = srim.generateCompoundAZUREEquation(sio2, p);
  check_contains("SiO2 average equation weights Si by 1/3", eq, "(0.333333*((a0*");
  check_contains("SiO2 average equation weights O by 2/3", eq, "+ 0.666667*((a8*");
  eq = srim.generateCompoundAZUREEquation(ta2o5, p, Ta);
  check_contains("Ta2O5/Ta equation weights Ta by 1", eq, "(1*((a0*");
  check_contains("Ta2O5/Ta equation weights O by 2.5", eq, "+ 2.5*((a8*");

  std::printf("%s: %d failure(s)\n", failures ? "FAILED" : "OK", failures);
  return failures ? 1 : 0;
}
