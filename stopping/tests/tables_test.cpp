// The compiled tables (docs/upgrade/stopping.md sections 3 and 4): complete, every row with its source, and the
// parsers refuse malformed text with the file and line.
#include <stdexcept>
#include <string>

#include "azure2/physics/Elements.h"
#include "azure2/stopping/Tables.h"
#include "check.h"

using namespace azure2::stopping;

int main(int, char **) {
  const auto &elements = ElementTable();
  check::True("elements.tsv: 92 elements", elements.size() == 92);
  bool symbols = true, positive = true;
  for (const auto &e : elements) {
    symbols = symbols && e.symbol == azure2::physics::ElementSymbol(e.z);
    positive = positive && e.atomic_weight_u > 0.0 && e.mean_excitation_eV > 0.0 && !e.atomic_weight_source.empty();
  }
  check::True("elements.tsv: symbols are physics/'s ElementSymbol(Z)", symbols);
  check::True("elements.tsv: A_r, I > 0 and a source on every row", positive);
  check::Rel("A_r(Au) (CIAAW 196.97)", Element(79).atomic_weight_u, 196.97, 1e-12);
  check::Rel("I(Al) (ICRU 37: 166 eV)", Element(13).mean_excitation_eV, 166.0, 1e-12);
  check::True("O is a gas, Ta a solid", Element(8).phase == Phase::kGas && Element(73).phase == Phase::kSolid);

  const auto &rows = CoefficientTable();
  check::True("electronic_coefficients.tsv: 184 rows (92 p, 92 He)", rows.size() == 184);
  bool complete = true;
  for (int z = 1; z <= 92; ++z) {
    for (int p = 1; p <= 2; ++p) {
      const CoefficientRow &r = Coefficients(z, p);
      complete = complete && !r.source.empty() && r.rel_uncertainty > 0.0 && r.origin.find(".tsv:") != std::string::npos;
    }
  }
  check::True("every element has a p and a He row with source, uncertainty and origin", complete);
  check::Rel("Al p row A3 (ICRU 49: 2766)", Coefficients(13, 1).a[2], 2766.0, 1e-12);
  check::Rel("Au He row A1 (Ziegler 1977: 2.892)", Coefficients(79, 2).a[0], 2.892, 1e-12);

  check::Throws<std::out_of_range>("Element(0) refused", [] { Element(0); });
  check::Throws<std::out_of_range>("Coefficients(93, 1) refused", [] { Coefficients(93, 1); });
  const std::string header = "Z\tprojectile\tform\tA1\tA2\tA3\tA4\tA5\tfit_max_MeV_u\trel_unc\tsource\tnote\n";
  check::Throws<std::invalid_argument>("a wrong header is refused", [] { ParseCoefficientTable("Z\tfoo\n", "t"); });
  check::Throws<std::invalid_argument>("a non-number is refused", [&] {
    ParseCoefficientTable(header + "1\tp\tAZ_p\t1\tx\t1\t1\t1\t1\t0.04\ts\t-\n", "t");
  });
  check::Throws<std::invalid_argument>("a He row in the proton form is refused", [&] {
    ParseCoefficientTable(header + "1\tHe\tAZ_p\t1\t1\t1\t1\t1\t1\t0.04\ts\t-\n", "t");
  });
  check::Throws<std::invalid_argument>("a missing column is refused", [&] {
    ParseCoefficientTable(header + "1\tp\tAZ_p\t1\t1\t1\t1\t1\t1\t0.04\ts\n", "t");
  });
  check::Throws<std::invalid_argument>("a gap in Z is refused", [] {
    ParseElementTable("Z\tsymbol\tA_r\tA_r_unc\tA_r_source\tI_eV\tphase\n2\tHe\t4\t0\ts\t41.8\tgas\n", "t");
  });
  check::Throws<std::invalid_argument>("an unknown phase is refused", [] {
    ParseElementTable("Z\tsymbol\tA_r\tA_r_unc\tA_r_source\tI_eV\tphase\n1\tH\t1\t0\ts\t19\tplasma\n", "t");
  });
  return check::Done();
}
