#include "azure2/stopping/Tables.h"

#include <cstdlib>
#include <sstream>
#include <stdexcept>

#include "TableText.h"

namespace azure2 {
namespace stopping {

namespace {

std::vector<std::string> SplitTabs(const std::string &line) {
  std::vector<std::string> fields;
  std::string::size_type start = 0;
  while (true) {
    const auto tab = line.find('\t', start);
    fields.push_back(line.substr(start, tab == std::string::npos ? std::string::npos : tab - start));
    if (tab == std::string::npos) break;
    start = tab + 1;
  }
  return fields;
}

[[noreturn]] void Fail(const std::string &where, const std::string &what) {
  throw std::invalid_argument(where + ": " + what);
}

double Number(const std::string &text, const std::string &where) {
  char *end = nullptr;
  const double v = std::strtod(text.c_str(), &end);
  if (text.empty() || end == text.c_str() || *end != '\0') Fail(where, "'" + text + "' is not a number");
  return v;
}

int Integer(const std::string &text, const std::string &where) {
  const double v = Number(text, where);
  if (v != static_cast<int>(v)) Fail(where, "'" + text + "' is not an integer");
  return static_cast<int>(v);
}

/// The data rows of a table: comment lines and empty lines skipped, the header checked and skipped.
std::vector<std::pair<std::string, std::vector<std::string>>> DataRows(const std::string &text, const std::string &origin,
                                                                       const std::string &header) {
  std::vector<std::pair<std::string, std::vector<std::string>>> rows;
  std::istringstream in(text);
  std::string line;
  int number = 0;
  bool header_seen = false;
  const std::size_t columns = SplitTabs(header).size();
  while (std::getline(in, line)) {
    ++number;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty() || line[0] == '#') continue;
    const std::string where = origin + ":" + std::to_string(number);
    if (!header_seen) {
      if (line != header) Fail(where, "header is not '" + header + "'");
      header_seen = true;
      continue;
    }
    auto fields = SplitTabs(line);
    if (fields.size() != columns) {
      Fail(where, std::to_string(fields.size()) + " columns, " + std::to_string(columns) + " expected");
    }
    rows.emplace_back(where, fields);
  }
  if (!header_seen) Fail(origin, "no header line");
  return rows;
}

}  // namespace

std::string PhaseName(Phase phase) {
  switch (phase) {
    case Phase::kSolid: return "solid";
    case Phase::kLiquid: return "liquid";
    case Phase::kGas: return "gas";
  }
  return "solid";
}

Phase PhaseFromName(const std::string &name) {
  if (name == "solid") return Phase::kSolid;
  if (name == "liquid") return Phase::kLiquid;
  if (name == "gas") return Phase::kGas;
  throw std::invalid_argument("phase '" + name + "': expected solid, liquid or gas");
}

std::vector<ElementRecord> ParseElementTable(const std::string &text, const std::string &origin) {
  std::vector<ElementRecord> table;
  for (const auto &[where, f] : DataRows(text, origin, "Z\tsymbol\tA_r\tA_r_unc\tA_r_source\tI_eV\tphase")) {
    ElementRecord r;
    r.z = Integer(f[0], where);
    r.symbol = f[1];
    r.atomic_weight_u = Number(f[2], where);
    r.atomic_weight_unc_u = Number(f[3], where);
    r.atomic_weight_source = f[4];
    r.mean_excitation_eV = Number(f[5], where);
    r.phase = PhaseFromName(f[6]);
    if (r.z != static_cast<int>(table.size()) + 1) Fail(where, "Z must run 1, 2, 3, ... without gaps");
    if (!(r.atomic_weight_u > 0.0) || !(r.mean_excitation_eV > 0.0)) Fail(where, "A_r and I must be positive");
    table.push_back(r);
  }
  return table;
}

std::vector<CoefficientRow> ParseCoefficientTable(const std::string &text, const std::string &origin) {
  std::vector<CoefficientRow> table;
  const std::string header = "Z\tprojectile\tform\tA1\tA2\tA3\tA4\tA5\tfit_max_MeV_u\trel_unc\tsource\tnote";
  for (const auto &[where, f] : DataRows(text, origin, header)) {
    CoefficientRow r;
    r.z = Integer(f[0], where);
    if (f[1] == "p") r.projectile_z = 1;
    else if (f[1] == "He") r.projectile_z = 2;
    else Fail(where, "projectile '" + f[1] + "': expected p or He");
    if (f[2] == "AZ_p") r.form = FitForm::kAndersenZieglerProton;
    else if (f[2] == "Z_He") r.form = FitForm::kZieglerHelium;
    else Fail(where, "form '" + f[2] + "': expected AZ_p or Z_He");
    if ((r.form == FitForm::kAndersenZieglerProton) != (r.projectile_z == 1)) Fail(where, "form does not fit projectile");
    for (int k = 0; k < 5; ++k) r.a[k] = Number(f[3 + k], where);
    r.fit_max_MeV_u = Number(f[8], where);
    r.rel_uncertainty = Number(f[9], where);
    r.source = f[10];
    r.note = f[11];
    r.origin = where;
    if (r.z < 1 || r.z > 92) Fail(where, "Z outside 1..92");
    if (!(r.fit_max_MeV_u > 0.0) || !(r.rel_uncertainty > 0.0)) Fail(where, "fit_max and rel_unc must be positive");
    table.push_back(r);
  }
  return table;
}

const std::vector<ElementRecord> &ElementTable() {
  static const std::vector<ElementRecord> table = ParseElementTable(ElementsTableText(), "stopping/data/elements.tsv");
  return table;
}

const std::vector<CoefficientRow> &CoefficientTable() {
  static const std::vector<CoefficientRow> table =
      ParseCoefficientTable(CoefficientsTableText(), "stopping/data/electronic_coefficients.tsv");
  return table;
}

const ElementRecord &Element(int z) {
  const auto &table = ElementTable();
  if (z < 1 || z > static_cast<int>(table.size())) {
    throw std::out_of_range("stopping: element Z = " + std::to_string(z) + " outside 1.." + std::to_string(table.size()));
  }
  return table[z - 1];
}

const CoefficientRow &Coefficients(int z, int projectile_z) {
  for (const auto &row : CoefficientTable()) {
    if (row.z == z && row.projectile_z == projectile_z) return row;
  }
  throw std::out_of_range("stopping: no electronic stopping coefficients for Z2 = " + std::to_string(z) +
                          " and projectile charge " + std::to_string(projectile_z));
}

}  // namespace stopping
}  // namespace azure2
