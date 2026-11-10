// The model data of a project and their rows of the objective's term T9 (ModelData.h).
#include "ModelData.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <istream>
#include <limits>
#include <set>
#include <sstream>

#include "AChannel.h"
#include "ALevel.h"
#include "AZUREGrad.h"
#include "CNuc.h"
#include "Config.h"
#include "ExperimentAdapter.h"
#include "JGroup.h"
#include "PPair.h"
#include "azure2/experiment/ProjectReader.h"
#include "azure2/experiment/Schema.h"
#include "azure2/io/RecordReader.h"
#include "azure2/physics/Elements.h"
#include "azure2/physics/LevelParameters.h"
#include "azure2/physics/ResonanceQuantities.h"

namespace ex = azure2::experiment;
namespace phys = azure2::physics;

int ModelDataSet::NumResiduals() const {
  int n = 0;
  for (const ModelDataRow &row : rows) n += row.isScatteringLength ? ScatteringLength::NumResiduals(row.constraint) : 1;
  return n;
}

int ModelDataSet::NumSectionResiduals() const {
  int n = 0;
  for (const ModelDataRow &row : rows)
    if (row.inSection) n += row.isScatteringLength ? ScatteringLength::NumResiduals(row.constraint) : 1;
  return n;
}

namespace ModelData {

namespace {

const char *kSection = "modelData";

/// Central-difference steps of the level quantities' derivatives: an energy (MeV) by kEnergyStep max(1, |E|), an
/// amplitude (MeV^1/2) by kAmplitudeStep max(|gamma|, kAmplitudeFloor).  The observed widths are smooth in both
/// on these scales (the penetrability and shift function vary on the scale of the channel energy, a width as
/// gamma^2); what limits the accuracy is the noise of the channel functions, about 1e-8 relative in dS/dE
/// (level_parameters.md section 1), and in the formal basis the 1e-6 MeV tolerance of the Brune iteration.  The
/// accuracy is that of AZURE_GRAD_CHECK on tests/model_data (the report of the branch states it).
constexpr double kEnergyStep = 1e-5;
constexpr double kAmplitudeStep = 1e-5;
constexpr double kAmplitudeFloor = 1e-3;

/// The levels a model datum names are those of the tables made from the <levels> rows (LevelTablesOfRows): their
/// indices must be the compound's.  "" when they are, else what differs (a defect of one of the two readers).
std::string SameLevel(CNuc &compound, const ex::ModelDatum &d, const ex::LevelTables &levels) {
  for (const ex::LevelInfo &l : levels.levels) {
    if (l.group != d.group || l.level != d.level) continue;
    if (d.group < 1 || d.group > compound.NumJGroups()) break;
    JGroup *group = compound.GetJGroup(d.group);
    if (d.level < 1 || d.level > group->NumLevels()) break;
    if ((int)std::lround(2.0 * group->GetJ()) == l.twoJ && group->GetPi() == l.parity &&
        group->GetLevel(d.level)->GetE() == l.energy && group->GetLevel(d.level)->IDs() == l.ids)
      return "";
    break;
  }
  return "level " + std::to_string(d.levelId) + " of the <levels> tables is not the compound nucleus's J^pi group " +
         std::to_string(d.group) + ", level " + std::to_string(d.level) + " (internal error: report it)";
}

/// The value of one level quantity from the observed slots (CNuc::ObservedOfFitValues), in the engine unit.
double LevelValue(CNuc *compound, const Config &configure, const ex::ModelDatum &d) {
  JGroup *group = compound->GetJGroup(d.group);
  ALevel *level = group->GetLevel(d.level);
  const unsigned char width = (unsigned char)phys::ObservedKind::Width;
  auto keyOf = [&](int ch) { return compound->GetPair(group->GetChannel(ch)->GetPairNum())->GetPairKey(); };
  // Gamma (MeV) to one pair (pairKey > 0) or in all, from the observed W of the channels that are widths.
  auto gamma = [&](int pairKey) {
    std::vector<double> w;
    for (int ch = 1; ch <= group->NumChannels(); ch++)
      if (level->GetObservedKind(ch) == width && (pairKey == 0 || keyOf(ch) == pairKey)) w.push_back(level->GetBigGamma(ch));
    return phys::total_width(w);
  };
  switch (d.quantity) {
    case ex::ModelQuantity::ResonanceStrength: {
      const PPair *in = compound->GetPairFromKey(d.entrancePair);
      const double omega = phys::statistical_factor((int)std::lround(2.0 * group->GetJ()), (int)std::lround(2.0 * in->GetJ(1)),
                                                    (int)std::lround(2.0 * in->GetJ(2)), in->IsIdentical());
      return phys::resonance_strength(omega, gamma(d.entrancePair), gamma(d.exitPair), gamma(0)).value;
    }
    case ex::ModelQuantity::BranchingRatio:
      return phys::branching_ratio(gamma(d.exitPair), gamma(0)).value;
    case ex::ModelQuantity::TotalWidth:
      return gamma(0);
    case ex::ModelQuantity::PartialWidth:
      return gamma(d.pair);
    case ex::ModelQuantity::ResonanceEnergy: {
      const double E = level->GetTransformE();  // the observed (Brune) energy, E_x
      if (d.frame == "E_x") return E;
      const PPair *p = compound->GetPairFromKey(d.pair);
      const double eCM = E - (p->GetSepE() + p->GetExE());
      if (d.frame == "E_cm") return eCM;
      if (configure.relativistic) return ex::LabFromCM(p->GetM(1) * uconv, p->GetM(2) * uconv, eCM);
      return eCM / ex::BeamCMFactor(PartitionOf(p));
    }
    case ex::ModelQuantity::Anc:
      for (int ch = 1; ch <= group->NumChannels(); ch++) {
        const AChannel *c = group->GetChannel(ch);
        if (keyOf(ch) != d.pair || c->GetL() != d.l || (int)std::lround(2.0 * c->GetS()) != d.twoS) continue;
        // C of a closed channel; 0 once the level has moved above the threshold (no ANC there).
        if (level->GetObservedKind(ch) == (unsigned char)phys::ObservedKind::ANC) return std::fabs(level->GetObservedValue(ch));
      }
      return 0.0;
    case ex::ModelQuantity::ScatteringLength:
      break;
  }
  return std::numeric_limits<double>::quiet_NaN();
}

bool HasLevelRows(const ModelDataSet &set) {
  for (const ModelDataRow &row : set.rows)
    if (!row.isScatteringLength) return true;
  return false;
}

/// The model values of the level rows (index = row; scattering-length rows NaN), from one transformation of the
/// current fit values; NaN where the transformation fails.
void LevelValues(CNuc *compound, const Config &configure, const ModelDataSet &set, std::vector<double> &q) {
  q.assign(set.rows.size(), std::numeric_limits<double>::quiet_NaN());
  if (!HasLevelRows(set) || !compound->ObservedOfFitValues(configure)) return;
  for (std::size_t k = 0; k < set.rows.size(); k++)
    if (!set.rows[k].isScatteringLength) q[k] = LevelValue(compound, configure, set.rows[k].datum);
}

/// The residual rows of row k, appended to r (q: LevelValues).
void RowResiduals(CNuc *compound, const Config &configure, const ModelDataRow &row, double q, vector_r &r) {
  if (row.isScatteringLength)
    ScatteringLength::ConstraintResiduals(compound, configure, row.constraint, r);
  else
    r.push_back((q - row.datum.value) / row.datum.error);
}

/// d q_k / d theta of the level rows for one level parameter theta (E of level la of group j for ch = 0, else its
/// amplitude in channel ch), by a central difference.
struct Column {
  int j = 0, la = 0, ch = 0;
  std::vector<double> dq;  ///< index = row (scattering-length rows 0)
};

/// The columns of every level parameter the level rows depend on that has a packed column in pmap: the parameters of
/// the rows' J^pi groups, of all groups with external capture (a photon width's external part depends on the final
/// state).  Each parameter is set back exactly (CNuc::SetFitLevelParameter); Park's shift functions, which the
/// transformation moved with the energies, are recomputed at the end.
std::vector<Column> LevelJacobian(CNuc *compound, const Config &configure, const ModelDataSet &set,
                                  const ParamIndexMap &pmap) {
  std::vector<Column> cols;
  std::set<int> groups;
  for (const ModelDataRow &row : set.rows)
    if (!row.isScatteringLength) groups.insert(row.datum.group);
  const bool allGroups = configure.paramMask & Config::USE_EXTERNAL_CAPTURE;
  std::vector<double> qp, qm;
  for (int j = 1; j <= compound->NumJGroups(); j++) {
    if (!allGroups && !groups.count(j)) continue;
    JGroup *group = compound->GetJGroup(j);
    for (int la = 1; la <= group->NumLevels(); la++)
      for (int ch = 0; ch <= group->NumChannels(); ch++) {
        const int f = ch == 0 ? pmap.EIndexFlat(j, la) : pmap.GammaIndexFlat(j, la, ch);
        if (f < 0 || f >= pmap.NumFull() || pmap.FullToPacked(f) < 0) continue;
        ALevel *level = group->GetLevel(la);
        const double v = ch == 0 ? level->GetFitE() : level->GetFitGamma(ch);
        const double h = ch == 0 ? kEnergyStep * std::max(1.0, std::fabs(v)) : kAmplitudeStep * std::max(std::fabs(v), kAmplitudeFloor);
        const double up = v + h, down = v - h;
        compound->SetFitLevelParameter(j, la, ch, up);
        LevelValues(compound, configure, set, qp);
        compound->SetFitLevelParameter(j, la, ch, down);
        LevelValues(compound, configure, set, qm);
        compound->SetFitLevelParameter(j, la, ch, v);
        Column c;
        c.j = j;
        c.la = la;
        c.ch = ch;
        c.dq.assign(set.rows.size(), 0.0);
        for (std::size_t k = 0; k < set.rows.size(); k++)
          if (!set.rows[k].isScatteringLength) c.dq[k] = (qp[k] - qm[k]) / (up - down);
        cols.push_back(c);
      }
  }
  if (configure.paramMask & Config::USE_PARK_FORMALISM) compound->CalcShiftFunctions(configure);
  return cols;
}

void AddToAccum(const Column &c, double v, GradAccum &accum) {
  if (c.ch == 0)
    accum.AddE(c.j, c.la, v);
  else
    accum.AddGamma(c.j, c.la, c.ch, v);
}

ex::ModelDatum ScatteringLengthLine(const ScatteringLength::Constraint &c) {
  ex::ModelDatum d;
  d.quantity = ex::ModelQuantity::ScatteringLength;
  d.quantityName = "scattering_length";
  d.name = "<scatteringLengths> line " + std::to_string(c.line);
  d.line = c.line;
  return d;
}

void Print(const Config &configure, const ex::Diagnostic &d) {
  configure.outStream << (d.severity == ex::Severity::Warning ? "WARNING: " : "ERROR: ") << "<modelData> line " << d.line
                      << " [" << d.code << "]: " << d.message << std::endl;
}

}  // namespace

std::shared_ptr<const ModelDataSet> Read(CNuc &compound, const Config &configure,
                                         const std::vector<azure2::io::LevelRow> &rows) {
  auto set = std::make_shared<ModelDataSet>();
  if (ScatteringLength::ReadConstraints(&compound, configure, set->scatteringLengthLines) == -1) return nullptr;
  for (const ScatteringLength::Constraint &c : set->scatteringLengthLines)
    if (ScatteringLength::NumResiduals(c) > 0) {
      ModelDataRow row;
      row.constraint = c;
      row.datum = ScatteringLengthLine(c);
      set->rows.push_back(row);
    }

  const ex::LevelTables levels = ex::LevelTablesOfRows(rows);
  // <modelData>: format-3 datum records, read and validated by the format-3 table, then found in <levels>.  Line
  // numbers count from the section's first line after its tag.
  std::unique_ptr<std::istream> section = configure.SectionStream(kSection);
  if (!section) return set;
  set->hasSection = true;
  std::vector<std::string> lines;
  std::string line;
  std::getline(*section, line);  // <modelData>
  while (std::getline(*section, line) && line.find("</modelData>") == std::string::npos) lines.push_back(line);
  const azure2::io::RecordReadResult read = azure2::io::ReadRecords(lines, 1, ex::Schema::KindsOfSection(kSection));
  std::vector<ex::Diagnostic> diagnostics;
  for (const azure2::io::RecordDiagnostic &d : read.diagnostics) diagnostics.push_back(ex::FromRecordDiagnostic(d, kSection));
  const ex::ValidationResult validated = ex::Schema::Validate({{kSection, read.set.records}});
  diagnostics.insert(diagnostics.end(), validated.diagnostics.begin(), validated.diagnostics.end());
  bool ok = true;
  for (const ex::Diagnostic &d : diagnostics) ok = ok && d.severity == ex::Severity::Warning;
  ex::ModelDataBuild built;
  if (ok) {
    built = ex::BuildModelData(validated, levels.levels, levels.pairs);
    diagnostics.insert(diagnostics.end(), built.diagnostics.begin(), built.diagnostics.end());
    ok = built.Ok();
  }
  for (const ex::Diagnostic &d : diagnostics) Print(configure, d);
  if (!ok) return nullptr;
  for (const ex::ModelDatum &d : built.data) {
    if (d.quantity != ex::ModelQuantity::ScatteringLength) {
      const std::string why = SameLevel(compound, d, levels);
      if (!why.empty()) {
        configure.outStream << "ERROR: <modelData> line " << d.line << ": datum '" << d.name << "': " << why << "."
                            << std::endl;
        return nullptr;
      }
    }
    ModelDataRow row;
    row.datum = d;
    row.inSection = true;
    row.isScatteringLength = d.quantity == ex::ModelQuantity::ScatteringLength;
    if (row.isScatteringLength) {
      ScatteringLength::Constraint &c = row.constraint;
      c.pairKey = d.pair;
      c.s = d.coherent ? -1.0 : 0.5 * d.twoS;
      c.value = d.value;
      c.sigma = d.error;
      if (d.hasImag) {
        c.imValue = d.imag;
        c.imSigma = d.imagError;
      }
      c.line = d.line;
      const std::string why = ScatteringLength::ResolveConstraint(&compound, c);
      if (!why.empty()) {
        configure.outStream << "ERROR: <modelData> line " << d.line << ": datum '" << d.name << "': " << why << "."
                            << std::endl;
        return nullptr;
      }
    }
    set->rows.push_back(row);
  }
  return set;
}

void Residuals(CNuc *compound, const Config &configure, const ModelDataSet &set, vector_r &r) {
  std::vector<double> q;
  LevelValues(compound, configure, set, q);
  for (std::size_t k = 0; k < set.rows.size(); k++) RowResiduals(compound, configure, set.rows[k], q[k], r);
}

double Chi2(CNuc *compound, const Config &configure, const ModelDataSet &set, double *sectionPart, double *linesPart) {
  std::vector<double> q;
  LevelValues(compound, configure, set, q);
  double chi2 = 0.0, section = 0.0, lines = 0.0;
  vector_r r;
  for (std::size_t k = 0; k < set.rows.size(); k++) {
    r.clear();
    RowResiduals(compound, configure, set.rows[k], q[k], r);
    for (double x : r) {
      chi2 += x * x;
      (set.rows[k].inSection ? section : lines) += x * x;
    }
  }
  if (sectionPart) *sectionPart = section;
  if (linesPart) *linesPart = lines;
  return chi2;
}

bool Gradient(CNuc *compound, const Config &configure, const ModelDataSet &set, const ParamIndexMap &pmap,
              const vector_matrix_r *shiftDeriv, GradAccum &accum) {
  for (const ModelDataRow &row : set.rows)
    if (row.isScatteringLength &&
        !ScatteringLength::ConstraintGradient(compound, configure, row.constraint, accum, shiftDeriv, -1))
      return false;
  if (!HasLevelRows(set)) return true;
  // d chi2 / d theta = sum_k 2 r_k / e_k d q_k / d theta over the level rows.
  std::vector<double> q;
  LevelValues(compound, configure, set, q);
  for (const Column &c : LevelJacobian(compound, configure, set, pmap)) {
    double g = 0.0;
    for (std::size_t k = 0; k < set.rows.size(); k++) {
      const ModelDataRow &row = set.rows[k];
      if (!row.isScatteringLength) g += 2.0 * (q[k] - row.datum.value) / (row.datum.error * row.datum.error) * c.dq[k];
    }
    AddToAccum(c, g, accum);
  }
  return true;
}

bool AppendRows(CNuc *compound, const Config &configure, const ModelDataSet &set, const ParamIndexMap &pmap,
                const vector_matrix_r *shiftDeriv, int nCols, vector_r &residuals, vector_r &jac) {
  if (set.rows.empty()) return true;
  if (configure.paramMask & Config::USE_BRUNE_FORMALISM) compound->CalcShiftFunctions(configure);
  // One Jacobian row from an accumulator of d r / d(E, gamma), scattered to the packed columns.
  auto appendRow = [&](double r, const GradAccum &acc) {
    vector_r rowFull(pmap.NumFull(), 0.0);
    acc.Scatter(pmap, rowFull);
    residuals.push_back(r);
    const size_t base = jac.size();
    jac.resize(base + (size_t)nCols, 0.0);
    for (int f = 0; f < pmap.NumFull(); f++) {
      const int a = pmap.FullToPacked(f);
      if (a >= 0 && a < nCols) jac[base + (size_t)a] = rowFull[f];
    }
  };
  std::vector<double> q;
  LevelValues(compound, configure, set, q);
  const std::vector<Column> cols = HasLevelRows(set) ? LevelJacobian(compound, configure, set, pmap) : std::vector<Column>();
  for (std::size_t k = 0; k < set.rows.size(); k++) {
    const ModelDataRow &row = set.rows[k];
    if (!row.isScatteringLength) {
      // d r_k / d theta = (d q_k / d theta) / e_k.
      GradAccum acc;
      acc.Init(compound);
      for (const Column &c : cols) AddToAccum(c, c.dq[k] / row.datum.error, acc);
      appendRow((q[k] - row.datum.value) / row.datum.error, acc);
      continue;
    }
    // Exact rows from the level-matrix adjoint at the pair's threshold.
    vector_r rs;
    ScatteringLength::ConstraintResiduals(compound, configure, row.constraint, rs);
    for (size_t r = 0; r < rs.size(); r++) {
      GradAccum acc;
      acc.Init(compound);
      if (!ScatteringLength::ConstraintGradient(compound, configure, row.constraint, acc, shiftDeriv, (int)r)) return false;
      appendRow(rs[r], acc);
    }
  }
  return true;
}

std::vector<Value> SectionValues(CNuc *compound, const Config &configure, const ModelDataSet &set) {
  std::vector<Value> out;
  std::vector<double> q;
  LevelValues(compound, configure, set, q);
  for (std::size_t k = 0; k < set.rows.size(); k++) {
    const ModelDataRow &row = set.rows[k];
    if (!row.inSection) continue;
    Value v;
    v.datum = &row.datum;
    if (row.isScatteringLength) {
      const complex a = ScatteringLength::ConstraintValue(compound, configure, row.constraint);
      v.model = std::real(a);
      v.modelImag = std::imag(a);
      v.residual = (v.model - row.datum.value) / row.datum.error;
      v.chi2 = v.residual * v.residual;
      if (row.datum.hasImag) {
        v.residualImag = (v.modelImag - row.datum.imag) / row.datum.imagError;
        v.chi2 += v.residualImag * v.residualImag;
      }
    } else {
      v.model = q[k];
      v.residual = (v.model - row.datum.value) / row.datum.error;
      v.chi2 = v.residual * v.residual;
    }
    out.push_back(v);
  }
  return out;
}

std::string Report(CNuc *compound, const Config &configure, const ModelDataSet &set) {
  if (!set.hasSection) return "";
  std::ostringstream out;
  out << "Model data (<modelData>) at the parameters of this run: observed energies and widths, as parameters.out "
         "gives them (docs/upgrade/data_kinds.md section 4.7)"
      << std::endl;
  out << std::left << std::setw(20) << "name" << std::setw(20) << "quantity" << std::right << std::setw(18) << "model"
      << std::setw(18) << "value" << std::setw(14) << "error" << "  " << std::left << std::setw(8) << "unit" << std::right
      << std::setw(14) << "residual" << std::setw(14) << "chi2" << std::endl;
  double total = 0.0;
  for (const Value &v : SectionValues(compound, configure, set)) {
    const ex::ModelDatum &d = *v.datum;
    auto line = [&](const std::string &quantity, double model, double value, double error, double residual, double chi2) {
      out << std::left << std::setw(20) << d.name << std::setw(20) << quantity << std::right << std::scientific
          << std::setprecision(9) << std::setw(18) << model / d.unitFactor << std::setw(18) << value / d.unitFactor
          << std::setw(14) << std::setprecision(4) << error / d.unitFactor << "  " << std::left << std::setw(8) << d.unit
          << std::right << std::setw(14) << residual << std::setw(14) << chi2 << std::defaultfloat << std::endl;
    };
    if (d.quantity == ex::ModelQuantity::ScatteringLength && d.hasImag) {
      line("Re " + d.quantityName, v.model, d.value, d.error, v.residual, v.residual * v.residual);
      line("Im " + d.quantityName, v.modelImag, d.imag, d.imagError, v.residualImag, v.residualImag * v.residualImag);
    } else {
      line(d.quantityName, v.model, d.value, d.error, v.residual, v.chi2);
    }
    total += v.chi2;
  }
  out << "Total model-data chi2 (<modelData>): " << std::setprecision(10) << total << std::endl;
  return out.str();
}

}  // namespace ModelData
