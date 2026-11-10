#ifndef MODELDATA_H
#define MODELDATA_H

#include <memory>
#include <string>
#include <vector>

#include "Constants.h"
#include "ScatteringLength.h"
#include "azure2/experiment/ModelDatum.h"

class CNuc;
class Config;
class ParamIndexMap;
struct GradAccum;

/*!
 * The model data of a project: values compared with a function of the level parameters rather than with a model at a
 * beam energy (docs/upgrade/data_model.md section 4, decision M10; docs/upgrade/project_format3.md section 3.10).
 * Each datum gives one or two residual rows (q_k(theta) - d_k)/e_k; together they are the one term T9 of the
 * objective (Objective.cpp, ModelDataTerm), which every path -- the fit, the gradient, the least-squares rows, the API,
 * the MCMC, the writer -- reaches through Objective::Engine().
 *
 * The rows, in this order:
 *  - the measured scattering lengths of the <scatteringLengths> section, line by line (the section and its grammar
 *    unchanged, include/ScatteringLength.h): Re then Im of each line, each when its sigma is > 0;
 *  - the `datum` records of the <modelData> section (format-3 records, io::ReadRecords and Schema::Validate, resolved
 *    against <levels> by experiment::BuildModelData), in file order: one row each, two for a scattering length with
 *    an imaginary part.
 * A scattering length is evaluated by ScatteringLength (exact gradient through the level-matrix adjoint).  A level
 * quantity (resonance strength, energy, widths, branching ratio, ANC) is evaluated from the observed energies and
 * widths of the current fit values (CNuc::ObservedOfFitValues, the quantities parameters.out reports), never the
 * formal ones (data_kinds.md section 4.7); its gradient is a central difference over the level parameters it depends
 * on (ModelData.cpp, kEnergyStep, kAmplitudeStep).
 *
 * Setup state: read once by CNuc::Fill (ModelData::Read), immutable afterwards, shared by the working copies.
 */
struct ModelDataRow {
  /// True: a measured scattering length, `constraint`; false: a level quantity of `datum`.
  bool isScatteringLength = true;
  ScatteringLength::Constraint constraint;
  /// The datum as read (a <scatteringLengths> line has its line number and quantity only).
  azure2::experiment::ModelDatum datum;
  /// From the <modelData> section (reported in modelData.out), not from <scatteringLengths>.
  bool inSection = false;
};

struct ModelDataSet {
  /// The <scatteringLengths> lines as read, for scatteringLengths.out (ScatteringLength::Report).
  std::vector<ScatteringLength::Constraint> scatteringLengthLines;
  /// The data that give residual rows, in row order.
  std::vector<ModelDataRow> rows;
  /// The project has a <modelData> section (modelData.out and its total are written only then).
  bool hasSection = false;
  /// Number of residual rows.
  int NumResiduals() const;
  /// Number of residual rows of the <modelData> section's data.
  int NumSectionResiduals() const;
};

namespace ModelData {

/// Read the model data of the project: the <scatteringLengths> and <modelData> sections, against the levels and
/// pairs of `compound` (after <levels> is read).  Prints the warnings; nullptr (with every error printed) on an error.
std::shared_ptr<const ModelDataSet> Read(CNuc &compound, const Config &configure);

/// The residual rows at the compound's current fit values, appended to r in row order.
void Residuals(CNuc *compound, const Config &configure, const ModelDataSet &set, vector_r &r);

/// The chi-squared of the rows, sum r_k^2 in row order (0 without rows); `sectionPart` and `linesPart`, when given,
/// receive the parts of the <modelData> rows and of the <scatteringLengths> lines, each summed in row order.
double Chi2(CNuc *compound, const Config &configure, const ModelDataSet &set, double *sectionPart = nullptr,
            double *linesPart = nullptr);

/*!
 * d chi2 / d(E, gamma) of the level parameters into accum: the scattering lengths exactly, through the level-matrix
 * adjoint at each threshold (ScatteringLength::Chi2Gradient); the level quantities by central differences over the
 * parameters of their J^pi groups (all groups with external capture, whose amplitudes depend on the final states),
 * those with a packed column in pmap.  False where the adjoint does not support the configuration (the caller then
 * differentiates numerically).
 */
bool Gradient(CNuc *compound, const Config &configure, const ModelDataSet &set, const ParamIndexMap &pmap,
              const vector_matrix_r *shiftDeriv, GradAccum &accum);

/*!
 * The residual rows with their Jacobian rows over nCols packed columns (pmap's), appended to residuals and jac.
 * Brune and Park: the shift functions are brought up to date first (CNuc::CalcShiftFunctions).  False, with the rows
 * appended so far, where Gradient would return false.
 */
bool AppendRows(CNuc *compound, const Config &configure, const ModelDataSet &set, const ParamIndexMap &pmap,
                const vector_matrix_r *shiftDeriv, int nCols, vector_r &residuals, vector_r &jac);

/// One <modelData> datum at the current fit values: model value (Re and Im for a scattering length) and residuals.
struct Value {
  const azure2::experiment::ModelDatum *datum = nullptr;
  double model = 0.0;      ///< in the engine unit (MeV, 1, fm, fm^-1/2); Re a0 for a scattering length
  double modelImag = 0.0;  ///< Im a0 (fm)
  double residual = 0.0;   ///< (model - value) / error
  double residualImag = 0.0;
  double chi2 = 0.0;       ///< the datum's share of the term
};
/// Every datum of the <modelData> section, in file order.
std::vector<Value> SectionValues(CNuc *compound, const Config &configure, const ModelDataSet &set);

/// The table of modelData.out and of the run log: one line per datum (name, quantity, model value, datum, error in
/// the unit written, residual, chi2), and the total.  Empty without a <modelData> section.
std::string Report(CNuc *compound, const Config &configure, const ModelDataSet &set);

}  // namespace ModelData

#endif
