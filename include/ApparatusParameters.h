// A number of an apparatus record that a fit may vary: the parameter kind `apparatus` (ParamKind::Apparatus,
// docs/upgrade/apparatus_parameters.md).  One table of rows, one per way a value reaches the calculation (the stopping
// scale: the targets' stopping; `kernel`: a table kernel's own parameter); every parameter of every kernel of the
// kernel table goes through the row `kernel`, so a new fitted kernel parameter needs no line here.
#ifndef APPARATUSPARAMETERS_H
#define APPARATUSPARAMETERS_H

#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

class Config;
class TargetEffect;
class WorkingCopy;
struct TargetEffectValues;

/// What one effect's record says of one of its fitted numbers.  Effects whose declarations have one `key` share one
/// parameter, which takes the first one's settings.
struct ApparatusDeclaration {
  int slot = 0;      ///< the row's place for the value (kernel: KernelParameters::values[slot])
  std::string key;   ///< what shares one parameter
  /// What it belongs to, in the order of its name (the values joined by '.'): ParameterRecord::owner.
  std::vector<std::pair<std::string, std::string>> owner;
  std::string unit;  ///< the value's unit
  double layoutBase = 0.0;  ///< what the limits and the layout are laid out from: the centre or the start (LayoutBase.h)
  double centre = 0.0;    ///< the prior's centre: as written, else the start (docs/upgrade/parameter_kinds.md section 4)
  bool vary = false;
  double sigma = 0.0;     ///< prior width, the value's unit; 0 = no prior
  double layout = 0.0;    ///< the largest value the sub-point grids span (the upper limit of a fitted one)
  std::string note;       ///< printed once at set-up when not empty (a default prior, say)
  std::string Name() const;
};

/// One row: how its effects declare their parameters, how a value reaches the working copy, what must follow.
struct ApparatusRow {
  const char *name;  ///< ParameterRecord::apparatus, "stopping_scale", "kernel"
  const char *label;     ///< the covariance files' label: "<label> <name> (<doc>)"
  const char *chiLabel;  ///< chiSquared.out: " Total-<chiLabel>-Chi-Squared: " when a parameter of the row has a prior
  const char *doc;
  /// The declarations of effect `te` (appended to `out`).
  void (*declare)(const TargetEffect &te, std::vector<ApparatusDeclaration> &out);
  /// Writes value v of slot `slot` into the working-copy values of effect `te`; true if the value changed.
  bool (*write)(const TargetEffect &te, int slot, double v, TargetEffectValues &values);
  /// The working-copy value of slot `slot` (the start value of the vector).
  double (*read)(const TargetEffect &te, int slot, const TargetEffectValues &values);
  /// After the values of the effects marked in `changed` (0-based) moved: what follows them (null: nothing).
  void (*follow)(WorkingCopy &, const std::vector<char> &changed, const Config &);
};

/// The rows, in the order of their parameters in the vector.  A constant table.
const std::vector<ApparatusRow> &ApparatusRows();

/// One fit parameter of kind `apparatus`: its row, the settings (its first effect's declaration), its effects.
struct ApparatusParam {
  int row = 0;                  ///< index in ApparatusRows()
  ApparatusDeclaration declaration;
  int lead = 0;                 ///< its first target effect (copy, 1-based)
  std::vector<int> effects;     ///< every target effect (copy, 1-based) of the parameter
  bool priorReplaced = false;   ///< the prior covariance holds its prior (Objective, ReplacedPriors)
  const ApparatusRow &Row() const { return ApparatusRows()[row]; }
  double Penalty(double value) const;
};

/// The rows' own declarations (their files): the stopping scale (StoppingScale.cpp), a table kernel's parameter.
ApparatusRow StoppingScaleRow();
ApparatusRow KernelParameterRow();

class AZUREParams;
class DataSet;
class ParamIndexMap;
/// The apparatus parameters' slots of the vector (WorkingCopy::FillMnParams), after every other kind.
void AddApparatusParams(WorkingCopy &data, AZUREParams &params);
/// The chiSquared.out lines of the rows with a fitted prior (WorkingCopy::WriteChiSquared).
void WriteApparatusChiSquared(const WorkingCopy &data, std::ostream &out);
/// The covariance files' label of apparatus parameter k (ParameterLabel.cpp).
std::string ApparatusLabel(const DataSet &data, int k);
/// The apparatus block of the index map, from full index \p i on (advanced past it); no analytic derivative: the
/// parameters are differenced on the whole chi^2, as a density (AZURECalc.cpp).
void AppendApparatusParams(ParamIndexMap &map, WorkingCopy *data, int &i);

#endif
