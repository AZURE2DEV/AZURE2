#ifndef TARGETSTOPPING_H
#define TARGETSTOPPING_H

#include <memory>
#include <string>
#include <vector>

#include "Equation.h"
#include "azure2/experiment/Diagnostics.h"
#include "azure2/experiment/SetupSpecs.h"
#include "azure2/stopping/Material.h"
#include "azure2/stopping/Tables.h"

class Config;

/*!
 * The stopping of one target, the one place the engine asks for it (docs/upgrade/energy_grid.md section 6.1):
 * epsilon_active(E), lab MeV cm^2 per atom of the active element, at lab energy E (MeV) of the projectile, from one of
 *   - formula:    the record's expression (an Equation), already per active nucleus;
 *   - evaluation: SRIM-2013's stopping of the composition (stopping/'s Stopping(): Bragg's rule with SRIM's phase factor
 *                 and compound correction; cite SRIM, stopping::kSrimCitation);
 *   - table:      a user table of the composition's stopping (stopping/'s UserStopping, read by ReadStoppingTable).
 * For the last two the convention of SetupSpecs.h StoppingSpec: epsilon_active = 1e-21 eps_atom / f_active, eps_atom in
 * eV/(1e15 atoms/cm^2) per atom of the material and f_active the atom fraction of the active element.  The stopping
 * scale s is not applied here: TargetEffect::StoppingAt multiplies it, once.
 */
class TargetStopping {
 public:
  using Source = azure2::experiment::StoppingSpec::Source;
  TargetStopping() = default;
  /// The stopping of a spec: \p formula and \p parameters are the formula's (stoppingPowerEquation, a0, a1, ...).  A
  /// table that cannot be read leaves Error() set; nothing here throws.
  TargetStopping(const azure2::experiment::StoppingSpec &, const std::string &formula, const std::vector<double> &parameters,
                 const Config &);
  /// The library part alone (evaluation, table), for --validate and --describe, which have no Config: a formula spec
  /// gives a stopping that is not to be evaluated (only Provenance and Check are meaningful).  A table file is
  /// relative to \p projectDirectory (azure2::experiment::PathToOpen).
  TargetStopping(const azure2::experiment::StoppingSpec &, const std::string &projectDirectory);
  Source GetSource() const { return source_; }
  bool IsLibrary() const { return source_ != Source::Formula; }
  /// epsilon_active at lab energy \p eLab (MeV), lab MeV cm^2.  Outside the library's range (or without SRIM's data) a
  /// quiet NaN, which the set-up refuses with its reason (Check).
  double PerActiveNucleus(double eLab, const Config &) const;
  /// The same for an evaluation or a table (no Config needed); NaN for a formula.
  double LibraryPerActiveNucleus(double eLab) const;
  /// The reason the stopping cannot be evaluated at \p eLab, "" when it can (the library's exception text).  Its
  /// diagnostic code: kStoppingData when \p dataMissing is set on return, else kStoppingRange.
  std::string Check(double eLab, bool &dataMissing) const;
  /// Why a table could not be read ("" when it was, and for a formula or an evaluation).
  const std::string &Error() const { return error_; }
  /// The atom fraction of the active element in the composition (the divisor of the convention).
  double ActiveFraction() const { return activeFraction_; }
  const azure2::stopping::Material &GetMaterial() const { return material_; }
  const azure2::stopping::Projectile &GetProjectile() const { return projectile_; }
  /// The relative 1-sigma uncertainty the library states for this projectile class (stopping.md section 7): the
  /// default prior width of a fitted scale of an evaluation; 0 for a formula or a table.
  double StatedRelUncertainty() const;
  /// Where the stopping comes from: the library's provenance text (SRIM-2013 and its citation, LunaSRIM, the data
  /// source, phase, compound correction, per element), the table's file, or the formula; with the active fraction.
  std::string Provenance() const;

 private:
  Source source_ = Source::Formula;
  Equation formula_;
  std::string formulaText_;
  azure2::stopping::Material material_;
  azure2::stopping::Projectile projectile_;
  double activeFraction_ = 1.0;
  std::string error_;
  /// SRIM's tables of the run's directory (Config::environment.srimDataDirectory: AZURE2_SRIM_DATA_DIR), for an
  /// evaluation made with a Config; null for the stopping of --validate and --describe, which take the library's
  /// SrimTables().  srimError_: why they could not be read.
  std::shared_ptr<const azure2::stopping::SrimData> srim_;
  std::string srimError_;
  void SetUpLibrary(const azure2::experiment::StoppingSpec &, const std::string &projectDirectory);
};

/*!
 * --validate and --describe (ProjectCommands.cpp): for every target effect of a built project whose stopping comes from
 * the library, its source (Provenance, with SRIM's citation) as text, and a diagnostic (StoppingChecks.h) when the
 * library cannot give it at all (SRIM's data files missing, a table unreadable).  \p effects are the built specs, their
 * table files relative to \p projectDirectory.
 */
struct StoppingReport {
  std::string text;
  std::vector<azure2::experiment::Diagnostic> refusals;
};
StoppingReport ReportTargetStoppings(const std::vector<azure2::experiment::TargetEffectSpec> &effects,
                                     const std::string &projectDirectory);

#endif
