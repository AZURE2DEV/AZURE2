// The stopping of a target: formula, SRIM-2013 evaluation or user table (TargetStopping.h).
#include "TargetStopping.h"

#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "Config.h"
#include "azure2/experiment/ProjectPaths.h"
#include "azure2/experiment/StoppingChecks.h"
#include "azure2/physics/Elements.h"
#include "azure2/stopping/Electronic.h"
#include "azure2/stopping/Stopping.h"
#include "azure2/stopping/Tables.h"

TargetStopping::TargetStopping(const azure2::experiment::StoppingSpec &spec, const std::string &formula,
                               const std::vector<double> &parameters, const Config &configure)
    : source_(spec.source), formulaText_(formula) {
  if (source_ == Source::Formula) {
    formula_.Initialize(formula, (int)parameters.size(), configure);
    for (std::size_t i = 0; i < parameters.size(); i++) formula_.SetParameter((unsigned)i, parameters[i], configure);
    return;
  }
  this->SetUpLibrary(spec, configure.baseDirectory);
}

TargetStopping::TargetStopping(const azure2::experiment::StoppingSpec &spec, const std::string &projectDirectory)
    : source_(spec.source) {
  if (source_ != Source::Formula) this->SetUpLibrary(spec, projectDirectory);
}

void TargetStopping::SetUpLibrary(const azure2::experiment::StoppingSpec &spec, const std::string &projectDirectory) {
  std::vector<azure2::stopping::Constituent> constituents;
  double total = 0.0, active = 0.0;
  for (const auto &zn : spec.composition) {
    constituents.push_back({zn.first, zn.second});
    total += zn.second;
    if (zn.first == spec.activeZ) active += zn.second;
  }
  try {
    material_ = azure2::stopping::MaterialFromFractions(spec.material, constituents);
  } catch (const std::exception &e) {
    error_ = e.what();
    return;
  }
  if (spec.phase == "gas") material_.phase = azure2::stopping::Phase::kGas;
  if (spec.phase == "solid") material_.phase = azure2::stopping::Phase::kSolid;
  material_.compound_correction = spec.compoundCorrection;
  // the atom fraction of the active element: n_active / sum n_i (the material's fractions are normalised the same way)
  activeFraction_ = (total > 0.0) ? active / total : 0.0;
  projectile_ = azure2::stopping::Ion(spec.projectileZ, spec.projectileMass,
                                      std::to_string(spec.projectileA) + azure2::physics::ElementSymbol(spec.projectileZ));
  if (source_ == Source::Table) {
    std::ifstream in(azure2::experiment::PathToOpen(projectDirectory, spec.tableFile));
    if (!in) {
      error_ = "cannot open the stopping table " + spec.tableFile;
      return;
    }
    std::ostringstream text;
    // the units are the record's; the file holds the two columns
    text << "E[" << spec.tableEnergyUnit << "] ";
    text << (spec.tableStoppingUnit == "MeV*cm2/g" ? "S[MeV*cm2/g]" : "eps[eV/1e15atoms/cm2]") << "\n" << in.rdbuf();
    try {
      material_.user_stopping.push_back(
          azure2::stopping::ReadStoppingTable(text.str(), material_, projectile_, 0.0, "the table " + spec.tableFile));
    } catch (const std::exception &e) {
      error_ = std::string(e.what()) + " (" + spec.tableFile + ")";
    }
  }
}

double TargetStopping::PerActiveNucleus(double eLab, const Config &configure) const {
  if (source_ == Source::Formula) return formula_.Evaluate(configure, eLab);
  return this->LibraryPerActiveNucleus(eLab);
}

double TargetStopping::LibraryPerActiveNucleus(double eLab) const {
  if (source_ == Source::Formula || !error_.empty()) return std::numeric_limits<double>::quiet_NaN();
  try {
    return azure2::stopping::kStoppingToMeVcm2 * azure2::stopping::Stopping(material_, projectile_, eLab).total / activeFraction_;
  } catch (const std::exception &) {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

std::string TargetStopping::Check(double eLab, bool &dataMissing) const {
  dataMissing = false;
  if (source_ == Source::Formula) return "";
  if (!error_.empty()) {
    dataMissing = true;
    return error_;
  }
  try {
    if (source_ == Source::Evaluation) azure2::stopping::SrimTables();
  } catch (const std::exception &e) {
    dataMissing = true;
    return e.what();
  }
  try {
    azure2::stopping::Stopping(material_, projectile_, eLab);
  } catch (const std::exception &e) {
    return e.what();
  }
  return "";
}

double TargetStopping::StatedRelUncertainty() const {
  return source_ == Source::Evaluation ? azure2::stopping::ElectronicRelUncertainty(projectile_.z) : 0.0;
}

std::string TargetStopping::Provenance() const {
  std::ostringstream o;
  if (source_ == Source::Formula) {
    o << "  user formula " << formulaText_ << " (lab MeV cm^2 per active nucleus)\n";
    return o.str();
  }
  o << "  per atom of the active element: the material's stopping per atom / " << activeFraction_
    << " (the active element's atom fraction)\n";
  try {
    std::istringstream text(azure2::stopping::StoppingProvenance(material_, projectile_));
    for (std::string line; std::getline(text, line);) o << "  " << line << "\n";
  } catch (const std::exception &e) {
    o << "  (no provenance: " << e.what() << ")\n";
  }
  return o.str();
}

StoppingReport ReportTargetStoppings(const std::vector<azure2::experiment::TargetEffectSpec> &effects,
                                     const std::string &projectDirectory) {
  StoppingReport r;
  std::vector<std::string> seen;
  for (const azure2::experiment::TargetEffectSpec &e : effects) {
    if (!e.targetIntegration || e.stopping.source == TargetStopping::Source::Formula) continue;
    const TargetStopping stopping(e.stopping, projectDirectory);
    const std::string label = e.targetName + " for " + stopping.GetProjectile().label;
    bool known = false;
    for (const std::string &s : seen) known = known || s == label;
    if (known) continue;
    seen.push_back(label);
    bool dataMissing = false;
    const std::string why = stopping.Check(1.0, dataMissing);
    if (dataMissing) {
      r.refusals.push_back(azure2::experiment::StoppingDataMissing(e.targetName, why));
      continue;
    }
    r.text += "Stopping of target " + label + ":\n" + stopping.Provenance();
  }
  return r;
}
