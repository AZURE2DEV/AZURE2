// AZURECalcMCMC, the priors: the diagonal Gaussian priors given by the caller or built from the project
// (SetPriors, BuildAutoPriors), the prior covariance (SetUpPriorCovariance), what a run needs prepared before its
// first evaluation (PrepareRun), and the log-prior itself with its normalisation (CalculateLogPrior,
// LogPriorNormalisation).  One unit because the set-up and the evaluation must agree on which parameter has
// which prior and on the posterior's support.
#include "Objective.h"
#include "AZURECalcMCMC.h"
#include "ModelPass.h"
#include "ParameterLabel.h"
#include "Config.h"
#include "CompoundSystem.h"
#include "LevelReport.h"
#include "WorkingCopy.h"
#include "ParameterLimitsManager.h"
#include "AZUREParams.h"
#include "PriorCovarianceFile.h"
#include "GSLException.h"
#include "ALevel.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include <vector>
#include <limits>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void AZURECalcMCMC::SetPriors(const std::vector<double> &priorMeans,
                              const std::vector<double> &priorStds,
                              const std::vector<bool> &usePriors) {
  priorMeans_ = priorMeans;
  priorStds_ = priorStds;
  usePriors_ = usePriors;
}

void AZURECalcMCMC::BuildAutoPriors(bool useRWA) {
  if (!parametersInitialized_) {
    configure().outStream << "Warning: BuildAutoPriors() called before parameters were "
                          << "initialized; automatic priors skipped.\n";
    return;
  }

  // Every parameter of the vector, in its length and order, by its kind (its record, written by CompoundSystem::FillMnParams
  // and WorkingCopy::FillMnParams).  (Re-walking the levels and segments, as this once did, gave an entry to a free THM
  // norm, which is profiled and not a parameter; matching names duplicated the fit's prior rules.)
  std::vector<ParamKind> allKinds = paramKinds_;
  std::vector<double> allAutoMean;
  std::vector<double> allAutoStd;  // <= 0 means "no automatic prior available"

  // The diagonal priors from the fit's own enumeration (ForEachDiagonalPrior, DiagonalPriors.h; docs/upgrade/objective.md
  // step 5: the MCMC shares the fit's configuration), by index, not by name.
  const ParamIndexMap pmap = BuildParamIndexMap(compound(), data(), std::vector<bool>());
  std::map<int, std::pair<double, double>> diagonal;  // full index -> (mean, sigma)
  ForEachDiagonalPrior(data(), pmap,
                       [&](int idx, double nominal, double sigma) { diagonal[idx] = std::make_pair(nominal, sigma); });
  for (size_t i = 0; i < allKinds.size(); ++i) {
    const ParamKind kind = allKinds[i];
    auto it = diagonal.find((int)i);
    allAutoMean.push_back(it == diagonal.end() ? 0.0 : it->second.first);
    // A prior kind without a prior (no quoted error) is 0, sampled uniformly; levels and cbkg -1.
    allAutoStd.push_back(it != diagonal.end() ? it->second.second : KindRow(kind).ownPrior ? 0.0 : -1.0);
  }

  if (allKinds.size() != fixed_.size()) {
    configure().outStream << "Warning: parameter classification produced "
                          << allKinds.size() << " entries for "
                          << fixed_.size() << " parameters";
    configure().outStream << "; automatic priors skipped to avoid mis-assigning them.\n";
    freeKinds_.clear();
    return;
  }

  // Compress to the varying parameters, which is how the priors and the
  // sampled vector are indexed.
  freeKinds_.clear();
  std::vector<double> freeAutoMean;
  std::vector<double> freeAutoStd;
  std::vector<int> freeAllIndex;  // varying index -> index among all parameters
  for (size_t i = 0; i < fixed_.size(); ++i) {
    if (!fixed_[i]) {
      freeKinds_.push_back(allKinds[i]);
      freeAutoMean.push_back(allAutoMean[i]);
      freeAutoStd.push_back(allAutoStd[i]);
      freeAllIndex.push_back((int)i);
    }
  }

  const size_t nfree = freeKinds_.size();

  // Any priors the caller set cover the varying parameters; grow them to the
  // full length so the loop below can index them unconditionally.
  priorMeans_.resize(nfree, 0.0);
  priorStds_.resize(nfree, 0.0);
  usePriors_.resize(nfree, false);

  int nAutoNorm = 0, nAutoShift = 0, nNoError = 0, nUserNeeded = 0, nUserSet = 0, nRecord = 0;
  int nCbkg = 0, nCbkgSet = 0;
  std::vector<std::string> noErrorLabels;
  std::vector<std::string> noPriorLabels;

  for (size_t i = 0; i < nfree; ++i) {
    const ParamKind kind = freeKinds_[i];

    if (KindRow(kind).ownPrior && IsRMatrixParameter(kind)) {
      // A number of a gamma_decay record, an isospin shift or mix: its prior is a term of the objective (the
      // record's; the isospin term, with the soft ties), which the likelihood evaluates (Objective::SamplerChi2); no
      // Gaussian log-prior here.
      usePriors_[i] = false;
      nRecord++;
    } else if (KindRow(kind).ownPrior) {
      // Derived from the data, so it always wins over whatever the GUI holds.
      if (freeAutoStd[i] > 0.0) {
        priorMeans_[i] = freeAutoMean[i];
        priorStds_[i] = freeAutoStd[i];
        usePriors_[i] = true;
        if (kind == ParamKind::Norm)
          nAutoNorm++;
        else
          nAutoShift++;
      } else {
        // No experimental error quoted: nothing to build a prior from, so the
        // parameter stays uniform rather than getting an invented width.
        usePriors_[i] = false;
        nNoError++;
        noErrorLabels.push_back(AZURELabel::Parameter(compound(), levels_.inputs, data(), freeAllIndex[i]));
      }
    } else if (kind == ParamKind::ThmCoherent) {
      // THM coherent background: nothing in the data quotes an error on it,
      // so no automatic prior; a prior the caller defined is kept.
      nCbkg++;
      if (usePriors_[i] && priorStds_[i] > 0.0) nCbkgSet++;
    } else {
      // Level energies and widths: the user's business.
      nUserNeeded++;
      if (usePriors_[i] && priorStds_[i] > 0.0)
        nUserSet++;
      else
        noPriorLabels.push_back(AZURELabel::Parameter(compound(), levels_.inputs, data(), freeAllIndex[i]));
    }
  }

  // Nuisance priors (<parameterSettings> useAsNuisance) on level parameters and THM coherent backgrounds, as the fit
  // applies them (docs/upgrade/objective.md D6, step 5): reduced amplitudes in an RWA chain, the physical width in a
  // physical-parameter chain.  They are part of the shared configuration, so they replace a prior the MCMC
  // configuration set on the same parameter (noted).
  {
    AZUREParams scratch;
    compound()->FillMnParams(scratch, levels_.state, levels_.inputs, &configure());
    data()->FillMnParams(scratch);
    ParameterLimitsManager limits(&configure(), compound(), levels_.state, levels_.inputs, data(), &scratch);
    int nNuisance = 0;
    std::string overridden;
    if (limits.ReadParameterSettings()) {
      limits.ApplyAllParameterSettings(scratch);
      for (size_t i = 0; i < nfree; ++i) {
        const ParamKind kind = freeKinds_[i];
        if (!IsLevelSlot(kind) && kind != ParamKind::ThmCoherent) continue;
        double nominal, sigma;
        const bool has = useRWA ? limits.NuisancePrior(freeAllIndex[i], nominal, sigma)
                                : limits.NuisancePriorPhysical(freeAllIndex[i], nominal, sigma);
        if (!has) continue;
        const std::string label = AZURELabel::Parameter(compound(), levels_.inputs, data(), freeAllIndex[i]);
        if (usePriors_[i] && priorStds_[i] > 0.0) {
          overridden += (overridden.empty() ? "" : ", ") + label;
        } else if (kind == ParamKind::ThmCoherent) {
          nCbkgSet++;
        } else {  // the summary below counted it as a level parameter without a prior
          nUserSet++;
          noPriorLabels.erase(std::remove(noPriorLabels.begin(), noPriorLabels.end(), label), noPriorLabels.end());
        }
        priorMeans_[i] = nominal;
        priorStds_[i] = sigma;
        usePriors_[i] = true;
        nNuisance++;
      }
    }
    if (nNuisance > 0)
      configure().outStream << "Nuisance priors from <parameterSettings>: " << nNuisance << " level/background parameter"
                            << (nNuisance == 1 ? "" : "s") << ", as in the fit\n";
    if (!overridden.empty())
      configure().Note("Note: MCMC: the <parameterSettings> nuisance prior replaces the MCMC prior set on " + overridden +
                       " (the MCMC shares the fit's configuration).");
  }

  configure().outStream << "Automatic priors: " << nAutoNorm << " normalization"
                        << (nAutoNorm == 1 ? "" : "s") << ", " << nAutoShift
                        << " energy shift" << (nAutoShift == 1 ? "" : "s")
                        << " taken from the quoted experimental errors\n";
  if (nRecord > 0)
    configure().outStream << "  " << nRecord << " number" << (nRecord == 1 ? "" : "s")
                          << " of gamma_decay or isospin records: their priors (the remainder's wall, the soft ties) "
                          << "are terms of the objective, sampled inside their limits\n";
  if (nNoError > 0) {
    configure().outStream << "  " << nNoError << " normalization/energy-shift "
                          << "parameter" << (nNoError == 1 ? " has" : "s have")
                          << " no quoted error and " << (nNoError == 1 ? "is" : "are")
                          << " sampled with a uniform prior:\n";
    for (size_t i = 0; i < noErrorLabels.size(); ++i)
      configure().outStream << "      " << noErrorLabels[i] << "\n";
  }
  configure().outStream << "  " << nUserSet << " of " << nUserNeeded
                        << " level-energy/width parameters have a user-defined prior\n";
  if (nUserSet < nUserNeeded) {
    configure().outStream << "  Warning: " << (nUserNeeded - nUserSet)
                          << " level-energy/width parameter"
                          << ((nUserNeeded - nUserSet) == 1 ? " is" : "s are")
                          << " sampled with an unbounded uniform prior. Define "
                          << "priors for them in the MCMC tab if the chain wanders:\n";
    // Cap the list so a large model cannot bury the rest of the log.
    const size_t maxListed = 10;
    for (size_t i = 0; i < noPriorLabels.size() && i < maxListed; ++i)
      configure().outStream << "      " << noPriorLabels[i] << "\n";
    if (noPriorLabels.size() > maxListed)
      configure().outStream << "      ...and " << (noPriorLabels.size() - maxListed)
                            << " more\n";
  }
  if (nCbkg > 0)
    configure().outStream << "  " << nCbkg << " THM coherent-background parameter"
                          << (nCbkg == 1 ? "" : "s") << ": " << nCbkgSet
                          << " with a user-defined prior, the others uniform "
                          << "(no automatic prior)\n";
  configure().outStream.flush();
}

bool AZURECalcMCMC::SetUpPriorCovariance(bool useRWA) {
  useRWA_ = useRWA;
  priorCov_.reset();
  if (!(configure().paramMask & Config::CALCULATE_WITH_DATA) || configure().priorCovariance.empty()) return true;
  const int basis = AZUREParams::BasisForConfig(configure());
  if (LoadPriorCovariance(configure(), paramNames_, fixed_, basis, priorCov_) < 0) return false;
  if (!useRWA) {
    // The prior is in the param.sav space: level energies and reduced-width amplitudes.  A physical-parameter chain
    // samples partial widths, so a prior on a level parameter would be applied to the wrong quantity.
    for (size_t a = 0; a < priorCov_->Names().size(); a++) {
      const std::string &name = priorCov_->Names()[a];
      const int f = priorCov_->FullIndex()[a];
      const bool level = f >= 0 && f < (int)paramKinds_.size() &&
                         IsLevelSlot(paramKinds_[f]);
      if (level) {
        configure().outStream << "**ERROR: priorCovariance: " << name << " is a level parameter and its prior is on "
                              << "the reduced-width (RWA) parameters of param.sav; sample with RWA parameters."
                              << std::endl;
        priorCov_.reset();
        return false;
      }
    }
  }
  // The multivariate Gaussian replaces the diagonal prior of each of its parameters.
  std::vector<int> packedOf(fixed_.size(), -1);
  for (int i = 0, k = 0; i < (int)fixed_.size(); i++)
    if (!fixed_[i]) packedOf[i] = k++;
  std::string replaced;
  for (size_t a = 0; a < priorCov_->FullIndex().size(); a++) {
    const int f = priorCov_->FullIndex()[a];
    const int k = (f >= 0 && f < (int)packedOf.size()) ? packedOf[f] : -1;
    if (k >= 0 && k < (int)usePriors_.size() && usePriors_[k] && k < (int)priorStds_.size() && priorStds_[k] > 0.0) {
      usePriors_[k] = false;
      replaced += (replaced.empty() ? "" : ", ") + priorCov_->Names()[a];
    }
  }
  if (!replaced.empty())
    configure().Note("Note: prior covariance: its multivariate Gaussian replaces the diagonal prior of " + replaced +
                     " in the MCMC run (a fit applies both).");
  return true;
}

bool AZURECalcMCMC::PrepareRun(const std::vector<double> &initialParams, bool useRWA) {
  // Isospin levels (IsospinTies.h): the ties act on the fit's basis, and the transformation from observed widths
  // mixes the levels of a J^pi group, so a chain on the observed parameters cannot hold them.
  if (!useRWA && compound()->GetIsospin().Active()) {
    configure().outStream << "**ERROR: a project with isospin levels is sampled in the fit's basis (the level energies and "
                          << "reduced-width amplitudes): switch the MCMC parameter type to RWA." << std::endl;
    return false;
  }
  UpdateParameterVectors(initialParams);
  // Classify the varying parameters and derive the normalization and energy-shift priors from the data. Must follow
  // UpdateParameterVectors() (it needs fixed_) and precede any log-probability evaluation.
  BuildAutoPriors(useRWA);
  // The prior parameter covariance (after the automatic priors, which it replaces for its parameters).
  return SetUpPriorCovariance(useRWA);
}

double AZURECalcMCMC::LogPriorNormalisation() const {
  double c = (priorCov_ && priorCov_->Size() > 0) ? priorCov_->LogNormalization() : 0.0;
  for (size_t i = 0; i < usePriors_.size() && i < priorStds_.size(); i++)
    if (usePriors_[i] && priorStds_[i] > 0.0) c -= log(priorStds_[i] * sqrt(2.0 * M_PI));
  return c;
}

double AZURECalcMCMC::CalculateLogPrior(const std::vector<double> &params) const {
  // Prior parameter covariance: -0.5 (p - mu)^T C^-1 (p - mu) - 0.5 log det(2 pi C).
  double covLogPrior = 0.0;
  if (priorCov_ && priorCov_->Size() > 0) {
    const vector_r full = useRWA_ ? ReconstructFullParameters(params) : ReconstructFullParametersPhysical(params);
    covLogPrior = -0.5 * priorCov_->Chi2(full) + priorCov_->LogNormalization();
  }
  if (priorMeans_.empty() || priorStds_.empty() || usePriors_.empty()) {
    return covLogPrior;  // No diagonal priors specified - uniform priors
  }

  double logPrior = covLogPrior;

  for (size_t i = 0; i < params.size() && i < usePriors_.size(); i++) {
    if (usePriors_[i] && i < priorMeans_.size() && i < priorStds_.size()) {
      double mean = priorMeans_[i];
      double std = priorStds_[i];

      if (std > 0.0) {
        // Gaussian prior: log(prior) = -0.5 * ((x - mean) / std)^2 - log(std * sqrt(2*pi))
        double deviation = (params[i] - mean) / std;
        logPrior += -0.5 * deviation * deviation - log(std * sqrt(2.0 * M_PI));
      }
    }
  }

  return logPrior;
}
