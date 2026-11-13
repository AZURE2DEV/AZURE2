/*!
 * The ensemble sampler samples the fit's objective and owns its state (docs/upgrade/objective.md section 6).
 *
 * (a) The 15O project of tests/gamma_decay (models/record.azr: the 14N(p,g)15O level at 259.7 keV written as a
 *     gamma_decay record), with Gaussian priors on Gamma_gamma, on the branching to 6.1763 MeV and on the remainder
 *     (the record's own priors, term T11 of the objective), in Brune's and Park's bases, sampled in the fit's basis
 *     (reduced-width amplitudes) and in the observed one: at three points the sampler's log-probability equals
 *     -1/2 x the fit's objective (AZURECalc::operator(), the limits manager as the CLI's fit) plus the sampler's
 *     documented constant (LogPriorNormalisation; rel. 1e-10); with a branching below 0 or above 1, or Gamma_gamma
 *     below 0, it is -infinity (three points: start, displaced, outside); a short chain (8 walkers, 60 steps) of the
 *     branching to 6.8594 MeV alone (the other numbers fixed) under a prior centred at -0.01 (width 0.00001), which
 *     pulls it below its limit 0, stays inside [0, 1].
 * (b) The sampler moves with a term of the objective it did not have before (T11): the log-probability of the
 *     project with the record's priors minus that of the project without them is -1/2 x the difference of the fit's
 *     objectives at the same point, and not 0.
 * (c) Two samplers in one process, on two projects (the record project with its priors, observed basis, and
 *     without them, fit basis), run at once in
 *     two threads, write the chains (samples.mcmc) of the same runs made one after the other, byte for byte.
 * And T2 of tests/upgrade (normalisation, shift, density, abundance priors and a prior covariance), fit basis: the
 * identity of (a) at three points (the constant is the Gaussian priors' normalisation).
 *
 * Run from the directory CMake prepares (data/ with both projects' data, output/, checks/):
 *   sampler_objective_test <record.azr> <T2.azr>   (ctest: sampler_objective)
 */
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <thread>
#include <vector>

#include "AZURECalc.h"
#include "AZURECalcMCMC.h"
#include "AZUREParams.h"
#include "CNuc.h"
#include "CaptureAyRequests.h"
#include "Config.h"
#include "Objective.h"
#include "ParameterLimitsManager.h"
#include "PriorCovariance.h"
#include "WorkingCopy.h"

namespace {

int failures = 0;

void check(const std::string &name, bool ok, const std::string &detail = "") {
  std::printf("  %s  %s%s%s\n", ok ? "ok  " : "FAIL", name.c_str(), ok || detail.empty() ? "" : "  -- ",
              ok ? "" : detail.c_str());
  if (!ok) failures++;
}

std::string num(double x) {
  std::ostringstream s;
  s.precision(14);
  s << x;
  return s.str();
}

std::string readFile(const std::string &path) {
  std::ifstream in(path.c_str());
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

// `text` with every `from` replaced by `to`; false (and text unchanged) when `from` is not there.
bool replaceAll(std::string &text, const std::string &from, const std::string &to) {
  size_t pos = text.find(from);
  if (pos == std::string::npos) return false;
  for (; pos != std::string::npos; pos = text.find(from, pos + to.size())) text.replace(pos, from.size(), to);
  return true;
}

// A project as the CLI sets it up for a fit and for the sampler; its output into `outdir`.
struct Model {
  std::ostringstream log;  // the engine's messages for this model (not checked), one stream per thread
  Config cfg;
  CNuc *compound = nullptr;
  WorkingCopy *data = nullptr;
  AZUREParams params;
  std::unique_ptr<ParameterLimitsManager> limits;
  std::shared_ptr<const PriorCovariance> priorCov;
  Model(const std::string &azr, bool park, const std::string &outdir) : cfg(log) {
    cfg.configfile = azr;
    ok = cfg.ReadConfigFile() >= 0;
    cfg.paramMask &= ~(Config::USE_PARK_FORMALISM | Config::USE_BRUNE_FORMALISM);
    cfg.paramMask |= Config::USE_BRUNE_FORMALISM;
    if (park) cfg.paramMask |= Config::USE_PARK_FORMALISM;
    mkdir(outdir.c_str(), 0755);
    cfg.outputdir = outdir + "/";
    if (!ok) return;
    compound = new CNuc();
    data = new WorkingCopy();
    if (compound->Fill(cfg) == -1 || data->Fill(cfg, compound) != 0) {
      ok = false;
      return;
    }
    RequestCaptureAyTerms(*compound, *data);
    compound->Initialize(cfg);
    ok = data->Initialize(compound, cfg) != -1;
    compound->FillMnParams(params, &cfg);
    data->FillMnParams(params);
    compound->FillCompoundFromParams(params.GetMinuitParams().Params());
    compound->CalcShiftFunctions(cfg);
  }
  ~Model() {
    delete compound;
    delete data;
  }
  bool ok = true;
  int N() const { return (int)params.GetMinuitParams().Params().size(); }
  bool Fixed(int i) const { return params.GetMinuitParams().Parameter(i).IsFixed(); }
  // The sampler's start, as AZURE2.cpp makes it: the varying parameters, amplitudes or observed values.
  std::vector<double> Start(bool useRWA) const {
    vector_r all;
    for (int i = 0; i < N(); i++) all.push_back(params.GetMinuitParams().Value(i));
    if (!useRWA) {
      CNuc *c = compound->Clone();
      c->FillCompoundFromParams(all);
      c->CalcShiftFunctions(cfg);
      c->TransformOut(cfg);
      vector_r phys = c->GetTransformParams(cfg);
      for (int i = (int)phys.size(); i < N(); i++) phys.push_back(all[i]);
      all = phys;
      delete c;
    }
    std::vector<double> x;
    for (int i = 0; i < N(); i++)
      if (!Fixed(i)) x.push_back(all[i]);
    return x;
  }
};

// The fit: AZURECalc::operator() with the limits manager and the prior covariance as AZUREMain sets them.
struct Fit {
  Model &m;
  azure2::Results results;
  std::unique_ptr<AZURECalc> calc;
  explicit Fit(Model &model) : m(model) {
    m.limits.reset(new ParameterLimitsManager(&m.cfg, m.compound, m.data, &m.params));
    m.limits->ReadParameterSettings();
    m.limits->ApplyAllParameterSettings(m.params);
    if (!m.cfg.priorCovariance.empty()) {
      std::vector<std::string> names;
      std::vector<bool> fixed;
      for (int i = 0; i < m.N(); i++) {
        names.push_back(m.params.GetMinuitParams().GetName(i));
        fixed.push_back(m.Fixed(i));
      }
      PriorCovariance::Load(m.cfg, names, fixed, AZUREParams::BasisForConfig(m.cfg), m.priorCov);
      m.data->SetPriorCovariance(m.priorCov);
      if (m.priorCov && m.priorCov->Size() > 0) ReplaceDiagonalPriors(m.compound, m.data, *m.priorCov);
    }
    calc.reset(new AZURECalc(m.data, results, m.compound, m.cfg, m.limits.get()));
    calc->SetErrorDef(1.0);
    std::vector<bool> fixed;
    for (int i = 0; i < m.N(); i++) fixed.push_back(m.Fixed(i));
    calc->SetFixedMask(fixed);
    m.data->SetFit(true);
  }
  // chi^2 of the fit at the sampler's point x (observed values turned into the fit's vector as the sampler does).
  double Chi2At(const AZURECalcMCMC &mc, const std::vector<double> &x, bool useRWA) {
    vector_r full = useRWA ? mc.ReconstructFullParameters(x) : mc.ReconstructFullParametersPhysical(x);
    if (!useRWA) {
      CNuc *c = m.compound->Clone();
      c->FillCompoundFromParamsPhysical(full);
      c->TransformIn(m.cfg);
      AZUREParams p;
      c->FillMnParams(p, &m.cfg);
      m.data->FillMnParams(p);
      vector_r fit = p.GetMinuitParams().Params();
      for (size_t i = (size_t)m.data->Data().GetNormParamOffset(); i < fit.size() && i < full.size(); i++) fit[i] = full[i];
      full = fit;
      delete c;
    }
    return (*calc)(full);
  }
};

// The sampler on a model, set up as RunMCMCSampling sets it up (parameters, automatic priors, prior covariance).
std::unique_ptr<AZURECalcMCMC> Sampler(Model &m, bool useRWA) {
  std::unique_ptr<AZURECalcMCMC> mc(new AZURECalcMCMC(m.data, m.compound, m.cfg));
  mc->PrepareRun(m.Start(useRWA), useRWA);
  return mc;
}

bool Close(double a, double b, double rel) { return std::fabs(a - b) <= rel * std::max(1.0, std::fabs(b)); }

// (a) and the T2 identity: ln P(x) = -chi2(x)/2 + C at three points x.
void Identity(const std::string &label, Model &m, bool useRWA, const std::vector<std::vector<double>> &points) {
  Fit fit(m);
  std::unique_ptr<AZURECalcMCMC> mc = Sampler(m, useRWA);
  const double C = mc->LogPriorNormalisation();
  for (size_t k = 0; k < points.size(); k++) {
    const double lp = mc->LogProbability(points[k], useRWA);
    const double expect = -0.5 * fit.Chi2At(*mc, points[k], useRWA) + C;
    check(label + ": point " + std::to_string(k) + ": ln P = -chi2/2 + C", Close(lp, expect, 1e-10),
          "ln P " + num(lp) + ", -chi2/2 + C " + num(expect) + " (C " + num(C) + ")");
  }
}

// The record project's points: the start, a displaced point, and the start with entry `out` set to `value`.
std::vector<double> Displaced(std::vector<double> x) {
  for (size_t k = 0; k < x.size(); k++) x[k] *= 1.0 + 0.02 * (k % 2 ? 1 : -1);
  return x;
}

std::string WriteVariant(const std::string &source, const std::string &target, bool priors, bool nearZero) {
  std::string text = readFile(source);
  bool ok = !text.empty();
  if (priors) {
    ok = ok && replaceAll(text, "    total_width  0.0385027 eV fit\n", "    total_width  0.0385027 eV fit prior 0.040 0.004\n");
    ok = ok && replaceAll(text, "    remainder    g+15O[5.183]\n", "    remainder    g+15O[5.183] prior 0.15 0.01\n");
    ok = ok && replaceAll(text, "    branching    g+15O[6.1763] 0.5843746023006179 fit\n",
                          "    branching    g+15O[6.1763] 0.5843746023006179 fit prior 0.6 0.01\n");
  }
  if (nearZero) {  // one free number, under a prior that pulls it below its limit
    ok = ok && replaceAll(text, "    total_width  0.0385027 eV fit prior", "    total_width  0.0385027 eV fixed prior");
    ok = ok && replaceAll(text, "    branching    g+15O 0.015030114771171893 fit\n",
                          "    branching    g+15O 0.015030114771171893 fixed\n");
    ok = ok && replaceAll(text, "    branching    g+15O[6.7931] 0.23037345432917689 fit\n",
                          "    branching    g+15O[6.7931] 0.23037345432917689 fixed\n");
    ok = ok && replaceAll(text, "    branching    g+15O[6.1763] 0.5843746023006179 fit prior",
                          "    branching    g+15O[6.1763] 0.5843746023006179 fixed prior");
    ok = ok && replaceAll(text, "    branching    g+15O[6.8594] 0.0014024990455214828 fit\n",
                          "    branching    g+15O[6.8594] 0.0014024990455214828 fit prior -0.01 0.00001\n");
  }
  std::ofstream(target.c_str()) << text;
  return ok ? target : "";
}

void RemoveChain(const std::string &dir) {
  std::remove((dir + "/samples.mcmc").c_str());
  std::remove((dir + "/walkers.mcmc").c_str());
}

// A short chain (8 walkers, `steps` steps, one thread) of model m into its output directory; the samples (the
// initial ensemble and every step).
std::vector<std::vector<double>> Chain(Model &m, bool useRWA, int steps, double spreadPercent = 1.0) {
  RemoveChain(m.cfg.outputdir);
  AZURECalcMCMC mc(m.data, m.compound, m.cfg);
  std::vector<std::vector<double>> samples;
  mc.RunMCMCSampling(8, steps, m.Start(useRWA), samples, spreadPercent, 1, useRWA, 0.5);
  return samples;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: sampler_objective_test <record.azr> <T2.azr>\n");
    return 2;
  }
  // Copies in the working directory: a project's relative paths (data, output) are its file's directory's.
  const std::string record = "record.azr", t2 = "T2.azr";
  std::ofstream(record.c_str()) << readFile(argv[1]);
  std::ofstream(t2.c_str()) << readFile(argv[2]);
  const std::string withPriors = WriteVariant(record, "record_priors.azr", true, false);
  const std::string nearZero = WriteVariant(record, "record_near_zero.azr", true, true);
  check("the record project's variants are made (every edit found)", !withPriors.empty() && !nearZero.empty());
  if (withPriors.empty() || nearZero.empty()) return 1;

  // The free parameters of the record project: Gamma_gamma, then four branchings (the remainder is not a parameter).
  for (int park = 0; park <= 1; park++)
    for (int useRWA = 1; useRWA >= 0; useRWA--) {
      const std::string label = std::string(park ? "Park" : "Brune") + (useRWA ? ", fit basis" : ", observed basis");
      std::printf("== (a) %s\n", label.c_str());
      Model m(withPriors, park, "output");
      if (!m.ok) {
        check(label + ": the project loads", false);
        continue;
      }
      const std::vector<double> x0 = m.Start(useRWA);
      check(label + ": five free parameters (the record's)", x0.size() == 5, std::to_string(x0.size()));
      if (x0.size() != 5) continue;
      Identity(label, m, useRWA, {x0, Displaced(x0)});
      std::unique_ptr<AZURECalcMCMC> mc = Sampler(m, useRWA);
      for (const auto &out : std::vector<std::pair<int, double>>{{1, -0.01}, {2, 1.01}, {0, -1e-4}}) {
        std::vector<double> x = x0;
        x[out.first] = out.second;
        const double lp = mc->LogProbability(x, useRWA);
        check(label + ": outside the limits (entry " + std::to_string(out.first) + " = " + num(out.second) +
                  "): ln P = -infinity",
              std::isinf(lp) && lp < 0, num(lp));
      }
      // (b) the record's priors (T11) reach the sampler: ln P with them minus ln P without = -1/2 x the fit's
      // difference, and not 0.
      Model bare(record, park, "output");
      Fit fitBare(bare), fitPriors(m);
      std::unique_ptr<AZURECalcMCMC> mcBare = Sampler(bare, useRWA);
      for (const std::vector<double> &x : {x0, Displaced(x0)}) {
        const double dSampler = mc->LogProbability(x, useRWA) - mcBare->LogProbability(x, useRWA);
        const double dFit = -0.5 * (fitPriors.Chi2At(*mc, x, useRWA) - fitBare.Chi2At(*mcBare, x, useRWA));
        check(label + ": (b) the record's priors move ln P by -1/2 x their term: " + num(dSampler),
              dSampler != 0.0 && Close(dSampler, dFit, 1e-9), "fit " + num(dFit));
      }
      // The chain stays inside the limits.
      Model z(nearZero, park, "output");
      const std::vector<std::vector<double>> chain = Chain(z, useRWA, 60, 30.0);  // walkers start at b 0.0014 +- 30 %
      int outside = 0;
      double minB = 1.0;
      for (const std::vector<double> &row : chain) {
        if (row.size() != 1 || row[0] < 0.0 || row[0] > 1.0) outside++;
        if (!row.empty()) minB = std::min(minB, row[0]);
      }
      check(label + ": a chain of " + std::to_string(chain.size()) + " rows stays inside the limits (smallest b " +
                num(minB) + ")",
            chain.size() >= 8 * 60 && outside == 0, std::to_string(outside) + " values outside");
    }

  std::printf("== T2 (fit basis): the identity with the Gaussian priors' constant\n");
  {
    Model m(t2, false, "output_t2");
    if (!m.ok) {
      check("T2: the project loads", false);
    } else {
      const std::vector<double> x0 = m.Start(true);
      std::vector<double> x1 = x0;
      for (size_t k = 0; k < x1.size(); k++) x1[k] = x1[k] * 1.01 + (x1[k] == 0.0 ? 1e-4 : 0.0);
      Identity("T2", m, true, {x0, x1});
    }
  }

  std::printf("== (c) two samplers in two threads give the chains of two runs one after the other\n");
  {
    auto chainFile = [](const std::string &azr, bool useRWA, int steps, const std::string &dir) {
      Model m(azr, false, dir);
      if (!m.ok) return std::string();
      Chain(m, useRWA, steps);
      const std::string chain = readFile(dir + "/samples.mcmc");
      if (chain.empty()) std::printf("  (no chain in %s; the run's log:)\n%s\n", dir.c_str(), m.log.str().c_str());
      return chain;
    };
    const std::string soloA = chainFile(withPriors, false, 60, "solo_a"), soloB = chainFile(record, true, 60, "solo_b");
    std::string pairA, pairB;
    std::thread ta([&] { pairA = chainFile(withPriors, false, 60, "pair_a"); });
    std::thread tb([&] { pairB = chainFile(record, true, 60, "pair_b"); });
    ta.join();
    tb.join();
    check("the two chains differ from each other (two projects)", !soloA.empty() && !soloB.empty() && soloA != soloB);
    check("the record project's chain in a thread beside another sampler = its chain alone (" +
              std::to_string(soloA.size()) + " bytes)",
          !soloA.empty() && pairA == soloA);
    check("the project without the priors: its chain in a thread beside another sampler = its chain alone (" +
              std::to_string(soloB.size()) + " bytes)",
          !soloB.empty() && pairB == soloB);
  }

  std::printf("%s\n", failures ? "FAILED" : "all passed");
  return failures ? 1 : 0;
}
