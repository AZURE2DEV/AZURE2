#include "ECIntegral.h"
#include "AngCoeff.h"
#include "EffectiveCharge.h"
#include "Config.h"
#include <math.h>
#include <gsl/gsl_integration.h>
#include <assert.h>
#include <cmath>
#include <map>
#include <vector>

extern double DoubleFactorial(int);

namespace {
/*!
 * Per-x memo of the two expensive quantities shared by the FW and GW
 * integrands: the Coulomb functions and the scaled bound-state Whittaker.
 *
 * FWIntegrand and GWIntegrand are identical apart from a final coul.F vs coul.G,
 * and GSL evaluates them at bit-identical abscissae (same Gauss-Kronrod nodes,
 * same variable transform).  So the first integrand to reach a given x computes
 * coul + whit once and stores them; the second reads them back.  Nothing about
 * the integration or the arithmetic changes -- the result is byte-identical to
 * two independent passes -- only the ~2x redundant special-function work is
 * removed.  The map is thread-local and cleared at the start of each FW/GW
 * Integrate().
 */
struct FGMemoEntry {
  double F;
  double G;
  double whit;
};
thread_local std::map<double, FGMemoEntry> g_fgMemo;
}  // namespace

void ECIntegral::ClearFGMemo() { g_fgMemo.clear(); }

/*!
 * Computes (or reuses) the shared FW/GW quantities at x: the Coulomb functions
 * and the scaled bound-state Whittaker.  File-scope helper; ECIntegral::Params
 * is public so it can be named here.
 */
static FGMemoEntry fgEvaluate(double x, ECIntegral::Params *params) {
  auto it = g_fgMemo.find(x);
  if (it != g_fgMemo.end()) return it->second;
  struct CoulWaves coul = params->coulFunc->operator()(params->liValue, x, params->pairEnergy);
  double sOut;
  double mOut = params->whitFunc->Scaled(params->lfValue, x, params->bindingEnergy, sOut);
  double exponent = sOut - params->log10RefOut;
  double whit = (exponent < -280.0) ? 0.0 : mOut * pow(10.0, exponent);
  FGMemoEntry e = {coul.F, coul.G, whit};
  g_fgMemo[x] = e;
  return e;
}

double ECIntegral::FWIntegrand(double x, void *p) {
  Params *params = (Params *)p;
  FGMemoEntry e = fgEvaluate(x, params);
  double returnValue = e.F * e.whit * pow(x, params->multLValue);
  return (!params->useLongWavelengthApprox) ? params->effectiveCharge->operator()(x) * returnValue : returnValue;
}

double ECIntegral::GWIntegrand(double x, void *p) {
  Params *params = (Params *)p;
  FGMemoEntry e = fgEvaluate(x, params);
  double returnValue = e.G * e.whit * pow(x, params->multLValue);
  return (!params->useLongWavelengthApprox) ? params->effectiveCharge->operator()(x) * returnValue : returnValue;
}

double ECIntegral::WWIntegrand(double x, void *p) {
  Params *params = (Params *)p;
  WhitFunc *theWhitFunc = (params->whitFunc.get());
  int liValue = (params->liValue);
  int lfValue = (params->lfValue);
  int multLValue = (params->multLValue);
  double pairEnergy = (params->pairEnergy);
  double bindingEnergy = (params->bindingEnergy);

  double returnValue;
  if (params->useScaledWhittaker) {
    // Both Whittaker factors can be ~1e-500 here.  Evaluate them with their
    // powers of ten kept separate and subtract the reference decades measured
    // at the channel radius, so the integrand is O(1) where it matters.  The
    // common scale is put back by the caller.
    double sIn, sOut;
    double mIn = theWhitFunc->Scaled(liValue, x, fabs(pairEnergy), sIn);
    double mOut = theWhitFunc->Scaled(lfValue, x, fabs(bindingEnergy), sOut);
    double exponent = sIn + sOut - params->log10RefIn - params->log10RefOut;
    if (exponent < -280.0) return 0.0;  // genuinely negligible
    returnValue = mIn * mOut * pow(10.0, exponent) * pow(x, multLValue);
  } else {
    double whitIn = theWhitFunc->operator()(liValue, x, fabs(pairEnergy));
    double whitOut = theWhitFunc->operator()(lfValue, x, fabs(bindingEnergy));
    returnValue = whitIn * whitOut * pow(x, multLValue);
  }
  return (!params->useLongWavelengthApprox) ? params->effectiveCharge->operator()(x) * returnValue : returnValue;
}

namespace {
/*!
 * Per-thread GSL integration workspace, reused across calls.  Allocating and
 * freeing a 1000-element workspace on every external-capture integral showed
 * up as ~8% of total runtime in a reaction-rate profile, and the workspace
 * carries no state between calls.
 */
struct ECWorkspace {
  gsl_integration_workspace *w;
  ECWorkspace() :
    w(gsl_integration_workspace_alloc(1000)) {}
  ~ECWorkspace() {
    if (w) gsl_integration_workspace_free(w);
  }
};
gsl_integration_workspace *ec_workspace() {
  static thread_local ECWorkspace holder;
  return holder.w;
}

/* ------------------------------------------------------------------------
 * Fixed radial grid + cached final-state Whittaker, used ONLY for the reaction
 * rate (FW/GW branch).
 *
 * In the rate the radial integral is evaluated at thousands of incident
 * energies.  The final-state Whittaker whit(lf, r, bindingEnergy) does NOT
 * depend on the incident energy (bindingEnergy is the bound level's energy,
 * fixed), yet it is the dominant cost (gsl_sf_hyperg_U, ~85% of runtime) and
 * the adaptive qagiu recomputes it from scratch every time because it samples
 * different radii at each energy.
 *
 * Replacing qagiu with a FIXED Gauss-Laguerre radial grid makes the sample
 * radii energy-independent, so the final Whittaker on that grid is computed
 * once and reused for every energy.  Only the entrance Coulomb function (which
 * does depend on energy, ~10% of cost) is recomputed per energy.
 *
 * Gauss-Laguerre integrates weight (x-a)^0 exp(-b(x-a)) exactly for a
 * polynomial * that weight; the integrand here is Coulomb * Whittaker * r^L,
 * which decays like exp(-b r) with b set to the Whittaker decay rate, so the
 * un-weighted part is smooth and the rule converges quickly.
 * ------------------------------------------------------------------------ */
struct RateGridKey {
  int z1, z2;
  double redmass;
  int lf;
  double bindingEnergy;
  double chanrad;
  bool operator<(const RateGridKey &o) const {
    if (z1 != o.z1) return z1 < o.z1;
    if (z2 != o.z2) return z2 < o.z2;
    if (redmass != o.redmass) return redmass < o.redmass;
    if (lf != o.lf) return lf < o.lf;
    if (bindingEnergy != o.bindingEnergy) return bindingEnergy < o.bindingEnergy;
    return chanrad < o.chanrad;
  }
};
struct RateGridEntry {
  std::vector<double> x;        // radial nodes
  std::vector<double> W;        // effective weights (Laguerre weight * e^{+b(x-a)})
  std::vector<double> whitOut;  // scaled final Whittaker at each node (ref subtracted)
  double refOut;                // log10 reference decades of the final Whittaker
  bool ok;
};
const int kRateRadialPoints = 150;

const RateGridEntry &getRateGrid(WhitFunc *wf, int lf, double bindingEnergy, double chanrad) {
  static thread_local std::map<RateGridKey, RateGridEntry> cache;
  RateGridKey key{wf->z1(), wf->z2(), wf->redmass(), lf, bindingEnergy, chanrad};
  auto it = cache.find(key);
  if (it != cache.end()) return it->second;

  RateGridEntry e;
  e.ok = false;
  e.refOut = 0.0;
  double b = sqrt(2.0 * uconv * wf->redmass() * bindingEnergy) / hbarc;  // Whittaker decay rate
  if (b > 1.0e-6 && std::isfinite(b)) {
    gsl_integration_fixed_workspace *ws = gsl_integration_fixed_alloc(
        gsl_integration_fixed_laguerre, kRateRadialPoints, chanrad, b, 0.0, 0.0);
    if (ws) {
      const double *xk = gsl_integration_fixed_nodes(ws);
      const double *wk = gsl_integration_fixed_weights(ws);
      double refOut;
      wf->Scaled(lf, chanrad, bindingEnergy, refOut);
      e.refOut = refOut;
      e.x.resize(kRateRadialPoints);
      e.W.resize(kRateRadialPoints);
      e.whitOut.resize(kRateRadialPoints);
      for (int k = 0; k < kRateRadialPoints; ++k) {
        e.x[k] = xk[k];
        e.W[k] = wk[k] * exp(b * (xk[k] - chanrad));  // undo the exp(-b(x-a)) weight
        double sk;
        double mk = wf->Scaled(lf, xk[k], bindingEnergy, sk);
        double expo = sk - refOut;
        e.whitOut[k] = (expo < -280.0) ? 0.0 : mk * pow(10.0, expo);
      }
      e.ok = true;
      gsl_integration_fixed_free(ws);
    }
  }
  auto res = cache.emplace(std::move(key), std::move(e));
  return res.first->second;
}
}  // namespace

void ECIntegral::Integrate(double chanrad) {
  gsl_integration_workspace *w = ec_workspace();

  if (params_.pairEnergy < 0.0) {
    gsl_function WW;
    WW.function = &WWIntegrand;
    WW.params = &params_;

    double wwintresult, wwinterror;

    gsl_integration_qagiu(&WW, chanrad, 0.0, 1e-4, 1000, w, &wwintresult, &wwinterror);
    GW_ = wwintresult;
    FW_ = 0.0;
  } else {
    // Reaction rate: use the fixed Gauss-Laguerre grid with the cached
    // energy-independent final Whittaker.  Everything else keeps the adaptive
    // qagiu, which is the trusted path for fits and extrapolations.
    bool rateMode = configure_ && (configure_->paramMask & Config::CALCULATE_REACTION_RATE);
    const RateGridEntry *grid = rateMode ? &getRateGrid(params_.whitFunc.get(), params_.lfValue, params_.bindingEnergy, chanrad)
                                         : nullptr;
    if (grid && grid->ok) {
      params_.log10RefOut = grid->refOut;  // CalculateDirect restores 10^refOut
      long double fw = 0.0L, gw = 0.0L;
      for (size_t k = 0; k < grid->x.size(); ++k) {
        double x = grid->x[k];
        struct CoulWaves coul = params_.coulFunc->operator()(params_.liValue, x, params_.pairEnergy);
        double eff = params_.useLongWavelengthApprox ? 1.0
                                                     : params_.effectiveCharge->operator()(x);
        double common = grid->W[k] * eff * grid->whitOut[k] * pow(x, params_.multLValue);
        fw += (long double)coul.F * common;
        gw += (long double)coul.G * common;
      }
      FW_ = (double)fw;
      GW_ = (double)gw;
    } else {
      gsl_function FW;
      FW.function = &FWIntegrand;
      FW.params = &params_;
      gsl_function GW;
      GW.function = &GWIntegrand;
      GW.params = &params_;
      double fwintresult, fwinterror, gwintresult, gwinterror;
      // FW and GW share their expensive per-x evaluations through g_fgMemo;
      // clear it so this call cannot see stale entries from a previous one.
      ClearFGMemo();
      gsl_integration_qagiu(&FW, chanrad, 0.0, 1e-4, 1000, w, &fwintresult, &fwinterror);
      gsl_integration_qagiu(&GW, chanrad, 0.0, 1e-4, 1000, w, &gwintresult, &gwinterror);
      FW_ = fwintresult;
      GW_ = gwintresult;
    }
  }
}

/*!
 * Overloaded operator to make the class instance callable as a function.  The intial and final orbital
 * angular momentum, the gamma multipolarity, and the incoming energy and final state energy in the
 * compound system are the dependent variables. The function returns an ECIntResult structure.
 */

complex ECIntegral::operator()(int theInitialLValue, int theFinalLValue,
                               double theInitialSValue, double theFinalSValue,
                               double theInitialJValue, double theFinalJValue,
                               int theLMult, char radType,
                               double inEnergy, double levelEnergy,
                               bool isChannelCapture, bool forceRecalc) {
  complex integral = CalculateDirect(theInitialLValue, theFinalLValue, theInitialSValue, theFinalSValue,
                                     theInitialJValue, theFinalJValue, theLMult, radType,
                                     inEnergy, levelEnergy, isChannelCapture);
  return integral;
}

complex ECIntegral::CalculateDirect(int theInitialLValue, int theFinalLValue,
                                    double theInitialSValue, double theFinalSValue,
                                    double theInitialJValue, double theFinalJValue,
                                    int theLMult, char radType,
                                    double inEnergy, double levelEnergy,
                                    bool isChannelCapture) {
  ResetIntegrals();

  double sepEnergy = pair()->GetSepE() + pair()->GetExE();
  double outEnergy = inEnergy - sepEnergy;
  double chanRad = pair()->GetChRad();
  double redMass = pair()->GetRedMass();

  EffectiveCharge effectiveChargeFunc(pair(), inEnergy - levelEnergy, theLMult);

  params_.effectiveCharge = &effectiveChargeFunc;
  params_.liValue = theInitialLValue;
  params_.lfValue = theFinalLValue;
  params_.multLValue = theLMult;
  params_.pairEnergy = outEnergy;
  params_.bindingEnergy = fabs(levelEnergy - sepEnergy);

  // Whittaker functions are evaluated in scaled form throughout: measure their
  // decades at the channel radius, integrate the rescaled integrand, and put
  // the scale back afterwards.  Both branches need this -- a sub-threshold
  // entrance channel makes the WW integrand ~1e-500 (which used to divide as
  // 0/0), and a final level close to the separation energy does the same to
  // the bound-state factor in FW/GW (which used to vanish silently).
  {
    double sOut;
    whitfunction()->Scaled(theFinalLValue, chanRad, params_.bindingEnergy, sOut);
    params_.log10RefOut = sOut;
    params_.log10RefIn = 0.0;
    params_.useScaledWhittaker = false;
    if (outEnergy <= 0.0) {
      double sIn;
      whitfunction()->Scaled(theInitialLValue, chanRad, fabs(outEnergy), sIn);
      params_.log10RefIn = sIn;
      params_.useScaledWhittaker = true;
    }
  }

  if (radType == 'E')
    params_.multLValue = theLMult;
  else
    params_.multLValue = 0;
  Integrate(chanRad);

  complex overlapIntegral(0., 0.);
  if (outEnergy > 0.0) {
    struct CoulWaves
        coul = coulfunction()->operator()(theInitialLValue, chanRad, outEnergy);
    if (isChannelCapture) {
      complex chanExpHSP(coul.G / sqrt(pow(coul.F, 2.0) + pow(coul.G, 2.0)),
                         -coul.F / sqrt(pow(coul.F, 2.0) + pow(coul.G, 2.0)));
      overlapIntegral = complex(0.0, -0.5) *
          sqrt(coulfunction()->Penetrability(theInitialLValue, chanRad, outEnergy)) *
          pow(redMass * uconv / 2. / fabs(outEnergy), 0.25) / sqrt(hbarc) *
          chanExpHSP * (GW() + complex(0.0, 1.0) * FW());
    } else
      overlapIntegral = (coul.G / sqrt(pow(coul.F, 2.0) + pow(coul.G, 2.0)) * FW() - coul.F / sqrt(pow(coul.F, 2.0) + pow(coul.G, 2.0)) * GW()) *
          pow(redMass * uconv / 2. / fabs(outEnergy), 0.25) / sqrt(hbarc);
    // FW_/GW_ were integrated with the bound-state Whittaker rescaled by
    // 10^-log10RefOut; restore it.  When the reference is far enough below the
    // double range the contribution really is negligible, so it becomes zero
    // rather than a denormal or an infinity.
    if (params_.log10RefOut != 0.0) {
      if (params_.log10RefOut < -290.0)
        overlapIntegral = complex(0.0, 0.0);
      else if (params_.log10RefOut > 290.0)
        overlapIntegral *= pow(10.0, 290.0);
      else
        overlapIntegral *= pow(10.0, params_.log10RefOut);
    }
  } else {
    assert(isChannelCapture);
    // GW() now holds the integral divided by 10^(log10RefIn+log10RefOut), and
    // the Whittaker in the denominator is mWhit*10^log10RefIn.  The RefIn
    // decades therefore cancel exactly -- which is the whole point, since both
    // quantities are ~1e-500 and each underflows to zero on its own.  What
    // survives is the bound-state scale, applied in log space so that a
    // genuinely negligible result becomes zero instead of a NaN.
    double sWhit;
    double mWhit = whitfunction()->Scaled(theInitialLValue, chanRad, fabs(outEnergy), sWhit);
    double ratio = 0.0;
    if (mWhit != 0.0) {
      double log10Ratio = log10(fabs(GW() / mWhit)) + params_.log10RefOut;
      if (log10Ratio > -290.0) {
        double sign = ((GW() / mWhit) < 0.0) ? -1.0 : 1.0;
        ratio = sign * pow(10.0, log10Ratio);
      }
    }
    overlapIntegral = complex(0.0, -0.5) * ratio *
        sqrt(redMass * uconv * chanRad) / hbarc;
  }


  double effectiveCharge;
  if (radType == 'E') {
    if (params_.useLongWavelengthApprox) {
      double totalM = pair()->GetM(1) + pair()->GetM(2);
      effectiveCharge = sqrt(fstruc * hbarc) * (pair()->GetZ(1) * pow(pair()->GetM(2) / totalM, theLMult) + pair()->GetZ(2) * pow(-pair()->GetM(1) / totalM, theLMult));
    } else
      effectiveCharge = 1.;
  } else {
    effectiveCharge = redMass * 1.00727638 *
        (pair()->GetZ(1) / pow(pair()->GetM(1), 2.) +
         pair()->GetZ(2) / pow(pair()->GetM(2), 2.));
  }

  complex ecAmplitude(0.0, 0.0);
  if (radType == 'E') {
    ecAmplitude = complex(0.0, -1.0) *
        effectiveCharge * sqrt((8. * (2. * theLMult + 1.) * (theLMult + 1.)) / theLMult) / DoubleFactorial(2 * theLMult + 1) *
        pow(complex(0., 1.0), theInitialLValue + theLMult - theFinalLValue) *
        AngCoeff::ClebGord(theInitialLValue, theLMult, theFinalLValue, 0, 0, 0) * sqrt(2. * theInitialLValue + 1.) * sqrt(2. * theFinalJValue + 1.) *
        AngCoeff::Racah(theLMult, theFinalLValue, theInitialJValue, theInitialSValue, theInitialLValue, theFinalJValue);
  } else {
    complex orbitalTerm = effectiveCharge *
        sqrt((2. * theInitialLValue + 1.) * (theInitialLValue + 1.) * theInitialLValue) *
        AngCoeff::Racah(1., theInitialLValue, theInitialJValue, theInitialSValue, theInitialLValue, theFinalJValue);
    complex tau = pow(std::complex<double>(-1., 0.), pair()->GetJ(1) + pair()->GetJ(2)) *
        (pow(complex(-1., 0.), theFinalSValue) *
             sqrt(pair()->GetJ(1) * (pair()->GetJ(1) + 1.) * (2. * pair()->GetJ(1) + 1.)) *
             AngCoeff::Racah(theFinalSValue, pair()->GetJ(1), theInitialSValue, pair()->GetJ(1), pair()->GetJ(2), 1.) *
             pair()->GetG(1) +
         pow(complex(-1., 0.), theInitialSValue) *
             sqrt(pair()->GetJ(2) * (pair()->GetJ(2) + 1.) * (2. * pair()->GetJ(2) + 1.)) *
             AngCoeff::Racah(theFinalSValue, pair()->GetJ(2), theInitialSValue, pair()->GetJ(2), pair()->GetJ(1), 1.) *
             pair()->GetG(2));
    complex spinTerm = -sqrt((2. * theInitialSValue + 1.) * (2. * theFinalSValue + 1.)) *
        AngCoeff::Racah(1, theInitialSValue, theFinalJValue, theInitialLValue, theFinalSValue, theInitialJValue) * tau;
    ecAmplitude = complex(0.0, 1.0) *
        sqrt(fstruc) * pow(hbarc, 1.5) / (2 * 1.00727638 * uconv) * sqrt(16 / 3) * sqrt(2 * theFinalJValue + 1.) *
        (orbitalTerm + spinTerm);
  }

  return ecAmplitude * overlapIntegral;
}
