#include "ECTable.h"
#include "ECIntegral.h"
#include "CoulFunc.h"
#include "Config.h"
#include "PPair.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <future>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <vector>

namespace {

/// The tabulated (reduced) integral at one node; ok = false when the node could
/// not be computed (non-finite value or a GSL error), which sends every
/// energy that needs it to the direct integral.
struct NodeValue {
  complex value;
  bool ok;
};

struct Term {
  std::mutex m;
  std::map<long, std::shared_future<NodeValue>> nodes;
  // AZURE2_EC_TABLE_CHECK: the worst relative difference and where, the
  // largest direct value and the largest absolute difference.
  double maxDev = -1.0, eAtMax = 0.0, absAtMax = 0.0, maxAbs = 0.0, maxAbsDiff = 0.0;
  double eMin = 1e300, eMax = -1e300;
};

std::mutex g_termsMutex;
std::map<std::vector<double>, std::unique_ptr<Term>> g_terms;

std::atomic<long long> g_nodes(0), g_lookups(0), g_direct(0);
std::atomic<long long> g_nodesReport(0);
std::mutex g_checkMutex;
double g_maxDeviation = -1.0;

const int kStencil = 8;  // Lagrange points; the stencil is nodes i0-3 .. i0+4

bool CheckMode() {
  static const bool on = [] {
    const char *v = std::getenv("AZURE2_EC_TABLE_CHECK");
    return v && *v && std::string(v) != "0";
  }();
  return on;
}

int NodesPerDecade(const Config &configure) {
  static const int fromEnv = [] {
    const char *v = std::getenv("AZURE2_EC_TABLE_NODES");
    return v ? std::atoi(v) : 0;
  }();
  int n = fromEnv > 0 ? fromEnv : configure.ecTableNodesPerDecade;
  return n > 0 ? n : 200;
}

bool RateMode(const Config &configure) {
  return (configure.paramMask & Config::CALCULATE_REACTION_RATE) != 0;
}

/// Everything the integral depends on besides the energy.
std::vector<double> TermKey(PPair *pair, const Config &configure, const CoulFunc &coul, int li, int lf,
                            double si, double sf, double ji, double jf, int multL, char radType,
                            double levelEnergy, bool isChannelCapture, int nodesPerDecade) {
  return {(double)pair->GetZ(1), (double)pair->GetZ(2), pair->GetM(1), pair->GetM(2), pair->GetRedMass(),
          pair->GetJ(1), pair->GetJ(2), pair->GetG(1), pair->GetG(2), pair->GetChRad(),
          pair->GetSepE() + pair->GetExE(),
          (double)li, (double)lf, si, sf, ji, jf, (double)multL, (double)radType, levelEnergy,
          isChannelCapture ? 1.0 : 0.0,
          (configure.paramMask & Config::USE_LONGWAVELENGTH_APPROX) ? 1.0 : 0.0,
          (configure.paramMask & Config::USE_GSL_COULOMB_FUNC) ? 1.0 : 0.0,
          RateMode(configure) ? 1.0 : 0.0,  // the rate integrates on a fixed radial grid
          (configure.ecConvectionCurrent && !(configure.paramMask & Config::USE_LONGWAVELENGTH_APPROX)) ? 1.0 : 0.0,
          configure.useCoulombTable ? 1.0 : 0.0,
          (double)coul.hybridTag(),  // the nuclear potential of a hybrid pair
          (double)nodesPerDecade};
}

/*!
 * The factor divided out of the integral before it is tabulated: sqrt(P_li) of
 * the final pair for an ordinary term (the integral falls like the Coulomb
 * penetration), its hard-sphere phase factor for channel capture (the integral
 * carries it, and it turns with energy).  Both are exact and cheap (the
 * Coulomb table answers them at the channel radius).
 */
bool Reduction(PPair *pair, CoulFunc &coul, int li, double outEnergy, bool isChannelCapture, complex &n) {
  if (isChannelCapture)
    n = coul.HardSpherePhase(li, pair->GetChRad(), outEnergy);
  else
    n = complex(sqrt(coul.Penetrability(li, pair->GetChRad(), outEnergy)), 0.0);
  return std::isfinite(n.real()) && std::isfinite(n.imag()) && std::abs(n) > 0.0;
}

/// Relative tolerance of the node integrals: the direct integrals use 1e-4,
/// which is the scatter (relative to a term's largest value) a table built from
/// them would inherit.  Nodes are few, so they are converged instead.
const double kNodeTolerance = 1e-10;

complex Direct(PPair *pair, const Config &configure, int li, int lf, double si, double sf, double ji,
               double jf, int multL, char radType, double inEnergy, double levelEnergy, bool isChannelCapture) {
  ECIntegral integral(pair, configure);
  return integral(li, lf, si, sf, ji, jf, multL, radType, inEnergy, levelEnergy, isChannelCapture);
}

/// The integral converged to kNodeTolerance, or at the direct integrals' own
/// tolerance where the quadrature cannot reach it (a GSL error or status).
complex Converged(PPair *pair, const Config &configure, int li, int lf, double si, double sf, double ji,
                  double jf, int multL, char radType, double inEnergy, double levelEnergy, bool isChannelCapture) {
  try {
    ECIntegral integral(pair, configure);
    integral.SetRelativeTolerance(kNodeTolerance);
    complex v = integral(li, lf, si, sf, ji, jf, multL, radType, inEnergy, levelEnergy, isChannelCapture);
    if (integral.Status() == 0) return v;
  } catch (...) {
  }
  return Direct(pair, configure, li, lf, si, sf, ji, jf, multL, radType, inEnergy, levelEnergy, isChannelCapture);
}

}  // namespace

complex ECTable::Integral(PPair *finalPair, const Config &configure,
                          int li, int lf, double si, double sf, double ji, double jf,
                          int multL, char radType, double inEnergy, double levelEnergy,
                          bool isChannelCapture) {
  const double threshold = finalPair->GetSepE() + finalPair->GetExE();
  const double outEnergy = inEnergy - threshold;
  CoulFunc coul(finalPair, !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC));
  if (RateMode(configure)) coul.SetUseGlobalCache(false);  // rate energies never recur
  complex reduction;
  if (!(outEnergy > 0.0) || !Reduction(finalPair, coul, li, outEnergy, isChannelCapture, reduction)) {
    ++g_direct;
    return Direct(finalPair, configure, li, lf, si, sf, ji, jf, multL, radType, inEnergy, levelEnergy,
                  isChannelCapture);
  }

  const int perDecade = NodesPerDecade(configure);
  const double h = log(10.0) / perDecade;
  Term *term;
  {
    std::vector<double> key = TermKey(finalPair, configure, coul, li, lf, si, sf, ji, jf, multL, radType,
                                      levelEnergy, isChannelCapture, perDecade);
    std::lock_guard<std::mutex> lock(g_termsMutex);
    std::unique_ptr<Term> &slot = g_terms[key];
    if (!slot) slot.reset(new Term);
    term = slot.get();
  }

  // The stencil: nodes first-3 .. first+4 around u = ln(E_out), t its position in it.
  const double u = log(outEnergy) / h;
  const long first = (long)floor(u) - (kStencil / 2 - 1);
  const double t = u - (double)first;

  // Claim the missing nodes under the lock, compute them outside it, and wait
  // for the ones another thread is computing.
  std::shared_future<NodeValue> futures[kStencil];
  std::unique_ptr<std::promise<NodeValue>> promises[kStencil];  // the nodes this call computes
  bool mine[kStencil];
  {
    std::lock_guard<std::mutex> lock(term->m);
    for (int s = 0; s < kStencil; ++s) {
      auto it = term->nodes.find(first + s);
      mine[s] = (it == term->nodes.end());
      if (mine[s]) {
        promises[s].reset(new std::promise<NodeValue>);
        futures[s] = promises[s]->get_future().share();
        term->nodes.emplace(first + s, futures[s]);
      } else
        futures[s] = it->second;
    }
  }
  for (int s = 0; s < kStencil; ++s) {
    if (!mine[s]) continue;
    NodeValue node{complex(0.0, 0.0), false};
    const double nodeOut = exp((double)(first + s) * h);
    try {
      complex n;
      if (Reduction(finalPair, coul, li, nodeOut, isChannelCapture, n)) {
        complex v = Converged(finalPair, configure, li, lf, si, sf, ji, jf, multL, radType, nodeOut + threshold,
                              levelEnergy, isChannelCapture) / n;
        node.ok = std::isfinite(v.real()) && std::isfinite(v.imag());
        node.value = v;
      }
    } catch (...) {
      node.ok = false;
    }
    promises[s]->set_value(node);
    ++g_nodes;
  }

  complex sum(0.0, 0.0);
  bool ok = true;
  for (int s = 0; s < kStencil; ++s) {
    const NodeValue &node = futures[s].get();
    if (!node.ok) {
      ok = false;
      break;
    }
    double w = 1.0;
    for (int q = 0; q < kStencil; ++q)
      if (q != s) w *= (t - (double)q) / (double)(s - q);
    sum += w * node.value;
  }
  if (!ok) {
    ++g_direct;
    return Direct(finalPair, configure, li, lf, si, sf, ji, jf, multL, radType, inEnergy, levelEnergy,
                  isChannelCapture);
  }
  ++g_lookups;
  complex result = sum * reduction;

  if (CheckMode()) {
    complex exact = Converged(finalPair, configure, li, lf, si, sf, ji, jf, multL, radType, inEnergy, levelEnergy,
                              isChannelCapture);
    if (std::abs(exact) > 0.0) {
      double dev = std::abs(result - exact) / std::abs(exact);
      {
        std::lock_guard<std::mutex> lock(term->m);
        if (dev > term->maxDev) {
          term->maxDev = dev;
          term->eAtMax = outEnergy;
          term->absAtMax = std::abs(exact);
        }
        term->maxAbs = std::max(term->maxAbs, std::abs(exact));
        term->maxAbsDiff = std::max(term->maxAbsDiff, std::abs(result - exact));
        term->eMin = std::min(term->eMin, outEnergy);
        term->eMax = std::max(term->eMax, outEnergy);
      }
      std::lock_guard<std::mutex> lock(g_checkMutex);
      if (dev > g_maxDeviation) g_maxDeviation = dev;
    }
  }
  return result;
}

long long ECTable::NodesComputed() { return g_nodes; }
long long ECTable::Lookups() { return g_lookups; }
long long ECTable::DirectCalls() { return g_direct; }

double ECTable::MaxCheckDeviation() {
  std::lock_guard<std::mutex> lock(g_checkMutex);
  return g_maxDeviation;
}

std::string ECTable::Report() {
  size_t terms;
  {
    std::lock_guard<std::mutex> lock(g_termsMutex);
    terms = g_terms.size();
  }
  std::ostringstream s;
  s << "External capture table: " << (g_lookups + g_direct) << " integrals, " << (g_nodes - g_nodesReport)
    << " computed at table nodes (" << terms << " terms)";
  if (g_direct > 0) s << ", " << g_direct << " directly";
  if (CheckMode()) {
    {
      std::lock_guard<std::mutex> lock(g_checkMutex);
      s << "; check: largest relative difference from the direct integrals " << g_maxDeviation;
    }
    // One line per term: li lf L type, E_out range, the worst point, and the
    // largest difference relative to the term's largest value.
    std::lock_guard<std::mutex> lock(g_termsMutex);
    for (auto &kv : g_terms) {
      const std::vector<double> &k = kv.first;
      Term &t = *kv.second;
      std::lock_guard<std::mutex> tl(t.m);
      if (t.maxDev < 0.0) continue;
      s << "\n\t  term li=" << k[11] << " lf=" << k[12] << " " << (char)k[18] << k[17]
        << (k[20] != 0.0 ? " chancap" : "") << " Ef=" << k[19] << ": E_out " << t.eMin << ".." << t.eMax
        << ", worst " << t.maxDev << " at " << t.eAtMax << " (|I| " << t.absAtMax << ", max |I| " << t.maxAbs
        << "), max|diff|/max|I| " << (t.maxAbs > 0.0 ? t.maxAbsDiff / t.maxAbs : 0.0);
    }
  }
  return s.str();
}

void ECTable::ResetCounters() {
  g_lookups = 0;
  g_direct = 0;
  g_nodesReport = g_nodes.load();
  std::lock_guard<std::mutex> lock(g_checkMutex);
  g_maxDeviation = -1.0;
}

void ECTable::Clear() {
  std::lock_guard<std::mutex> lock(g_termsMutex);
  g_terms.clear();
}
