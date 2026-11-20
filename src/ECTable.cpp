#include "RunScratch.h"
#include "ECTable.h"
#include "azure2/support/Environment.h"
#include "ChannelTable.h"
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

using NodeValue = ECTermNodes::Value;
using Term = ECTermNodes;

const int kStencil = 8;  // Lagrange points; the stencil is nodes i0-3 .. i0+4

bool CheckMode() { return azure2::support::ProcessEnvironment().On(azure2::support::Switch::kEcTableCheck); }

int NodesPerDecade(const Config &configure) {
  const int fromEnv = configure.environment.ecTableNodes;  // AZURE2_EC_TABLE_NODES
  int n = fromEnv > 0 ? fromEnv : configure.ecTableNodesPerDecade;
  return n > 0 ? n : 200;
}

bool RateMode(const Config &configure) {
  return (configure.paramMask & Config::CALCULATE_REACTION_RATE) != 0;
}

/// Everything the integral depends on besides the energy.
std::vector<double> TermKey(const PPair *pair, double radius, const Config &configure, const CoulFunc &coul, int li, int lf,
                            double si, double sf, double ji, double jf, int multL, char radType,
                            double levelEnergy, bool isChannelCapture, int nodesPerDecade) {
  std::vector<double> key{(double)pair->GetZ(1), (double)pair->GetZ(2), pair->GetM(1), pair->GetM(2), pair->GetRedMass(),
          pair->GetJ(1), pair->GetJ(2), pair->GetG(1), pair->GetG(2), radius,
          pair->GetSepE() + pair->GetExE(),
          (double)li, (double)lf, si, sf, ji, jf, (double)multL, (double)radType, levelEnergy,
          isChannelCapture ? 1.0 : 0.0,
          (configure.paramMask & Config::USE_LONGWAVELENGTH_APPROX) ? 1.0 : 0.0,
          RateMode(configure) ? 1.0 : 0.0,  // the rate integrates on a fixed radial grid
          (configure.ecConvectionCurrent && !(configure.paramMask & Config::USE_LONGWAVELENGTH_APPROX)) ? 1.0 : 0.0,
          configure.useCoulombTable ? 1.0 : 0.0,
          (double)nodesPerDecade,
          // relativistic kinematics change eta, rho and the flux factor of every node (checklist t3-25)
          configure.relativistic ? 1.0 : 0.0};
  // the nuclear potential of a hybrid pair, by value
  const std::vector<double> potential = azure2::physics::PotentialKey(coul.HybridPotential());
  key.insert(key.end(), potential.begin(), potential.end());
  return key;
}

/*!
 * The factor divided out of the integral before it is tabulated: sqrt(P_li) of
 * the final pair for an ordinary term (the integral falls like the Coulomb
 * penetration), its hard-sphere phase factor for channel capture (the integral
 * carries it, and it turns with energy).  Both are exact and cheap (the
 * Coulomb table answers them at the channel radius).
 */
bool Reduction(double radius, CoulFunc &coul, int li, int twoJi, double outEnergy, bool isChannelCapture, complex &n) {
  if (isChannelCapture)
    n = coul.HardSpherePhase(li, twoJi, radius, outEnergy);
  else
    n = complex(sqrt(coul.Penetrability(li, twoJi, radius, outEnergy)), 0.0);
  return std::isfinite(n.real()) && std::isfinite(n.imag()) && std::abs(n) > 0.0;
}

/// Relative tolerance of the node integrals: the direct integrals use 1e-4,
/// which is the scatter (relative to a term's largest value) a table built from
/// them would inherit.  Nodes are few, so they are converged instead.
const double kNodeTolerance = 1e-10;

complex Direct(const PPair *pair, double radius, const Config &configure, int li, int lf, double si, double sf, double ji,
               double jf, int multL, char radType, double inEnergy, double levelEnergy, bool isChannelCapture) {
  ECIntegral integral(pair, configure, radius);
  return integral(li, lf, si, sf, ji, jf, multL, radType, inEnergy, levelEnergy, isChannelCapture);
}

/// The integral converged to kNodeTolerance, or at the direct integrals' own
/// tolerance where the quadrature cannot reach it (a GSL error or status).
complex Converged(const PPair *pair, double radius, const Config &configure, int li, int lf, double si, double sf, double ji,
                  double jf, int multL, char radType, double inEnergy, double levelEnergy, bool isChannelCapture) {
  try {
    ECIntegral integral(pair, configure, radius);
    integral.SetRelativeTolerance(kNodeTolerance);
    complex v = integral(li, lf, si, sf, ji, jf, multL, radType, inEnergy, levelEnergy, isChannelCapture);
    if (integral.Status() == 0) return v;
  } catch (...) {
  }
  return Direct(pair, radius, configure, li, lf, si, sf, ji, jf, multL, radType, inEnergy, levelEnergy, isChannelCapture);
}

}  // namespace

complex ECTable::Integral(const PPair *finalPair, double radius, const Config &configure,
                          int li, int lf, double si, double sf, double ji, double jf,
                          int multL, char radType, double inEnergy, double levelEnergy,
                          bool isChannelCapture, const ChannelTable *table, int finalPairNum) {
  ECTermStore &counts = configure.scratch->ecTerms;  // the run's counters (Report)
  const double threshold = finalPair->GetSepE() + finalPair->GetExE();
  // The channel energy by the channel table's expression, (E - S) - e (docs/upgrade/channel_table.md).
  const double outEnergy = (inEnergy - finalPair->GetSepE()) - finalPair->GetExE();
  CoulFunc coul(finalPair);
  // sqrt(P) at a_c: one more radius of the (l_i, E) block of the integrals (not in the rate mode: energies never recur)
  if (!RateMode(configure)) coul.SetRadialMemo(&configure.scratch->ForThisThread().ec.radial);
  complex reduction;
  // sqrt(P) or exp(-i phi) at a_c from the point's channel table where it holds the open channel (same function).
  // the entrance channel's 2J key: its waves depend on J_i when the pair has a spin-orbit term
  const int twoJi = finalPair->SpinOrbitTwoJ(ji);
  const ChannelTable::Row *row = table ? table->Find({finalPairNum, li, twoJi, radius, outEnergy}) : nullptr;
  // An unbound final state (levelEnergy above this pair's threshold, ECIntegral::IntegrateUnbound) is integrated
  // directly: its integral grows like (k_i - k_f)^-2 towards the photon threshold E = levelEnergy and is zero below
  // it, which a Lagrange stencil in ln(E_out) cannot follow.
  const bool unboundFinal = levelEnergy > threshold;
  // An entrance wave with an interior background is integrated directly, never from the table
  // (potential_background.md section 5.8): its amplitude follows the potential's resonances and the poles of R0,
  // whose width can be far below the spacing of the nodes (a narrow s-wave resonance of 3He + 4He at 0.29 MeV: up to
  // 21 % off with the table), and the dressed factors the table divided by vanish or jump there.  With the table
  // off the two routes are the same integral, so a run with or without --no-ec-table gives the same numbers for
  // such a wave.
  if (coul.HybridPotential() && azure2::physics::BackgroundOf(*coul.HybridPotential(), li, twoJi)) {
    ++counts.direct;
    return Direct(finalPair, radius, configure, li, lf, si, sf, ji, jf, multL, radType, inEnergy, levelEnergy,
                  isChannelCapture);
  }
  if (!unboundFinal && row && row->P > 0.0)
    reduction = isChannelCapture ? row->expHSP : complex(sqrt(row->P), 0.0);
  else if (unboundFinal || !(outEnergy > 0.0) ||
           !Reduction(radius, coul, li, twoJi, outEnergy, isChannelCapture, reduction)) {
    ++counts.direct;
    return Direct(finalPair, radius, configure, li, lf, si, sf, ji, jf, multL, radType, inEnergy, levelEnergy,
                  isChannelCapture);
  }

  const int perDecade = NodesPerDecade(configure);
  const double h = log(10.0) / perDecade;
  Term *term;
  {
    std::vector<double> key = TermKey(finalPair, radius, configure, coul, li, lf, si, sf, ji, jf, multL, radType,
                                      levelEnergy, isChannelCapture, perDecade);
    ECTermStore &store = configure.scratch->ecTerms;
    std::lock_guard<std::mutex> lock(store.mutex);
    std::unique_ptr<Term> &slot = store.terms[key];
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
      if (Reduction(radius, coul, li, twoJi, nodeOut, isChannelCapture, n)) {
        complex v = Converged(finalPair, radius, configure, li, lf, si, sf, ji, jf, multL, radType, nodeOut + threshold,
                              levelEnergy, isChannelCapture) / n;
        node.ok = std::isfinite(v.real()) && std::isfinite(v.imag());
        node.value = v;
      }
    } catch (...) {
      node.ok = false;
    }
    promises[s]->set_value(node);
    ++counts.nodes;
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
    ++counts.direct;
    return Direct(finalPair, radius, configure, li, lf, si, sf, ji, jf, multL, radType, inEnergy, levelEnergy,
                  isChannelCapture);
  }
  ++counts.lookups;
  complex result = sum * reduction;

  if (CheckMode()) {
    complex exact = Converged(finalPair, radius, configure, li, lf, si, sf, ji, jf, multL, radType, inEnergy, levelEnergy,
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
      std::lock_guard<std::mutex> lock(counts.checkMutex);
      if (dev > counts.maxDeviation) counts.maxDeviation = dev;
    }
  }
  return result;
}

long long ECTable::NodesComputed(const Config &configure) { return configure.scratch->ecTerms.nodes; }
long long ECTable::Lookups(const Config &configure) { return configure.scratch->ecTerms.lookups; }
long long ECTable::DirectCalls(const Config &configure) { return configure.scratch->ecTerms.direct; }

double ECTable::MaxCheckDeviation(const Config &configure) {
  ECTermStore &store = configure.scratch->ecTerms;
  std::lock_guard<std::mutex> lock(store.checkMutex);
  return store.maxDeviation;
}

std::string ECTable::Report(const Config &configure) {
  ECTermStore &store = configure.scratch->ecTerms;
  size_t terms;
  {
    std::lock_guard<std::mutex> lock(store.mutex);
    terms = store.terms.size();
  }
  std::ostringstream s;
  s << "External capture table: " << (store.lookups + store.direct) << " integrals, " << (store.nodes - store.nodesAtReset)
    << " computed at table nodes (" << terms << " terms)";
  if (store.direct > 0) s << ", " << store.direct << " directly";
  if (CheckMode()) {
    {
      std::lock_guard<std::mutex> lock(store.checkMutex);
      s << "; check: largest relative difference from the direct integrals " << store.maxDeviation;
    }
    // One line per term: li lf L type, E_out range, the worst point, and the
    // largest difference relative to the term's largest value.
    std::lock_guard<std::mutex> lock(store.mutex);
    for (auto &kv : store.terms) {
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

void ECTable::ResetCounters(const Config &configure) {
  ECTermStore &store = configure.scratch->ecTerms;
  store.lookups = 0;
  store.direct = 0;
  store.nodesAtReset = store.nodes.load();
  std::lock_guard<std::mutex> lock(store.checkMutex);
  store.maxDeviation = -1.0;
}

void ECTable::Clear(const Config &configure) {
  ECTermStore &store = configure.scratch->ecTerms;
  std::lock_guard<std::mutex> lock(store.mutex);
  store.terms.clear();
}
