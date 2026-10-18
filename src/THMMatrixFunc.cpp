#include "THMMatrixFunc.h"
#include "AZUREGrad.h"
#include "ChannelFunc.h"
#include "CNuc.h"
#include "Config.h"
#include "EPoint.h"
#include "JGroup.h"
#include "ALevel.h"
#include "AChannel.h"
#include "PPair.h"
#include "Constants.h"
#include "ThmLineshape.h"
#include "ThmDistortion.h"
#include "ThmDwVertex.h"
#include "ThmAngular.h"
#include "ThmExperiment.h"
#include "ThmVertexBoundary.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

namespace {

// An entrance bucket: channel spin s and orbital l -- (s, 0) for
// entranceL=coherent, (s, 2 l + k) for the two components of the DW vertex.
typedef std::pair<double, int> Bucket;
// Per bucket the entrance vertex of every level, v[la] = sum_{c_in} gamma_{la,c_in} M_l.
typedef std::map<Bucket, std::vector<complex>> BucketVertices;
// cbackground=: per bucket its entrance channel and the vertex without a width.
typedef std::map<Bucket, std::pair<int, complex>> BucketUnitVertices;

/// The kinematic factor of the model, fixed by how the HOES data were divided
/// out of the triple cross section (Config::ThmOptions::kinematics).
double ThmFluxFactor(const Config &configure, const PPair *entrancePair, const EPoint *point, double kf,
                     double muf) {
  double fluxFactor = 1.0;
  switch (configure.thm.kinematics) {
    case Config::ThmOptions::LA_COGNATA:  // k_f/mu_f (mrmpy)
      fluxFactor = kf / muf;
      break;
    case Config::ThmOptions::TRIPLE:
      break;
    case Config::ThmOptions::KF_THREE_BODY:
      fluxFactor = 1.0 / (muf * kf);
      break;
    case Config::ThmOptions::LAMBDA32: {
      // On-shell entrance momentum, as Pizzone et al. take it.  Undefined
      // below the entrance threshold; |E| keeps sub-threshold tails finite.
      double mui = entrancePair->GetRedMass() * uconv;
      double ki = std::sqrt(2.0 * mui * std::max(std::fabs(point->GetCMEnergy()), 1.0e-6)) / hbarc;
      fluxFactor = 1.0 / ki;
      break;
    }
  }
  return fluxFactor;
}

/*!
 * The entrance vertices of J group j for one pass (spectator-window node
 * `node`, -1 without a window; with the DW vertex `pass` indexes its nodes):
 * one bucket per entrance channel spin s and orbital l, summed incoherently
 * later.  With the quantization axis along p_xA only m_l = 0 enters, and once
 * the exit direction is integrated and the spin projections summed,
 * sum_{M,m_s} <s m_s l 0|J M><s m_s l' 0|J M> = delta_{ll'} (2J+1)/(2l+1)
 * removes the l cross terms.  `entranceL=coherent` keeps one bucket per s, as
 * mrmpy does.
 *
 * Vertex boundary B: ThmVertexBoundary (Config::ThmOptions::vertex), with the
 * on-shell L_c(E) = S_c(E) + i P_c(E) recovered from L_o = L_c - B_c.  With
 * coherentJ the vertex without a width of every bucket goes to unitVertex.
 */
void ThmEntranceVertices(const Config &configure, CNuc *compound, EPoint *point, JGroup *jg, int j, int aa,
                         int pass, int node, const ThmDwVertex *dw, const ThmDwVertex::At &dwAt, bool coherentJ,
                         BucketVertices &vbys, BucketUnitVertices &unitVertex) {
  const int numLevels = jg->NumLevels();
  const bool coherentL = configure.thm.coherentL;
  for (int ch = 1; ch <= jg->NumChannels(); ch++) {
    AChannel *c = jg->GetChannel(ch);
    if (c->GetPairNum() != aa) continue;
    std::vector<complex> *pwVertex = nullptr;
    if (!dw) {
      pwVertex = &vbys[std::make_pair(c->GetS(), coherentL ? 0 : c->GetL())];
      if (pwVertex->empty()) pwVertex->assign(numLevels + 1, complex(0.0, 0.0));
    }
    const ThmVertexBoundary boundaryOf(configure, compound->GetPair(aa), jg, ch);
    const complex onShellL = point->GetLoElement(j, ch) + c->GetBoundaryCondition();
    if (dw) {
      // Two components per (s, l): buckets 2 l + k.
      const int li = dw->LIndex(c->GetL());
      std::vector<complex> &vertex1 = vbys[std::make_pair(c->GetS(), 2 * c->GetL() + 1)];
      if (vertex1.empty()) vertex1.assign(numLevels + 1, complex(0.0, 0.0));
      std::vector<complex> &vertex0 = vbys[std::make_pair(c->GetS(), 2 * c->GetL())];
      if (vertex0.empty()) vertex0.assign(numLevels + 1, complex(0.0, 0.0));
      if (li < 0) continue;
      const complex *ak = &dwAt.a[((size_t)pass * dwAt.nl + li) * 2];
      const complex *dk = &dwAt.d[((size_t)pass * dwAt.nl + li) * 2];
      if (coherentJ) {
        complex b = boundaryOf.At(nullptr, onShellL);
        unitVertex[std::make_pair(c->GetS(), 2 * c->GetL())] = std::make_pair(ch, ak[0] * (b - 1.0) - dk[0]);
        unitVertex[std::make_pair(c->GetS(), 2 * c->GetL() + 1)] = std::make_pair(ch, ak[1] * (b - 1.0) - dk[1]);
      }
      for (int la = 1; la <= numLevels; la++) {
        ALevel *level = jg->GetLevel(la);
        if (!level->IsInRMatrix()) continue;
        complex boundary = boundaryOf.At(level, onShellL);
        double g = level->GetFitGamma(ch);
        vertex0[la] += g * (ak[0] * (boundary - 1.0) - dk[0]);
        vertex1[la] += g * (ak[1] * (boundary - 1.0) - dk[1]);
      }
      continue;
    }
    if (coherentJ) {
      complex b = boundaryOf.At(nullptr, onShellL);
      unitVertex[std::make_pair(c->GetS(), c->GetL())] =
          std::make_pair(ch, node < 0 ? point->GetThmFormFactor(j, ch, b) : point->GetThmFormFactor(j, ch, b, node));
    }
    for (int la = 1; la <= numLevels; la++) {
      ALevel *level = jg->GetLevel(la);
      if (!level->IsInRMatrix()) continue;
      complex boundary = boundaryOf.At(level, onShellL);
      (*pwVertex)[la] += level->GetFitGamma(ch) * (node < 0 ? point->GetThmFormFactor(j, ch, boundary)
                                                            : point->GetThmFormFactor(j, ch, boundary, node));
    }
  }
}

/// Widths below this are taken as a closed exit channel of the level.
const double kNegligibleExitWidth = 1.0e-12;

/*!
 * The resonant exit amplitude of one bucket into exit channel ch of J group
 * j: sum_{la,lap} gamma_{la,ch} [N_C(la)] A_{la,lap} v[lap] (act: the level's
 * index among the R-matrix levels; nc: the line-shape factors, or null).
 */
complex ThmExitAmplitude(const AMatrixFunc &func, JGroup *jg, int j, int ch, const std::vector<int> &act,
                         const std::vector<complex> &vertex, const std::vector<complex> *nc) {
  const int numLevels = jg->NumLevels();
  complex amp(0.0, 0.0);
  for (int la = 1; la <= numLevels; la++) {
    if (!jg->GetLevel(la)->IsInRMatrix()) continue;
    double gEx = jg->GetLevel(la)->GetFitGamma(ch);
    if (std::fabs(gEx) < kNegligibleExitWidth) continue;
    if (nc) {
      complex gExC = gEx * (*nc)[la];
      for (int lap = 1; lap <= numLevels; lap++) {
        if (!jg->GetLevel(lap)->IsInRMatrix()) continue;
        amp += gExC * func.GetAMatrixElement(j, act[la], act[lap]) * vertex[lap];
      }
      continue;
    }
    for (int lap = 1; lap <= numLevels; lap++) {
      if (!jg->GetLevel(lap)->IsInRMatrix()) continue;
      amp += gEx * func.GetAMatrixElement(j, act[la], act[lap]) * vertex[lap];
    }
  }
  return amp;
}

/*!
 * The energy-dependent factors of the point's model, applied in turn before
 * the folding: the segment's weight w(E) (<thm> weight[k]=, e.g. the
 * Coulomb-distortion factor R(E) of Mukhamedzhanov & Pang, PRC 99 (2019)
 * 064618, from a table), then the distortion factor R(E) of the point's THM
 * experiment (distortion=...; ThmDistortion.h).  Both apply to every point and
 * folding sub-point.
 */
void ThmApplyWeights(const Config &configure, const EPoint *point, double &sigma) {
  if (const ThmWeightTable *weight = point->GetThmWeight()) {
    bool outside = false;
    sigma *= (*weight)(point->GetCMEnergy(), &outside);
    if (outside && !weight->warned.exchange(true))
      configure.outStream << "WARNING: <thm> weight '" << weight->name << "' is evaluated at E = "
                          << point->GetCMEnergy() << " MeV, outside its table ["
                          << weight->e.front() << ", " << weight->e.back()
                          << "] MeV; the end value is used there (reported once)." << std::endl;
  }
  if (const ThmDistortion *distortion = point->GetThmDistortion()) {
    bool outside = false;
    sigma *= distortion->Weight(point->GetCMEnergy(), &outside);
    if (outside && !distortion->warned.exchange(true))
      configure.outStream << "WARNING: <thm> experiment[" << distortion->experiment
                          << "]: the distortion factor is evaluated at E = " << point->GetCMEnergy()
                          << " MeV, outside its grid; the end value is used there (reported once)." << std::endl;
  }
}

}  // namespace

THMMatrixFunc::THMMatrixFunc(CNuc *compound, const Config &configure) :
  AMatrixFunc(compound, configure) {}

/*!
 * HOES cross section of the modified R-matrix formalism (arbitrary units).
 * Assumes ClearMatrices/FillMatrices/InvertMatrices have already produced the
 * level matrix A (shared interior). Uses the entrance transfer form factors
 * M_l and exit penetrabilities stored on the point by CalcEDependentValues.
 */
void THMMatrixFunc::CalculateTHMCrossSection(EPoint *point) {
  int aa = compound()->GetPairNumFromKey(point->GetEntranceKey());
  int exitPairNum = compound()->GetPairNumFromKey(point->GetExitKey());
  PPair *entrancePair = compound()->GetPair(aa);
  PPair *exitPair = compound()->GetPair(exitPairNum);

  // The entrance c.m. energy sets the compound-system energy; the exit
  // channel energy is that minus the exit pair threshold. Below the exit
  // threshold there is no flux.
  double inEnergy = point->GetCMEnergy() + entrancePair->GetSepE() + entrancePair->GetExE();
  double exitEnergy = inEnergy - exitPair->GetSepE() - exitPair->GetExE();
  double muf = exitPair->GetRedMass() * uconv;  // MeV/c^2
  double kf = (exitEnergy > 0.0) ? std::sqrt(2.0 * muf * exitEnergy) / hbarc
                                 : 0.0;  // fm^-1
  if (kf == 0.0) {
    point->SetFitCrossSection(0.0);
    return;
  }
  const double fluxFactor = ThmFluxFactor(configure(), entrancePair, point, kf, muf);

  // Coulomb line shape of the spectator (<thm> experiment[...] lineshape=on,
  // ThmLineshape.h): the exit amplitude of level lambda -- the level that
  // decays to b + B -- carries N_C(zeta, E_lambda - E, Gamma_lambda), inside
  // the coherent level sum.  zeta depends on the energy and the exit pair only.
  const ThmLineshape *lineshape = point->GetThmLineshape();
  double zeta = 0.0;
  if (lineshape) {
    int light = exitPair->GetM(1) <= exitPair->GetM(2) ? 1 : 2;
    zeta = lineshape->Zeta(point->GetCMEnergy(), exitPair->GetZ(3 - light), exitPair->GetM(3 - light));
    if (!(lineshape->EsF(point->GetCMEnergy()) > 0.0) && !lineshape->warned.exchange(true))
      configure().outStream << "WARNING: <thm> experiment[" << lineshape->experiment << "] lineshape: at E = "
                            << point->GetCMEnergy() << " MeV the spectator has no energy left (E_sF <= 0); "
                            << "E_sF = 1 keV is used there (reported once)." << std::endl;
  }

  const ThmSpectatorWindow *window = point->GetThmSpectatorWindow();
  int numNodes = window ? point->NumThmPsNodes() : 0;

  // Distorted-wave entrance vertex (vertexModel=dw, ThmDwVertex.h): per node
  // (one: the spectatorAngle direction or the acceptance average) and entrance l, two
  // incoherent components M^(k)(B) = a_k (B - 1) - d_k replace M_l.
  const ThmDwVertex *dw = point->GetThmDwVertex();
  ThmDwVertex::At dwAt;
  if (dw) {
    dw->Evaluate(point->GetCMEnergy(), dwAt);
    if (dwAt.outside && !dw->warned.exchange(true))
      configure().outStream << "WARNING: <thm> experiment[" << dw->experiment
                            << "]: the DW vertex is evaluated at E = " << point->GetCMEnergy()
                            << " MeV, outside its grid; the end value is used there (reported once)." << std::endl;
    window = nullptr;
    numNodes = dwAt.nodes;
  }

  // Angular window of the exit pair (theta=thmin-thmax; ThmAngular.h): at a
  // fixed angle the J^pi groups and the entrance l of one channel spin
  // interfere, so the partial amplitudes of every J group are collected --
  // x = sqrt(K 2P) e^{i(omega - phi)} (amplitude of the exit channel), the exit
  // phase as in the T matrix of AMatrixFunc::CalculateTMatrix, with the
  // entrance vertex M_l in place of e^{i(omega-phi)} sqrt(P) -- and combined in
  // the Legendre sum after the loop, one set per spectator-window node.
  const ThmAngleWindow *angle = point->GetThmAngleWindow();
  std::vector<ThmPartialWave> waves;
  std::vector<std::vector<complex>> partial(angle ? std::max(numNodes, 1) : 0);

  // Coherent background of the experiment (cbackground=,
  // ThmCoherentBackground): c(E) times the bucket's vertex M_l (no width),
  // added to the resonant amplitude of its combinations before squaring.
  const ThmCoherentBackground *coherent = point->GetThmCoherent();
  const std::vector<double> &coherentValues = compound()->ThmCoherentValues();

  double sigma = 0.0;
  for (int j = 1; j <= compound()->NumJGroups(); j++) {
    JGroup *jg = compound()->GetJGroup(j);
    if (!jg->IsInRMatrix()) continue;
    int numLevels = jg->NumLevels();
    int numChannels = jg->NumChannels();
    double spinWeight = 2.0 * jg->GetJ() + 1.0;

    // Active-level index map (a_matrices_ / GetAMatrixElement are indexed by the
    // compacted in-R-matrix levels, so a level not in the R matrix is skipped).
    std::vector<int> act(numLevels + 1, 0);
    int na = 0;
    for (int la = 1; la <= numLevels; la++)
      if (jg->GetLevel(la)->IsInRMatrix()) act[la] = ++na;
    if (na == 0) continue;

    bool hasEntrance = false;
    for (int ch = 1; ch <= numChannels; ch++)
      if (jg->GetChannel(ch)->GetPairNum() == aa) hasEntrance = true;
    if (!hasEntrance) continue;  // this J group does not couple the entrance pair

    // N_C of each level (lineshape=on); the pole is the observed energy and
    // total width of the level (Brune).
    std::vector<complex> nc;
    if (lineshape) {
      nc.assign(numLevels + 1, complex(1.0, 0.0));
      for (int la = 1; la <= numLevels; la++) {
        ALevel *level = jg->GetLevel(la);
        if (!level->IsInRMatrix()) continue;
        nc[la] = ThmLineshapeFactor(zeta, level->GetFitE() - inEnergy,
                                    ThmLevelWidth(compound(), jg, level, configure()));
      }
    }

    // Spectator-momentum window (ps=..., ThmLineshape.h ThmSpectatorWindow):
    // the cross section -- not the amplitude -- is averaged over the accepted
    // directions at this energy (nodes p_k) with the normalized event weights
    // w_k (|phi|^2 d cos theta_cm), since each direction is a distinct final
    // state (Mukhamedzhanov et al. 2017 eq. 34).  Without a window one
    // pass with the stored single-node vertex (node -1), as before.
    bool coherentJ = false;
    if (coherent)
      for (const ThmCoherentBackground::Combo &cb : coherent->combos) coherentJ = coherentJ || cb.jGroup == j;

    const int passes = numNodes > 0 ? numNodes : 1;
    for (int pass = 0; pass < passes; pass++) {
      const int node = numNodes > 0 ? pass : -1;
      BucketVertices vbys;
      BucketUnitVertices unitVertex;
      ThmEntranceVertices(configure(), compound(), point, jg, j, aa, pass, node, dw, dwAt, coherentJ, vbys,
                          unitVertex);

      // Exit channels (incoherent), each with sqrt(2 P) folded in as 2 P outside.
      for (int ch = 1; ch <= numChannels; ch++) {
        AChannel *c = jg->GetChannel(ch);
        if (c->GetPairNum() != exitPairNum) continue;
        double sqrtPen = point->GetSqrtPenetrability(j, ch);
        double pex = sqrtPen * sqrtPen;  // P_l(E_exit); 0 if closed
        if (pex == 0.0) continue;

        double term = 0.0;
        for (BucketVertices::iterator it = vbys.begin(); it != vbys.end(); ++it) {
          complex amp = ThmExitAmplitude(*this, jg, j, ch, act, it->second, lineshape ? &nc : nullptr);
          if (coherentJ) {
            BucketUnitVertices::const_iterator u = unitVertex.find(it->first);
            const ThmCoherentBackground::Combo *cb = u != unitVertex.end() ? coherent->Find(j, u->second.first, ch)
                                                                           : nullptr;
            if (cb && cb->index + 2 * cb->form <= (int)coherentValues.size()) {
              const double *v = &coherentValues[cb->index];
              complex cE(v[0], v[1]);
              if (cb->form == 2) cE += complex(v[2], v[3]) * point->GetCMEnergy();
              amp += cE * u->second.second;
            }
          }
          term += std::norm(amp);  // |amp|^2, incoherent over (s, l)
          if (angle) {
            if (pass == 0) {
              ThmPartialWave w;
              w.J = jg->GetJ();
              w.s = it->first.first;
              w.l = it->first.second;
              w.sp = c->GetS();
              w.lp = c->GetL();
              waves.push_back(w);
            }
            partial[pass].push_back(std::sqrt(fluxFactor * 2.0 * pex) * point->GetExpCoulombPhase(j, ch) *
                                    point->GetExpHardSpherePhase(j, ch) * amp);
          }
        }
        if (angle) continue;
        if (dw)
          sigma += dwAt.weight[pass] * (spinWeight * fluxFactor * 2.0 * pex * term);
        else if (node < 0)
          sigma += spinWeight * fluxFactor * 2.0 * pex * term;
        else
          sigma += point->GetThmPsWeight(node) * (spinWeight * fluxFactor * 2.0 * pex * term);
      }
    }
  }

  if (angle) {
    // sum_spins |F|^2 = (1/pi) sum_L b_L P_L(cos theta), averaged over the
    // window (and over the spectator-window nodes with their weights).
    // (theta= with vertexModel=dw is refused, CheckThmExperiments: no DW weights here.)
    std::vector<double> b;
    for (size_t k = 0; k < partial.size(); k++) {
      ThmLegendreCoefficients(waves, partial[k], b);
      sigma += (numNodes > 0 ? point->GetThmPsWeight((int)k) : 1.0) * angle->Mean(b);
    }
  }

  ThmApplyWeights(configure(), point, sigma);
  point->SetFitCrossSection(sigma);
}

// ===========================================================================
// Reverse-mode (adjoint) derivative of the HOES cross section.
// ===========================================================================

ThmAdjointContext BuildThmAdjointContext(CNuc *compound, const Config &configure, const vector_matrix_r *shiftDeriv,
                                         bool lineshape) {
  ThmAdjointContext ctx;
  ctx.shiftDeriv = shiftDeriv;
  const bool useGSL = !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC);
  const int nJ = compound->NumJGroups();
  ctx.lowestShiftDeriv.assign(nJ, std::vector<double>());
  if (lineshape) ctx.widthDeriv.assign(nJ, std::vector<std::vector<double>>());
  for (int j = 1; j <= nJ; j++) {
    JGroup *jg = compound->GetJGroup(j);
    if (!jg->IsInRMatrix()) continue;
    const int numChannels = jg->NumChannels();
    // vertex=constant: B = S_c(E_min) of the lowest level (ThmVertexBoundary).
    if (ThmVertexBoundary::RuleOf(configure) == ThmVertexBoundary::LOWEST_LEVEL_SHIFT) {
      ctx.lowestShiftDeriv[j - 1].assign(numChannels, 0.0);
      int lowest = 0;
      double eMin = 0.0;
      for (int la = 1; la <= jg->NumLevels(); la++) {
        ALevel *level = jg->GetLevel(la);
        if (!level->IsInRMatrix()) continue;
        if (!lowest || level->GetFitE() < eMin) {
          eMin = level->GetFitE();
          lowest = la;
        }
      }
      if (lowest)
        for (int ch = 1; ch <= numChannels; ch++) {
          AChannel *c = jg->GetChannel(ch);
          PPair *pair = compound->GetPair(c->GetPairNum());
          if (pair->GetPType() != 0) continue;  // only a particle pair has an entrance vertex
          ctx.lowestShiftDeriv[j - 1][ch - 1] =
              ChannelFunc(pair, useGSL).ShiftDerivative(c->GetL(), eMin - (pair->GetSepE() + pair->GetExE()));
        }
    }
    // lineshape=on: dGamma/dE and dGamma/dgamma_c of every level, by central
    // differences of ThmLevelWidth on a copy of the level (it depends on that
    // level's own energy and widths only).
    if (lineshape) {
      ctx.widthDeriv[j - 1].assign(jg->NumLevels(), std::vector<double>());
      for (int la = 1; la <= jg->NumLevels(); la++) {
        ALevel *level = jg->GetLevel(la);
        if (!level->IsInRMatrix()) continue;
        std::vector<double> &d = ctx.widthDeriv[j - 1][la - 1];
        d.assign(numChannels + 1, 0.0);
        ALevel copy(*level);
        {
          const double x0 = level->GetFitE(), h = 1.0e-6 * (std::fabs(x0) + 1.0);
          copy.SetFitE(x0 + h);
          double wp = ThmLevelWidth(compound, jg, &copy, configure);
          copy.SetFitE(x0 - h);
          double wm = ThmLevelWidth(compound, jg, &copy, configure);
          copy.SetFitE(x0);
          d[0] = (wp - wm) / (2.0 * h);
        }
        for (int ch = 1; ch <= numChannels; ch++) {
          const double x0 = level->GetFitGamma(ch), h = 1.0e-6 * (std::fabs(x0) + 1.0);
          copy.SetFitGamma(ch, x0 + h);
          double wp = ThmLevelWidth(compound, jg, &copy, configure);
          copy.SetFitGamma(ch, x0 - h);
          double wm = ThmLevelWidth(compound, jg, &copy, configure);
          copy.SetFitGamma(ch, x0);
          d[ch] = (wp - wm) / (2.0 * h);
        }
      }
    }
  }
  return ctx;
}

/*!
 * The adjoint mirrors CalculateTHMCrossSection step by step (see there for
 * the physics); a change to one must be made to the other, and
 * tests/pyazr/thm_adjoint_test.py checks the derivatives against differences.  In order:
 *
 *  1. forward: per J group and pass (spectator-window node / DW node) the
 *     bucket vertices v, per exit channel and bucket the amplitude
 *       amp = sum_{la,lap} gEx_la [N_C(la)] A_{la,lap} v_lap  [+ c(E) u],
 *     and sigma (as the forward, the same sums in the same order);
 *  2. ampBar = d sigma/d amp (cotangent convention of PointAdjoint:
 *     zbar = d/dRe z + i d/dIm z): 2 W amp for an incoherent |amp|^2 term of
 *     weight W, or through the Legendre sum of the angular window;
 *  3. amp -> gEx (explicit, also for a width under the cutoff the forward
 *     skips, whose derivative is not zero), A (-> LevelMatrixAdjoint), v,
 *     N_C, the coherent background c(E) and its unit vertex u;
 *  4. v_la = sum_c gamma_{la,c} M(B) -> gamma_{la,c}, and through the
 *     boundary B (ThmVertexBoundary): S_c(E_la) under vertex=perlevel with
 *     Brune (dS/dE from ctx.shiftDeriv), S_c(E_min) under vertex=constant
 *     (the lowest level's energy, ctx.lowestShiftDeriv); M is linear in B;
 *  5. N_C(zeta, E_la - E, Gamma_la) -> E_la, and Gamma_la -> E_la, gamma_la,c
 *     (ctx.widthDeriv).
 * The kinematic factor, the exit penetrability, the weights w(E) and R(E),
 * and the vertex pieces j_l, rho j_l', the Coulomb term and the DW a_k, d_k
 * depend on the energy only.
 */
bool THMMatrixFunc::THMPointAdjoint(EPoint *point, double fitBar, GradAccum &accum, const ThmAdjointContext &ctx,
                                    std::vector<double> *coherentBar, double *sigmaOut) {
  if (configure().paramMask & Config::USE_RMC_FORMALISM) return false;
  int aa = compound()->GetPairNumFromKey(point->GetEntranceKey());
  int exitPairNum = compound()->GetPairNumFromKey(point->GetExitKey());
  PPair *entrancePair = compound()->GetPair(aa);
  PPair *exitPair = compound()->GetPair(exitPairNum);

  double inEnergy = point->GetCMEnergy() + entrancePair->GetSepE() + entrancePair->GetExE();
  double exitEnergy = inEnergy - exitPair->GetSepE() - exitPair->GetExE();
  double muf = exitPair->GetRedMass() * uconv;
  double kf = (exitEnergy > 0.0) ? std::sqrt(2.0 * muf * exitEnergy) / hbarc : 0.0;
  if (sigmaOut) *sigmaOut = 0.0;
  if (kf == 0.0) return true;
  const double fluxFactor = ThmFluxFactor(configure(), entrancePair, point, kf, muf);
  // The energy-only weights w(E) R(E), as one factor (ThmApplyWeights).
  double weightFactor = 1.0;
  ThmApplyWeights(configure(), point, weightFactor);

  const ThmLineshape *lineshape = point->GetThmLineshape();
  double zeta = 0.0;
  if (lineshape) {
    int light = exitPair->GetM(1) <= exitPair->GetM(2) ? 1 : 2;
    zeta = lineshape->Zeta(point->GetCMEnergy(), exitPair->GetZ(3 - light), exitPair->GetM(3 - light));
    if (zeta != 0.0 && ctx.widthDeriv.empty()) return false;  // no width derivatives to hand
  }

  const ThmSpectatorWindow *window = point->GetThmSpectatorWindow();
  int numNodes = window ? point->NumThmPsNodes() : 0;
  const ThmDwVertex *dw = point->GetThmDwVertex();
  ThmDwVertex::At dwAt;
  if (dw) {
    dw->Evaluate(point->GetCMEnergy(), dwAt);
    window = nullptr;
    numNodes = dwAt.nodes;
  }
  const ThmAngleWindow *angle = point->GetThmAngleWindow();
  std::vector<ThmPartialWave> waves;
  std::vector<std::vector<complex>> partial(angle ? std::max(numNodes, 1) : 0);
  const ThmCoherentBackground *coherent = point->GetThmCoherent();
  const std::vector<double> &coherentValues = compound()->ThmCoherentValues();
  const ThmVertexBoundary::Rule rule = ThmVertexBoundary::RuleOf(configure());

  // One block per (J group, pass): what the backward sweep needs again.
  struct Block {
    int j, pass, node;
    std::vector<int> act, levelOf;
    std::vector<complex> nc;  // N_C(la), lineshape=on
    BucketVertices vbys, vBar;
    BucketUnitVertices unit;
    std::map<Bucket, complex> uBar;
    std::vector<complex> ncBar;
  };
  // One record per (block, exit channel, bucket) amplitude.
  struct Rec {
    int block, ch;
    Bucket key;
    complex amp;
    double w;        // weight of |amp|^2 in sigma (no angular window)
    complex pref;    // x = pref amp (angular window)
    int slot;        // index of x in partial[pass]
    const ThmCoherentBackground::Combo *cb;
    complex unitVertex, cE;
  };
  std::vector<Block> blocks;
  std::vector<Rec> recs;

  double sigma = 0.0;
  for (int j = 1; j <= compound()->NumJGroups(); j++) {
    JGroup *jg = compound()->GetJGroup(j);
    if (!jg->IsInRMatrix()) continue;
    int numLevels = jg->NumLevels();
    int numChannels = jg->NumChannels();
    double spinWeight = 2.0 * jg->GetJ() + 1.0;
    std::vector<int> act(numLevels + 1, 0), levelOf(1, 0);
    int na = 0;
    for (int la = 1; la <= numLevels; la++)
      if (jg->GetLevel(la)->IsInRMatrix()) {
        act[la] = ++na;
        levelOf.push_back(la);
      }
    if (na == 0) continue;
    bool hasEntrance = false;
    for (int ch = 1; ch <= numChannels; ch++)
      if (jg->GetChannel(ch)->GetPairNum() == aa) hasEntrance = true;
    if (!hasEntrance) continue;

    std::vector<complex> nc;
    if (lineshape) {
      nc.assign(numLevels + 1, complex(1.0, 0.0));
      for (int la = 1; la <= numLevels; la++) {
        ALevel *level = jg->GetLevel(la);
        if (!level->IsInRMatrix()) continue;
        nc[la] = ThmLineshapeFactor(zeta, level->GetFitE() - inEnergy,
                                    ThmLevelWidth(compound(), jg, level, configure()));
      }
    }
    bool coherentJ = false;
    if (coherent)
      for (const ThmCoherentBackground::Combo &cb : coherent->combos) coherentJ = coherentJ || cb.jGroup == j;

    const int passes = numNodes > 0 ? numNodes : 1;
    for (int pass = 0; pass < passes; pass++) {
      const int node = numNodes > 0 ? pass : -1;
      blocks.push_back(Block());
      Block &b = blocks.back();
      b.j = j;
      b.pass = pass;
      b.node = node;
      b.act = act;
      b.levelOf = levelOf;
      b.nc = nc;
      ThmEntranceVertices(configure(), compound(), point, jg, j, aa, pass, node, dw, dwAt, coherentJ, b.vbys, b.unit);

      for (int ch = 1; ch <= numChannels; ch++) {
        AChannel *c = jg->GetChannel(ch);
        if (c->GetPairNum() != exitPairNum) continue;
        double sqrtPen = point->GetSqrtPenetrability(j, ch);
        double pex = sqrtPen * sqrtPen;
        if (pex == 0.0) continue;
        double term = 0.0;
        double passWeight = dw ? dwAt.weight[pass] : node < 0 ? 1.0 : point->GetThmPsWeight(node);
        for (BucketVertices::iterator it = b.vbys.begin(); it != b.vbys.end(); ++it) {
          Rec r;
          r.block = (int)blocks.size() - 1;
          r.ch = ch;
          r.key = it->first;
          r.cb = nullptr;
          r.slot = -1;
          r.w = 0.0;
          complex amp = ThmExitAmplitude(*this, jg, j, ch, act, it->second, lineshape ? &b.nc : nullptr);
          if (coherentJ) {
            BucketUnitVertices::const_iterator u = b.unit.find(it->first);
            const ThmCoherentBackground::Combo *cb = u != b.unit.end() ? coherent->Find(j, u->second.first, ch)
                                                                       : nullptr;
            if (cb && cb->index + 2 * cb->form <= (int)coherentValues.size()) {
              const double *v = &coherentValues[cb->index];
              complex cE(v[0], v[1]);
              if (cb->form == 2) cE += complex(v[2], v[3]) * point->GetCMEnergy();
              amp += cE * u->second.second;
              r.cb = cb;
              r.cE = cE;
              r.unitVertex = u->second.second;
            }
          }
          r.amp = amp;
          term += std::norm(amp);
          if (angle) {
            if (pass == 0) {
              ThmPartialWave w;
              w.J = jg->GetJ();
              w.s = it->first.first;
              w.l = it->first.second;
              w.sp = c->GetS();
              w.lp = c->GetL();
              waves.push_back(w);
            }
            r.pref = std::sqrt(fluxFactor * 2.0 * pex) * point->GetExpCoulombPhase(j, ch) *
                     point->GetExpHardSpherePhase(j, ch);
            r.slot = (int)partial[pass].size();
            partial[pass].push_back(r.pref * amp);
          } else {
            r.w = passWeight * spinWeight * fluxFactor * 2.0 * pex;
          }
          recs.push_back(r);
        }
        if (angle) continue;
        if (dw)
          sigma += dwAt.weight[pass] * (spinWeight * fluxFactor * 2.0 * pex * term);
        else if (node < 0)
          sigma += spinWeight * fluxFactor * 2.0 * pex * term;
        else
          sigma += point->GetThmPsWeight(node) * (spinWeight * fluxFactor * 2.0 * pex * term);
      }
    }
  }
  std::vector<std::vector<complex>> xBar(partial.size());
  if (angle) {
    std::vector<double> bL;
    for (size_t k = 0; k < partial.size(); k++) {
      const double pw = numNodes > 0 ? point->GetThmPsWeight((int)k) : 1.0;
      ThmLegendreCoefficients(waves, partial[k], bL);
      sigma += pw * angle->Mean(bL);
      // sigma_k = (1/pi) sum_L b_L meanP[L]  (ThmAngleWindow::Mean)
      std::vector<double> wL(angle->meanP.size());
      for (size_t L = 0; L < wL.size(); L++) wL[L] = fitBar * weightFactor * pw * angle->meanP[L] / M_PI;
      ThmLegendreAdjoint(waves, partial[k], wL, xBar[k]);
    }
  }
  sigma *= weightFactor;
  if (sigmaOut) *sigmaOut = sigma;
  if (fitBar == 0.0) return true;

  // --- Backward: amplitudes -> gEx, A, v, N_C, coherent background. ---
  std::vector<matrix_c> aBar(compound()->NumJGroups());
  for (Block &b : blocks) {
    const int na = (int)b.levelOf.size() - 1;
    if (aBar[b.j - 1].empty()) aBar[b.j - 1].assign(na, vector_c(na, complex(0.0, 0.0)));
    for (BucketVertices::iterator it = b.vbys.begin(); it != b.vbys.end(); ++it)
      b.vBar[it->first].assign(it->second.size(), complex(0.0, 0.0));
    if (!b.nc.empty()) b.ncBar.assign(b.nc.size(), complex(0.0, 0.0));
  }
  std::vector<complex> av;
  for (const Rec &r : recs) {
    Block &b = blocks[r.block];
    const int j = b.j;
    JGroup *jg = compound()->GetJGroup(j);
    const int numLevels = jg->NumLevels();
    const matrix_c &A = this->GetAMatrix(j);
    if (A.empty()) continue;
    complex ampBar = r.slot >= 0 ? std::conj(r.pref) * xBar[b.pass][r.slot]
                                 : fitBar * weightFactor * r.w * 2.0 * r.amp;
    if (ampBar == complex(0.0, 0.0)) continue;
    const std::vector<complex> &v = b.vbys[r.key];
    std::vector<complex> &vb = b.vBar[r.key];
    matrix_c &Ab = aBar[j - 1];
    // (A v)_la over the active levels.
    av.assign(numLevels + 1, complex(0.0, 0.0));
    for (int la = 1; la <= numLevels; la++) {
      if (!b.act[la]) continue;
      const vector_c &Arow = A[b.act[la] - 1];
      complex s(0.0, 0.0);
      for (int lap = 1; lap <= numLevels; lap++)
        if (b.act[lap]) s += Arow[b.act[lap] - 1] * v[lap];
      av[la] = s;
    }
    for (int la = 1; la <= numLevels; la++) {
      if (!b.act[la]) continue;
      double gEx = jg->GetLevel(la)->GetFitGamma(r.ch);
      complex ncl = b.nc.empty() ? complex(1.0, 0.0) : b.nc[la];
      accum.AddGamma(j, la, r.ch, std::real(std::conj(ampBar) * ncl * av[la]));
      if (std::fabs(gEx) < kNegligibleExitWidth) continue;  // out of the forward sum
      if (!b.nc.empty()) b.ncBar[la] += std::conj(gEx * av[la]) * ampBar;
      const complex cb = std::conj(gEx * ncl) * ampBar;  // cotangent of (A v)_la
      const vector_c &Arow = A[b.act[la] - 1];
      vector_c &Abrow = Ab[b.act[la] - 1];
      for (int lap = 1; lap <= numLevels; lap++) {
        if (!b.act[lap]) continue;
        Abrow[b.act[lap] - 1] += cb * std::conj(v[lap]);
        vb[lap] += std::conj(Arow[b.act[lap] - 1]) * cb;
      }
    }
    if (r.cb) {
      // amp += c(E) u:  c(E) = (v0 + i v1) [+ (v2 + i v3) E].
      const complex q = std::conj(r.unitVertex) * ampBar;
      if (coherentBar && r.cb->index + 2 * r.cb->form <= (int)coherentBar->size()) {
        (*coherentBar)[r.cb->index] += std::real(q);
        (*coherentBar)[r.cb->index + 1] += std::imag(q);
        if (r.cb->form == 2) {
          (*coherentBar)[r.cb->index + 2] += std::real(q) * point->GetCMEnergy();
          (*coherentBar)[r.cb->index + 3] += std::imag(q) * point->GetCMEnergy();
        }
      }
      b.uBar[r.key] += std::conj(r.cE) * ampBar;
    }
  }

  // --- Vertices -> gamma_in, and the boundary B -> level energies. ---
  const bool coherentL = configure().thm.coherentL;
  for (Block &b : blocks) {
    const int j = b.j;
    JGroup *jg = compound()->GetJGroup(j);
    const int numLevels = jg->NumLevels();
    // Bbar (real B) to the level energy it depends on.
    auto routeB = [&](int ch, const ALevel *level, int la, int lowest, double bBar) {
      if (bBar == 0.0) return;
      if (rule == ThmVertexBoundary::LEVEL_SHIFT && level) {
        if (ctx.shiftDeriv) accum.AddE(j, la, bBar * (*ctx.shiftDeriv)[j - 1][la - 1][ch - 1]);
      } else if (rule == ThmVertexBoundary::LOWEST_LEVEL_SHIFT && lowest > 0) {
        if (j - 1 < (int)ctx.lowestShiftDeriv.size() && ch - 1 < (int)ctx.lowestShiftDeriv[j - 1].size())
          accum.AddE(j, lowest, bBar * ctx.lowestShiftDeriv[j - 1][ch - 1]);
      }
    };
    for (int ch = 1; ch <= jg->NumChannels(); ch++) {
      AChannel *c = jg->GetChannel(ch);
      if (c->GetPairNum() != aa) continue;
      const ThmVertexBoundary boundaryOf(configure(), compound()->GetPair(aa), jg, ch);
      const int lowest = boundaryOf.LowestLevel();
      const complex onShellL = point->GetLoElement(j, ch) + c->GetBoundaryCondition();
      if (dw) {
        const int li = dw->LIndex(c->GetL());
        if (li < 0) continue;
        const complex *ak = &dwAt.a[((size_t)b.pass * dwAt.nl + li) * 2];
        const complex *dk = &dwAt.d[((size_t)b.pass * dwAt.nl + li) * 2];
        for (int k = 0; k < 2; k++) {
          const Bucket key = std::make_pair(c->GetS(), 2 * c->GetL() + k);
          std::map<Bucket, complex>::const_iterator ub = b.uBar.find(key);
          if (ub != b.uBar.end() && b.unit[key].first == ch && !boundaryOf.OnShell())
            routeB(ch, nullptr, 0, lowest, std::real(std::conj(ub->second) * ak[k]));
          BucketVertices::iterator vbIt = b.vBar.find(key);
          if (vbIt == b.vBar.end()) continue;
          const std::vector<complex> &vb = vbIt->second;
          for (int la = 1; la <= numLevels; la++) {
            ALevel *level = jg->GetLevel(la);
            if (!level->IsInRMatrix()) continue;
            const complex boundary = boundaryOf.At(level, onShellL);
            const double g = level->GetFitGamma(ch);
            accum.AddGamma(j, la, ch, std::real(std::conj(vb[la]) * (ak[k] * (boundary - 1.0) - dk[k])));
            if (!boundaryOf.OnShell()) routeB(ch, level, la, lowest, std::real(std::conj(vb[la]) * (g * ak[k])));
          }
        }
        continue;
      }
      const Bucket key = std::make_pair(c->GetS(), coherentL ? 0 : c->GetL());
      // M(B) = (B - 1) j_l - rho j_l' + C_l is linear in B: dM/dB = M(B + 1) - M(B).
      auto formFactor = [&](complex boundary) {
        return b.node < 0 ? point->GetThmFormFactor(j, ch, boundary) : point->GetThmFormFactor(j, ch, boundary, b.node);
      };
      std::map<Bucket, complex>::const_iterator ub = b.uBar.find(std::make_pair(c->GetS(), c->GetL()));
      if (ub != b.uBar.end() && b.unit[ub->first].first == ch && !boundaryOf.OnShell()) {
        const complex bf = boundaryOf.At(nullptr, onShellL);
        routeB(ch, nullptr, 0, lowest, std::real(std::conj(ub->second) * (formFactor(bf + 1.0) - formFactor(bf))));
      }
      BucketVertices::iterator vbIt = b.vBar.find(key);
      if (vbIt == b.vBar.end()) continue;
      const std::vector<complex> &vb = vbIt->second;
      for (int la = 1; la <= numLevels; la++) {
        ALevel *level = jg->GetLevel(la);
        if (!level->IsInRMatrix()) continue;
        const complex boundary = boundaryOf.At(level, onShellL);
        const complex m = formFactor(boundary);
        accum.AddGamma(j, la, ch, std::real(std::conj(vb[la]) * m));
        if (!boundaryOf.OnShell())
          routeB(ch, level, la, lowest,
                 std::real(std::conj(vb[la]) * (level->GetFitGamma(ch) * (formFactor(boundary + 1.0) - m))));
      }
    }

    // N_C(zeta, x = E_la - E, Gamma_la): d ln N/dx = zeta (G/2 - i x)/r^2,
    // d ln N/dG = -zeta (x/2 + i G/4)/r^2, r^2 = x^2 + G^2/4 (G > 0).
    if (!b.ncBar.empty() && zeta != 0.0) {
      for (int la = 1; la <= numLevels; la++) {
        ALevel *level = jg->GetLevel(la);
        if (!level->IsInRMatrix() || b.ncBar[la] == complex(0.0, 0.0)) continue;
        const double x = level->GetFitE() - inEnergy;
        const double G = ThmLevelWidth(compound(), jg, level, configure());
        complex dlnx, dlnG(0.0, 0.0);
        if (G > 0.0) {
          const double r2 = x * x + 0.25 * G * G;
          dlnx = zeta * complex(0.5 * G, -x) / r2;
          dlnG = -zeta * complex(0.5 * x, 0.25 * G) / r2;
        } else {
          dlnx = x != 0.0 ? complex(0.0, -zeta / x) : complex(0.0, 0.0);
        }
        const complex n = b.nc[la];
        const double dx = std::real(std::conj(b.ncBar[la]) * (n * dlnx));
        const double dG = std::real(std::conj(b.ncBar[la]) * (n * dlnG));
        accum.AddE(j, la, dx);
        const std::vector<double> &wd = ctx.widthDeriv[j - 1][la - 1];
        if (dG != 0.0 && !wd.empty()) {
          accum.AddE(j, la, dG * wd[0]);
          for (int ch = 1; ch <= jg->NumChannels(); ch++) accum.AddGamma(j, la, ch, dG * wd[ch]);
        }
      }
    }
  }

  // --- Level matrices -> E, gamma. ---
  for (Block &b : blocks) {
    matrix_c &Ab = aBar[b.j - 1];
    if (Ab.empty()) continue;
    LevelMatrixAdjoint(point, b.j, Ab, inEnergy, ctx.shiftDeriv, accum, &b.levelOf);
    Ab.clear();  // once per J group (the blocks of its passes share it)
  }
  return true;
}
