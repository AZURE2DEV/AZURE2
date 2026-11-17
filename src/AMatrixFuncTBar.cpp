// AMatrixFunc, cross section -> T: the T-matrix cotangents of one point's observable (PointTBar), by kind of
// observable (angle-integrated, total, beta-delayed, A_y of particle and capture exits, phase shift, differential and
// Legendre coefficients), and the per-pathway Legendre products of the differential sum (TBarProducts).  The next
// step, T -> level energies and reduced widths, is AMatrixFuncAdjoint.cpp; the forward T matrix AMatrixFunc.cpp.
#include "AMatrixFunc.h"
#include "RunDiagnostics.h"
#include <mutex>
#include "azure2/support/Environment.h"
#include "Profile.h"
#include "AngCoeff.h"
#include "PolarizationFunc.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "Unit.h"
#include "AZUREGrad.h"
#include "azure2/physics/ExternalCapture.h"
#include <stdexcept>
#include "JGroup.h"
#include "ALevel.h"
#include "AChannel.h"
#include "PPair.h"
#include <assert.h>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif

/*!
 * Reverse mode of one unit's observable: PointTBar (cross section -> T) and TBarAdjoint (T -> level energies and
 * reduced widths); the gradient over the energy tasks (GradEnergyTasks.cpp) sums the T cotangents of the units of a
 * T-matrix group between the two.
 *
 * Adjoint convention (consistent throughout, see PLAN.md "Gotchas"):
 *   For a complex intermediate z and the real scalar lnL we store the cotangent
 *     zbar := dlnL/dRe(z) + i dlnL/dIm(z).
 *   - A linear map z = c * a (c complex constant) back-propagates as
 *         abar += conj(c) * zbar.
 *   - The contribution of z to a *real* parameter p is
 *         dlnL/dp += Re( conj(zbar) * dz/dp ).
 *   - The holomorphic inverse A = M^{-1} back-propagates as
 *         Mbar = -A^H * Abar * A^H.
 *
 * Forward chain (per point, standard A-matrix):
 *   T(k,m) = [tphase] - uphase ( [1] + umatrix )         ([..] only if chNum==chpNum)
 *   umatrix = sum_{l,l'} 2i sqrtP_en sqrtP_ex gamma_en,l gamma_ex,l' A_j[l][l']
 *   A_j = M_j^{-1},   M_j[r][c] = (E_r - inE) delta_rc - sum_ch gamma_r,ch gamma_c,ch L_ch
 * with the cross-section -> T step depending on whether the point is
 * angle-integrated (quadratic |T_temp|^2 sum) or differential (Legendre
 * interference sum + elastic Coulomb interference).
 */

/// A gate of PointTBar no project reaches (why in the message): a broken invariant, never a fallback point.
[[noreturn]] static void UnreachableGate(const char *w) { throw std::logic_error(std::string("PointTBar: ") + w); }
/*!
 * Cross section -> T: the T-matrix cotangents of one point's observable
 * (AMatrixFunc.h).  Reads the T-matrix held by this object and the
 * angle-dependent quantities of `point`.  False: the finite-difference fallback (park_chi2_fallback_test.cpp).
 */
bool AMatrixFunc::PointTBar(Site *point, const azure2::Results &results, double fitBar, TMatrixCotangent &bar,
                            int xsComponent, bool ayCrossSection, const Site &rows, int rowSlot) {
  AZUREProfile::Timer profileTimer(configure().profile, AZUREProfile::kGradTBar);
  const Site *angles = &rows;  // the angular row (Legendre rows, coefficient order) is the reader's; the rest the site's

  // --- Common gating. ---
  if (configure().paramMask & Config::USE_RMC_FORMALISM) return false;
  if (point->IsAngularDist() && !angles->IsCoefficient()) UnreachableGate("an extrapolation's computed coefficients");
  // xsComponent < 0 asks for A_0 = sigma_int / 4 pi of an angular-distribution point: the angle-integrated
  // cotangent scaled by 1 / 4 pi (GradTargetEffectAdjoint, normalised coefficients under a target effect).
  const bool a0Mode = (xsComponent < 0);
  if (a0Mode) xsComponent = 0;
  if (xsComponent != 0 && (point->IsDifferential() || point->IsPhase()))  // field multipole: Quantities.cpp
    UnreachableGate("an E1/E2 part exists for sigma, sfactor and modified_sfactor only, not differential or phase");

  int aa = compound()->GetPairNumFromKey(point->GetEntranceKey());
  int exitPairNum = compound()->GetPairNumFromKey(point->GetExitKey());
  // Beta-delayed particle emission is supported below; its observable is a plain
  // incoherent sum over pathways, so only the T cotangents differ from the
  // cross-section case. E1/E2 component selection has no meaning for it.
  const bool isBetaDelayed = (compound()->GetPair(aa)->GetPType() == 20);
  if (isBetaDelayed && xsComponent != 0) return false;

  // Active-level compaction unsupported here (original idx == active; ClearMatrices records it).
  if (!all_levels_active_) UnreachableGate("a level outside the R matrix: every level read from <levels> is in it");

  // The exit of the point's exit pair (CalculateTMatrix, called for the point's group first, throws without one).
  const azure2::EntrancePathways &P = compound()->Pathways().Entrance(aa);
  const int ir = P.ExitOf(exitPairNum);
  if (!ir) UnreachableGate("no exit for the point's pair, after CalculateTMatrix of its group");
  const azure2::ExitRow &theDecay = P.exits[ir - 1];

  const double geom = point->GetGeometricalFactor();
  const double i1i2 = compound()->GetPair(aa)->GetI1I2Factor();
  const bool isDiff = point->IsDifferential();
  const bool isPhase = point->IsPhase();
  if (isPhase && aa != exitPairNum) UnreachableGate("an inelastic phase shift (validator and Build refuse it)");
  // The cross section carried by an analyzing-power point is the differential
  // one at the point's angle (GenMatrixFunc computes it before A_y and a
  // target-effect sub-point keeps it); anything else is not supported.
  if (ayCrossSection && (!point->IsAnalyzingPower() || !isDiff || isPhase || xsComponent != 0)) return false;

  // Electron screening (GenMatrixFunc::CalculateCrossSection): a screened observable is f times the bare
  // one, so the bare one's cotangent is f fitBar, and dY/dU_e = k Y per MeV, 1e-6 k Y per eV, with Y the
  // screened value this point holds from its forward calculation.  A phase shift and an A_y value are
  // not screened (an A_y sub-point's cross section is: ayCrossSection).
  bar.screen.clear();
  if (point->GetScreeningParam() > 0 && point->GetScreeningK() != 0.0 && !isPhase && !angles->IsCoefficient() &&
      !(point->IsAnalyzingPower() && !ayCrossSection)) {
    double y = point->ValueFor(*angles, rowSlot, results);
    if (xsComponent == 1) y = point->GetFitE1CrossSection(results);
    else if (xsComponent == 2) y = point->GetFitE2CrossSection(results);
    bar.screen.push_back({point->GetScreeningParam(), fitBar * 1.0e-6 * point->GetScreeningK() * y});
    fitBar *= point->ScreeningFactor();
  }

  // === Cross-section -> T:  build the T-matrix cotangents. ===
  // tBar  : cotangents of the internal T(k,m).
  // ecBar : cotangents of the external-capture T(k,ecm).
  int nK = theDecay.num_groups;
  std::vector<vector_c> &tBar = bar.t, &ecBar = bar.ec;
  tBar.assign(nK, vector_c());
  ecBar.assign(nK, vector_c());
  for (int k = 1; k <= nK; k++) {
    tBar[k - 1].assign(P.Group(theDecay, k).num_internal, complex(0.0, 0.0));
    ecBar[k - 1].assign(P.Group(theDecay, k).num_external, complex(0.0, 0.0));
  }

  // Keep only the contributions feeding the requested cross-section component
  // (0 = full, 1 = E1, 2 = E2); E1/E2 select electric channels of that L.
  auto includeInternal = [&](const azure2::InternalRow *mg) -> bool {
    if (xsComponent == 0) return true;
    const AChannel *ex = compound()->GetJGroup(mg->group)->GetChannel(mg->chp);
    return ex->GetRadType() == 'E' && ex->GetL() == xsComponent;
  };
  auto includeEC = [&](const azure2::ExternalRow *ecg) -> bool {
    if (xsComponent == 0) return true;
    return ecg->rad_type == 'E' && ecg->mult == xsComponent;
  };

  // The angle-integrated cotangent as a function of its own cotangent `bar`: the chain below calls it with
  // fitBar; a normalised Legendre coefficient a_k = N_k / sigma_int calls it with -fitBar a_k / sigma_int.
  // symmetryFactor: the (1+delta_12) of an identical entrance pair where the forward value carries it
  // (the reported cross section), 1 where it does not (sigma_int of the angular coefficients).
  auto integratedBar = [&](double bar, double symmetryFactor) {
    // ---- Angle-integrated: model = (1/100) sum_k sum_temp w_temp |T_temp|^2,
    //      w_temp = geom (2J+1) I1I2 (1+delta_12),  T_temp = sum_m T(k,m) sharing (J,l,l'). ----
    for (int k = 1; k <= nK; k++) {
      // Mirror the angle-integrated UPOS over-counting guard.
      if (aa != exitPairNum && compound()->GetPair(exitPairNum)->GetPType() == 0) {
        if (P.Group(theDecay, k).sp != P.Group(theDecay, k).sp2) continue;
      }
      int numM = P.Group(theDecay, k).num_internal;
      int numEC = P.Group(theDecay, k).num_external;
      struct Group {
        double j;
        int l;
        int lp;
        complex tsum;
      };
      std::vector<Group> groups;
      std::vector<int> mGroupIdx(numM + 1, -1);
      std::vector<int> ecGroupIdx(numEC + 1, -1);
      auto findGroup = [&](double jv, int lv, int lpv) -> int {
        for (int g = 0; g < (int)groups.size(); g++)
          if (groups[g].j == jv && groups[g].l == lv && groups[g].lp == lpv) return g;
        groups.push_back(Group{jv, lv, lpv, complex(0.0, 0.0)});
        return (int)groups.size() - 1;
      };
      // Internal pathways T(k,m), grouped by (J,l,l').
      for (int m = 1; m <= numM; m++) {
        const azure2::InternalRow *mg = &P.Internal(P.Group(theDecay, k), m);
        if (!includeInternal(mg)) continue;
        const JGroup *jg = compound()->GetJGroup(mg->group);
        int gi = findGroup(jg->GetJ(),
                           jg->GetChannel(mg->ch)->GetL(),
                           jg->GetChannel(mg->chp)->GetL());
        groups[gi].tsum += this->GetTMatrixElement(k, m);
        mGroupIdx[m] = gi;
      }
      // External-capture pathways T(k,ecm), grouped by (J, L, mult) and merged
      // with the matching internal group (forward: same temp T-matrix list).
      for (int m = 1; m <= numEC; m++) {
        const azure2::ExternalRow *ecg = &P.External(P.Group(theDecay, k), m);
        if (!includeEC(ecg)) continue;
        int gi = findGroup(ecg->J, ecg->l, ecg->mult);
        groups[gi].tsum += this->GetECTMatrixElement(k, m);
        ecGroupIdx[m] = gi;
      }
      std::vector<complex> tBarGroup(groups.size());
      for (int g = 0; g < (int)groups.size(); g++) {
        tBarGroup[g] = azure2::physics::partial_wave_cotangent(bar, geom, groups[g].j, i1i2, symmetryFactor, groups[g].tsum);
      }
      // accumulate (not assign): a normalised Legendre coefficient adds this to its order-k cotangents
      for (int m = 1; m <= numM; m++)
        if (mGroupIdx[m] >= 0) tBar[k - 1][m - 1] += tBarGroup[mGroupIdx[m]];
      for (int m = 1; m <= numEC; m++)
        if (ecGroupIdx[m] >= 0) ecBar[k - 1][m - 1] += tBarGroup[ecGroupIdx[m]];
    }
  };

  if (point->IsTotal()) {
    // ---- Total cross section by the optical theorem (GenMatrixFunc::CalculateBareCrossSection):
    //        model = (geom i1i2 / 100) sum (2J+1) 2 Re T(k,m) over the diagonal pathways,
    //      so d(model)/dRe T = (geom i1i2 / 100) 2 (2J+1), d(model)/dIm T = 0, and the cotangent (dL/dRe + i dL/dIm)
    //      is real. ----
    if (xsComponent != 0 || isDiff) UnreachableGate("sigma_total is angle-integrated and takes no multipole");
    for (int k = 1; k <= nK; k++) {
      if (P.Group(theDecay, k).CrossSpin()) continue;
      for (int m = 1; m <= P.Group(theDecay, k).num_internal; m++) {
        const azure2::InternalRow *mg = &P.Internal(P.Group(theDecay, k), m);
        if (mg->ch != mg->chp) continue;
        const double jValue = compound()->GetJGroup(mg->group)->GetJ();
        tBar[k - 1][m - 1] += azure2::physics::optical_theorem_cotangent(fitBar, geom, i1i2, jValue);
      }
    }
  } else if (isBetaDelayed) {
    // ---- Beta-delayed particle emission. GenMatrixFunc forms
    //        sum = sum_{k,m} 25 |T(k,m)|^2 ,   model = Re(sum)/100
    //      with no geometrical or spin factor and no external-capture term --
    //      an incoherent sum over pathways, with no interference between them.
    //      Hence d(model)/dT*(k,m) = (25/100) T(k,m), and the factor of two is
    //      the cotangent convention used throughout this function. ----
    for (int k = 1; k <= nK; k++)
      for (int m = 1; m <= P.Group(theDecay, k).num_internal; m++)
        tBar[k - 1][m - 1] += azure2::physics::beta_delayed_cotangent(fitBar, this->GetTMatrixElement(k, m));
  } else if (point->IsAnalyzingPower() && !ayCrossSection) {
    // ---- Vector analyzing power. A_y is built from the channel-spin amplitude
    //      matrix, which is linear in T(k,m) with the coefficients AddPathway
    //      applies, so the whole derivative is a contraction against those same
    //      coefficients. The Coulomb amplitude carries no parameter dependence
    //      and drops out. ----
    if (compound()->GetPair(aa)->GetPType() != 0) UnreachableGate("an A_y needs a particle entrance (Build refuses)");

    if (compound()->GetPair(exitPairNum)->GetPType() == 10) {
      // ---- Capture analyzing power (Seyler and Weller, Phys. Rev. C 20
      //      (1979) 453). The observable is a ratio of two forms that are each
      //      bilinear in the T-matrix elements,
      //        A_y = N/D,  N = sum bk Re(i T1 T2*) P_k^1,
      //                    D = sum ak Re(T1 T2*)   P_k,
      //      so with d Re(T1 T2*)/dT1* = T2/2 and d Re(i T1 T2*)/dT1* = -i T2/2
      //      (and their partners with 1 <-> 2, i -> -i),
      //        dA_y/dT1* = [ -i bk P_k^1 T2 - A_y ak P_k T2 ] / (2 D),
      //        dA_y/dT2* = [ +i bk P_k^1 T1 - A_y ak P_k T1 ] / (2 D),
      //      and the cotangent convention of this function supplies the 2. ----
      const int maxL = angles->GetMaxLOrder();
      if (!theDecay.has_capture_ay) return false;  // terms not requested: the forward pass refused the point
      const int nTerms = theDecay.num_capture_ay;
      if (nTerms == 0) return false;

      const double xCos = cos(angles->SlotAngle(rowSlot) * pi / 180.0);  // the reader's angle at the slot
      std::vector<double> assocP(maxL + 1, 0.0);
      for (int k = 1; k <= maxL; k++) assocP[k] = AngCoeff::LegendreP1(k, xCos);

      auto element = [&](int kg, int p, bool isEC) -> complex {
        return isEC ? this->GetECTMatrixElement(kg, p)
                    : this->GetTMatrixElement(kg, p);
      };
      double num = 0.0, den = 0.0;
      for (int i = 0; i < nTerms; i++) {
        const azure2::CaptureAyTerm *t = &P.capture_ay[theDecay.first_capture_ay + i];
        if (t->kOrder > maxL) continue;
        azure2::physics::add_capture_ay_sums(t->ak, t->bk, angles->RowP(rowSlot, t->kOrder), assocP[t->kOrder],
                                             element(t->kGroup1, t->path1, t->isEC1),
                                             element(t->kGroup2, t->path2, t->isEC2), den, num);
      }
      if (fabs(den) < 1.e-300) return false;
      const double ay = num / den;

      for (int i = 0; i < nTerms; i++) {
        const azure2::CaptureAyTerm *t = &P.capture_ay[theDecay.first_capture_ay + i];
        if (t->kOrder > maxL) continue;
        const complex T1 = element(t->kGroup1, t->path1, t->isEC1);
        const complex T2 = element(t->kGroup2, t->path2, t->isEC2);
        const double pk = angles->RowP(rowSlot, t->kOrder);
        const double pk1 = assocP[t->kOrder];
        complex d1, d2;
        azure2::physics::capture_ay_cotangents(t->ak, t->bk, pk, pk1, ay, den, T1, T2, d1, d2);
        if (t->isEC1)
          ecBar[t->kGroup1 - 1][t->path1 - 1] += fitBar * d1;
        else
          tBar[t->kGroup1 - 1][t->path1 - 1] += fitBar * d1;
        if (t->isEC2)
          ecBar[t->kGroup2 - 1][t->path2 - 1] += fitBar * d2;
        else
          tBar[t->kGroup2 - 1][t->path2 - 1] += fitBar * d2;
      }
    } else if (compound()->GetPair(exitPairNum)->GetPType() != 0) {
      return false;
    } else {
      Polarization::AmplitudeMatrix M(compound(), angles->SlotAngle(rowSlot), aa, exitPairNum);  // the reader's, at the slot
      for (int k = 1; k <= nK; k++) {
        if (P.Group(theDecay, k).CrossSpin()) continue;  // as the forward sum
        for (int m = 1; m <= P.Group(theDecay, k).num_internal; m++) {
          const azure2::InternalRow *mg = &P.Internal(P.Group(theDecay, k), m);
          M.AddPathway(mg->group, mg->ch, mg->chp,
                       this->GetTMatrixElement(k, m));
        }
      }
      complex fdAy(0.0, 0.0), fxAy(0.0, 0.0);  // the reader's Coulomb amplitudes at the slot's energy
      if (aa == exitPairNum) angles->CoulombAt(rowSlot, fdAy, fxAy), M.AddCoulomb(fdAy, fxAy);
      if (M.size() == 0) return false;

      const std::vector<complex> bar = M.AnalyzingPowerBar();
      for (int k = 1; k <= nK; k++) {
        if (P.Group(theDecay, k).CrossSpin()) continue;
        for (int m = 1; m <= P.Group(theDecay, k).num_internal; m++) {
          const azure2::InternalRow *mg = &P.Internal(P.Group(theDecay, k), m);
          if (!includeInternal(mg)) continue;
          tBar[k - 1][m - 1] += fitBar * M.PathwayAdjoint(mg->group, mg->ch, mg->chp, bar);
        }
      }
    }
  } else if (isPhase) {
    // ---- Phase shift: model = (90/pi) * arg(U) [+ const for identical pairs],
    //      U = sum_{matching k,m} (expCP^2 - T(k,m)) / expCP^2, matched by
    //      (J, l, elastic channel) = (segment J, segment L). ----
    double segJ = point->GetJ();
    int segL = point->GetL();
    struct PM {
      int k;
      int m;
      complex expCP;
    };
    std::vector<PM> matches;
    complex U(0.0, 0.0);
    for (int k = 1; k <= nK; k++) {
      for (int m = 1; m <= P.Group(theDecay, k).num_internal; m++) {
        const azure2::InternalRow *mg = &P.Internal(P.Group(theDecay, k), m);
        const JGroup *jg = compound()->GetJGroup(mg->group);
        const AChannel *en = jg->GetChannel(mg->ch);
        const AChannel *ex = jg->GetChannel(mg->chp);
        if (jg->GetJ() == segJ && en->GetL() == segL && en == ex) {
          complex expCP = point->GetExpCoulombPhase(mg->group, mg->ch);
          U += azure2::physics::collision_from_t_matrix(expCP, this->GetTMatrixElement(k, m));
          matches.push_back(PM{k, m, expCP});
        }
      }
    }
    double absU2 = std::norm(U);
    if (absU2 > 0.0) {
      // phase = (90/pi) arg(U) -> Ubar = fitBar (90/pi)(-Im U + i Re U)/|U|^2.
      const complex Ubar = azure2::physics::phase_shift_cotangent(fitBar, U);
      // U = sum [1 - T/expCP^2]  ->  dU/dT = -1/expCP^2.
      for (const PM &pm : matches)
        tBar[pm.k - 1][pm.m - 1] += azure2::physics::collision_cotangent_to_t_matrix(pm.expCP, Ubar);
    }
  } else if (!isDiff && (!angles->IsCoefficient() || a0Mode)) {
    // The reported angle-integrated cross section carries (1+delta_12) of an identical entrance pair
    // (GenMatrixFunc); the A_0 of an angular-coefficient point is the unscaled integral.
    integratedBar(a0Mode ? fitBar / (4.0 * pi) : fitBar,
                  a0Mode ? 1.0 : compound()->GetPair(aa)->GetEntranceSymmetryFactor());
  } else {
    // ---- Differential: model = (Re(CT)+Re(RT)+Re(IT))/100. ----
    //   RT = (geom I1I2 rtFactor/pi) Re( sum_{kL,inter} weight T1 conj(T2) P_L )
    //   IT = (geom itFactor/sqrt(pi)) Re( i sum_{k,m: en==ex} statSpin coulombAmp conj(T) P_l )  (elastic)
    //   CT = const (no T dependence).
    // The standard case (entrance==exit pair, or gamma exit) with external-capture interferences, and the
    // inelastic-particle UPOS / normal-angular-distribution case.
    const PPair *exitPairDiff = compound()->GetPair(exitPairNum);
    // Identical entrance pair: x4 / x2 for elastic, (1+delta_12) otherwise,
    // mirroring GenMatrixFunc::CalculateCrossSection.
    double rtFactor = 1.0, itFactor = 1.0;
    if (compound()->GetPair(aa)->IsIdentical()) {
      if (aa == exitPairNum) {
        rtFactor = 4.0;
        itFactor = 2.0;
      } else {
        rtFactor = compound()->GetPair(aa)->GetEntranceSymmetryFactor();
      }
    }
    double cR = azure2::physics::legendre_sum_factor(geom, i1i2, rtFactor);
    // An angular-coefficient point (Segment isDiff 9): the forward value is
    //   a_k = (4 geom i1i2 / (100 sigma_int)) sum_{L=k} z T1 T2*   (normalised)   or   A_k = a_k sigma_int / 4 pi,
    // i.e. the differential model's Legendre coefficient of order k, coeff_k = cR sum_{L=k} z T1 T2*, times
    // 4 pi / (rtFactor sigma_int) (normalised) or 1 / rtFactor (absolute).  So the differential machinery below
    // with a Kronecker delta in place of P_l(cos theta) gives the N_k part; the normalised one then adds the
    // sigma_int cotangent -fitBar a_k / sigma_int (integratedBar).  No Coulomb terms: elastic pairs are refused.
    const bool coef = angles->IsCoefficient();
    const int kOrd = coef ? angles->GetCoefficientOrder() : -1;
    const double sigmaInt = coef ? 4.0 * pi * point->GetAngDistA0(results) : 0.0;
    if (coef && point->IsUPOS()) UnreachableGate("UPOS coefficients: quantity legendre takes no secondary_gamma");
    if (coef && !(sigmaInt > 0.0)) return false;  // a vanishing integrated cross section (parameter values)
    const double dfac =
        coef ? azure2::physics::legendre_coefficient_factor(angles->GetCoefficientAbsolute(), rtFactor, sigmaInt) : 0.0;
    auto pl = [&](int l) -> double { return coef ? ((l == kOrd) ? dfac : 0.0) : angles->RowP(rowSlot, l); };

    // Fetch a T value (internal 'R' or external-capture 'E') and route a
    // cotangent contribution to the right accumulator.
    auto Tval = [&](int K, int idx, char rad) -> complex {
      return (rad == 'E') ? this->GetECTMatrixElement(K, idx)
                          : this->GetTMatrixElement(K, idx);
    };
    auto Tadd = [&](int K, int idx, char rad, complex v) {
      if (rad == 'E')
        ecBar[K - 1][idx - 1] += v;
      else
        tBar[K - 1][idx - 1] += v;
    };
    // model += W Re(T1 conj T2)  ->  T1bar += fitBar W T2 ;  T2bar += fitBar W T1.
    auto addQuad = [&](int K, int m1, char r1, int m2, char r2, double W) {
      complex T1 = Tval(K, m1, r1), T2 = Tval(K, m2, r2);
      Tadd(K, m1, r1, azure2::physics::interference_pair_cotangent(fitBar, W, T2));
      Tadd(K, m2, r2, azure2::physics::interference_pair_cotangent(fitBar, W, T1));
    };

    const bool inelasticParticle = (aa != exitPairNum && exitPairDiff->GetPType() == 0);
    // Away from UPOS the sum is linear in the Legendre polynomials with
    // angle-independent coefficients: the cotangent of each pathway is a sum
    // over orders of products taken once per T-matrix (TBarProducts).
    if (!(inelasticParticle && point->IsUPOS()) && this->TBarProducts(P, theDecay, inelasticParticle)) {
      const complex scale = fitBar * cR;
      const azure2::InterferencePathway *paths = P.interference_pathways.data() + theDecay.first_interference_pathway;
      const size_t numPaths = theDecay.num_interference_pathways;
      const int nl = theDecay.num_orders;
      tbar_legendre_.assign(nl, 0.0);
      const int nP = std::min(nl, angles->GetMaxLOrder() + 1);
      for (int l = 0; l < nP; l++) tbar_legendre_[l] = pl(l);
      for (size_t p = 0; p < numPaths; p++) {
        const complex s = azure2::physics::legendre_contraction(&tbar_q_[p * nl], tbar_legendre_.data(), nl);
        if (s != complex(0.0, 0.0)) Tadd(paths[p].k, paths[p].m, paths[p].rad, scale * s);
      }
    } else
    for (int inter = 0; inter < theDecay.num_interferences; inter++) {
      const azure2::FlatInterference &itf = P.interferences[theDecay.first_interference + inter];
      int K = itf.k;
      int lOrder = itf.l;
      double legP = pl(lOrder);
      const azure2::SpinGroupRow &kg = P.Group(theDecay, K);
      {
        const std::string type{itf.r1, itf.r2};  // "RR","ER","RE","EE": BuildPathwayTable writes 'R' or 'E' only
        char rad1 = type[0], rad2 = type[1];  // 'R' internal, 'E' external
        int m1 = itf.m1, m2 = itf.m2;

        if (inelasticParticle) {
          // Inelastic particle exit: pick the normal or UPOS sub-case.
          double sp1 = kg.sp, sp2 = kg.sp2;
          int lp1 = compound()->GetJGroup(P.Internal(kg, m1).group)->GetChannel(P.Internal(kg, m1).chp)->GetL();
          int lp2 = compound()->GetJGroup(P.Internal(kg, m2).group)->GetChannel(P.Internal(kg, m2).chp)->GetL();
          if (sp1 == sp2 && !point->IsUPOS()) {
            addQuad(K, m1, rad1, m2, rad2, azure2::physics::interference_weight(cR, itf.z, legP));
          }
          if (lp1 == lp2 && point->IsUPOS() && type == "RR") {
            // R_L coefficient (constant w.r.t. parameters), even lOrder only.
            double finalL = (double)point->GetSecondaryDecayL();
            double Ic = point->GetIc();
            double j2f = compound()->GetPair(exitPairNum)->GetJ(2);  // the exit pair, not the decay index (GenMatrixFunc)
            double delta = point->GetDelta();
            double R_L = 0.0;
            if ((int)lOrder % 2 == 0) {
              R_L = this->GetRk(j2f, finalL, finalL, Ic, lOrder);
              if (delta != 0.0) {
                double finalLp = finalL + 1.0;
                double R_LLp = this->GetRk(j2f, finalL, finalLp, Ic, lOrder);
                double R_LpLp = this->GetRk(j2f, finalLp, finalLp, Ic, lOrder);
                // unconditionally, as in GenMatrixFunc (checklist t3-24)
                R_L = azure2::physics::upos_mixed_rk(R_L, R_LLp, R_LpLp, delta);
              }
            }
            double W = azure2::physics::upos_interference_weight(cR, itf.z_upos, lOrder, R_L, legP);
            addQuad(K, m1, 'R', m2, 'R', W);
          }
        } else {
          // Standard case (entrance==exit pair, or gamma exit).
          addQuad(K, m1, rad1, m2, rad2, azure2::physics::interference_weight(cR, itf.z, legP));
        }
      }
    }

    if (coef && !angles->GetCoefficientAbsolute()) {
      // a_k = N_k / sigma_int: the sigma_int cotangent through the angle-integrated form
      integratedBar(-fitBar * point->ValueFor(*angles, rowSlot, results) / sigmaInt, 1.0);
    }

    // f_C of channel spin s, the reader's at its slot (Site::CoulombAt): f_d + (-1)^s f_x (f_x = 0 but with spin exchange)
    complex fd(0.0, 0.0), fx(0.0, 0.0);
    if (!coef && aa == exitPairNum) angles->CoulombAt(rowSlot, fd, fx);
    auto readerCoulomb = [&](double s) { return azure2::physics::channel_spin_coulomb_amplitude(s, fd, fx); };
    if (!coef && aa == exitPairNum && theDecay.has_interferences) {  // elastic Coulomb-nuclear interference
      const complex cI = azure2::physics::coulomb_interference_factor(geom, itFactor);
      for (int n = 0; n < theDecay.num_diagonal; n++) {
        const azure2::DiagonalPathway &d = P.diagonal[theDecay.first_diagonal + n];
        tBar[d.k - 1][d.m - 1] += azure2::physics::coulomb_interference_cotangent(fitBar, cI, d.statSpin, readerCoulomb(d.s),
                                                                                 angles->RowP(rowSlot, d.l));
      }
    } else if (!coef && aa == exitPairNum) {  // elastic Coulomb-nuclear interference
      complex cI = azure2::physics::coulomb_interference_factor(geom, itFactor);
      for (int k = 1; k <= nK; k++) {
        for (int m = 1; m <= P.Group(theDecay, k).num_internal; m++) {
          const azure2::InternalRow *mg = &P.Internal(P.Group(theDecay, k), m);
          const AChannel *en = compound()->GetJGroup(mg->group)->GetChannel(mg->ch);
          const AChannel *ex = compound()->GetJGroup(mg->group)->GetChannel(mg->chp);
          if (en == ex) {
            int l = en->GetL();
            // model += Re(B conj T)  ->  Tbar += fitBar B.
            tBar[k - 1][m - 1] += azure2::physics::coulomb_interference_cotangent(
                fitBar, cI, mg->stat_spin_factor, readerCoulomb(en->GetS()), angles->RowP(rowSlot, l));
          }
        }
      }
    }
  }

  return true;
}

/*!
 * The cotangent coefficients of the differential sum of `theDecay` for the
 * T-matrix held now, by pathway and Legendre order: model += cR Re(sum_i
 * Z1Z2_i T1_i conj(T2_i) P_{L_i}) gives T1bar += cR Z1Z2 P_L T2 and
 * T2bar += cR Z1Z2 P_L T1, so each pathway's cotangent is cR sum_L q_L P_L with
 * q_L summed over the interferences once.  `particleExit` keeps only the
 * k-groups of equal channel spins, as the forward sum.  False (nothing built)
 * for an interference type the adjoint does not know.
 */
bool AMatrixFunc::TBarProducts(const azure2::EntrancePathways &P, const azure2::ExitRow &theDecay, bool particleExit) {
  if (tbar_decay_ == &theDecay && tbar_particle_ == particleExit) return true;
  if (!theDecay.has_interferences) return false;
  const azure2::FlatInterference *flat = P.interferences.data() + theDecay.first_interference;
  const int nl = theDecay.num_orders;
  const size_t np = theDecay.num_interference_pathways;
  tbar_q_.assign(np * nl, complex(0.0, 0.0));
  const std::vector<vector_c> *T = tmatrix_.empty() ? nullptr : &tmatrix_[0];
  for (const azure2::FlatInterference *fp = flat; fp != flat + theDecay.num_interferences; fp++) {
    const azure2::FlatInterference &f = *fp;
    if (particleExit && !f.spEqual) continue;
    if (!f.r1 || !f.r2) {  // an interference type the adjoint does not know
      tbar_decay_ = nullptr;
      return false;
    }
    const complex T1 = (f.r1 == 'E') ? ec_tmatrix_[f.k - 1][f.m1 - 1] : (T ? (*T)[f.k - 1][f.m1 - 1] : complex(0.0, 0.0));
    const complex T2 = (f.r2 == 'E') ? ec_tmatrix_[f.k - 1][f.m2 - 1] : (T ? (*T)[f.k - 1][f.m2 - 1] : complex(0.0, 0.0));
    azure2::physics::add_interference_cotangents(f.z, T1, T2, tbar_q_[(size_t)f.s1 * nl + f.l],
                                                 tbar_q_[(size_t)f.s2 * nl + f.l]);
  }
  tbar_decay_ = &theDecay;
  tbar_particle_ = particleExit;
  return true;
}
