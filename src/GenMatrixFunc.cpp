#include "AngCoeff.h"
#include "Profile.h"
#include "AChannel.h"
#include "CompoundSystem.h"
#include "Unit.h"
#include "GenMatrixFunc.h"
#include "PolarizationFunc.h"
#include <assert.h>
#include <iostream>
#include <stdexcept>

/*!
 * The child classes AMatrixFunc or RMatrixFunc contain functions
 * to calculate the T-Matrix from the fitted R-Matrix parameters.  This function
 * then calculates the cross section from the T-Matrix elements.
 */

/*!
 * Electron screening multiplies the bare cross section by f = exp(pi eta U_e / E) (Assenbaum et al.
 * 1987).  It acts on cross sections only: a phase shift, the angular-distribution coefficients and an
 * analyzing power are left as they are, while an analyzing-power sub-point keeps its (screened) cross
 * section as the weight of the target average.  The E1 and E2 components are scaled with the total.
 */
void GenMatrixFunc::CalculateCrossSection(Site *point, azure2::Results &results) {
  AZUREProfile::Timer profileTimer(configure().profile, AZUREProfile::kCrossSection);
  this->CalculateBareCrossSection(point, results);
  if (point->GetScreeningParam() == 0) return;
  if (point->IsAngularDist() || point->IsPhase()) return;
  if (point->IsAnalyzingPower() && !point->IsSubPoint()) return;
  const double f = point->ScreeningFactor();
  if (f == 1.0) return;
  if (point->AngleFree()) point->ScaleObservable(results, f);
  point->SetFitCrossSection(results, f * point->GetFitCrossSection(results));
  point->SetFitE1CrossSection(results, f * point->GetFitE1CrossSection(results));
  point->SetFitE2CrossSection(results, f * point->GetFitE2CrossSection(results));
}

void GenMatrixFunc::CalculateBareCrossSection(Site *point, azure2::Results &results) {
  complex sum(0., 0.);
  complex sumE1(0., 0.);
  complex sumE2(0., 0.);
  int aa = compound()->GetPairNumFromKey(point->GetEntranceKey());
  // `ir` below is the DECAY INDEX within the entrance pair's list; whether the
  // reaction is elastic is a question about pair numbers, so compare those.
  // (aa == ir happened to hold for the first pair of a file, which is why a
  // second compound nucleus's elastic segments came out as inelastic.)
  const int exitPairNum = compound()->GetPairNumFromKey(point->GetExitKey());
  const azure2::EntrancePathways &P = compound()->Pathways().Entrance(aa);
  if (compound()->GetPair(compound()->GetPairNumFromKey(point->GetExitKey()))->GetPType() == 10 &&
      (configure().paramMask & Config::USE_RMC_FORMALISM)) {
    const int decayNum = P.ExitOf(aa);
    if (!decayNum) throw std::runtime_error(NoPathwayMessage(aa, aa));
    const azure2::ExitRow &elastic = P.exits[decayNum - 1];
    for (int k = 1; k <= elastic.num_groups; k++) {
      const azure2::SpinGroupRow &g = P.Group(elastic, k);
      for (int m = 1; m <= g.num_internal; m++) {
        const azure2::InternalRow &theMGroup = P.Internal(g, m);
        if (theMGroup.ch == theMGroup.chp) {
          double jValue = compound()->GetJGroup(theMGroup.group)->GetJ();
          sum += 2. * point->GetGeometricalFactor() *
              (2. * jValue + 1.) * compound()->GetPair(aa)->GetI1I2Factor() *
              imag(this->GetTMatrixElement(k, m, decayNum));
        }
      }
    }
    for (int dp = 1; dp <= (int)P.exits.size(); dp++) {
      if (compound()->GetPair(P.exits[dp - 1].pair)->GetPType() == 0) {
        for (int k = 1; k <= P.exits[dp - 1].num_groups; k++) {
          const azure2::SpinGroupRow &g = P.Group(P.exits[dp - 1], k);
          // A cross-spin UPOS group repeats pathways of the equal-spin groups.
          if (g.CrossSpin()) continue;
          this->ClearTempTMatrices();
          for (int m = 1; m <= g.num_internal; m++) {
            const azure2::InternalRow &theMGroup = P.Internal(g, m);
            int lValue = compound()->GetJGroup(theMGroup.group)->GetChannel(theMGroup.ch)->GetL();
            int lpValue = compound()->GetJGroup(theMGroup.group)->GetChannel(theMGroup.chp)->GetL();
            double jValue = compound()->GetJGroup(theMGroup.group)->GetJ();
            int tempTNum = this->IsTempTMatrix(jValue, lValue, lpValue);
            if (!tempTNum) {
              TempTMatrix temptmatrix = {jValue, lValue, lpValue, this->GetTMatrixElement(k, m, dp)};
              this->NewTempTMatrix(temptmatrix);
            } else
              this->AddToTempTMatrix(tempTNum, this->GetTMatrixElement(k, m, dp));
          }
          for (int temp = 1; temp <= this->NumTempTMatrices(); temp++) {
            sum -= point->GetGeometricalFactor() *
                (2. * this->GetTempTMatrix(temp)->jValue + 1.) *
                compound()->GetPair(aa)->GetI1I2Factor() *
                (this->GetTempTMatrix(temp)->TMatrix) * conj(this->GetTempTMatrix(temp)->TMatrix);
          }
        }
      }
    }
    // (1+delta_12) of an identical entrance pair: see the angle-integrated
    // branch below.
    point->SetFitCrossSection(results, real(sum) / 100. * compound()->GetPair(aa)->GetEntranceSymmetryFactor());
    // This formalism has no angular dependence: a site holding Legendre coefficients holds the value as B_0.
    if (point->AngleFree()) point->SetIsotropicObservable(results, point->GetFitCrossSection(results));
  } else {
    const int ir = P.ExitOf(exitPairNum);
    if (!ir) throw std::runtime_error(NoPathwayMessage(aa, exitPairNum));
    const azure2::ExitRow &theDecay = P.exits[ir - 1];
    if (point->IsTotal() && !point->IsPhase()) {
      // The total cross section of the entrance pair by the optical theorem: with every channel of the J-group in its
      // level matrix, U is unitary, so sum_c' |T_cc'|^2 over all exit channels is 2 Re T_cc (T = 1 - U for a neutral
      // pair), and sigma_tot = (pi/k^2) g sum_J (2J+1) sum_c 2 Re T_cc over the entrance channels -- one T-matrix (the
      // elastic one) instead of one per exit pair.  Photon channels carry no width in the level matrix and external
      // capture is outside it, so capture is not included (add it as a SUM component).
      double total = 0.;
      for (int k = 1; k <= theDecay.num_groups; k++) {
        const azure2::SpinGroupRow &g = P.Group(theDecay, k);
        if (g.CrossSpin()) continue;
        for (int m = 1; m <= g.num_internal; m++) {
          const azure2::InternalRow &theMGroup = P.Internal(g, m);
          if (theMGroup.ch != theMGroup.chp) continue;
          const double jValue = compound()->GetJGroup(theMGroup.group)->GetJ();
          total += (2. * jValue + 1.) * 2. * real(this->GetTMatrixElement(k, m));
        }
      }
      point->SetFitCrossSection(results, point->GetGeometricalFactor() * compound()->GetPair(aa)->GetI1I2Factor() * total / 100.);
      point->SetFitE1CrossSection(results, 0.);
      point->SetFitE2CrossSection(results, 0.);
      return;
    }
    if (!point->IsPhase()) {
      double angleIntegratedXS = 0.;
      double angleIntegratedE1XS = 0.;
      double angleIntegratedE2XS = 0.;
      if (!point->IsDifferential()) {
        for (int k = 1; k <= theDecay.num_groups; k++) {
          const azure2::SpinGroupRow &g = P.Group(theDecay, k);
          // For inelastic particle reactions (aa!=ir, GetPType==0) with UPOS 3-param KGroups,
          // skip KGroups where sp!=sp2 to avoid over-counting in angle-integrated XS.
          if (aa != exitPairNum && compound()->GetPair(compound()->GetPairNumFromKey(point->GetExitKey()))->GetPType() == 0) {
            if (g.sp != g.sp2) continue;
          }
          this->ClearTempTMatrices();
          this->ClearTempTMatricesE1();
          this->ClearTempTMatricesE2();
          for (int m = 1; m <= g.num_internal; m++) {
            const azure2::InternalRow &theMGroup = P.Internal(g, m);
            if (compound()->GetPair(aa)->GetPType() == 20) {
              sum += 25. * this->GetTMatrixElement(k, m) * conj(this->GetTMatrixElement(k, m));
            } else {
              int lValue = compound()->GetJGroup(theMGroup.group)->GetChannel(theMGroup.ch)->GetL();
              int lpValue = compound()->GetJGroup(theMGroup.group)->GetChannel(theMGroup.chp)->GetL();
              double jValue = compound()->GetJGroup(theMGroup.group)->GetJ();
              int tempTNum = this->IsTempTMatrix(jValue, lValue, lpValue);
              if (!tempTNum) {
                TempTMatrix temptmatrix = {jValue, lValue, lpValue, this->GetTMatrixElement(k, m)};
                this->NewTempTMatrix(temptmatrix);
              } else
                this->AddToTempTMatrix(tempTNum, this->GetTMatrixElement(k, m));
              if (compound()->GetJGroup(theMGroup.group)->GetChannel(theMGroup.chp)->GetRadType() == 'E' &&
                  compound()->GetJGroup(theMGroup.group)->GetChannel(theMGroup.chp)->GetL() == 1) {
                int tempTNumE1 = this->IsTempTMatrixE1(jValue, lValue, lpValue);
                if (!tempTNumE1) {
                  TempTMatrix temptmatrixE1 = {jValue, lValue, lpValue, this->GetTMatrixElement(k, m)};
                  this->NewTempTMatrixE1(temptmatrixE1);
                } else
                  this->AddToTempTMatrixE1(tempTNumE1, this->GetTMatrixElement(k, m));
              }
              if (compound()->GetJGroup(theMGroup.group)->GetChannel(theMGroup.chp)->GetRadType() == 'E' &&
                  compound()->GetJGroup(theMGroup.group)->GetChannel(theMGroup.chp)->GetL() == 2) {
                int tempTNumE2 = this->IsTempTMatrixE2(jValue, lValue, lpValue);
                if (!tempTNumE2) {
                  TempTMatrix temptmatrixE2 = {jValue, lValue, lpValue, this->GetTMatrixElement(k, m)};
                  this->NewTempTMatrixE2(temptmatrixE2);
                } else
                  this->AddToTempTMatrixE2(tempTNumE2, this->GetTMatrixElement(k, m));
              }
            }
          }
          if (compound()->GetPair(aa)->GetPType() == 20) continue;
          for (int m = 1; m <= g.num_external; m++) {
            const azure2::ExternalRow &theECMGroup = P.External(g, m);
            int lValue = theECMGroup.l;
            int lpValue = theECMGroup.mult;
            double jValue = theECMGroup.J;
            int tempTNum = this->IsTempTMatrix(jValue, lValue, lpValue);
            if (!tempTNum) {
              TempTMatrix temptmatrix = {jValue, lValue, lpValue, this->GetECTMatrixElement(k, m)};
              this->NewTempTMatrix(temptmatrix);
            } else
              this->AddToTempTMatrix(tempTNum, this->GetECTMatrixElement(k, m));
            if (theECMGroup.rad_type == 'E' && theECMGroup.mult == 1) {
              int tempTNumE1 = this->IsTempTMatrixE1(jValue, lValue, lpValue);
              if (!tempTNumE1) {
                TempTMatrix temptmatrixE1 = {jValue, lValue, lpValue, this->GetECTMatrixElement(k, m)};
                this->NewTempTMatrixE1(temptmatrixE1);
              } else
                this->AddToTempTMatrixE1(tempTNumE1, this->GetECTMatrixElement(k, m));
            }
            if (theECMGroup.rad_type == 'E' && theECMGroup.mult == 2) {
              int tempTNumE2 = this->IsTempTMatrixE2(jValue, lValue, lpValue);
              if (!tempTNumE2) {
                TempTMatrix temptmatrixE2 = {jValue, lValue, lpValue, this->GetECTMatrixElement(k, m)};
                this->NewTempTMatrixE2(temptmatrixE2);
              } else
                this->AddToTempTMatrixE2(tempTNumE2, this->GetECTMatrixElement(k, m));
            }
          }
          for (int temp = 1; temp <= this->NumTempTMatrices(); temp++) {
            sum += point->GetGeometricalFactor() *
                (2. * this->GetTempTMatrix(temp)->jValue + 1.) *
                compound()->GetPair(aa)->GetI1I2Factor() *
                (this->GetTempTMatrix(temp)->TMatrix) * conj(this->GetTempTMatrix(temp)->TMatrix);
          }
          for (int temp = 1; temp <= this->NumTempTMatricesE1(); temp++) {
            sumE1 += point->GetGeometricalFactor() *
                (2. * this->GetTempTMatrixE1(temp)->jValue + 1.) *
                compound()->GetPair(aa)->GetI1I2Factor() *
                (this->GetTempTMatrixE1(temp)->TMatrix) * conj(this->GetTempTMatrixE1(temp)->TMatrix);
          }
          for (int temp = 1; temp <= this->NumTempTMatricesE2(); temp++) {
            sumE2 += point->GetGeometricalFactor() *
                (2. * this->GetTempTMatrixE2(temp)->jValue + 1.) *
                compound()->GetPair(aa)->GetI1I2Factor() *
                (this->GetTempTMatrixE2(temp)->TMatrix) * conj(this->GetTempTMatrixE2(temp)->TMatrix);
          }
        }
        angleIntegratedXS = real(sum) / 100.;
        angleIntegratedE1XS = real(sumE1) / 100.;
        angleIntegratedE2XS = real(sumE2) / 100.;
        if (!point->IsAngularDist()) {
          // Identical particles in the ENTRANCE channel: the reported total
          // is sigma = (pi/k^2) sum_J g_J (1+delta_12) sum |T|^2 over the
          // channels symmetry allows (l+s even; CompoundSystem refuses or warns about
          // the others) -- for spin-0 bosons the familiar
          // (pi/k^2) sum_l (2l+1)[1+(-1)^l] T_l.  The sum above is the
          // distinguishable-particle one, so an identical entrance pair takes
          // x2 whatever the exit: elastic (4He+4He), reaction (12C+12C ->
          // 20Ne+a, d+d -> t+p) and capture (d+d -> 4He+g) alike.  It is the
          // convention of the reciprocity theorem with a symmetric U,
          //   w_12 k_12^2 sigma(12->34)/(1+delta_12)
          //     = w_34 k_34^2 sigma(34->12)/(1+delta_34),  w = (2i_1+1)(2i_2+1),
          // and counts each reaction once.  An identical EXIT pair (7Li+p ->
          // a+a) takes no factor: sum |U_cc'|^2 is already the probability of
          // the reaction, so the total counts events, not alphas.  For elastic
          // scattering x2 is half the full 0-pi integral of the symmetrized
          // |F(theta)|^2 (x4), which counts both particles of every event.
          // Checked against JINA ReacLib's triple-alpha rate (elastic, x2),
          // and in tests/identical_entrance_reaction against the one-level
          // (1+delta) Breit-Wigner for 12C+12C -> a+20Ne and d+d -> p+t.
          //
          // Scoped to the reported cross section only -- angleIntegratedXS
          // itself (and E1/E2) is left unscaled below, since it is also the
          // normalization denominator for angular-distribution coefficients
          // (angularCoeff), which are ratios the factor cancels from.
          double reportedXS = angleIntegratedXS;
          double reportedE1XS = angleIntegratedE1XS;
          double reportedE2XS = angleIntegratedE2XS;
          const double symmetryFactor = compound()->GetPair(aa)->GetEntranceSymmetryFactor();
          reportedXS *= symmetryFactor;
          reportedE1XS *= symmetryFactor;
          reportedE2XS *= symmetryFactor;
          point->SetFitCrossSection(results, reportedXS);
          point->SetFitE1CrossSection(results, reportedE1XS);
          point->SetFitE2CrossSection(results, reportedE2XS);
          return;
        }
      }
      // The observables of the reaction at this energy (SiteObservable.h).  A site of an angle-free kind holds them
      // in its record and its readers apply their own rows (Site::ValueFor), a site that is a row of its data (not a
      // slot) included.
      const SiteObservable kind = point->Observable();
      if (kind == SiteObservable::kLegendre) {
        point->SetLegendreCoefficients(results, this->LegendreCoefficientsOf(point, P, theDecay, aa, exitPairNum));
        if (point->IsAngularDist()) point->SetAngDistA0(results, angleIntegratedXS / (4.0 * pi));
        if (!point->IsSubPoint() && (point->IsDifferential() || point->IsCoefficient()))
          point->SetFitCrossSection(results, point->ValueFor(*point, 0, results));
        return;
      }
      if (kind == SiteObservable::kElastic) {
        point->SetElasticCoefficients(results, this->ElasticCoefficientsOf(point, P, theDecay, aa));
        if (!point->IsSubPoint()) point->SetFitCrossSection(results, point->ValueFor(*point, 0, results));
        return;
      }
      if (kind == SiteObservable::kAnalyzingPower) {
        // An analyzing-power segment reports A_y in place of the cross section; a slot keeps the coefficients, and
        // each reader contracts the cross section that weights its target average and A_y (Site::SlotValue,
        // SlotAnalyzingPower).
        point->SetAnalyzingPowerCoefficients(results, this->AnalyzingPowerCoefficientsOf(point, P, theDecay, aa,
                                                                                           exitPairNum));
        if (!point->IsSubPoint()) {
          const double ay = point->AnalyzingPowerFor(*point, 0, results);
          point->SetAnalyzingPower(results, ay);
          point->SetFitCrossSection(results, ay);
        }
        return;
      }
      // Every differential kind has an angle-free kind (ObservableOf); THM is THMMatrixFunc's.
      throw std::logic_error("GenMatrixFunc::CalculateBareCrossSection: a differential site of kind " +
                             std::string(kSiteObservables[(int)kind].name) + " (data set " +
                             std::to_string(point->GetSegmentKey()) + ")");
    } else if (aa == exitPairNum) {
      double segmentJ = point->GetJ();
      int segmentL = point->GetL();
      // The channel spin the phase belongs to ("spin <s>" on the segment line); without it the
      // segment must have a single spin for its (J, l), which WorkingCopy::Initialize enforces: summing
      // the U-matrix elements of two channel spins (as happened up to 2026-10-27) is not a phase.
      const double segmentS = point->GetPhaseS();
      const PPair *entrancePair = compound()->GetPair(aa);
      this->ClearTempTMatrices();
      for (int k = 1; k <= theDecay.num_groups; k++) {
        const azure2::SpinGroupRow &g = P.Group(theDecay, k);
        for (int m = 1; m <= g.num_internal; m++) {
          const azure2::InternalRow &theMGroup = P.Internal(g, m);
          double jValue = compound()->GetJGroup(theMGroup.group)->GetJ();
          const AChannel *entranceChannel = compound()->GetJGroup(theMGroup.group)->GetChannel(theMGroup.ch);
          int lValue = entranceChannel->GetL();
          const AChannel *exitChannel = compound()->GetJGroup(theMGroup.group)->GetChannel(theMGroup.chp);
          if (jValue == segmentJ && lValue == segmentL && entranceChannel == exitChannel &&
              (segmentS < 0.0 || std::fabs(entranceChannel->GetS() - segmentS) < 1.0e-6)) {
            complex expCoulPhaseSquared = point->GetExpCoulombPhase(theMGroup.group, theMGroup.ch) *
                point->GetExpCoulombPhase(theMGroup.group, theMGroup.ch);
            complex theUMatrix = (expCoulPhaseSquared - this->GetTMatrixElement(k, m)) / expCoulPhaseSquared;
            int tempTNum = this->IsTempTMatrix(jValue, lValue, lValue);
            if (!tempTNum) {
              TempTMatrix temptmatrix = {jValue, lValue, lValue, theUMatrix};
              this->NewTempTMatrix(temptmatrix);
            } else
              this->AddToTempTMatrix(tempTNum, theUMatrix);
          }
        }
      }
      assert(this->NumTempTMatrices() <= 1);
      // The U-matrix accumulated above is S_L = exp(2 i delta_L); the nuclear
      // phase shift delta_L is half its argument. For identical particles the
      // physical partial-wave amplitude is symmetrised by [1 + eps (-1)^(L+S)],
      // which doubles allowed waves and removes forbidden ones, but delta_L
      // per allowed wave has the same meaning as in the distinguishable case.
      double phase = 0.0;
      if (this->NumTempTMatrices() == 1) {
        phase = 180.0 / pi / 2.0 *
            atan2(imag(this->GetTempTMatrix(1)->TMatrix),
                  real(this->GetTempTMatrix(1)->TMatrix));
        // atan2 returns (-pi, pi], so phase lives in (-90, 90]. For
        // identical-particle pairs we wrap to [0, 180), the standard
        // convention in partial-wave phase-shift analyses (e^{2 i delta}
        // is pi-periodic, and delta_L grows monotonically through a
        // resonance in this range). Left untouched for distinguishable
        // pairs to preserve existing fit conventions.
        if (entrancePair->IsIdentical() && phase < 0.0) phase += 180.0;
      }
      point->SetFitCrossSection(results, phase);
    }
  }
}

void GenMatrixFunc::BuildOrderProducts(const azure2::EntrancePathways &P, const azure2::ExitRow &theDecay,
                                       bool particleExit, bool elastic) {
  if (rt_decay_ == &theDecay) return;
  rt_wl_.clear();
  // An exit without interference lists (a beta entrance) has an empty list here.
  const azure2::FlatInterference *flat = P.interferences.data() + theDecay.first_interference;
  const std::vector<vector_c> *T = tmatrix_.empty() ? nullptr : &tmatrix_[0];
  for (const azure2::FlatInterference *fp = flat; fp != flat + theDecay.num_interferences; fp++) {
    const azure2::FlatInterference &f = *fp;
    if (particleExit && !f.spEqual) continue;
    if ((int)rt_wl_.size() <= f.l) rt_wl_.resize(f.l + 1, complex(0.0, 0.0));
    if (!f.r1) continue;  // unknown type: no contribution
    const complex T1 = (f.r1 == 'E') ? ec_tmatrix_[f.k - 1][f.m1 - 1] : (T ? (*T)[f.k - 1][f.m1 - 1] : complex(0.0, 0.0));
    const complex T2 = (f.r2 == 'E') ? ec_tmatrix_[f.k - 1][f.m2 - 1] : (T ? (*T)[f.k - 1][f.m2 - 1] : complex(0.0, 0.0));
    rt_wl_[f.l] += f.z * T1 * conj(T2);
  }
  this->BuildCoulombInterference(P, theDecay, elastic);
  rt_decay_ = &theDecay;
}

/*!
 * dsigma/dOmega_c.m. = sum_L B_L P_L(cos theta) of a reaction (Blatt and Biedenharn, Rev. Mod. Phys. 24 (1952) 258,
 * in the coupling of the pathway table; for elastic scattering the nuclear amplitude's part): B_L = (geom I1I2 s /
 * (100 pi)) sum Re(Z T1 T2*) over the interferences of order L (b/sr; geom = pi/k^2 in fm^2, 1/100 to b; s the
 * identical-particle factor, LegendreSymmetryFactor).  A particle exit
 * keeps the k-groups of equal channel spins; an unobserved primary with an observed secondary (UPOS) weighs the
 * internal pairs of equal exit l by Z_upos sqrt(2L+1)/4 R_L of the secondary decay (Brune and deBoer, PRC 102, 024628
 * (2020), Eq. (7); R_L = 0 for odd L).
 */
std::vector<double> GenMatrixFunc::LegendreCoefficientsOf(const Site *point, const azure2::EntrancePathways &P,
                                                          const azure2::ExitRow &theDecay, int aa, int exitPairNum) {
  const bool particleExit = aa != exitPairNum && compound()->GetPair(exitPairNum)->GetPType() == 0;
  std::vector<double> c;
  if (particleExit && point->IsUPOS()) {
    c.assign(theDecay.num_orders, 0.0);
    const double finalL = (double)point->GetSecondaryDecayL();
    const double Ic = point->GetIc();
    // the residual nucleus of the exit pair (tests/upos_pair_index)
    const double j2f = compound()->GetPair(exitPairNum)->GetJ(2);
    const double delta = point->GetDelta();
    std::vector<double> rL(theDecay.num_orders, 0.0);
    for (int lOrder = 0; lOrder < theDecay.num_orders; lOrder += 2) {
      double R_L = this->GetRk(j2f, finalL, finalL, Ic, lOrder);
      if (delta != 0.) {
        // Eq. (7) unconditionally: every R_k is 0 where its coupling is not allowed (checklist t3-24).
        const double finalLp = finalL + 1.;
        const double R_LLp = this->GetRk(j2f, finalL, finalLp, Ic, lOrder);
        const double R_LpLp = this->GetRk(j2f, finalLp, finalLp, Ic, lOrder);
        R_L = (R_L + 2. * delta * R_LLp + delta * delta * R_LpLp) / (1. + delta * delta);
      }
      rL[lOrder] = R_L;
    }
    for (int inter = 0; inter < theDecay.num_interferences; inter++) {
      const azure2::FlatInterference &f = P.interferences[theDecay.first_interference + inter];
      if (f.r1 != 'R' || f.r2 != 'R') continue;  // only internal pathways interfere here
      const azure2::SpinGroupRow &g = P.Group(theDecay, f.k);
      const azure2::InternalRow &a = P.Internal(g, f.m1), &b = P.Internal(g, f.m2);
      if (compound()->GetJGroup(a.group)->GetChannel(a.chp)->GetL() != compound()->GetJGroup(b.group)->GetChannel(b.chp)->GetL())
        continue;
      const complex T1 = this->GetTMatrixElement(f.k, f.m1), T2 = this->GetTMatrixElement(f.k, f.m2);
      c[f.l] += real(f.z_upos * T1 * conj(T2)) * pow(2. * f.l + 1., 0.5) / (4.) * rL[f.l];
    }
  } else {
    this->BuildOrderProducts(P, theDecay, particleExit, aa == exitPairNum);
    c.resize(rt_wl_.size());
    for (std::size_t l = 0; l < rt_wl_.size(); l++) c[l] = real(rt_wl_[l]);
  }
  const double geom = point->GetGeometricalFactor(), i1i2 = compound()->GetPair(aa)->GetI1I2Factor();
  const double symmetry = LegendreSymmetryFactor(*compound()->GetPair(aa), aa == exitPairNum);
  for (double &b : c) b = b / pi * geom * i1i2 * symmetry / 100.;
  return c;
}

/*!
 * The observables of elastic scattering at the site's energy (SiteObservable::kElastic):
 *   dsigma/dOmega = |F_C|^2 + |F_N|^2 + 2 Re[F_C* F_N]  (geom = pi/k^2 times the spin average),
 * B_L^NN of |F_N|^2 from LegendreCoefficientsOf, and the interference
 *   2 Re[F_C* F_N] / 100 = Re sum_l (f_d D_l + f_x X_l) P_l,  D_l = (i geom f_it / (100 sqrt(pi))) sum_{s} c_{s,l},
 *   X_l = (i geom f_it / (100 sqrt(pi))) sum_s (-1)^s c_{s,l},
 * c_{s,l} = sum statSpin conj(T) over the diagonal pathways of channel spin s and orbital l (BuildCoulombInterference).
 * Identical particles: in channel spin s, F_X = f_X(theta) + (-1)^s f_X(pi - theta) (exchange multiplies |s v> by
 * (-1)^(2j-s), the state must pick up (-1)^(2j), and Y_l(pi - theta, pi) = (-1)^l Y_l(theta, 0)); the nuclear term
 * carries 1 + (-1)^(l'+s'), 2 on every channel CompoundSystem admits (l + s even), so |F_N|^2 = 4 |f_N|^2 (the factor 4 of
 * LegendreSymmetryFactor) and 2 Re[F_C* F_N] = 4 Re[F_C,s* f_N] (f_it = 2; 1 otherwise).  Only the Coulomb amplitude
 * depends on s: for j = 0 there is s = 0 alone and f_d is the symmetrized F_C (f_x = 0, X_l = 0); for j != 0,
 * F_C,s = f_d + (-1)^s f_x, and the spin average of |F_C,s|^2, sum_s (2s+1) |F_C,s|^2 / (2j+1)^2, is
 * |f_d|^2 + |f_x|^2 + 2 (-1)^(2j) / (2j+1) Re(f_d f_x*): c_0 = geom / 100, c_x = 2 c_0 (-1)^(2j) / (2j+1).
 */
Site::ElasticCoefficients GenMatrixFunc::ElasticCoefficientsOf(const Site *point, const azure2::EntrancePathways &P,
                                                              const azure2::ExitRow &theDecay, int aa) {
  Site::ElasticCoefficients e;
  e.nuclear = this->LegendreCoefficientsOf(point, P, theDecay, aa, aa);  // builds rt_it_ too (BuildOrderProducts)
  const PPair *pair = compound()->GetPair(aa);
  const double geom = point->GetGeometricalFactor();
  e.coulomb = geom / 100.;
  const bool spinDependent = pair->IsIdentical() && pair->HasSpinDependentExchange();
  if (spinDependent) {
    const double j = pair->GetJ(1);
    e.cross = 2.0 * e.coulomb * ((std::lround(2.0 * j) % 2 == 0) ? 1.0 : -1.0) / (2.0 * j + 1.0);
  }
  const double itFactor = pair->IsIdentical() ? 2.0 : 1.0;
  const complex cI = complex(0., 1.) / sqrt(pi) * geom * itFactor / 100.;
  for (const CoulombTerm &t : rt_it_) {
    if (t.l > Config::maxLOrder)
      throw std::runtime_error("elastic scattering of pair " + std::to_string(aa) + ": a channel with l = " +
                               std::to_string(t.l) + " above the highest Legendre order " +
                               std::to_string(Config::maxLOrder));
    if ((int)e.direct.size() <= t.l) e.direct.resize(t.l + 1, complex(0.0, 0.0));
    e.direct[t.l] += cI * t.c;
    if (!spinDependent) continue;
    if ((int)e.exchange.size() <= t.l) e.exchange.resize(t.l + 1, complex(0.0, 0.0));
    // channel spins of two identical particles are integers (GetCoulombAmplitude(s))
    e.exchange[t.l] += (((std::lround(2.0 * t.s) / 2) % 2 == 0) ? 1.0 : -1.0) * cI * t.c;
  }
  return e;
}

/*!
 * Creates a new temporary T-Matrix element.
 */

void GenMatrixFunc::NewTempTMatrix(TempTMatrix tempTMatrix) {
  temp_t_matrices_.push_back(tempTMatrix);
}

void GenMatrixFunc::NewTempTMatrixE1(TempTMatrix tempTMatrix) {
  temp_t_matrices_E1_.push_back(tempTMatrix);
}

void GenMatrixFunc::NewTempTMatrixE2(TempTMatrix tempTMatrix) {
  temp_t_matrices_E2_.push_back(tempTMatrix);
}

/*!
 * Adds a value to the temporary T-Matrix element specified by its position in the TempTMatrix vector.
 */

void GenMatrixFunc::AddToTempTMatrix(int tempTMatrixNum, complex tempValue) {
  this->GetTempTMatrix(tempTMatrixNum)->TMatrix += tempValue;
}

void GenMatrixFunc::AddToTempTMatrixE1(int tempTMatrixNum, complex tempValue) {
  this->GetTempTMatrixE1(tempTMatrixNum)->TMatrix += tempValue;
}

void GenMatrixFunc::AddToTempTMatrixE2(int tempTMatrixNum, complex tempValue) {
  this->GetTempTMatrixE2(tempTMatrixNum)->TMatrix += tempValue;
}

/*!
 * Clears the temporary T-Matrices.
 */

void GenMatrixFunc::ClearTempTMatrices() {
  temp_t_matrices_.clear();
}

void GenMatrixFunc::ClearTempTMatricesE1() {
  temp_t_matrices_E1_.clear();
}

void GenMatrixFunc::ClearTempTMatricesE2() {
  temp_t_matrices_E2_.clear();
}

/*!
 * Adds an internal T-Matrix element to the vector of internal T-matrix elements
 * corresponding to a specified internal reaction pathway.
 */

void GenMatrixFunc::AddTMatrixElement(int kGroupNum, int mGroupNum, complex tMatrixElement, int decayNum) {
  rt_decay_ = nullptr;
  matrix_c d;
  vector_c e;
  while (decayNum > tmatrix_.size()) tmatrix_.push_back(d);
  while (kGroupNum > tmatrix_[decayNum - 1].size()) tmatrix_[decayNum - 1].push_back(e);
  tmatrix_[decayNum - 1][kGroupNum - 1].push_back(tMatrixElement);
  assert(kGroupNum == tmatrix_[decayNum - 1].size());
  assert(mGroupNum == tmatrix_[decayNum - 1][kGroupNum - 1].size());
}

/*!
 * Adds an external T-Matrix element to the vector of external T-matrix elements
 * corresponding to a specified external reaction pathway.
 */

void GenMatrixFunc::AddECTMatrixElement(int kGroupNum, int mGroupNum, complex tMatrixElement) {
  rt_decay_ = nullptr;
  vector_c d;
  while (kGroupNum > ec_tmatrix_.size()) ec_tmatrix_.push_back(d);
  ec_tmatrix_[kGroupNum - 1].push_back(tMatrixElement);
  assert(mGroupNum == ec_tmatrix_[kGroupNum - 1].size());
}

/*!
 * Tests if a temporary T-Matrix element already exists for a given \f$ J,l,l' \f$ combination.
 * If the element exists, returns the position in the TempTMatrix vector, otherwise returns 0.
 */

int GenMatrixFunc::IsTempTMatrix(double jValue, int lValue, int lPrimeValue) {
  int d = 0;
  bool e = false;
  while (!e && d < this->NumTempTMatrices()) {
    if (jValue == this->GetTempTMatrix(d + 1)->jValue &&
        lValue == this->GetTempTMatrix(d + 1)->lValue &&
        lPrimeValue == this->GetTempTMatrix(d + 1)->lpValue) e = true;
    d++;
  }
  if (!e)
    return 0;
  else
    return d;
}

int GenMatrixFunc::IsTempTMatrixE1(double jValue, int lValue, int lPrimeValue) {
  int d = 0;
  bool e = false;
  while (!e && d < this->NumTempTMatricesE1()) {
    if (jValue == this->GetTempTMatrixE1(d + 1)->jValue &&
        lValue == this->GetTempTMatrixE1(d + 1)->lValue &&
        lPrimeValue == this->GetTempTMatrixE1(d + 1)->lpValue) e = true;
    d++;
  }
  if (!e)
    return 0;
  else
    return d;
}

int GenMatrixFunc::IsTempTMatrixE2(double jValue, int lValue, int lPrimeValue) {
  int d = 0;
  bool e = false;
  while (!e && d < this->NumTempTMatricesE2()) {
    if (jValue == this->GetTempTMatrixE2(d + 1)->jValue &&
        lValue == this->GetTempTMatrixE2(d + 1)->lValue &&
        lPrimeValue == this->GetTempTMatrixE2(d + 1)->lpValue) e = true;
    d++;
  }
  if (!e)
    return 0;
  else
    return d;
}

/*!
 * Returns the number of temporary T-Matrix elements in the TempTMatrix vector.
 */

int GenMatrixFunc::NumTempTMatrices() const {
  return temp_t_matrices_.size();
}

int GenMatrixFunc::NumTempTMatricesE1() const {
  return temp_t_matrices_E1_.size();
}

int GenMatrixFunc::NumTempTMatricesE2() const {
  return temp_t_matrices_E2_.size();
}

/*!
 * Returns a pointer to the temporary T-Matrix element specified by a position in the TempTMatrix vector.
 */

TempTMatrix *GenMatrixFunc::GetTempTMatrix(int tempTMatrixNum) {
  TempTMatrix *b = &temp_t_matrices_[tempTMatrixNum - 1];
  return b;
}

TempTMatrix *GenMatrixFunc::GetTempTMatrixE1(int tempTMatrixNum) {
  TempTMatrix *b = &temp_t_matrices_E1_[tempTMatrixNum - 1];
  return b;
}

TempTMatrix *GenMatrixFunc::GetTempTMatrixE2(int tempTMatrixNum) {
  TempTMatrix *b = &temp_t_matrices_E2_[tempTMatrixNum - 1];
  return b;
}

/*!
 * Returns the value of the internal T-Matrix element specified by an internal reaction pathway.
 */

complex GenMatrixFunc::GetTMatrixElement(int kGroupNum, int mGroupNum, int decayNum) const {
  return tmatrix_[decayNum - 1][kGroupNum - 1][mGroupNum - 1];
}

/*!
 * Returns the value of the external T-Matrix element specified by an external reaction pathway.
 */

complex GenMatrixFunc::GetECTMatrixElement(int kGroupNum, int ecMGroupNum) const {
  return ec_tmatrix_[kGroupNum - 1][ecMGroupNum - 1];
}

/*!
 * Calculates the R_k coefficient for the unobserved primary, observed secondary (UPOS)
 * angular distribution calculation.
 *
 * R_k = sqrt(2j2f+1)*sqrt(2L+1)*sqrt(2L'+1)*(-1)^(j2f-Ic+L-L'+k+1)
 *       * ClebGord(L',L,k,1,-1,0) * Racah(L,L',j2f,j2f,k,Ic)
 */

double GenMatrixFunc::GetRk(double j2f, double finalL, double finalLp, double Ic, int lOrder) {
  AngCoeff angCoeff;
  return sqrt(2. * j2f + 1.) * sqrt(2. * finalL + 1.) * sqrt(2. * finalLp + 1.) *
      pow(-1., j2f - Ic + finalL - finalLp + lOrder + 1) *
      angCoeff.ClebGord(finalLp, finalL, lOrder, 1., -1., 0.) *
      angCoeff.Racah(finalL, finalLp, j2f, j2f, lOrder, Ic);
}


std::string GenMatrixFunc::NoPathwayMessage(int entrancePair, int exitPair) const {
  return "no reaction pathway leads from pair " + std::to_string(compound()->GetPair(entrancePair)->GetPairKey()) +
         " to pair " + std::to_string(compound()->GetPair(exitPair)->GetPairKey()) +
         " (no R-matrix level couples them, and no external capture leads there)";
}

std::string GenMatrixFunc::CaptureAyNotBuiltMessage(const Site *point) const {
  return "the capture analyzing power of data set " + std::to_string(point->GetSegmentKey()) + " (pair " +
         std::to_string(point->GetEntranceKey()) + " to pair " + std::to_string(point->GetExitKey()) +
         ") has no terms: they are built at set-up for the reactions the data sets name (RequestCaptureAyTerms), and "
         "this one was not named";
}

Site::AnalyzingPowerCoefficients GenMatrixFunc::AnalyzingPowerCoefficientsOf(const Site *point,
                                                                             const azure2::EntrancePathways &P,
                                                                             const azure2::ExitRow &theDecay, int aa,
                                                                             int exitPairNum) {
  Site::AnalyzingPowerCoefficients a;
  // the cross section that weights A_y: the differential one, by the route of every other point (no angle here)
  a.elastic = aa == exitPairNum;
  if (a.elastic)
    a.weightElastic = this->ElasticCoefficientsOf(point, P, theDecay, aa);
  else
    a.weightLegendre = this->LegendreCoefficientsOf(point, P, theDecay, aa, exitPairNum);
  const PPair *in = compound()->GetPair(aa), *out = compound()->GetPair(exitPairNum);
  a.capture = true;  // without a particle or capture form (no projectile in), the empty capture form: A_y = 0
  if (in->GetPType() != 0) return a;
  if (out->GetPType() == 10) {
    // Seyler and Weller (1979): the coefficients of P_k and P_k^1, summed over the terms per order
    if (!theDecay.has_capture_ay) throw std::runtime_error(CaptureAyNotBuiltMessage(point));
    int orders = 0;
    for (int i = 0; i < theDecay.num_capture_ay; i++)
      orders = std::max(orders, P.capture_ay[theDecay.first_capture_ay + i].kOrder + 1);
    a.denominator.assign(orders, 0.0);
    a.numerator.assign(orders, 0.0);
    for (int i = 0; i < theDecay.num_capture_ay; i++) {
      const azure2::CaptureAyTerm *term = &P.capture_ay[theDecay.first_capture_ay + i];
      complex T1 = term->isEC1 ? this->GetECTMatrixElement(term->kGroup1, term->path1)
                               : this->GetTMatrixElement(term->kGroup1, term->path1);
      complex T2 = term->isEC2 ? this->GetECTMatrixElement(term->kGroup2, term->path2)
                               : this->GetTMatrixElement(term->kGroup2, term->path2);
      const complex prod = T1 * conj(T2);
      if (term->ak != 0.0) a.denominator[term->kOrder] += term->ak * real(prod);
      if (term->bk != 0.0) a.numerator[term->kOrder] += term->bk * (-imag(prod));  // Re(i T T'*)
    }
    return a;
  }
  if (out->GetPType() != 0) return a;
  // Seyler (1969): the amplitude matrix's coefficients per exit orbital momentum
  a.capture = false;
  a.spins = {in->GetJ(1), in->GetJ(2), out->GetJ(1), out->GetJ(2)};
  Polarization::AmplitudeMatrix M(compound(), aa, exitPairNum);
  a.amplitudes = (int)M.size();
  a.orders = 1;
  for (int k = 1; k <= theDecay.num_groups; k++) {
    const azure2::SpinGroupRow &g = P.Group(theDecay, k);
    if (g.CrossSpin()) continue;  // a cross-spin UPOS group repeats pathways of the equal-spin groups
    for (int m = 1; m <= g.num_internal; m++) {
      const azure2::InternalRow &r = P.Internal(g, m);
      a.orders = std::max(a.orders, compound()->GetJGroup(r.group)->GetChannel(r.chp)->GetL() + 1);
    }
  }
  a.coefficients.assign((std::size_t)a.amplitudes * a.orders, complex(0.0, 0.0));
  for (int k = 1; k <= theDecay.num_groups; k++) {
    const azure2::SpinGroupRow &g = P.Group(theDecay, k);
    if (g.CrossSpin()) continue;
    for (int m = 1; m <= g.num_internal; m++) {
      const azure2::InternalRow &r = P.Internal(g, m);
      M.AddPathwayCoefficients(r.group, r.ch, r.chp, this->GetTMatrixElement(k, m), a.coefficients.data(), a.orders);
    }
  }
  return a;
}

/*!
 * The Coulomb-nuclear interference of the elastic differential cross section,
 * sum over pathways with equal entrance and exit channel of
 * statSpin f_C(s) conj(T) P_l, summed per (s, l) without the angle-dependent
 * f_C(s) P_l, for the T-matrix held now.  Only for an elastic decay.
 */
void GenMatrixFunc::BuildCoulombInterference(const azure2::EntrancePathways &P, const azure2::ExitRow &exit,
                                             bool elastic) {
  rt_it_.clear();
  if (!elastic) return;
  // The pathways with equal entrance and exit channel; an exit without interference lists (a beta entrance) has none
  // and is never elastic.
  for (int n = 0; n < exit.num_diagonal; n++) {
    const azure2::DiagonalPathway &d = P.diagonal[exit.first_diagonal + n];
    size_t i = 0;
    while (i < rt_it_.size() && !(rt_it_[i].s == d.s && rt_it_[i].l == d.l)) i++;
    if (i == rt_it_.size()) rt_it_.push_back(CoulombTerm{d.s, d.l, complex(0.0, 0.0)});
    rt_it_[i].c += d.statSpin * conj(this->GetTMatrixElement(d.k, d.m));
  }
}
