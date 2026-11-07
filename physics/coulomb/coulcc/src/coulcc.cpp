// COULCC driver: the main routine of coulcc.f (I. J. Thompson and A. R.
// Barnett, Comput. Phys. Commun. 36 (1985) 363) with structured control flow.
//
// The local variables of the Fortran routine, its COMMON blocks /STEED/ and
// /RCFCM1/ and the work array XRCF are the members of `Steed`, one object per
// call.  The labelled sections of the Fortran are member functions that return
// where the computation continues (`Flow`); the Fortran labels are given in
// the comments.  The arithmetic is that of the Fortran, statement by
// statement (fortran_complex.h).  DESIGN.md section 3.
#include "coulcc.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>

#include "kernels.h"

namespace coulcc {

namespace {

using namespace detail;

// The public warning bits are the internal ones.
static_assert(unsigned(Warning::ReflectionNotUsedL) == kWarnReflectionNotUsedL, "warning bit");
static_assert(unsigned(Warning::ReflectionNotUsedX) == kWarnReflectionNotUsedX, "warning bit");
static_assert(unsigned(Warning::LinearIndependenceLost) == kWarnLinearIndependenceLost, "warning bit");
static_assert(unsigned(Warning::IrregularOverflow) == kWarnIrregularOverflow, "warning bit");
static_assert(unsigned(Warning::IrregularUnderflow) == kWarnIrregularUnderflow, "warning bit");
static_assert(unsigned(Warning::RegularOverflow) == kWarnRegularOverflow, "warning bit");
static_assert(unsigned(Warning::RegularUnderflow) == kWarnRegularUnderflow, "warning bit");
static_assert(unsigned(Warning::Cf1aAndCf1Unstable) == kWarnCf1aAndCf1Unstable, "warning bit");
static_assert(unsigned(Warning::F11Overflow) == kWarnF11Overflow, "warning bit");
static_assert(unsigned(Warning::PolesAndInstabilities) == kWarnPolesAndInstabilities, "warning bit");
static_assert(unsigned(Warning::Cf1AccuracyLoss) == kWarnCf1AccuracyLoss, "warning bit");
static_assert(unsigned(Warning::Cf1NotConverged) == kWarnCf1NotConverged, "warning bit");
static_assert(unsigned(Warning::Cf2NotConverged) == kWarnCf2NotConverged, "warning bit");
static_assert(unsigned(Warning::LogGammaPole) == kWarnLogGammaPole, "warning bit");
static_assert(unsigned(Warning::RcfRestart) == kWarnRcfRestart, "warning bit");
static_assert(unsigned(Warning::IrregularTooSmallRho) == kWarnIrregularTooSmallRho, "warning bit");
static_assert(unsigned(Warning::RoundoffWarning) == kWarnRoundoff, "warning bit");

// Parameters determining the region of each method (DATA in coulcc.f).
constexpr double kR20 = 3.0;    // estimate of (2F0 iterations)/(CF2 iterations)
constexpr double kAsym = 3.0;   // minimum X/(ETA**2+L) for CF1A to converge easily
constexpr double kXnear = 0.5;  // minimum ABS(X) for CF2 to converge accurately
constexpr int kNdrop = 5;       // number of successive decrements that define instability

// Where the computation continues after a section.
enum class Flow {
  next,      // the following section
  reverse,   // label 60: redo all lambda values in the reverse direction
  soft,      // label 280: IFAIL from the number of lambda values completed
  hard,      // label 390: IFAIL = -1
};

class Steed {
 public:
  Steed(complex xx, complex eta1, complex zlmin, int nl, complex* fc, complex* gc, complex* fcp, complex* gcp,
        complex* sig, int mode1, int kfn)
      : xx_(xx), eta1_(eta1), zlmin_(zlmin), nl_(nl), mode1_(mode1), kfn_(kfn), fc_(fc), gc_(gc), fcp_(fcp),
        gcp_(gcp), sig_(sig) {}

  int run();
  void report(Info& info) const;

 private:
  // 1-based access to the output arrays, as in the Fortran.
  complex& FC(int l) { return fc_[l - 1]; }
  complex& GC(int l) { return gc_[l - 1]; }
  complex& FCP(int l) { return fcp_[l - 1]; }
  complex& GCP(int l) { return gcp_[l - 1]; }
  complex& SIG(int l) { return sig_[l - 1]; }

  void start_block();
  Flow regular_at_top();
  Flow recur_regular();
  Flow irregular_at_bottom();
  Flow irregular_by_kase();
  Flow from_series();
  Flow normalise_and_recur();

  // Arguments.
  fcomplex xx_, eta1_, zlmin_;
  int nl_, mode1_, kfn_;
  complex *fc_, *gc_, *fcp_, *gcp_, *sig_;

  // Constants set at entry.
  double accur_ = 0, acct_ = 0, acch_ = 0, accb_ = 0, hpi_ = 0, tlog_ = 0;
  int mode_ = 0;
  bool ifcp_ = false;

  // COMMON /STEED/.
  double rerr_ = 0;
  int nfp_ = 0, n11_ = 0, npq_[3] = {0, 0, 0}, n20_ = 0, kas_[3] = {0, 0, 0};  // 1-based NPQ, KAS
  unsigned warnings_ = 0;

  // Quantities fixed for the call.
  fcomplex cik_, x_, eta_, etai_, zm1_, xi_, xlog_;
  bool etane0_ = false, rlel_ = false, axial_ = false;
  double dell_ = 0, scale_ = 0, absx_ = 0;
  int m1_ = 1, l1_ = 1;

  // State of the block (lambda range) being computed.
  int ifail_ = -2, id_ = 1, lf_ = 1, last_ = 1;
  bool donem_ = false, unstab_ = false, zlneg_ = false, nocf2_ = false;
  fcomplex zlm_, zll_, z11_, p11_, sigma_, cll_, theta_, thetam_, first_;

  // COMMON /RCFCM1/ and the remaining locals of COULCC.
  fcomplex pk_, ek_, clgaa_, clgab_, clgbb_, dsig_, tpk1_, w_, rl_, fcl1_, q_, gam_, hcl_, hpl_, fcm_, hcl1_,
      alpha_, beta_, pl_;
  fcomplex aa_, ab_, bb_, p_, f_, fest_, fesl_, fcl_, fpl_, zl_, sl_, pm_, etap_, pq1_, pq2_, f11v_, f20v_, chi_,
      zlog_, df_;
  double err_ = 0, ta_ = 0, rk_ = 0, sf_ = 0, off_ = 0, omega_ = 0, sfsh_ = 0, eps_ = 0, paccq_ = 0;
  int i_ = 0, l_ = 0, mono_ = 0, lh_ = 1, kase_ = 0, n_ = 0, ih_ = 0;
  PadeWork work_;  // XRCF(JMAX,4)
};

// ---------------------------------------------------------------------------
int Steed::run() {
  mode_ = std::abs(mode1_) % 10;
  ifcp_ = mode_ % 2 == 1;
  ifail_ = -2;
  hpi_ = 2.0 * std::atan(1.0);
  tlog_ = std::log(2.0);
  accur_ = std::max(1e-14, 50 * kAcc8);
  acct_ = accur_ * 0.5;
  acch_ = std::sqrt(accur_);
  accb_ = std::sqrt(acch_);
  rerr_ = acct_;

  cik_ = 1.0;
  if (kfn_ >= 3) cik_ = CI * sign(1.0, kAcc8 - imag(xx_));
  x_ = xx_ * cik_;
  eta_ = eta1_;
  if (kfn_ > 0) eta_ = 0.0;
  etane0_ = absc(eta_) > kAcc8;
  etai_ = eta_ * CI;
  dell_ = 0.0;
  if (kfn_ >= 2) dell_ = 0.5;
  zm1_ = zlmin_ - dell_;
  scale_ = 0.0;
  if (mode1_ < 0) scale_ = imag(x_);

  m1_ = 1;
  l1_ = m1_ + nl_ - 1;
  rlel_ = std::fabs(imag(eta_)) + std::fabs(imag(zm1_)) < kAcc8;
  absx_ = abs(x_);
  axial_ = rlel_ && std::fabs(imag(x_)) < kAcc8 * absx_;
  if (mode_ <= 2 && absx_ < kFpmin) {  // 310
    warnings_ |= kWarnIrregularTooSmallRho;
    return ifail_;
  }
  xi_ = 1.0 / x_;
  xlog_ = log(x_);  // cut along the negative real axis
  id_ = 1;
  donem_ = false;
  unstab_ = false;
  lf_ = m1_;
  ifail_ = -1;

  // One block per pass: lambda from ZLL down to ZLM (ID = 1), or the reverse.
  for (;;) {
    start_block();  // 10
    Flow flow = regular_at_top();
    if (flow == Flow::next) flow = recur_regular();
    if (flow == Flow::reverse) {  // 60
      id_ = -1;
      lf_ = l1_;
      l1_ = m1_;
      rerr_ = acct_;
      continue;
    }
    if (flow == Flow::next) flow = irregular_at_bottom();
    if (flow == Flow::next) flow = normalise_and_recur();
    // flow is now soft (280, also the normal end) or hard (390).
    if (flow == Flow::hard) {
      ifail_ = -1;  // 390
    } else {
      if (id_ > 0 || last_ == 0) ifail_ = last_;  // 280
      if (id_ < 0 && last_ != 0) ifail_ = -3;
    }
    if (id_ > 0 && lf_ != m1_) {  // 290 -> 300
      // On the first block F started decreasing monotonically, or hit bound
      // states for low lambda: redo M1 to LF-1 in the reverse direction.
      id_ = -1;
      if (!unstab_) lf_ = lf_ - 1;
      donem_ = unstab_;
      lf_ = std::min(lf_, l1_);
      l1_ = m1_;
      continue;
    }
    break;
  }
  if (ifail_ < 0) return ifail_;
  if (rerr_ > accb_) warnings_ |= kWarnRoundoff;
  if (rerr_ > double(0.1f)) ifail_ = -4;  // RERR.GT.0.1: a single-precision constant
  return ifail_;
}

// 10: the lambda range of this block; phase shift and Gamow factor at ZLL.
void Steed::start_block() {
  zlm_ = zm1_ + lf_ - m1_;
  zll_ = zm1_ + l1_ - m1_;  // final lambda value, or 0.5 smaller for J, Y Bessels
  z11_ = zll_;
  if (id_ < 0) z11_ = zlm_;
  p11_ = CI * sign(1.0, kAcc8 - imag(eta_));
  last_ = l1_;

  pk_ = zll_ + 1.0;
  aa_ = pk_ - etai_;
  ab_ = pk_ + etai_;
  bb_ = 2.0 * pk_;
  zlneg_ = npint(bb_, accb_);
  clgaa_ = clogam(aa_, warnings_);
  clgab_ = clgaa_;
  if (etane0_ && !rlel_) clgab_ = clogam(ab_, warnings_);
  if (etane0_ && rlel_) clgab_ = conj(clgaa_);
  sigma_ = (clgaa_ - clgab_) * CI * 0.5;
  if (kfn_ == 0) SIG(l1_) = sigma_.std();
  if (!zlneg_) cll_ = zll_ * tlog_ - hpi_ * eta_ - clogam(bb_, warnings_) + (clgaa_ + clgab_) * 0.5;
  theta_ = x_ - eta_ * (xlog_ + tlog_) - zll_ * hpi_ + sigma_;
}

// 216-275: f = F'/F at ZLL by CF1A or CF1, then F(ZLL) with the right sign.
Flow Steed::regular_at_top() {
  ta_ = (imag(aa_) * imag(aa_) + imag(ab_) * imag(ab_) + std::fabs(real(aa_)) + std::fabs(real(ab_))) * 0.5;
  bool cf1 = id_ > 0 && absx_ < ta_ * kAsym && !zlneg_;
  bool skip_to_40 = false;
  if (!cf1) {
    // CF1A treats negative lambda and recurrence-unstable cases properly.
    rk_ = sign(1.0, real(x_) + kAcc8);
    p_ = theta_;
    if (rk_ < 0) p_ = -x_ + eta_ * (log(-x_) + tlog_) - zll_ * hpi_ - sigma_;
    f_ = rk_ * cf1a(x_ * rk_, eta_ * rk_, zll_, p_, acct_, kJmax, nfp_, fest_, err_, kFpmax, work_, warnings_);
    fesl_ = log(fest_) + std::fabs(imag(x_));
    nfp_ = -nfp_;
    if (nfp_ < 0 || (unstab_ && err_ < accb_)) {
      skip_to_40 = true;
    } else if (!zlneg_ || (unstab_ && err_ > accb_)) {
      cf1 = true;  // 20
    } else {
      warnings_ |= kWarnReflectionNotUsedL;
      if (err_ > accb_) return Flow::soft;
      skip_to_40 = true;
    }
  }
  if (!skip_to_40) {
    // 20: CF1 = f = F'(ZLL,ETA,X)/F(ZLL,ETA,X).
    if (axial_) {
      f_ = cf1r(real(x_), real(eta_), real(zll_), kAcc8, sf_, rk_, etane0_, kLimit, err_, nfp_, acch_, kFpmin,
                kFpmax, warnings_);
      fcl_ = sf_;
      tpk1_ = rk_;
    } else {
      f_ = cf1c(x_, eta_, zll_, kAcc8, fcl_, tpk1_, etane0_, kLimit, err_, nfp_, acch_, kFpmin, kFpmax, warnings_);
    }
    if (err_ > 1.0) return Flow::hard;
    // A simple check for CF1 being badly unstable.
    if (id_ >= 0) {
      unstab_ = real((1.0 - eta_ * xi_) * CI * imag(theta_) / f_) > 0.0 && !axial_ &&
                std::fabs(imag(theta_)) > -(std::log(kAcc8) * 0.5) && absc(eta_) + absc(zll_) < absc(x_);
      if (unstab_) return Flow::reverse;
    }
    // 30: compare the accumulated phase FCL with the asymptotic phase for
    // G(k+1) to estimate F(ZLL), with the correct sign, to start the recursion.
    w_ = x_ * x_ * (0.5 / tpk1_ + 1.0 / (tpk1_ * tpk1_)) + eta_ * (eta_ - 2.0 * x_) / tpk1_;
    fesl_ = (zll_ + 1.0) * xlog_ + cll_ - w_ - log(fcl_);
  }
  fesl_ = fesl_ - std::fabs(scale_);  // 40
  rk_ = std::max(real(fesl_), kFplmin * 0.5);
  fesl_ = fcomplex(std::min(rk_, kFplmax * 0.5), imag(fesl_));
  fest_ = exp(fesl_);
  rerr_ = std::max(std::max(rerr_, err_), kAcc8 * std::fabs(real(theta_)));
  fcl_ = fest_;
  fpl_ = fcl_ * f_;
  if (ifcp_) FCP(l1_) = fpl_.std();
  FC(l1_) = fcl_.std();
  return Flow::next;
}

// 276-353: recurrence of F, F' to lambda = ZLM; GC (and GCP) store the
// coefficients RL (and PL) for the irregular recurrence.
Flow Steed::recur_regular() {
  i_ = std::max(-id_, 0);
  zl_ = zll_ + i_;
  mono_ = 0;
  off_ = abs(fcl_);
  ta_ = absc(sigma_);
  bool stop = false;  // label 50: cannot or should not recur below this ZL
  for (l_ = l1_ - id_; id_ > 0 ? l_ >= lf_ : l_ <= lf_; l_ -= id_) {
    if (etane0_) {
      if (rlel_) {
        dsig_ = std::atan2(real(eta_), real(zl_));
        rl_ = std::sqrt(real(zl_) * real(zl_) + real(eta_) * real(eta_));
      } else {
        aa_ = zl_ - etai_;
        bb_ = zl_ + etai_;
        if (absc(aa_) < acch_ || absc(bb_) < acch_) {
          stop = true;
          break;
        }
        dsig_ = (log(aa_) - log(bb_)) * CI * 0.5;
        rl_ = aa_ * exp(CI * dsig_);
      }
      if (absc(sigma_) < ta_ * 0.5) {
        // Recalculate SIGMA because of accumulating roundoff.
        sl_ = (clogam(zl_ + i_ - etai_, warnings_) - clogam(zl_ + i_ + etai_, warnings_)) * CI * 0.5;
        rl_ = (zl_ - etai_) * exp(CI * id_ * (sigma_ - sl_));
        sigma_ = sl_;
        ta_ = 0.0;
      } else {
        sigma_ = sigma_ - dsig_ * id_;
      }
      ta_ = std::max(ta_, absc(sigma_));
      sl_ = eta_ + zl_ * zl_ * xi_;
      pl_ = 0.0;
      if (absc(zl_) > acch_) pl_ = (sl_ * sl_ - rl_ * rl_) / zl_;
      fcl1_ = (fcl_ * sl_ + id_ * zl_ * fpl_) / rl_;
      sf_ = abs(fcl1_);
      if (sf_ > kFpmax) {  // 350
        warnings_ |= kWarnRegularOverflow;
        return Flow::soft;
      }
      fpl_ = (fpl_ * sl_ + id_ * pl_ * fcl_) / rl_;
      if (mode_ <= 1) GCP(l_ + id_) = (pl_ * id_).std();
    } else {
      // ETA = 0, including Bessels; RL = SL.
      rl_ = zl_ * xi_;
      fcl1_ = fcl_ * rl_ + fpl_ * id_;
      sf_ = abs(fcl1_);
      if (sf_ > kFpmax) {  // 350
        warnings_ |= kWarnRegularOverflow;
        return Flow::soft;
      }
      fpl_ = (fcl1_ * rl_ - fcl_) * id_;
    }
    if (sf_ < off_) {
      mono_ = mono_ + 1;
    } else {
      mono_ = 0;
    }
    fcl_ = fcl1_;
    off_ = sf_;
    FC(l_) = fcl_.std();
    if (ifcp_) FCP(l_) = fpl_.std();
    if (kfn_ == 0) SIG(l_) = sigma_.std();
    if (mode_ <= 2) GC(l_ + id_) = rl_.std();
    zl_ = zl_ - id_;
    if (mono_ < kNdrop) continue;
    if (axial_ || (real(zlm_) * id_ > -kNdrop && !etane0_)) continue;
    unstab_ = true;
    stop = true;
    break;
  }
  if (!stop) return Flow::next;  // 80
  zlm_ = zl_;  // 50
  lf_ = l_;
  if (id_ < 0) {  // 380
    warnings_ |= kWarnPolesAndInstabilities;
    return Flow::hard;
  }
  if (!unstab_) lf_ = l_ + 1;
  if (l_ + mono_ < l1_ - 2 || id_ < 0 || !unstab_) return Flow::next;  // 80
  return Flow::reverse;  // all lambda values (for stability) in the reverse direction
}

// 80-130: f at ZLM; choice of omega and of the method (KASE) for the
// irregular solution H(omega) at ZLM.
Flow Steed::irregular_at_bottom() {
  if (fcl_ == 0.0) fcl_ = kAcc8;  // 80
  f_ = fpl_ / fcl_;
  // Second time around: check that the f values agree.
  if (id_ > 0) first_ = f_;
  if (donem_) rerr_ = std::max(rerr_, absc(f_ - first_) / absc(f_));
  if (donem_) {  // 90 with L = 1, then 130
    lh_ = 1;
    if (omega_ < 0.0) lh_ = 2;
    pm_ = CI * omega_;
    etap_ = eta_ * pm_;
    return irregular_by_kase();
  }
  nocf2_ = false;
  thetam_ = x_ - eta_ * (xlog_ + tlog_) - zlm_ * hpi_ + sigma_;
  // On the left x-plane, OMEGA puts the cut on the negative x axis; on the
  // right, OMEGA (estimated from THETAM) makes H(omega) the smaller solution,
  // which recurs upwards accurately.  The boundary is shifted to give CF2 a
  // chance to converge.
  omega_ = sign(1.0, imag(x_) + kAcc8);
  if (real(x_) >= kXnear) omega_ = sign(1.0, imag(thetam_) + kAcc8);
  sfsh_ = std::exp(omega_ * scale_ - std::fabs(scale_));
  off_ = std::exp(std::min(
      2.0 * std::max(std::max(std::fabs(imag(x_)), std::fabs(imag(thetam_))), std::fabs(imag(zlm_)) * 3), kFplmax));
  eps_ = std::max(kAcc8, acct_ * 0.5 / off_);

  // 90: try the estimated omega, then its opposite, to find the H(omega)
  // linearly independent of F, i.e. maximise CF1-CF2 = 1/(F H(omega)).
  enum class Then { accept, choose, series } then = Then::accept;
  bool found = false;
  for (int trial = 1; trial <= 2 && !found; ++trial) {
    lh_ = 1;
    if (omega_ < 0.0) lh_ = 2;
    pm_ = CI * omega_;
    etap_ = eta_ * pm_;
    pq1_ = 0.0;
    paccq_ = 1.0;
    kase_ = 0;
    // Small X: avoid CF2.
    if (mode_ >= 3 && absx_ < 1.0) {
      then = Then::series;  // 190
      found = true;
      break;
    }
    if (mode_ < 3 && (nocf2_ || (absx_ < kXnear && absc(eta_) * absx_ < 5 && absc(zlm_) < 4))) {
      kase_ = 5;
      then = Then::choose;  // 120
      found = true;
      break;
    }
    // CF2: PQ1 = p + i omega q at lambda = ZLM.
    pq1_ = cf2(x_, eta_, zlm_, pm_, eps_, kLimit, err_, npq_[lh_], kAcc8, acch_, accur_, warnings_);
    err_ = err_ * std::max(1.0, absc(pq1_) / std::max(absc(f_ - pq1_), kAcc8));
    if (err_ < acch_) {
      then = Then::accept;  // 110
      found = true;
      break;
    }
    // Impossible to get F-PQ accurately because of cancellation?
    nocf2_ = real(x_) < kXnear && std::fabs(imag(x_)) < -std::log(kAcc8);
    omega_ = -omega_;  // 100: the first guess of omega was wrong
  }
  if (!found) {
    if (unstab_) {  // 360
      warnings_ |= kWarnCf1aAndCf1Unstable;
      return Flow::soft;
    }
    if (real(x_) < -kXnear) warnings_ |= kWarnReflectionNotUsedX;
  }
  if (then == Then::series) {
    // 190
    f11v_ = f11(x_, eta_, z11_, p11_, acct_, kLimit, 0, err_, n11_, kFpmax, kAcc8, kAcc16, warnings_);
    return from_series();
  }
  if (then == Then::accept) rerr_ = std::max(rerr_, err_);  // 110
  // 120: the case for the irregular solution.
  if (kase_ < 5) {
    if (real(x_) > kXnear) {
      // Estimate the errors if KASE 2 or 3 were used.
      paccq_ = eps_ * off_ * absc(pq1_) / std::max(std::fabs(imag(pq1_)), kAcc8);
    }
    if (paccq_ < accur_) {
      kase_ = 2;
      if (axial_) kase_ = 3;
    } else {
      kase_ = 1;
      if (npq_[1] * kR20 < kJmax) kase_ = 4;  // the 2F0 is predicted to converge
    }
  }
  return irregular_by_kase();
}

// 130-180: H(omega) and F at ZLM by the method KASE.
Flow Steed::irregular_by_kase() {
  for (;;) {
    const int kase = std::abs(kase_);  // 130
    if (kase == 2 || kase == 3) {
      if (kase == 2) {  // 140: PQ2 = p - i omega q
        if (!donem_)
          pq2_ = cf2(x_, eta_, zlm_, -pm_, eps_, kLimit, err_, npq_[3 - lh_], kAcc8, acch_, accur_, warnings_);
        p_ = (pq2_ + pq1_) * 0.5;
        q_ = (pq2_ - pq1_) * 0.5 * pm_;
      } else {  // 150: on the real axis p and q are real and PQ2 = PQ1*
        p_ = real(pq1_);
        q_ = imag(pq1_);
        pq2_ = conj(pq1_);
      }
      // 160: solve for FCM = F at ZLM, then the norm factor W = FCM/FCL.
      w_ = (pq1_ - f_) * (pq2_ - f_);
      sf_ = std::exp(-std::fabs(scale_));
      fcm_ = sqrt(q_ / w_) * sf_;
      // The sign of the square root: FCM nearest to the phase of FCL.
      if (real(fcm_ / fcl_) < 0.0) fcm_ = -fcm_;
      gam_ = (f_ - p_) / q_;
      ta_ = absc(gam_ + pm_);
      paccq_ = eps_ * std::max(ta_, 1.0 / ta_);
      hcl_ = fcm_ * (gam_ + pm_) * (sfsh_ / (sf_ * sf_));
      if (paccq_ > accur_ && kase_ > 0) {
        // Consider a KASE = 1 calculation.
        f11v_ = f11(x_, eta_, z11_, p11_, acct_, kLimit, 0, err_, n11_, kFpmax, kAcc8, kAcc16, warnings_);
        if (err_ < paccq_) return from_series();  // 200
      }
      rerr_ = std::max(rerr_, paccq_);
      return Flow::next;  // 230
    }
    if (kase == 4) {
      // 170: the exponentially decreasing H(LH) directly from 2F0.
      bool series = false;
      if (!donem_) {
        aa_ = etap_ - zlm_;
        bb_ = etap_ + zlm_ + 1.0;
        f20v_ = f20(aa_, bb_, -(0.5 * pm_ * xi_), acct_, kJmax, err_, kFpmax, n20_, work_, warnings_);
        if (n20_ <= 0) {
          series = true;  // 190
        } else {
          rerr_ = std::max(rerr_, err_);
          hcl_ = kFpmin;
          if (std::fabs(real(pm_ * thetam_) + omega_ * scale_) > kFplmax) {  // 330
            warnings_ |= kWarnIrregularUnderflow;
            return Flow::soft;
          }
        }
      }
      if (!series) {
        hcl_ = f20v_ * exp(pm_ * thetam_ + omega_ * scale_);  // 180
        fcm_ = sfsh_ / ((f_ - pq1_) * hcl_);
        return Flow::next;  // 230
      }
    }
    // 190 (KASE 1, 5, 6, or 2F0 tried and failed): F at Z11 from the 1F1
    // series, where Z11 = ZLL if ID > 0 and ZLM if ID < 0.
    f11v_ = f11(x_, eta_, z11_, p11_, acct_, kLimit, 0, err_, n11_, kFpmax, kAcc8, kAcc16, warnings_);
    const int before = kase_;
    const Flow flow = from_series();
    if (kase_ < 0 && kase_ != before) continue;  // 200 asked for KASE 2 or 3 (-2, -3): back to 130
    return flow;
  }
}

// 200-220: F (and for KASE 5, 6 the irregular solution) from the 1F1 series.
Flow Steed::from_series() {
  if (n11_ < 0) {  // F11 failed: BB a negative integer
    warnings_ |= kWarnReflectionNotUsedL;
    return Flow::hard;  // 390
  }
  if (err_ > paccq_ && paccq_ < accb_) {
    // Consider a KASE 2 or 3 calculation.
    kase_ = -2;
    if (axial_) kase_ = -3;
    return Flow::next;  // the caller returns to 130
  }
  rerr_ = std::max(rerr_, err_);
  if (err_ > kFpmax) {  // 370
    warnings_ |= kWarnF11Overflow;
    return Flow::hard;
  }
  if (id_ < 0)
    cll_ = z11_ * tlog_ - hpi_ * eta_ - clogam(bb_, warnings_) + clogam(z11_ + 1.0 + p11_ * eta_, warnings_) -
           p11_ * sigma_;
  ek_ = (z11_ + 1.0) * xlog_ - p11_ * x_ + cll_ - std::fabs(scale_);
  if (id_ > 0) ek_ = ek_ - fesl_ + log(fcl_);
  if (real(ek_) > kFplmax) {  // 350
    warnings_ |= kWarnRegularOverflow;
    return Flow::soft;
  }
  if (real(ek_) < kFplmin) {  // 340
    warnings_ |= kWarnRegularUnderflow;
    return Flow::soft;
  }
  fcm_ = f11v_ * exp(ek_);

  if (kase_ < 5) {
    if (mode_ <= 2) hcl_ = sfsh_ / ((f_ - pq1_) * fcm_);
    kase_ = 1;
    return Flow::next;  // 230
  }
  // |X| < XNEAR: CF2 may not converge accurately, so use an expansion of the
  // irregular solution from the origin.
  if (absc(zlm_ + zlm_ - nintc(zlm_ + zlm_)) < acch_) kase_ = 6;
  sl_ = zlm_;
  zlneg_ = real(zlm_) < -1.0 + accb_;
  if (kase_ == 5 || zlneg_) sl_ = -zlm_ - 1.0;
  pk_ = sl_ + 1.0;
  aa_ = pk_ - etap_;
  ab_ = pk_ + etap_;
  bb_ = 2.0 * pk_;
  clgaa_ = clogam(aa_, warnings_);
  clgab_ = clgaa_;
  if (etane0_) clgab_ = clogam(ab_, warnings_);
  clgbb_ = clogam(bb_, warnings_);
  if (kase_ == 6 && !zlneg_) {
    if (npint(aa_, accur_)) clgaa_ = clgab_ - 2.0 * pm_ * sigma_;
    if (npint(ab_, accur_)) clgab_ = clgaa_ + 2.0 * pm_ * sigma_;
  }
  cll_ = sl_ * tlog_ - hpi_ * eta_ - clgbb_ + (clgaa_ + clgab_) * 0.5;
  dsig_ = (clgaa_ - clgab_) * pm_ * 0.5;
  if (kase_ == 6) p11_ = -pm_;
  ek_ = pk_ * xlog_ - p11_ * x_ + cll_ - std::fabs(scale_);
  sf_ = std::exp(-std::fabs(scale_));
  chi_ = 0.0;
  if (kase_ == 5 || zlneg_) {
    // G(l) = (cos(CHI) F(l) - F(-l-1)) / sin(CHI), CHI = sig(l) - sig(-l-1) - (2l+1) pi/2.
    chi_ = sigma_ - dsig_ - (zlm_ - sl_) * hpi_;
    f11v_ = f11(x_, eta_, sl_, p11_, acct_, kLimit, 0, err_, npq_[1], kFpmax, kAcc8, kAcc16, warnings_);
    rerr_ = std::max(rerr_, err_);
    if (kase_ != 6) {
      fesl_ = f11v_ * exp(ek_);
      fcl1_ = exp(pm_ * chi_) * fcm_;
      hcl_ = fcl1_ - fesl_;
      rerr_ = std::max(rerr_, acct_ * std::max(absc(fcl1_), absc(fesl_)) / absc(hcl_));
      hcl_ = hcl_ / sin(chi_) * (sfsh_ / (sf_ * sf_));
      pq1_ = f_ - sfsh_ / (fcm_ * hcl_);  // 220
      return Flow::next;
    }
  }
  // 210: the logarithmic expansion of the irregular solution (KASE 6), for
  // integral BB where sin(CHI) would be zero.
  rl_ = bb_ - 1.0;
  n_ = nintc(rl_);
  zlog_ = xlog_ + tlog_ - pm_ * hpi_;
  chi_ = chi_ + pm_ * thetam_ + omega_ * scale_ + ab_ * zlog_;
  aa_ = 1.0 - aa_;
  if (npint(aa_, accur_)) {
    hcl_ = 0.0;
  } else {
    if (id_ > 0 && !zlneg_) f11v_ = fcm_ * exp(-ek_);
    const fcomplex psi_series =
        f11(x_, eta_, sl_, -pm_, acct_, kLimit, 2, err_, npq_[2], kFpmax, kAcc8, kAcc16, warnings_);
    const int sign_n = ((n_ + 1) & 1) ? -1 : 1;  // (-1)**(N+1)
    hcl_ = exp(chi_ - clgbb_ - clogam(aa_, warnings_)) * sign_n * (f11v_ * zlog_ + psi_series);
    rerr_ = std::max(rerr_, err_);
  }
  if (n_ > 0) {
    ek_ = chi_ + clogam(rl_, warnings_) - clgab_ - rl_ * zlog_;
    int nits = 0;  // the Fortran writes it into L, which is not used afterwards
    df_ = f11(x_, eta_, -sl_ - 1.0, -pm_, 0.0, n_, 0, err_, nits, kFpmax, kAcc8, kAcc16, warnings_);
    hcl_ = hcl_ + exp(ek_) * df_;
  }
  pq1_ = f_ - sfsh_ / (fcm_ * hcl_);  // 220
  return Flow::next;
}

// 230-270: absolute normalisation at ZLM, the recurrence of the irregular
// solution, and the transformation to the functions requested.
Flow Steed::normalise_and_recur() {
  ih_ = std::abs(mode1_) / 10;
  if (kfn_ == 3) ih_ = int((3 - imag(cik_)) / 2 + 0.5);
  p11_ = 1.0;
  if (ih_ == 1) p11_ = CI;
  if (ih_ == 2) p11_ = -CI;
  df_ = -pm_;
  if (ih_ >= 1) df_ = -pm_ + p11_;
  if (absc(df_) < acch_) df_ = 0.0;
  // Normalisations of spherical and cylindrical Bessel functions.
  alpha_ = 0.0;
  if (kfn_ == 1) alpha_ = xi_;
  if (kfn_ >= 2) alpha_ = xi_ * 0.5;
  beta_ = 1.0;
  if (kfn_ == 1) beta_ = xi_;
  if (kfn_ >= 2) beta_ = sqrt(xi_ / hpi_);
  if (kfn_ >= 2 && real(beta_) < 0.0) beta_ = -beta_;

  aa_ = 1.0;
  if (kfn_ > 0) aa_ = -(p11_ * beta_);
  if (kfn_ >= 3) {
    // Rescaling factors for I and K.
    p_ = exp((zlm_ + dell_) * hpi_ * cik_);
    aa_ = beta_ * hpi_ * p_;
    beta_ = beta_ / p_;
    q_ = cik_ * id_;
  }
  // Rescaling factors for the GC output.
  if (ih_ == 0) {
    ta_ = std::fabs(scale_) + imag(pm_) * scale_;
    rk_ = 0.0;
    if (ta_ < kFplmax) rk_ = std::exp(-ta_);
  } else {
    ta_ = std::fabs(scale_) + imag(p11_) * scale_;
    if (absc(df_) > acch_ && ta_ > kFplmax) {  // 320
      warnings_ |= kWarnIrregularOverflow;
      return Flow::soft;
    }
    if (absc(df_) > acch_) df_ = df_ * std::exp(ta_);
    sf_ = 2.0 * (lh_ - ih_) * scale_;
    rk_ = 0.0;
    if (sf_ > kFplmax) {  // 320
      warnings_ |= kWarnIrregularOverflow;
      return Flow::soft;
    }
    if (sf_ > kFplmin) rk_ = std::exp(sf_);
  }

  kas_[(3 - id_) / 2] = kase_;
  w_ = fcm_ / fcl_;
  if (std::log(absc(w_)) + std::log(absc(fcomplex(FC(lf_)))) < kFplmin) {  // 340
    warnings_ |= kWarnRegularUnderflow;
    return Flow::soft;
  }
  if (mode_ < 3) {
    if (absc(f_ - pq1_) < acch_ * absc(f_)) warnings_ |= kWarnLinearIndependenceLost;
    hpl_ = hcl_ * pq1_;
    if (absc(hpl_) < kFpmin || absc(hcl_) < kFpmin) {  // 330
      warnings_ |= kWarnIrregularUnderflow;
      return Flow::soft;
    }
  }
  // 240: recurrence in the direction ID from HCL, HPL at LF; renormalise FC,
  // FCP at each lambda (ZL = ZLM - MIN(ID,0) here).
  for (l_ = lf_; id_ > 0 ? l_ <= l1_ : l_ >= l1_; l_ += id_) {
    fcl_ = w_ * FC(l_);
    if (absc(fcl_) < kFpmin) {  // 340
      warnings_ |= kWarnRegularUnderflow;
      return Flow::soft;
    }
    if (ifcp_) fpl_ = w_ * FCP(l_);
    fcomplex fc = beta_ * fcl_;
    FC(l_) = fc.std();
    if (ifcp_) FCP(l_) = (beta_ * (fpl_ - alpha_ * fcl_) * cik_).std();
    FC(l_) = tidy(FC(l_), accur_).std();
    if (ifcp_) FCP(l_) = tidy(FCP(l_), accur_).std();
    if (mode_ < 3) {
      if (l_ != lf_) {
        zl_ = zl_ + id_;
        const fcomplex zid = zl_ * id_;
        rl_ = GC(l_);
        if (etane0_) {
          sl_ = eta_ + zl_ * zl_ * xi_;
          if (mode_ == 1) {
            pl_ = GCP(l_);
          } else {
            pl_ = 0.0;
            if (absc(zl_) > acch_) pl_ = (sl_ * sl_ - rl_ * rl_) / zid;
          }
          hcl1_ = (sl_ * hcl_ - zid * hpl_) / rl_;
          hpl_ = (sl_ * hpl_ - pl_ * hcl_) / rl_;
        } else {
          hcl1_ = rl_ * hcl_ - hpl_ * id_;
          hpl_ = (hcl_ - rl_ * hcl1_) * id_;
        }
        hcl_ = hcl1_;
        if (absc(hcl_) > kFpmax) {  // 320
          warnings_ |= kWarnIrregularOverflow;
          return Flow::soft;
        }
      }
      GC(l_) = (aa_ * (rk_ * hcl_ + df_ * fcl_)).std();  // 250
      if (mode_ == 1) GCP(l_) = ((aa_ * (rk_ * hpl_ + df_ * fpl_) - alpha_ * fcomplex(GC(l_))) * cik_).std();
      GC(l_) = tidy(GC(l_), accur_).std();
      if (mode_ == 1) GCP(l_) = tidy(GCP(l_), accur_).std();
      if (kfn_ >= 3) aa_ = aa_ * q_;
    }
    if (kfn_ >= 3) beta_ = -(beta_ * q_);  // 260
    last_ = std::min(last_, (l1_ - l_) * id_);  // 270
  }
  return Flow::soft;  // 280
}

void Steed::report(Info& info) const {
  info.rerr = rerr_;
  info.warnings = warnings_;
  info.diag.nfp = nfp_;
  info.diag.n11 = n11_;
  info.diag.npq[0] = npq_[1];
  info.diag.npq[1] = npq_[2];
  info.diag.n20 = n20_;
  info.diag.kase[0] = kas_[1];
  info.diag.kase[1] = kas_[2];
}

bool valid_arguments(int nl, int mode1, int kfn) {
  const int mode = std::abs(mode1) % 10, ih = std::abs(mode1) / 10;
  return nl >= 1 && mode >= 1 && mode <= 4 && ih <= 2 && kfn >= -1 && kfn <= 3;
}

}  // namespace

// ---------------------------------------------------------------------------
int Options::mode1() const {
  const int mode = irregular ? (derivatives ? 1 : 2) : (derivatives ? 3 : 4);
  const int m = 10 * static_cast<int>(second) + mode;
  return scaled ? -m : m;
}

Options Options::from_fortran(int mode1, int kfn) {
  Options o;
  const int mode = std::abs(mode1) % 10;
  o.kind = static_cast<Functions>(kfn);
  o.second = static_cast<Irregular>(std::abs(mode1) / 10);
  o.irregular = mode <= 2;
  o.derivatives = mode % 2 == 1;
  o.scaled = mode1 < 0;
  return o;
}

int coulcc(complex xx, complex eta1, complex zlmin, int nl, complex* fc, complex* gc, complex* fcp, complex* gcp,
           complex* sig, int mode1, int kfn, Info& info) {
  info = Info();
  if (!valid_arguments(nl, mode1, kfn)) return static_cast<int>(Status::InvalidArgument);
  Steed steed(xx, eta1, zlmin, nl, fc, gc, fcp, gcp, sig, mode1, kfn);
  const int ifail = steed.run();
  steed.report(info);
  return ifail;
}

Result coulcc(complex rho, complex eta, complex lambda_min, int nl, const Options& options) {
  Result r;
  const int mode1 = options.mode1(), kfn = static_cast<int>(options.kind);
  if (!valid_arguments(nl, mode1, kfn)) {
    r.ifail = static_cast<int>(Status::InvalidArgument);
    r.status = Status::InvalidArgument;
    return r;
  }
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const complex blank(nan, nan);
  const int mode = std::abs(mode1) % 10;
  r.F.assign(nl, blank);
  if (mode <= 2) r.G.assign(nl, blank);
  if (mode % 2 == 1) r.Fp.assign(nl, blank);
  if (mode == 1) r.Gp.assign(nl, blank);
  if (kfn == 0) r.sigma.assign(nl, blank);
  Info info;
  r.ifail = coulcc(rho, eta, lambda_min, nl, r.F.data(), r.G.empty() ? nullptr : r.G.data(),
                   r.Fp.empty() ? nullptr : r.Fp.data(), r.Gp.empty() ? nullptr : r.Gp.data(),
                   r.sigma.empty() ? nullptr : r.sigma.data(), mode1, kfn, info);
  r.status = r.ifail > 0 ? Status::Partial : static_cast<Status>(r.ifail);
  r.n_valid = r.ifail >= 0 ? nl - r.ifail : 0;
  r.rerr = info.rerr;
  r.warnings = info.warnings;
  r.diag = info.diag;
  for (auto* v : {&r.F, &r.G, &r.Fp, &r.Gp, &r.sigma})
    for (int i = r.n_valid; i < static_cast<int>(v->size()); ++i) (*v)[i] = blank;
  return r;
}

}  // namespace coulcc
