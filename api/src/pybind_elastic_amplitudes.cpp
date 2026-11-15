// The door of the nuclear amplitude expansion of charged-particle elastic scattering for pyazr.endf (ENDF-102, 2023
// edition, LAW=5 LTP=1, pp. 145-147; docs/upgrade/endf.md section 8.2, row "nuclear amplitude coefficients"): the
// coefficients b_l and a_l of eqs. 6.13 and 6.14 at an elastic point, read from the record the engine wrote at the
// point's energy.  Nothing is computed again: the record is what the evaluation contracts into dsigma/dOmega.
//
// Where the numbers come from.  The collision matrix of the point's energy is the engine's T-matrix (the A-matrix
// route, AMatrixFunc::CalculateTMatrix, or RMatrixFunc::CalculateTMatrix; GenMatrixFunc::GetTMatrixElement).
// GenMatrixFunc::ElasticCoefficientsOf (src/GenMatrixFunc.cpp) forms from it the kElastic record of the site
// (include/SiteObservable.h, elastic_record): the nuclear B_L of |F_N|^2 and the complex Coulomb-nuclear D_l (and
// X_l for identical particles with spin), so that the engine's
//   dsigma/dOmega = c_0 (|f_d|^2 + |f_x|^2) + c_x Re(f_d f_x*) + sum_L B_L P_L + Re sum_l (f_d D_l + f_x X_l) P_l
// (ContractElastic, src/SiteObservables.cpp) with the engine's Coulomb amplitude (Site::CoulombAmplitudesAt,
// src/ChannelTables.cpp)
//   f_d = eta / (sqrt(pi) (1 - mu)) exp(-i eta ln((1 - mu)/2))   (distinguishable particles; + f_C(pi - theta) for
//                                                                  identical particles without spin)
//   f_x = eta / (sqrt(pi) (1 + mu)) exp(-i eta ln((1 + mu)/2))   (identical particles with spin), c_0 = pi/(100 k^2).
//
// The manual's form (eq. 6.13, p. 145; eta, k of eqs. 6.11-6.12, the c.m. cosine mu):
//   sigma_e = sigma_c - (2 eta / (1 - mu)) Re{exp(i eta ln((1 - mu)/2)) sum_l (2l+1)/2 a_l P_l}
//             + sum_{l=0}^{2NL} (2l+1)/2 b_l P_l
// Its interference term equals the engine's Re sum_l f_d D_l P_l = Re sum_l conj(f_d) conj(D_l) P_l term by term when
//   a_l = -conj(D_l) / (sqrt(pi) (2l + 1)),          and its nuclear term the engine's sum B_L P_L when
//   b_l = 2 B_l / (2l + 1).
// Identical particles (eq. 6.14): the interference term is the same sum with (1 + mu) exp(i eta ln((1-mu)/2)) +
// (-1)^l (1 - mu) exp(i eta ln((1+mu)/2)) over 1 - mu^2, which is the engine's f_d + f_x term with the same a_l where
// X_l = (-1)^l D_l: every channel of two identical particles has l + s even (CompoundSystem admits no other), so the
// exchange sign (-1)^s of X_l is (-1)^l (checked here; a record that breaks it is refused).  The nuclear term
// sum_{l=0}^{NL} (4l+1)/2 b_l P_{2l} is sum_L B_L P_L with b_l = 2 B_{2l} / (4l + 1), the odd B_L zero (checked).
// The spins of the projectile and the target enter only through the engine's T-matrix and its spin average (the
// record); the manual's a_l expand the trace of the amplitude matrix over the spin states divided by their number, the
// Coulomb amplitude being diagonal in them, which is what D_l is (statSpin of the diagonal pathways).
//
// What is returned: the values of the LIST record of LTP=1 in the manual's order (p. 147), barns/sr:
//   LIDP=0: b_0 ... b_{2NL}, Re a_0, Im a_0, ..., Re a_NL, Im a_NL   (NW = 4 NL + 3)
//   LIDP=1: b_0 ... b_NL,    Re a_0, Im a_0, ..., Re a_NL, Im a_NL   (NW = 3 NL + 3)
// NL the highest l of the record's D_l (or of its B_L, whichever asks more).  The reader's eta and k must be the
// engine's for the sum to be the engine's cross section: the writer refuses relativistic projects.
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "Constants.h"
#include "Results.h"
#include "SiteObservable.h"
#include "Unit.h"

namespace {

// The largest departure from the identities a record of identical particles must meet, relative to its largest
// coefficient (the record's own rounding is about 1e-16 of it).
constexpr double kIdentityTolerance = 1e-10;

}  // namespace

/// The LIST values of LAW=5 LTP=1 (the file's text) at `point`'s energy from the kElastic record it holds in
/// `results`; empty where the point holds no elastic record of its own (another kind of observable, or a point whose
/// value is a fold of sub-points).
std::vector<double> NuclearAmplitudeCoefficients(const Unit &point, const azure2::Results &results) {
  const Site &site = point.OwnSite();
  if (site.Observable() != SiteObservable::kElastic || site.ResultAngularSlots() == 0 || point.NumSubPoints() > 0)
    return {};
  // The record's angular part: LegendreCoefficients gives its first slot after the count, kNumB's place + 1.
  const double *record = site.LegendreCoefficients(results) - 1;
  using namespace elastic_record;
  if (std::isnan(record[kNumB]) || std::isnan(record[kNumL])) return {};
  const int n = (int)record[kNumB], m = (int)record[kNumL];
  const double *B = record + kFirst, *D = B + n, *X = D + 2 * m;
  const bool identical = site.EntranceSymmetry() == 4.0;  // LegendreSymmetryFactor of identical elastic scattering
  const bool exchange = record[kCross] != 0.0;
  double scale = 0.0;
  for (int L = 0; L < n; L++) scale = std::max(scale, std::fabs(B[L]));
  for (int l = 0; l < 2 * m; l++) scale = std::max(scale, std::fabs(D[l]));
  if (exchange) {
    for (int l = 0; l < m; l++) {
      const double sign = (l % 2 == 0) ? 1.0 : -1.0;
      if (std::fabs(X[2 * l] - sign * D[2 * l]) > kIdentityTolerance * scale ||
          std::fabs(X[2 * l + 1] - sign * D[2 * l + 1]) > kIdentityTolerance * scale)
        throw std::domain_error("elastic scattering of identical particles: the exchange coefficient X_" +
                                std::to_string(l) + " is not (-1)^l D_l, which eq. 6.14 of ENDF-102 (p. 145) needs");
    }
  }
  if (identical)
    for (int L = 1; L < n; L += 2)
      if (std::fabs(B[L]) > kIdentityTolerance * scale)
        throw std::domain_error("elastic scattering of identical particles: an odd Legendre coefficient B_" +
                                std::to_string(L) + " of the nuclear cross section, which eq. 6.14 has no place for");
  // NL: the highest l of the interference sum, at least as high as the nuclear sum's orders ask (2 NL or NL)
  int NL = m > 0 ? m - 1 : 0;
  const int highestB = n - 1;
  if (identical)
    NL = std::max(NL, highestB / 2);
  else
    NL = std::max(NL, (highestB + 1) / 2);
  std::vector<double> out;
  const int nb = identical ? NL + 1 : 2 * NL + 1;
  for (int l = 0; l < nb; l++) {
    const int L = identical ? 2 * l : l;
    out.push_back(L < n ? 2.0 * B[L] / (2.0 * L + 1.0) : 0.0);
  }
  for (int l = 0; l <= NL; l++) {
    const complex d = l < m ? complex(D[2 * l], D[2 * l + 1]) : complex(0.0, 0.0);
    const complex a = -std::conj(d) / (std::sqrt(pi) * (2.0 * l + 1.0));
    out.push_back(a.real());
    out.push_back(a.imag());
  }
  return out;
}
