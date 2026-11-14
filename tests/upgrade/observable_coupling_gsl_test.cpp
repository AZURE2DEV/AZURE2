// The observable coupling coefficients of physics/ (ObservableCoupling.h) against the engine's expressions before
// step 5 (CompoundSystem::CalcAngularDists, CompoundSystem::CalcCaptureAnalyzingPower, ThmAngular's Zbar at 9be284f0), written here
// with AngCoeff's Clebsch-Gordan and Racah coefficients and GSL's 9j, an independent implementation that physics/
// does not link (docs/upgrade/cnuc_split.md section 9 (b)).  For every coefficient of each kind: the difference is
// at most 1e-14 of the kind's largest magnitude, and the engine's keep test (|c| > 1e-10) gives the same answer.
// Quantum numbers: l, L <= 4, 2j, 2s <= 8, k <= 6.
#include <gsl/gsl_sf_coupling.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "AngCoeff.h"
#include "azure2/physics/ObservableCoupling.h"

using namespace azure2::physics;

namespace {

int failures = 0;

struct Kind {
  std::string name;
  std::vector<double> got, want;
};

void Report(const Kind &kind) {
  double largest = 0.0, diff = 0.0;
  for (double w : kind.want) largest = std::fmax(largest, std::fabs(w));
  long keepMismatch = 0, nonzero = 0;
  for (size_t i = 0; i < kind.got.size(); i++) {
    diff = std::fmax(diff, std::fabs(kind.got[i] - kind.want[i]));
    if ((std::fabs(kind.got[i]) > 1e-10) != (std::fabs(kind.want[i]) > 1e-10)) keepMismatch++;
    if (std::fabs(kind.want[i]) > 1e-10) nonzero++;
  }
  const bool ok = diff <= 1e-14 * largest && keepMismatch == 0;
  if (!ok) failures++;
  std::printf("  %s %-24s %8zu coefficients (%ld kept), largest %.3g, largest difference %.2g (%.2g relative), "
              "keep mismatches %ld\n",
              ok ? "ok  " : "FAIL", kind.name.c_str(), kind.got.size(), nonzero, largest, diff,
              largest > 0 ? diff / largest : 0.0, keepMismatch);
}

// The engine's expressions before step 5, verbatim but for the variable names.
double OldZ(double l1, double j1, double l2, double j2, double s, double L) {
  return sqrt(2. * l1 + 1.) * sqrt(2. * l2 + 1.) * sqrt(2. * j1 + 1.) * sqrt(2. * j2 + 1.) *
         AngCoeff::ClebGord(l1, l2, L, 0., 0., 0.) * AngCoeff::Racah(l1, j1, l2, j2, s, L);
}

double OldZPhoton(double L1, double j1, double L2, double j2, double jf, double L) {
  return sqrt(2. * L1 + 1.) * sqrt(2. * L2 + 1.) * sqrt(2. * j1 + 1.) * sqrt(2. * j2 + 1.) *
         AngCoeff::ClebGord(L1, L2, L, 1., -1., 0) * AngCoeff::Racah(L1, j1, L2, j2, jf, L);
}

double OldUpos(double l1, double l1p, double j1, double l2, double j2, double s, double sp, double sp2, double j1f,
               double j2f, double lOrder) {
  return pow(-1., sp - sp2) * pow(-1., lOrder + sp2 - sp) * (2. * j1 + 1.) * (2. * j2 + 1.) *
         pow((2. * l1 + 1.) * (2. * j2f + 1.) * (2. * sp + 1.) * (2. * sp2 + 1.), 0.5) *
         AngCoeff::ClebGord(lOrder, l1, l2, 0., 0., 0.) * AngCoeff::Racah(lOrder, j2f, sp2, j1f, j2f, sp) *
         AngCoeff::Racah(lOrder, sp, j2, l1p, sp2, j1) * AngCoeff::Racah(lOrder, j1, l2, s, j2, l1);
}

double PhaseM1(double n) { return (std::labs(std::lround(n)) % 2) ? -1.0 : 1.0; }
double Hat(double j) { return sqrt(2.0 * j + 1.0); }

CaptureAyCoefficients OldAy(const CapturePathway &p1, const CapturePathway &p2, int k, double xSpin, double aSpin,
                            double cSpin) {
  CaptureAyCoefficients out;
  const double l1 = p1.l, L1 = p1.L, b1 = p1.two_b / 2.0, s1 = p1.two_s / 2.0;
  const double l2 = p2.l, L2 = p2.L, b2 = p2.two_b / 2.0, s2 = p2.two_s / 2.0;
  if ((int)(L1 + p1.p + L2 + p2.p + k) % 2 != 0) return out;
  double cgL = AngCoeff::ClebGord(l1, l2, (double)k, 0., 0., 0.);
  if (fabs(cgL) < 1.e-12) return out;
  double cgG = AngCoeff::ClebGord(L1, L2, (double)k, 1., -1., 0.);
  if (fabs(cgG) < 1.e-12) return out;
  double wG = AngCoeff::Racah(L1, b1, L2, b2, cSpin, (double)k);
  if (fabs(wG) < 1.e-12) return out;
  double common = cgL * cgG * wG * Hat(l1) * Hat(l2) * Hat(L1) * Hat(L2) * (2. * b1 + 1.) * (2. * b2 + 1.);
  if (fabs(s1 - s2) < 1.e-6) {
    double wL = AngCoeff::Racah(l1, b1, l2, b2, s1, (double)k);
    out.a = PhaseM1(s1 - cSpin + 1.0) * common * wL;
  }
  if (k >= 1) {
    const double pre = 3.0 * sqrt(xSpin) * Hat(xSpin) * Hat((double)k) / sqrt((xSpin + 1.0) * (double)k * ((double)k + 1.0));
    double wS = AngCoeff::Racah(xSpin, s1, xSpin, s2, aSpin, 1.0);
    if (fabs(wS) >= 1.e-12) {
      double x9j = gsl_sf_coupling_9j(2 * p1.l, p1.two_s, p1.two_b, 2 * p2.l, p2.two_s, p2.two_b, 2 * k, 2, 2 * k);
      if (fabs(x9j) >= 1.e-14)
        out.b = pre * common * Hat(s1) * Hat(s2) * wS * x9j * PhaseM1(aSpin - xSpin + cSpin - b1 - s1 + l1) *
                PhaseM1(s1 - s2);
    }
  }
  return out;
}

}  // namespace

int main() {
  const int maxl = 4, maxTwo = 8, maxk = 6;
  Kind z{"Zbar"}, zp{"Zbar photon"}, part{"particle interference"}, phot{"photon interference"}, upos{"UPOS"},
      ayA{"capture a_k"}, ayB{"capture b_k"};

  for (int l1 = 0; l1 <= maxl; l1++)
    for (int l2 = 0; l2 <= maxl; l2++)
      for (int tj1 = 0; tj1 <= maxTwo; tj1++)
        for (int tj2 = 0; tj2 <= maxTwo; tj2++)
          for (int ts = 0; ts <= maxTwo; ts++)
            for (int k = 0; k <= maxk; k++) {
              z.got.push_back(z_bar(l1, tj1, l2, tj2, ts, k));
              z.want.push_back(OldZ(l1, tj1 / 2., l2, tj2 / 2., ts / 2., k));
              if (l1 >= 1 && l2 >= 1) {
                zp.got.push_back(z_bar_photon(l1, tj1, l2, tj2, ts, k));
                zp.want.push_back(OldZPhoton(l1, tj1 / 2., l2, tj2 / 2., ts / 2., k));
              }
            }

  // Interference coefficients: channel spins and level spins of matching parity of 2j (one pair of spins).
  for (int l1 = 0; l1 <= 3; l1++)
    for (int l1p = 0; l1p <= 3; l1p++)
      for (int l2 = 0; l2 <= 3; l2++)
        for (int l2p = 0; l2p <= 3; l2p++)
          for (int tj1 = 0; tj1 <= 7; tj1++)
            for (int tj2 = tj1 % 2; tj2 <= 7; tj2 += 2)
              for (int ts = tj1 % 2; ts <= 6; ts += 2)
                for (int tsp = tj1 % 2; tsp <= 6; tsp += 2)
                  for (int k = 0; k <= maxk; k++) {
                    part.got.push_back(particle_interference(l1, l1p, tj1, l2, l2p, tj2, ts, tsp, k));
                    part.want.push_back(pow(-1.0, tsp / 2. - ts / 2.) / 4. * OldZ(l1, tj1 / 2., l2, tj2 / 2., ts / 2., k) *
                                        OldZ(l1p, tj1 / 2., l2p, tj2 / 2., tsp / 2., k));
                    if (l1p >= 1 && l2p >= 1 && l1p <= 2 && l2p <= 2) {
                      const int tjf = tsp;
                      phot.got.push_back(photon_interference(l1, l1p, tj1, l2, l2p, tj2, ts, tjf, k));
                      phot.want.push_back(pow(-1., 1. + ts / 2. - tjf / 2.) / 4. * OldZ(l1, tj1 / 2., l2, tj2 / 2., ts / 2., k) *
                                          OldZPhoton(l1p, tj1 / 2., l2p, tj2 / 2., tjf / 2., k));
                    }
                  }

  for (int l1 = 0; l1 <= 3; l1++)
    for (int l1p = 0; l1p <= 3; l1p++)
      for (int l2 = 0; l2 <= 3; l2++)
        for (int tj1 = 0; tj1 <= 7; tj1++)
          for (int tj2 = tj1 % 2; tj2 <= 7; tj2 += 2)
            for (int ts = tj1 % 2; ts <= 4; ts += 2)
              for (int tjb : {1, 2})
                for (int tjB : {0, 1, 3})
                  for (int tsp = std::abs(tjb - tjB); tsp <= tjb + tjB; tsp += 2)
                    for (int tsp2 = std::abs(tjb - tjB); tsp2 <= tjb + tjB; tsp2 += 2)
                      for (int k = 0; k <= maxk; k++) {
                        upos.got.push_back(upos_interference(l1, l1p, tj1, l2, tj2, ts, tsp, tsp2, tjb, tjB, k));
                        upos.want.push_back(OldUpos(l1, l1p, tj1 / 2., l2, tj2 / 2., ts / 2., tsp / 2., tsp2 / 2.,
                                                    tjb / 2., tjB / 2., k));
                      }

  for (int tx : {1, 2, 3})
    for (int ta : {0, 1, 2, 3})
      for (int tc : {0, 1, 2, 3, 4}) {
        std::vector<CapturePathway> paths;
        for (int ts = std::abs(tx - ta); ts <= tx + ta; ts += 2)
          for (int l = 0; l <= 3; l++)
            for (int tb = std::abs(2 * l - ts); tb <= 2 * l + ts; tb += 2)
              for (int L = 1; L <= 2; L++)
                for (int p = 0; p <= 1; p++)
                  if (tb <= maxTwo && (tb + tc) % 2 == 0 && std::abs(tb - tc) <= 2 * L && tb + tc >= 2 * L)
                    paths.push_back(CapturePathway{l, L, tb, ts, p});
        for (const auto &t1 : paths)
          for (const auto &t2 : paths)
            for (int k = 0; k <= maxk; k++) {
              const CaptureAyCoefficients got = capture_ay_coefficients(t1, t2, k, tx, ta, tc);
              const CaptureAyCoefficients want = OldAy(t1, t2, k, tx / 2., ta / 2., tc / 2.);
              ayA.got.push_back(got.a);
              ayA.want.push_back(want.a);
              ayB.got.push_back(got.b);
              ayB.want.push_back(want.b);
            }
      }

  for (const Kind *kind : {&z, &zp, &part, &phot, &upos, &ayA, &ayB}) Report(*kind);
  std::printf("%s: %d kind(s) failed\n", failures ? "FAIL" : "PASS", failures);
  return failures ? 1 : 0;
}
