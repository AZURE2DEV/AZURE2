// The engine's Clebsch-Gordan and Racah coefficients (AngCoeff, since 2026-11-16 those of physics/) against the
// expressions AngCoeff held before: GSL's 3j and 6j symbols with the phases written out, an implementation that
// physics/ does not link.  Every argument set with 2j <= 12 (Clebsch-Gordan, all projections with m1 + m2 = m3 and
// some without) and 2j <= 10 (Racah): the difference is at most 1e-14, and the engine's keep test (|c| > 1e-10)
// gives the same answer.  Where the old expression was not a number (pow(-1, x) with x not an integer, which only
// arguments of a vanishing coefficient give) the new value must be 0.
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_coupling.h>

#include <cmath>
#include <cstdio>

#include "AngCoeff.h"

namespace {

// AngCoeff::ClebGord and AngCoeff::Racah before 2026-11-16, verbatim.
double OldClebGord(double j1, double j2, double j3, double m1, double m2, double m3) {
  m3 = -m3;
  int j1x2 = (int)(2 * j1);
  int j2x2 = (int)(2 * j2);
  int j3x2 = (int)(2 * j3);
  int m1x2 = (int)(2 * m1);
  int m2x2 = (int)(2 * m2);
  int m3x2 = (int)(2 * m3);
  double w3j = gsl_sf_coupling_3j(j1x2, j2x2, j3x2, m1x2, m2x2, m3x2);
  return pow(-1.0, j1 - j2 - m3) * sqrt(2.0 * j3 + 1.) * w3j;
}

double OldRacah(double j1, double j2, double l2, double l1, double j3, double l3) {
  int j1x2 = (int)(2 * j1);
  int j2x2 = (int)(2 * j2);
  int j3x2 = (int)(2 * j3);
  int l1x2 = (int)(2 * l1);
  int l2x2 = (int)(2 * l2);
  int l3x2 = (int)(2 * l3);
  double w6j = gsl_sf_coupling_6j(j1x2, j2x2, j3x2, l1x2, l2x2, l3x2);
  return pow(-1.0, j1 + j2 + l2 + l1) * w6j;
}

struct Tally {
  long n = 0, kept = 0, notANumber = 0, keepMismatch = 0, bad = 0;
  double largest = 0.0, diff = 0.0;
  void Add(double got, double want) {
    n++;
    if (std::isnan(want)) {
      notANumber++;
      if (got != 0.0) bad++;
      return;
    }
    if (!std::isfinite(got)) bad++;
    largest = std::fmax(largest, std::fabs(want));
    diff = std::fmax(diff, std::fabs(got - want));
    if (std::fabs(want) > 1e-10) kept++;
    if ((std::fabs(got) > 1e-10) != (std::fabs(want) > 1e-10)) keepMismatch++;
  }
  bool Report(const char *name) const {
    const bool ok = n > 0 && kept > 0 && diff <= 1e-14 && keepMismatch == 0 && bad == 0;
    std::printf("  %s %-16s %9ld argument sets (%ld kept, %ld where the old expression was not a number), largest "
                "%.3g, largest difference %.2g, keep mismatches %ld, bad %ld\n",
                ok ? "ok  " : "FAIL", name, n, kept, notANumber, largest, diff, keepMismatch, bad);
    return ok;
  }
};

}  // namespace

int main() {
  gsl_set_error_handler_off();
  Tally cg, racah;
  const int maxCg = 12, maxRacah = 10;
  for (int a = 0; a <= maxCg; a++)
    for (int b = 0; b <= maxCg; b++)
      for (int c = 0; c <= maxCg; c++)
        for (int ma = -a; ma <= a; ma += 2)
          for (int mb = -b; mb <= b; mb += 2)
            for (int mc = ma + mb - 2; mc <= ma + mb + 2; mc += 2) {  // m3 = m1 + m2 and one step to each side
              const double j1 = 0.5 * a, j2 = 0.5 * b, j3 = 0.5 * c, m1 = 0.5 * ma, m2 = 0.5 * mb, m3 = 0.5 * mc;
              cg.Add(AngCoeff::ClebGord(j1, j2, j3, m1, m2, m3), OldClebGord(j1, j2, j3, m1, m2, m3));
            }
  for (int a = 0; a <= maxRacah; a++)
    for (int b = 0; b <= maxRacah; b++)
      for (int c = 0; c <= maxRacah; c++)
        for (int d = 0; d <= maxRacah; d++)
          for (int e = 0; e <= maxRacah; e++)
            for (int f = 0; f <= maxRacah; f++)
              racah.Add(AngCoeff::Racah(0.5 * a, 0.5 * b, 0.5 * c, 0.5 * d, 0.5 * e, 0.5 * f),
                        OldRacah(0.5 * a, 0.5 * b, 0.5 * c, 0.5 * d, 0.5 * e, 0.5 * f));
  const bool okCg = cg.Report("Clebsch-Gordan");
  const bool okRacah = racah.Report("Racah");
  if (okCg && okRacah) {
    std::printf("engine coupling coefficients: physics/ and the former GSL expressions agree\n");
    return 0;
  }
  return 1;
}
