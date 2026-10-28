// AChannel::SetWignerLimit: gamma_W^2 = 3 hbar^2 / (2 mu a^2) for a particle channel, and no limit
// (0, which CNuc::FillMnParams skips) for a photon or beta channel or a channel without a radius.
// A photon pair has mu = 0, so until 2026-10-27 a photon channel with a nonzero radius got an
// infinite limit, and FillMnParams set Minuit limits (-inf, inf), which its sine transform turns
// into NaN (audit 2026-10-25, item 32; the archive's photon pairs all have radius 0, so none hit it).
#include "AChannel.h"
#include "Constants.h"
#include <cmath>
#include <cstdio>
static int failures = 0;
static void check(bool ok, const char *what) {
  std::printf("  %s  %s\n", ok ? "ok  " : "FAIL", what);
  if (!ok) failures++;
}
int main() {
  // AChannel(l, s, pairNum, radType)
  AChannel particle(1, 0.5, 1, 'P');
  const double mu = 1.00728 * 12.0 / 13.00728, a = 4.0;
  particle.SetWignerLimit(mu, a);
  const double expected = 3.0 * hbarc * hbarc / (2.0 * mu * uconv * a * a);
  check(std::fabs(particle.GetWignerLimit() - expected) < 1e-9 * expected, "particle channel: 3 hbar^2 / (2 mu a^2)");
  AChannel noRadius(1, 0.5, 1, 'P');
  noRadius.SetWignerLimit(mu, 0.0);
  check(noRadius.GetWignerLimit() == 0.0, "no radius: no limit");
  AChannel photon(1, 0.5, 2, 'E');
  photon.SetWignerLimit(0.0, 4.0);  // a photon pair's reduced mass is 0
  check(photon.GetWignerLimit() == 0.0 && std::isfinite(photon.GetWignerLimit()), "photon channel with a radius: no limit (was inf)");
  AChannel magnetic(1, 0.5, 2, 'M');
  magnetic.SetWignerLimit(0.0, 4.0);
  check(magnetic.GetWignerLimit() == 0.0, "M photon channel: no limit");
  AChannel beta(0, 0.5, 3, 'F');
  beta.SetWignerLimit(0.000548, 4.0);  // a beta pair: a finite but meaningless number before
  check(beta.GetWignerLimit() == 0.0, "beta channel: no limit");
  std::printf(failures ? "FAILED: %d\n" : "PASS\n", failures);
  return failures ? 1 : 0;
}
