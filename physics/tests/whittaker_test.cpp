// physics/coulomb Whittaker functions (docs/upgrade/coulomb.md section 1.4) against mpmath at 40 digits
// (data/whittaker.csv, from data/make_whittaker_reference.py): l = 0-10, eta_b = 0-1000, rho = 1e-4-200 with
// eta_b rho <= 2000.  Checked: whittaker_scaled's log10 W and rho W'/W at every point; whittaker's W and W' where W
// is inside the double range, and out_of_range where it is not; the argument checks; the same bits from 8 threads
// walking the grid in different orders.
//
// Error measures: W relative (|10^(log10 W - ref) - 1|); rho W'/W relative where |rho W'/W| >= 1 and absolute
// below (it enters the R matrix as S - B, so an absolute error is what matters near 0).
//
// What it detects: a wrong normalization (W ~ e^-rho (2rho)^-eta_b), a wrong sign or factor in the derivative, a
// quadrature that stops early (accuracy regression of more than a decade), state shared between calls.
// Usage: physics_whittaker_test <path to whittaker.csv>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "azure2/physics/coulomb/Whittaker.h"

using namespace azure2::physics::coulomb;

namespace {

struct Row {
  int l;
  double eta, rho, log10W, rhoDlog;
};

std::vector<Row> Read(const std::string &path) {
  std::vector<Row> rows;
  std::ifstream in(path);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    for (char &c : line)
      if (c == ',') c = ' ';
    std::istringstream ss(line);
    Row r{};
    std::string t[4];
    ss >> r.l >> t[0] >> t[1] >> t[2] >> t[3];
    r.eta = std::strtod(t[0].c_str(), nullptr);
    r.rho = std::strtod(t[1].c_str(), nullptr);
    r.log10W = std::strtod(t[2].c_str(), nullptr);
    r.rhoDlog = std::strtod(t[3].c_str(), nullptr);
    rows.push_back(r);
  }
  return rows;
}

int failures = 0;
void check(const std::string &what, bool ok, const std::string &detail = "") {
  std::printf("  %s  %s%s%s\n", ok ? "ok  " : "FAIL", what.c_str(), ok || detail.empty() ? "" : "  -- ", ok ? "" : detail.c_str());
  if (!ok) failures++;
}

// Measured on GCC 11.5, x86-64, -O3 (2026-11-08): W 2.1e-12 (scaled form at eta_b = 1000, rho = 0.01, where |ln W|
// is in the thousands and epsilon |ln W| is the limit of any log-scaled double; 1.2e-13 inside the double range),
// rho W'/W 1.6e-13.  The tolerances are a decade above.
const double kTolW = 2e-11, kTolDlog = 2e-12;

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s whittaker.csv\n", argv[0]);
    return 2;
  }
  const std::vector<Row> rows = Read(argv[1]);
  check("reference read (" + std::to_string(rows.size()) + " points)", rows.size() > 1000);

  // 1. whittaker_scaled everywhere.
  double worstW = 0.0, worstD = 0.0;
  const Row *whereW = nullptr, *whereD = nullptr;
  bool allOk = true;
  std::vector<WhittakerScaled> serial;
  for (const Row &r : rows) {
    const WhittakerScaled s = whittaker_scaled(r.l, r.eta, r.rho);
    serial.push_back(s);
    if (s.status != Status::ok) {
      allOk = false;
      continue;
    }
    const double log10W = std::log10(s.mantissa) + (double)s.log10_scale;
    const double eW = std::fabs(std::expm1((log10W - r.log10W) * std::log(10.0)));
    const double eD = std::fabs(s.rho_dlog - r.rhoDlog) / std::fmax(std::fabs(r.rhoDlog), 1.0);
    if (eW > worstW) worstW = eW, whereW = &r;
    if (eD > worstD) worstD = eD, whereD = &r;
  }
  char buf[256];
  check("whittaker_scaled: status ok at every point", allOk);
  std::snprintf(buf, sizeof buf, "whittaker_scaled: W relative error %.2e (l %d, eta_b %g, rho %g), tolerance %.0e", worstW,
                whereW ? whereW->l : -1, whereW ? whereW->eta : 0.0, whereW ? whereW->rho : 0.0, kTolW);
  check(buf, worstW <= kTolW);
  std::snprintf(buf, sizeof buf, "whittaker_scaled: rho W'/W error %.2e (l %d, eta_b %g, rho %g), tolerance %.0e", worstD,
                whereD ? whereD->l : -1, whereD ? whereD->eta : 0.0, whereD ? whereD->rho : 0.0, kTolDlog);
  check(buf, worstD <= kTolDlog);
  bool mantissas = true;
  for (const WhittakerScaled &s : serial) mantissas = mantissas && s.mantissa >= 1.0 && s.mantissa < 10.0;
  check("whittaker_scaled: mantissa in [1, 10)", mantissas);

  // 2. whittaker: W and W' inside the double range, out_of_range outside.
  double worstPlain = 0.0;
  int inRange = 0, outRange = 0;
  bool statuses = true;
  const double logMin = std::log10(std::numeric_limits<double>::min()), logMax = std::log10(std::numeric_limits<double>::max());
  for (const Row &r : rows) {
    const WhittakerWave w = whittaker(r.l, r.eta, r.rho);
    // well inside the range (two decades of margin for W' = W rho_dlog / rho)
    const double lwp = r.log10W + std::log10(std::fabs(r.rhoDlog) / r.rho + 1e-300);
    const bool inside = r.log10W > logMin + 2 && r.log10W < logMax - 2 && lwp > logMin + 2 && lwp < logMax - 2;
    const bool outside = r.log10W < logMin - 2 || r.log10W > logMax + 2;
    if (inside) {
      inRange++;
      if (w.status != Status::ok) {
        statuses = false;
        continue;
      }
      const double W = std::pow(10.0, r.log10W), Wp = W * r.rhoDlog / r.rho;
      worstPlain = std::fmax(worstPlain, std::fabs(w.W - W) / W);
      worstPlain = std::fmax(worstPlain, std::fabs(w.Wp - Wp) / std::fmax(std::fabs(Wp), W / r.rho));
    } else if (outside) {
      outRange++;
      statuses = statuses && w.status == Status::out_of_range && std::isnan(w.W) && std::isnan(w.Wp);
    }
  }
  std::snprintf(buf, sizeof buf, "whittaker: W, W' relative error %.2e on %d points inside the double range", worstPlain, inRange);
  check(buf, worstPlain <= 10 * kTolW);
  std::snprintf(buf, sizeof buf, "whittaker: ok inside the range, out_of_range (NaN) on the %d points outside", outRange);
  check(buf, statuses && outRange > 0);

  // 3. Argument checks.
  const double nan = std::numeric_limits<double>::quiet_NaN(), inf = std::numeric_limits<double>::infinity();
  bool invalid = true;
  for (const auto &a : std::vector<std::vector<double>>{{-1, 1, 1}, {10001, 1, 1}, {0, -0.1, 1}, {0, nan, 1}, {0, inf, 1}, {0, 1, 0}, {0, 1, -1}, {0, 1, nan}, {0, 1, inf}}) {
    const WhittakerScaled s = whittaker_scaled((int)a[0], a[1], a[2]);
    const WhittakerWave w = whittaker((int)a[0], a[1], a[2]);
    invalid = invalid && s.status == Status::invalid_argument && std::isnan(s.mantissa) && std::isnan(s.rho_dlog) &&
              w.status == Status::invalid_argument && std::isnan(w.W);
  }
  check("invalid arguments: invalid_argument, NaN", invalid);

  // 4. Closed form at eta_b = 0, l = 0: W = e^-rho exactly, rho W'/W = -rho.
  double worstClosed = 0.0;
  for (double rho : {1e-3, 0.1, 1.0, 7.5, 40.0}) {
    const WhittakerWave w = whittaker(0, 0.0, rho);
    worstClosed = std::fmax(worstClosed, std::fabs(w.W / std::exp(-rho) - 1.0));
    worstClosed = std::fmax(worstClosed, std::fabs(w.Wp / -std::exp(-rho) - 1.0) * std::fmin(1.0, rho));
  }
  std::snprintf(buf, sizeof buf, "closed form eta_b = 0, l = 0 (W = e^-rho): %.2e", worstClosed);
  check(buf, worstClosed <= kTolW);

  // 5. Threads: 8 threads, each walking the grid from a different start, against the serial results.
  std::vector<int> same(8, 1);
  std::vector<std::thread> pool;
  for (int t = 0; t < 8; t++)
    pool.emplace_back([&, t]() {
      const std::size_t n = rows.size();
      for (std::size_t k = 0; k < n; k++) {
        const std::size_t i = (k + t * n / 8) % n;
        const WhittakerScaled s = whittaker_scaled(rows[i].l, rows[i].eta, rows[i].rho);
        if (s.status != serial[i].status || s.log10_scale != serial[i].log10_scale ||
            std::memcmp(&s.mantissa, &serial[i].mantissa, sizeof(double)) != 0 ||
            std::memcmp(&s.rho_dlog, &serial[i].rho_dlog, sizeof(double)) != 0)
          same[t] = 0;
      }
    });
  for (std::thread &th : pool) th.join();
  bool bits = true;
  for (int v : same) bits = bits && v;
  check("8 threads, different orders: the same bits as one", bits);

  std::printf("\n%s\n", failures ? "FAIL" : "all checks passed");
  return failures ? 1 : 0;
}
