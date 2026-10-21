// azure2-coulomb-table: build, validate and describe the universal Coulomb tables
// (CoulombTable).  AZURE2 itself builds the default table with --build-coulomb-table
// (and the GUI from Configure > Build Coulomb Table...); this tool adds the build
// options, validation and info.
//
//   azure2-coulomb-table build FILE [--lmax 10] [--eta 1e-3 100] [--rho 1e-3 100]
//                              [--per-decade 100] [--eta-rho-max 150] [--threads N]
//   azure2-coulomb-table validate FILE [--n 200000] [--seed 1] [--order 4|6|8] [--eta-rho-max X]
//                                 [--limit 1e-7] [--threads N]
//                                 (--eta-rho-max: test only that part of the domain;
//                                  --limit: the interpolation error that passes)
//   azure2-coulomb-table info FILE
//
// build writes FILE.new and renames it over FILE, so jobs that have FILE open
// keep reading the old one.  Without FILE: the engine's default file
// ($AZURE2_COULOMB_TABLE, else ~/.azure2/coulomb_table.bin).

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <random>
#include <string>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_coulomb.h>

#include "Config.h"
#include "CoulombTable.h"
#include "ShftFunc.h"

#ifdef _OPENMP
#include <omp.h>
#endif

Config *g_config = nullptr;  // the engine library expects the executable to define it

namespace {

int Usage() {
  std::fprintf(stderr,
               "FILE defaults to $AZURE2_COULOMB_TABLE, else ~/.azure2/coulomb_table.bin\n"
               "usage: azure2-coulomb-table build [FILE] [--lmax L] [--eta LO HI] [--rho LO HI] [--per-decade N] [--eta-rho-max X] [--threads N]\n"
               "       azure2-coulomb-table validate [FILE] [--n N] [--seed S] [--order 4|6|8] [--eta-rho-max X] [--limit E] [--threads N]\n"
               "       azure2-coulomb-table info [FILE]\n");
  return 2;
}

double Wrap(double d) {
  d = std::fmod(d, 2 * M_PI);
  if (d > M_PI) d -= 2 * M_PI;
  if (d < -M_PI) d += 2 * M_PI;
  return d;
}

bool LoadOrSay(CoulombTable &t, const std::string &path) {
  std::string err;
  if (!t.Load(path, err)) {
    std::fprintf(stderr, "ERROR: %s\n", err.c_str());
    return false;
  }
  return true;
}

void Info(const CoulombTable &t, const std::string &path) {
  const auto &e = t.EtaAxis(), &r = t.RhoAxis();
  std::printf("%s\n  generator: %s\n  l = 0..%d\n  eta: 10^%.4f .. 10^%.4f, %d nodes (%.2f per decade)\n"
              "  rho: 10^%.4f .. 10^%.4f, %d nodes (%.2f per decade)\n",
              path.c_str(), t.Generator().c_str(), t.LMax(), e.x0, e.Last(), e.n, 1.0 / e.dx, r.x0, r.Last(), r.n,
              1.0 / r.dx);
  if (t.EtaRhoMax() > 0) std::printf("  domain limited to eta * rho <= %g\n", t.EtaRhoMax());
}

// Errors of one quantity, overall and on a map of (eta decade, rho decade).
struct ErrorMap {
  const char *name;
  int ex0, nex, rx0, nrx;
  std::vector<double> cell;
  std::vector<double> perL = std::vector<double>(64, 0.0);
  double worst = 0.0;
  double wEta = 0, wRho = 0;
  int wL = 0;
  long n = 0, bad = 0;
  void Init(const char *nm, const CoulombTable &t) {
    name = nm;
    ex0 = int(std::floor(t.EtaAxis().x0 + 1e-9));
    nex = int(std::ceil(t.EtaAxis().Last() - 1e-9)) - ex0;
    rx0 = int(std::floor(t.RhoAxis().x0 + 1e-9));
    nrx = int(std::ceil(t.RhoAxis().Last() - 1e-9)) - rx0;
    cell.assign(nex * nrx, 0.0);
  }
  void Add(int l, double eta, double rho, double e) {
    n++;
    if (!std::isfinite(e)) {
      bad++;
      return;
    }
    int a = std::min(nex - 1, std::max(0, int(std::floor(std::log10(eta))) - ex0));
    int b = std::min(nrx - 1, std::max(0, int(std::floor(std::log10(rho))) - rx0));
    cell[a * nrx + b] = std::max(cell[a * nrx + b], e);
    if (l < 64) perL[l] = std::max(perL[l], e);
    if (e > worst) {
      worst = e;
      wEta = eta;
      wRho = rho;
      wL = l;
    }
  }
  void Merge(const ErrorMap &o) {
    for (size_t k = 0; k < cell.size(); k++) cell[k] = std::max(cell[k], o.cell[k]);
    for (size_t k = 0; k < perL.size(); k++) perL[k] = std::max(perL[k], o.perL[k]);
    if (o.worst > worst) {
      worst = o.worst;
      wEta = o.wEta;
      wRho = o.wRho;
      wL = o.wL;
    }
    n += o.n;
    bad += o.bad;
  }
  void Print(double limit, int lMax) const {
    std::printf("\n%s: worst %.2e at l=%d eta=%.4e rho=%.4e (%ld points, %ld without a reference value)  %s\n", name,
                worst, wL, wEta, wRho, n, bad, worst < limit ? "ok" : "ABOVE LIMIT");
    std::printf("  worst per l:");
    for (int l = 0; l <= lMax && l < 64; l++) std::printf(" %d:%.1e", l, perL[l]);
    std::printf("\n");
    std::printf("  log10 eta \\ log10 rho:");
    for (int b = 0; b < nrx; b++) std::printf(" %8d", rx0 + b);
    std::printf("\n");
    for (int a = 0; a < nex; a++) {
      std::printf("  %21d", ex0 + a);
      for (int b = 0; b < nrx; b++) std::printf(" %8.1e", cell[a * nrx + b]);
      std::printf("\n");
    }
  }
};

int Validate(const CoulombTable &t, long n, unsigned seed, int threads, double etaRhoLimit, double limit) {
  const auto &E = t.EtaAxis(), &R = t.RhoAxis();
  ErrorMap mP, mS, mPhi, mB;
  mP.Init("P (relative)", t);
  mS.Init("S (absolute, relative above |S| = 1)", t);
  mPhi.Init("phi (absolute, rad)", t);
  mB.Init("bound S (absolute, relative above |S| = 1)", t);
  const long nTurn = n / 20, nEdge = n / 20;
  const long total = n + nTurn + nEdge;
  double tTable = 0, tExact = 0;
#pragma omp parallel num_threads(threads)
  {
    ErrorMap lP = mP, lS = mS, lPhi = mPhi, lB = mB;
#ifdef _OPENMP
    std::mt19937_64 rng(seed + 7919 * omp_get_thread_num());
#else
    std::mt19937_64 rng(seed);
#endif
    std::uniform_real_distribution<double> U(0.0, 1.0);
    double myTable = 0, myExact = 0;
#pragma omp for schedule(dynamic, 1000)
    for (long q = 0; q < total; q++) {
      const int l = std::min(t.LMax(), int(U(rng) * (t.LMax() + 1)));
      double x = E.x0 + U(rng) * (E.Last() - E.x0), y = R.x0 + U(rng) * (R.Last() - R.x0);
      if (q >= n && q < n + nTurn) {  // near the turning point
        const double eta = std::pow(10.0, x);
        y = std::log10(eta + std::sqrt(eta * eta + l * (l + 1.0))) + 0.1 * (U(rng) - 0.5);
        if (y < R.x0 || y > R.Last()) continue;
      } else if (q >= n + nTurn) {  // within two cells of an edge
        const double w = 2.0 * U(rng);
        switch (int(U(rng) * 4)) {
          case 0: x = E.x0 + w * E.dx; break;
          case 1: x = E.Last() - w * E.dx; break;
          case 2: y = R.x0 + w * R.dx; break;
          default: y = R.Last() - w * R.dx; break;
        }
      }
      const double eta = std::pow(10.0, x), rho = std::pow(10.0, y);
      if (!t.Inside(eta, rho) || (etaRhoLimit > 0 && eta * rho > etaRhoLimit)) continue;
      auto t0 = std::chrono::steady_clock::now();
      CoulombTable::Open to;
      double tb = 0;
      bool okO = t.LookupOpen(l, eta, rho, to);
      bool okB = t.LookupBound(l, eta, rho, tb);
      auto t1 = std::chrono::steady_clock::now();
      CoulombTable::Open eo = CoulombTable::ExactOpen(l, eta, rho);
      double eb = CoulombTable::ExactBound(l, eta, rho);
      auto t2 = std::chrono::steady_clock::now();
      myTable += std::chrono::duration<double>(t1 - t0).count();
      myExact += std::chrono::duration<double>(t2 - t1).count();
      if (!okO) {
        lP.Add(l, eta, rho, NAN);
        continue;
      }
      // P below ~1e-250 is at the edge of the double range (l = 10, eta ~ 100, rho ~ 1e-3:
      // the exact value underflows); no width can use it
      if (eo.P > 1e-250) lP.Add(l, eta, rho, std::fabs(to.P / eo.P - 1.0));
      lS.Add(l, eta, rho, std::fabs(to.S - eo.S) / std::max(1.0, std::fabs(eo.S)));
      lPhi.Add(l, eta, rho, std::fabs(Wrap(to.phi - eo.phi)));
      if (std::isfinite(eb)) lB.Add(l, eta, rho, okB ? std::fabs(tb - eb) / std::max(1.0, std::fabs(eb)) : NAN);
    }
#pragma omp critical(validate_merge)
    {
      mP.Merge(lP);
      mS.Merge(lS);
      mPhi.Merge(lPhi);
      mB.Merge(lB);
      tTable += myTable;
      tExact += myExact;
    }
  }
  std::printf("%ld points: %ld uniform in (log eta, log rho), %ld near the turning point, %ld at the edges\n", total, n,
              nTurn, nEdge);
  std::printf("time per point: table %.3f us (open + bound), exact %.3f us (x%.0f)\n", 1e6 * tTable / total,
              1e6 * tExact / total, tExact / tTable);
  mP.Print(limit, t.LMax());
  mS.Print(limit, t.LMax());
  mPhi.Print(limit, t.LMax());
  mB.Print(limit, t.LMax());

  // the references themselves: Michel vs GSL's Coulomb functions, and the bound
  // route at eta -> 0 vs the closed neutral form
  std::mt19937_64 rng(seed + 1);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  double wP = 0, wS = 0, wN = 0;
  int refused = 0;
  gsl_set_error_handler_off();
  for (int q = 0; q < 5000; q++) {
    const int l = std::min(t.LMax(), int(U(rng) * (t.LMax() + 1)));
    const double eta = std::pow(10.0, E.x0 + U(rng) * (E.Last() - E.x0));
    const double rho = std::pow(10.0, R.x0 + U(rng) * (R.Last() - R.x0));
    gsl_sf_result F, Fp, G, Gp;
    double eF, eG;
    if (gsl_sf_coulomb_wave_FG_e(eta, rho, l, 0, &F, &Fp, &G, &Gp, &eF, &eG) != GSL_SUCCESS || eF != 0 || eG != 0) {
      refused++;
    } else {
      const double d = F.val * F.val + G.val * G.val;
      CoulombTable::Open eo = CoulombTable::ExactOpen(l, eta, rho);
      if (std::isfinite(d) && d > 0) {
        wP = std::max(wP, std::fabs(rho / d / eo.P - 1.0));
        wS = std::max(wS, std::fabs(rho * (F.val * Fp.val + G.val * Gp.val) / d - eo.S) / std::max(1.0, std::fabs(eo.S)));
      }
    }
    wN = std::max(wN, std::fabs(CoulombTable::ExactBound(l, 0.0, rho) - ShftFunc::NeutralBoundShift(l, rho)) /
                          std::max(1.0, std::fabs(ShftFunc::NeutralBoundShift(l, rho))));
  }
  std::printf("\nreferences: Michel vs GSL Coulomb P %.2e, S %.2e (%d of 5000 points GSL refused or scaled); "
              "bound route at eta = 0 vs the neutral closed form %.2e\n",
              wP, wS, refused, wN);

  const bool pass = mP.worst < limit && mS.worst < limit && mPhi.worst < limit && mB.worst < limit && mP.bad == 0;
  std::printf("%s: interpolation error limit %.0e\n", pass ? "PASS" : "FAIL", limit);
  return pass ? 0 : 1;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) return Usage();
  const std::string cmd = argv[1];
  // FILE may be omitted: the file the engine reads by default
  const bool hasPath = argc >= 3 && argv[2][0] != '-';
  const std::string path = hasPath ? argv[2] : CoulombTable::DefaultPath();
  const int firstOption = hasPath ? 3 : 2;
  CoulombTable::BuildSpec spec;
  long n = 200000;
  unsigned seed = 1;
  int order = 8;
  double limit = 1e-7;
#ifdef _OPENMP
  int threads = omp_get_max_threads();
#else
  int threads = 1;
#endif
  for (int a = firstOption; a < argc; a++) {
    const std::string o = argv[a];
    auto need = [&](int k) {
      if (a + k >= argc) {
        std::fprintf(stderr, "ERROR: %s needs %d value(s)\n", o.c_str(), k);
        std::exit(2);
      }
    };
    if (o == "--lmax") { need(1); spec.lMax = std::atoi(argv[++a]); }
    else if (o == "--eta") { need(2); spec.etaMin = std::atof(argv[++a]); spec.etaMax = std::atof(argv[++a]); }
    else if (o == "--rho") { need(2); spec.rhoMin = std::atof(argv[++a]); spec.rhoMax = std::atof(argv[++a]); }
    else if (o == "--eta-rho-max") { need(1); spec.etaRhoMax = std::atof(argv[++a]); }
    else if (o == "--limit") { need(1); limit = std::atof(argv[++a]); }
    else if (o == "--order") { need(1); order = std::atoi(argv[++a]); }
    else if (o == "--per-decade") { need(1); spec.perDecade = std::atof(argv[++a]); }
    else if (o == "--threads") { need(1); threads = std::atoi(argv[++a]); }
    else if (o == "--n") { need(1); n = std::atol(argv[++a]); }
    else if (o == "--seed") { need(1); seed = unsigned(std::atol(argv[++a])); }
    else return Usage();
  }
  spec.threads = threads;
  if (cmd == "build") {
    auto t0 = std::chrono::steady_clock::now();
    std::string error;
    if (!CoulombTable::BuildFile(path, spec, error)) {
      std::fprintf(stderr, "ERROR: %s\n", error.c_str());
      return 1;
    }
    auto t1 = std::chrono::steady_clock::now();
    std::printf("built in %.1f s on %d threads\n", std::chrono::duration<double>(t1 - t0).count(), threads);
    CoulombTable t;
    if (!LoadOrSay(t, path)) return 1;
    Info(t, path);
    return 0;
  }
  CoulombTable t;
  if (!LoadOrSay(t, path)) return 1;
  if (cmd == "info") {
    Info(t, path);
    return 0;
  }
  if (cmd == "validate") {
    t.SetOrder(order);
    Info(t, path);
    std::printf("  interpolation order %d\n", t.Order());
    if (spec.etaRhoMax != CoulombTable::BuildSpec().etaRhoMax)
      std::printf("  validating only eta * rho <= %g\n", spec.etaRhoMax);
    return Validate(t, n, seed, threads, spec.etaRhoMax != CoulombTable::BuildSpec().etaRhoMax ? spec.etaRhoMax : 0.0, limit);
  }
  return Usage();
}
