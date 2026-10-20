#include "CoulombTable.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <set>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_hyperg.h>

#include "cwfcomp_accurate.H"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

const char kMagic[8] = {'A', 'Z', 'C', 'T', 'A', 'B', '1', '\0'};
const std::int32_t kFormat = 1;
// Identifies the exact routine the values came from; a table from another
// generator is still readable, the engine decides what to do with it.
const char kGenerator[] = "Michel cwfcomp_accurate (open), GSL hyperg_U ratio (bound)";

// Weights of the first derivative at node k of the 7 equally spaced nodes
// 0..6 (unit spacing), from Lagrange interpolation: exact to degree 6.
void FDWeights(int k, double w[7]) {
  for (int m = 0; m < 7; m++) {
    // d/dx of the Lagrange basis polynomial L_m at x = k
    double denom = 1.0;
    for (int q = 0; q < 7; q++)
      if (q != m) denom *= (m - q);
    double sum = 0.0;
    for (int r = 0; r < 7; r++) {
      if (r == m) continue;
      double prod = 1.0;
      for (int q = 0; q < 7; q++)
        if (q != m && q != r) prod *= (k - q);
      sum += prod;
    }
    w[m] = sum / denom;
  }
}

// First derivative (per unit of the axis variable) of n values spaced by h.
void Derivative(const double *f, int stride, int n, double h, double *d, int dstride) {
  if (n < 7) {  // too short for the stencil: central differences
    for (int i = 0; i < n; i++) {
      int a = std::max(0, i - 1), b = std::min(n - 1, i + 1);
      d[i * dstride] = (b > a) ? (f[b * stride] - f[a * stride]) / ((b - a) * h) : 0.0;
    }
    return;
  }
  struct Weights {
    double w[7][7];
    Weights() {
      for (int k = 0; k < 7; k++) FDWeights(k, w[k]);
    }
  };
  static const Weights W;  // thread-safe initialisation
  const auto &w = W.w;
  for (int i = 0; i < n; i++) {
    int start = std::min(std::max(i - 3, 0), n - 7);
    int k = i - start;
    double s = 0.0;
    for (int m = 0; m < 7; m++) s += w[k][m] * f[(start + m) * stride];
    d[i * dstride] = s / h;
  }
}

inline double H00(double t) { return (2 * t - 3) * t * t + 1; }
inline double H10(double t) { return ((t - 2) * t + 1) * t; }
inline double H01(double t) { return (3 - 2 * t) * t * t; }
inline double H11(double t) { return (t - 1) * t * t; }

}  // namespace

std::atomic<long> CoulombTable::hits_{0}, CoulombTable::misses_{0};

namespace {
// The tables Shared() has loaded and the files it could not read, for the whole process;
// ForgetMissing() clears the second and moves the generation on, so the per-thread cache
// in Shared() looks again.
std::mutex &SharedMutex() {
  static std::mutex m;
  return m;
}
std::map<std::string, std::unique_ptr<CoulombTable>> &SharedTables() {
  static std::map<std::string, std::unique_ptr<CoulombTable>> t;
  return t;
}
std::set<std::string> &SharedFailed() {
  static std::set<std::string> f;
  return f;
}
std::atomic<int> g_sharedGeneration{0};
}  // namespace

namespace {
// Prints how much the table was used when the process ends, if it was used at all.
struct UsageReport {
  std::atomic<long> *hits, *misses;
  ~UsageReport() {
    const long h = hits->load(), m = misses->load();
    if (h + m > 0)
      std::cerr << "Coulomb table: " << h << " lookups, " << m
                << " outside the table (exact routine used)" << std::endl;
  }
};
}  // namespace

std::string CoulombTable::DefaultPath() {
  if (const char *e = std::getenv("AZURE2_COULOMB_TABLE"))
    if (*e) return e;
  const char *home = std::getenv("HOME");
  if (!home || !*home) home = std::getenv("USERPROFILE");
  return std::string(home ? home : ".") + "/.azure2/coulomb_table.bin";
}

const CoulombTable *CoulombTable::Shared(const std::string &path) {
  // CoulFunc objects are built per channel and point, so the common case -- the
  // same path as last time on this thread -- must not take the lock.
  thread_local std::string lastPath;
  thread_local const CoulombTable *last = nullptr;
  thread_local bool lastSet = false;
  thread_local int lastGeneration = -1;
  const int generation = g_sharedGeneration.load(std::memory_order_acquire);
  if (lastSet && path == lastPath && generation == lastGeneration) return last;

  std::lock_guard<std::mutex> lock(SharedMutex());
  auto &tables = SharedTables();
  auto &failed = SharedFailed();
  static UsageReport report{&hits_, &misses_};
  const std::string file = path.empty() ? DefaultPath() : path;
  const CoulombTable *result = nullptr;
  {
    auto it = tables.find(file);
    if (it != tables.end()) {
      result = it->second.get();
    } else if (!failed.count(file)) {
      auto t = std::make_unique<CoulombTable>();
      std::string error;
      if (t->Load(file, error)) {
        std::cerr << "Coulomb table: " << file << " (l = 0.." << t->LMax() << ", "
                  << 1.0 / t->EtaAxis().dx << " nodes per decade, eta*rho <= " << t->EtaRhoMax() << ")"
                  << std::endl;
        result = t.get();
        tables[file] = std::move(t);
      } else {
        failed.insert(file);
        // The table is on by default, so a missing default file is the normal
        // state of a new installation: say so once, without alarm.  A file that
        // was named, or exists and cannot be read, is a warning.
        if (path.empty() && !std::ifstream(file.c_str()).good())
          std::cerr << "Coulomb table: none at " << file
                    << "; using the exact Coulomb functions (build one with "
                       "'AZURE2 --build-coulomb-table' or, in the GUI, Configure > Build Coulomb "
                       "Table..., about 5 min on 4 threads)." << std::endl;
        else
          std::cerr << "WARNING: Coulomb table not used: " << error
                    << "; using the exact Coulomb functions. Build one with "
                       "'AZURE2 --build-coulomb-table="
                    << file << "'." << std::endl;
      }
    }
  }
  lastPath = path;
  last = result;
  lastSet = true;
  lastGeneration = generation;
  return result;
}

void CoulombTable::ForgetMissing() {
  std::lock_guard<std::mutex> lock(SharedMutex());
  SharedFailed().clear();
  g_sharedGeneration.fetch_add(1, std::memory_order_release);
}

double CoulombTable::LogGamow(double eta) {
  const double x = 2.0 * M_PI * eta;
  if (x < 1e-8) return -0.5 * x;  // x / (e^x - 1) = 1 - x/2 + ...
  if (x < 1.0) return std::log(x / std::expm1(x));
  return std::log(x) - x - std::log1p(-std::exp(-x));
}

CoulombTable::Open CoulombTable::ExactOpen(int l, double eta, double rho) {
  Coulomb_wave_functions_accurate coul(true, std::complex<double>(l, 0.), std::complex<double>(eta, 0.));
  std::complex<double> cF, cdF, cG, cdG;
  coul.compute(std::complex<double>(rho, 0.), cF, cdF, cG, cdG);
  const double F = std::real(cF), dF = std::real(cdF), G = std::real(cG), dG = std::real(cdG);
  Open o;
  double d = F * F + G * G;
  if (std::isfinite(d) && d > 0) {
    o.P = rho / d;
    o.S = rho * (F * dF + G * dG) / d;
  } else {  // the CoulFunc::PEShift overflow guard
    const double h = std::hypot(F, G);
    o.P = rho / h / h;
    o.S = rho * ((F / h) * (dF / h) + (G / h) * (dG / h));
  }
  o.phi = std::atan2(F, G);
  return o;
}

double CoulombTable::ExactBound(int l, double etaB, double rhoB) {
  const double z = 2.0 * rhoB;
  const double a = l + 1.0 + etaB, b = 2.0 * l + 2.0;
  gsl_sf_result_e10 u0, u1;
  gsl_error_handler_t *old = gsl_set_error_handler_off();
  int s0 = gsl_sf_hyperg_U_e10_e(a, b, z, &u0);
  int s1 = gsl_sf_hyperg_U_e10_e(a + 1.0, b + 1.0, z, &u1);
  gsl_set_error_handler(old);
  if (s0 != GSL_SUCCESS || s1 != GSL_SUCCESS || u0.val == 0.0) return std::nan("");
  const double ratio = (u1.val / u0.val) * std::pow(10.0, double(u1.e10 - u0.e10));
  // S = z dlnW/dz with W = exp(-z/2) z^(l+1) U(a, b, z), U' = -a U(a+1, b+1, z)
  return -rhoB + (l + 1.0) - a * z * ratio;
}

bool CoulombTable::Build(const BuildSpec &spec) {
  lMax_ = spec.lMax;
  auto axis = [&](double lo, double hi) {
    Axis a;
    a.x0 = std::log10(lo);
    const double span = std::log10(hi) - a.x0;
    a.n = std::max(2, int(std::lround(span * spec.perDecade)) + 1);
    a.dx = span / (a.n - 1);
    return a;
  };
  eta_ = axis(spec.etaMin, spec.etaMax);
  rho_ = axis(spec.rhoMin, spec.rhoMax);
  generator_ = kGenerator;
  etaRhoMax_ = spec.etaRhoMax;
  const int nx = eta_.n, ny = rho_.n;
  sheets_.assign((lMax_ + 1) * NKIND, Sheet());
  for (Sheet &s : sheets_) s.f.assign(nx * ny, 0.0);

  // rows (l, eta) are independent
  const int nrows = (lMax_ + 1) * nx;
  int threads = std::max(1, spec.threads);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
  for (int r = 0; r < nrows; r++) {
    if (spec.cancel && spec.cancel->load(std::memory_order_relaxed)) continue;
    const int l = r / nx, i = r % nx;
    const double eta = std::pow(10.0, eta_.x0 + i * eta_.dx);
    const double lg = LogGamow(eta);
    Sheet &sp = At(l, LOGP), &ss = At(l, SHIFT), &sf = At(l, PHASE), &sb = At(l, BOUND);
    double phiPrev = 0.0;
    for (int j = 0; j < ny; j++) {
      const double y = rho_.x0 + j * rho_.dx, rho = std::pow(10.0, y);
      Open o = ExactOpen(l, eta, rho);
      double phi = o.phi;
      if (j == 0) {
        // below the turning point G > 0 and phi = arg(G + iF) is the continuous value
        if (phi < 0) phi += 2 * M_PI;
      } else {
        // dphi/dln(rho) = P: Simpson over the step picks phi's branch
        const int m = 8;
        const double t0 = (y - rho_.dx) * M_LN10, t1 = y * M_LN10, ht = (t1 - t0) / m;
        double sum = 0.0;
        for (int q = 0; q <= m; q++) {
          double p = (q == m) ? o.P : ExactOpen(l, eta, std::exp(t0 + q * ht)).P;
          sum += p * ((q == 0 || q == m) ? 1.0 : (q % 2 ? 4.0 : 2.0));
        }
        const double est = phiPrev + sum * ht / 3.0;
        phi += 2 * M_PI * std::round((est - phi) / (2 * M_PI));
      }
      phiPrev = phi;
      const int k = i * ny + j;
      sp.f[k] = std::log(o.P) - lg;
      ss.f[k] = o.S;
      sf.f[k] = phi - rho;
      sb.f[k] = ExactBound(l, eta, rho) + rho;  // S_b -> -rho_b at large rho_b
    }
    if (spec.rowsDone) spec.rowsDone->fetch_add(1, std::memory_order_relaxed);
  }
  if (spec.cancel && spec.cancel->load()) {
    sheets_.clear();
    return false;
  }
  if (order_ == 4)
    for (Sheet &s : sheets_) Derivatives(s);
  return true;
}

int CoulombTable::BuildRows(const BuildSpec &spec) {
  const double span = std::log10(spec.etaMax) - std::log10(spec.etaMin);
  return (spec.lMax + 1) * std::max(2, int(std::lround(span * spec.perDecade)) + 1);
}

bool CoulombTable::BuildFile(const std::string &path, const BuildSpec &spec, std::string &error) {
  const std::string file = path.empty() ? DefaultPath() : path;
  CoulombTable t;
  if (!t.Build(spec)) {
    error = "cancelled";
    return false;
  }
  std::error_code ec;
  const std::filesystem::path parent = std::filesystem::path(file).parent_path();
  if (!parent.empty()) std::filesystem::create_directories(parent, ec);
  const std::string tmp = file + ".new";
  if (!t.Save(tmp) || std::rename(tmp.c_str(), file.c_str()) != 0) {
    std::remove(tmp.c_str());
    error = "could not write " + file;
    return false;
  }
  return true;
}

void CoulombTable::SetOrder(int order) {
  order_ = (order == 6 || order == 8) ? order : 4;
  for (Sheet &s : sheets_) {
    if (order_ == 4 && s.fx.empty()) Derivatives(s);
    if (order_ != 4) s.fx.clear(), s.fy.clear(), s.fxy.clear();
  }
}

bool CoulombTable::Inside(double eta, double rho) const {
  if (!(eta > 0.0) || !(rho > 0.0)) return false;
  if (etaRhoMax_ > 0.0 && eta * rho > etaRhoMax_) return false;
  const double fx = (std::log10(eta) - eta_.x0) / eta_.dx, fy = (std::log10(rho) - rho_.x0) / rho_.dx;
  return fx >= 0.0 && fy >= 0.0 && fx <= eta_.n - 1 && fy <= rho_.n - 1;
}

CoulombTable::Stencil CoulombTable::MakeStencil(double x, double y) const {
  Stencil st;
  const int p = order_;
  const double fx = (x - eta_.x0) / eta_.dx, fy = (y - rho_.x0) / rho_.dx;
  st.i0 = std::min(std::max(int(fx) - (p / 2 - 1), 0), eta_.n - p);
  st.j0 = std::min(std::max(int(fy) - (p / 2 - 1), 0), rho_.n - p);
  auto weights = [p](double t, double *w) {
    for (int m = 0; m < p; m++) {
      double num = 1.0, den = 1.0;
      for (int q = 0; q < p; q++)
        if (q != m) num *= (t - q), den *= (m - q);
      w[m] = num / den;
    }
  };
  weights(fx - st.i0, st.wx);
  weights(fy - st.j0, st.wy);
  return st;
}

double CoulombTable::Apply(const Sheet &s, const Stencil &st) const {
  const int p = order_, ny = rho_.n;
  double sum = 0.0;
  for (int a = 0; a < p; a++) {
    const double *row = &s.f[(st.i0 + a) * ny + st.j0];
    double r = 0.0;
    for (int b = 0; b < p; b++) r += st.wy[b] * row[b];
    sum += st.wx[a] * r;
  }
  return sum;
}

void CoulombTable::Derivatives(Sheet &s) const {
  const int nx = eta_.n, ny = rho_.n;
  s.fx.assign(nx * ny, 0.0);
  s.fy.assign(nx * ny, 0.0);
  s.fxy.assign(nx * ny, 0.0);
  for (int i = 0; i < nx; i++) Derivative(&s.f[i * ny], 1, ny, rho_.dx, &s.fy[i * ny], 1);
  for (int j = 0; j < ny; j++) {
    Derivative(&s.f[j], ny, nx, eta_.dx, &s.fx[j], ny);
    Derivative(&s.fy[j], ny, nx, eta_.dx, &s.fxy[j], ny);
  }
}

bool CoulombTable::Cell(double x, double y, int &i, int &j, double &u, double &v) const {
  const double fx = (x - eta_.x0) / eta_.dx, fy = (y - rho_.x0) / rho_.dx;
  if (!(fx >= 0.0 && fy >= 0.0 && fx <= eta_.n - 1 && fy <= rho_.n - 1)) return false;
  i = std::min(int(fx), eta_.n - 2);
  j = std::min(int(fy), rho_.n - 2);
  u = fx - i;
  v = fy - j;
  return true;
}

double CoulombTable::Eval(const Sheet &s, double x, double y) const {
  if (order_ != 4) return Apply(s, MakeStencil(x, y));
  int i, j;
  double u, v;
  Cell(x, y, i, j, u, v);
  const int ny = rho_.n;
  const int k00 = i * ny + j, k10 = k00 + ny, k01 = k00 + 1, k11 = k10 + 1;
  const double hx = eta_.dx, hy = rho_.dx;
  const double a0 = H00(u), a1 = H01(u), b0 = H00(v), b1 = H01(v);
  const double c0 = H10(u) * hx, c1 = H11(u) * hx, d0 = H10(v) * hy, d1 = H11(v) * hy;
  return s.f[k00] * a0 * b0 + s.f[k10] * a1 * b0 + s.f[k01] * a0 * b1 + s.f[k11] * a1 * b1 +
      s.fx[k00] * c0 * b0 + s.fx[k10] * c1 * b0 + s.fx[k01] * c0 * b1 + s.fx[k11] * c1 * b1 +
      s.fy[k00] * a0 * d0 + s.fy[k10] * a1 * d0 + s.fy[k01] * a0 * d1 + s.fy[k11] * a1 * d1 +
      s.fxy[k00] * c0 * d0 + s.fxy[k10] * c1 * d0 + s.fxy[k01] * c0 * d1 + s.fxy[k11] * c1 * d1;
}

bool CoulombTable::LookupOpen(int l, double eta, double rho, Open &out) const {
  if (l < 0 || l > lMax_ || !Inside(eta, rho)) {
    misses_.fetch_add(1, std::memory_order_relaxed);
    return false;
  }
  hits_.fetch_add(1, std::memory_order_relaxed);
  const double x = std::log10(eta), y = std::log10(rho);
  double lp, sv, ph;
  if (order_ == 4) {
    lp = Eval(At(l, LOGP), x, y);
    sv = Eval(At(l, SHIFT), x, y);
    ph = Eval(At(l, PHASE), x, y);
  } else {  // one set of weights for the three sheets
    const Stencil st = MakeStencil(x, y);
    lp = Apply(At(l, LOGP), st);
    sv = Apply(At(l, SHIFT), st);
    ph = Apply(At(l, PHASE), st);
  }
  out.P = std::exp(lp + LogGamow(eta));
  out.S = sv;
  out.phi = ph + rho;
  return true;
}

bool CoulombTable::LookupBound(int l, double etaB, double rhoB, double &S) const {
  if (l < 0 || l > lMax_ || !Inside(etaB, rhoB)) {
    misses_.fetch_add(1, std::memory_order_relaxed);
    return false;
  }
  hits_.fetch_add(1, std::memory_order_relaxed);
  const double x = std::log10(etaB), y = std::log10(rhoB);
  S = Eval(At(l, BOUND), x, y) - rhoB;
  return std::isfinite(S);
}

bool CoulombTable::Save(const std::string &path) const {
  std::ofstream out(path, std::ios::binary);
  if (!out) return false;
  const double one = 1.0;  // byte-order sentinel
  out.write(kMagic, 8);
  out.write(reinterpret_cast<const char *>(&kFormat), sizeof kFormat);
  out.write(reinterpret_cast<const char *>(&one), sizeof one);
  std::int32_t lmax = lMax_;
  out.write(reinterpret_cast<const char *>(&lmax), sizeof lmax);
  for (const Axis *a : {&eta_, &rho_}) {
    std::int32_t n = a->n;
    out.write(reinterpret_cast<const char *>(&a->x0), sizeof a->x0);
    out.write(reinterpret_cast<const char *>(&a->dx), sizeof a->dx);
    out.write(reinterpret_cast<const char *>(&n), sizeof n);
  }
  char gen[128] = {0};
  std::strncpy(gen, generator_.c_str(), sizeof gen - 1);
  out.write(gen, sizeof gen);
  out.write(reinterpret_cast<const char *>(&etaRhoMax_), sizeof etaRhoMax_);
  for (const Sheet &s : sheets_) out.write(reinterpret_cast<const char *>(s.f.data()), s.f.size() * sizeof(double));
  return bool(out);
}

bool CoulombTable::Load(const std::string &path, std::string &error) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    error = "cannot open " + path;
    return false;
  }
  char magic[8];
  std::int32_t format, lmax;
  double one;
  in.read(magic, 8);
  in.read(reinterpret_cast<char *>(&format), sizeof format);
  in.read(reinterpret_cast<char *>(&one), sizeof one);
  in.read(reinterpret_cast<char *>(&lmax), sizeof lmax);
  if (!in || std::memcmp(magic, kMagic, 8) != 0) {
    error = path + " is not a Coulomb table";
    return false;
  }
  if (format != kFormat) {
    error = path + ": table format " + std::to_string(format) + ", this build reads " + std::to_string(kFormat);
    return false;
  }
  if (one != 1.0) {
    error = path + ": written with another byte order";
    return false;
  }
  Axis axes[2];
  for (Axis &a : axes) {
    std::int32_t n;
    in.read(reinterpret_cast<char *>(&a.x0), sizeof a.x0);
    in.read(reinterpret_cast<char *>(&a.dx), sizeof a.dx);
    in.read(reinterpret_cast<char *>(&n), sizeof n);
    a.n = n;
  }
  char gen[128];
  in.read(gen, sizeof gen);
  gen[sizeof gen - 1] = '\0';
  in.read(reinterpret_cast<char *>(&etaRhoMax_), sizeof etaRhoMax_);
  if (!in || lmax < 0 || axes[0].n < 2 || axes[1].n < 2) {
    error = path + ": damaged header";
    return false;
  }
  lMax_ = lmax;
  eta_ = axes[0];
  rho_ = axes[1];
  generator_ = gen;
  sheets_.assign((lMax_ + 1) * NKIND, Sheet());
  for (Sheet &s : sheets_) {
    s.f.resize(eta_.n * rho_.n);
    in.read(reinterpret_cast<char *>(s.f.data()), s.f.size() * sizeof(double));
  }
  if (!in) {
    error = path + ": truncated";
    sheets_.clear();
    return false;
  }
  if (order_ == 4)
    for (Sheet &s : sheets_) Derivatives(s);
  return true;
}
