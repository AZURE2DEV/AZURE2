#include "CoulombTable.h"

#include "azure2/support/Environment.h"
#include "azure2/physics/coulomb/Whittaker.h"
#include "GSLException.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <set>

#include <gsl/gsl_errno.h>

#include "CoulFunc.h"
#include "azure2/support/FloatingPointPrecision.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

const char kMagic[8] = {'A', 'Z', 'C', 'T', 'A', 'B', '1', '\0'};
// 2: the bound sheet from physics/coulomb's whittaker_scaled, not GSL's U ratio (coulomb.md section 1.4); a format-1
// file is refused, not reused.
const std::int32_t kFormat = 2;
// Identifies the exact routine the values came from; a table from another
// generator is still readable, the engine decides what to do with it.
const char kGenerator[] = "Michel cwfcomp_accurate (open), physics whittaker_scaled (bound)";

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


namespace {
// The tables Shared() has loaded (immutable once loaded) and the files it could not read, for the whole process;
// ForgetMissing() clears the second, so Shared() looks again.
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
}  // namespace

std::string CoulombTable::DefaultPathIn(const azure2::support::Environment &environment) {
  const std::string &named = environment.Text(azure2::support::Switch::kCoulombTable);
  if (!named.empty()) return named;
  return environment.Home() + "/.azure2/coulomb_table.bin";
}

std::string CoulombTable::DefaultPath() { return DefaultPathIn(azure2::support::ProcessEnvironment()); }

const CoulombTable *CoulombTable::Shared(const std::string &path, const std::string &defaultPath) {
  std::lock_guard<std::mutex> lock(SharedMutex());
  auto &tables = SharedTables();
  auto &failed = SharedFailed();
  const std::string file = path.empty() ? defaultPath : path;
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
                       "'AZURE2 --build-coulomb-table' or, in the setup program, Run > Build the "
                       "Coulomb Table..., about 5 min on 4 threads)." << std::endl;
        else
          std::cerr << "WARNING: Coulomb table not used: " << error
                    << "; using the exact Coulomb functions. Build one with "
                       "'AZURE2 --build-coulomb-table="
                    << file << "'." << std::endl;
      }
    }
  }
  return result;
}

void CoulombTable::ForgetMissing() {
  std::lock_guard<std::mutex> lock(SharedMutex());
  SharedFailed().clear();
}

double CoulombTable::LogGamow(double eta) {
  const double x = 2.0 * M_PI * eta;
  if (x < 1e-8) return -0.5 * x;  // x / (e^x - 1) = 1 - x/2 + ...
  if (x < 1.0) return std::log(x / std::expm1(x));
  return std::log(x) - x - std::log1p(-std::exp(-x));
}

double CoulombTable::SmallPhase(int l, double eta, double rho, double logP) {
  // F_l = C_l rho^(l+1) sum_n a_n, a_0 = 1, a_1 = eta rho / (l + 1), n (n + 2l + 1) a_n = 2 eta rho a_(n-1) -
  // rho^2 a_(n-2); C_0^2 = 2 pi eta / (exp(2 pi eta) - 1), C_l = C_(l-1) sqrt(l^2 + eta^2) / (l (2l + 1)).  Below
  // the turning point the terms are positive.  phi = atan(F / G) and P = rho / (F^2 + G^2) give sin(phi) =
  // F sqrt(P / rho); in logarithms, since C_l and P leave the double range separately.
  double logC = 0.5 * LogGamow(eta);
  for (int k = 1; k <= l; k++) logC += 0.5 * std::log(k * k + eta * eta) - std::log(k * (2.0 * k + 1.0));
  double before = 1.0, last = eta * rho / (l + 1.0), sum = before + last;
  bool converged = false;
  for (int n = 2; n < 2000 && !converged; n++) {
    const double term = (2.0 * eta * rho * last - rho * rho * before) / (n * (n + 2.0 * l + 1.0));
    sum += term;
    converged = std::abs(term) < 1e-17 * std::abs(sum) && std::abs(last) < 1e-17 * std::abs(sum);
    before = last;
    last = term;
  }
  if (!converged || !(sum > 0.0) || !std::isfinite(sum)) return std::numeric_limits<double>::quiet_NaN();
  const double sine = std::exp(logC + (l + 1) * std::log(rho) + std::log(sum) + 0.5 * (logP - std::log(rho)));
  return sine < 1.0 ? std::asin(sine) : std::numeric_limits<double>::quiet_NaN();
}

CoulombTable::Open CoulombTable::ExactOpen(int l, double eta, double rho) {
  const CoulWaves w = CoulombWavesAt(l, eta, rho);
  const double F = w.F, dF = w.dF, G = w.G, dG = w.dG;
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
  // S = rho W'/W of W_l(rho) = W_{-eta_b, l+1/2}(2 rho), from physics/coulomb (docs/upgrade/coulomb.md section 1.4):
  // finite where W itself under- or overflows; NaN if the arguments are invalid.
  namespace co = azure2::physics::coulomb;
  const co::WhittakerScaled ws = co::whittaker_scaled(l, etaB, rhoB);
  return ws.status == co::Status::ok ? ws.rho_dlog : std::nan("");
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
    const azure2::support::ExtendedPrecisionScope fpPrecision;
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

int CoulombTable::DefaultBuildThreads() {
#ifdef _OPENMP
  return std::max(1, omp_get_max_threads());
#else
  return 1;
#endif
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

bool CoulombTable::LookupOpen(int l, double eta, double rho, Open &out, CoulombTableUse *use) const {
  if (l < 0 || l > lMax_ || !Inside(eta, rho)) {
    if (use) use->outside.fetch_add(1, std::memory_order_relaxed);
    return false;
  }
  if (use) use->lookups.fetch_add(1, std::memory_order_relaxed);
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
  // The sheet holds phi - rho, so the sum is good to the rounding of rho (about 1e-16 rad) and no better.  Deep
  // below the barrier phi itself is far smaller than that (phi ~ P), and the elastic amplitude exp(-2i phi) - 1
  // would be rounding noise: there the phase comes from SmallPhase (the exact routine there would cost about 100
  // lookups, and a heavy target has most of its lookups there).  phi is monotonic from phi(0) = 0, so a small
  // value is met nowhere else, and below the limit the continuous phase and atan2's are the same.
  if (std::abs(out.phi) < kSmallPhase) {
    const double small = SmallPhase(l, eta, rho, lp + LogGamow(eta));
    out.phi = std::isfinite(small) ? small : ExactOpen(l, eta, rho).phi;
  }
  return true;
}

bool CoulombTable::LookupBound(int l, double etaB, double rhoB, double &S, CoulombTableUse *use) const {
  if (l < 0 || l > lMax_ || !Inside(etaB, rhoB)) {
    if (use) use->outside.fetch_add(1, std::memory_order_relaxed);
    return false;
  }
  if (use) use->lookups.fetch_add(1, std::memory_order_relaxed);
  const double x = std::log10(etaB), y = std::log10(rhoB);
  S = Eval(At(l, BOUND), x, y) - rhoB;
  return std::isfinite(S);
}

namespace {
// 64-bit hash of a stream of 8-byte words: MurmurHash3's 64-bit mixing step per word and its finaliser (Appleby,
// public domain), so that a change of any bit of any word moves every bit of the result.
class WordHash {
 public:
  void Add(std::uint64_t word) {
    word *= 0x87c37b91114253d5ULL;
    word = (word << 31) | (word >> 33);
    word *= 0x4cf5ad432745937fULL;
    h_ ^= word;
    h_ = (h_ << 27) | (h_ >> 37);
    h_ = h_ * 5 + 0x52dce729;
    ++count_;
  }
  std::uint64_t Final() const {
    std::uint64_t h = h_ ^ count_;
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdULL;
    h ^= h >> 33;
    h *= 0xc4ceb9fe1a85ec53ULL;
    h ^= h >> 33;
    return h;
  }

 private:
  std::uint64_t h_ = 0x9e3779b97f4a7c15ULL;
  std::uint64_t count_ = 0;
};
}  // namespace

const std::string &CoulombTable::ContentIdentity() const {
  std::call_once(identity_->once, [this] {
    WordHash hash;
    for (const Sheet &s : sheets_) {
      hash.Add(s.f.size());
      for (double v : s.f) {
        std::uint64_t word;
        std::memcpy(&word, &v, sizeof word);
        hash.Add(word);
      }
    }
    char hex[32];
    std::snprintf(hex, sizeof hex, "%016llx", static_cast<unsigned long long>(hash.Final()));
    char head[256];
    std::snprintf(head, sizeof head, "format %d order %d lmax %d eta %.17g %.17g %d rho %.17g %.17g %d etarho %.17g ",
                  static_cast<int>(kFormat), order_, lMax_, eta_.x0, eta_.dx, eta_.n, rho_.x0, rho_.dx, rho_.n,
                  etaRhoMax_);
    identity_->text = std::string(head) + "generator '" + generator_ + "' values " + hex;
  });
  return identity_->text;
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
