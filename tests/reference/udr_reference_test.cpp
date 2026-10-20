// Reference checks for the user-defined numerical resolution function.
//
// The kernel is built from a temporary UDR file through the same TargetEffect
// parser and loader the engine uses, and checked against closed forms:
//   * a tabulated Gaussian delay distribution comes back with unit integral
//     and the right second moment;
//   * centring moves the centroid to zero and leaves the shape alone;
//   * the burst and channel widths add in quadrature to the variance
//     (Gaussian: sigma^2, rectangle of width c: c^2/12);
//   * the energy interpolation is linear between the tabulated energies;
//   * the time-of-flight conversions invert each other.
// It needs no model and runs in milliseconds.

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "TargetEffect.h"
#include "Config.h"

// Each consumer of the engine defines this itself (src/AZURE2.cpp for the
// CLI, api/src/pybind_api.cpp for the Python API).
Config *g_config = nullptr;

static int failures = 0;

static void check(const char *what, double got, double want, double tol) {
  bool ok = std::fabs(got - want) <= tol;
  std::printf("%-60s %s  (%.8g vs %.8g)\n", what, ok ? "PASS" : "FAIL", got, want);
  if (!ok) failures++;
}

static void moments(const UdrKernel &k, double &norm, double &mean, double &var) {
  norm = mean = 0.0;
  double second = 0.0;
  for (size_t j = 0; j + 1 < k.tau.size(); j++) {
    double dt = k.tau[j + 1] - k.tau[j];
    double a = k.tau[j], b = k.tau[j + 1], ra = k.r[j], rb = k.r[j + 1];
    norm += 0.5 * (ra + rb) * dt;
    mean += dt * (ra * (2 * a + b) + rb * (a + 2 * b)) / 6.0;
    second += dt * (ra * (3 * a * a + 2 * a * b + b * b) + rb * (a * a + 2 * a * b + 3 * b * b)) / 12.0;
  }
  mean /= norm;
  var = second / norm - mean * mean;
}

// A Gaussian delay distribution of width s (us) centred at mu, with an
// amplitude that is deliberately not normalised, at the given energies.
static std::string writeUdr(double s, double mu, double amp, const std::vector<std::pair<double, double> > &energyScale) {
  std::string path = "udr_reference_test_tmp.txt";
  std::ofstream out(path.c_str());
  out << "# synthetic\n---------------------------------------------------\n";
  for (size_t e = 0; e < energyScale.size(); e++) {
    out << energyScale[e].first << "   # eV\n";
    double w = s * energyScale[e].second;
    for (int i = 0; i <= 800; i++) {
      double t = mu + (-6.0 + 12.0 * i / 800.0) * w;
      out << t << " " << amp * std::exp(-0.5 * (t - mu) * (t - mu) / (w * w)) << "\n";
    }
    out << "\n";
  }
  return path;
}

static TargetEffect makeEffect(const Config &config, const std::string &file, double burst, double channel, int centred) {
  std::ostringstream line;
  line << "1 \"1\" 200 0 0 0 0 \"\" 0 0 0 0 \"\" 0 0 0.04 20 50 udr \"" << file << "\" 100 "
       << burst << " " << channel << " " << centred;
  std::istringstream in(line.str());
  return TargetEffect(in, config);
}

int main() {
  std::ostringstream sink;
  Config config(sink);
  const double s = 2.0e-3;  // 2 ns
  const double mn = 939.56542;
  std::vector<std::pair<double, double> > one(1, std::make_pair(1.0e6, 1.0));
  std::string file = writeUdr(s, 0.0, 7.0, one);

  // Plain, uncentred, no widths: unit integral, Gaussian moments, parser fields.
  {
    TargetEffect te = makeEffect(config, file, 0.0, 0.0, 0);
    check("parser: udr flag", te.IsUdr() ? 1 : 0, 1, 0);
    check("parser: flight path", te.GetUdrFlightPath(), 100.0, 0);
    check("parser: centred flag off", te.IsUdrCentred() ? 1 : 0, 0, 0);
    std::ostringstream log;
    check("loader: file read", te.LoadUdrTable(0.5, 2.0, log) ? 1 : 0, 1, 0);
    std::shared_ptr<const UdrKernel> k = te.BuildUdrKernel(1.0, mn);
    double norm, mean, var;
    moments(*k, norm, mean, var);
    check("kernel integral is one", norm, 1.0, 1e-9);
    check("kernel mean (uncentred Gaussian at 0)", mean, 0.0, 1e-9);
    // The tabulation is truncated at +-6 s and sampled at s/67, so the
    // trapezoid moments of the piecewise-linear density are off by ~4e-5.
    check("kernel variance = s^2", var, s * s, 1e-4 * s * s);
    check("support lower edge ~ -6 s", k->tauMin, -6.0 * s, 0.1 * s);
    check("support upper edge ~ +6 s", k->tauMax, 6.0 * s, 0.1 * s);
  }
  // Shifted Gaussian: centring removes the mean, keeps the variance; without
  // centring the mean survives.
  {
    std::string shifted = writeUdr(s, 5.0e-3, 1.0, one);
    TargetEffect c = makeEffect(config, shifted, 0.0, 0.0, 1);
    std::ostringstream log;
    c.LoadUdrTable(0.5, 2.0, log);
    std::shared_ptr<const UdrKernel> k = c.BuildUdrKernel(1.0, mn);
    double norm, mean, var;
    moments(*k, norm, mean, var);
    check("centred: mean is zero", mean, 0.0, 1e-9);
    check("centred: variance unchanged", var, s * s, 1e-4 * s * s);
    TargetEffect u = makeEffect(config, shifted, 0.0, 0.0, 0);
    u.LoadUdrTable(0.5, 2.0, log);
    k = u.BuildUdrKernel(1.0, mn);
    moments(*k, norm, mean, var);
    check("uncentred: mean is the tabulated shift", mean, 5.0e-3, 1e-9);
  }
  // Burst and channel width add in quadrature.
  {
    const double burstNs = 3.0, channelNs = 4.0;
    TargetEffect b = makeEffect(config, file, burstNs, 0.0, 1);
    std::ostringstream log;
    b.LoadUdrTable(0.5, 2.0, log);
    double norm, mean, var;
    moments(*b.BuildUdrKernel(1.0, mn), norm, mean, var);
    double sb = burstNs * 1e-3 / 2.3548200450309493;
    check("burst: integral still one", norm, 1.0, 1e-9);
    check("burst: variance = s^2 + sigma_b^2", var, s * s + sb * sb, 2e-3 * (s * s + sb * sb));
    TargetEffect c = makeEffect(config, file, 0.0, channelNs, 1);
    c.LoadUdrTable(0.5, 2.0, log);
    moments(*c.BuildUdrKernel(1.0, mn), norm, mean, var);
    double vc = channelNs * 1e-3 * channelNs * 1e-3 / 12.0;
    check("channel: variance = s^2 + c^2/12", var, s * s + vc, 2e-3 * (s * s + vc));
    TargetEffect bc = makeEffect(config, file, burstNs, channelNs, 1);
    bc.LoadUdrTable(0.5, 2.0, log);
    moments(*bc.BuildUdrKernel(1.0, mn), norm, mean, var);
    check("burst+channel: variance adds in quadrature", var, s * s + sb * sb + vc, 2e-3 * (s * s + sb * sb + vc));
  }
  // Energy interpolation: widths 1 s at 1 MeV and 3 s at 3 MeV give 2 s at
  // 2 MeV for the linearly interpolated tabulation (the same tau grid is
  // interpolated too, so the Gaussian stays a Gaussian).
  {
    std::vector<std::pair<double, double> > two;
    two.push_back(std::make_pair(1.0e6, 1.0));
    two.push_back(std::make_pair(3.0e6, 3.0));
    std::string f2 = writeUdr(s, 0.0, 1.0, two);
    TargetEffect te = makeEffect(config, f2, 0.0, 0.0, 1);
    std::ostringstream log;
    te.LoadUdrTable(1.0, 3.0, log);
    double norm, mean, var;
    moments(*te.BuildUdrKernel(2.0, mn), norm, mean, var);
    check("interpolated width at the midpoint energy", std::sqrt(var), 2.0 * s, 1e-4 * s);
    moments(*te.BuildUdrKernel(10.0, mn), norm, mean, var);
    check("above the table: last tabulation used", std::sqrt(var), 3.0 * s, 1e-4 * s);
  }
  // Delay cap: a tabulation reaching delays close to the flight time (true
  // energy above 100 x nominal) is cut there and renormalised, so a badly
  // scaled or uncentred table cannot send the sub-point window to infinity.
  {
    const double t0 = UdrKinematics::TimeOfFlight(1.0, mn, 100.0);
    const double cap = t0 - UdrKinematics::TimeOfFlight(100.0, mn, 100.0);
    std::string late = writeUdr(0.1, cap + 0.15, 1.0, one);  // Gaussian, 0.1 us wide, centred 0.15 us past the cap
    TargetEffect te = makeEffect(config, late, 0.0, 0.0, 0);
    std::ostringstream log;
    te.LoadUdrTable(0.5, 2.0, log);
    std::shared_ptr<const UdrKernel> k = te.BuildUdrKernel(1.0, mn);
    double norm, mean, var;
    moments(*k, norm, mean, var);
    check("delay cap: support ends at the cap", k->tauMax, cap, 3e-3);
    check("delay cap: integral still one", norm, 1.0, 1e-9);
    check("delay cap: mean below the cap", mean < cap ? 1 : 0, 1, 0);
  }
  // Kinematics: a 1 MeV neutron over 100 m takes 7.2356 us (7.2296 us
  // non-relativistically; beta = 0.046101), and the inverse conversion
  // recovers the energy.
  {
    double t = UdrKinematics::TimeOfFlight(1.0, mn, 100.0);
    check("1 MeV neutron, 100 m: 7.2356 us", t, 7.2356, 2e-4);
    check("energy from time inverts", UdrKinematics::EnergyFromTime(t, mn, 100.0), 1.0, 1e-12);
  }
  std::remove("udr_reference_test_tmp.txt");
  std::printf("%s: %d failure(s)\n", failures ? "FAILED" : "OK", failures);
  return failures ? 1 : 0;
}
