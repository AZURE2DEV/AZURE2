#include "azure2/physics/transfer/Continuum.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>

#include "azure2/physics/Constants.h"
#include "azure2/physics/HybridPotential.h"
#include "azure2/physics/transfer/DistortedWaves.h"

namespace azure2 {
namespace physics {
namespace transfer {

namespace {

using cplx = std::complex<double>;
const double kPi = std::acos(-1.0);

void CheckRequest(const ContinuumRequest &q) {
  if (q.l < 0 || q.two_s < 0 || q.two_s > 2) throw std::invalid_argument("continuum: l >= 0, spin 0, 1/2 or 1");
  if (q.two_s == 0 && q.two_j != 2 * q.l) throw std::invalid_argument("continuum: j = l without spin");
  if (std::abs(q.two_j - 2 * q.l) > q.two_s || (q.two_j + q.two_s) % 2 != 0)
    throw std::invalid_argument("continuum: j does not couple l and s");
  if (!(q.step_fm > 0.0) || !(q.r_max_fm > 40.0 * q.step_fm)) throw std::invalid_argument("continuum: bad grid");
}

void CheckRequest(const OpticalPotential &p, const ContinuumRequest &q) {
  if (p.imaginary_volume.depth_MeV != 0.0 || p.imaginary_surface.depth_MeV != 0.0 ||
      p.imaginary_spin_orbit.depth_MeV != 0.0)
    throw std::invalid_argument("continuum: the potential must be real");
  CheckRequest(q);
}

double TwoMuOverHbar2(double mu_u) { return 2.0 * mu_u * kAtomicMassUnit_MeV / (kHbarC_MeV_fm * kHbarC_MeV_fm); }

}  // namespace

double NormalisationFactor(const ScatteringState &s, double mu_u, ContinuumNormalisation n) {
  switch (n) {
    case ContinuumNormalisation::AmplitudeOne: return 1.0;
    case ContinuumNormalisation::DeltaK: return std::sqrt(2.0 / kPi);
    case ContinuumNormalisation::DeltaE: return std::sqrt(TwoMuOverHbar2(mu_u) / (kPi * s.k_fm));
  }
  return 1.0;
}

RealPotential RealPotentialOf(const OpticalPotential &potential, const ContinuumRequest &q) {
  CheckRequest(potential, q);
  const int z1z2 = q.pair.z1 * q.pair.z2;
  const double two_ls = q.two_s > 0 ? TwoLDotI(q.l, q.two_j, q.two_s) : 0.0;
  RealPotential v;
  v.total_MeV = [potential, z1z2, two_ls](double r) { return PotentialAt(potential, z1z2, two_ls, r).real(); };
  v.coulomb_MeV = [potential, z1z2](double r) { return TermsAt(potential, z1z2, r).coulomb_MeV; };
  v.point_charge_MeV_fm = RadiusOf(potential, potential.coulomb_r0_fm) == 0.0 ? z1z2 * kFineStructure * kHbarC_MeV_fm : 0.0;
  return v;
}

ScatteringState SolveScatteringStateIn(const RealPotential &v, const ContinuumRequest &q, double energy) {
  CheckRequest(q);
  if (!(energy > 0.0)) throw std::invalid_argument("continuum: E must be positive");
  ScatteringRequest sq;
  sq.pair = q.pair;
  sq.two_s = 0;  // V holds the spin-orbit term of (l, j): one wave
  sq.energy_MeV = energy;
  sq.step_fm = q.step_fm;
  sq.matching_radius_fm = q.r_max_fm;
  sq.l_max = q.l;
  auto central = [&v](double r) { return std::complex<double>(v.total_MeV(r), 0.0); };
  auto none = [](double) { return std::complex<double>(0.0, 0.0); };
  const DistortedWaves w = SolveDistortedWavesInPotential(central, none, false, v.point_charge_MeV_fm, sq);
  const PartialWave &pw = WaveOf(w, q.l, 2 * q.l);
  ScatteringState s;
  s.l = q.l;
  s.two_j = q.two_j;
  s.two_s = q.two_s;
  s.energy_MeV = energy;
  s.k_fm = w.k_fm;
  s.eta = w.eta;
  s.step_fm = q.step_fm;
  s.phase_shift_rad = 0.5 * std::arg(pw.S);  // (-pi/2, pi/2]
  // Eq. (4.2) is e^{i delta} u_E: u_E = e^{-i delta} u, real up to the round-off of |S| = 1.
  const cplx phase = std::exp(cplx(0.0, -s.phase_shift_rad));
  s.u.resize(pw.u.size());
  s.potential_MeV.assign(s.u.size(), 0.0);
  s.coulomb_MeV.assign(s.u.size(), 0.0);
  for (std::size_t i = 0; i < pw.u.size(); ++i) {
    s.u[i] = (phase * pw.u[i]).real();
    if (i == 0) continue;
    s.potential_MeV[i] = v.total_MeV(i * q.step_fm);
    s.coulomb_MeV[i] = v.coulomb_MeV(i * q.step_fm);
  }
  return s;
}

ScatteringState SolveScatteringState(const OpticalPotential &potential, const ContinuumRequest &q, double energy) {
  return SolveScatteringStateIn(RealPotentialOf(potential, q), q, energy);
}

ContinuumBin MakeBin(const OpticalPotential &potential, const ContinuumRequest &q, double centre, double width,
                     int nodes, BinWeight weight, BinRule rule) {
  return MakeBinIn(RealPotentialOf(potential, q), q, centre, width, nodes, weight, rule);
}

ContinuumBin MakeBinIn(const RealPotential &potential, const ContinuumRequest &q, double centre, double width,
                       int nodes, BinWeight weight, BinRule rule) {
  CheckRequest(q);
  if (!(width > 0.0) || !(centre > 0.0)) throw std::invalid_argument("bin: centre and width must be positive");
  if (nodes < 3 || (rule == BinRule::Simpson && nodes % 2 == 0))
    throw std::invalid_argument("bin: at least 3 nodes (an odd number for Simpson's rule)");
  ContinuumBin b;
  b.l = q.l;
  b.two_j = q.two_j;
  b.two_s = q.two_s;
  b.step_fm = q.step_fm;
  b.energy_min_MeV = std::max(centre - 0.5 * width, 1e-5);  // FRESCO's lower limit
  b.energy_max_MeV = b.energy_min_MeV + width;
  const double c = TwoMuOverHbar2(q.pair.mu_u);
  const double k_min = std::sqrt(c * b.energy_min_MeV), k_max = std::sqrt(c * b.energy_max_MeV);
  const double dk = (k_max - k_min) / (nodes - 1);
  double norm = 0.0, previous = 0.0, added = 0.0;
  for (int n = 0; n < nodes; ++n) {
    const double k = k_min + n * dk;
    const ScatteringState s = SolveScatteringStateIn(potential, q, k * k / c);
    // The phase continued across the bin; u_E changes sign with each pi added.
    double delta = s.phase_shift_rad;
    if (n > 0) {
      if (delta + added < previous - 0.5 * kPi) added += kPi;
      if (delta + added > previous + 0.5 * kPi) added -= kPi;
    }
    delta += added;
    previous = delta;
    const double sign = std::fmod(std::fabs(added / kPi) + 0.5, 2.0) < 1.0 ? 1.0 : -1.0;
    double w = (weight == BinWeight::TMatrix || weight == BinWeight::TMatrixTimesK) ? std::sin(delta) : 1.0;
    if (weight == BinWeight::FlatTimesK || weight == BinWeight::TMatrixTimesK) w *= k;
    double qw = dk;
    if (rule == BinRule::Simpson) qw = dk / 3.0 * ((n == 0 || n == nodes - 1) ? 1.0 : (n % 2 == 1 ? 4.0 : 2.0));
    if (n == 0) {
      b.u.assign(s.u.size(), 0.0);
      b.potential_MeV = s.potential_MeV;
      b.coulomb_MeV = s.coulomb_MeV;
    }
    for (std::size_t i = 0; i < s.u.size(); ++i) b.u[i] += qw * w * sign * s.u[i];
    norm += qw * w * w;
    b.k_fm.push_back(k);
    b.phase_shift_rad.push_back(delta);
  }
  const double scale = std::sqrt(2.0 / (kPi * norm));
  for (double &x : b.u) x *= scale;
  return b;
}

BinListing ListBin(const ContinuumBin &b) {
  BinListing out;
  const double h = b.step_fm;
  double n2 = 0.0, r2 = 0.0, d0 = 0.0;
  for (std::size_t i = 1; i + 1 < b.u.size(); ++i) {
    const double r = i * h;
    n2 += b.u[i] * b.u[i] * h;
    r2 += r * r * b.u[i] * b.u[i] * h;
    d0 += r * b.potential_MeV[i] * b.u[i] * h;
  }
  out.norm = std::sqrt(n2);
  out.rms_fm = std::sqrt(r2) / out.norm;
  out.d0_MeV_fm32 = b.l == 0 ? std::sqrt(4.0 * kPi) * d0 : 0.0;
  return out;
}

}  // namespace transfer
}  // namespace physics
}  // namespace azure2
