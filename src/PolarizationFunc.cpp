#include "PolarizationFunc.h"

#include <cmath>
#include <cstdio>

#include "AChannel.h"
#include "AngCoeff.h"
#include "CNuc.h"
#include "EPoint.h"
#include "JGroup.h"
#include "PPair.h"

namespace Polarization {
namespace {

//! Half-integer-safe equality for spins and projections.
inline bool Same(double a, double b) { return std::fabs(a - b) < 1.e-6; }

}  // namespace

// ---------------------------------------------------------------------------

AmplitudeMatrix::AmplitudeMatrix(CNuc* compound, EPoint* point, int aa, int ir)
    : compound_(compound), point_(point), aa_(aa), ir_(ir),
      theta_(point->GetCMAngle() * pi / 180.0) {
  // Channel spins available in the entrance and exit pairs: |j1-j2| .. j1+j2.
  PPair* entrance = compound_->GetPair(aa_);
  PPair* exit = compound_->GetPair(ir_);
  for (double s = std::fabs(entrance->GetJ(1) - entrance->GetJ(2));
       s <= entrance->GetJ(1) + entrance->GetJ(2) + 1.e-6; s += 1.0)
    entranceSpins_.push_back(s);
  for (double s = std::fabs(exit->GetJ(1) - exit->GetJ(2));
       s <= exit->GetJ(1) + exit->GetJ(2) + 1.e-6; s += 1.0)
    exitSpins_.push_back(s);

  // One amplitude per (s, v, s', v'); v and v' run over their projections.
  for (std::size_t i = 0; i < entranceSpins_.size(); i++) {
    const double s = entranceSpins_[i];
    for (double v = -s; v <= s + 1.e-6; v += 1.0) {
      for (std::size_t j = 0; j < exitSpins_.size(); j++) {
        const double sp = exitSpins_[j];
        for (double vp = -sp; vp <= sp + 1.e-6; vp += 1.0) {
          Amplitude a = {s, v, sp, vp, complex(0.0, 0.0)};
          amplitudes_.push_back(a);
        }
      }
    }
  }
}

complex& AmplitudeMatrix::At(double s, double v, double sp, double vp) {
  for (std::size_t i = 0; i < amplitudes_.size(); i++) {
    Amplitude& a = amplitudes_[i];
    if (Same(a.s, s) && Same(a.v, v) && Same(a.sp, sp) && Same(a.vp, vp))
      return a.value;
  }
  static complex dummy(0.0, 0.0);
  dummy = complex(0.0, 0.0);
  return dummy;
}

complex AmplitudeMatrix::Get(double s, double v, double sp, double vp) const {
  for (std::size_t i = 0; i < amplitudes_.size(); i++) {
    const Amplitude& a = amplitudes_[i];
    if (Same(a.s, s) && Same(a.v, v) && Same(a.sp, sp) && Same(a.vp, vp))
      return a.value;
  }
  return complex(0.0, 0.0);
}

void AmplitudeMatrix::AddPathway(int jNum, int chNum, int chpNum,
                                 complex tMatrixElement) {
  JGroup* jgroup = compound_->GetJGroup(jNum);
  AChannel* entrance = jgroup->GetChannel(chNum);
  AChannel* exitCh = jgroup->GetChannel(chpNum);

  const double jValue = jgroup->GetJ();
  const int l = entrance->GetL();
  const int lp = exitCh->GetL();
  const double s = entrance->GetS();
  const double sp = exitCh->GetS();

  // Seyler Eq. (4): the sum over J, l, l' of
  //   sqrt(2l+1) (s l v 0|J v) (s' l' v' v-v'|J v) [bracket] Y_{l'}^{v-v'}
  // The incident wave travels along z, so the entrance orbital projection is
  // zero and the entrance Clebsch-Gordan fixes the total projection to v.
  for (double v = -s; v <= s + 1.e-6; v += 1.0) {
    const double cg1 = AngCoeff::ClebGord(s, (double)l, jValue, v, 0.0, v);
    if (std::fabs(cg1) < 1.e-12) continue;
    for (double vp = -sp; vp <= sp + 1.e-6; vp += 1.0) {
      const double mu = v - vp;                 // outgoing orbital projection
      if (std::fabs(mu) > lp + 1.e-6) continue;
      const double cg2 = AngCoeff::ClebGord(sp, (double)lp, jValue, vp, mu, v);
      if (std::fabs(cg2) < 1.e-12) continue;

      const complex y = AngCoeff::SphericalHarmonic(lp, (int)std::lround(mu), theta_);
      At(s, v, sp, vp) += complex(0.0, 1.0) * std::sqrt(2.0 * l + 1.0) *
                          cg1 * cg2 * tMatrixElement * y;
    }
  }
}

void AmplitudeMatrix::AddCoulomb(complex coulombAmplitude) {
  // Coulomb scattering is diagonal in channel spin and its projection, and
  // only exists when entrance and exit pairs are the same.
  if (aa_ != ir_) return;
  for (std::size_t i = 0; i < entranceSpins_.size(); i++) {
    const double s = entranceSpins_[i];
    for (double v = -s; v <= s + 1.e-6; v += 1.0)
      At(s, v, s, v) += -coulombAmplitude;
  }
}

double AmplitudeMatrix::UnpolarizedCrossSection() const {
  // Average over entrance projections, sum over exit ones. The (2s+1) weights
  // and the wave-number factor are supplied by the caller, which knows the
  // conventions the rest of AZURE2 works in; this returns the bare spin sum.
  double total = 0.0;
  int nEntrance = 0;
  for (std::size_t i = 0; i < entranceSpins_.size(); i++)
    nEntrance += (int)std::lround(2.0 * entranceSpins_[i] + 1.0);
  if (nEntrance == 0) return 0.0;

  for (std::size_t i = 0; i < amplitudes_.size(); i++)
    total += std::norm(amplitudes_[i].value);
  return total / nEntrance;
}

int AmplitudeMatrix::IndexOf(double s, double v, double sp, double vp) const {
  for (std::size_t i = 0; i < amplitudes_.size(); i++) {
    const Amplitude& a = amplitudes_[i];
    if (Same(a.s, s) && Same(a.v, v) && Same(a.sp, sp) && Same(a.vp, vp))
      return (int)i;
  }
  return -1;
}

std::vector<complex> AmplitudeMatrix::AnalyzingPowerBar() const {
  std::vector<complex> bar(amplitudes_.size(), complex(0.0, 0.0));

  // A_y = N/D with N = 2 sum Im[u d*] and D = sum (|u|^2 + |d|^2), where u and
  // d are the amplitudes reached from the two entrance projections. Both are
  // needed before any slot can be differentiated, so N and D come first.
  bool haveHalf = false;
  for (std::size_t i = 0; i < entranceSpins_.size(); i++)
    if (Same(entranceSpins_[i], 0.5)) haveHalf = true;
  if (!haveHalf) return bar;

  complex interference(0.0, 0.0);
  double denominator = 0.0;
  for (std::size_t j = 0; j < exitSpins_.size(); j++) {
    const double sp = exitSpins_[j];
    for (double vp = -sp; vp <= sp + 1.e-6; vp += 1.0) {
      const complex up = Get(0.5, 0.5, sp, vp);
      const complex down = Get(0.5, -0.5, sp, vp);
      interference += up * std::conj(down);
      denominator += std::norm(up) + std::norm(down);
    }
  }
  if (denominator <= 0.0) return bar;
  const double N = 2.0 * std::imag(interference);
  const double D = denominator;

  // dN/du* = i d,  dN/dd* = -i u,  dD/du* = u,  dD/dd* = d, so that
  // dA/du* = i d/D - (N/D^2) u and dA/dd* = -i u/D - (N/D^2) d. The factor of
  // two is the cotangent convention, not part of the derivative.
  const complex I(0.0, 1.0);
  for (std::size_t j = 0; j < exitSpins_.size(); j++) {
    const double sp = exitSpins_[j];
    for (double vp = -sp; vp <= sp + 1.e-6; vp += 1.0) {
      const complex up = Get(0.5, 0.5, sp, vp);
      const complex down = Get(0.5, -0.5, sp, vp);
      const int iu = IndexOf(0.5, 0.5, sp, vp);
      const int id = IndexOf(0.5, -0.5, sp, vp);
      if (iu >= 0) bar[iu] = 2.0 * (I * down / D - (N / (D * D)) * up);
      if (id >= 0) bar[id] = 2.0 * (-I * up / D - (N / (D * D)) * down);
    }
  }
  return bar;
}

complex AmplitudeMatrix::PathwayAdjoint(int jNum, int chNum, int chpNum,
                                        const std::vector<complex>& bar) const {
  JGroup* jgroup = compound_->GetJGroup(jNum);
  AChannel* entrance = jgroup->GetChannel(chNum);
  AChannel* exitCh = jgroup->GetChannel(chpNum);

  const double jValue = jgroup->GetJ();
  const int l = entrance->GetL();
  const int lp = exitCh->GetL();
  const double s = entrance->GetS();
  const double sp = exitCh->GetS();

  // The mirror of AddPathway: same loop, same coefficients, contracted against
  // the cotangents instead of multiplied by T. M is linear in T, so the
  // coefficient is the entire derivative and nothing has to be re-derived.
  complex tbar(0.0, 0.0);
  for (double v = -s; v <= s + 1.e-6; v += 1.0) {
    const double cg1 = AngCoeff::ClebGord(s, (double)l, jValue, v, 0.0, v);
    if (std::fabs(cg1) < 1.e-12) continue;
    for (double vp = -sp; vp <= sp + 1.e-6; vp += 1.0) {
      const double mu = v - vp;
      if (std::fabs(mu) > lp + 1.e-6) continue;
      const double cg2 = AngCoeff::ClebGord(sp, (double)lp, jValue, vp, mu, v);
      if (std::fabs(cg2) < 1.e-12) continue;
      const int idx = IndexOf(s, v, sp, vp);
      if (idx < 0) continue;
      const complex y = AngCoeff::SphericalHarmonic(lp, (int)std::lround(mu), theta_);
      const complex coeff = complex(0.0, 1.0) * std::sqrt(2.0 * l + 1.0) *
                            cg1 * cg2 * y;
      tbar += std::conj(coeff) * bar[idx];
    }
  }
  return tbar;
}

void AmplitudeMatrix::DumpSpinHalf() const {
  // For spin-1/2 on spin-0 the matrix should read
  //   M = g + h sigma.n  =  [[g, -i h], [i h, g]]
  // so the two flip elements must be equal and opposite in phase. If instead
  // they come out equal, the sigma.n structure is absent and A_y vanishes by
  // construction rather than by physics.
  const complex pp = Get(0.5,  0.5, 0.5,  0.5);
  const complex pm = Get(0.5, -0.5, 0.5,  0.5);   // exit +1/2 from entrance -1/2
  const complex mp = Get(0.5,  0.5, 0.5, -0.5);   // exit -1/2 from entrance +1/2
  const complex mm = Get(0.5, -0.5, 0.5, -0.5);
  std::printf("MDUMP nonflip(++)=(%.4e,%.4e) nonflip(--)=(%.4e,%.4e) "
              "flip(+-)=(%.4e,%.4e) flip(-+)=(%.4e,%.4e)\n",
              pp.real(), pp.imag(), mm.real(), mm.imag(),
              pm.real(), pm.imag(), mp.real(), mp.imag());
}

double AmplitudeMatrix::MaxSpinFlip() const {
  double m = 0.0;
  for (std::size_t i = 0; i < amplitudes_.size(); i++) {
    const Amplitude& a = amplitudes_[i];
    if (!Same(a.v, a.vp)) m = std::max(m, std::abs(a.value));
  }
  return m;
}

double AmplitudeMatrix::AnalyzingPowerAy() const {
  // For a spin-1/2 projectile on a spin-0 target the channel spin is 1/2 and
  // the projections are the projectile's own. A_y then follows from the
  // interference between the two projections:
  //
  //   A_y = 2 Im[ sum_{s'v'} M(v=+1/2) conj(M(v=-1/2)) ] / sum |M|^2
  //
  // with y along k_in x k_out, which is the Madison convention.
  bool haveHalf = false;
  for (std::size_t i = 0; i < entranceSpins_.size(); i++)
    if (Same(entranceSpins_[i], 0.5)) haveHalf = true;
  if (!haveHalf) return 0.0;

  complex interference(0.0, 0.0);
  double denominator = 0.0;
  for (std::size_t j = 0; j < exitSpins_.size(); j++) {
    const double sp = exitSpins_[j];
    for (double vp = -sp; vp <= sp + 1.e-6; vp += 1.0) {
      const complex up = Get(0.5, 0.5, sp, vp);
      const complex down = Get(0.5, -0.5, sp, vp);
      interference += up * std::conj(down);
      denominator += std::norm(up) + std::norm(down);
    }
  }
  if (denominator <= 0.0) return 0.0;
  return 2.0 * std::imag(interference) / denominator;
}

}  // namespace Polarization
