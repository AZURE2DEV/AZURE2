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

  PPair* entrance = compound_->GetPair(aa_);
  const double j1 = entrance->GetJ(1);
  const double j2 = entrance->GetJ(2);
  if (std::fabs(j1 - 0.5) > 1.e-6) return bar;

  // N and D are needed in full before any slot can be differentiated, so the
  // projectile/target decomposition is walked twice: once to accumulate them,
  // once to distribute the derivative back onto the channel-spin amplitudes.
  complex interference(0.0, 0.0);
  double denominator = 0.0;
  for (std::size_t j = 0; j < exitSpins_.size(); j++) {
    const double sp = exitSpins_[j];
    for (double vp = -sp; vp <= sp + 1.e-6; vp += 1.0) {
      for (double m2 = -j2; m2 <= j2 + 1.e-6; m2 += 1.0) {
        complex up(0.0, 0.0), down(0.0, 0.0);
        for (std::size_t i = 0; i < entranceSpins_.size(); i++) {
          const double s = entranceSpins_[i];
          const double nuUp = 0.5 + m2, nuDn = -0.5 + m2;
          if (std::fabs(nuUp) <= s + 1.e-6)
            up   += AngCoeff::ClebGord(j1, j2, s,  0.5, m2, nuUp) * Get(s, nuUp, sp, vp);
          if (std::fabs(nuDn) <= s + 1.e-6)
            down += AngCoeff::ClebGord(j1, j2, s, -0.5, m2, nuDn) * Get(s, nuDn, sp, vp);
        }
        interference += up * std::conj(down);
        denominator  += std::norm(up) + std::norm(down);
      }
    }
  }
  if (denominator <= 0.0) return bar;
  const double N = 2.0 * std::imag(interference);
  const double D = denominator;
  const complex I(0.0, 1.0);

  for (std::size_t j = 0; j < exitSpins_.size(); j++) {
    const double sp = exitSpins_[j];
    for (double vp = -sp; vp <= sp + 1.e-6; vp += 1.0) {
      for (double m2 = -j2; m2 <= j2 + 1.e-6; m2 += 1.0) {
        const double nuUp = 0.5 + m2, nuDn = -0.5 + m2;
        complex up(0.0, 0.0), down(0.0, 0.0);
        for (std::size_t i = 0; i < entranceSpins_.size(); i++) {
          const double s = entranceSpins_[i];
          if (std::fabs(nuUp) <= s + 1.e-6)
            up   += AngCoeff::ClebGord(j1, j2, s,  0.5, m2, nuUp) * Get(s, nuUp, sp, vp);
          if (std::fabs(nuDn) <= s + 1.e-6)
            down += AngCoeff::ClebGord(j1, j2, s, -0.5, m2, nuDn) * Get(s, nuDn, sp, vp);
        }
        // dA/du* and dA/dd*, with u and d the decomposed amplitudes.
        const complex dA_du =  I * down / D - (N / (D * D)) * up;
        const complex dA_dd = -I * up   / D - (N / (D * D)) * down;
        // u and d are holomorphic in M, so du*/dM* is just the (real)
        // Clebsch-Gordan coefficient. The factor of two is the cotangent
        // convention AMatrixFunc uses, not part of the derivative.
        for (std::size_t i = 0; i < entranceSpins_.size(); i++) {
          const double s = entranceSpins_[i];
          if (std::fabs(nuUp) <= s + 1.e-6) {
            const int idx = IndexOf(s, nuUp, sp, vp);
            if (idx >= 0)
              bar[idx] += 2.0 * AngCoeff::ClebGord(j1, j2, s, 0.5, m2, nuUp) * dA_du;
          }
          if (std::fabs(nuDn) <= s + 1.e-6) {
            const int idx = IndexOf(s, nuDn, sp, vp);
            if (idx >= 0)
              bar[idx] += 2.0 * AngCoeff::ClebGord(j1, j2, s, -0.5, m2, nuDn) * dA_dd;
          }
        }
      }
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
  // A_y is defined with respect to the polarization of the *projectile*, so the
  // Pauli matrix acts on the projectile spin alone and the target spin is
  // traced over. The amplitudes are held in the channel-spin basis, in which
  // projectile and target spin are coupled, so the entrance index has to be
  // decomposed before the projectile can be addressed on its own:
  //
  //   M_{out; m1 m2} = sum_s <j1 m1 j2 m2 | s, m1+m2> M_{out; s, m1+m2}
  //
  // For a spin-0 target this collapses to the channel spin being the
  // projectile's own, which is the only case that needs no decomposition -- and
  // was the only case the first implementation handled. With a spin-1/2 target
  // such as 15N the channel spins are 0 and 1, never 1/2, and looking for a
  // channel spin of 1/2 finds nothing and returns zero for what is a perfectly
  // well defined and non-zero observable.
  PPair* entrance = compound_->GetPair(aa_);
  const double j1 = entrance->GetJ(1);   // projectile
  const double j2 = entrance->GetJ(2);   // target
  // The vector analyzing power is a spin-1/2 beam observable.
  if (std::fabs(j1 - 0.5) > 1.e-6) return 0.0;

  complex interference(0.0, 0.0);
  double denominator = 0.0;
  for (std::size_t j = 0; j < exitSpins_.size(); j++) {
    const double sp = exitSpins_[j];
    for (double vp = -sp; vp <= sp + 1.e-6; vp += 1.0) {
      for (double m2 = -j2; m2 <= j2 + 1.e-6; m2 += 1.0) {
        complex up(0.0, 0.0), down(0.0, 0.0);
        for (std::size_t i = 0; i < entranceSpins_.size(); i++) {
          const double s = entranceSpins_[i];
          const double nuUp = 0.5 + m2, nuDn = -0.5 + m2;
          if (std::fabs(nuUp) <= s + 1.e-6)
            up   += AngCoeff::ClebGord(j1, j2, s,  0.5, m2, nuUp) * Get(s, nuUp, sp, vp);
          if (std::fabs(nuDn) <= s + 1.e-6)
            down += AngCoeff::ClebGord(j1, j2, s, -0.5, m2, nuDn) * Get(s, nuDn, sp, vp);
        }
        interference += up * std::conj(down);
        denominator  += std::norm(up) + std::norm(down);
      }
    }
  }
  if (denominator <= 0.0) return 0.0;
  return 2.0 * std::imag(interference) / denominator;
}

}  // namespace Polarization
