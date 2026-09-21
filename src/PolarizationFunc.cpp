#include "PolarizationFunc.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

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

/*!
 * <j1 m1 j2 m2 | s nu>: decomposes a channel spin s into its two particles,
 * particle 1 (the light one) first, as Lane and Thomas couple them.
 *
 * AZURE2_SWAP_COUPLING_ORDER=1 couples them the other way round,
 * <j2 m2 j1 m1 | s nu> = (-1)^(j1+j2-s) <j1 m1 j2 m2 | s nu>.  It exists solely
 * to test convention invariance: the swap alone must move A_y on a target with
 * spin, and the swap together with a sign flip of every reduced-width amplitude
 * whose channel spin has (-1)^(j1+j2-s) = -1 must reproduce the original A_y to
 * round-off.  It is read once per process and is not a physics option.
 */
inline double ParticleCG(double j1, double j2, double s, double m1, double m2, double nu) {
  static const bool swapOrder = [] {
    const char *value = std::getenv("AZURE2_SWAP_COUPLING_ORDER");
    return value != nullptr && *value != '\0' && std::strcmp(value, "0") != 0;
  }();
  return swapOrder ? AngCoeff::ClebGord(j2, j1, s, m2, m1, nu)
                   : AngCoeff::ClebGord(j1, j2, s, m1, m2, nu);
}

}  // namespace

// ---------------------------------------------------------------------------

AmplitudeMatrix::AmplitudeMatrix(CNuc *compound, EPoint *point, int aa, int ir) :
  compound_(compound),
  point_(point),
  aa_(aa),
  ir_(ir),
  theta_(point->GetCMAngle() * pi / 180.0) {
  // Channel spins available in the entrance and exit pairs: |j1-j2| .. j1+j2.
  PPair *entrance = compound_->GetPair(aa_);
  PPair *exit = compound_->GetPair(ir_);
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

complex &AmplitudeMatrix::At(double s, double v, double sp, double vp) {
  for (std::size_t i = 0; i < amplitudes_.size(); i++) {
    Amplitude &a = amplitudes_[i];
    if (Same(a.s, s) && Same(a.v, v) && Same(a.sp, sp) && Same(a.vp, vp))
      return a.value;
  }
  static complex dummy(0.0, 0.0);
  dummy = complex(0.0, 0.0);
  return dummy;
}

complex AmplitudeMatrix::Get(double s, double v, double sp, double vp) const {
  for (std::size_t i = 0; i < amplitudes_.size(); i++) {
    const Amplitude &a = amplitudes_[i];
    if (Same(a.s, s) && Same(a.v, v) && Same(a.sp, sp) && Same(a.vp, vp))
      return a.value;
  }
  return complex(0.0, 0.0);
}

void AmplitudeMatrix::AddPathway(int jNum, int chNum, int chpNum,
                                 complex tMatrixElement) {
  JGroup *jgroup = compound_->GetJGroup(jNum);
  AChannel *entrance = jgroup->GetChannel(chNum);
  AChannel *exitCh = jgroup->GetChannel(chpNum);

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
      const double mu = v - vp;  // outgoing orbital projection
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
    const Amplitude &a = amplitudes_[i];
    if (Same(a.s, s) && Same(a.v, v) && Same(a.sp, sp) && Same(a.vp, vp))
      return (int)i;
  }
  return -1;
}

std::vector<complex> AmplitudeMatrix::AnalyzingPowerBar() const {
  std::vector<complex> bar(amplitudes_.size(), complex(0.0, 0.0));

  PPair *entrance = compound_->GetPair(aa_);
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
            up += ParticleCG(j1, j2, s, 0.5, m2, nuUp) * Get(s, nuUp, sp, vp);
          if (std::fabs(nuDn) <= s + 1.e-6)
            down += ParticleCG(j1, j2, s, -0.5, m2, nuDn) * Get(s, nuDn, sp, vp);
        }
        interference += up * std::conj(down);
        denominator += std::norm(up) + std::norm(down);
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
            up += ParticleCG(j1, j2, s, 0.5, m2, nuUp) * Get(s, nuUp, sp, vp);
          if (std::fabs(nuDn) <= s + 1.e-6)
            down += ParticleCG(j1, j2, s, -0.5, m2, nuDn) * Get(s, nuDn, sp, vp);
        }
        // dA/du* and dA/dd*, with u and d the decomposed amplitudes.
        const complex dA_du = I * down / D - (N / (D * D)) * up;
        const complex dA_dd = -I * up / D - (N / (D * D)) * down;
        // u and d are holomorphic in M, so du*/dM* is just the (real)
        // Clebsch-Gordan coefficient. The factor of two is the cotangent
        // convention AMatrixFunc uses, not part of the derivative.
        for (std::size_t i = 0; i < entranceSpins_.size(); i++) {
          const double s = entranceSpins_[i];
          if (std::fabs(nuUp) <= s + 1.e-6) {
            const int idx = IndexOf(s, nuUp, sp, vp);
            if (idx >= 0)
              bar[idx] += 2.0 * ParticleCG(j1, j2, s, 0.5, m2, nuUp) * dA_du;
          }
          if (std::fabs(nuDn) <= s + 1.e-6) {
            const int idx = IndexOf(s, nuDn, sp, vp);
            if (idx >= 0)
              bar[idx] += 2.0 * ParticleCG(j1, j2, s, -0.5, m2, nuDn) * dA_dd;
          }
        }
      }
    }
  }
  return bar;
}

complex AmplitudeMatrix::PathwayAdjoint(int jNum, int chNum, int chpNum,
                                        const std::vector<complex> &bar) const {
  JGroup *jgroup = compound_->GetJGroup(jNum);
  AChannel *entrance = jgroup->GetChannel(chNum);
  AChannel *exitCh = jgroup->GetChannel(chpNum);

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
  const complex pp = Get(0.5, 0.5, 0.5, 0.5);
  const complex pm = Get(0.5, -0.5, 0.5, 0.5);  // exit +1/2 from entrance -1/2
  const complex mp = Get(0.5, 0.5, 0.5, -0.5);  // exit -1/2 from entrance +1/2
  const complex mm = Get(0.5, -0.5, 0.5, -0.5);
  std::printf("MDUMP nonflip(++)=(%.4e,%.4e) nonflip(--)=(%.4e,%.4e) "
              "flip(+-)=(%.4e,%.4e) flip(-+)=(%.4e,%.4e)\n",
              pp.real(), pp.imag(), mm.real(), mm.imag(),
              pm.real(), pm.imag(), mp.real(), mp.imag());
}

double AmplitudeMatrix::MaxSpinFlip() const {
  double m = 0.0;
  for (std::size_t i = 0; i < amplitudes_.size(); i++) {
    const Amplitude &a = amplitudes_[i];
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
  // The coupling order follows Lane and Thomas, on whose formalism AZURE2 is
  // built: "this channel spin s is formed by coupling I1 and I2 together:
  // s = I1 + I2", with coefficients "(I1 I2 i1 i2 | s v) ... from the
  // (I1 i1, I2 i2) scheme to the (I1 I2, s v) scheme ... as discussed by
  // Condon and Shortley" (RMP 30 (1958) 257, sec. III.2a). So particle 1 comes
  // first and the phases are Condon-Shortley, which is what AngCoeff::ClebGord
  // provides. Nothing else in AZURE2 fixes this order -- the unpolarized cross
  // section adds channel spins incoherently and is blind to it -- so it is
  // recorded here rather than left implicit.  ParticleCG supplies the
  // coefficients, and can reverse the order for a convention-invariance test.
  PPair *entrance = compound_->GetPair(aa_);
  const double j1 = entrance->GetJ(1);  // particle 1, the light one
  const double j2 = entrance->GetJ(2);  // particle 2, the heavy one
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
            up += ParticleCG(j1, j2, s, 0.5, m2, nuUp) * Get(s, nuUp, sp, vp);
          if (std::fabs(nuDn) <= s + 1.e-6)
            down += ParticleCG(j1, j2, s, -0.5, m2, nuDn) * Get(s, nuDn, sp, vp);
        }
        interference += up * std::conj(down);
        denominator += std::norm(up) + std::norm(down);
      }
    }
  }
  if (denominator <= 0.0) return 0.0;
  return 2.0 * std::imag(interference) / denominator;
}


/*!
 * Vector polarization of the outgoing particle, produced with an unpolarized
 * beam. See the header for why this is the exit-index counterpart of
 * AnalyzingPowerAy and where the sign comes from.
 */

double AmplitudeMatrix::OutgoingPolarizationPy() const {
  // The Pauli matrix acts on the ejectile alone, so the exit channel spin has
  // to be decomposed into ejectile and residual just as AnalyzingPowerAy
  // decomposes the entrance channel spin into projectile and target:
  //
  //   M_{s' , m1' m2' ; in} = sum_{s'} <j1' m1' j2' m2' | s' m1'+m2'> M_{s' v' ; in}
  //
  // with the same Lane and Thomas coupling order (particle 1 first) and
  // Condon-Shortley phases that AngCoeff::ClebGord supplies.
  PPair *exit = compound_->GetPair(ir_);
  const double j1p = exit->GetJ(1);  // the ejectile
  const double j2p = exit->GetJ(2);  // the residual nucleus
  // A vector polarization is a spin-1/2 ejectile observable.
  if (std::fabs(j1p - 0.5) > 1.e-6) return 0.0;

  complex interference(0.0, 0.0);
  double denominator = 0.0;
  for (std::size_t i = 0; i < entranceSpins_.size(); i++) {
    const double s = entranceSpins_[i];
    for (double v = -s; v <= s + 1.e-6; v += 1.0) {
      for (double m2p = -j2p; m2p <= j2p + 1.e-6; m2p += 1.0) {
        complex up(0.0, 0.0), down(0.0, 0.0);
        for (std::size_t j = 0; j < exitSpins_.size(); j++) {
          const double sp = exitSpins_[j];
          const double nuUp = 0.5 + m2p, nuDn = -0.5 + m2p;
          if (std::fabs(nuUp) <= sp + 1.e-6)
            up += ParticleCG(j1p, j2p, sp, 0.5, m2p, nuUp) * Get(s, v, sp, nuUp);
          if (std::fabs(nuDn) <= sp + 1.e-6)
            down += ParticleCG(j1p, j2p, sp, -0.5, m2p, nuDn) * Get(s, v, sp, nuDn);
        }
        interference += up * std::conj(down);
        denominator += std::norm(up) + std::norm(down);
      }
    }
  }
  if (denominator <= 0.0) return 0.0;
  // Exit-index trace: Tr(sigma_y M M+) = -2 Im(sum up conj(down)).
  return -2.0 * std::imag(interference) / denominator;
}


/*!
 * Numerator of P_y: N = -2 Im(sum over the exit decomposition of u' conj(d')).
 */

double AmplitudeMatrix::OutgoingPolarizationNumerator() const {
  PPair *exit = compound_->GetPair(ir_);
  const double j1p = exit->GetJ(1);
  const double j2p = exit->GetJ(2);
  if (std::fabs(j1p - 0.5) > 1.e-6) return 0.0;
  complex interference(0.0, 0.0);
  for (std::size_t i = 0; i < entranceSpins_.size(); i++) {
    const double s = entranceSpins_[i];
    for (double v = -s; v <= s + 1.e-6; v += 1.0) {
      for (double m2p = -j2p; m2p <= j2p + 1.e-6; m2p += 1.0) {
        complex up(0.0, 0.0), down(0.0, 0.0);
        for (std::size_t j = 0; j < exitSpins_.size(); j++) {
          const double sp = exitSpins_[j];
          const double nuUp = 0.5 + m2p, nuDn = -0.5 + m2p;
          if (std::fabs(nuUp) <= sp + 1.e-6)
            up += ParticleCG(j1p, j2p, sp, 0.5, m2p, nuUp) * Get(s, v, sp, nuUp);
          if (std::fabs(nuDn) <= sp + 1.e-6)
            down += ParticleCG(j1p, j2p, sp, -0.5, m2p, nuDn) * Get(s, v, sp, nuDn);
        }
        interference += up * std::conj(down);
      }
    }
  }
  return -2.0 * std::imag(interference);
}

/*!
 * Reverse mode for the numerator. N is bilinear in M, so the derivative is
 * exact and has no denominator: dN/du'* = -i d', dN/dd'* = +i u'.
 */

std::vector<complex> AmplitudeMatrix::OutgoingPolarizationNumeratorBar() const {
  std::vector<complex> bar(amplitudes_.size(), complex(0.0, 0.0));
  PPair *exit = compound_->GetPair(ir_);
  const double j1p = exit->GetJ(1);
  const double j2p = exit->GetJ(2);
  if (std::fabs(j1p - 0.5) > 1.e-6) return bar;
  const complex I(0.0, 1.0);
  for (std::size_t i = 0; i < entranceSpins_.size(); i++) {
    const double s = entranceSpins_[i];
    for (double v = -s; v <= s + 1.e-6; v += 1.0) {
      for (double m2p = -j2p; m2p <= j2p + 1.e-6; m2p += 1.0) {
        const double nuUp = 0.5 + m2p, nuDn = -0.5 + m2p;
        complex up(0.0, 0.0), down(0.0, 0.0);
        for (std::size_t j = 0; j < exitSpins_.size(); j++) {
          const double sp = exitSpins_[j];
          if (std::fabs(nuUp) <= sp + 1.e-6)
            up += ParticleCG(j1p, j2p, sp, 0.5, m2p, nuUp) * Get(s, v, sp, nuUp);
          if (std::fabs(nuDn) <= sp + 1.e-6)
            down += ParticleCG(j1p, j2p, sp, -0.5, m2p, nuDn) * Get(s, v, sp, nuDn);
        }
        const complex dN_du = -I * down;
        const complex dN_dd = I * up;
        for (std::size_t j = 0; j < exitSpins_.size(); j++) {
          const double sp = exitSpins_[j];
          if (std::fabs(nuUp) <= sp + 1.e-6) {
            const int idx = IndexOf(s, v, sp, nuUp);
            if (idx >= 0)
              bar[idx] += 2.0 * ParticleCG(j1p, j2p, sp, 0.5, m2p, nuUp) * dN_du;
          }
          if (std::fabs(nuDn) <= sp + 1.e-6) {
            const int idx = IndexOf(s, v, sp, nuDn);
            if (idx >= 0)
              bar[idx] += 2.0 * ParticleCG(j1p, j2p, sp, -0.5, m2p, nuDn) * dN_dd;
          }
        }
      }
    }
  }
  return bar;
}

/*!
 * Reverse mode for P_y itself (the ratio), for completeness and for anyone
 * fitting the bare polarization rather than the published product.
 */

std::vector<complex> AmplitudeMatrix::OutgoingPolarizationBar() const {
  std::vector<complex> bar(amplitudes_.size(), complex(0.0, 0.0));
  PPair *exit = compound_->GetPair(ir_);
  const double j1p = exit->GetJ(1);
  const double j2p = exit->GetJ(2);
  if (std::fabs(j1p - 0.5) > 1.e-6) return bar;

  complex interference(0.0, 0.0);
  double denominator = 0.0;
  for (std::size_t i = 0; i < entranceSpins_.size(); i++) {
    const double s = entranceSpins_[i];
    for (double v = -s; v <= s + 1.e-6; v += 1.0) {
      for (double m2p = -j2p; m2p <= j2p + 1.e-6; m2p += 1.0) {
        complex up(0.0, 0.0), down(0.0, 0.0);
        for (std::size_t j = 0; j < exitSpins_.size(); j++) {
          const double sp = exitSpins_[j];
          const double nuUp = 0.5 + m2p, nuDn = -0.5 + m2p;
          if (std::fabs(nuUp) <= sp + 1.e-6)
            up += ParticleCG(j1p, j2p, sp, 0.5, m2p, nuUp) * Get(s, v, sp, nuUp);
          if (std::fabs(nuDn) <= sp + 1.e-6)
            down += ParticleCG(j1p, j2p, sp, -0.5, m2p, nuDn) * Get(s, v, sp, nuDn);
        }
        interference += up * std::conj(down);
        denominator += std::norm(up) + std::norm(down);
      }
    }
  }
  if (denominator <= 0.0) return bar;
  const double N = -2.0 * std::imag(interference);
  const double D = denominator;
  const complex I(0.0, 1.0);

  for (std::size_t i = 0; i < entranceSpins_.size(); i++) {
    const double s = entranceSpins_[i];
    for (double v = -s; v <= s + 1.e-6; v += 1.0) {
      for (double m2p = -j2p; m2p <= j2p + 1.e-6; m2p += 1.0) {
        const double nuUp = 0.5 + m2p, nuDn = -0.5 + m2p;
        complex up(0.0, 0.0), down(0.0, 0.0);
        for (std::size_t j = 0; j < exitSpins_.size(); j++) {
          const double sp = exitSpins_[j];
          if (std::fabs(nuUp) <= sp + 1.e-6)
            up += ParticleCG(j1p, j2p, sp, 0.5, m2p, nuUp) * Get(s, v, sp, nuUp);
          if (std::fabs(nuDn) <= sp + 1.e-6)
            down += ParticleCG(j1p, j2p, sp, -0.5, m2p, nuDn) * Get(s, v, sp, nuDn);
        }
        // Sign mirrors the entrance-index case with N = -2 Im(...).
        const complex dP_du = -I * down / D - (N / (D * D)) * up;
        const complex dP_dd = I * up / D - (N / (D * D)) * down;
        for (std::size_t j = 0; j < exitSpins_.size(); j++) {
          const double sp = exitSpins_[j];
          if (std::fabs(nuUp) <= sp + 1.e-6) {
            const int idx = IndexOf(s, v, sp, nuUp);
            if (idx >= 0)
              bar[idx] += 2.0 * ParticleCG(j1p, j2p, sp, 0.5, m2p, nuUp) * dP_du;
          }
          if (std::fabs(nuDn) <= sp + 1.e-6) {
            const int idx = IndexOf(s, v, sp, nuDn);
            if (idx >= 0)
              bar[idx] += 2.0 * ParticleCG(j1p, j2p, sp, -0.5, m2p, nuDn) * dP_dd;
          }
        }
      }
    }
  }
  return bar;
}


double AmplitudeMatrix::SelfCheckNumeratorBar(double h) const {
  const std::vector<complex> bar = OutgoingPolarizationNumeratorBar();
  AmplitudeMatrix probe(*this);
  double worst = 0.0;
  for (std::size_t i = 0; i < amplitudes_.size(); i++) {
    const complex saved = amplitudes_[i].value;
    const double scale = std::max(std::abs(saved), 1.e-12) * h;
    // d/dRe
    probe.amplitudes_[i].value = saved + complex(scale, 0.0);
    const double np = probe.OutgoingPolarizationNumerator();
    probe.amplitudes_[i].value = saved - complex(scale, 0.0);
    const double nm = probe.OutgoingPolarizationNumerator();
    const double dRe = (np - nm) / (2.0 * scale);
    // d/dIm
    probe.amplitudes_[i].value = saved + complex(0.0, scale);
    const double ip = probe.OutgoingPolarizationNumerator();
    probe.amplitudes_[i].value = saved - complex(0.0, scale);
    const double im = probe.OutgoingPolarizationNumerator();
    const double dIm = (ip - im) / (2.0 * scale);
    probe.amplitudes_[i].value = saved;
    const double refR = std::max(std::fabs(std::real(bar[i])), std::fabs(dRe));
    const double refI = std::max(std::fabs(std::imag(bar[i])), std::fabs(dIm));
    if (refR > 1.e-14)
      worst = std::max(worst, std::fabs(std::real(bar[i]) - dRe) / refR);
    if (refI > 1.e-14)
      worst = std::max(worst, std::fabs(std::imag(bar[i]) - dIm) / refI);
  }
  return worst;
}

}  // namespace Polarization
