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

  const double identical = IdenticalFactor(lp, sp);
  if (identical == 0.0) return;

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
          cg1 * cg2 * tMatrixElement * y * identical;
    }
  }
}

double AmplitudeMatrix::IdenticalFactor(int lp, double sp) const {
  // Two identical particles in the exit channel: the physical amplitude is
  // A(k') + eps P12 A(-k'), eps = (-1)^(2j). P12 on the exit channel-spin
  // state gives (-1)^(2j-s'), and Y_l'(pi-theta, pi) = (-1)^l' Y_l'(theta, 0),
  // so each pathway is multiplied by 1 + (-1)^(l'+s'): 2 on allowed channels,
  // 0 on forbidden ones. Applied to elastic scattering only, which is where
  // GenMatrixFunc applies the same symmetrization (its factor 4 on |f_N|^2),
  // so the two routes stay equal; a reaction into or out of an identical pair
  // keeps the unsymmetrized normalization the cross-section route uses.
  if (aa_ != ir_ || !compound_->GetPair(aa_)->IsIdentical()) return 1.0;
  const int parity = (lp + (int)std::lround(sp)) % 2;
  return (parity == 0) ? 2.0 : 0.0;
}

void AmplitudeMatrix::AddCoulomb() {
  // Coulomb scattering is diagonal in channel spin and its projection, and
  // only exists when entrance and exit pairs are the same.
  if (aa_ != ir_) return;
  for (std::size_t i = 0; i < entranceSpins_.size(); i++) {
    const double s = entranceSpins_[i];
    const complex coulombAmplitude = point_->GetCoulombAmplitude(s);
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

  const double identical = IdenticalFactor(lp, sp);
  if (identical == 0.0) return complex(0.0, 0.0);

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
          cg1 * cg2 * y * identical;
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

AmplitudeMatrix::Product AmplitudeMatrix::BuildProduct(const SpinObservable &obs) const {
  // Both channel-spin indices decomposed into the particles' projections with the same Lane-Thomas
  // coefficients A_y uses (particle 1, the light one, first):
  //   F(m3 m4; m1 m2) = sum_{s s'} (j1 m1 j2 m2|s nu)(j3 m3 j4 m4|s' nu') M(s nu -> s' nu').
  // Projections are ordered m = j, j-1, ..., -j, as in SpinObservable::TermMatrix, and a two-particle
  // index is (row of particle 1) * dim2 + (row of particle 2).
  PPair *entrance = compound_->GetPair(aa_);
  PPair *exit = compound_->GetPair(ir_);
  const double j[5] = {0.0, entrance->GetJ(1), entrance->GetJ(2), exit->GetJ(1), exit->GetJ(2)};
  int d[5];
  for (int p = 1; p <= 4; p++) d[p] = (int)std::lround(2.0 * j[p] + 1.0);

  Product P;
  P.din = d[1] * d[2];
  P.dout = d[3] * d[4];
  P.F.assign(P.dout * P.din, complex(0.0, 0.0));
  P.slots.assign(amplitudes_.size(), {});
  for (std::size_t a = 0; a < amplitudes_.size(); a++) {
    const Amplitude &am = amplitudes_[a];
    for (int r1 = 0; r1 < d[1]; r1++)
      for (int r2 = 0; r2 < d[2]; r2++) {
        const double m1 = j[1] - r1, m2 = j[2] - r2;
        if (!Same(m1 + m2, am.v)) continue;
        const double cin = ParticleCG(j[1], j[2], am.s, m1, m2, am.v);
        if (std::fabs(cin) < 1.e-14) continue;
        for (int r3 = 0; r3 < d[3]; r3++)
          for (int r4 = 0; r4 < d[4]; r4++) {
            const double m3 = j[3] - r3, m4 = j[4] - r4;
            if (!Same(m3 + m4, am.vp)) continue;
            const double cout = ParticleCG(j[3], j[4], am.sp, m3, m4, am.vp);
            if (std::fabs(cout) < 1.e-14) continue;
            const int idx = (r3 * d[4] + r4) * P.din + (r1 * d[2] + r2);
            P.F[idx] += cin * cout * am.value;
            P.slots[a].push_back(std::make_pair(idx, cin * cout));
          }
      }
  }

  // Entrance operator A = (beam term) x (target term), exit operator B = (e1 term) x (e2 term); a particle
  // without a term contributes the unit matrix.
  auto kron = [&](int pa, int pb) {
    std::vector<complex> ma(d[pa] * d[pa], complex(0.0, 0.0)), mb(d[pb] * d[pb], complex(0.0, 0.0));
    for (int r = 0; r < d[pa]; r++) ma[r * d[pa] + r] = 1.0;
    for (int r = 0; r < d[pb]; r++) mb[r * d[pb] + r] = 1.0;
    for (const SpinTerm &t : obs.terms) {
      if (t.who == pa) ma = SpinObservable::TermMatrix(t, j[pa]);
      if (t.who == pb) mb = SpinObservable::TermMatrix(t, j[pb]);
    }
    const int n = d[pa] * d[pb];
    std::vector<complex> m(n * n, complex(0.0, 0.0));
    for (int a1 = 0; a1 < d[pa]; a1++)
      for (int b1 = 0; b1 < d[pb]; b1++)
        for (int a2 = 0; a2 < d[pa]; a2++)
          for (int b2 = 0; b2 < d[pb]; b2++)
            m[(a1 * d[pb] + b1) * n + (a2 * d[pb] + b2)] = ma[a1 * d[pa] + a2] * mb[b1 * d[pb] + b2];
    return m;
  };
  P.A = kron(1, 2);
  P.B = kron(3, 4);
  return P;
}

namespace {

//! B F A for F (dout x din), A (din x din), B (dout x dout), all row-major.
std::vector<complex> Sandwich(const std::vector<complex> &B, const std::vector<complex> &F,
                              const std::vector<complex> &A, int dout, int din) {
  std::vector<complex> FA(dout * din, complex(0.0, 0.0)), BFA(dout * din, complex(0.0, 0.0));
  for (int o = 0; o < dout; o++)
    for (int i = 0; i < din; i++) {
      const complex f = F[o * din + i];
      if (f == complex(0.0, 0.0)) continue;
      for (int k = 0; k < din; k++) FA[o * din + k] += f * A[i * din + k];
    }
  for (int o = 0; o < dout; o++)
    for (int p = 0; p < dout; p++) {
      const complex b = B[o * dout + p];
      if (b == complex(0.0, 0.0)) continue;
      for (int k = 0; k < din; k++) BFA[o * din + k] += b * FA[p * din + k];
    }
  return BFA;
}

std::vector<complex> Adjoint(const std::vector<complex> &m, int n) {
  std::vector<complex> a(n * n);
  for (int r = 0; r < n; r++)
    for (int c = 0; c < n; c++) a[r * n + c] = std::conj(m[c * n + r]);
  return a;
}

}  // namespace

double AmplitudeMatrix::SpinObservableValue(const SpinObservable &obs) const {
  const Product P = BuildProduct(obs);
  // N = Tr(F A F^dagger B) = sum_{o,i} conj(F_oi) (B F A)_oi,  D = sum |F_oi|^2.
  const std::vector<complex> BFA = Sandwich(P.B, P.F, P.A, P.dout, P.din);
  complex N(0.0, 0.0);
  double D = 0.0;
  for (std::size_t x = 0; x < P.F.size(); x++) {
    N += std::conj(P.F[x]) * BFA[x];
    D += std::norm(P.F[x]);
  }
  if (D <= 0.0) return 0.0;
  return std::real(obs.Prefactor() * N) / D;
}

std::vector<complex> AmplitudeMatrix::SpinObservableBar(const SpinObservable &obs) const {
  std::vector<complex> bar(amplitudes_.size(), complex(0.0, 0.0));
  const Product P = BuildProduct(obs);
  const std::vector<complex> BFA = Sandwich(P.B, P.F, P.A, P.dout, P.din);
  const std::vector<complex> BhFAh = Sandwich(Adjoint(P.B, P.dout), P.F, Adjoint(P.A, P.din), P.dout, P.din);
  complex N(0.0, 0.0);
  double D = 0.0;
  for (std::size_t x = 0; x < P.F.size(); x++) {
    N += std::conj(P.F[x]) * BFA[x];
    D += std::norm(P.F[x]);
  }
  if (D <= 0.0) return bar;
  const complex c = obs.Prefactor();
  const double O = std::real(c * N) / D;
  // O = Re(c N)/D with N = Tr(F A F^+ B):  dN/dF* = B F A,  d(N*)/dF* = B^+ F A^+,  dD/dF* = F, so
  //   dO/dF* = [ (c B F A + c* B^+ F A^+)/2 - O F ] / D.
  // F is linear in M with real coefficients, so a slot's cotangent is the coefficient-weighted sum; the
  // factor of two is the cotangent convention AMatrixFunc uses (as in AnalyzingPowerBar).
  for (std::size_t a = 0; a < amplitudes_.size(); a++)
    for (const std::pair<int, double> &sl : P.slots[a]) {
      const int x = sl.first;
      const complex g = (0.5 * (c * BFA[x] + std::conj(c) * BhFAh[x]) - O * P.F[x]) / D;
      bar[a] += 2.0 * sl.second * g;
    }
  return bar;
}

}  // namespace Polarization
