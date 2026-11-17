#include "azure2/physics/CollisionMatrix.h"

#include "azure2/physics/LevelPairShift.h"

#include <cmath>

namespace azure2 {
namespace physics {

namespace {
typedef std::complex<double> complex;

// The shift function of channel ch between levels la and lb (E, S, S', S'' of both levels).
inline LevelPairShift<double> pair_shift(const LevelMatrixTables &t, int la, int lb, int ch) {
  auto at = [&](int l) {
    const std::size_t i = (std::size_t)l * t.nc + ch;
    return ShiftAtLevel<double>{t.E[l], t.shift[i], t.shift_d1[i], t.shift_d2[i]};
  };
  return level_pair_shift(at(la), at(lb));
}

inline bool decoupled(const LevelMatrixTables &t, int la) {
  const double *g = &t.gamma[(std::size_t)la * t.nc];
  for (int ch = 0; ch < t.nc; ch++)
    if (std::fabs(g[ch]) >= kWidthCutoff) return false;
  return true;
}
}  // namespace

// Moved from AMatrixFunc::FillLevelMatrix with its order of operations.
void fill_level_matrix(const LevelMatrixTables &t, double inEnergy, LevelMatrixSolver::matrix_c &M) {
  const int numLevels = t.n, numChannels = t.nc;
  const bool park = t.park;
  M.resize(numLevels);
  for (int r = 0; r < numLevels; ++r) M[r].assign(numLevels, complex(0.0, 0.0));
  for (int la = 0; la < numLevels; la++) {
    const double *gammaLa = &t.gamma[(std::size_t)la * numChannels];
    const double *shiftLa = &t.shift[(std::size_t)la * numChannels];
    for (int lap = 0; lap < numLevels; lap++) {
      const double *gammaLap = &t.gamma[(std::size_t)lap * numChannels];
      complex sum(0.0, 0.0);
      double parkOverlap = 0.0;
      for (int ch = 0; ch < numChannels; ch++) {
        double gammaCh = gammaLa[ch];
        double gammaChp = gammaLap[ch];

        // Early termination for effectively zero gamma values
        if (std::fabs(gammaCh) < kWidthCutoff || std::fabs(gammaChp) < kWidthCutoff) continue;

        complex loElement = t.lo[ch];
        sum += gammaCh * gammaChp * loElement;

        char radType = t.rad_type[ch];

        if ((radType == 'M' || radType == 'E') && la == lap && t.rmc) sum += complex(0.0, 1.0) * gammaCh * gammaChp;
        if (park && radType == 'P') {
          // Park (2021) Eq. (22) with B_{lambda c} = S_c(E_lambda), Eq. (26):
          // the kernel is S_c(E) - S_c(E_lambda) + i P_c, where loElement
          // is S_c(E) - B_c + i P_c.  Eq. (9) gives the off-diagonal overlap.
          sum += gammaCh * gammaChp * (t.boundary[ch] - shiftLa[ch]);
          if (la != lap) {
            const LevelPairShift<double> chord = pair_shift(t, la, lap, ch);
            parkOverlap -= gammaCh * gammaChp * chord.slope() / chord.gap;
          }
        } else if (t.brune && radType == 'P') {
          sum += gammaCh * gammaChp * t.boundary[ch];
          if (la == lap) {
            sum -= gammaCh * gammaChp * shiftLa[ch];
          } else {
            const LevelPairShift<double> chord = pair_shift(t, la, lap, ch);
            sum -= gammaCh * gammaChp * chord.chord(inEnergy) / chord.gap;
          }
        }
      }
      complex element;
      if (park) {
        // Park (2021) Eqs. (21)-(22): A^-1 = (E_lambda - E) J - gamma (S - B_lambda + iP) gamma^T,
        // with the diagonal overlap J_{lambda lambda} of Eq. (28).
        if (la == lap) parkOverlap = t.park_norm[la];
        element = (t.E[la] - inEnergy) * parkOverlap - sum;
      } else if (la == lap) {
        double resenergy = t.E[la];
        element = resenergy - inEnergy - sum;
      } else
        element = -sum;
      // A decoupled level (see the header): a unit diagonal keeps the matrix regular.  The failed factorization
      // used to drop the whole block at that point (a zero-width level on a data energy).
      if (la == lap && element == complex(0.0, 0.0) && decoupled(t, la)) element = 1.0;
      if (!std::isfinite(element.real()) || !std::isfinite(element.imag())) continue;
      M[la][lap] = element;
    }
  }
}

CollisionMatrix collision_matrix(const LevelMatrixTables &t, double E_MeV, const ChannelTablesAtE &channels) {
  CollisionMatrix out;
  const int n = out.n = t.n, nc = out.nc = t.nc;
  out.A.assign((std::size_t)n * n, complex(0.0, 0.0));
  out.U.assign((std::size_t)nc * nc, complex(0.0, 0.0));
  out.T.assign((std::size_t)nc * nc, complex(0.0, 0.0));
  LevelMatrixSolver::matrix_c M;
  fill_level_matrix(t, E_MeV, M);
  LevelMatrixSolver solver;
  solver.Decompose(M);
  if (solver.singular()) {
    out.singular = true;
    return out;
  }
  if (n > 0) {
    const LevelMatrixSolver::matrix_c &A = solver.Inverse();
    for (int l = 0; l < n; l++)
      for (int m = 0; m < n; m++) out.A[(std::size_t)l * n + m] = A[l][m];
  }
  // The widths of each channel over the levels, with the cutoff of the level matrix; one solve per exit channel.
  std::vector<double> g((std::size_t)nc * n, 0.0);
  for (int c = 0; c < nc; c++)
    for (int l = 0; l < n; l++) {
      const double gamma = t.gamma[(std::size_t)l * nc + c];
      if (std::fabs(gamma) < kWidthCutoff) continue;
      g[(std::size_t)c * n + l] = gamma;
    }
  std::vector<complex> solve(n);
  for (int cp = 0; cp < nc; cp++) {
    for (int l = 0; l < n; l++) solve[l] = g[(std::size_t)cp * n + l];
    if (n > 0) solver.Solve(solve.data());
    const ChannelPhases exitPhases{channels.coulomb_phase[cp], channels.hard_sphere_phase[cp]};
    for (int c = 0; c < nc; c++) {
      const ChannelPhases entrancePhases{channels.coulomb_phase[c], channels.hard_sphere_phase[c]};
      const complex umatrix = pathway_amplitude(channels.sqrt_penetrability[c], channels.sqrt_penetrability[cp],
                                                width_bilinear(&g[(std::size_t)c * n], solve.data(), (std::size_t)n));
      out.U[(std::size_t)c * nc + cp] = collision_element(c == cp, entrancePhases, exitPhases, umatrix);
      out.T[(std::size_t)c * nc + cp] = t_matrix_element(c == cp, entrancePhases, exitPhases, umatrix);
    }
  }
  return out;
}

}  // namespace physics
}  // namespace azure2
