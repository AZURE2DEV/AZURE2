#ifndef AZURE2_PHYSICS_SPINOBSERVABLE_H
#define AZURE2_PHYSICS_SPINOBSERVABLE_H

#include <complex>
#include <string>
#include <vector>

namespace azure2 {
namespace physics {

/*!
 * A spin observable of a two-body reaction, named by an operator token.
 *
 * Grammar (one token without white space):
 *
 *   token := side [ '>' side ]
 *   side  := term { '*' term }
 *   term  := who ':' op
 *   who   := b | t | e1 | e2            beam, target, exit particle 1 (light) or 2 (heavy)
 *   op    := x | y | z                  vector: S_i / j (Pauli matrices for spin 1/2, Ohlsen's S_i for spin 1)
 *          | xx | yy | zz | xy | xz | yz  Cartesian tensor, spin 1 only: (3/2)(S_i S_j + S_j S_i) - 2 delta_ij
 *          | T<k><q> | iT<k><q>         spherical tensor tau_kq (Madison normalization), q >= 0; "i" multiplies by i
 *
 * Entrance terms (b, t) form the operator on the incident spins, exit terms (e1, e2) the operator on the outgoing
 * ones.  Without '>' a side made only of exit terms is an outgoing polarization; one made of entrance terms is an
 * analyzing power or a spin correlation.  With '>' the left side is the entrance and the right side the exit
 * (polarization transfer).  Examples: b:y (A_y), t:y (A_0y), b:y*t:y (C_yy), b:yy, b:iT11, b:T20, e1:y, b:y>e1:y.
 *
 * The value is  Re( c Tr(F A F^dagger B) ) / Tr(F F^dagger),  with F the amplitude in the product basis of the
 * particles' spin projections, A the entrance operator, B the exit operator and c = i for an iT_kq term (1
 * otherwise).  Frame: Madison (z along k_in, y along k_in x k_out); the exit operators use the same y.
 */
struct SpinTerm {
  int who = 0;          ///< 1 beam, 2 target, 3 exit particle 1, 4 exit particle 2
  int kind = 0;         ///< 1 vector, 2 Cartesian tensor, 3 spherical
  int a = 0, b = 0;     ///< Cartesian axes (0 x, 1 y, 2 z); b unused for a vector
  int k = 0, q = 0;     ///< spherical rank and projection
  bool timesI = false;  ///< "iT_kq"
};

struct SpinObservable {
  typedef std::complex<double> complex;

  std::string token;
  std::vector<SpinTerm> terms;

  /// Parses `token`; on failure returns false and fills `error`.
  bool Parse(const std::string &token, std::string &error);

  /// Checks every term against the particle spins (entrance j1, j2; exit j3, j4); on failure returns false and
  /// fills `error`.
  bool Check(double j1, double j2, double j3, double j4, std::string &error) const;

  /// The complex prefactor c of the trace (i for an iT_kq term, else 1).
  complex Prefactor() const;

  /// The (2j+1)x(2j+1) matrix of one term, row-major, rows and columns ordered m = j, j-1, ..., -j.
  static std::vector<complex> TermMatrix(const SpinTerm &t, double j);

  /// The operator on two particles `first` and `second` (who numbers; spins jFirst, jSecond): the Kronecker product
  /// of their terms, the unit matrix for a particle without one; a two-particle index is (row of the first) x
  /// (2 jSecond + 1) + (row of the second).  PairOperator(1, 2, j1, j2) is A, PairOperator(3, 4, j3, j4) is B.
  std::vector<complex> PairOperator(int first, int second, double jFirst, double jSecond) const;
};

/*!
 * Re(c Tr(F A F^dagger B)) / Tr(F F^dagger) for the amplitude F (dout x din, row-major: exit index, entrance index;
 * din = (2 j1 + 1)(2 j2 + 1), dout = (2 j3 + 1)(2 j4 + 1)), 0 when F vanishes.  With `gradient` not null it also
 * gives dO/dF* per element of F (Wirtinger: [(c B F A + c* B^dagger F A^dagger)/2 - O F] / Tr(F F^dagger)).
 */
double SpinObservableValue(const SpinObservable &observable, double j1, double j2, double j3, double j4,
                           const std::vector<std::complex<double>> &F,
                           std::vector<std::complex<double>> *gradient = nullptr);

}  // namespace physics
}  // namespace azure2

#endif
