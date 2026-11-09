#ifndef AZURE2_PHYSICS_ANGULARCOUPLING_H
#define AZURE2_PHYSICS_ANGULARCOUPLING_H

namespace azure2 {
namespace physics {

/*!
 * Angular-momentum coupling coefficients for physics/, which may not use GSL (docs/upgrade/RULES.md section 1,
 * repo_consistency).  All angular momenta and projections are passed doubled (two_j = 2j), as everywhere in
 * physics/, so half-integers are exact.  The phases are those of Condon and Shortley, the same as the engine's
 * AngCoeff.cpp (GSL) uses (docs/source/theory/conventions.rst, "Clebsch-Gordan and Racah coefficients").
 *
 * The sums are Racah's closed forms (A. R. Edmonds, Angular Momentum in Quantum Mechanics (1957), Eqs. 3.7.3 and
 * 6.3.7), evaluated in long double with factorials from a table.  A symbol whose arguments violate a triangle
 * condition, a projection sum or |m| <= j, or whose j + m is not an integer, is 0.  Arguments are limited to
 * 2j <= kMaxTwiceJ (j <= 40); larger arguments throw std::out_of_range.  At j <= 10 the results agree
 * with exact values to 1e-15 (physics/tests/angular_coupling_test.cpp).
 */

/// Largest 2j accepted by the functions below.
constexpr int kMaxTwiceJ = 80;

/// Wigner 3j symbol (j1 j2 j3; m1 m2 m3).
double wigner_3j(int two_j1, int two_j2, int two_j3, int two_m1, int two_m2, int two_m3);

/// Clebsch-Gordan coefficient (j1 m1 j2 m2 | J M) = (-1)^(j1 - j2 + M) sqrt(2J + 1) (j1 j2 J; m1 m2 -M).
double clebsch_gordan(int two_j1, int two_m1, int two_j2, int two_m2, int two_J, int two_M);

/// Wigner 6j symbol {j1 j2 j3; j4 j5 j6}.
double wigner_6j(int two_j1, int two_j2, int two_j3, int two_j4, int two_j5, int two_j6);

/// Racah W(a b c d; e f) = (-1)^(a + b + c + d) {a b e; d c f} (Edmonds 6.2.13); 0 when the 6j is.
double racah_w(int two_a, int two_b, int two_c, int two_d, int two_e, int two_f);

/// Wigner 9j symbol {j1 j2 j3; j4 j5 j6; j7 j8 j9}, as the sum over x of three 6j (Edmonds 6.4.3).
double wigner_9j(int two_j1, int two_j2, int two_j3, int two_j4, int two_j5, int two_j6, int two_j7, int two_j8,
                 int two_j9);

/// True when a, b, c satisfy the triangle condition |a - b| <= c <= a + b with a + b + c an integer.
bool triangle(int two_a, int two_b, int two_c);

}  // namespace physics
}  // namespace azure2

#endif
