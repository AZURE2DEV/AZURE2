# reference — checks against closed forms, not against ourselves

`tests/run_tests.sh` and the projects beside it pin a chi-squared against a
number this code produced earlier. That catches a *change*, which is most of
what a regression suite is for, but it can never say the number was right:
the reference and the code under test came from the same place.

Everything in this directory is the other kind of test. It computes the same
quantity a second way — from a closed form where one exists, otherwise from an
independent implementation of the documented formula — and compares.

## beam_profile_reference_test

The beam-profile experimental effect (`docs/source/user_guide/experimental_effects.rst`),
which averages a cross section over an absolute beam energy profile rather than
a Gaussian centred on each data point. It drives
`EPoint::IntegrateTargetEffect` directly, with sub-point cross sections chosen
by hand, so it needs no model, no data files and no GUI, and runs in
milliseconds.

What it checks:

* **Parsing.** The `beamprofile` token block, its components, the detector
  resolution, truncation and the detailed-balance flag.
* **The profile.** `BeamProfileWeight` against the closed-form skewed Gaussian
  over 2–5 MeV, for one component and for two with unequal weights; truncation
  keeps the centre and zeroes the tail; the support brackets the profile; the
  lab → c.m. conversion is applied once however many segments share the effect.
* **The fold, by identity.** A constant cross section must come back *exactly*,
  whatever the profile, window and detailed-balance weight are doing, because
  numerator and denominator share the kernel. Held to 1e-12, with and without
  the window, with and without detailed balance. This is the check that would
  catch a normalisation mistake or a numerator and denominator evaluated on
  different nodes.
* **The fold, against a second implementation.** A linear and then a curved
  (narrow Gaussian bump) cross section against a dense Simpson quadrature of
  `K(E) sigma(E)` written here from the documentation, evaluating sigma exactly
  rather than interpolating it. The curved case also asserts that refining the
  sub-point grid moves the engine's answer *toward* the reference, which is the
  statement that the quadrature converges.
* **The energy shift.** A segment energy shift must carry the resolution window
  with the point; checked against the reference fold with the window moved by
  hand.
* **Detailed balance.** Matches the reference, and is not a no-op.

Both halves are mutation-tested: forcing `delta = 0` in the kernel fails the
shift check and nothing else, and dropping the `dbWeight` factor fails the
detailed-balance check and nothing else.

### One number worth knowing

The profile weight agrees with the closed form only to 5.7e-10, not to machine
precision, and the test says so in a named tolerance. `include/Constants.h`
defines `pi = 3.141592650` — the true value truncated at ten digits, a relative
error of 1.14e-9 — and the weight carries `1/(omega sqrt(2 pi))`, so it inherits
half of that. It cancels in the fold itself, where numerator and denominator
share it, which is why the identity checks hold to 1e-12. Giving the constant
its full precision would let the tolerance go to 1e-12 too, at the cost of
moving every recorded chi-squared in the suite very slightly.

## thm_coulomb_term_test

The external Coulomb term of the THM entrance vertex, `ThmCoulombTerm`
(`src/ThmFunc.cpp`, contract in `include/ThmFunc.h`):

    C_l(E) = 2 eta k Int_a^inf dr O_l(k r) / O_l(k a) j_l(p r)

(Tribble et al. 2014 eq. 2.79; Mukhamedzhanov et al. 2017 eq. 27; Typel & Baur
2003 eq. A.4). For E > 0 the integral converges only conditionally (1/r times
oscillation), which is what makes it worth an independent check.

The expected values come from `thm_coulomb_term_reference.py` (mpmath, 25
digits). It shares no numerics with the engine: O_l is `mpmath.whitw`
(DLMF 33.2.7, H+ ~ W_{-i eta, l+1/2}(-2 i rho)) rather than the library Coulomb
functions, and the whole integral is taken along two rays into the complex r
plane starting at the channel radius (j_l = (h1 + h2)/2, each half on the ray
where it decays) — no real-axis quadrature and no asymptotic expansion. Its
`--check` option recomputes every value with a real-axis segment a..a+15 fm
followed by rays, and the two agree to ~1e-14 or better. Regenerating the table takes
about half an hour:

    python3 tests/reference/thm_coulomb_term_reference.py

The grid: 7Li+p (a = 4 fm, B = 2.2246 MeV), 6Li+d (4.5 fm, 1.4735 MeV),
12C+12C (6 fm, 10.27 MeV; eta up to 62 at 50 keV, deep below the barrier)
at E = -1, -0.1, 0, 0.05, 0.3, 1, 3, 6 MeV and l = 0..3, plus 17O+n, which must
return exactly zero. Tolerance 1e-5 relative + 1e-8 absolute; the observed
agreement is ~1e-13.

Every case is run with `useGSL` both off and on, which must not matter: the
engine builds O_l from its asymptotic expansion by inward integration and calls
no library Coulomb function. (Were it to use GSL's, C_l would inherit their
errors — G_l and F_l from GSL are off by up to ~2% for l >= 2, eta ~ 2,
rho ~ 1-4, checked against mpmath.)

It also checks that E <= 0 gives a real result and that C_l is continuous
through threshold (E = +-1e-4 MeV against the E = 0 limit). `-v` prints every
value and the time per call (~0.05-1 ms).

## thm_distortion_test

The THM distortion factor (`src/ThmDistortion.cpp`,
docs/source/theory/thm_implementation.rst, "Distortion factor R(E)"): the
zero-range DWBA amplitude M = <χ⁻_sF φ_sx χ⁺_aA(βr)> and its plane-wave limit.
`thm_distortion_reference.py` evaluates both independently -- mpmath's
Coulomb functions and adaptive quadrature for point Coulomb, scipy's DOP853
with mpmath matching for a complex Woods-Saxon potential -- and prints the
table hard-coded in the test (a few minutes; the 18O(3He,d) cases need
l ≈ 80 and take longest).

What it checks: |M|² and M_PW for 12C(14N,d) at 30 MeV (forward, 90° with a
Yukawa tail, a 3 fm cutoff, a Woods-Saxon in both channels) and a made-up
18O(3He,d) at 115 MeV with the Trojan horse as target (forward and 60°);
plane waves in both channels against 4π/(κ² + q²) (R ≡ 1); an optical
potential with the nuclear part off against point Coulomb; the interpolated
weight against R between the grid nodes.

## thm_fixed_angle_test

The fixed-angle THM observable (`include/ThmAngular.h`): the spin-summed
angular distribution of the exit pair for an entrance with m_l = 0 along
p_xA, in the Blatt–Biedenharn form the engine uses (GSL 3j/6j symbols,
Gauss–Legendre window averages), against `thm_fixed_angle_reference.py`,
which does the literal M-sum
F = Σ √(2l+1) ⟨s ν l 0|J ν⟩⟨s′ ν′ l′ m′|J ν⟩ x Y_l′^m′(θ, 0) with sympy's
exact Clebsch–Gordan coefficients and mpmath's spherical harmonics and
quadrature. The amplitudes x come from a toy two-level R-matrix with two
entrance l per channel spin (a 7Li+p → α+α-like case and a spin-1/2 case with
odd L), assembled as `THMMatrixFunc` does.

What it checks: the Legendre coefficients, dσ/dΩ at five angles and five
window averages (agreement 5e-16, tolerance 1e-12); the 0–180 window equals
Σ(2J+1)|x|²/4π (the angle-integrated observable); θ = 0 equals the m_l = 0
amplitude along the axis and small windows approach it as t²; for an α+α
exit θ = 0 is the Clebsch–Gordan-weighted coherent sum, not the
`entranceL=coherent` recipe; no odd L and symmetry about 90° for identical
bosons, asymmetry for the spin-1/2 case.

