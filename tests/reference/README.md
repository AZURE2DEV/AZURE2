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
