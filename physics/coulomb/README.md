# physics/coulomb: sources and provenance

The Coulomb functions F_l, F′_l, G_l, G′_l and the Coulomb phase σ_l for real η, real ρ > 0 and integer l, behind
one interface (`include/azure2/physics/coulomb/Coulomb.h`; design in `docs/upgrade/coulomb.md`). The algorithms are
those of Jakub Skowroński's pycoulomb (github.com/skowrons94/pycoulomb), copied from two of its branches; nothing
here is fetched at build time.

## Source commits

| Here | pycoulomb branch, commit | path there |
|---|---|---|
| `coulsm/`, `coul/`, `coulfg/` | `azure2-prereqs`, 7adef98ab956e74f1b8ddee117ce4e8108939552 (2026-11-07) | `libcoulomb/coulsm/`, `libcoulomb/coul/`, `libcoulomb/coulfg/`, byte for byte |
| `src/Coulomb.cpp`, `src/methods.h`, `src/method_coulsm.cpp`, `src/method_coul.cpp`, `src/method_coulfg.cpp` | same | `libcoulomb/src/coulomb.cpp`, `methods.h`, `method_*.cpp`, adapted (below) |
| `coulcc/include/coulcc.h`, `coulcc/src/` | `coulcc-cpp`, 37f0f9b4d40c7dcfaf1827890d16cd14d50b37d1 (2026-11-07) | `coulcc/include/`, `coulcc/src/`, byte for byte |
| `src/method_coulcc.cpp` | — | written here (adaptor) |
| `include/azure2/physics/coulomb/Whittaker.h`, `src/Whittaker.cpp` | — | written here (2026-11-08; `docs/upgrade/coulomb.md` §1.4) |
| `../tests/data/whittaker.csv`, `../tests/data/make_whittaker_reference.py` | — | written here (mpmath 1.4.1, 40 digits) |
| `../tests/data/coulomb.csv`, `../tests/data/make_reference.py` | `azure2-prereqs`, as above | `tests/data/coulomb.csv`, `tests/make_reference.py` (mpmath 1.4.1, 40 digits) |

Not taken: pycoulomb's GSL method and its Whittaker functions (both need GSL, which physics/ does not link;
`docs/upgrade/coulomb.md` §1.3–1.4; the Whittaker functions here are written anew from the integral representation
of U, with pycoulomb's interface), `coulomb_range`, the Python module, the Numerov hybrid model, `coulcc.f`.

## Local changes

- Namespace `azure2::physics::coulomb` instead of `pycoulomb`; no `Method::gsl`, no `not_available` status, and
  `Method::coulcc` added; `coulomb_range` left out (the engine asks for one l at a time).
- The adaptors take one l (`Status <method>(int l, double eta, double rho, CoulombWaves &)`) instead of a range.
- The third-party sources are unchanged. They are compiled with warnings off (`-w`; they are not ours to edit), COULCC without FMA
  contraction (`-ffp-contract=off`, as in pycoulomb), and on Windows with `_USE_MATH_DEFINES` (the `M_*` constants under
  `-std=c++17`, which AZURE2's `external/coul` copy had as `#ifndef` guards in its sources).

`coulsm/` is AZURE2's former `external/coul` (in it since 17fd9ded, 2025-09-05) with the namespace `coulsm`, comments
and the guards above as the only differences; its numbers are the same bits (checked when `external/coul` was
deleted, `docs/upgrade/cnuc_split.md` §6, step 3).

## Accuracy (physics_coulomb test)

Against mpmath at 40 digits on 1296 points (l = 0–8; η = 0 and 0.01–20; ρ = 10⁻³–50), largest error per region
(relative inside the outer turning point ρ_tp = η + [η² + l(l+1)]^½, relative to the amplitude beyond it), GCC 11.5,
x86-64, `-O3`:

| Method | η = 0 | η > 0, ρ < ρ_tp | η > 0, ρ ≥ ρ_tp | F′G − FG′ − 1 (η > 0) |
|---|---|---|---|---|
| `coulsm` | G 5.4e-16; F wrong for l ≥ 2 at small ρ (upward recurrence; enters P, S, φ only through F/G, ≤ 1.4e-16) | 6.9e-10 | 8.3e-10 | 4.8e-11 |
| `coul` | F 3.2e-02 (l ≥ 6, 1 ≤ ρ ≤ 2), G 3.2e-11 | 6.9e-10 | 8.3e-10 | 4.8e-11 |
| `coulfg` | 1.9e-01 | 3.9e+00, 7 points fail (JWKB inside the turning point) | 6.3e-14 | 2.4e-15 |
| `coulcc` | 4.5e-13 | 1.5e-12 | 3.0e-13 | 4.4e-16 |

Where all are accurate (ρ ≥ ρ_tp) the four agree to 8.3 × 10⁻¹⁰, the precision of `cwfcomp`; `coulfg` and `coulcc`
agree to 3.0 × 10⁻¹³. σ_l agrees with mpmath to 7.1 × 10⁻¹⁵ rad. Eight threads calling every method in different
orders give the same bits as one thread. COULCC's own limits (found when it was translated; HANDOFF 2026-11-07): for
non-integer 2λ and ρ < ½, G loses accuracy (not reachable here, λ = l), and for integer 2λ and ρ < ½ the origin series
is good to 8 × 10⁻¹¹.

## Accuracy of the Whittaker functions (physics_whittaker, whittaker_gsl tests)

W_l(ρ) = W_{−η_b, l+½}(2ρ) and ρW′/W against mpmath at 40 digits on 1162 points (l = 0–10; η_b = 0–1000;
ρ = 10⁻⁴–200 with η_b ρ ≤ 2000), GCC 11.5, x86-64, `-O3`: W 2.1 × 10⁻¹² (scaled form at η_b = 1000, where ε |ln W|
is the limit; 1.4 × 10⁻¹³ on the 1057 points inside the double range), ρW′/W 3.3 × 10⁻¹³. GSL's U as the engine
calls it (`gsl_sf_hyperg_U_e10_e`): 6.9 × 10⁻⁹, no value at 10 points. Time per call on l = 0–3, η_b = 0–20,
ρ = 0.1–20: 4.8 µs (W and ρW′/W), GSL 0.36 µs (W only).

## Licensing (facts, from pycoulomb's `LICENSING_NOTES.md`)

`coulsm/` and `coul/` (N. Michel, Comput. Phys. Commun. 176 (2007) 232; the SAMMY-side patch) carry no licence text;
both were in AZURE2 before. `coulfg/` (A. R. Barnett's COULFG, C translation © 2022 Lawrence Livermore National
Security, LLC) carries `SPDX-License-Identifier: BSD-3-Clause` and refers to a "top-level COPYRIGHT file" that is not
in pycoulomb; the full BSD-3 text should be added next to `coulfg2.c` (open, owner's decision). `coulcc/` is a
translation of COULCC (I. J. Thompson and A. R. Barnett, Comput. Phys. Commun. 36 (1985) 363), no licence text. The
owner states that pycoulomb's C++ sources may be taken into AZURE2.
