# thm_fixed_angle — the fixed-angle (differential) HOES observable

`check.sh` runs copies of `tests/7Li_p_a` (7Li+p → α+α: two entrance l per
channel spin, identical-boson exit, 30 keV folding), `tests/17O`
(17O(n,α)14C: Jπ groups of both parities), `tests/6Li_d` (6Li(d,α)α) and
`tests/18O_p_a_thm` (one 1/2+ group; two segments, the 2H(18O,α15N)n
kinematics) with `experiment[<name>] ... theta=...` in the `<thm>` block
(docs/source/theory/thm_implementation.rst, "Fixed-angle observable"). With a
window the model is dσ/dΩ of the exit pair averaged over it, the
Blatt–Biedenharn sum over the HOES partial amplitudes of all Jπ groups.

| case | expectation |
|---|---|
| (a) `theta=all` | output files byte-identical to the same line without the key; model and chi2 byte-identical to no `<thm>` block |
| (b) `theta=0-180` | 4π × model = angle-integrated model (measured 4e-11, the output file's 11 digits; tol 1.2e-10) and the same chi2, for 7Li(p,α) and 17O(n,α); 18O(p,α) (isotropic 1/2+) the same for `theta=50-70` together with `ps=hulthen:0-40` and `background=linear` |
| (c) θ → 0 | windows 0–0.5° and 0–1° approach `theta=0-0` as t² (ratio of the differences 3.998–4.003); the 0° chi2 (1747.26) is not the `entranceL=coherent` one (3196.77) |
| (e) identical exit | 30–60° == 120–150° for 7Li(p,α)α, 10–40° == 140–170° for 6Li(d,α)α (to the last digit); 17O(n,α)14C differs by up to 27 % between 20–60° and 120–160° (odd L) |
| (f) refusals | a word, a single number, a reversed, negative or beyond-180° window, an open end, `theta` twice, a window with `entranceL=coherent` |
| (size) | Tumino et al. EPJA 27 (2006) 243, θ_cm = 50–70°, parameters not refitted: shape 14 % rms, 45 % max against the angle-integrated model; chi2 2138.48 → 2349.86 |

The formula itself is checked against an independent M-sum (sympy
Clebsch–Gordan, mpmath quadrature) in `tests/reference` (ctest
`thm_fixed_angle`), the Python side in `tests/pyazr/thm_fixed_angle_test.py`.
