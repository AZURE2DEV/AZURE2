# shift_fn — energy-shift functions (`shiftfn`)

A segment line may end with

    shiftfn "<expression in x>" N v0 s0 vary0 ... v(N-1) s(N-1) vary(N-1)

which adds f(x; a0..a(N-1)) to the segment's energy shift,
dE(x) = a + f(x), where x is the point's original lab energy in MeV and a the
ordinary constant shift. Each aj is a fit parameter with a starting value (the
nominal value of its penalty), a penalty width and a vary flag; its name in
`param.sav` is `segment_<key>_energy_shift_fn<j>`. The expression is read by the
same Equation class as `<targetInt>` equations. On a shift group (`shiftgroup
<label>`) the block goes on the group's first segment and the function applies to
every member, with one set of parameters. Introduced 2026-10 (design note
`11B+a/10-13-26_sqrtshift/shift_function_design.html`): different calibration
errors give different forms (magnet field offset ~ sqrt(E), magnet constant or
flight path ~ E, time zero ~ E^1.5), and adding each as a fixed keyword was
expensive. Write f so that it vanishes inside the data, e.g.
`a0*(sqrt(x)-sqrt(1.7))`: then a is the shift at that energy and its penalty
means what it says.

Model: the p+p model of tests/energy_shift_sqrt (itself tests/identical_pp_res),
E_cm 0.2-4.8 MeV. `xs_40.dat` and `xs_90.dat` are the model x 1.01 with 1 %
errors (chi2 = 24 unshifted). Generated files:

* `xs_40_x15_preshifted.dat`: `xs_40.dat` written at E + 0.002 + 0.003*E^1.5;
* `xs_40_lin_mislabelled.dat`, `xs_90_lin_mislabelled.dat`: labels L with
  L - 0.004 + 0.006*(L - 1) = E exactly, and the model's exact values as yields,
  so a = -0.004, a0 = 0.006 in `a0*(x-1)` is the chi2 = 0 truth.

`shift_fn.azr` (mode 1, pinned by `expected/chiSquared.out`):

1. `xs_40` with a = 0.004 and `sqrtshift 0.003`;
2. the same with `shiftfn "a0*sqrt(x)" 1 0.003 0 0`: must equal segment 1;
3. `xs_40` with a = 0.002 and `shiftfn "a0*x^1.5" 1 0.003 0 0` and
4. `xs_40_x15_preshifted.dat` without shift: must equal segment 3;
5. and 6. the two mislabelled files as group G, the function `a0*(x-1)` on the
   leader at a0 = 0.005 (the truth is 0.006, so chi2 is not zero);
7. `xs_90`, untouched: 24.

`check.sh` adds, with tolerances on the model values in `AZUREOut` (10 digits):
the sqrt(E) equivalence to 1e-10 with and without target integration (sub-points
are shifted by f at their own energies), the centred form, the pre-shifted twin,
a follower against the same function given directly, fits of a and a0 from 0 on
the mislabelled group (numerical MIGRAD, `--use-gradient`, `--use-lm`), the
group-local shift gradient against full-dataset finite differences on both
`--use-gradient` paths, the start-up errors (a follower declaring a function, a
function together with `sqrtshift`, an undeclared or unused parameter, a short
block, a term with no effect at the group's energies), the degeneracy warning,
and that a `.sav` value does not free a fixed function parameter.
