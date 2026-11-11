# spin_orbit_potential — the spin-orbit term of a hybrid pair

`hybrid_potential` (p + ¹²C, Woods–Saxon 20 MeV / 3.6 fm / 0.6 fm, R_C = 3.5 fm, a_c = 3.6 fm) written in format 3,
with a spin-orbit term on the proton:

```
<potentials>
potential p12C
    pair           p+12C
    coulomb_radius 3.5 fm
    central        woods_saxon depth 20 MeV radius 3.6 fm diffuseness 0.6 fm
    spin_orbit     woods_saxon depth 6 MeV radius 3.6 fm diffuseness 0.6 fm on light operator 2l.s
    ...
</potentials>
```

The term is v_so(r) ⟨2 l·i⟩ with v_so = V_so (2.000 fm²) (1/r) df/dr and ⟨2 l·i⟩ = j(j+1) − l(l+1) − i(i+1), FRESCO's
convention (`docs/source/theory/hybrid_potential.rst`); with one spinning particle j = J, so the channel functions of
the pair depend on (l, J).

| file | what |
|---|---|
| `spin_orbit_potential.azr` | the Meyer elastic data at 144.1° of `hybrid_potential` (no target integration); χ² pinned in `expected/chiSquared.out` (2117.49; 2259.87 without the term) |
| `models/phase.azr` | p₁/₂, p₃/₂ and s₁/₂ phase shifts at E_cm = 0.5–3 MeV, levels without widths: δ = −φ of the potential |
| `models/capture.azr` | `13N` cut to the LUNA HPGe ¹²C(p,γ)¹³N cross section, with external capture into the ground state |
| `phase_check.py` | δ against `physics/tests/data/spin_orbit.csv` (`make_spin_orbit_reference.py`: scipy DOP853 + mpmath, independent of the engine), 1e-6 rad |
| `unitarity.py` | n + ¹⁶O with a term on the neutron: σ_tot (optical theorem) = σ_el |
| `ec_reference.py` | the external-capture integrals with / without the term for the s- and d-wave E1 pathways into the ground state, by an independent integration (`expected/ec_ratio.csv`); compared with `pyazr.ec_integrals` in `tests/pyazr/spin_orbit_potential_test.py` |
| `check.sh` | the phase shifts, unitarity, 1 against 3 threads (three models), the gradient against finite differences, the term reaching the capture path |

Injected errors that the checks were shown to catch (each once, then removed): the channel row of a point read with
J = l + 1/2 for every channel (phase shifts off by 2.8e-2 rad); the bound state of external capture without its
J (`ec_integrals` ratio off by 1.7e-2); in `physics_spin_orbit`, l·s instead of 2 l·s, the sign of df/dr, ⟨2 l·i⟩ halved.
