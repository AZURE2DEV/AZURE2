# The full form factor of the THM transfer vertex (`vertexModel=ff`)

Design note, 2026-11-17. Branch `feature/upgrade-5-thm-ff`.

## What it is

The THM model multiplies each level's reduced width by a transfer vertex. Until now the vertex was the plane-wave
surface term M_l(pa; B) = (B − 1) j_l(pa) − pa j_l′(pa) (`vertexModel=pw`), or the surface term of the prior-form
DWBA (`dw`).

The Catania groups convert peak areas with the full form factor of a single-particle level function: La Cognata et
al., ApJ 708 (2010) 796, Eqs. 20–22, W = ⟨φ_λ|V|j_l⟩. M. La Cognata's review (Nov 2026, items T4 and T5, annotations
A071, A077, A143, A148 and A222) asks for that form, with the dependence on the node count.

Green's theorem splits W exactly into the surface term and a volume term:

W = (ħ²/2μ) u(a) [M_l + Δ_λ(E)],  Δ_λ = (2μ/ħ²)(e_λ − ε)/u(a) ∫₀ᵃ u r j_l(pr) dr.

`vertexModel=ff` adds Δ_λ to M_l for every level and entrance channel.

## Where it plugs in

- **Physics** (`physics/transfer`, a leaf, no engine headers):
  - `azure2/physics/transfer/FormFactor.h` holds the Woods–Saxon + Coulomb level function with a given node count and
    boundary value, the overlaps, the volume term, M_l, the Wildermuth rule and a Chebyshev series.
  - It is exposed in the binding as `pyazr.physics.thm_full_form_factor` (R, Δ, M_l, V₀), so the web UI can plot it.
- **Keys** (`src/ThmKeys.cpp`; one table drives the parser, `pyazr.kinds()["thm"]` and `reference/thm_keys.rst`):
  - `vertexModel` gains `ff`;
  - `ffWell=<r0>,<a0>`, `ffNodes=<J><±>[,k<k>][,l<l>]:<n>;...` and `ffVolume=on|off`;
  - an `ff*` key without `vertexModel=ff` is refused.
- **Set-up** (`src/ThmFormFactor.cpp`, `ThmBuildFormFactor`, called from `WorkingCopy::BuildThmGroups`):
  - one constant `ThmFormFactorTable` per experiment, held like the DW vertex (group → segments → `TermSettings`);
  - per (J group, level, entrance channel): the node count, E_λ, B_λ = S_c(E_λ), the level function and the Chebyshev
    series of g(p) = ∫ u r j_l(pr) dr / u(a).
- **Evaluation** (`src/THMMatrixFunc.cpp`, `ThmEntranceVertices`):
  - the plane-wave vertex of a level becomes γ (M_l + Δ_λ(ε));
  - ε = E + B_xs + T_s of the point, or of the spectator-window node, is stored beside the vertex pieces
    (`ThmVertexRows::eps`, `ThmPsTable::eps`);
  - no cache and no state: the table is built once and read `const`.

## The one restriction, and why

**Fixed energies.** The level energies of the levels the table covers must be fixed (refused at set-up otherwise).

- The level functions then do not change with the fit: the table is part of the problem, built once, and the
  analytic gradient (through the widths) is exact without the derivatives of Δ in E_λ and B.
- The review prescribes literature resonance energies for THM analyses, so this costs nothing in practice.
- A guard in the evaluation solves a level function again, with a once-only warning, if a level's energy or boundary
  ever differs from the table's.

**Basis and boundary.** `ff` also needs `vertex=perlevel` in the Brune basis: B_λ is the level function's own
boundary value, the one the surface term uses.

## What it replaces

Nothing is removed. `pw` and `dw` stay, and `ff` with `ffVolume=off` is `pw` with `vertex=perlevel`, bit for bit.

The `thm-review` branch's implementation (4882ad04, never pushed) is the source of the physics code. Its per-thread
cache keyed by E_λ and B, and its finite-difference derivatives of Δ, are not ported: both are unnecessary with fixed
energies, and the cache is excluded by the "no hidden state" rule.

## Tests

- `physics_transfer_form_factor`:
  - Green's identity to 1e−8 (W direct = surface + volume);
  - boundary value, node count and norm of every level function;
  - R = W/surface for the n + ¹⁷O levels against an independent NumPy implementation to 1e−4;
  - j_l, the Wildermuth rule and the Chebyshev series against closed forms.
- `pyazr_physics_table`: the binding's row against the same independent values.
- `tests/thm_full_form_factor/check.sh`:
  - `ffVolume=off` identical to `pw`;
  - `ff` moves the χ² of a THM project;
  - the refusals.

## Open (for the Catania groups)

- Their potential and node counts.
- Whether they normalize the level function inside a and stop the integral at a.
- Whether their volume term is the one computed here (the reference values are illustrations with our well).
