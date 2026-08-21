---
name: azr-project-builder
description: Create a brand-new AZURE2 .azr project from scratch without the Qt GUI — declare particle pairs, photon (capture) pairs, level schemes and channels, and emit a valid <levels> block. Use whenever an evaluation needs a compound nucleus that has no existing .azr to edit, when AzrModel.add_level refuses because the particle pair does not exist yet, or when the 31-field NucLine format has to be read or written by hand.
---

# Building an AZURE2 project from nothing

`pyazr.AzrModel` edits an existing `.azr`: it can add levels, channels, data
segments and extrapolation grids, but it **cannot create a particle pair**, so
a compound nucleus that has never been set up cannot be reached from it. The
GUI is the documented way to create one. This skill is the headless way.

The key fact that makes it possible: **there is no `<pairs>` section**. Every
channel line in `<levels>` repeats the full particle-pair description, and
`CNuc::Fill` (`src/CNuc.cpp:130`) rebuilds the pair list by de-duplicating
those repeats with `CNuc::IsPair`. Writing correct level lines therefore
creates the pairs as a side effect.

A ready-made implementation lives in `evaluations/_tools/azrbuild.py` of the
IAEA/AI-R evaluation tree: `Particle`, `Pair`, `Channel`, `Level`,
`AzrProject`, plus `particle_pair()`, `gamma_pair()`, `make_level()` and a
`validate()` that catches the model errors AZURE2 reports only by refusing to
start. Prefer extending it to re-deriving the format.

## The 31 fields of a level line

Stream order from `include/NucLine.h`, one line per **channel** of each level,
blank line between levels:

```
levelJ levelPi levelE levelFix aa ir 2S 2L levelID isActive channelFix
gamma j1 pi1 j2 pi2 e2 m1 m2 z1 z2 entranceSepE sepE j3 pi3 e3
pType chRad g1 g2 ecMultMask
```

| field | meaning | trap |
|---|---|---|
| `levelE` | compound-nucleus **excitation** energy, MeV | not the resonance energy in the channel |
| `levelFix` | 1 = energy fixed, 0 = fit parameter | per level, repeated on every one of its lines |
| `aa` | legacy entrance key | GUI always writes `1`; `CNuc` never reads it |
| `ir` | this channel's pair key | what `entrance`/`exit` in a segment refers to |
| `2S`, `2L` | doubled channel spin and orbital momentum | integers; `NucLine` halves them on read |
| `channelFix` | 1 = width fixed | per channel |
| `gamma` | Γ in eV (open), ANC in fm^-1/2 (closed), Γ_γ in eV (photon) | **not** a reduced-width amplitude |
| `j1 pi1` | light partner spin/parity | for a photon pair the GUI writes `1  1` |
| `j2 pi2 e2` | heavy partner spin/parity/**excitation** | `e2` carries the final-state Ex of a capture pair |
| `m1 m2` | **nuclear** masses in u | 13N uses 1.00728 for the proton, i.e. atomic − Z·mₑ |
| `entranceSepE` | separation energy of the **first** pair | same value on every line of the file |
| `sepE` | this pair's separation energy | 0 for a photon pair |
| `j3 pi3 e3` | unused | GUI writes the literal `0    0          0.0` |
| `pType` | 0 particle+particle, 10 particle+gamma, 20 beta decay | `src/CNuc.cpp:432` |
| `chRad` | channel radius, fm | 0 for a photon pair |
| `g1 g2` | nuclear g-factors | only used for M1 external capture |
| `ecMultMask` | external-capture bitmask, per pair | bit0 E1, bit1 M1, bit2 E2 (`Constants.h:15`); 5 = E1+E2 |

## Rules the code enforces silently

**A pair's threshold sits at `Ex = sepE + exE`** (`AMatrixFunc.cpp:144`,
`AdaptiveIntegrationGrid.cpp:228`). Two consequences:

- A **photon pair feeding an excited final state** carries `sepE = 0` and
  `e2 = Ex(final)`. Capture to the ⁷Be 1/2⁻ state is `sepE=0, e2=0.42908`, not
  `sepE=-0.42908`.
- A compound nucleus **unbound** against a pair takes a **negative** `sepE`:
  α+α in ⁸Be is −0.0918, n+α in ⁵He is −0.7348, p+α in ⁵Li is −1.9650. AZURE2
  accepts this; `E_cm = Ex − sepE` then exceeds `Ex`, which is correct.

**The radiation type of a photon channel is not stored** — it is derived
(`src/AChannel.cpp:19`):

```
radType = 'E'  if  levelPi * pi2 == (-1)^L   else  'M'
```

so `L` is the multipolarity and `S` is **the spin of the final state**. One
`(L, S)` pair therefore yields exactly one of E*L* or M*L*, never both; to get
both E1 and M1 you need two different final states or two different level
parities. `S` for a photon channel is not a channel spin in the usual sense.

**All levels of one J^π must carry an identical channel set.** This is a
J-group requirement, not a suggestion — a mismatched set makes AZURE2 refuse
the file with no useful message. `AzrProject.validate()` checks it.

**Identical-particle pairs are auto-detected** (`src/PPair.cpp:33`) from equal
Z, mass, spin, parity and `e2`, with the boson/fermion sign taken from `2j`.
Declaring d+d, α+α, p+p or ³He+³He needs nothing extra, but the symmetrisation
it triggers restricts which `(L, S)` survive, so enumerate channels and let
AZURE2 drop what it must.

**`<parameterSettings>` may be empty.** It only carries limits, errors and
nuisance flags (`src/ParameterLimitsManager.cpp:116`); free-versus-fixed comes
from `levelFix` and `channelFix`. A new project can ship the block with its two
comment lines and nothing else.

**Level ordering.** The GUI writes levels sorted by ascending J, then parity
(−1 before +1), then energy, and assigns `levelID` sequentially over that
order. AZURE2 does not require the order, but matching it keeps `levelID`,
`param.sav` line numbers and the GUI's own display consistent.

## Channel enumeration

For a particle pair, `(L, S)` is allowed when

```
|J1 − J2| <= S <= J1 + J2 ,   |J − S| <= L <= J + S ,
pi1 * pi2 * (−1)^L == levelPi
```

For a photon pair, `S = J_final` and `1 <= L <= |J − J_final| … J + J_final`.
`azrbuild.allowed_channels()` implements both and is the function to reuse.

## Minimum viable project

```python
from azrbuild import (AzrProject, TestSegment, E1, E2,
                      particle_pair, gamma_pair, make_level)

p1 = particle_pair(1, "3He", "4He", "7Be", radius=4.2)   # sepE from AME2020
p2 = particle_pair(2, "p", "6Li", "7Be", radius=4.0)
p3 = gamma_pair(3, "7Be", 1.5, -1, 0.0,     ec_mult_mask=E1 | E2)
p4 = gamma_pair(4, "7Be", 0.5, -1, 0.42908, ec_mult_mask=E1 | E2)
pairs = [p1, p2, p3, p4]

levels = [make_level(1.5, -1, 0.0,   pairs, seed=1.0),
          make_level(3.5, -1, 4.57,  pairs, seed=1e3, free=(1,))]

AzrProject("7Be", pairs, levels, data=[], tests=[...]).write("7Be/7Be.azr")
```

`write()` runs `validate()` and creates `output/`, `checks/` and `data/`
alongside the file.

## Verifying a newly built file

1. **CLI mode 3 first.** `printf '3\n\nn\n' | AZURE2 --no-gui new.azr` needs no
   data and prints the compound-nucleus construction, the boundary conditions
   and the parameter transformation. It is the fastest way to see the real
   error message; the pyazr exception only says "could not initialize".
2. **pyazr needs at least one `<segmentsData>` line.** A project with test
   segments only initialises fine under CLI mode 3 and fails under
   `azure2(...)`, whose default is data mode. Add a data segment before loading
   it from Python.
3. **Read the transformation warnings.** *"Denominator less than zero while
   transforming"* means the seed width exceeds what the level shift can absorb
   — the seed, not the model, is wrong. *"radiative width larger than 10% of
   the particle width"* means a seed γ width is unphysical for the same reason.
4. Then `cat output/chiSquared.out` as usual.

## What this does not cover

- **Beta-delayed pairs (`pType = 20`)** are read by the same line format but
  add `j3/pi3/e3` semantics this skill has not exercised.
- **`<targetInt>`** (target integration, straggling, stopping-power formulae)
  is emitted empty. Copy a block from `tests/13N/13N.azr` and edit it when a
  thick-target dataset needs it.
- Channel radii are an input, not a derived quantity. Start from
  `1.4 (A1^(1/3) + A2^(1/3))` fm and treat the radius as a model parameter to
  be scanned, as `azure2-eval` describes.

## See also

- `azure2-eval` — running, fitting and decomposing once the project exists.
- `nds-explorer` — filling `data/` from EXFOR and the level scheme from ENSDF.
- `evaluations/_tools/azrbuild.py` — the implementation.
- `gui/src/LevelsTab.cpp:726` — the GUI's own writer, the authority on field
  order and widths.
