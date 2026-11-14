# The `<levels>` block: pairs, levels and channels

This reference gives the line format of the `<levels>` block of an `.azr` project and the rules the engine applies to
it, so that a project can be written or audited by hand. The authoritative column table is
`docs/source/reference/data_formats.rst` ("Project Files"); `pyazr.Project` edits rows through the engine's writer
(`io/LevelRowText.h`).

## 1. There is no pairs section

Every channel line repeats the full particle-pair description. `CompoundSystem::Fill` rebuilds the pair list by
de-duplicating those repeats (`CompoundSystem::IsPair`). Writing correct level lines therefore creates the pairs.
`pyazr.Project` edits an existing format-3 project (levels, channels, pair columns, records) but cannot create a
particle pair; a compound nucleus that has never been set up is reached by writing level lines directly.

An external helper module, `azrbuild.py` (kept in the evaluation tree, not in this repository), implements
`Particle`, `Pair`, `Channel`, `Level`, `AzrProject`, `particle_pair()`, `gamma_pair()`, `make_level()`,
`allowed_channels()` and a `validate()` that catches the model errors AZURE2 reports only by refusing to start.
Prefer extending it to re-deriving the format.

## 2. The 31 fields of a level line

One line per **channel** of each level, a blank line between levels (`io/include/azure2/io/LevelRows.h`):

```
levelJ levelPi levelE levelFix aa ir 2S 2L levelID isActive channelFix
gamma j1 pi1 j2 pi2 e2 m1 m2 z1 z2 entranceSepE sepE j3 pi3 e3
pType chRad g1 g2 ecMultMask
```

| field | meaning | trap |
|---|---|---|
| `levelJ` | J itself (not 2J) | reading it as 2J turns a 3⁻ level into "3/2⁻" |
| `levelE` | compound-nucleus excitation energy E_x (MeV) | not the resonance energy in a channel |
| `levelFix` | 1 = energy fixed, 0 = fit parameter | per level, repeated on each of its lines |
| `aa` | legacy entrance key | the Qt GUI wrote `1`; `CompoundSystem` never reads it |
| `ir` | this channel's pair key | what `entrance`/`exit` of a segment refer to |
| `2S`, `2L` | doubled channel spin and orbital momentum | integers; the engine halves them on read |
| `channelFix` | 1 = width fixed | per channel |
| `gamma` | Γ (eV) for an open particle channel, ANC (fm^-1/2) for a closed one, Γ_γ (eV) for a photon channel | not a reduced-width amplitude (§4) |
| `j1 pi1` | light partner spin and parity | a photon pair carries `1 1` |
| `j2 pi2 e2` | heavy partner spin, parity, excitation (MeV) | `e2` carries the final-state E_x of a capture pair |
| `m1 m2` | **atomic** masses (u), AME | a bare proton (1.00728 u) among atomic masses draws the electron-mass hint |
| `entranceSepE` | separation energy of the first pair (MeV) | the same value on every line of the file |
| `sepE` | this pair's separation energy (MeV) | 0 for a photon pair |
| `j3 pi3 e3` | unused for particle and photon pairs | the Qt GUI wrote `0 0 0.0` |
| `pType` | 0 particle+particle, 10 particle+γ, 20 β decay | `src/CompoundSystemBuild.cpp` |
| `chRad` | channel radius a_c (fm) | 0 for a photon pair |
| `g1 g2` | nuclear g-factors | used only for M1 external capture |
| `ecMultMask` | external-capture multipoles of this pair | bit0 E1, bit1 M1, bit2 E2 (`include/Constants.h`); 5 = E1+E2 |

Optional fields: **32** is `e1`, the light particle's excitation energy (MeV), written only when nonzero; the pair's
total excitation is e1 + e2. On a THM entrance line field 32 is the binding energy B (MeV) and field 33 the
`gammaIsRWA` flag (azure2-thm). Engines that read 31 columns ignore column 32 silently and open the pair at the
wrong threshold.

Token map for hand parsing (0-based, whitespace split): `t0` J, `t1` parity, `t2` E_x, `t3` energy-fixed flag, `t5` pair
key (the file's key, not pyazr's engine pair number), `t6` 2S, `t7` 2L, `t8` levelID, `t10` channel-fixed flag, `t11`
gamma. Prefer `Project.levels()` (`lv.J`, `lv.parity`, `lv.energy`, `c.pair`, `c.l`, `c.s`, `c.width`) and `pairs()`.

## 3. Rules the engine applies without telling you

- **Thresholds.** A pair opens at E_x = sepE + e1 + e2. A photon pair feeding an excited final state has
  `sepE = 0` and `e2 = E_x(final)`: capture to the ⁷Be 1/2⁻ state is `sepE = 0, e2 = 0.42908`, not
  `sepE = -0.42908`. A compound nucleus unbound against a pair takes a negative `sepE` (α+α in ⁸Be −0.0918 MeV,
  n+α in ⁵He −0.7348 MeV, p+α in ⁵Li −1.9650 MeV); E_c.m. = E_x − sepE then exceeds E_x, which is correct.
- **Photon multipole type is derived, not stored** (`src/AChannel.cpp`): E if levelPi · pi2 = (−1)^L, else M.
  For a photon channel L is the multipolarity and S is the spin of the final state. One (L, S) yields exactly one
  of EL or ML; E1 and M1 together need two final states or two level parities.
- **All levels of one J^π carry an identical channel set** (one J-group). A mismatched set makes AZURE2 refuse the
  file without a useful message.
- **Identical particles are auto-detected** (`src/PPair.cpp`) from equal Z, mass, spin, parity, g-factor and
  excitation (e1 = e2), with the boson/fermion sign from 2j. α + α* is not identical. The symmetrisation restricts
  which (L, S) survive: enumerate all channels and let AZURE2 drop what it must.
- **`<parameterSettings>` may be empty.** It carries limits, errors and nuisance flags only; free versus fixed comes
  from `levelFix` and `channelFix`.
- **Level ordering.** The Qt GUI sorted by J, then parity (− before +), then energy, and numbered `levelID` in that
  order. AZURE2 does not require it; matching it keeps `levelID` and `param.sav` names consistent with files it saved.
- **Energies written by a program**: round E_x to at most 8 decimals; a 16-significant-digit token in an added
  level line has made the engine misread the model.

## 4. The `gamma` field is a physical width, not an amplitude

It holds what `parameters.out` prints: Γ (eV) for an open particle channel, the ANC (fm^-1/2) for a closed
(sub-threshold) channel, Γ_γ (eV) for a photon channel, i.e. exactly `s.transform(x)` at the width rows. In the ⁷Be
model the reduced-width amplitude and the physical value differ by factors 10² to 10⁷; a file written with amplitudes
loads without complaint and is wrong. (THM lines with `gammaIsRWA = 1` are the one exception; azure2-thm.)

## 5. Channel enumeration

For a particle pair (L, S) is allowed when

```
|J1 − J2| <= S <= J1 + J2,   |J − S| <= L <= J + S,   pi1 · pi2 · (−1)^L = levelPi
```

For a photon pair S = J_final and max(1, |J − J_final|) <= L <= J + J_final. Drop high-L and strongly closed
channels only with a stated reason.

## 6. Partial waves without levels: `includeAllJPi`

The elastic T-matrix diagonal is e^{2iω} − e^{2i(ω−φ)}(1 + U): a J^π group with channels and no level still
scatters through its hard-sphere phase φ. A group exists only if an active level line creates it, unless the
`<config>` block carries `includeAllJPi=1` (and optionally `maxOrbitalMomentum=N`, default the largest particle
l in the file); then a level-less group is added for every J^π the pairs allow. Without either, the p+¹²C elastic
cross section is off by up to 14 % (`tests/hard_sphere`). Consequences:

- A file without the key computes as before: check the Levels tab status line ("Not in the calculation: …") or
  `grep includeAllJPi` before trusting an elastic curve. Old models used a fixed zero-width "dummy" level per J^π
  for the same purpose; such dummies stay harmless and are not background poles.
- With the key on, deactivating a level removes the resonance and keeps the partial wave; all levels off gives the
  pure hard-sphere + Coulomb curve. With it off, a J^π disappears with its levels.
- Capture and β pairs get nothing from it.

## 7. External-capture final states

- The final state of every γ pair that a segment (or the rate) exits into must be a level of the project, carrying
  the ANC in its particle channels and a γ line. Otherwise the run stops with
  `ERROR: external capture into pair K needs its final state ...`; add the level at the pair's E_x. A pair no
  segment uses only draws a warning. Test: `tests/ec_final_state_listed`.
- The final state is the level of its J^π within 1 keV of the pair's E_x. Two such levels stop the run
  (`tests/ec_final_state_ambiguous`): move the extra level, or set the pair's energy to the level it decays to.
- **Does a transition need external capture?** Direct (external) capture runs from the channel radius outward with an
  amplitude proportional to the final state's ANC; it dominates non-resonant capture to a loosely bound state, and
  matters little for capture to a deeply bound state with a small ANC. Switch it on per photon pair with `ecMultMask`
  (E1, M1, E2 bits; `ignore_externals` is another switch, of the resonant external widths); take the ANC from a transfer
  measurement or a structure calculation with its error, never from the fit alone; then measure its share: the `bare`
  column of azure2-calculate's `decompose.py` (every resonance's widths zeroed, the final states' ANCs kept) is the
  direct capture alone. The capture operator (long-wavelength or exact, convection current) and the hybrid nuclear
  potential of the external wave: azure2-calculate (its physics-options reference).

## 8. Verifying a newly built file

1. `AZURE2 --validate new.azr` (exit 2 with an error) and `scripts/audit_project.py new.azr`.
2. A calculation of the extrapolations needs no data: `AZURE2 extrapolate new.azr`
   prints the compound-nucleus construction, boundary conditions and the parameter transformation with the engine's
   own messages.
3. Read the transformation warnings. "Denominator less than zero while transforming" means an observed width larger
   than the level can have at that radius (the message gives 2P/(dS/dE), the largest); "radiative width larger than
   10% of the particle width" an unphysical seed Γ_γ. The seed is wrong, not the model.
4. Then a calculation with data, and `output/chiSquared.out`.

## 9. Not covered here

- β pairs (`pType = 20`) use `j3 pi3 e3`; `tests/upgrade/T3` has one (a Gamow-Teller channel).
- Targets, beams and detectors are records of `<targets>` (azure2-targets).
