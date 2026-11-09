# AZURE2 upgrade: the kinds of data a measurement can supply (Phase 4 survey and proposal)

This note lists the quantities an R-matrix evaluation of light nuclei takes from experiments and from EXFOR, says
which of them AZURE2 accepts today and exactly what it predicts for each, gives the defining formula of every kind in
terms of the collision matrix or the level parameters, and proposes how each is written as a model expression
(`ARCHITECTURE.md` §3, decision 12) and as an explicit, named quantity in the `.azr`. It answers the maintainer's
request: "The data input could be more flexible, as for example yield data that links to a target effect, or
Rutherford ratio data or other type of data as omega gamma or other that can be found on EXFOR." It is input to the
Phase 4 design note (`ROADMAP.md`, Phase 4); it changes no code.

Marks used throughout: **[V]** verified in the code or in the cited document on 2026-11-09; **[I]** inferred, not
checked by running anything. Code references are to `feature/upgrade-4-data` at `eb79d95a`. EXFOR references are to

- LEXFOR, the EXFOR Compiler's Manual, IAEA-NDS-208 Rev. 2015/08, ed. O. Schwerer
  (`https://nds.iaea.org/nrdc/nrdc_doc/iaea-nds-0208-201508.pdf`), cited as "LEXFOR p. X.n";
- the EXFOR/CINDA dictionaries in the IAEA's JSON form (`github.com/IAEA-NDS/exfor_dictionary`,
  `src/exfor_dictionary/latest.json` = `trans.9133`, fetched 2026-11-09), cited as "Dict. 236" (quantities), "Dict. 34"
  (modifiers), "Dict. 25" (units);
- NRDC working papers wp2012-28 and wp2021-08 for the correlation flags of `ERR-ANALYS` (the EXFOR Formats Manual,
  IAEA-NDS-207, could not be downloaded: the IAEA server answered 403 behind a Cloudflare challenge).

## 1. The comparison every kind feeds

A data set is a list of points i with a datum d_i and an uncertainty e_i. The engine predicts m_i(θ, ν) from the level
parameters θ and the nuisances ν, and adds

```
χ²_data = Σ_i ((m_i − n d_i) / (n e_i))²                                             (1)
```

to the objective, n being the normalization of the data set, which multiplies the data, not the model
(`src/AZURECalc.cpp:110-117`; objective term T1 in `objective.md` §1) **[V]**. A normalization with a quoted relative
uncertainty σ_n (per cent) adds ((n − n₀)/(n₀ σ_n/100))² (`src/AZURECalc.cpp:901`, term T3) **[V]**. A covariance file
replaces the sum of squares by rᵀC⁻¹r and can carry rank-one terms `vector` (absolute) and `relvector` (fractions of
the model) for correlated systematics (`docs/source/reference/data_formats.rst:68-120`) **[V]**.

"A kind of data" is therefore three things: what the datum d_i is (quantity, unit, frame), the function m_i of the
collision matrix U(E) or of the level parameters, and the apparatus fold that turns the first into the second.

## 2. What AZURE2 accepts today

### 2.1 The data file and the segment line

A data file has four columns, `E  θ  y  δy`: lab energy (MeV, forward kinematics: light particle on heavy target), lab
angle (degrees; a dummy for angle-integrated data), the value and its uncertainty in b or b/sr
(`docs/source/reference/data_formats.rst:7-66`) **[V]**. A `<segmentsData>` line reads
`active entrance exit Emin Emax θmin θmax isDiff [J l] norm varyNorm normError [shift shiftError vary] file
[composite block] [UPOS block] [keyword blocks]` (`include/SegLine.h:103-106, 146-210`) **[V]**. The observable is the
integer `isDiff`, decoded in `src/ESegment.cpp:38-83` **[V]**:

| `isDiff` (data) | kind | value column | frame of E, θ, value | source |
|---|---|---|---|---|
| 0 | angle-integrated σ (any exit; with a γ exit the sum of all multipoles) | b | E lab; no angle | `ESegment.cpp:50-57` |
| 1 | dσ/dΩ | b/sr | E, θ, value lab; converted to c.m. at load | `EPoint.cpp:1023-1060` |
| 2 | phase shift δ of one (J, l, s) wave | degrees | E lab; δ = ½ arg(U_cc e^{−2iω_c}) | `data_formats.rst:122-141` |
| 3 | total capture (sum over all γ pairs; exit key −1) | b | E lab | `data_formats.rst:366-371` |
| 4 | dσ/dΩ | b/sr | E lab; θ and value c.m. | `ESegment.cpp:53-57` |
| 5, 6 | angle-integrated capture, E1 part only, E2 part only | b | E lab | `ESegment.cpp:69-76` |
| 7 | vector analyzing power A_y (Madison), spin-½ projectile only | — | E lab; θ c.m. | `ESegment.cpp:44-49`; `EData.cpp:57-99` |
| 8 | refused (former P·dσ/dΩ) with advice to use A_y of the inverse reaction | — | — | `EData.cpp:397`; `segments.rst:177-186` |
| 9 | Legendre coefficient a_k = A_k/A_0 of dσ/dΩ_cm, or A_k (b/sr) with `absolute 1` | — / b/sr | E lab; k in the angle column | `data_formats.rst:143-154`; `SegLine.h:351-358` |
| +10 | the same observable as a THM (half-off-shell) measurement | arbitrary | `THM.md` | `ESegment.cpp:38-43` |

Test segments (`<segmentsTest>`) use another numbering (3 = angular distribution, 4 = total capture, 5 = c.m.
differential: `src/ESegment.cpp:166-196`; `pyazr/datasets.py:90-98`) **[V]**. The same integer therefore means
different observables in the two tables.

Further forms, all on the same line **[V]**:

- **Composites.** After the file, `1 <op> <components>`: op 0 = SUM, op 1 = RATIO; components are (entrance, exit,
  angle, scaling) (`include/ComponentList.h`, `SegLine.h:188-198`). SUM adds the segment's own term and every
  component times its scaling, each in the base point's frame (`ESegment.cpp:1748-1773`). RATIO is the segment's own
  term over the **first** component, with the lab Jacobian of each term's own pair (`ESegment.cpp:1775-1815`). Both
  terms are R-matrix observables: there is no analytic denominator such as the Rutherford cross section.
- **Natural and enriched targets.** `sample <1|name>` weights the term and each component by the isotopic abundance
  of its compound nucleus in an `<abundances>` sample, abundances being optional fit parameters with priors
  (`SegLine.h:316-320`; `data_formats.rst:325-436`). Refused for a RATIO (`EData.cpp:209`).
- **UPOS.** `isUPOS L Ic δ` after the file: the secondary γ of B → C + γ instead of the undetected particle
  (`SegLine.h:199-210`; `segments.rst:512-552`).
- **Neutron total, transmission, self-shielded yield.** `total 1`: the optical-theorem total cross section;
  `transmission n`: e^{−nσ}; `yield n`: (σ_x/σ_t)(1 − e^{−nσ_t}), n in atoms/b, optionally fitted and shared with
  `thicknessgroup` (`SegLine.h:328-350`; `ESegment.h:366-418`). Neutral entrance pair only (`segments.rst:727`).
- **β-delayed spectra.** A segment whose entrance pair is a β pair (type 20): the energy column is the light
  particle's energy from a decay at rest, converted by E_cm = E₁(m₁ + m₂)/m₂ (`EPoint.cpp:801-812`), and the model is
  Σ 25|T|² in units fixed by the β-feeding parameters (`GenMatrixFunc.cpp:140-141`). The scale is in practice arbitrary
  and needs a free norm **[I]**.
- **Measured scattering lengths.** Not a segment: the `<scatteringLengths>` block, `active pair s value σ [Im σ_Im]`
  in fm, each line a term ((Re a₀ − value)/σ)² of the objective (`data_formats.rst:438-469`;
  `include/ScatteringLength.h`; term T9, `objective.md` §1). This is today's only datum without a beam energy.

**S-factor input.** The engine has no S-factor data type: the value column is a cross section
(`data_formats.rst:20-23`). S factors are outputs only (`EPoint::GetSFactorConversion`, `include/EPoint.h:396`).
pyazr converts an EXFOR `SFC` data set to barns before writing the file, σ = S exp(−2πη)/E with the engine's
constant 0.157488 (`pyazr/nds.py:798-809, 965-984`) **[V]**.

### 2.2 How a target effect is attached today

A `<targetInt>` line begins `active segmentList nPoints isConvolution σ_beam isTargetIntegration density
stoppingEq nPar par…` (`src/TargetEffect.cpp:61`) **[V]**. `segmentList` is a string of **segment indices**
(`3,6,7-14`) (`docs/source/user_guide/experimental_effects.rst:66-78`); a line listing several segments gives each its
own copy of the effect (`include/TargetEffect.h:57-63`). The target has no name. What one physical target means is
spread over three places:

| property | where today | how shared |
|---|---|---|
| areal density N (atoms/cm²), stopping ε(E), straggling, beam spread | `<targetInt>` line | `densityfit vary σ` and `densitygroup <label>` (`TargetEffect.h:89-123`) |
| composition seen by the reactions (isotopic or stoichiometric weights) | `<abundances>` sample | `sample <name>` on segment lines |
| neutron sample thickness n (atoms/b) | segment line | `thicknessgroup <label>` |

The stopping cross section is a user formula in x = E_lab (MeV), "per active atom" being the user's responsibility
(`experimental_effects.rst:263-320`) **[V]**. Nothing ties the stoichiometry used for ε to the one used for the
abundance weights.

### 2.3 What a "yield" predicts today

With `isTargetIntegration = 1` the model of a point at beam energy E₀ is, with every energy in the c.m. frame inside
the code and dE_cm/ε_cm = dE_lab/ε_lab (`src/EPoint.cpp:2816-2830`, 2850-2890) **[V]**:

- **thin target** (density N > 0):

  ```
  m(E₀) = (1/N) ∫_{E₀−Δ}^{E₀} dE ⟨σ⟩_g(E) / ε(E),     ∫_{E₀−Δ}^{E₀} dE/ε = N              (2)
  ```

  where ⟨σ⟩_g(E) is σ averaged over a normalized Gaussian of width (σ_beam² + σ_strag²(E))^{1/2} (two-point
  Gauss–Legendre per sub-interval, `EPoint.cpp:2862-2887`). The unit of m is **b** (b/sr for differential data): the
  "yield" of a thin target is a target-averaged cross section, divided by N·10⁻²⁴ (`EPoint.cpp:2890`;
  `experimental_effects.rst:217-231`).
- **thick target** (density 0):

  ```
  m(E₀) = 10⁻²⁴ ∫_{1 keV}^{E₀} dE ⟨σ⟩_g(E) / ε(E)                                          (3)
  ```

  with σ in b and ε in MeV cm² per active atom: reactions **per incident particle** (b/sr ⇒ per incident particle per
  sr) (`EPoint.cpp:2890`; `TargetEffect.h:150-157`; `experimental_effects.rst:243-256`).

The data file must already hold these units: a target-averaged cross section, or reactions per incident projectile.
Counts per µC, a charge state, a detection efficiency, a γ branching or an arbitrary scale are the user's to remove
(the last through a free norm). A beam-current-integrated yield in µC⁻¹ needs Y_part = Y_µC · q e / 10⁻⁶ C with q the
charge state; nothing in the engine knows q **[V]** (no such field exists in `TargetEffect.h` or `SegLine.h`).
Finite detector solid angle enters through Q-coefficients on the Legendre terms (`TargetEffect.h:32-34`;
`experimental_effects.rst:912-940`).

### 2.4 Ratio-to-Rutherford data today

The engine cannot express it **[V]**: no observable or composite has the Rutherford (or Mott) cross section as a term
(`grep -i rutherford src include` finds only the Coulomb amplitude in `EPoint.h:441`), and a RATIO's denominator is an
R-matrix observable (`ESegment.cpp:1780-1812`), i.e. the full Coulomb + nuclear elastic cross section, not σ_R.
The only route is to multiply the data beforehand. pyazr does this on import: `to_azr(rutherford=True)` multiplies the
ratio and its error by the c.m. Rutherford cross section at the c.m. angle and divides by the lab-to-c.m. Jacobian for
lab angles (`pyazr/nds.py:232-247, 985-1004`) **[V]**; for identical particles it warns that the reference should be
Mott and writes ratio × Rutherford anyway (`nds.py:997-1000`) **[V]**. The masses used are mass numbers
(`nds.py:186-198`) **[V]**, not the atomic masses the engine uses.

### 2.5 The present EXFOR client (pyazr/nds.py)

What it does **[V]**: fetch by dataset ID from the Web-API or a local `exfor_master` copy (`nds.py:463-490`); convert
`EN-CM` to lab energy with mass numbers (`nds.py:834-863`); convert units by prefix (`nds.py:347-382`); convert `SFC`
to σ (`nds.py:965-984`); convert Rutherford ratios on request (above); take the point error from the first of
`ERR-S`, `DATA-ERR`, `ERR-T` that exists (`nds.py:894-910`); take the first row of `ERR-SYS` in per cent as the norm
error (`nds.py:284-292`); guess the observable (`nds.py:294-316`).

What it does not do **[V]**:

1. It refuses no quantity. `TTY`, `WID/STR`, `SGV`, `LEG`, `REL`, ratios and sums of reactions arrive as numbers in a
   b column. `NO-DIM` and `ARB-UNITS` values pass through unconverted (`nds.py:1005-1006`) and no free norm is set.
2. A point without an error gets 5 % with no warning (`nds.py:85, 909-910`): a silent fallback.
3. Any (x,γ) reaction is guessed to be total capture (`nds.py:304`), including a `PAR,SIG` to one final level; checked
   by calling `_guess_observable` on `6-C-12(P,G)7-N-13,PAR,SIG`, which returns `total-capture`.
4. The tautology form of inverse kinematics, `(1-H-1(6-C-12,EL)1-H-1,,DA)=(6-C-12(P,EL)6-C-12,,DA)`, is misparsed:
   projectile `'1-H-1(6-C-12'`, target `''` (same check); the mass lookup then fails.
5. An analyzing power with a lab `ANG` column is written with its lab angle as if c.m. (code 7 takes c.m. angles):
   `_angle_column` passes the angle unchanged (`nds.py:866-875`) and `_guess_observable` returns `analyzing-power`
   without looking at the frame (`nds.py:307-308`). Wrong by the lab-to-c.m. angle shift **[I]** (not run on a real
   entry).
6. `ERR-ANALYS` (the meaning and correlation of each error) and the partial errors `ERR-1…`, `MONIT-ERR` are ignored;
   `ERR-SYS` given per point is reduced to its first value.

## 3. What EXFOR holds for light-nucleus R-matrix work

### 3.1 Frames, errors, relative data, combinations

- **Frame.** All quantities are lab unless the data heading carries `-CM` (`EN-CM`, `ANG-CM`, `E-CM`, `NUMBER-CM` for
  Legendre coefficients fitted in c.m. angle, `DATA-CM`). For the cross-section-like ratios — "cross section,
  Rutherford ratio, analyzing power" — the heading is `DATA` whatever the angle frame, because the value is the same
  at corresponding angles (LEXFOR p. C.1) **[V]**. (For a cross section this holds only for a ratio-type value; a
  dσ/dΩ in b/sr is frame dependent, and LEXFOR's sentence must be read as applying to the ratios.)
- **Inverse kinematics.** When the energy is `EN-CM` or per nucleon, A_targ ≤ 4, A_proj ≥ 5 and swapping target and
  projectile does not change the numbers, the reaction is coded with the tautology `((exp. form)=(forward form))`
  (LEXFOR p. I.2) **[V]**.
- **Errors** (LEXFOR p. E.2–E.4) **[V]**: `ERR-S` statistical; `ERR-SYS` total systematic; `ERR-T` total, whose
  components are also given; `ERR-1, ERR-2, …` partial; `MONIT-ERR` monitor; `DATA-ERR` of unknown type or total
  without components; `ERR-DIG` digitizing. Only 1σ values are compiled. The correlation property is a fourth field of
  each `ERR-ANALYS` code: U uncorrelated, F fully correlated, P partially correlated, C correlated (F or P unknown),
  blank unknown (wp2012-28; wp2021-08) **[V, from the working papers, not the Formats Manual]**.
- **Relative data.** Modifier `REL` in SF8 with units `ARB-UNITS`: "Relative data. To be combined with arbitrary
  units. Otherwise use 'FCT'" (Dict. 34; LEXFOR p. D.13, l. 2390–2392) **[V]**.
- **Ratios and sums of reactions.** Separators `/` (ratio, units `NO-DIM`), `//` (ratio at different values of an
  independent variable, with `-NM`/`-DN` headings) and `+` (sum), with `A` or `RAB` in SF8 for natural abundances
  (LEXFOR p. R.1, S.23) **[V]**.
- **Natural elements** are coded with A = 0 (`6-C-0`); the nearly monoisotopic list includes 1-H-1, 2-HE-4, 6-C-12,
  7-N-14, 8-O-16 (LEXFOR p. E.1) **[V]**.
- **Charged-particle elastic σ is forbidden**: `(P,EL),,SIG` (LEXFOR p. S.2) **[V]**.

### 3.2 Quantity codes (Dict. 236, verified in `latest.json`; dimension code after the slash)

| code (SF5,SF6,SF7,SF8) | meaning (Dict. 236 text unless noted) | dim. |
|---|---|---|
| `,SIG` / `PAR,SIG` | cross section / partial (to a level or γ line, `E-LVL`, `E`) | B |
| `,SIG,,SFC` / `PAR,SIG,,SFC` / `EP/PAR,SIG,,SFC` | S factor / partial / electric multipole part | B*E |
| `,DA,,SFC` | S factor as a function of angle (LEXFOR p. A.7) | EDA |
| `,DA` / `PAR,DA` | dσ/dΩ / partial | DA |
| `,DA,,RTH` / `PAR,DA,,RTH` / `,DA,*,RTH` | ratio to Rutherford | NO |
| `,DA,,MOT` | ratio to Mott scattering | NO |
| `,DA,,LEG` | dσ/dΩ = a₀ + Σ a_L P_L | DA |
| `,DA,,LEG/RS` | (4π/σ) dσ/dΩ = 1 + Σ W_L P_L | NO |
| `,DA,,LEG/RSL` | (4π/σ) dσ/dΩ = Σ (2L+1) B_L P_L | NO |
| `,DA,,LEG/2L2` | dσ/dΩ = ½ Σ (2L+1) a_L P_L | DA |
| `,DA,,LEG/RS0`, `,DA,,LEG/RSD` | relative to dσ/dΩ at 0° / 90° | NO |
| `,DA,,LEG/1K2` | k² dσ/dΩ = Σ a_L P_L | NO |
| `,DA,,LEG/4PI` | 4π dσ/dΩ = a₀ + Σ a_L P_L | B |
| `,DA,,LEG/RTH` | Legendre coefficients of the Rutherford ratio | NO |
| `,DA,,COS…` | the same forms in powers of cos θ | |
| `,DA,,RS`, `,DA,,RSD`, `,DA,,RS0` | (4π/σ)dσ/dΩ, dσ/dΩ / dσ/dΩ(90°), / dσ/dΩ(0°) | NO |
| `,POL/DA,,ANA` | analyzing power (beam spin normal to the plane) | NO |
| `,POL/DA,,VAP`, `20,POL/DA,,TAP`, `,POL/DA,,AYY` (inactive) | iT₁₁, T₂₀, A_yy | NO |
| `,POL/DA` | polarization of the outgoing particle | NO |
| `,POL/DA,,ANA/LEG` | Legendre coefficients of the analyzing power | NO |
| `,SIG,,TT` / `,SIG,,TTA` / `,DA,,TT` | cross section measured on a thick target, not corrected / thin-target approximation | B, DA |
| `,MLT,,TT`, `,MLT,,TT/CH`, `PAR,MLT,G,TT/CH` | particles (γ) per incident projectile / per charge | YLD, YLDC |
| `,PY,,TT`, `,PY,,TT/CH` | products per incident projectile / per charge | YLD, YLDC |
| `,TTY,,PHY`, `,TTY,,EOB`, `,TTY,,SAT` | **activation** yields: decay rate per charge or per current | TTT, TTY |
| `,WID/STR`, `PAR,WID/STR`, `,WID/STR,,RG` | resonance strength ωγ; partial; (2J+1)Γ_iΓ_j/Γ | E |
| `,EN` | resonance energy | E |
| `,WID`, `PAR,WID`, `,WID,,G`, `,WID,,2G` | total, partial, g·Γ, 2g·Γ | E |
| `,WID/RED`, `,WID/RED,,RMT`, `,WID/RED,,RMT/AMP` | reduced width; R-matrix γ² (Γ = 2γ²P); amplitude γ | E, RE |
| `,ARE` | resonance area ∫σ dE | B*E |
| `,AMP`, `COH,AMP` | scattering length | L |
| `,SGV`, `PAR,SGV` | ⟨σv⟩ at temperature `KT` | B*V |

Three findings from this check **[V]**:

1. `TTY` is not a prompt reaction yield. It is the activity of a produced nuclide per unit charge or current (units
   MBq/µA, MBq/C; LEXFOR p. T.8–T.9). A prompt γ or particle thick-target yield is `MLT…,TT` or `PY…,TT` (per incident
   particle, `PRT/INC`, `PRD/INC`) or per charge (`TT/CH`, `PRT/MUCOUL`), or a cross section not corrected for the
   target, `SIG,,TT` (LEXFOR p. T.9). The bare `,TTY` code is marked inactive.
2. Ratio-to-reference-scattering codes are only `RTH` and `MOT`. Dict. 34 defines `MOT` as "relative to Mott
   scattering", LEXFOR p. S.4 as the identical-particle Rutherford cross section symmetric about 90°, while Dict. 236's
   text for `,DA,,MOT` says "relativistic scattering by point nucleus" — two different "Mott" cross sections. No code
   for a ratio to "Coulomb plus hard sphere" or to an optical-model cross section exists; such data can only be coded
   with `MSC` and free text. The JSON text of `,DA,,RTH` reads "coeff. for a sum in power of sine**2", which
   contradicts Dict. 34 and LEXFOR; it is an error in the JSON description and the code means the Rutherford ratio.
3. The Legendre normalizations disagree between sources: Dict. 34 writes `RSD` as "(dσ/dΩ)/(dσ/dΩ at 90°) = 1 + Σ
   a_L P_L", LEXFOR p. F.18 as "= Σ a_L P_L", Dict. 236 `,DA,,LEG/RSD` as "4π/σ(90°) dσ/dΩ = 1 + Σ a_L P_L". A reader
   cannot convert `RSD`/`RS0` coefficients without the entry's own text.

### 3.3 How often each occurs

The EXFOR Web-API `x4list` was queried on 2026-11-09 for the targets ²H, ³H, ³He, ⁴He, ⁶Li, ⁷Li, natLi, ⁹Be, ¹⁰B, ¹¹B,
natB, ¹²C, ¹³C, natC, ¹⁴N, ¹⁵N, natN, ¹⁶O, ¹⁷O, ¹⁸O, natO, ¹⁹F, ²⁰Ne, ²¹Ne, ²²Ne, natNe with the projectiles p, d, t,
³He, α (any exit), keeping data sets with a lowest lab energy below 20 MeV (or none quoted). Script:
`scripts/survey.py`, `scripts/classify.py` in the session scratchpad (not committed). Counts are of data sets
(EXFOR subentries or pointers), not of points.

16 989 distinct data sets were returned; 10 132 pass the energy cut. Classified by SF5–SF8 of the reaction code (for a
tautology, the forward form), families are matched in the order listed, so each data set is counted once:

| family | data sets | points |
|---|---|---|
| dσ/dΩ, absolute (`DA`, `PAR,DA`, no ratio modifier) | 4129 | 662 266 |
| σ, absolute (`SIG`, `PAR,SIG`) | 1365 | 54 014 |
| dσ/dΩ, relative (`DA … REL`) | 730 | 60 537 |
| analyzing powers and polarizations (`POL…`; codes containing `ANA` 214, `VAP` 76, `TAP` 163) | 596 | 27 556 |
| resonance energy, width, J, parity (`EN`, `WID`, `PAR,WID`, `J`, `PTY`) | 563 | 5379 |
| thick-target yields of particles or γ (`MLT…TT`, `PY…TT`, `TTY/MLT`, `TTY/PY`; mostly (x,X) neutron production) | 506 | 75 201 |
| S factor (`SFC`, all forms) | 350 | 10 703 |
| double-differential and spectra (`DA/DE`, `DE`, `DA/DA`, `SPC`) | 339 | 38 867 |
| Legendre and cosine coefficients (`LEG…`, `COS…`; of these `LEG/RS` 146, `LEG` 62) | 312 | 29 369 |
| ratio to Rutherford (`RTH`, all forms; `RTH/REL` 5, `LEG/RTH` 1) | 293 | 23 706 |
| dσ/dΩ normalized to σ or to 0°/90° (`RS`, `RSD`, `RS0`, `4PI`) | 189 | 5897 |
| ratios of reactions `/` | 99 | 933 |
| σ, relative (`SIG … REL`) | 107 | 9896 |
| resonance strength (`WID/STR`, `PAR,WID/STR`, `WID/STR,,RG`) | 101 | 446 |
| activation yields (`TTY,,SAT/PHY/EOB`) | 82 | 2223 |
| ratios at different independent variables `//` | 53 | 6463 |
| thermonuclear rate (`SGV`) | 30 | 6299 |
| sums of reactions `+` | 25 | 683 |
| reduced widths (`WID/RED`) | 7 | 33 |
| ratio to Mott (`MOT`) | 0 | 0 |
| other (`IPA`, `SIG,,RES`, `AV`, …) | 256 | — |

Across all families: `REL` in SF8 on 1013 data sets, process `X` (production) on 792, a natural target (A = 0) on 233,
a tautology on 98. α + α elastic data (target ⁴He, beam α) are coded `,DA` (5 data sets), not `RTH` or `MOT`. These
numbers order the proposal (§5.3): after plain σ and dσ/dΩ, the frequent kinds are relative data, polarizations,
resonance parameters, S factors, Legendre coefficients and Rutherford ratios; ωγ is rarer in data sets but each is a
precise number with a large weight in a fit.

## 4. The physics of each kind

Notation. Entrance channel α = a + A with spins j_a, j_A, reduced mass μ, c.m. energy E (MeV), wave number k =
(2μE)^{1/2}/ħ (fm⁻¹), Sommerfeld parameter η = Z_a Z_A e²/(ħv). U is the collision matrix, T = U − 1 with the Coulomb
phases removed as in `docs/source/theory/conventions.rst`. g_J = (2J+1)/((2j_a+1)(2j_A+1)). Lab quantities carry the
subscript lab. Every observable below is computed by the engine in c.m.; the frame conversions are load-time functions
(ARCHITECTURE §3, `Measurement`).

### 4.1 Cross sections, S factors, partial cross sections

σ_{αβ}(E) = (π/k²) Σ_J g_J Σ_{c∈α, c′∈β} |T_{c′c}^J|² (b, c.m. = lab: angle-integrated cross sections are frame
independent). The S factor S(E) = E σ(E) exp(2πη) (MeV b or keV b; LEXFOR p. A.6–A.7, with the actual reduced mass,
not mass numbers) is a diagonal rescaling: a `Scale` node with a factor that depends on E only. Exact if the engine's
η is used. pyazr's conversion uses mass numbers for μ (`nds.py:186-198, 975-982`) **[V]**, which changes 2πη and
hence σ by exp(2πΔη) − 1: +2.9 % for ³He + α at E = 93 keV, +6.7 % for p + ¹²C at 100 keV, +0.9 % for α + ¹²C at
300 keV (computed for this note with the engine's constant 0.157488 and AME masses); the engine then converts back
with atomic masses, so the S factor it reports differs from the one imported by the same amount. A wrong Z is a
factor of 5000 for ³He(α,γ) at 93 keV (`nds.py:66-71`). A partial cross section to a final level or γ line
(`PAR,SIG`, `E-LVL`, `E`) is the leaf of one exit pair; a γ line fed by several levels is a `Sum` of capture leaves to each feeding level times its branching (a
`Scale`); the branching comes from ENSDF with its uncertainty as a nuisance. Apparatus: none for σ itself; target and
resolution folds when the data are not deconvolved.

### 4.2 Differential cross sections and Jacobians

dσ/dΩ_cm(E, θ) = |f_C(θ)δ_{αβ} + f_N(θ)|² summed over spin projections; the Coulomb amplitude f_C is not a finite
Legendre series, so the engine keeps nuclear amplitudes and adds f_C at the point (ARCHITECTURE §3, `observables`).
Lab data: dσ/dΩ_cm = K(θ) dσ/dΩ_lab with K = dΩ_lab/dΩ_cm; for two-body kinematics at θ_lab ≠ 0, 180°
K = sin²θ_lab cos(θ_cm − θ_lab)/sin²θ_cm (`EPoint.cpp:1053-1058`), for photons the aberration Jacobian
(1 − β cos θ)²/(1 − β²) (`data_formats.rst:56-60`). A point has one Jacobian per exit pair, so a lab-frame Sum or
Ratio of different exit pairs needs one per term, as the engine does since October 2026 (`segments.rst`, "Frames")
**[V]**. Inverse kinematics (heavy beam on a light target) changes the lab angles and the lab energy, not the c.m.
quantities: E_cm = E_beam m_light/(m_light + m_heavy); an inverse-kinematics data set is converted to forward lab or
to c.m. at load, with the detected particle and its branch (two c.m. angles per lab angle) stated **[I]**.

### 4.3 Ratio to Rutherford or Mott (treated fully)

**Definition.** For elastic scattering of a + A,

```
R(E, θ) = (dσ/dΩ)(E, θ) / (dσ_R/dΩ)(E, θ),     dσ_R/dΩ_cm = (Z_a Z_A e² / 4E)² / sin⁴(θ_cm/2)      (4)
```

with E the c.m. energy and e² = αħc = 1.43996 MeV fm, giving fm²/sr (×10⁻² for b/sr) (LEXFOR p. S.4; the same
formula as `nds.py:386-396`) **[V]**.

**Frame.** Numerator and denominator are both differential in the same solid-angle element. In elastic scattering
both refer to the same two-body kinematics, so the lab-to-c.m. Jacobian is the same function of angle for both and
cancels: R_lab(θ_lab) = R_cm(θ_cm(θ_lab)). R is frame independent only in this sense — numerator and denominator in
the same frame and at corresponding angles. Three ways to break it: (i) a lab cross section divided by the c.m.
Rutherford value (the pre-October-2026 pyazr bug, 1.13× for p + ¹⁶O and 1.8–2.2× for α + ¹²C, `nds.py:245-247`);
(ii) the c.m. Rutherford value evaluated at the lab angle; (iii) inelastic scattering (Q ≠ 0), where the Jacobians
of the reaction and of elastic Rutherford differ, so a "ratio to Rutherford" of an inelastic group is defined only
together with the angle and frame at which σ_R was taken. The angle still has to be converted for the energy
dependence and for the frame of the stored angle column.

**Identical particles** (α + α, ¹²C + ¹²C, p + p, ¹⁶O + ¹⁶O). The Coulomb amplitude is symmetrized; for two identical
particles of spin s,

```
dσ_Mott/dΩ_cm = (Z²e²/4E)² [ 1/sin⁴(θ/2) + 1/cos⁴(θ/2)
                 + (2(−1)^{2s}/(2s+1)) cos(η ln tan²(θ/2)) / (sin²(θ/2) cos²(θ/2)) ]                 (5)
```

(θ = θ_cm), symmetric about 90°. The interference coefficient 2(−1)^{2s}/(2s+1) is +2 for spin-0 bosons (α + α,
¹²C + ¹²C) and −1 for spin-½ fermions (p + p) **[I: standard result, not checked against a reference in this
session]**. EXFOR's `MOT` is this reference
(LEXFOR p. S.4), but see §3.2 finding 2. A ratio coded `RTH` for identical particles is ambiguous and must be
asked about, not guessed.

**θ → 0.** dσ_R/dΩ diverges as θ⁻⁴ and R → 1: the nuclear amplitude is negligible against Coulomb. Forward-angle
points carry the normalization of the data set (many experiments normalize there), and essentially no information on
θ; their uncertainties are the normalization uncertainty. A model ratio near θ = 0 is a ratio of two large numbers
and must be formed from amplitudes, R = |f_C + f_N|²/|f_C|² = 1 + 2 Re(f_C* f_N)/|f_C|² + |f_N|²/|f_C|², not from two
separately rounded cross sections **[I]**.

**Neutral channels.** η = 0: no Rutherford cross section; a ratio to Rutherford is undefined and must be refused. A
ratio between charged-particle elastic and inelastic groups is a RATIO of two observables, not a Rutherford ratio.

**Why fit the ratio rather than multiply the data back.** With fixed E, θ and no fold, (1) gives the same χ² for
(m/σ_R vs nR, error n δR) and (m vs n R σ_R, error n δR σ_R): the two are algebraically identical. They differ when

1. the energy moves (a fitted shift ΔE): the ratio form re-evaluates σ_R(E + ΔE, θ), the multiplied form keeps
   σ_R(E, θ) frozen into the data, so a shift of 1 % changes the converted cross section by about 2 % (σ_R ∝ E⁻²)
   and the χ² is computed against data that should have moved with the shift;
2. a fold acts (beam spread, target thickness, finite angular acceptance Δθ): the measured quantity is a ratio of
   averages, ⟨dσ/dΩ⟩/⟨dσ_R/dΩ⟩, not ⟨R⟩ and not ⟨dσ/dΩ⟩/σ_R(⟨θ⟩). Because σ_R varies as sin⁻⁴(θ/2), a uniform
   ±1° acceptance raises ⟨σ_R⟩ over σ_R(θ) by 3.4 % at 10°, 0.83 % at 20° and 0.36 % at 30° (numerical average of
   sin⁻⁴(θ/2), computed for this note); the model must fold numerator and denominator separately;
3. the masses or constants of the conversion differ from the engine's (mass numbers in `nds.py`), which moves the
   converted cross sections by a constant fraction per point;
4. the plotted and reported residuals of a ratio stay O(1) over the whole angular range, while multiplied data span
   several decades (forward-angle points dominate any plot and any badly scaled covariance).

**Apparatus.** Beam energy (and its spread), target thickness (energy loss), angular acceptance of the detector;
not the absolute charge or solid angle, which cancel in the ratio, and that is why such data are published as
ratios. The normalization n stays: it is the uncertainty of the experiment's own normalization of R (often to R = 1
at forward angles).

### 4.4 Legendre and cosine coefficients

dσ/dΩ_cm = Σ_L A_L P_L(cos θ_cm) with A₀ = σ/4π. AZURE2's normalized coefficient a_L = A_L/A₀ (`isDiff` 9) is EXFOR's
`LEG/RS` W_L; its `absolute 1` A_L is EXFOR's `LEG` a_L (b/sr). The other forms are fixed linear maps: `LEG/RSL`
B_L = W_L/(2L+1); `LEG/2L2` a_L = 2A_L/(2L+1); `LEG/4PI` a_L = 4πA_L; `LEG/RS0`, `LEG/RSD` divide by
Σ_L A_L P_L(1) or Σ_L A_L P_L(0) with the normalization of §3.2 finding 3; `LEG/1K2` multiplies by k² (energy
dependent). Cosine forms are a triangular linear map of the Legendre ones. Only c.m. coefficients (`NUMBER-CM`) are
a linear function of the c.m. A_L; lab-frame Legendre fits of particle distributions are not, and must be refused.
Elastic scattering of charged particles has no finite expansion; `LEG/RTH` coefficients of the Rutherford ratio are
the only form that can exist there, and they need the ratio of §4.3 projected on P_L **[I]**.

### 4.5 Polarization observables

A_y(θ) = Tr(U σ_y U†)/Tr(U U†) for a spin-½ beam, Madison convention (`docs/source/theory/polarization_theory.rst`);
frame independent at corresponding angles (LEXFOR p. C.1), so a lab angle only needs converting. The outgoing
polarization P of A(a,b)B equals A_y of B(b,a)A (time reversal; `segments.rst:177-186`). Tensor analyzing powers (T₂₀,
T₂₁, T₂₂, A_yy) for a spin-1 beam need the rank-2 spin trace: a new leaf, not a combination. The Legendre
coefficients of A_y (`ANA/LEG`) are coefficients of the product A_y dσ/dΩ in the usual convention, which LEXFOR's
description (`d-Sig/d-Omega = a(0) + …`) does not settle **[I]**; refuse until an entry's text says which.

### 4.6 Yield linked to a target (treated fully)

**Thin and thick target.** A beam of lab energy E₀ (MeV), Gaussian spread σ_b, loses energy in a target of areal
density N (active nuclei cm⁻²) with stopping cross section per active nucleus ε(E) (MeV cm²), effective for a compound
A_xB_y: ε_eff = ε_A + (y/x) ε_B (Bragg additivity; `experimental_effects.rst:263-300`). The number of reactions per
incident particle is

```
Y(E₀) = ∫_0^{E₀} dE′ ∫ dE g(E′; E₀) P(E, E′) σ(E) / ε(E)                                      (6)
```

with g the beam energy distribution and P(E, E′) the probability that a particle entering at E′ reaches the depth
where its mean energy is E, broadened by straggling σ_s²(x) = S(E(x))² ∫₀^x (dΩ²/dx′)/S(E(x′))² dx′ (Tschalär;
ARCHITECTURE §7). Limits:

- thin target, Δ = ε N ≪ the scale on which σ varies: Y = N σ(E₀ − Δ/2) · 10⁻²⁴ (σ in b);
- thick target, N → ∞ (the beam stops): Y(E₀) = 10⁻²⁴ ∫_0^{E₀} σ/ε dE, eq. (3);
- an isolated narrow resonance (Breit–Wigner, constant ε, λ and partial widths across it; all energies c.m.):

  ```
  Y(E₀) = (λ_r² ωγ / (2π ε_r)) [ arctan((E₀ − E_r)/(Γ/2)) − arctan((E₀ − Δ − E_r)/(Γ/2)) ] · 10⁻²⁴ (7)
  ```

  with λ_r = 2π/k_r the c.m. de Broglie wavelength at the resonance and ε_r the stopping in the c.m. frame,
  ε_cm = ε_lab m_A/(m_a + m_A). Its plateau for Δ ≫ Γ is the thick-target step height
  Y_∞ = λ_r² ωγ / (2 ε_r^cm), which is LEXFOR's ωγ = (2ε/λ²)(M_t/(M_t + M_i)) Y_thick with ε in the lab
  (LEXFOR p. S.11) **[V for the LEXFOR formula; eq. (7) is the standard Breit–Wigner integral, I]**.

**Units.** The model predicts reactions per incident particle (dimensionless) or per particle and steradian. Data
appear as (a) per incident particle (`PRT/INC`); (b) per unit charge, Y_µC = Y · 10⁻⁶ C/(q e) with q the charge state
of the beam (a ¹²C²⁺ beam has half as many particles per µC as ¹²C⁺); (c) per µA·h (×3600); (d) arbitrary units
(free n, no prior); (e) target-averaged cross sections, eq. (2), when the experiment has already divided by N;
(f) counts of one γ line, which needs the branching (a `Scale`) and the detector's efficiency and solid angle, which
belong to the data reduction and must already be removed, with their uncertainty in the normalization. The engine
must never guess q: a per-charge yield without q is refused.

**The target is a named object.** One physical target is often seen by several data sets: several γ lines or
detectors at one beam energy, an excitation function and angular distributions, a resonance scan and an off-resonance
run. They share the areal density, the stoichiometry (which sets both ε_eff and the number of active nuclei), the
stopping evaluation and its scale nuisance, the depth profile (implanted targets) and their priors. Today this is
`densitygroup`, `<abundances>` samples and `thicknessgroup` (§2.2), three mechanisms that do not know about each
other. A `Target` (ARCHITECTURE §3, `Material`) with a name, layers, stoichiometry, areal density (fitted or fixed,
prior), profile, and a stopping model with a scale nuisance, referred to by name from each data set, removes the
copies: the same N enters the yield kernel and the same stoichiometry enters ε_eff and the abundance weights.

**A target that changes during the measurement.** Carbon build-up, sputtering or loss of an implanted species make N
and the profile functions of the accumulated charge Q. Two representations: (i) **epochs**: the target has one
density per epoch (groups of runs), each a nuisance with a prior about the nominal and tied to the others by the
fixed stoichiometry, the data set naming `target T@epoch`; (ii) a **law** N(Q) = N₀ exp(−Q/Q₀) or N₀(1 − Q/Q₀) with Q
per point in a column of the data file. (i) needs no model of the degradation and no new column; (ii) is the physics when the reference runs (a resonance scan repeated during the campaign) are in the fit. Recommend (i)
first, (ii) when a project needs it.

**Thick-target step heights** are data with a beam energy range rather than a beam energy: Y_∞ of one resonance, eq.
(7) at E₀ − E_r ≫ Γ, minus the non-resonant yield. Entered as a yield datum at E₀ above the step, they are a
thick-target yield and need nothing new; entered as an ωγ, they are §4.7.

**Approximations.** Eq. (2) evaluates the beam spread and straggling with two-point Gauss–Legendre per sub-interval;
eq. (6) neglects multiple scattering and the angular spread of the beam, and the energy dependence of the detection
(angular distribution, Q-coefficients) across the target; screening (`screening U_e`) multiplies σ before the fold.
Each needs its tolerance stated in the Phase 5 kernel note.

### 4.7 Resonance strength, energy and widths (treated fully)

**ωγ from the level parameters.** For a level λ of spin J in the compound nucleus, entrance a + A and exit b + B,

```
ωγ = ω Γ_a Γ_b / Γ,    ω = (2J+1)(1 + δ_aA) / ((2j_a+1)(2j_A+1)),    Γ = Σ_c Γ_c                     (8)
```

(eV), with Γ_c the **observed** partial widths at the observed resonance energy E_r: Γ_c = 2 P_c(E_r) γ̃²_c, γ̃ the
observed (Brune) reduced width, P_c the penetrability at E_r, for a γ channel P_γ ∝ k_γ^{2L+1}. AZURE2's `<levels>`
already stores open-channel widths in this form ("partial width (eV) of an open channel",
`data_formats.rst:273-274`), so ωγ is an algebraic function of the parameters θ with an analytic derivative. EXFOR
`WID/STR` is ωγ; `WID/STR,,RG` is (2J+1)Γ_iΓ_j/Γ (without the spins of the entrance particles); `PAR,WID/STR` uses the
partial width to one final level, which in AZURE2 is one γ pair (Dict. 236; LEXFOR p. S.11) **[V]**.

**Formal versus observed widths.** The Lane–Thomas formal width Γ_c^f = 2 P_c(E_λ) γ²_λc is evaluated at the formal
energy with the boundary condition B_c; for one level the observed width is Γ_c = Γ_c^f / (1 + Σ_c γ²_c dS_c/dE|_{E_r}).
For several levels of one Jπ the map is a matrix transformation (Brune, Phys. Rev. C 66, 044611 (2002)), which
`level_parameters` performs (ARCHITECTURE §3, §5). An ωγ compared with formal widths is wrong by the factor
(1 + Σγ² dS/dE)², of order 10 % for broad levels **[I]**; it must always use observed widths at E_r.

**ωγ from the integral of the model.** For an isolated narrow resonance,

```
∫ σ_{ab}(E) dE = (λ_r²/2) ωγ = (2π²/k_r²) ωγ,                                                   (9)
```

which is how the experiment determines ωγ (eq. (7), or the area of a thin-target scan, LEXFOR p. S.11). Eq. (8) and the
model integral agree only when (a) Γ ≪ E_r, so λ², P_c and the γ-ray energy factor are constant across the
resonance; (b) no other level of the same Jπ interferes within a few Γ; (c) there is no non-resonant (direct capture,
background-pole, hard-sphere) amplitude in the same channel, whose interference changes the area; (d) the integration
window of the experiment and of the model are the same. Where these fail, ωγ is a property of a fit to the yield
curve, not a model-independent datum: the experiment's own definition (step height after subtracting a linear
non-resonant yield) is what must be computed, and that is the yield of §4.6 with a declared energy window. A
definition based on eq. (8) is then compared with a number that does not mean the same thing.

**Resonance energy and widths quoted by an experiment.** E_r appears as a lab beam energy, a c.m. energy or an
excitation energy (EXFOR `EN` with lab/`EN-CM` headings, or `E-LVL`); the comparison quantity is the observed energy
of the level (Brune's E_λ in AZURE2's observed basis), which is neither the S-matrix pole nor the peak of a broad,
interfering line. Γ quoted from a Breit–Wigner fit is the observed total width; Γ_p, Γ_α from an analysis assume a
penetrability and a radius that may differ from the project's — a reduced width γ² (`WID/RED`, `WID/RED,,RMT`) is
meaningful only with its channel radius and boundary condition, which EXFOR does not code (LEXFOR p. M.7 gives
Γ = 2γ²P only) **[V]**.

**Sub-threshold states.** E_r < 0 in the entrance channel: Γ_a = 0, ωγ is undefined; the entrance-channel strength
is an ANC (fm⁻¹ᐟ²), which AZURE2 stores for closed channels (`data_formats.rst:273-274`). An "ωγ" of a sub-threshold
level is refused; an ANC datum is the analogue (a constraint on γ̃ through the Whittaker function at a_c).

**Where it plugs in.** These data have no beam energy: they are terms of the objective on functions of θ, like the
measured scattering length (`include/ScatteringLength.h`; T9 in `objective.md`). The value is computed by
`level_parameters` from θ and its derivative is the chain rule through the Brune transformation. One `ObjectiveTerm`
per kind (resonance strength, energy, total width, partial width, branching ratio Γ_b/Γ, ANC), each
((f(θ) − d)/e)², in `src/Objective.cpp` as the recipe in `CLAUDE.md` ("A term of the objective") requires.

### 4.8 Reaction rates

N_A⟨σv⟩(T) is an integral of σ over a Maxwellian (LEXFOR p. T.6). A measured rate does not exist: an `SGV` value in
EXFOR is derived from cross sections or resonance parameters, often the same ones in the fit. Refuse as data;
compare rates as outputs.

## 5. Proposal

### 5.1 One explicit quantity per data set

Each data set states, in words, what its numbers are:

```
quantity <name>  unit <unit>  frame <lab|cm>  [reference rutherford|mott]  [target <name>[@epoch]]
         [per particle|charge q=<n>|arbitrary]  [norm free|fixed|prior <%>]
```

as keyword blocks on the segment line (parsed in `SegLine.h`, listed in `_KEYWORDS`, recipe "A segment or
extrapolation keyword" in `CLAUDE.md`). `io/` reads an old line without `quantity` and builds the same expression tree
from `isDiff`, the composite block, `sample`, `total`, `transmission`, `yield` and the `<targetInt>` index list; a line
with both a `quantity` and an `isDiff` that disagree is refused. The writer keeps old lines byte for byte when unedited
(`project.md`).

Unit strings are a closed list with conversion factors to the engine's units, checked at load (b, mb, µb, nb; b/sr …;
MeV b, keV b, eV b; deg, rad; eV, keV, MeV; fm; per-particle, per-µC, per-µAh; arbitrary); an unknown unit is an
error. The frame applies to E, θ and value together, except where the quantity is frame independent (ratios, A_y),
where it applies to E and θ only.

Named targets go in a `<targets>` section read by `Config::SectionStream("targets")` (recipe "A `<config>` key or a
new `.azr` section"):

```
<targets>
# name    layers / stoichiometry   active   N (atoms/cm2)  fit  prior   stopping            straggling
TaO       Ta2O5                    O-16     4.1e17         1    5%      azure2 scale 1 3%   bohr-chu
Ne22gas   Ne gas enr(22Ne=0.999)   Ne-22    2.0e17         0    0       azure2 scale 1 4%   bohr-chu
C13@run1  C 13C=0.99               C-13     1.0e18         1    10%     azure2 scale 1 3%   bohr-chu
C13@run2  same C13@run1                     8.5e17         1    10%
</targets>
```

The stoichiometry gives both ε_eff (Bragg) and the abundance weights of `Weighted`; the beam spread stays with the
beam (a property of the data set or of a `beam` record), not with the target. Today's `<targetInt>` lines remain
readable: each becomes an anonymous target bound to its segment index list.

### 5.2 Kind by kind

Columns: (a) does the engine express it today; (b) model expression (leaf ∘ nodes; folds written as kernels K);
(c) `.azr` form; (d) the test with an independent answer.

1. **Angle-integrated σ, partial σ to a level.** (a) yes, code 0. (b) `σ(α→β)`. (c) `quantity sigma unit mb`.
   (d) existing projects bit-identical; a mb file equals the b file ×10⁻³.
2. **S factor.** (a) no (pyazr converts). (b) `Scale_SF(E) ∘ σ`, S = E σ e^{2πη} with the engine's η. (c)
   `quantity sfactor unit keV*b frame cm`. (d) χ² of a fit to S data equals, to 10⁻¹², the χ² of the same data
   converted to σ by the engine's own formula; a narrow-resonance-free S(E) of a single hard-sphere channel against
   mpmath.
3. **dσ/dΩ, lab or c.m.** (a) yes, codes 1, 4. (b) `dσ/dΩ(α→β, θ)` with the Jacobian in the load-time conversion.
   (c) `quantity dsigma/domega unit mb/sr frame lab`. (d) the same data in lab and in c.m. (converted by an independent
   script) give the same χ² to 10⁻⁹; existing photon-aberration test.
4. **Ratio to Rutherford / Mott.** (a) no. (b) `Ratio(K[dσ/dΩ(α→α, θ)], K[dσ_C/dΩ(α, θ)])` with a new analytic leaf
   `Coulomb(α)`: Rutherford, or Mott for identical pairs (spin from the pair), folded by the same kernel K. A ratio of
   two folds, not a fold of a ratio. The leaf exists only for charged pairs; a neutral pair is refused. (c)
   `quantity ratio unit 1 reference rutherford frame cm` (`reference mott` required for identical pairs, `rutherford`
   refused there). (d) with every width set to zero the model is R = |f_C + f_HS|²/|f_C|², f_HS the hard-sphere
   amplitude, which a short mpmath script computes from Coulomb functions at a_c: agreement to 10⁻⁹; for α + α the
   same with eq. (5) as denominator; without shift or fold the χ² equals that of the multiplied data to 10⁻¹²; with a
   5-keV shift it differs by the amount the script computes.
5. **Legendre coefficients, every EXFOR normalization.** (a) partly (code 9: `LEG/RS` and `LEG`). (b) linear maps
   (`Scale`, `Ratio` of linear combinations) of the coefficient leaves A_L. (c)
   `quantity legendre unit 1 norm rs|rsl|2l2|4pi|abs frame cm`; `rs0`, `rsd`, `1k2` only after K3. (d) the identities
   of §4.4 on T2's (p,α₀) coefficients; lab-frame coefficients refused with a message.
6. **A_y.** (a) yes, code 7 (c.m. angles). (b) `A_y(α→β, θ)`, folded as numerator A_y σ and denominator σ
   (ARCHITECTURE §3). (c) `quantity analyzing-power frame lab|cm`; a lab angle is converted at load, the value
   unchanged. (d) existing polarization tests; a lab-angle file equals the c.m. file after conversion.
7. **Tensor analyzing powers, outgoing polarization.** (a) no / via inverse reaction. (b) a new leaf (rank-2 spin
   trace) — later; P via `A_y` of the inverse channel. (c) refused now with the advice the engine gives for code 8.
8. **Phase shifts, total, transmission, β-delayed, UPOS, THM, total capture, E1/E2 parts.** (a) yes. (b) the existing
   observables as leaves; E1/E2 as amplitude filters (ARCHITECTURE §3). (c) `quantity phase-shift J l s`,
   `quantity sigma exit all-gamma`, `filter E1`. (d) existing tests bit-identical (Phase 4 gate).
9. **Sums and ratios of reactions, natural targets.** (a) yes (SUM, RATIO, `sample`). (b) `Sum`, `Ratio`,
   `Weighted(abundances of the named target)`. (c) components named by reaction (`sum (p,a0) (p,a1)`) rather than by
   pair keys; `target natNe`. (d) `tests/natural_neon` bit-identical; a sum of two partial σ equals their separately
   fitted sum.
10. **Thin-target yield (target-averaged σ).** (a) yes (eq. (2), unit b). (b) `K_thin(target) ∘ σ`. (c)
    `quantity sigma unit mb target TaO` — the data are a cross section, averaged over the target. (d) constant σ
    survives the fold exactly; a linear σ gives σ(Ē) exactly (`TargetEffect.h:163-168`).
11. **Thick-target yield.** (a) per incident particle only (eq. (3)). (b) `K_thick(target) ∘ σ`, then `Scale(1/(q e))`
    for per-charge data and `Scale(branching)` for one γ line. (c) `quantity yield unit 1/uC per charge q=2 target
    C13@run1`; `per arbitrary` = free norm without prior. (d) σ = σ₀ constant and ε = a/E (closed form
    Y = 10⁻²⁴ σ₀ (E₀² − E_min²)/(2a), E_min = 1 keV the engine's lower limit); a narrow Breit–Wigner against
    eq. (7) with its stated O(Γ/Δ) and O(Δ/E) corrections;
    a per-µC file with q = 2 equals the per-particle file × 10⁻⁶/(2e).
12. **Shared and evolving targets.** (a) partly (`densitygroup`, `thicknessgroup`, samples). (b) one `Target`
    object whose parameters (N per epoch, stoichiometry, stopping scale) are nuisances ν of every kernel that names
    it. (c) `<targets>` above. (d) two segments naming one target give the χ² and gradient of today's `densitygroup`
    to bit identity; an epoch split with equal densities equals the unsplit target.
13. **ωγ, Γ, Γ_b/Γ, E_r (observed), ANC.** (a) no (only a₀). (b) objective terms on functions of θ from
    `level_parameters`, like T9. (c) a `<levelData>` section: `strength level <id> entrance <pair> exit <pair> value
    eV sigma eV`, `energy level <id> frame cm|ex|lab value σ`, `width level <id> [channel <pair>] value σ`,
    `branching …`, `anc level <id> pair <pair> s l value σ`. A warning at load when the level's Γ/E_r > 1 % or another
    level of the same Jπ lies within 10 Γ (the conditions under which eq. (8) and eq. (9) differ); refused for
    E_r < 0. (d) single isolated level: eq. (8) against the numerical integral (9) of the model σ to O(Γ/E_r); the
    gradient against finite differences; Brune ↔ formal round trip leaves the term unchanged.
14. **Relative data (`REL`, arbitrary units).** (a) yes, by a free norm. (b) `Scale(n)` with n free and no prior.
    (c) `norm free`; a prior on such a norm is refused. (d) χ² invariant under multiplying the file by 10⁶.

### 5.3 Priority

1. **The named quantity with unit and frame (5.1), read from today's codes.** It is the Phase 4 item itself
   ("Model expressions … today's segment types read into them"), it removes the two tables of `isDiff`, and every other
   kind below is a new value of `quantity`. Gate: bit-identical.
2. **Ratio to Rutherford/Mott.** The maintainer asked for it; about a fifth of the elastic angular distributions in
   the survey are Rutherford ratios (286 `RTH` against 1055 `DA` and 34 `DA … REL` elastic data sets, §3.3); today it can only be multiplied in, which is wrong under shifts and folds (§4.3) and uses
   mass numbers. Needs one analytic leaf and the existing `Ratio`.
3. **S factor as input.** Frequent (§3.3), trivial (a diagonal `Scale`), and removes the second place where η is
   computed (pyazr), which today uses mass numbers and is off by several per cent at low energy (§4.1).
4. **Named targets and thick-target yields with explicit units (per particle, per charge with q, arbitrary).** The
   maintainer's first example. It needs the `Target` object of ARCHITECTURE §3 and fits the Phase 5 kernel work;
   the `.azr` grammar and the object should be designed in Phase 4 so Phase 5 only fills the kernel.
5. **Resonance-parameter data (ωγ, E_r, Γ, branching, ANC).** Independent of the data path (objective terms), so it can
   be built in parallel; it needs `level_parameters` (Phase 1/2 work) to give observed widths with derivatives.
6. **Legendre normalizations beyond `RS`/`LEG`; relative data made explicit.** Small, after 1.
7. **Tensor analyzing powers, activation yields (with decay data), target drift laws.** When a project needs them.

### 5.4 Refused, with a message

- `SGV` (rates: derived, not measured) — "compare rates as outputs".
- `TTY…` activation yields (decay rate per charge; depend on half-life, irradiation history, and the product, not a
  channel) — until a decay-data node exists; message names `MLT…,TT` / `PY…,TT` as the prompt-yield codes.
- Process `X` (production of a nuclide or particle from all channels), `SCT`, `NON`, and sums whose competing reactions
  are not listed (LEXFOR p. S.23, "Sum reactions … are not coded using the form above").
- Spectrum-averaged quantities (`MXW`, `SPA`, `AV`, `BRA`, `FIS`), `RAW`, `MSC`, `FCT` without a stated factor.
- Lab-frame Legendre or cosine coefficients of particle distributions; `LEG/RSD`, `LEG/RS0` until the normalization is
  settled per entry (K3); `ANA/LEG`.
- `WID/RED` and `WID/RED,,RMT(/AMP)` without the channel radius and boundary condition in the entry's text.
- A Rutherford ratio for a neutral pair; `RTH` for identical particles (ask for `MOT`); `MOT` whose text means the
  relativistic point-nucleus cross section.
- ωγ of a sub-threshold level.
- Per-charge yields without the charge state.
- Double-differential, spectra (`SPC`), multiplicities of unidentified particles, `DA/DA` correlations.

## 6. What an EXFOR reader must record

The reader lives in pyazr `sources/` (ARCHITECTURE §3, `io/`) and writes the data file plus the `quantity` line. For
every data set it records:

1. **Identity and provenance**: dataset ID with pointer (`O2599004.1`), EXFOR entry version or the master-file commit
   and checksum, the full `REACTION` string including any combination and tautology, `REFERENCE`, DOI, `AUTHOR`, year,
   `STATUS` (Dict. 16: `DEP` dependent data, which must not be fitted together with its parent; `SPSDD` superseded or
   withdrawn).
2. **Quantity**: SF5–SF8 parsed against Dict. 236 (the dictionary file and its `trans` number stored), mapped to one
   `quantity` value or refused (§5.4); SF7 particle; `E-LVL`/`E-EXC` or γ energy `E` for partial data, mapped to an
   exit pair, with refusal if no pair matches.
3. **Reaction**: target and projectile Z, A from the code; natural target (A = 0) mapped to a named target with
   natural abundances; inverse kinematics (tautology) resolved to the forward form, with the experiment's form stored.
4. **Units and frame**: each column's heading and unit from Dict. 24 and Dict. 25, converted with the dictionary's
   factor (not by prefix matching), `-CM` headings recorded as frame; energy as lab or c.m.; masses from AME2020
   (`pyazr.masses`), never mass numbers.
5. **Errors**: every error column with its heading and the `ERR-ANALYS` text and correlation flag. Mapping:
   point error e_i = quadrature sum of `ERR-S` and partial errors flagged U; normalization uncertainty = quadrature sum
   of partials flagged F (and `ERR-SYS` when no partials are given and the text says systematic); partials flagged P or
   C go to a covariance file as `relvector` terms; `ERR-T` alone becomes e_i with no normalization uncertainty, and a
   warning; `DATA-ERR` of unknown type likewise. An error column that is missing for a point is an error (no default),
   unless the user states a value.
6. **Normalization status**: `REL`/`ARB-UNITS` → `norm free`; `MONIT`/`MONIT-REF` stored, and data normalized to a
   monitor flagged so the user can renormalize to a newer monitor value.
7. **Conversions made**: each one listed in a header comment of the written file (factor, from which column), so the
   file can be re-derived.

Refused with the reason: every code in §5.4, a quantity not in Dict. 236, a unit not in Dict. 25, a `DEP` data set
whose parent is also imported, a data set with only a `MISC` energy.

## 7. Open questions for the maintainer

**Director's working choices, 2026-11-09 (the recommendations below, taken so that implementation can proceed; each is
what users write or what a number means, so each stays open to the maintainer until the format is released):** K1
keyword blocks on the segment line; K2 refuse the ambiguous case; K3 refuse `LEG/RSD`, `LEG/RS0` in the reader; K4 the
level-parameter form with the load-time warning; K5 epochs; K6 the charge state with the data set; K7 `--convert`
rewrites `<targetInt>` once `<targets>` exists; K8 the stated error mapping with an override; K9 refuse `TTY`.
Model-only data go in one section, `<modelData>` (`data_model.md` M10), not a separate `<levelData>`: the entries of
§5.2 item 13 plus the scattering lengths. The pyazr S-factor conversion with mass numbers (§4.1, finding of §2.5) is
a defect of the shipped importer and is fixed with priority 3 (S factor as an engine quantity).

- **K1.** Should the named quantity live on the segment line (keyword blocks) or should data sets move to a new
  section with one record per measurement (quantity, file, target, normalization, shift, covariance)?
  Recommendation: keyword blocks now (old files and tools keep working, one parser), a new section only if Phase 9's
  web UI needs it.
- **K2.** For identical particles, is a `RTH` ratio refused or interpreted as Mott? Recommendation: refuse, ask the
  user to state `reference mott` or `reference rutherford`, since both occur in the literature.
- **K3.** `LEG/RSD` and `LEG/RS0`: accept with the normalization the user states, or refuse? Recommendation: refuse
  in the EXFOR reader; accept in the `.azr` with an explicit `norm rsd-1plus|rsd-sum`.
- **K4.** ωγ: compare the parameter form (8) only, or also offer the integral form (9) over a stated window?
  Recommendation: (8) as the datum, with the load-time warning; for broad or interfering levels the user fits the
  yield curve instead. No integral form.
- **K5.** Target that changes during a measurement: epochs only, or also a drift law in the accumulated charge?
  Recommendation: epochs now; a law when a project with repeated resonance scans asks for it.
- **K6.** Charge state of per-charge yields: a property of the data set (`per charge q=2`) or of a `beam` record
  shared by data sets? Recommendation: of the data set; one number, stated where the unit is.
- **K7.** Should `<targetInt>` lines be converted to `<targets>` by `--convert` (format 2), or stay as they are?
  Recommendation: convert, with the round trip test; the anonymous-target reading stays for old files.
- **K8.** EXFOR errors: is the U/F/P mapping of §6 acceptable as the default, with ERR-T-only data sets entering
  without a normalization uncertainty? Recommendation: yes, with the warning and the user's override.
- **K9.** Activation (`TTY`) yields: refuse permanently, or plan a decay node (λ, irradiation history)?
  Recommendation: refuse; they are rarely used for R-matrix fits of light nuclei (§3.3).

## 8. Not verified

- The EXFOR Formats Manual (IAEA-NDS-207) could not be downloaded (403); the correlation flags come from NRDC working
  papers wp2012-28 and wp2021-08, the heading rules from LEXFOR 2015/08. LEXFOR is from 2015; Dict. 236 is current
  (trans 9133). Where LEXFOR and the dictionary disagree (§3.2), both are quoted.
- Eq. (5) (Mott, with the spin factor), eq. (7) and the formal-to-observed width relation for one level are standard
  results written here without a check against a reference text in this session.
- The occurrence counts (§3.3) rely on `x4list`'s reaction codes and its `enMin`; data sets with energies in EN-CM may
  be filtered on a c.m. value; a data set listed under two targets is counted once.
- The β-delayed normalization (§2.1) and the lab-angle A_y import (§2.5 item 5) were read in the code, not run.
- No test project was run with the built engine.
