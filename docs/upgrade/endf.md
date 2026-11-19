# ENDF for the evaluators: design, with a measured spike

This note designs what the PI asked for on 2026-11-14 ("an ENDF file reader or exporter ... for the evaluators ... if
the ENDF utilities codes are installed, it could be able to launch it to check the ENDF file and correct it ... an
easy way of checking the evaluations"): a pyazr module that reads and writes ENDF-6 and runs the checking codes when
they are found, and screens of the page over it. Nothing of it is in the tree yet. Every statement about a record
cites the ENDF-6 formats manual, ENDF-102 (2023 edition, BNL, 431 pages, fetched from the NNDC), by its **printed**
page number. Every number of section 6 comes from a script of `work/endf/scripts/` (outside the repository; list in
section 6.4), run on the base build `work/b-trial71` (start commit 9e209ace) at one thread. "Not checked" says that
neither the manual nor a run shows it.

## 0. In one page

- **What maps.** An AZURE2 nucleus is an R-matrix in the sense of ENDF's "R-Matrix Limited" format (MF=2, MT=151,
  LRF=7, KRM=4: full R-matrix, every channel on an equal footing; pp. 71-77, 354-364). Pairs, spin groups, particle
  channels (l, s, radius), level energies and reduced-width amplitudes go into it after three conversions the
  engine already does or nearly does: Brune (or Park) to formal parameters with the engine's boundary values
  B_c = S_c(E_1); E_x to the lab energy of the incident particle; gamma^2 scaled by the same lab/c.m. factor, so that
  the R-matrix is unchanged. A project run with `--use-rmc` maps to KRM=3 (Reich-Moore) exactly.
- **What does not.** External capture, the photon channels' energy dependence E_gamma^(2L+1) and the rule that
  radiative widths stay out of the total width have no place in MF=2 (an LRF=7 photon channel has P = 1 and enters
  the level matrix, p. 356, p. 74). The design writes the capture cross sections pointwise in MF=3/MT=102 and the
  primary photons in MF=6 from the engine, and keeps MF=2 to the particle channels: exact for AZURE2's model. THM
  records, target effects, normalizations, the data and the hybrid potential have no ENDF place (they stay in the
  project; the MF=1 text names them). ENDF's background R-matrix (KBK>0), tabulated phases (KPS>0), separate
  effective and true radii, SHF=0 and fission channels cannot be read into AZURE2: refused with what, allowed, fix.
- **Measured: digits.** T1 and T2 written as LRF=7 text and read back into the engine: the classic ENDF number form
  (7 significant digits) moves chi^2 by -0.031 of 482 (T1) and -0.69 of 10905 (T2), and the model by up to 4e-4;
  the form with the letter E (6 digits) by +0.17 and +4.4; the widest form 11 columns hold (fixed point where it
  fits, as FUDGE writes; 8-10 digits) by 1.2e-5 and 4.4e-4, the model by 3e-7 and 4e-8. The separation energy must
  not come from the file's masses: T1's photon pair gives S_p 1.1 keV wrong in the classic form, chi^2 +710.
- **Measured: a real evaluation.** ENDF/B-VIII.1 3He + 4He (LLNL/INDEN 2022, LRF=7 KRM=4, boundary B_c = -l) read by
  the spike, its B changed to the engine's, put into a format-3 project: the engine reproduces the file's own
  MF=3/MT=600 (4He(3He,p)6Li, 7.2-20 MeV) to 1.1e-6 (median 7e-8). Read the other way (widths taken as c.m.), the
  same file is off by up to 143 %: the lab reading of ER and GAM is confirmed by the file itself. The same run found
  three things a checker should report: the file's 6Li mass is 21.7 keV below AME2020 while its Q value agrees with
  AME, the MF=3 points just above threshold are 37 orders above what its own MF=2 gives, and the header text says
  LTP=2 where the record holds LTP=1.
- **The checking codes** (CHECKR, FIZCON, PSYCHE, INTER, STANEF; PREPRO's LINEAR, RECENT, FIXUP, DICTIN; none
  installed here) were designed from their manuals and read in their sources: for LRF=7, CHECKR tests three flags
  and the list lengths, FIZCON and PSYCHE read the records and skip them (FIZCON's charge test is commented out).
  "Correct" from them can only mean the bookkeeping STANEF, DICTIN and FIXUP rewrite, and STANEF's standard form
  rounds every number it touches to 7 digits. The honest check of an evaluation's resonance parameters is the
  engine's reconstruction compared with its MF=3, as the spike did, and RECENT's when it is found.
- **Where it lives.** `pyazr.endf`: records and number forms, interpolation laws, the format's own laws (MF=3, MF=6
  LAW=2 and LAW=5), the LRF=7 frame and unit conversions, the checkers. Two doors in `api/src/`: the formal levels of
  a session with their B_c, and formal levels with a caller's B_c to the observed form. No physics in Python.
- **Order.** Six stages, each mergeable; the first (an evaluation's curves over the fit and the data) needs no
  engine change and answers "check the evaluations" on its own.

## 1. What was read, fetched and run

| what | where (outside the repository) | how |
|---|---|---|
| ENDF-102 (2023), 431 pp. | `work/endf/manual/ENDF-102-2023.pdf`, text in `work/endf/docs/endf102.txt` | NNDC `endfdocs/`, read with `pdftotext` |
| ENDF Utility Codes 7.01/02 manual (Dunford 2005, 17 pp.) | `work/endf/manual_util/endf_util7p01.pdf` | IAEA NDS |
| CHECKR 8.24, FIZCON 8.20, PSYCHE 8.11, STANEF 8.10: sources, sample decks and reports | `work/endf/util_src/` | github.com/IAEA-NDS/ENDF-utility-codes; read, not built |
| PREPRO 2023: README; RECENT, FIXUP, DICTIN, LINEAR, ENDF2C sources | `work/endf/manual_prepro/`, `work/endf/prepro_src/` | github.com/IAEA-NDS/PREPRO; read, not built (the IAEA's PREPRO pages answer 403 to scripts) |
| ENDF/B-VIII.1 release metadata; material lists, READMEs, change logs of the p, d, t, 3He, alpha sublibraries | `work/endf/manual/meta81.json`, `work/endf/listings/` | NNDC `endf-releases/releases/B-VIII.1/` |
| one evaluation: ENDF/B-VIII.1 3He + 4He, MAT 228 (and the archive's LICENSE) | `work/endf/evals/h_He4/` | the helium3s archive (190 109 bytes, MD5 as published), one file extracted, archive deleted |
| the spike's scripts and outputs | `work/endf/scripts/`, `work/endf/run/` | section 6.4 |

The IAEA mirror named in `endf-comparison.md` (`www-nds.iaea.org/public/download-endf/`) answered 403 behind a
browser challenge; the NNDC served everything. Nothing was installed, built or run of what was fetched.

## 2. The format, from the manual

### 2.1 Records, numbers, units

- Every number sits in an 11-column field. A float is an optional sign, digits with a decimal point, and an optional
  exponent written either as a sign and digits or as `e`/`E`, sign and digits; blank is zero; `d`/`D` is not allowed
  (p. 31). The manual adds that ENDF floats are "mid-way between single and double precision". Three forms follow:
  the classic `-1.234567+6` (7 significant digits; 6 with a two-digit exponent), the E form `-1.23456E+6` (6; 5), and
  a fixed-point form where it fits (`-9546287.6`, `1159.42866`: up to 10 digits), which is what FUDGE wrote into the
  3He + 4He file (section 6.2).
- Units: energies in eV, cross sections in b, masses in neutron masses, angles as cosines (Table 15, p. 23). The
  neutron mass is 1.00866491595 u (Appendix H, p. 394).
- Masses: AWR of a material and of heavy products in atomic units, the projectile AWI and light products in nuclear
  units (section 0.5.1.1, p. 24). The LRF=7 pair masses MA, MB are "in units of neutron mass" with no statement of
  atomic or nuclear (p. 72); the 3He + 4He file writes atomic ones (3He 3.0160301 u, p 1.0078250 u) and a nuclear AWI
  in MF=1 (2.989033 m_n = 3.01493 u). Q values are from atomic masses (section 3.3.2, p. 105).
- Interpolation laws INT 1-5 (constant, lin-lin, lin-log, log-lin, log-log) and INT=6, the charged-particle law
  sigma = (A/E) exp(-B/sqrt(E - T)) between two points, T = 0 for Q > 0 and the threshold otherwise (Table 16 and
  eqs. 1-6, pp. 25-26).
- A material: MF=1/MT=451 header and directory (pp. 38-42), then files in MF, MT order, each section HEAD ... SEND
  (p. 32).

### 2.2 MF=2, MT=151 and LRF=7

- File 2 is required only for incident neutrons; LRP in MF=1 says whether File 2 is given and used: -1 none,
  0 radius only, 1 parameters whose cross sections are added to File 3, 2 parameters given but not to be added
  (pp. 38-39, 59-60); "the LRP=2 option is not allowed in primary evaluations" (p. 42).
- LRF=7 "contains all the generality of LRF=3 plus unlimited numbers and types of channels" (p. 62); competitive
  reactions are channels of their own (p. 60).
- The records (pp. 73-74): a CONT `0, 0, IFG, KRM, NJS, KRL`; a LIST of the NPP particle pairs, 12 numbers each:
  MA, MB (masses), ZA, ZB (charges), IA, IB (spins; the sign of IA carries a parity), Q (eV), PNT (penetrability
  computed or not; default -1 for capture), SHF (shift computed or not), MT, PA, PB; then per spin group a LIST
  `AJ, PJ, KBK, KPS, 6*NCH, NCH` with six numbers per channel: PPI (pair), L, SCH (channel spin), BND (boundary
  condition), APE (radius for the phase), APT (radius for P and S), both in 1e-12 cm; and a LIST `0, 0, 0, NRS,
  6*NX, NX` with, per level, ER and NCH values GAM (pp. 72-74).
- KRM: 1 SLBW, 2 MLBW, 3 Reich-Moore, 4 full R-matrix ("identical to ... D.1.7 with the elimination of the imaginary
  term in the denominator ... and the inclusion of each gamma-channel on an equal basis with all other channels",
  p. 74). KRL: 0 non-relativistic, 1 relativistic; the manual gives the flag only (p. 71), no equations.
- IFG=1: GAM is the reduced-width amplitude in eV^1/2; IFG=0: the width in eV at E_lambda, with Gamma = 2 gamma^2 P
  (eq. D.45), at |E| for a negative-energy level, the sign of GAM the sign of gamma (p. 73).
- ER is "in eV" (p. 73) and, for every resolved format, "in the laboratory system" (p. 66). The R-matrix is
  R_cc' = sum gamma gamma / (E_lambda - E) (eq. D.44, p. 356) with E the energy of eq. D.39,
  k^2 = 2 m_a m_b^2 E / (m_a + m_b)^2 (p. 355), the lab energy of the incident particle; for a charged pair
  hbar^2 k^2 = 2 mu (E M/(m+M) + Q) (eq. D.82, p. 364). So gamma^2 in an LRF=7 file is on the lab energy scale: the
  spike confirms this on a library file (section 6.2).
- Optional extensions (pp. 74-77): a background R-matrix for KBK channels (LBK=0 none, 1 tabulated complex, 2 SAMMY's
  logarithmic form, 3 Froehner's form; "currently only allowed for KRM=3", p. 74; eqs. D.76-D.79, p. 363) and
  tabulated phase shifts (KPS=1). The incident channel is found by MA = 1 with Q = 0, or by MT = 2 (p. 71). For
  KRM = 1-3 the first pair is gamma + compound nucleus (p. 73).
- Formulas: U = Omega W Omega, W = P^1/2 (1 - RL)^-1 (1 - RL*) P^-1/2, L = S - B + iP (eqs. D.40-D.43, p. 356); for
  "non-eliminated capture channels, the penetrability is unity" (p. 356); P, S and phi of charged pairs from the
  Coulomb functions (eqs. D.80-D.85, p. 364); the angle-integrated cross sections D.49-D.53 are stated as "accurate
  only for cases in which one of particles in alpha is a neutron" (p. 358).
- Range record: EL, EH, LRU, LRF, NRO, NAPS (pp. 63-65); NRO=1 an energy-dependent radius table.

### 2.3 The other sections an evaluator hands over

- **MF=1/MT=451** (pp. 38-42): ZA, AWR, LRP, LFI, NLIB, NMOD; ELIS, STA, LIS, LISO, NFOR=6; AWI, EMAX, LREL, NSUB
  (10010 protons, 10020 deuterons, 10030 tritons, 20030 helions, 20040 alphas: Table 3, p. 7), NVER; TEMP, LDRV,
  NWD, NXC; three standard text lines (ZSYMAM, ALAB, EDATE, AUTH; REF, DDATE, RDATE, ENDATE; HSUB) and free text
  (66 columns); the directory MF, MT, NC (lines), MOD.
- **MF=3** (pp. 103-111): one TAB1 per MT, E the lab energy of the incident particle, QM and QI (atomic masses).
  Incident charged particles: no total cross section (p. 110); elastic MT=2 holds 1.0 or the nuclear-plus-
  interference integral sigma_NI for a cut-off angle (p. 110, and p. 150); inelastic MT 601-649 (protons out) ...
  801-849 (alphas), totals MT 103-107 (p. 111).
- **MF=4** is for incident neutrons and photons only (p. 112); a charged-particle file gives angular distributions
  in **MF=6**: LCT (frame), product subsections with yields and a LAW (pp. 133-137); LAW=2 two-body Legendre
  (LANG=0) or tabulated (12, 14) coefficients in the c.m., only for MT 50-90, 600-648, ... and for primary photons of
  MT=102 (pp. 142-143); LAW=5 charged-particle elastic: LTP=1 nuclear amplitude expansion (complex a_l, real b_l,
  eqs. 6.13-6.14), LTP=2 residual Legendre expansion (6.15-6.18), LTP=12, 14, 15 tabulated P_NI(mu) with sigma_NI in
  MF=3 (6.19-6.20); Rutherford and Mott terms eqs. 6.9-6.12 (pp. 144-147). "Whenever possible, the nuclear amplitude
  expansion should be used with LAW=5" (p. 150). (The 2023 text calls both LTP=12 and LTP=15 "linear in mu", pp. 146-147.)
- **Photon production**: MF=12 multiplicities (LO=1; LP=2 marks a primary photon) or transition probability arrays
  (LO=2) (pp. 199-203); MF=13 photon production cross sections (p. 205); MF=14 angular distributions and MF=15
  spectra (table of contents only, pp. 209-216); or MF=6 MT=102 with a primary photon to each final state as LAW=2,
  its recoil LAW=4, then the cascade (ordering rule, p. 136).
- **MF=32** (pp. 251-259, 267-268): covariances of the File 2 parameters, LCOMP=1 for LRF=7: blocks of spin groups,
  per group NCH, NRB, the parameters ER, GAM, then the upper triangle of the covariance of the (NCH+1) NRB
  parameters; LCOMP=2 a compact form with NDIGIT-digit correlations. Covariances are on the parameters as File 2
  lists them (GAM with its sign, p. 253). **MF=33** (pp. 274-276): covariances of File 3 cross sections by energy
  groups (NC- and NI-type sub-subsections).

## 3. The mapping

### 3.1 AZURE2 project → LRF=7, quantity by quantity

| AZURE2 | LRF=7 field | conversion | engine function that does it today |
|---|---|---|---|
| particle pair: masses (atomic u) | MA, MB | / m_n (p. 394) | none needed (pair records, `Session.pairs`) |
| charges, spins, parities | ZA, ZB, IA, IB, PA, PB | IA signed by the light particle's parity | none needed |
| separation energy S, residual excitation | QI of each pair, relative to the entrance pair | Q = S_ent - S_pair - E_x,residual (eV) | none needed |
| light-particle excitation (column 32) | QI and IA of the excited light particle | as above | none needed |
| entrance pair | MT=2 | the incident channel is found by MT=2 (p. 71) | `Pair.is_entrance` |
| exit particle pairs | MT=600+k (p), 650+k (d), 700+k (t), 750+k (3He), 800+k (alpha) | k = order of the residual's excitation | none needed |
| photon pair (final state J^pi, E_x) | a pair MT=102, PNT=-1, IB, PB the final state, Q = S_ent - E_x,f | only if photon channels go into MF=2 (decision 1) | none needed |
| beta-decay pairs (F, GT) | none | refused | — |
| J^pi group | AJ (sign = parity), PJ | — | parameter records (`J`, `parity`) |
| particle channel l, s | L, SCH | — | parameter records (`L`, `S`) |
| channel radius (per pair) | APE = APT (per channel) | fm / 10 | pair record `channel_radius` |
| boundary condition: Brune and Park need none; formal: B_c = S_c(E_1), E_1 the lowest level of the group | BND, SHF=1 | written as the engine's formal B_c | `boundary_values` (physics/), not exposed: door 1 (section 8.2) |
| level energy (Brune/Park: observed E_x) | ER | formal E_lambda, then (E_lambda - S_ent) (m_A + m_B)/m_B x 1e6 eV | `formal_from_brune` (physics/), reached today only by a session with `use_brune=False` (`Parameters.x0`) |
| widths, ANCs, amplitudes (any input kind, `width_is_amplitude`) | GAM, IFG=1 | formal gamma, x sqrt((m_A + m_B)/m_B) x 1e3 eV^1/2 | as above (TransformIn, then formal_from_brune) |
| Park amplitudes | GAM | Park -> Brune -> formal | `brune_from_park`, `formal_from_brune` |
| sub-threshold level given by an ANC | negative ER, GAM (IFG=1) | as any level; IFG=1 avoids the width-at-|E| rule of IFG=0 | as above |
| background poles | ordinary levels (the "dummy resonances" of LBK=0, p. 75) | — | — |
| photon channels, `gamma_decay` records (Gamma_gamma, branching, mixing inputs) | channels of the photon pair (decision 1) | gamma_ENDF = gamma_formal (E_gamma/hbar c)^(L+1/2) at E_lambda (P = 1 in LRF=7, p. 356), lab-scaled; E or M follows from the parities and L | `photon_penetrability` (physics row) |
| radiative widths outside the total width (the PI's rule) | none: KRM=4 puts photon channels in the level matrix, KRM=3 puts Gamma_gamma in the denominator | written in MF=3/MF=6 instead (decision 1) | — |
| external capture | none in MF=2 | only in MF=3/MT=102 and MF=6 (computed) | the engine's capture cross sections |
| `--use-rmc` (Reich-Moore, P = 1, -i sum gamma^2 on the diagonal) | KRM=3, the photon pair first, its channels the eliminated widths | exact | — |
| relativistic=1 | KRL=1 | the engine follows EDA's prescription (LA-UR-19-31084); ENDF-102 gives no equations for KRL=1 (p. 71); `relativistic.out` maps level parameters exactly at E_lambda only | decision 4 |
| several compound nuclei / a natural target | NIS isotopes (ZAI, ABN) in one MF=2 | one LRF=7 subsection per isotope, ABN from the abundance | not checked further |
| identical particles (alpha + alpha, 12C + 12C) | nothing in MF=2; LIDP=1 in MF=6 LAW=5 | named in the MF=1 text | — |
| hybrid nuclear potential (`potential` record) | none (KPS=1 would carry phases, not P and S) | refused | — |
| THM records, targets, beams, detectors, normalizations, energy shifts, data sets | none | not written; named in the MF=1 text with chi^2 per data set | — |

**Summary.** *Maps cleanly*: pairs, spins, parities, charges, masses, Q values, J^pi groups, channels, radii,
background poles, Reich-Moore projects. *Converts* (engine functions exist): the basis (Brune or Park to formal, with
B_c), the frame (E_x to lab ER, gamma^2 to the lab scale), the input kinds (widths, ANCs, amplitudes), photon
amplitudes to P = 1 at E_lambda. *Lost with a statement*: external capture and the photon energy dependence in MF=2
(kept exact in MF=3/MF=6), the radiative-width rule if photon channels go into MF=2, the apparatus and the data,
identical-particle symmetry in MF=2, KRL=1 away from the levels. *Refused*: beta-decay channels, hybrid potentials, a
Brune set with no formal equivalent (`M` not positive definite: the engine's own refusal, `tests/7Li_p_ay`).

### 3.2 LRF=7 → AZURE2, what the reader does

| LRF=7 holds | reader |
|---|---|
| KRM=4, IFG=1 or 0, KBK=0, KPS=0, SHF=1, APE=APT, one radius per pair | converts: ER, GAM to E_x and c.m. gamma; BND to the engine's B (door 2); a project in the Brune basis |
| KRM=3 (eliminated capture) | converts to `use_rmc=1`: one photon pair carries the eliminated widths (only their sum matters) |
| KRM=1, 2 (SLBW, MLBW inside LRF=7) | refused |
| KBK>0: background R-matrix (LBK 1, 2, 3) | refused: what (the option), allowed (KBK=0), fix (distant poles in the evaluation, or SAMMY) |
| KPS>0: tabulated phase shifts | refused |
| APE != APT, or channels of one pair with different radii, NRO=1 | refused (AZURE2 has one radius per pair for P, S and phi) |
| SHF=0 on a particle channel, PNT=-1 on a particle channel | refused |
| a fission pair (MT=18) | refused |
| KRL=1 | decision 4 |
| NIS>1, NER>1 | the user picks the isotope and range; one project each |
| the compound nucleus' separation energies | not in the file: from AME2020 (`pyazr.sources.masses`), the file's Q values kept, the file's masses for the kinematics; the difference is reported (section 6.2) |

## 4. The other sections: what the engine fills, what the evaluator adds

| section | filled from | grid and frame | the evaluator adds |
|---|---|---|---|
| MF=1/MT=451 | the project (pairs, levels, options), the run (chi^2 per data set, norms), `pyazr.describe()` | — | MAT, NLIB, LREL, NVER, ALAB, AUTH, REF, dates, EMAX; the evaluation's story |
| MF=2/MT=151 | formal levels (door 1) | lab ER, eV | EL, EH (default: the data's range) |
| MF=3 MT=2 | 1.0 (LTP=1, 2) or sigma_NI (LTP=12-15) from the engine's dsigma/dOmega minus Rutherford, integrated to mu_max | lab E; an adaptive extrapolation grid (`energies adaptive`) to a tolerance, INT=2 | mu_max |
| MF=3 MT=600.., 800.. | the engine's `sigma` of each particle exit | as above | — |
| MF=3 MT=102 | the engine's capture `sigma`, summed over final states, external capture included | as above | — |
| MF=6 MT=2 LAW=5 | LTP=12 first: P_NI(mu) from the engine's dsigma/dOmega on an angle grid; LTP=1 needs the nuclear amplitudes (a new door, stage 5) | c.m. mu, lab E | LIDP for identical particles is set from the pair |
| MF=6 MT=600.. LAW=2 | the engine's Legendre coefficients (`legendre`, normalised) | c.m. (LCT=2, p. 143) | — |
| MF=6 MT=102 | primary photon per final state: yield sigma_f/sigma, LAW=2 angular distribution from the engine; recoil LAW=4 | c.m. | the cascade below each final state (ENSDF), MF=12 LO=2 if preferred |
| MF=32 LCOMP=1 | the covariance of the formal parameters: the fit's covariance carried through the transformation at the written B_c, J C J^T (section 13.13), scaled to eV and eV^1/2 | — | which levels to include |
| MF=33 | sensitivities and the parameter covariance on an energy grid (`pyazr.bands`), NI-type sub-subsections | groups | the long-range (normalisation) components |

The engine computes every number on a grid it is given; the module writes the TAB1 and LIST records. Nothing is
reconstructed in Python. Above the data, an evaluation's MF=3 is extrapolation: EMAX is the evaluator's statement.

## 5. The library: which charged-particle evaluations carry LRF=7

Looked at (section 1): the ENDF/B-VIII.1 release metadata (15 sublibraries; the charged-particle archives are 0.2 MB
each, protons 12.3 MB), the material lists (protons 49 materials, deuterons 6, tritons 5, helions 4, alphas 5), the
READMEs and change logs of those five (no word of resonance formats; the helion log reports "GNDS to ENDF
translation issues" fixed in VIII.1), and one file: helions on 4He (MAT 228, "Thompson and INDEN-Light-elements",
EVAL-OCT22): **LRF=7, KRM=4, IFG=1, KRL=0, NAPS=1, LRP=2**, 2 pairs (3He + 4He, p + 6Li), 10 spin groups, 15 levels
(9 of them background poles at ER = 35.07 MeV lab), B_c = -l, MF=3 MT=2 (placeholder 1.0) and MT=600 (pointwise,
7.05-20 MeV), MF=6 MT=2 (LAW=5, LTP=1) and MT=600. "Translated via GNDS to ENDF6 by FUDGE" (its MF=1 text).
p + 12C (VIII.1) carries MF=3/MF=6 only, LTP=12/14, from 1 MeV (`azure2-exchange/references/endf-comparison.md`;
not rerun here). **Not checked**: the other 67 materials. Candidates by author (LANL EDA evaluations of light
systems; the LLNL helion and triton files) need a look at their MF=2 records, which stage 1's reader prints in one
line per file (`pyazr endf summary`).

## 6. The spike (measured)

### 6.1 Round trips of T1 and T2

`roundtrip.py`: a session in the formal basis (`use_brune=False`) gives each level's formal E and gamma (the
parameter records' `start`), B_c = S_c(E_1) from the session's own channel functions; the nucleus is written as
MF=2 text in each number form (`lrf7.py`, with the conversions of section 3.1), read back by the same script,
converted back, and given to the engine as a `#parametrization 0` parameter file (`Session(..., use_brune=False,
parameter_file=...)`). Photon channels were written as P = 1 channels at E_lambda and converted back with the same
factor, so the round trip measures the numbers, not the photon-channel loss. Separation energy of the entrance pair:
from the project, or from the file's masses (the photon pair's MB is the only compound-nucleus mass in the file).

| project | chi^2 (formal; Brune) | number form | S from | chi^2 back | Delta chi^2 | max rel. model | max rel. formal parameter | max dE (eV) | max rel. physical (Brune) | max dB |
|---|---|---|---|---|---|---|---|---|---|---|
| T1 | 482.2752396172; 482.2752396173 | classic (7 digits) | project | 482.2442689 | -3.1e-2 | 4.1e-4 | 4.3e-7 | 3.5 | 1.7e-6 | 4.0e-7 |
| | | E form (6) | project | 482.4453717 | +1.7e-1 | 8.1e-3 | 4.9e-6 | 42 | 1.9e-5 | 3.6e-6 |
| | | widest (8-10) | project | 482.2752275 | -1.2e-5 | 2.6e-7 | 6.9e-9 | 0.0044 | 2.0e-8 | 3.2e-9 |
| | | classic | file masses | 1191.948 | +709.7 | 0.78 | 1.1e-3 | 1126 | 4.3e-3 | — |
| | | E form | file masses | 1237.225 | +754.9 | 0.49 | 1.4e-3 | 1408 | 5.3e-3 | — |
| | | widest | file masses | 482.5043 | +0.23 | 4.7e-3 | 2.0e-6 | 2.0 | 8.3e-6 | — |
| T2 | 10905.2313021; 10905.2313021 | classic | project | 10904.5446 | -0.69 | 5.5e-5 | 3.4e-7 | 0.46 | 1.4e-6 | 4.4e-7 |
| | | E form | project | 10909.6214 | +4.4 | 2.5e-4 | 2.2e-6 | 3.0 | 9.9e-6 | 4.8e-6 |
| | | widest | project | 10905.2317463 | +4.4e-4 | 3.8e-8 | 6.0e-10 | 0.00042 | 1.4e-9 | 4.5e-9 |

(dE is in the formal E_x scale; "physical" are `Session.transform` values; dB is the BND field's rounding: the engine
recomputes its own B, so it does not enter these chi^2; an outside code would use the rounded value, not measured.)
The formal and Brune chi^2 agree to 6e-12 relative. Each T1 session prints the engine's note that the 1/2+ level at
2.365 MeV is transformed by a bracketed root, the original as well as the read-back ones.

Conclusions: the writer uses the widest form (decision 2); the classic form costs about one unit of chi^2 per 10^4 on
these projects; a reader never takes a separation energy from the file's masses.

### 6.2 An ENDF/B-VIII.1 evaluation read into the engine

`inden.py` (run with `python -I` on the fetched file): MF=2 of 3He + 4He read by the same reader; ER and GAM converted
with the lab factor of the 3He + 4He pair; the file's B_c = -l changed to the engine's convention by the exact change
of boundary condition (N = diag(E) - sum_c gamma_c gamma_c^T (B'_c - B_c), diagonalised; B'_c = S_c at the lowest
observed level, found as the lowest root of lambda_k(N(S(e))) = e by a scan and bisection with the session's own
channel functions); a format-3 project written with the formal levels read without transformation
(`transform=False`); an extrapolation `4He(h,p)6Li`, quantity `sigma`, at the 95 energies of MF=3/MT=600.
Separation energies: AME2020 for 3He + 4He (1.587137 MeV), the file's Q (-4.0198 MeV; AME gives -4.019717) for
p + 6Li; masses from the file. 11.5 s.

- **MF=3/MT=600 reproduced**: above 7.2 MeV lab, the engine's sigma(E) agrees with the file to 1.13e-6 at most and
  7.4e-8 in the median (the file prints 7-9 digits). At 7.15 MeV: 2.8e-6; at 7.10 MeV: -1.3e-3; at 7.05 MeV (0.7
  keV c.m. above the p + 6Li threshold) the file has 4.0e-9 b, its own parameters give 3.4e-46 b (2 pi eta = 105).
  The first two points of the file do not follow from its MF=2; the cause is not established (a threshold from
  other masses would move it by tens of keV, which the file's other points exclude).
- **The frame**: the same run with GAM read as a c.m. amplitude is off by up to 143 % (median 54 %), and the 7/2-
  level moves from E_x = 4.5625 to 4.98 MeV. ER and GAM of an LRF=7 file are on the lab scale of the incident pair
  (pp. 66, 355-356), confirmed by the file's own MF=3.
- **Observed levels** from the file (the roots above, in E_x): 3/2- 0.0049, 1/2- 0.4523, 7/2- 4.5625, 5/2- 6.6204 and
  7.1804, 3/2+ 9.3788, then broad states and background; the engine's `Session.transform` of the built project gives
  the same energies.
- **Findings a checker should report**: the 6Li mass in the file (5.963427 m_n = 6.0150996 u) is 2.33e-5 u
  (21.7 keV) below AME2020, while the file's Q agrees with AME to 83 eV; the engine says so on opening the project
  ("separation energies and masses disagree ... +22.6 keV"). The MF=1 text says "Elastic Coulomb scattering is given
  with LTP=2 Legendre terms", the MF=6 record holds LTP=1 (NW = 4 NL + 3 = 27). The file uses LRP=2 (p. 42 reserves
  it for derived files). MF=2 has no gamma pair with KRM=4 (p. 71's "two particle-pairs will always be present:
  gamma + compound nucleus, and neutron + target" is written for neutrons).
- **What the engine needed** beyond today's API: (1) formal levels with a caller's B_c turned into the engine's basis
  (done in the spike by 30 lines of numpy and a root scan: belongs in the engine, door 2); (2) the frame and unit
  conversions (format work, the module); (3) the separation energies (AME through `pyazr.sources.masses`); (4) a
  project made from pairs and levels (`Project` exists; the spike wrote the text). Elastic (MF=6 LAW=5 LTP=1) was not
  reconstructed (not run).

### 6.3 Not run

The LRF=7 files written by the spike were not given to any checking code or to RECENT (none installed); the effect
of BND rounding on an outside reconstruction; MF=3/MF=6/MF=32 writing; the photon-channel loss of option (a) of
decision 1 (KRM=4 photon channels with P = 1) as a number; elastic reconstruction from LAW=5; KRL=1; a second
library file.

### 6.4 The scripts (`work/endf/`, outside the repository)

- `scripts/env.sh`: the environment of the runs (base build `work/b-trial71`, one thread).
- `scripts/lrf7.py`: number forms, ENDF records, the writer and reader of an LRF=7 subsection, the nucleus of a
  formal-basis session.
- `scripts/roundtrip.py <project.azr> <outdir>`: section 6.1 (`run/out_T1/`, `run/out_T2/`: the MF=2 texts, the
  parameter files, `result.json`).
- `scripts/inden.py <mf2.txt> <mf3_mt600.txt> <workdir>`: section 6.2 (`run/inden_work/`: the project, `compare.txt`;
  `GAMMA_FRAME=cm` for the c.m. reading, `run/inden_work_cm/`).
- `scripts/pg.sh <first> [last]`: prints printed pages of the manual's text.

## 7. Reading

### 7.1 An evaluation's resonance parameters as a start model

`pyazr.endf.read(path)` → a `Material` (sections, header); `material.resonances()` → an `LRF7` record (pairs, spin
groups, channels with BND, APE, APT, levels with ER and GAM, in the file's frame and units; refusals of section 3.2
raised here with the three parts); `pyazr.endf.to_project(lrf7, path, separation_energies="ame2020" | {pair: MeV},
basis="brune")` → a new `Project` file, never the open one. It converts ER and GAM (lab factor of the entrance pair,
eV to MeV), asks door 2 for the observed levels (Brune energies, widths or ANCs with signs) and writes them with
`Project.add_pair` / `add_level` / `set_channel`, every level fixed. What the file does not hold and is asked of the
user (the page shows the defaults): the separation energy of the compound nucleus (default AME2020; the difference
from the file's masses is shown), which isotope and energy range (NIS, NER > 1), the J^pi groups and levels to keep
(background poles marked by energy above EH), the data sets (none come from the file), the options of the
calculation (none in the file but KRL), the photon pairs' multipolarities (from the parities and L), and the
reaction labels of the pairs (from MT: 2 elastic, 600+k p_k, 800+k alpha_k, 102 capture). A report lists every
conversion and every refusal.

### 7.2 An evaluation's cross sections as curves: "checking the evaluations"

The plain use: draw an ENDF/B (or JENDL, TENDL, ...) evaluation over the fit and the data in the Results tab, and in
a script.

```python
m = pyazr.endf.read("p-006_C_012.endf")
c = m.curves(s.evaluate(x).data[i])          # the evaluation at the data set's own points and observable
c.energies, c.values, c.coverage, c.notes    # c.m. MeV as the data set, the data set's unit; where the file has data
```

- **Which section**: by the data set's reaction and quantity (`Evaluation.data[i].quantity`, the pairs of its
  reaction): `sigma` of an exit pair → MF=3/MT (2, 600+k, 800+k, 102); `dsigma_domega` of a particle exit → MF=3 x
  MF=6 LAW=2 (sum (2l+1)/2 a_l P_l, eq. 6.8, p. 143, times sigma/2pi, c.m.); elastic `dsigma_domega` and
  `ratio_to_rutherford` → MF=6/MT=2 LAW=5 by its LTP (eqs. 6.9-6.20) with sigma_NI from MF=3/MT=2 for LTP > 2;
  `sfactor` → sigma times the engine's `sfactor_conversion` row; capture to one final state → MF=6/MT=102 yields or
  MF=12 multiplicities times MF=3/MT=102. A quantity the file cannot give (A_y, phase shifts, THM) is said so, not
  drawn.
- **Interpolation**: the TAB1's own laws, INT 1-6 (INT=6 with the reaction's threshold), and for MF=6 the TAB2 law
  between incident energies; values only inside the tabulated range; `coverage` gives it (p + 12C: nothing below
  1 MeV lab).
- **Frames**: MF=3 and MF=6 incident energies are lab energies of the projectile (pp. 103, 134); the data set's
  energies are converted with the engine's `cm_to_lab_energy` row on the pair's masses (and its relativistic flag);
  MF=6 LAW=2 and LAW=5 angles are c.m. cosines (pp. 143-144); a lab data set's angles go through the engine's
  `elastic_lab_to_cm_angle` / `reaction_lab_to_cm_angle` rows and the Jacobian.
- **The Coulomb parts** (eta, k, the Rutherford and Mott cross sections) are the engine's physics rows
  (`sommerfeld_eta`, `wave_number`, `rutherford_cross_section`, `mott_cross_section`); the module evaluates only the
  format's own sums. Units: b, b/sr, MeV.
- **LRP=1 and MF=3 backgrounds**: when MF=2 holds parameters to be added (LRP=1) the curve is MF=3 plus the
  reconstruction: the reader builds the project of section 7.1 in a temporary directory, the engine computes the
  reaction on the data's grid (an extrapolation record, as the spike did), and the module adds MF=3. With LRP=2 (the
  3He + 4He file) MF=3 alone is the evaluation and the reconstruction is a check (below).
- **The check of an evaluation** (`pyazr.endf.check_resonances(material)`): reconstruct from MF=2 with the engine and
  compare with the file's MF=3/MF=6 point by point; report the relative difference, the points outside a tolerance,
  the threshold behaviour, the masses against AME and Q, the header text against the records (section 6.2 shows all
  four on one library file). This is the check the utility codes do not do for LRF=7 (section 9).

## 8. Where the code lives

### 8.1 Rules

One implementation of each formula; no physics in Python (`pyazr3.md` 2.1); an extension is one named place. The
format is format work: records, number forms, interpolation laws, the LAW=2 and LAW=5 sums of coefficients the file
states, the lab factor m_B/(m_A + m_B) of a pair, eV ↔ MeV, fm ↔ 1e-12 cm. That is the module's. Everything that
needs the R-matrix, a Coulomb function or a basis conversion is the engine's, reached through a door.

### 8.2 Engine doors

| door | what | exists | new |
|---|---|---|---|
| formal levels of a session | per J^pi group: channels (pair, l, s, B_c, radius), levels (formal E, gamma) at a fit vector x, in any basis of the session | `formal_from_brune`, `boundary_values`, `brune_from_park` in `physics/`; a formal-basis session's `Parameters.x0` | `Session.formal_levels(x)`: `AZUREAPI::FormalLevels`, one binding (about 120 lines C++, a test against a `use_brune=False` session); refuses as the standard basis refuses (no formal equivalent) |
| formal levels with a caller's B_c → observed | the reader's conversion: a System (pairs), a block's channels, E, gamma, B → Brune E, widths or ANCs (and the status per level) | `brune_from_formal(system, b, E, gamma, B, cf)` and `observed_from_amplitudes` in `physics/` take B as an argument | a free function of the binding (`pyazr.physics`-like, on pair records; about 150 lines C++ with the System built from the pairs, a test: formal parameters of T1 with the engine's B give T1's levels; with B = -l, the 3He + 4He levels of section 6.2) |
| cross sections on a given grid | the reaction and quantity at energies (and angles) | an extrapolation record and `Session(extrapolations=True)`; `pyazr3.md`'s planned `Session.evaluate` on given energies | none for stage 1-3 |
| nuclear amplitude coefficients a_l, b_l (LAW=5 LTP=1) | from the collision matrix | no | stage 5, a session function row |
| covariance of formal parameters | the fit's covariance through the Jacobian of the transformation at the written B_c (section 13.13; HESSE in a formal-basis session until 2026-11-18) | door 2 and iminuit's HESSE | none |

### 8.3 The module `pyazr.endf`

| file | holds | about |
|---|---|---|
| `records.py` | TEXT, CONT, HEAD, LIST, TAB1, TAB2, INTG read and write; the three number forms; line numbering, SEND/FEND/MEND/TEND | 250 lines |
| `interpolation.py` | INT 1-6 (Table 16, eqs. 1-6, pp. 25-26) | 80 |
| `material.py` | `Material`: sections by (MF, MT), MF=1/MT=451 header and directory (built, never typed), `summary()` | 220 |
| `resonances.py` | `LRF7`, `SpinGroup`, `Channel`, `Level`, `Pair` (plain records in ENDF units and frame); read, write; the refusals of section 3.2 | 280 |
| `distributions.py` | MF=3 curves, MF=6 LAW=2 and LAW=5 (LTP 1, 2, 12, 14, 15) evaluated at given (E, mu) | 300 |
| `model.py` | `from_session(session, x)` → `LRF7` (door 1, frame, units, photon channels by decision 1); `to_project(...)` (door 2); the losses list | 350 |
| `curves.py` | `Material.curves(data)`, `check_resonances(material)` | 200 |
| `checkers.py` | find, decks, run, parse reports, candidate files and their differences (section 9) | 300 |
| `writer.py` | `write(path, session, x, header, include=...)`: MF=1, 2, 3, 6, (12), 32 from the engine | 250 |

The calls:

```python
pyazr.endf.read(path) -> Material                         # any ENDF-6 file; sections parsed on demand
Material.summary() -> dict                               # MAT, NSUB, LRP, sections, LRF/KRM/IFG/KRL, energy ranges
Material.resonances(isotope=0, range=0) -> LRF7          # refusals raise pyazr.endf.FormatRefused (what, allowed, fix)
Material.curves(data_or_spec) -> Curve                   # section 7.2
pyazr.endf.check_resonances(material, tolerance=1e-4) -> Report
pyazr.endf.to_project(lrf7, path, separation_energies="ame2020", basis="brune") -> (Project, Report)
pyazr.endf.from_session(session, x) -> (LRF7, Losses)    # the nucleus as LRF=7, and what it cannot carry
pyazr.endf.write(path, session, x, header, number_form="widest", include=("mf2", "mf3", "mf6"),
                 grid=None, covariance=None) -> (path, Report)
pyazr.endf.checkers(search=None) -> list[Checker]         # found programs: name, path, version, kind
pyazr.endf.check(path, programs=("checkr", "fizcon", "psyche")) -> list[CheckReport]
pyazr.endf.standardize(path, program="stanef") -> Candidate  # a new file beside the old; .diff(), .digits_lost
Candidate.accept(path=None)                                 # writes; never over the original unless asked
```

Errors are pyazr's (`FormatError` for a malformed record with file:line:column, a new `FormatRefused`
(ValueError) for a valid file AZURE2 cannot read, `EngineError` from the doors), each with what, allowed, fix.

### 8.4 Ferdinand, FUDGE, GNDS

Ferdinand reads `.azr` files of the earlier formats and writes GNDS and ENDF; FUDGE reads and writes GNDS and
translates GNDS ↔ ENDF-6 (the 3He + 4He file was made that way). As the route of AZURE2 they would be a dependency
(none is installed here; not installable by this work) and they do not read format 3 (`azure2-exchange`
SKILL.md section 3). Ferdinand's known problems are the ones this design measures and avoids: the separation energy
from the masses (its problem 1; section 6.1: +710 in chi^2 from a 1.1 keV error), photon channels treated as particle
channels (problem 3; decision 1 keeps them out of MF=2). The module does not make their use harder: its files are
plain ENDF-6 that FUDGE translates to GNDS (`endf2gnds`), and its `LRF7` record holds exactly GNDS's
`resonances/resolved/RMatrix` content (spin groups, channels with boundary condition and radii, a table of energies
and reduced-width amplitudes, `projectileFrame="lab"`), so a GNDS writer is a later file of the same module if the PI
wants it (decision 6). An evaluator who has FUDGE writes ENDF with the module and converts; one who has Ferdinand
keeps using it on a format-2 copy.

## 9. The checking codes

### 9.1 What exists, how it is started, what it writes

| code | what (manual) | started | writes | rewrites the evaluation |
|---|---|---|---|---|
| CHECKR 8.24 (NNDC, IAEA) | format: "any file that passes ... without error messages fully conforms to the ENDF format" (manual p. 6) | batch deck on standard input: input file, message file, MAT range (blank: all), `DONE`; or `CHECKR in%out%Y` (pp. 6-7) | a message file: header, "ERROR(S) FOUND IN MAT=..., MF=..., MT=..." blocks or "No problems to report" | no |
| FIZCON 8.20 | procedures and simple physics: energy order, widths adding up, Q values, missing sections, normalisation; optional SUMUP and deviant-point tests with a fractional tolerance (pp. 8-9) | deck: input, messages, `MAT1,MAT2,deviant,sumup,tolerance` | message file, ends "Encountered N errors, M warnings" | no |
| PSYCHE 8.11 | physics: energy conservation, Wick's limit, resonance statistics, Q against mass tables (p. 12) | deck as CHECKR | message file of tables and warnings | no |
| INTER | thermal values, resonance integrals, 14 MeV (pp. 14-16); needs pointwise linearised MF=3 | deck | results file | no |
| STANEF 8.10 | bookkeeping: tape label, MT=451 directory and HSUB lines, resequencing, "numeric fields to standard form", binary (pp. 4-5) | deck: input, output, label, `Y` standard options (or the four option fields) | a new ENDF file | **yes** |
| PREPRO LINEAR, RECENT, SIGMA1, ... (Cullen, 2023) | linearise; reconstruct MF=3 from MF=2 (RECENT: LRF=7, Coulomb penetrabilities since 2004, "RML charged particles" since 2020, SHIFT=1 since 2021, from its source history) | fixed names in the working directory: `RECENT.INP` (cards), `ENDFB.IN` → `ENDFB.OUT`, listing `RECENT.LST` | derived files (LRP=2), not corrections | derived file |
| PREPRO FIXUP, DICTIN | FIXUP: redundant cross sections redefined as sums, ZA/AWR made consistent, sequencing (its header); DICTIN: the MF=1 directory | `FIXUP.INP`, `ENDFB.IN` → `ENDFB.OUT` | a new ENDF file | **yes** |
| ENDF2C (PREPRO) | numbers rewritten in a form C and Fortran read alike | `ENDFB.IN` → `ENDFB.OUT` | a new file | numbers' form |

What they check in LRF=7, read in their sources (not run): CHECKR (`CHECK_RL`) tests IFG in 0-1, KRM in 1-4, KRL in
0-1, NPP >= 1 and the two list lengths, then reads two LIST records per spin group (a file with KBK or KPS records
would be misread); FIZCON (`CHKRL`) reads the records, its test of the pairs' charges is commented out; PSYCHE reads
and skips them. So for an AZURE2 file the codes check the material's bookkeeping, MF=1, MF=3, MF=6 (FIZCON's
normalisation of distributions, energy order, Q values) and the structure of MF=2, not the R-matrix.

STANEF's standard form is the classic one: its `NORMAL` routine writes `F9.6` mantissa and a one-digit exponent
(7 significant digits). It is applied to every numeric field by default, MF=2 lists only on request (manual pp. 4-5).
On a widest-form file it rounds MF=3/MF=6 to 7 digits, and with the MF=2 option it costs what section 6.1's classic
row shows.

### 9.2 Finding them

`pyazr.endf.checkers(search=None)`: an explicit list of paths, else `$AZURE2_ENDF_UTILITIES` and `$AZURE2_PREPRO`
(directories), else `PATH`; names `checkr`, `fizcon`, `psyche`, `inter`, `stanef`, `linear`, `recent`, `fixup`,
`dictin`, `endf2c` (case and `.exe` variants). The version is read from the banner a run prints (CHECKR: "PROGRAM
CHECKR VERSION 8.24"); no program is run to find it. A machine setting, never a project key.

### 9.3 Running and reading the report

Each program runs in a fresh temporary directory (the PREPRO codes read and write fixed names), with a copy of the
file, a deck written by the module on standard input, one process, a timeout, its exit status, standard output and
the message file kept. The report is parsed into messages: program, MAT, MF, MT, severity (ERROR, WARNING, the
"SHOULD BE" lines of PSYCHE as notes), text, and the line range of the section, so the page links a message to the
section's text and, for files the module wrote, to the project item that made it (a level, a pair, a grid).

### 9.4 What "correct it" can mean

- **Reported, not changed**: everything CHECKR, FIZCON and PSYCHE say. The module shows the messages; the user
  fixes the project (a level, a mass, a grid) and writes again. The module's own writer makes the bookkeeping right
  by construction (directory counts, sequence numbers, MOD, HSUB), so those messages should not appear on its files
  (a test).
- **Rewritten by a program, shown before it is kept**: STANEF, DICTIN, FIXUP (and ENDF2C) produce a candidate file.
  `Candidate.diff()` lists the sections that changed and, per section, whether only bookkeeping changed (directory,
  sequence, text lines) or numbers; for numbers, how many and the largest relative change, and the digits lost
  (STANEF on a widest-form file). Nothing replaces a file until the user accepts, and the default target is a new
  file beside the original.
- **Never**: a silent rewrite of numbers, a change of a parameter to make a check pass, or a run of a program the user
  did not choose.

### 9.5 The stand-in program for the tests

`tests/endf/standin_checker.py`: a small Python program that reads the deck of any of the codes from standard input
(or `in%out%Y`), and writes a message file in the program's layout from a table given by the test (no message;
"ERROR(S) FOUND IN MAT=9999, MF=3, MT=600" with a line; a crash with a non-zero exit; a candidate file with the
directory changed, or with numbers rounded to 7 digits as STANEF does). The tests point `AZURE2_ENDF_UTILITIES` at a
directory of links named `checkr`, `stanef`, ... to it, and check finding, the deck, parsing, the candidate's diff and
that no file is replaced without `accept`. The real codes' sample reports (`util_src/*/*.out-0`) are the models of
the layouts (the test keeps a few lines of each, written by hand).

## 10. The page and the program

### 10.1 Screens

- **Results tab: Evaluations.** "Add an evaluation" (a file in the project's root, or one fetched by URL into the
  root): its summary (MAT, library, author, date, sections, energy range). For each data set, the curve is drawn
  over the fit and the data when the file has the observable; a legend entry per file; a note where the file has
  nothing (below its first energy, an observable it lacks). A table: chi^2 of the data against the evaluation per
  data set, beside the fit's. "Check the resonance parameters" (when MF=2 is there): the reconstruction against MF=3
  as a ratio plot and the list of findings (section 6.2's kind).
- **Project tab: Export ENDF.** The header fields asked (MAT, NLIB, LREL, NVER, lab, authors, reference, dates, EMAX;
  defaults from the project and the pairs: ZA, AWR, AWI, NSUB), what to include (MF=2, MF=3, MF=6, MF=12, MF=32 when
  the fit has a covariance), the energy range and grid tolerance, the number form (widest by default). Before
  writing: the list of losses (decision 1, external capture, the apparatus not carried, ...) and refusals, the MF=1
  text to edit. Writing is a job; the result offers "Check".
- **Project tab: New project from ENDF.** A file; its summary; the isotope and range; the separation energies (AME
  by default, the file's masses shown with their difference); the levels kept; the new project's name. The result
  opens as a new project; the conversion report stays beside it.
- **Check.** The programs found (path, version, kind), the ones to run, the file; each run a job; the report as a
  table of messages grouped by MF/MT, each linked to the section's text; for STANEF, FIXUP, DICTIN the candidate's
  difference and "Keep as new file" / "Replace" (the second asks again).

### 10.2 Routes and jobs

Rows of the route table (`web_ui.md` 4.2), each one pyazr call: `POST /api/endf/summary` → `read(path).summary()`;
`POST /api/endf/curves` → `Material.curves`; `POST /api/endf/check_resonances` → `check_resonances`;
`POST /api/endf/import` → `to_project`; `GET /api/endf/checkers` → `checkers()`. Jobs (`web_ui.md` 4.4), new kinds:
`endf_write` → `write`, `endf_check` → `check` (one job per program), `endf_standardize` → `standardize`, with
`accept` → `Candidate.accept`. The jobs use one thread each.

### 10.3 The program

No `AZURE2 endf` command: the format work is Python, the engine has nothing to add beyond the doors, and a second
implementation in C++ would break the rule of one. The command line is pyazr's (`pyazr3.md` 2.8 table):
`pyazr endf summary FILE`, `pyazr endf curves FILE PROJECT`, `pyazr endf import FILE NEW.azr`,
`pyazr endf write PROJECT OUT.endf [--header header.json]`, `pyazr endf check FILE [--program checkr ...]`, each with
`--json` (decision 5).

## 11. Decisions for the PI

Status (2026-11-14): decisions 1 (photon channels out of MF=2, capture pointwise in MF=3/MF=6), 3 (LRP=2 with full
pointwise MF=3), 7 (the MAT 228 extract committed as test data) and 10 (MF=32 from a formal-basis HESSE) are
**decided by the PI**, as recommended. Decisions 2, 4, 5, 6, 8 and 9 are the director's interim choices, as
recommended (HANDOFF, 2026-11-14 15:00). As built (sections 13.4, 13.5, 13.8): 1 (no photon channel in MF=2;
capture in MF=3/MT=102 and the primary photons in MF=6/MT=102), 2 (widest form by default), 3 (LRP=2, one named
value `writer.LRP`; MF=3 and MF=6 from the engine; MF=2 checked against MF=3 before writing), 4 (KRL=1 refused on
reading and writing), 7, 8, 9 (AME2020 for the incident
pair, the file's Q values for the others) are in the code.

1. **Photon channels.** (a) In MF=2 as channels of photon pairs (KRM=4, P = 1, Gamma_gamma kept at E_lambda): the
   radiative widths then enter the level matrix of every code that reads the file, the energy dependence
   E_gamma^(2L+1) and external capture are lost, and the file's capture is not AZURE2's. (b) MF=2 holds the particle
   channels only (exact); capture is written pointwise in MF=3/MT=102 and MF=6/MT=102 from the engine (exact on the
   grid, external capture included), and the photon widths are listed in the MF=1 text. **Recommended: (b)**; it
   keeps the PI's rule (radiative widths out of the total width) in what a processing code computes. A project with
   `use_rmc=1` is written as KRM=3 either way.
2. **Number form**: **widest** (fixed point where it fits, else the exponent without E; 8-10 digits) by default, the
   classic form on request; never the E form. Reason: section 6.1 (Delta chi^2 1e-5 against 3e-2 to 0.7).
3. **LRP and what MF=3 holds.** (a) LRP=1 with zero MF=3 backgrounds for the particle reactions: the file is minimal
   but every user depends on a processing code reconstructing charged-particle LRF=7. (b) LRP=2 with full pointwise
   MF=3/MF=6 from the engine and MF=2 for reuse, as ENDF/B-VIII.1 does for 3He + 4He, though p. 42 reserves LRP=2 for
   derived files. (c) LRP=1 with MF=3 pointwise and MF=2 particle channels and backgrounds zero only where MF=2
   contributes. **Recommended: (b)**, with a sentence in MF=1 that MF=2 reproduces MF=3 (the module checks it before
   writing); it is what the library ships and what a transport user reads without reconstruction.
4. **Relativistic projects (KRL=1)**: refuse until the convention is fixed (the manual has no equations; EDA's
   E(s) convention and AZURE2's linear E_x agree only at the levels), or write KRL=1 with `relativistic.out`'s
   mapping and a stated loss. **Recommended: refuse**, with the fix "write with relativistic=0 or ask for the KRL=1
   convention to be decided".
5. **A command of the program**: **none; `pyazr endf ...`** (section 10.3).
6. **GNDS**: not now; the `LRF7` record is GNDS's R-matrix content, so a writer is one file later. **Recommended: not
   now.**
7. **Test data from the library.** Commit an extract of ENDF/B-VIII.1 MAT 228 (MF=1 header, MF=2 and MF=3/MT=600:
   about 140 lines, 11 kB) as the reader's and the reconstruction's reference (section 6.2: 1.1e-6)? Source:
   NNDC, ENDF/B-VIII.1 helion sublibrary, archive MD5 402f3675ca8fb2b4dc02624a1de7ca3a. Terms: the archive's LICENSE
   asks only that the release be cited (Brown et al., Nucl. Data Sheets 148, 1 (2018), as the file states). Not
   committed by this work. **Recommended: yes, the extract with its source and citation in a README.**
8. **The parser**: the module's own (records are a few hundred lines; no dependency) or ENDFtk / the `endf` package
   (both not installable here; the second does not interpret charged-particle sublibraries, `endf-comparison.md`).
   **Recommended: the module's own.**
9. **Separation energies on import**: AME2020 by default with the file's Q values, the file's masses for the
   kinematics, and the differences reported. **Recommended: so** (section 6.1: the file's masses alone gave 1.1 keV).
10. **Covariances (MF=32)**: from HESSE in a formal-basis session at the transformed best fit (no new Jacobian), or
    the fit's covariance propagated through the Brune → formal Jacobian (a new door). **Recommended: the first**,
    stage 5. *The owner's answer (2026-11-18): the second, the fit's covariance carried through the transformation;
    built (section 13.13).*

## 12. The order of work

Each stage merges alone, with its tests and its user-guide page section (`docs/source/user_guide/endf.rst`, new), and
changes no output number of the engine.

| stage | what | needs decided | files, lines (estimate) | tests |
|---|---|---|---|---|
| 1. Curves | `records.py`, `interpolation.py`, `material.py`, `distributions.py`, `curves.py`: read any ENDF-6 file, summary, MF=3 and MF=6 (LAW=2, LAW=5 all LTP) at a data set's points; Results-tab overlay; `pyazr endf summary/curves` | 8 | 5 Python files, 1 100; page 250; route rows | a hand-written file per law (INT 1-6, LAW=2, LAW=5 LTP 1/2/12/14/15) against hand values; the classic, E and widest forms read back; the p + 12C data of a test project against an evaluation curve made by hand |
| 2. LRF=7 reader and the start model | `resonances.py`, `model.py` (to_project), door 2, refusals of section 3.2, "New project from ENDF" | 7, 9 | 2 Python files 450; C++ 150 + binding; page 150 | write-then-read of T1 and T2 (widest form, model to 1e-6, Delta chi^2 < 1e-3: section 6.1 gives 3e-7 and 4.4e-4); a hand-made LRF=7 with B_c != S(E_1) against the engine; the MAT 228 extract (decision 7) to 1e-5; each refusal with its three parts |
| 3. The writer | `from_session` (door 1), `writer.py`: MF=1, MF=2, MF=3, MF=6 from the engine; the losses list; Export screen; `pyazr endf write` | 1, 2, 3, 4 | 2 Python files 450; C++ 120 + binding; page 200 | T1, T2, 13N written and read back; MF=3 of the written file against the engine at its points; the directory and sequence by construction; a project with `use_rmc=1` as KRM=3 |
| 4. Checkers | `checkers.py`, stand-in program, Check screen, jobs, candidates and their difference | — | 1 Python file 300; stand-in 80; page 200 | the stand-in in each role (clean, errors, crash, a STANEF candidate with digits lost); no file replaced without accept |
| 5. Covariances and LTP=1 | MF=32 LCOMP=1 from formal HESSE; MF=6 LAW=5 LTP=1 (door: nuclear amplitude coefficients) | 10 | Python 200, C++ 150 | covariance read back equals the written one; LTP=1 against the engine's dsigma/dOmega |
| 6. The check of an evaluation | `check_resonances` (reconstruction against MF=3/MF=6, masses against AME and Q, header text against records); RECENT as a second reconstruction when found | — | Python 150; page 100 | the MAT 228 extract gives the three findings of section 6.2 |

Stage 1 alone answers "an easy way of checking the evaluations" for every library file with MF=3/MF=6 (all 69
charged-particle materials of VIII.1); stage 2 adds the evaluations that carry LRF=7. The class scheme
(`docs/source/developer/figures/classes.tsv`) gets the module's records as rows when they are built, not before.

## 13. As built

### 13.1 Stage 1: read, summarise, curves

`pyazr/endf/` holds `records.py` (lines, the three number forms, CONT, LIST, TAB1, TAB2 and a cursor over a
section), `interpolation.py` (INT 1-6), `names.py` (the manual's names of files, reactions, sublibraries, formats),
`material.py` (`read`, `Material`, `Header`), `distributions.py` (MF=3, MF=6 LAW=2 and LAW=5), `summary.py`,
`curves.py` (`CurveRequest`, `Curve`, `CoulombTerms`, `request_of`) and `format_check.py` (the own check of
stage 4, built here because the reader and the check share the cursor). `pyazr.FormatRefused` (`ValueError`) is a
class of `pyazr/errors.py`; a malformed file raises pyazr's `FormatError` with one diagnostic (file, line, columns,
found, allowed, the page). Commands: `pyazr endf summary FILE`, `pyazr endf curves FILE PROJECT [--record]`.

What differs from the design, and why:
- `Material.curves(data, session=s)`: an `EvaluatedData` does not carry its pairs, so the session is passed (its
  `datasets` and `pairs` give them); `CurveRequest` is the plain form (`request_of(session, data)` makes it).
- The lab energy of a point is the engine's own (`EvaluatedData.extra["energy_lab"]`, after the energy shift);
  the `cm_to_lab_energy` row is used only when the file's projectile is the pair's second particle. The data set's
  angles are already c.m. (the engine's), so no lab-angle row is needed.
- The MT of a particle exit is found by the file's own Q values (QM - QI = the residual's excitation, p. 105).
- `check_resonances` (stage 6) and the LRP=1 sum of MF=2 and MF=3 are not built: an LRP=1 file gives MF=3 alone.
- Capture to one final state (MF=6/MT=102 yields, MF=12) is refused; capture to every final state is MF=3/MT=102.

What the manual leaves open, and how it is read (each in the module's text):
- INT=6 threshold T for Q <= 0 (eq. 2, p. 26: "the threshold", no formula): the lab energy at E_cm = -QI, the
  engine's `cm_to_lab_energy` on the file's masses (AWI of MF=1, AWR of the section).
- INT=4, 5 with a zero end: 0 inside the interval (the law's limit); a change of sign: NaN.
- Between the incident energies of MF=6: the value at fixed mu (Cartesian, eq. 7 p. 26); unit-base laws 21-25 taken
  as 1-5 (equal for distributions over the same mu range, p. 27); corresponding-point laws refused. For LAW=5 the
  interpolated quantity is sigma_e - sigma_c for every LTP (the manual tabulates it so for LTP 12-15 and is silent
  for LTP 1 and 2); the Coulomb part is computed at the energy asked.
- LTP=12 and 15 both "linear in mu" (pp. 146-147): read alike; LTP=14 ln P linear in mu.
- LAW=5 has no yield factor (MF=3 holds 1.0 or sigma_NI, p. 145); LAW=2 multiplies by the product's yield.
- A tabulated distribution outside its cosines is zero (stated for LAW=1, p. 142; taken for LAW=2 and LAW=5).
- The TPID line may carry MAT 0 (the library file does); it is the first line with MF = MT = 0 and not an end record.
- Fields after a record's last value: blanks or zeros (FUDGE pads the interpolation line with zeros) are accepted;
  anything else is a warning of the check (a count smaller than what was written, which a reader cannot see).

Tested (`tests/pyazr/endf_read_test.py`, 82 checks): against hand values from the manual's equations -- the number
forms, INT 1-6 (INT=6 in both of the manual's forms), MF=3 on a table with five laws, LAW=2 LANG 0, 12, 14 at and
between panels, LAW=5 LTP 1, 2 (distinguishable and identical), 12, 14, 15, the Coulomb term against eqs. 6.9-6.12
(2e-6), the curves of sigma, dsigma/dOmega, sfactor and ratio_to_rutherford with the lab energy by hand; against
itself only -- the session route (`curves(data, session=s)` equals the request built by hand on T1). The column
rules: a line beyond column 80 (read; a warning), a wrong sequence number (read; a warning), a count larger than
the lines (refused; an error), lines left after the last record (refused; an error), a count one smaller on the
same lines (read; a warning: the reader cannot know), a number with D (refused with file:line:columns and p. 31), a
line without MAT, MF, MT (refused), a wrong directory NC and a missing SEND (errors of the check). Two injected errors
(the Legendre factor (2l+1)/2 made 2l+1; INT=3 made linear) failed nine checks; removed.

The own check on the library file `h-002_He_004.endf` (not in the repository) finds no error: one warning (LRP=2 in
a primary evaluation, p. 42) and one note (no sequence numbers). On the committed extract it finds in addition the
three sections the directory lists and the extract left out (the test expects exactly these).

### 13.2 Stage 4: the checking codes

`pyazr/endf/checkers.py`: `PROGRAMS` (one row per code: family, kind check / rewrite / derive, what it does),
`checkers(search=None)`, `check(path, codes=None, timeout=300)` (the own format check first, then CHECKR, FIZCON,
PSYCHE), `standardize(path, program="stanef", mf2=False)` returning a `Candidate` (`diff()`, `summary()`,
`accept(to=None, overwrite=False)`), `compare(a, b)`. Commands `pyazr endf check FILE [--program]` and
`pyazr endf standardize FILE [--program] [--mf2] [--write NEW [--overwrite]]`. The stand-in of the tests:
`tests/pyazr/endf_standin/standin.py`; the real codes run in `tests/pyazr/endf_check_test.py` when
`AZURE2_ENDF_UTILITIES` and `AZURE2_PREPRO` are set.

The codes were built by the director under `opt/` with the PI's yes (`opt/README-endf-codes.md`: utility codes at
ea2638b2, PREPRO 2023 at cc5fd087) and run here. Where the manuals, the sources and the built programs differ, the
programs' behaviour was taken:
- **The deck ends with DONE** (checkr.f lines 255-270, 806-861). Without it the program reaches the end of its input
  after a complete report and stops with a Fortran runtime error, exit status 2. With it the codes exit 0 whatever
  they found: the exit status says nothing about findings (the director's first runs had read status 2 as
  "findings"; it was the end of input).
- **CHECKR counts "SECTION CANNOT BE CHECKED" as an error**; FIZCON wraps one message over two lines (seen twice);
  PSYCHE prints tables and comments under "FILE n / SECTION mt" and has no count. The parser reads them so and
  compares its count with the report's "Encountered N errors, M warnings" (a finding when they differ).
- **STANEF's standard options rewrite MF=2 too** (the manual's pp. 4-5 say MF=2 on request; the built 8.10 rounds
  the 67 changed numbers of the 3He + 4He MF=2 with blank options as with `0,Y,Y,Y`).
- **STANEF corrupts a float written without a decimal point**: `35070343` and `12377043` (ER of background poles,
  eV) became `3.507034+3` and `1.237704+3`. The manual allows such a field (p. 31: "digits optionally containing a
  decimal symbol"). The cause (the director, read in the sources and reproduced with STANEF): the four utility
  codes read floats with `E11.4` (`stanef.f` line 1131; 7 READ statements in STANEF, 10 in CHECKR, 8 in FIZCON, 9
  in PSYCHE) where the manual says `F11.0`; a field without a decimal point then takes its last four digits as the
  fraction. The whole library file has ten such lines (four values) in MF=2/MT=151. That CHECKR, FIZCON and PSYCHE
  test the misread value is read from their sources, not run.
- **FIXUP 2023-1 refuses to run without option 5** ("Cross Section Reconstruction is MANDATORY"): the deck sets
  column 1 (ZA, AWR) and column 5 = 2 (built-in summation rules), every other test off.
- **CHECKR 8.24 does not know LTP=15** (p. 147 of the 2023 manual defines it).

Measured on the ENDF/B-VIII.1 3He + 4He file (`standardize` on a copy; digits lost = significant digits written
before minus after, the largest in the section):

- STANEF (standard options, and with MF=2): MF=1 bookkeeping; MF=2/151 67 numbers, 2 digits lost, two values 1e4
  too small; MF=3/2 1 number, 2 digits; MF=3/600 30 numbers, 2 digits (4.9e-7); MF=6/2 349, 2 digits (4.7e-7);
  MF=6/600 2, 2 digits.
- FIXUP (columns 1 and 5): MF=1 bookkeeping (51 -> 66 lines); MF=2, MF=6 numbers written another way, values equal;
  MF=3/2 other fields; MF=3/600 4 numbers, 1 digit (2.5e-6).
- DICTIN: MF=1 bookkeeping (51 -> 52 lines); every other section unchanged.
- ENDF2C: MF=1 bookkeeping; MF=2, MF=3 form only; MF=6/2 1441 numbers, 1 digit (3.6e-5); MF=6/600 540, 1 digit
  (2.2e-5).

So the design's reading of the sources holds for STANEF's 7 digits and DICTIN's bookkeeping; ENDF2C loses one digit
in exponent fields (`-1.38025-12` -> `-1.3803E-12`), FIXUP one in a few MF=3 fields; STANEF's handling of MF=2 and of
integer-looking floats is worse than the sources suggested.

The real codes on the files of the repository (exit status 0 everywhere; every finding the reports count was read):

- The 3He + 4He extract: own check 3 errors, 1 warning, 1 note; CHECKR 4 errors (LRP=2, three sections missing);
  FIZCON 1 (MF=3/600 needs MF=6); PSYCHE 1 warning (AP: LRF=7 misread).
- The library file itself (not committed): own 0/1/1; CHECKR 5 errors (LRP=2, NE > 200 twice and "cannot be checked"
  twice); FIZCON 1 (MF=6/600: energy range of the LIST records against the TAB1); PSYCHE 1 warning (AP), 8 notes.
- `p-C012_hand.endf` (the example): own 0/0/0; CHECKR, FIZCON none; PSYCHE 3 notes.
- The read test's hand files (MT=800/801; three number forms): own none (a note when unsequenced); CHECKR none;
  FIZCON 4 (EMAX not reached twice, LIST against TAB1 energy range twice); PSYCHE 2 warnings (a negative LAW=2
  distribution).
- The LAW=5 files (LTP 1, 2, 12, 14, 15): own none; CHECKR none, but 3 for LTP=15 (not known to it); FIZCON EMAX and,
  for LTP 12-15, the normalisation of P_NI; PSYCHE notes. Identical-particle files: CHECKR 1 warning (ZSYMAM
  text), FIZCON AWI atomic and EMAX.

The findings on the hand-written files are true of them (made-up coefficients, a 6 MeV table under EMAX = 20 MeV);
the own check does not test physics or procedures (FIZCON's and PSYCHE's domain): it found none of them, and no code
found a format error the own check missed on these files.

### 13.3 What the page needs (routes and screens): as built

Built as listed below (`web_ui.md` section 11.13: the ENDF tab, `pyazr/ui/evaluations.py`). Differences: the curves
route takes several files at once and answers every data set of the session with the data times its norm, the
model and two ratios formed in the server's Python (`ratio_to_data` with its error from the datum's, and
`ratio_to_model`); the codes' list adds, per program, what it does to numbers (a fourth column of `PROGRAMS`); the
accept is `POST /api/endf/accept` (job id, `to`, `overwrite`), not a route of the job table; the two jobs build no
session. Two rows of `MEANINGS` were added for CHECKR's LTP=15 messages ("INVALID ANGULAR REPRESENTATION LTP=15",
"NW SHOULD BE SET TO ****"), seen with the built CHECKR 8.24 on a LAW=5 LTP=15 file. Not built here: "New project
from ENDF" and the export (built since, `web_ui.md` section 11.16: the jobs `endf_import` and `endf_write`, the
routes `POST /api/endf/resonances`, `POST /api/endf/keep_project`, `GET /api/endf/export`), extrapolations' curves
(only data sets; the curves route also gives each file's values on the grid of a data set's model curve).

The list as designed:

- Results tab, "Evaluations", add a file: `pyazr.endf.read(path).summary().to_dict()` on a path in the project's
  root; answer: MAT, target, projectile, library, lab, date, authors, LRP, EMAX, the sections (MF, MT, reaction,
  lines, what) and the text.
- The curve of a data set: `Material.curves(data, session=s)` (or `pyazr endf curves FILE PROJECT --record NAME
  --json`) for the file and a data set of the open session; answer: energies, angles, values (NaN where none),
  coverage, quantity, unit, source, notes; or a refusal's what, allowed, fix.
- "Check ENDF", the codes found: `pyazr.endf.checkers()`; answer: name, path, family, kind, what it does.
- Run a check (a job, one thread): `pyazr.endf.check(path, codes=[...]).to_dict()`; answer: the own findings (line,
  MF, MT, message, page) and per code its version, exit status, findings (MF, MT, record, text, meaning) and report.
- "Standardize" (a job): `pyazr.endf.standardize(path, program=, mf2=)`; answer: the candidate's path and its
  `diff()` rows; "Keep as new file" is `accept(to=)`, "Replace" `accept(to=path, overwrite=True)` after a second
  question.
- "New project from ENDF" (a job, one thread): `m = pyazr.endf.read(path)`, `m.resonances(isotope=, energy_range=)`
  (a refusal's what, allowed, fix shown as it is), then `pyazr.endf.to_project(record, new_path,
  separation_energy="ame2020" | MeV, material=m)` (or `pyazr endf to-project FILE NEW --json`); answer: the new
  project's path (open it as the project), the separation energy and its source, each nuclide's mass in the file
  against AME2020 (keV), the pairs (key, MT, Q, S, excitation), the levels (spin group, ER, formal and Brune E_x,
  the widths and ANCs written, the iteration's status, background when ER > EH), what was not written, the notes.
  The page asks: the file, the isotope and range when MF=2 has several, the separation energy (AME2020 by
  default), the new project's name.
- "Export ENDF" (a job, one thread): `pyazr.endf.write(session, new_path, x, header={...}, energy_range_MeV=(lo,
  hi), tolerance=, number_form=)` (or `pyazr endf write PROJECT NEW --json`); answer: the file, its sections (MF,
  MT, lines), the grid of each MF=3 section (points, range), the check that MF=2 reproduces MF=3 (largest relative
  difference per MT), the losses, the sections not written. The page asks the MF=1 fields (MAT, NLIB, NVER, LREL,
  NMOD, ALAB, EDATE, AUTH, REF, DDATE, RDATE, the library's name in HSUB, free text), the energy range (default:
  the data's), the grid tolerance (1e-3) and the number form (widest); before writing it shows
  `pyazr.endf.from_session(s, x)[1]` (the losses and the photon widths). The result offers "Check" on the file.

Routes, one per call: `POST /api/endf/summary`, `POST /api/endf/curves`, `GET /api/endf/checkers`; jobs
`endf_check`, `endf_standardize` (with `accept`), `endf_import` (`to_project`) and `endf_write` (`write`).

### 13.4 Stage 2: the start model

`pyazr/endf/resonances.py`: `ParticlePair`, `Channel`, `Level`, `SpinGroup`, `LRF7` (plain records in the file's
units and frame, `print()` lists them, `LRF7.rows(number_form)` writes MF=2/MT=151), `read_resonances` behind
`Material.resonances(isotope=0, energy_range=0)` with every refusal of section 3.2 (three parts). `pyazr/endf/
model.py`: `to_project(record, path, separation_energy="ame2020" | MeV, material=None)` -> `(Project, Conversion)`,
`from_session(session, x)` -> `(LRF7, Losses)`. Command `pyazr endf to-project FILE NEW`.

The doors (`api/src/pybind_level_doors.cpp`, registered in `BindPhysics`): `formal_levels_of_block` (door 1:
a block's levels in the Brune or Park basis, as amplitudes or in the observed form, -> formal E, gamma and the
engine's B_c, through `physics::level_parameters`) and `observed_levels_of_block` (door 2: formal levels with the
caller's B_c -> Brune E, amplitudes, widths and ANCs, through `brune_from_formal` and `observed_levels`; its
iteration tolerance an argument, the engine's 1e-6 MeV by default). Both take plain pair records and build the
`System` and the engine's `PPair` objects themselves, with the engine's `ChannelFunc` / `WhitFunc` and the session
option `use_coulomb_table` (and `use_gsl_coul` until it was retired, 2026-11-16). What differs from the design, and why:
- Door 1 is a free function of pair records and a block's numbers, not `Session.formal_levels(x)`: the session's
  object is in files another agent rewrites (`api/src/pybind_api.cpp`, `AZUREAPI*.cpp`). `from_session` gives it
  the session's fit vector (Brune or Park amplitudes). A session in the formal basis is refused (its B_c belongs to
  the energies it was opened with).
- The particle channel functions of the doors repeat the 12 lines of `ChannelFunctionsFromCNuc`
  (`src/CNucLevelParameters.cpp`, another agent's file) for a `PPair`; that function could call the door's.
- `to_project` asks door 2 for a Brune tolerance of 1e-11 MeV (the engine's report uses 1e-6 MeV): with 1e-6 the
  write-then-read of T1 moved the model by 7.6e-5.
- IFG=0 is converted by `to_project` (section 13.7); KRM=3 is refused (section 13.7: the door it needs).
- The coupling order of SCH is not in the manual; read as AZURE2's (light particle first). No ENDF-6 section can
  settle it: unpolarized cross sections and distributions do not depend on it (section 13.9, measured).

Measured (`tests/pyazr/endf_model_test.py`, build `work/b-endf3`): door 1 on a Brune session's start equals the
formal session's levels bit for bit, door 2 on those with door 1's B_c equals the formal session's transform bit
for bit (T1, T2, 13N). The ENDF/B-VIII.1 extract: 14 levels (the design's 15 was a miscount), the observed energies of
section 6.2, MF=3/MT=600 above 7.2 MeV reproduced to 1.13e-6 (median 6.6e-8), 6Li 21.7 keV below AME2020 in the
report. Write then read (from_session, MF=2 in the widest form, resonances(), to_project, the levels put back into
the original session): T1 parameters to 5.8e-9, model to 1.7e-6 (one analysing-power point beside the narrow 3/2-
level; the 11-column fields carry BND to 1e-9 and ER to 1e-3 eV), Delta chi^2 -2.2e-4; T2 parameters to 1.2e-8,
model 1.7e-7, Delta chi^2 +6.2e-4. Ten refusals, each with its three parts.

### 13.5 Stage 3: the writer

`pyazr/endf/writer.py`: `write(session, path, x=None, header=, energy_range_MeV=, tolerance=1e-3, max_points=4000,
number_form="widest")` -> `Written` (sections, grids, the MF=2 check, losses, not written); `records.RecordWriter`
writes the records. Command `pyazr endf write PROJECT NEW`. MF=1/MT=451 (LRP=2; AWI the projectile's nuclear mass
from AME2020, p. 24; NSUB; the text: what wrote it, the MF=2/MF=3 sentence, the photon widths, the losses; the
directory), MF=2/MT=151 (door 1), MF=3 of every particle exit and of the total capture.

The grid: `pyazr.refine_grid` (the engine's refinement of an adaptive extrapolation grid; the engine's `energies
adaptive` record is refused by a pyazr session) on all reactions at once, each sweep one extrapolation session at
the sweep's lab energies. Its rule: |midpoint - linear interpolation| <= tolerance x max(|sigma|, 1e-3 max sigma).
Measured (`tests/pyazr/endf_write_test.py`): MF=3 read back with the module's reader at the data's lab energies and
at the grid's midpoints equals the engine to 9.2e-4 - 1.0e-3 by that rule (T1 MT=102 381 points, T2 MT=800/801 427
points, 13N MT=102 488 points). Before writing, MF=2 is read back (to_project) and its particle reactions computed
at the grid: T2 7.9e-8 (T1 and 13N have no particle exit: nothing to compare). The own format check: no error, the
LRP=2 warning.

MF=6 and MF=3/MT=2: section 13.8. The real codes on the written files: sections 13.6 and 13.8.

### 13.6 The checking codes on the written files

CHECKR 8.24, FIZCON 8.19, PSYCHE 8.11 (`opt/`), on T1, T2, 13N as written:
- CHECKR: LRP=2 (error; the PI's decision), MAT should be 625 / 728 (warning; MAT 9999 is the default, the
  evaluator's field).
- FIZCON: MF=3/MT needs MF=6 (one per MF=3 section: MF=6 not written); on 13N before AWI was taken from AME2020,
  "AWI SHOULD BE SET TO 9.98623E-01" (13N gives the proton's mass as 1.00728 u, nuclear; repaired).
- PSYCHE: AP not in range and, on T2, "SPIN OF ISOTOPE 7015. DISAGREES WITH 0.5" (both its reading of the record
  after the range as SPI, AP of LRF 1-3, section 13.2: that record holds 0, 0, IFG, KRM, NJS, KRL); for MT=102
  "calculated Q 6.67452E+06 disagrees with the given Q 1.94349E+06": PSYCHE's `QTEST` (psyche.f lines 2700-2820)
  takes the masses from its function `AMASS` (lines 8874-9046), a semi-empirical mass formula with shell terms and
  a table for n, p, d, t, 3He and alpha only; that formula, coded in `work/endf4/s/amass.py`, gives 6.6745 MeV for
  12C + p -> 13N. The file's QM = QI is AME2020's S_p of 13N (1.94349 MeV): the code's number, not the file's (the
  test of PSYCHE is 5 % and above 0.25 MeV, so any light capture target is reported); a negative thermal cross
  section and a resonance integral (neutron quantities, below the table).

### 13.7 The writer's loose ends (as built)

- **The start grid** (`pyazr/endf/grid.py`, `start_grid`): 41 even points over the range, every exit threshold and
  every level energy with E_r -+ Gamma/2 inside it. The lists are pyazr's: the physical values of the session's
  parameters at x (`Session.transform`: energies and the widths of the channels open at the level) and the pair
  records (`Session.pairs`), E_x to lab by the engine's `cm_to_lab_energy`. `grid.refine` runs `pyazr.refine_grid`
  and then checks the quarter points of every interval by the same rule, refining again from the failing ones (up to
  four rounds; the model is cached, so a point is computed once). Without the check, the seeded grid of 13N reached
  3.6e-3 on a flank (44 of 20001 points above 2e-3); with it T1 MT=102, T2 MT=800/801 and 13N MT=102 are within
  1.0e-3 on 20001 points each. Test: a copy of T2 with its 3- level at 13.0734 MeV made 0.2 eV wide, the level
  2.6 keV from an even point: the grid follows the engine to 1.8e-4 of the peak; with the seeds removed it misses it
  (0.92, the injected error). A file of that copy is not written: its MF=2 reproduces MF=3 only to 3.7e-4 near the
  level, since ER holds 1e-3 eV at 1 MeV (the 11-column field) against a width of 0.2 eV; the writer's 1e-5 check
  refuses it (a question for the PI: a narrow level cannot be written).
- **The loss list**: `from_session`'s losses name only what the project holds (photon pairs, data sets with their
  records, THM experiments, beta pairs, excited light particles); `Written.not_written` names every section (MF=3,
  MF=6 of elastic MT=2, of each exit, of MT=102) a reaction of the project needs that the file lacks, with the
  reason; both are in the MF=1 text. The test compares them with the file's sections and the project.
- **Refusal texts**: the fix of `Candidate.accept`, `write`, `to_project` and `from_session` is in words (no Python
  call), so that the page, the command line and a script read the same.
- **IFG=0** (p. 73): converted by `to_project`, gamma = sign(GAM) sqrt(|GAM| / 2P) (the manual prints
  "sign(GAM) x |GAM|/2P": its square root is meant, Gamma = 2 gamma^2 P, eq. D.45 p. 356), P the engine's
  `penetrability` row at the channel's c.m. energy |ER| m_B/(m_A+m_B) + Q (eq. D.82, p. 364; |ER| for a negative
  ER, p. 73); a width on a channel closed there is refused. T2 written as IFG=0 gives the levels of the IFG=1 file to
  6.9e-7. **KRM=3** stays refused: AZURE2's `use_rmc` holds eliminated photon widths, but a start model is in the
  Brune basis, and the Brune transformation of the eliminated widths needs a door the engine lacks: door 2
  (`observed_levels_of_block`) with, per level, the formal eliminated width Gamma_gamma (MeV, P = 1) as an
  argument, returning the observed Gamma_gamma of each level (the eigenvector transformation of the particle
  channels applied to the photon amplitudes). No library file of `work/endf/evals` exercises either: the only one,
  ENDF/B-VIII.1 3He + 4He, is IFG=1, KRM=4.
- **PSYCHE's Q of MT=102**: its mass formula (section 13.6); nothing to repair in the file.

### 13.8 MF=6 (as built)

Files: `pyazr/endf/engine_values.py` (the engine's values on the grids: `Reaction`, `Sweeps.energy_grid`,
`Sweeps.elastic_cosines`), `pyazr/endf/sections.py` (the records: MF=3, LAW=2 with its LAW=4 recoil, LAW=5 with
MF=3/MT=2, MF=6/MT=102), `pyazr/endf/grid.py` (`refine`, `thin`), `writer.py` (`_reactions`, the assembly). The
signatures of `write` and `to_project` are unchanged; `write` has one keyword more, `elastic_max_cosine` (0.98).

- **LAW=2** of every particle exit (p. 143): LCT=2, the emitted particle (LAW=2, LANG=0, yield 1) and the residual
  (LAW=4, yield 1); a_l per energy are the engine's `legendre` quantity with the expansion `rsl`, A_L/((2L+1)A_0),
  which is LAW=2's normalisation (eq. 6.8): no arithmetic in Python. Orders asked up to 2 l_max + 4; the last must
  come back zero (else refused). NL per panel is the highest non-zero order.
- **LAW=5** of charged-particle elastic scattering, LTP=12 (the design's choice: no new door; pp. 144-147) with
  MF=3/MT=2 holding sigma_NI (p. 150: "With LPT=12 and LPT=14 representation, MT=2 in File 3 contains
  sigma_NI(E) in barns"). sigma_e from the engine's `dsigma_domega` at c.m. angles, sigma_c by the reader's own
  terms (`CoulombTerms`, the engine's rows on the file's AWI and AWR) so that the reader gets the engine's sigma_e
  back at the nodes. sigma_NI = the trapezoid integral of sigma_e - sigma_c over the panel's cosines (exact for
  P_NI linear in mu), P_NI = (sigma_e - sigma_c)/sigma_NI. mu from -1 (0 for identical particles, LIDP=1, SPI the
  projectile's spin) to 0.98; above it the file gives the Coulomb term alone (listed in MF=1).
- **MF=6/MT=102** (p. 136): per final state (photon pair) in order of excitation, the primary photon (ZAP=0, LAW=2,
  yield sigma_f/sigma of MF=3/MT=102, a_l of the engine's capture distribution to that state) and the compound's
  recoil (LAW=4, the same yield). Not written: the cascade below the final states (MF=1 says so).
- **Masses** (p. 24): AWP of light products (A <= 4) nuclear, of heavier ones atomic, AME2020; AWI the same function.
- **The grid**: one energy grid for every section (MF=3 and MF=6), `grid.refine` from `grid.start_grid`, the rule
  applied to sigma of each exit and of capture, the yield of each final state, p(mu) of each LAW=2 distribution at
  9 cosines and sigma_e at 9 cosines of the LAW=5 range. The cosines of LAW=5: `grid.refine` on sigma_e - sigma_c
  plus a constant per energy (the smallest sigma_e of the start cosines: it cancels in the chord, so the rule holds
  the error against the size of the cross section; on sigma_e - sigma_c alone T2 needed 2418 cosines per energy and
  an 81 MB file), then each energy keeps the cosines its own values need (`grid.thin`, the same rule; T2: up to 203).

What the manual leaves open, and the reading taken:
- AWP of a primary photon "can contain the energy of a primary photon" (p. 135), not at which incident energy: EG of
  MF=12 LP=2 (p. 201), the photon energy at zero incident energy, S - E_x,f.
- The yield of a primary photon and of its recoil: the multiplicity of a discrete photon, sigma_f/sigma (eq. 12.1,
  p. 200); LIP=0 for a recoil in an excited state (LIP marks isomers, p. 135).
- Eq. 6.19 integrates sigma_e - sigma_c over mu, not over the solid angle, while p. 146 says sigma_NI is "in barns";
  the equations are taken as written (sigma_e = sigma_c + sigma_NI P_NI), as the reader computes it.
- Between incident energies the reader follows the manual (p. 146): sigma_NI(E) of MF=3 times P_NI(mu, E) of
  MF=6, each interpolated under its own law; the product reading (sigma_NI(E_i) P_NI(mu, E_i) interpolated) is a
  named option. Where sigma_NI changes sign the two differ by about sigma_e - sigma_c; the writer brackets every
  sign change (section 13.10).
- The recoil of LAW=5 for distinguishable particles: LAW=4 (two-body, p. 144); the library file writes LAW=0 there.

Measured (`tests/pyazr/endf_write_test.py`, build `work/b-endf4`, tolerance 1e-3): `Material.curves` of the written
file against the engine's dsigma/dOmega, relative to max(|engine|, 1e-3 of its largest value), (i) at the points
of the project's differential data sets inside the file's range and LAW=5 cosines, (ii) at 60 midpoints of the
energy grid x 16 c.m. angles 20-170 deg:

| project | elastic LAW=5 (i), (ii) | exits LAW=2 (i), (ii) | capture photons (i), (ii) | grid |
|---|---|---|---|---|
| T1 | 1.8e-3 (207 points), 5.6e-4 | - | 5.5e-4 (21), 8.3e-4 to the ground state; 3.0e-3 to the 1.0 MeV state | 736 energies, up to 267 cosines |
| T2 | no data, 2.5e-4 | 1.0e-3 (60), 1.1e-3 (MT=800); 9.9e-4 (MT=801) | - | 655 energies, up to 203 cosines |
| 13N | 3.4e-3 (136), 7.8e-4 | - | 1.1e-3 (67), 1.4e-3 | 740 energies |

So the stated tolerance: each tabulated quantity to the rule (1e-3); dsigma/dOmega, a product of two tables
interpolated each, to 3e-3 for LAW=2 and the photons and 5e-3 for elastic scattering (the test's bounds). The grid
grows with MF=6 (T1 736 energies against 381 with MF=3 alone; T2 655 against 446): the elastic cross section at 9
cosines is the column that asks most. The write takes 47 s for T2, about 2 min for T1 (the example).

The real codes on the written files (CHECKR 8.24, FIZCON 8.19, PSYCHE 8.11; T1, T2, 13N), every finding:
- CHECKR "LRP=2 ONLY ALLOWED IN DERIVED FILES" (the PI's decision) and "MATERIAL NUMBER SHOULD BE 625/728" (MAT
  9999, the evaluator's field): as before.
- CHECKR "NE = 655 (736, 740) OUT OF RANGE 1 - 200" in MF=6, then "SECTION CANNOT BE CHECKED", then "OUT OF
  SEQUENCE AT SEQUENCE NUMBER 1" of the next MF=6 section: CHECKR applies MF=5's limit of 200 incident energies
  (`NES5MAX`, checkr.f lines 337, 3419, 3459) to LAW=2 and LAW=5, where Appendix G (p. 391) gives MF=6 NE 2,000;
  the third follows from the second (a file written over 0.6-0.9 MeV, 48 energies, has none of them). The code's.
  With 48 energies CHECKR reports "NL = 208 OUT OF RANGE 1 - 201" for a LAW=5 panel: MF=4's 201 angular points
  (p. 391 gives no limit for the cosines of MF=6 LAW=5); the code's. The library file has the NE finding too. Whether
  to keep MF=6 under 200 energies and 201 cosines for the codes that copied these limits is a question for the PI.
- CHECKR "EMAX = 9.00E+05 OUT OF RANGE 1.00E+06 - 3.00E+09" on that 0.6-0.9 MeV file: CHECKR's bound; EMAX is the
  file's upper energy.
- FIZCON "DISCRETE 2-BODY LAW NOT PERMITTED FOR MT= 102" (T1, 13N): fizcon.f lines 4343-4349 allow LAW=2 only for MT
  2, 50-90 and 600-848; p. 143 allows it for the primary photons of MT=102. The code's.
- FIZCON "CHECK NORMALIZATION= 0.998889 AT E= 2.7275E+05" (T1 MF=6/MT=2): at that energy sigma_NI = -6.4e-5 b and
  P_NI reaches 2e4; the panel as written integrates to 1.0000001 (the module); FIZCON sums the intervals in single
  precision (fizcon.f lines 4418-4436). The code's arithmetic; the panel is the fragile case of the reading above.
- FIZCON "MF=3/MT REQUIRES MF=6": gone (MF=6 written).
- PSYCHE: "AP NOT IN RANGE", "SPIN OF ISOTOPE" (its reading of LRF=7, section 13.6); the calculated Q of MT=102 (its
  mass formula, section 13.6); thermal cross sections, resonance integrals and boundary tests of MT=2 and MT=102
  (neutron quantities); "NO ENERGY-BALANCE TEST FOR ..." and "NO GAMMAS" (notes of what PSYCHE does not test).

### 13.9 The channel-spin coupling order (measured)

The full ENDF/B-VIII.1 3He + 4He file (`work/endf/evals/h_He4/h-002_He_004.endf`, read with `python -I`;
`work/endf4/s/sch_order.py`): `to_project` of its MF=2, then the engine's 4He(3He,p)6Li at 13 panel energies of
MF=6/MT=600 from 7.2 to 19.2 MeV and 9 angles.
1. The file's LAW=2 coefficients are not normalised as eq. 6.8 says (p. 143: "the angular distribution p_i is
   normalized"): they are the engine's a_l times sigma(E)/(2 pi) of its own MF=3/MT=600, to 5 digits, at every
   energy up to 17.2 MeV (a defect of the file, or of its translation by FUDGE; a reader following the manual draws
   the distribution wrong by that factor). Divided by sigma/(2 pi), p(mu) of the file equals the engine's to 4.2e-4
   of its maximum (7.2-17.2 MeV; below 1e-7 at 6 of the 11 energies); at 18.2 and 19.2 MeV it differs by 52 % and 36 %
   (not explained: the coefficients' ratios to the engine's differ by order). The l-dependent signs of the start
   model are right where the file is consistent.
2. Every amplitude of channel spin 1/2 of p + 6Li flipped in sign in all levels (14 amplitudes, a copy of the
   project): sigma, dsigma/dOmega and every Legendre coefficient unchanged to the last bit (0.0). A control, one
   amplitude of one level flipped, changes dsigma/dOmega by 17 %. Exchanging the two spins in the coupling changes a
   channel's state by the constant sign (-1)^(I1+I2-s), the collision matrix by a sign per channel pair; unpolarized
   cross sections sum over the channel spins incoherently, so no MF=3 or MF=6 number depends on the order, and no
   library file can settle it. Analysing powers mix channel spins coherently and do depend on it: a project made by
   `to_project` and fitted to analysing powers rests on AZURE2's order (the light particle first), which the
   conversion report says when a pair has two non-zero spins.
3. The reader reads MF=6/MT=600 (LAW=2, LANG=0) and MF=6/MT=2 (LAW=5, LTP=1) of this file.

### 13.10 The elastic section under the manual's reading

**The two readings.** ENDF-102 p. 146 (eqs. 6.19, 6.20): "The angular distribution p_NI is given in File 6 as a
tabulated function of mu, and sigma_NI(E) in barns is given in File 3"; each record has its own interpolation, so
a processing code interpolates sigma_NI(E) in MF=3 and P_NI(mu, E) between the panels of MF=6 at fixed mu, and
multiplies (the *manual* reading). The manual states no other reading. The module had interpolated
sigma_NI(E_i) P_NI(mu, E_i), the tabulated sigma_e - sigma_c (the *product* reading). Where sigma_NI changes sign
between two energies, P_NI is large and of opposite sign on the two sides and the manual reading misses
sigma_e - sigma_c by about its own size, however close the energies; next to a zero of sigma_NI, where its ratio
between neighbours r is far from 1, the error is about (r - 1)^2 / (4 r) of sigma_e - sigma_c.

**Measured on the files of the previous writer** (`work/endf5/s/readings.py`: the file's elastic dsigma/dOmega
against the engine at the quarter points of every interval of MF=6/MT=2's energies, 12 c.m. angles 20-175 deg; the
deviation relative to max(|engine|, 1e-3 of its largest value), as the writer's test; build `work/b-endf5`):

| file | energies | sign changes of sigma_NI | product reading: max, median | manual reading: max (where), median | width off by > 1e-3 / 5e-3 / 1e-2 / 0.1 under the manual reading |
|---|---|---|---|---|---|
| T1 | 736 | 4 | 8.9e-4, 2.0e-5 | 0.79 (1.7361 MeV, 20 deg), 2.5e-5 | 62.1 / 24.7 / 19.6 / 4.8 keV (53 intervals) |
| T2 | 655 | 3 | 2.4e-4, 7.7e-7 | 2.9e-3 (0.3847 MeV, 30 deg), 9.0e-7 | 13.0 keV / 0 / 0 / 0 |
| 13N | 740 | 6 | 8.3e-4, 8.4e-6 | 0.56 (1.7338 MeV, 20 deg), 1.1e-5 | 42.4 / 20.3 / 14.8 / 2.9 keV (38 intervals) |

(Widths: the summed lab-energy spacing of the quarter points off by more than the bound.) Off by more than the
writer's tolerance are the intervals that hold a sign change and their neighbours, where sigma_NI is small.

**The library.** `work/endf/evals` holds one charged-particle file, ENDF/B-VIII.1 3He + 4He: MF=6/MT=2 is LAW=5
LTP=1 at 400 energies and MF=3/MT=2 holds 1.0 (p. 145), so it has no sigma_NI to change sign and says nothing
on the reading of LTP=12 (`work/endf5/s/lib_mt2.py`, `python -I`). No library file contradicts the manual's
text; none confirms it either.

**What is built.** `pyazr.endf.distributions.elastic_cross_section(..., reading="manual" | "product")` and
`Material.curves(..., elastic_reading=...)`, the manual's by default; a curve's notes say where the other reading
differs by more than `curves.READINGS_BOUND` (1e-3). The writer (`pyazr/endf/elastic.py`) gives elastic
scattering its own energies, chosen so that the manual reading holds the grid's rule:
1. every sign change of sigma_NI on the grid is closed in by ladders of energies E_0 -+ (w/2) 1.05^k around the
   chord's zero until the two energies that hold it are at most `BRACKET_WIDTH` = 1e-5 x E apart (the bracket);
2. the energies written are chosen from all the engine computed by `grid.thin_to` with the manual reading between
   two kept energies instead of the chord, the two ends of every bracket kept;
3. the quarter points of every kept interval that holds no computed energy are computed and the choice made again,
   until they all pass or the energies computed reach twice `max_points`.
MF=1, `Written.elastic` and its `print()` state the number of energies, the tolerance reached, whether every
interval was checked, and each bracket (inside which the manual reading is wrong by up to sigma_e - sigma_c).

**Measured on the files written so** (same script, quarter points of the new intervals):

| file | energies (computed) | brackets (eV wide) | checked | manual reading off by > 1e-3 / 5e-3 / 1e-2 / 0.1 | product reading max | write |
|---|---|---|---|---|---|---|
| T1 | 3551 (7987) | 4 (0.33, 0.14, 1.5, 15.7) | budget reached | 567 / 22.6 / 13.3 / 13.0 eV | 9.2e-4 | 293 s |
| T2 | 415 (1908) | 5 (0.4-3.8) | yes | 10.2 / 4.3 / 0 / 0 eV | 9.3e-4 | 74 s |
| 13N | 3174 (7965) | 6 (0.3-17) | budget reached | 337 / 14.5 / 14.5 / 14.3 eV | 9.2e-4 | 330 s |

Outside the brackets and a few eV around them the file holds 5e-3 under the manual reading; the regions off by
more than 1e-3 lie within 2.6 keV of a zero of sigma_NI. The cost: T1 and 13N need 3000-3500 elastic energies
(about 20 MB files) where the product reading needed 736. The format allows it (Appendix G, p. 391: NE 2,000 in
MF=6, exceeded here); the limits the PI chose for written MF=6 are the next section's.

### 13.11 CHECKR's limits on written MF=6 (the PI's decision)

The PI decided that written MF=6 stays within CHECKR's limits: at most 200 incident energies per section
(`checkr.f` NES5MAX, File 5's limit applied to MF=6) and 201 cosines per LAW=5 panel (NANGMAX, File 4's). The
manual's own limits (Appendix G, p. 391) are 2,000 energies in MF=6, no MF=6 limit on cosines (the general mesh
limit NP 10,000), 50,000 points in an MF=3 TAB1; CHECKR has NPTSMAX = 50,000 for a TAB1 too, so MF=3 is not
limited (the written MF=3 sections have 655-740 points).

`write(..., mf6_limits="checkr")` (default) or `"manual"` (`writer.MF6_LIMITS`; another value is refused with
its three parts). A section over the limit is thinned by the grid's rule to the limit (`grid.thin_to` with a
limit: the interval with the worst dropped point is split first), and the tolerance it reaches is stated: MF=1
("Accuracy:"), `Written.accuracy` (lines; they also end `Written.losses`, which the export screen shows without a
change of its code) and `Written.reached` ({MT: tolerance}). The exits (LAW=2) and MT=102 keep MF=3's grid up to
the limit and a subset above it; between MF=6 energies a reader interpolates p(mu) at fixed mu linearly (TAB2
INT=2) and multiplies by MF=3's sigma at the energy asked, so the error is the stated one of p(mu) plus MF=3's.
The elastic section chooses its 200 among the energies it computed by the manual's reading (13.10); the brackets
of the sign changes are kept first (a file whose brackets alone exceed the limit is refused); each panel keeps at
most 201 cosines.

Measured (build `work/b-endf5`, one thread each, three writes side by side on a loaded machine;
`work/endf5/s/readings.py` for the elastic column, as in 13.10):

| file | write | MF=6 sections thinned: tolerance reached | elastic: 200 energies; stated tolerance; measured max / median (12 angles, quarter points) | cosines |
|---|---|---|---|---|
| T1 | 100 s | MT=102 2.7e-3 | 5.5; 0.84 / 1.5e-3, off by > 1e-2 over 497 keV | 1.6e-3 |
| T2 | 51 s | MT=800 3.4e-3, MT=801 (144 energies) none | 0.12; 1.9e-2 / 2.2e-5, > 1e-2 over 8 eV | 1.0e-3 |
| 13N | 133 s | MT=102 (147 energies) none | 1.3; 0.47 / 6.7e-4, > 1e-2 over 350 keV | 1.5e-3 |

The stated tolerance is the largest deviation by the grid's rule at the energies the engine computed and at the
quarter points and midpoint of every written interval (one more session), cosines included: a single angle where
the cross section is small near a zero of sigma_NI sets it. So under CHECKR's limits the elastic section of a
proton on 12C (T1, 13N) is not usable near the zeros of sigma_NI under the manual's reading (the 200 energies
cannot follow sigma_NI's ratio near four to six zeros; the median stays at 1e-3) and MF=1 says so. With `mf6_limits="manual"` the section takes what 13.10 measured. T2
meets its tolerance outside the brackets.


### 13.12 Elastic scattering as nuclear amplitudes (LAW=5 LTP=1)

**The door.** The engine's elastic cross section at a point is already of the manual's LTP=1 form. The collision
matrix of the point's energy (the T-matrix of `AMatrixFunc::CalculateTMatrix` or `RMatrixFunc::CalculateTMatrix`,
read through `GenMatrixFunc::GetTMatrixElement`) is turned by `GenMatrixFunc::ElasticCoefficientsOf`
(`src/GenMatrixFunc.cpp`) into the site's kElastic record (`include/SiteObservable.h`, `elastic_record`): the
nuclear B_L of |F_N|^2, the complex Coulomb-nuclear D_l (X_l for identical particles with spin) and the Coulomb
factors, which `ContractElastic` (`src/SiteObservables.cpp`) contracts with the engine's Coulomb amplitude
f_d = eta / (sqrt(pi) (1 - mu)) exp(-i eta ln((1 - mu)/2)) (`Site::CoulombAmplitudesAt`). Eq. 6.13 (p. 145) is the
same sum term by term with

  a_l = -conj(D_l) / (sqrt(pi) (2l + 1)),   b_l = 2 B_l / (2l + 1);

eq. 6.14 (identical particles) with the same a_l where X_l = (-1)^l D_l (every channel of two identical particles
has l + s even, so the exchange sign (-1)^s is (-1)^l) and b_l = 2 B_{2l} / (4l + 1), the odd B_L zero (both
checked by the door; a record that breaks them is refused). The spins enter only through the engine's T-matrix and
its spin average: the a_l expand the trace of the amplitude matrix over the spin states divided by their number
(the Coulomb amplitude is diagonal in them), which is what D_l holds. The door is a point output of an elastic
data set, `nuclear_amplitudes` (`api/src/pybind_elastic_amplitudes.cpp`, one row of `PointOutputs()` in
`api/src/pybind_evaluation.cpp`): the LIST values of p. 147 per point, read from the record the evaluation wrote;
nothing is computed again and nothing of the evaluation changes. It is a data set's output, so the writer asks it
of a data session (`Sweeps.nuclear_amplitudes`: one data set at 90 deg c.m. with placeholder values, every energy a
row); the extrapolation path (`calculated_*` of `AZUREAPI`) is another agent's file.

**The check of the door** (`work/endf6/s/door_check.py`, eqs. 6.9-6.14 written from the page, CODATA 2018
constants; 50 lab energies over the data's range x 12 c.m. angles 10-170 deg): the cross section rebuilt from the
coefficients against the engine's dsigma/dOmega at the same points, largest relative difference T1 7.8e-12, T2
4.6e-12, 13N 7.4e-12, p + p (`identical_pp_res`, identical fermions) 1.3e-11, d + d (`identical_dd_mott`, spin 1)
4.5e-12. The constants of eqs. 6.11-6.12 differ from the engine's in the last digits; nothing else is left. The test
`tests/pyazr/endf_amplitudes_test.py` (20 energies) holds it to 1e-9; the mapping injected wrong (conj dropped)
fails it on all three projects.

**Against the library.** ENDF/B-VIII.1 3He + 4He (LTP=1 at 400 energies; `work/endf6/s/library_ltp1.py`, read
with `python -I`): the file's MF=2 read into a project (`to_project`), the door's coefficients at 16 of the file's
panel energies against the file's own. At 1 MeV they agree to 4e-6 of the largest coefficient (a_0 =
-1.443711e-2 + 4.89646e-4 i in both), which confirms the convention (sign, sqrt(pi), the conjugate). The agreement
falls off smoothly with energy: 9e-5 at 2.25 MeV, 2e-3 at 4.75, 3e-2 at 9.75, 2.5e-1 at 19.75 MeV (the file's
cross section by eqs. 6.9-6.13 against the engine's from its own MF=2: 3e-6 at 1 MeV, 4e-2 at 7.25, 6e-1 at 14.75).
Not explained: the file's elastic section and its MF=2 agree only at low energy (as its MF=6/MT=600 above 17 MeV,
13.9); the file states NL = 6 where its MF=2 has l <= 3 (the extra coefficients are zero or small).

**The writer.** `write(..., elastic_form="amplitudes" | "tabulated")`, the default the named constant
`writer.ELASTIC_FORM` ("amplitudes"; the PI decides at the merge). "tabulated" is LTP=12 as before (13.10,
13.11). "amplitudes" (`pyazr/endf/amplitudes.py`, `sections.elastic_amplitudes`): MF=3/MT=2 holds 1.0 at the first
and last energy (p. 145: "a value of 1.0 is used in File 3"; the library file does the same), MF=6/MT=2 LAW=5 with
per energy a LIST `0, E, LTP=1, 0, NW, NL` of the door's values. The energies: the coefficients at every energy of
the one grid, then `grid.thin_to` on the cross section at 17 cosines over the LAW=5 range with the reader's own
interpolation between two kept energies (sigma_e - sigma_c linear at fixed mu, the Coulomb term at the energy),
quarter points of the kept intervals computed and the choice made again while the limit allows, and the quarter
points and midpoints of every written interval computed once more for the tolerance stated (MF=1, `Written.elastic`,
an `AmplitudeReport`). The rule holds the engine's cross section (the Coulomb terms on the project's masses); the
reader's Coulomb terms use the file's masses (AWI nuclear, p. 24), which is in the measured deviation below.

What the manual leaves open, and the reading taken:
- The quantity a reader interpolates between the incident energies of LTP=1 (and 2): the manual says nothing; the
  module interpolates sigma_e - sigma_c at fixed mu, the Coulomb term computed at the energy asked (as 13.1).
- NL of a LIST of LTP=1: "the highest Legendre order of nuclear partial waves used" (p. 147): the highest l of the
  record's D_l, raised where the nuclear B_L ask more (2 NL + 1 values of b for distinguishable particles).
- The a_l normalisation: "complex coefficients for expanding the trace of the nuclear scattering amplitude matrix"
  (p. 145) is not more precise; the reading is fixed by eqs. 6.13-6.14 themselves (the coefficients that make the
  formula the cross section) and confirmed on the library file at 1 MeV.

**Measured** (`work/endf6/s/measure.py`: the file read back with the module's reader, the manual's reading for
LTP=12, against the engine at the quarter points and midpoint of every interval of the MF=6/MT=2 energies x 12 c.m.
angles 20-175 deg, relative to max(|engine|, 1e-3 of its largest value); build `work/b-endf6`, one thread each,
three runs side by side):

| file | form, limits | energies | stated | measured max (where) / median | width off by > 1e-3 / 5e-3 | write | file size |
|---|---|---|---|---|---|---|---|
| T1 | amplitudes, checkr | 200 | 6.4e-3 | 4.3e-3 (0.494 MeV, 60 deg) / 8.6e-5 | 763 keV / 0 | 28 s | 166 kB |
| T1 | tabulated, checkr | 200 | 5.5 | 0.84 (0.426 MeV, 90 deg) / 1.5e-3 | 1.26 MeV / 940 keV | 100 s | 1.1 MB |
| T1 | amplitudes, manual | 491 | 1.1e-3 | 8.5e-4 / 2.6e-5 | 0 / 0 | 29 s | 468 kB |
| T2 | amplitudes, checkr | 117 | 1.1e-3 | 1.8e-4 / 3.3e-6 | 0 / 0 | 1.5 s | 150 kB |
| T2 | tabulated, checkr | 200 | 0.12 | 1.9e-2 (0.687 MeV, 45 deg) / 2.2e-5 | 169 keV / 13 eV | 51 s | 845 kB |
| 13N | amplitudes, checkr | 200 | 6.0e-3 | 3.4e-3 (1.742 MeV, 30 deg) / 3.1e-5 | 478 keV / 0 | 37 s | 145 kB |
| 13N | tabulated, checkr | 200 | 1.3 | 0.47 (1.734 MeV, 20 deg) / 6.7e-4 | 1.28 MeV / 544 keV | 137 s | 1.0 MB |
| 13N | amplitudes, manual | 465 | 1.1e-3 | 6.2e-4 / 5.2e-6 | 0 / 0 | 37 s | 380 kB |

(T2 amplitudes needs 117 energies under either limit. The write times include MF=3 and the other MF=6 sections; T1
and 13N spend most of them on MT=102.) The real codes (CHECKR 8.24, FIZCON 8.19, PSYCHE 8.11) on the amplitude
files report only what 13.8 lists for every written file (LRP=2, MAT, FIZCON's LAW=2 for MT=102, PSYCHE's reading of
LRF=7 and its notes); with the manual's limits CHECKR adds its NE > 200 for MF=6/MT=2 and MT=102. On the tabulated
files FIZCON adds "CHECK NORMALIZATION" for 81 (T1), 62 (T2) and 70 (13N) LAW=5 panels: the panels near the zeros of
sigma_NI.

**What is left for "tabulated".** On these three projects the amplitude form is better on every number: its error
is the interpolation between energies alone, with no product of two tables and no sign change to bracket, and it
holds every angle (LTP=12 stops at mu = 0.98). What "tabulated" still offers: a file for a code that reads only
tabulated LAW=5 (none is known here; the checking codes read both, and the library's own charged-particle LRF=7
file uses LTP=1), sigma_NI in MF=3/MT=2 as a number a user may want, and a form that holds a cross section the
engine computes but whose amplitudes are not the elastic record's (none today). Not measured: a processing code
(NJOY, the PREPRO codes) reconstructing either form.

**sigma_NI's sign from the thinning** (the tabulated form): sigma_NI of a panel is the trapezoid integral over the
cosines the panel keeps; where the true sigma_NI is near zero, its sign is decided by those cosines (T2: three sign
changes within 22 eV near 385.93 keV). A repair (sigma_NI from the full cosine grid) would break the normalisation
of the written P_NI (the trapezoid over the kept cosines must be 1, eq. 6.20), so it is not a few lines; not done.
The amplitude form has no sigma_NI.

### 13.13 Covariances: MF=32 LCOMP=1 for LRF=7

`pyazr/endf/covariance.py`: `formal_covariance(session, x, covariance=None)` -> `FormalCovariance` (with `jacobian`
and `fit_covariance`), `jacobian(session, x, record)`, `rows(record, covariance)` (the records),
`read_covariance(material)` -> `ResonanceCovariance` behind `Material.covariance()`; `write(...,
covariance=None | "hesse" | the fit's covariance (array) | FormalCovariance)`; the own format check counts the records of MF=32 (`_read_mf32`).

**What is computed** (the owner's answer of 2026-11-18 to decision 10, replacing HESSE in a formal-basis session):
C_32 = J C J^T. C is the fit's covariance of the session's free parameters in its own basis (Brune or Park) with
its own fixed parameters: HESSE of the session's objective at x by iminuit (errordef 1 from the objective,
symmetrised, refused when not positive definite), or the matrix the caller passes (`EngineFitResult.covariance`).
J = d(ER, GAM of MF=2)/d(free parameters) at the B_c the file writes (BND): a reader varies ER and GAM with BND held,
so the derivative holds B_c too. Per J^pi group J = D2^-1 D_b, D2 the Jacobian of door 2 (formal levels at B_c ->
Brune levels) at MF=2's levels, D_b that of the group's Brune parameters in the session's (identity in the Brune
basis; door 2 after door 1 in the Park basis), central differences with steps 1e-5 (|E| + 1) MeV and 1e-5 max(|gamma|,
1e-3) MeV^1/2, door 2 iterated to 1e-12 MeV; then the file's units (ER: lab factor x 1e6; GAM: the amplitude
scale). Columns of norms, shifts, target parameters and photon widths are zero (marginalised). A resonance is listed
when a parameter of it has a non-zero row of J; the other parameters of a listed resonance have zero rows (on T1:
the formal ER of the lowest level of the two groups whose lowest Brune energy is fixed: at B_c = S_c(E_1) it equals
that energy).

**What it changed on T1** (fitted by MIGRAD; the method of `98a4cf60` run on the same build against this one, both
with iminuit's HESSE): 14 parameters listed instead of 10 (the two background poles' resonances, whose formal ER and
GAM move with the free parameters of their groups, are now listed); on the 10 parameters both list, the elements
differ by up to 3.89 sigma_i sigma_j of the new matrix; the sigmas move by factors 0.33 (a background level's ER, 512
-> 1565 eV) to 2.2 (a GAM, 17.8 -> 8.0 eV^1/2), five within 5 %. No MF=32 file is committed as an expected result
(the test writes and reads its own), so no file was regenerated.

**The check** (`tests/pyazr/endf_covariance_test.py`): J against central differences of the transformation at the
written B_c (the formal levels at theta +- 3h found by Newton's iteration on door 2, to 1e-11): 1.1e-8 of each row's
largest element; J C J^T taken back to the Brune basis through door 2's Jacobian equals the fit's own covariance of
its free level parameters to 1.2e-6 sigma_i sigma_j (the formal-basis HESSE differed from it by 3.9). Door 1's own
derivative (B_c = S_c(E_1) following E_1, what `from_session` does at each x) gives a matrix 0.93 sigma_i sigma_j
from this one (printed by the test). iminuit's HESSE against the engine fit's covariance as C: 8.5e-3 sigma_i sigma_j.

**The records** (pp. 253-259): HEAD, isotope, range (LRU=1, LRF=7, NRO=0); `0, 0, 0, LCOMP=1, 0, ISR=0`;
`AWRI, 0, 0, 0, NSRS=1, NLRS=0`; one block: `0, 0, NJSX, 0, 0, 0`, per spin group with a listed resonance a LIST
`0, 0, NCH, NRB, 6*NX, NX` of ER and GAM per resonance, each resonance starting a new line as in MF=2; then the
upper triangle, row by row, LIST `0, 0, 0, 0, N, NPARB`. A resonance whose every parameter is fixed in the fit is not
listed (p. 255: "not all resonances need to be included"); a fixed parameter of a listed resonance has a zero row
and column, and MF=1 names it ("zero rows").

What the manual leaves open, and the reading taken:
- L1 of the LCOMP=1 flags record: p. 258 prints `0.0, 0.0, 0, LCOMP, 0, ISR`, p. 255 says "IFG is relevant only for
  LRF=7 and its meaning is the same as in File 2"; the LCOMP=2 record (p. 268) puts IFG there. Written as p. 258
  prints it (0); GAM is the amplitude of MF=2's IFG=1.
- NX of a spin group's LIST: "the number of lines per resonance" and "6 times the number of lines" (p. 258): taken
  as MF=2's LIST (p. 74), NX the lines of the whole list, each resonance starting a new line.
- A spin group carries no J in MF=32: the reader identifies it by NCH and by the listed ER and GAM against MF=2
  (p. 256: the parameters are listed "for positive identification"); names `"<J^pi>#<k>.ER"`, `"...GAM<c>"`.
- The covariance LIST header is printed `[MAT,2,151]` on p. 259: taken as MF=32 (a misprint).

**Tests** (`tests/pyazr/endf_covariance_test.py`, T1 fitted by MIGRAD, 45 s): J and the door-2 round trip (above);
symmetric, positive semi-definite, zero rows exactly where J's row is zero; written and read back
(`Material.covariance`) to 1e-6 of each row's largest element (measured 2.7e-10; the number form keeps 7-10
digits), the values exactly MF=2's; the own format check without error; the fit's covariance as `covariance=`
written and read back; a matrix of the wrong shape, a word other than "hesse" and HESSE of a session of the
extrapolations refused with three parts. The off-diagonal layout injected wrong (1.001 on every off-diagonal
element) fails the read-back check. The real codes on the written T1 file: CHECKR, FIZCON and PSYCHE report nothing
on MF=32 (their findings are those of 13.12). No library file of `work/endf/evals` has MF=32 (the only one, 3He +
4He, holds MF=1, 2, 3, 6): no library file was read.

**The finding that led to the change.** HESSE in a formal-basis session (built until 2026-11-18) fixed the
parameters the project marks fixed *in the formal basis*: a fixed observed energy or a fixed width became a fixed
formal E_lambda or gamma, another constraint than the fit's, so the covariance was not the fit's uncertainty mapped to
the formal basis (on T1 up to 3.9 sigma_i sigma_j from it, through door 2). The owner chose the fit's covariance
carried through the transformation.

**For `pyazr.priors`** (not changed): a part kind `"endf"` (fields: the file, the isotope) would read
`Material.covariance()`; its parameters are formal ER and GAM on the lab scale, while `to_project` makes a Brune
project, so the part needs the Jacobian of door 2 at the file's levels (as the test computes it) and the units
(MeV, MeV^1/2 c.m.), and names that map `"<J^pi>#<k>.ER"` to the new project's parameter names (the order of
`to_project`'s levels is MF=2's). The lines: `PART_KINDS` (one row), a `_endf_part` beside `_fit_part`.

### 13.14 The check of an evaluation

`pyazr/endf/evaluation_check.py`: `check_resonances(material_or_path, tolerance=1e-4, separation_energy="ame2020")`
-> `ResonanceCheck` (findings, parts with the numbers, notes; `print()`, `to_dict()`). The findings are the own
format check's `Finding` with `source` "pyazr.endf resonances", so that `pyazr.endf.check(path, resonances=True)`
puts them among the own findings and the page's Check screen lists them without a change of its code; the page's
job passes no `resonances` today: the one line is `pyazr/ui/evaluations.py` line 169 (`endf.check(...,
timeout=settings["timeout"])` would add `resonances=True`, or a setting of the job; another agent's file, not
changed). Command: an option of the existing check, `pyazr endf check FILE --resonances` (one check of a file, one
report; a separate command would print the own findings twice). The three parts:

- **reconstruction**: `to_project` of MF=2 (AME2020 separation energy of the incident pair, the file's Q values),
  the engine's sigma of each particle exit at the energies of its MF=3 section, the relative difference per point
  (`RECONSTRUCTION`: largest, median, points off by more than the tolerance, a warning when any is), MF=3 above the
  reconstruction by more than 3 orders of magnitude (`THRESHOLD`). LRP=1 (MF=3 holds backgrounds) is not compared
  (a note). Not built: MF=6 against the reconstruction (the LAW=2 distributions and, with the door of 13.12, LAW=5
  LTP=1: `work/endf6/s/library_ltp1.py` does it by hand for the library file).
- **masses**: each nuclide of MF=2 against AME2020 (`MASS` above 1 keV: the masses' frame is not stated, p. 72,
  so the atomic mass is the reference), each pair's Q against AME2020's (`Q`, a note within 1 keV) and against the
  Q its own masses give (`Q_MASSES` above 1 keV).
- **text**: every `LRF=`, `KRM=`, `IFG=`, `LTP=`, `LRP=` the MF=1 text names against the record that holds it
  (`TEXT`); a flag whose section is absent is said so.

On the committed extract (`tests/pyazr/endf_evaluation_check_test.py`): RECONSTRUCTION (97 energies, 2 off: 7.05
and 7.10 MeV, median 6.9e-8), THRESHOLD (37 orders at 7.05 MeV: 4.01e-9 b against 3.36e-46 b), MASS (6Li -21.7 keV),
Q note (-0.083 keV), Q_MASSES (-22.6 keV), TEXT (the text's LTP=2; the extract has no MF=6/MT=2, which its
directory lists): the design's three findings. On the full library file (`work/endf6/s/check_full.py`, `python
-I`): the same, and TEXT says "the record of MF=6/MT=2 holds LTP=1". The mass tolerance injected at 100 keV fails
the test (MASS and Q_MASSES gone).

### 13.15 Memory, progress, the project's name, the export's time, the warnings of a sweep

- **Memory.** Peak resident size per engine session of a tabulated T2 write (`work/endf6/s/mem_tab.py`, a wrapper
  of `Sweeps.evaluate`): 0.17-0.96 GB for the first 26 sessions, then one session of 1 858 890 points (the cosine
  refinement of LAW=5, `Sweeps.elastic_cosines`: every energy of the grid x every cosine in one session, about
  2.6 kB per point) at 4.85 GB. That session is now cut into sessions of at most `SESSION_POINTS` (200 000) points,
  as `elastic_values` already was: the tabulated T2 write peaks at 0.78 GB with every section but MF=1 identical to
  the write before (MF=1 differs by the header fields given); the tabulated T1 file is byte for byte the file the
  start commit's code writes (same arguments; peak 0.85 GB, 2:03 min). The quick writer test (T1 in the default
  amplitude form, and the narrow level) peaks at 0.22 GB (was 7.3 GB), 70 s.
- **progress**: `write(..., progress=f)` calls `f(done, total, words)` after each of `writer.WRITE_STEPS` (six:
  the levels, the grid, elastic scattering, the covariance, MF=2 against MF=3, the file); a step the file does not
  need is counted all the same. The page's job row would use it where it builds the call (`pyazr/ui/`, another
  agent's; not changed).
- **project_name**: the quick writer test writes T1 as "T1.azr for the tests" and finds the name in the tape's first
  line and in MF=1's text.
- **The export's time.** `pyazr/evaluate.py` `_legendre_rows` (the Legendre frame of an extrapolation, read for
  every grid) looped over the points with a `round()` each: 12.8 s for 5.6 million points without coefficients
  (the elastic grids). A frame whose points all have the same count is now cut by one reshape (0.8 s); the arrays
  are the loop's, value for value (`tests/pyazr/endf_export_support_test.py`: equal, varying, empty frames).
  Every caller reads the same arrays: the callers are `_forward_pass` (the only one).
- **The warnings of a sweep.** `Project.new(nuclear_part=)` (`pyazr/project.py`, the classmethod `new`) emptied the
  record sections and kept `<parameterSettings>` whole, so every row of a removed data set's parameter (13N: its
  16 shift offsets) made the engine warn W_PARAMETER_SETTING_UNKNOWN in each session of the writer (208 lines in a
  13N write). The rows whose name starts with a removed record's name and a dot are dropped with the records
  (`<parameterSettings>` and `<mcmc>`); the door check of 13N opens with none (16 before). The test checks the rows.
