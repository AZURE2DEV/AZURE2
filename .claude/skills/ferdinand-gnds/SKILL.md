---
name: ferdinand-gnds
description: Translate R-matrix evaluations between AZURE2 (.azr), GNDS (XML), ENDF, FRESCO/SFRESCO, EDA, AMUR, RAC and LaTeX with LLNL's Ferdinand, and read or inspect GNDS files with FUDGE — e.g. export an AZURE2 fit to GNDS for an evaluation, import a GNDS/ENDF R-matrix evaluation into AZURE2, or compare an AZURE2 fit with another code's parameters. Use whenever the task mentions Ferdinand, ferdinand.py, FUDGE, GNDS, a .gnds/.xml R-matrix evaluation, or converting an .azr to or from another R-matrix code's format.
---

# Ferdinand + FUDGE (R-matrix format translation, GNDS)

| | |
|---|---|
| Ferdinand | `/groups/rdeboer1/user/rdeboer1/R-matrix/ferdinand/ferdinand.py` (git clone of github.com/llnl/ferdinand, v0.50, checked out 2026-09-02 commit) |
| FUDGE | `/groups/rdeboer1/user/rdeboer1/R-matrix/fudge` (github.com/LLNL/fudge, release 6.13.0) |
| Python | `python/3.12.15` module; FUDGE's compiled pieces (`numericalFunctions`, `crossSectionAdjustForHeatedTarget`) are in the user site |
| environment | `~/.login` appends the fudge directory to `PYTHONPATH`; in a shell without it: `setenv PYTHONPATH /groups/rdeboer1/user/rdeboer1/R-matrix/fudge` (csh) or `export PYTHONPATH=...` (bash) |
| full option list | `python3 ferdinand.py -h`, or `ferdinand/README.md` |

Installed 2026-09-07; first real use (12C+p round trip) 2026-10-03. Only the
AZURE2 <-> GNDS path has been exercised. Extend this file when a real task
uses another path or option, not speculatively.

## Running

```bash
F=/groups/rdeboer1/user/rdeboer1/R-matrix/ferdinand/ferdinand.py
python3 $F fit.azr gnds                       # -> fit.azr+.gnds
python3 $F fit.azr+.gnds azure -o back.azr    # GNDS -> AZURE2
python3 $F fit.azr tex -o fit.tex             # standalone RevTeX document with particle and parameter tables
```

- Input format comes from the file suffix (`.azr` -> azure) or `-i <format>`.
- Default output name is the input name **plus** `+.<format>`
  (`12C+p.azr+.gnds`); use `-o` for a clean name.
- Output goes next to the input file (or to `-o`). Work on a **copy** of an
  `.azr`, never inside a live fit directory.
- Fast: a 30-line-level 12C+p model converts in ~4 s.
- Reads the `.azr` `<levels>` block. To export fitted values, update
  `<levels>` from `param.sav` first (see the azure2-eval skill).
- The log reports boundary condition, elastic pair, and "Setting CN zero
  energy from eSepE = ..." when writing AZURE2. **Check that eSepE equals the
  separation energy in the original `.azr`.**

## Known problems (found on 12C+p, 2026-10-03)

### 1. Mixed nuclear/atomic masses shift every level energy (GNDS -> AZURE2)

GNDS -> AZURE2 recomputes the separation energy from the pair masses instead
of keeping the `.azr` value. The 12C+p fits use a **nuclear** proton mass
(`1.00728`) with an **atomic** 12C mass (`12`). Ferdinand then got
Sp = 1.43580 MeV instead of 1.94351, and every level energy came back
**0.5077 MeV low**. The GNDS file itself was correct (pole energies matched;
ground state at E_lab = -2.1066 MeV). With both masses atomic (p = `1.007825`)
the round trip returned the energies to within ~50 eV.

As of 2026-10-03, only the 12C+p directories use `1.00728` (35 `.azr` files).
The other 2026 fits use atomic masses (n `1.00866`, α `4.0026`, p
`1.00782/1.00783`). Before converting, check token `t17` (projectile mass)
of the level lines for consistency.

### 2. AZURE2-only settings are lost or changed, even with consistent masses

Token positions in the level lines follow `AZURE2/include/NucLine.h`:
`levelJ levelPi levelE levelFix aa ir s l levelID isActive channelFix gamma
j1 pi1 j2 pi2 e2 m1 m2 z1 z2 entranceSepE sepE j3 pi3 e3 pType chRad g1 g2 ecMultMask`.

After `.azr -> GNDS -> .azr` on 12C+p:
- `levelFix` (t3) and `channelFix` (t10) are all reset to 0. **Every parameter
  becomes free.**
- On all 18 p+γ (capture-pair) lines: `pType` (t26) 10 -> 0, `ecMultMask` (t30)
  5 -> 0, `chRad` (t27) 0 -> 1.0, and `e2`/`sepE` (t16/t22) 0 -> ±5.345. AZURE2
  treats `pType` 10 as a capture pair and needs `ecMultMask != 0` for external
  capture. **The γ channel is no longer a capture channel and external capture
  is off.**

So a round-tripped `.azr` is **not** a drop-in replacement for the original.
Use the GNDS export for exchange and archiving. When importing into AZURE2,
restore fix flags and capture-pair columns from a known-good `.azr`, then
recalculate and compare χ² before trusting the result.

Check any round trip token by token, not by eye:

```python
def levels(f):
    b = open(f).read().split("<levels>")[1].split("</levels>")[0]
    return [l.split() for l in b.strip().splitlines() if l.strip()]
from collections import Counter
diff = Counter()
for x, y in zip(levels("orig.azr"), levels("back.azr")):
    for i, (p, q) in enumerate(zip(x, y)):
        try: same = abs(float(p) - float(q)) <= 1e-4 * max(1, abs(float(p)))
        except ValueError: same = p == q
        if not same: diff[(i, p, q)] += 1
for (i, p, q), n in sorted(diff.items()): print(f"t{i} {p} -> {q}: {n} lines")
```

## Reading GNDS with FUDGE

```python
from fudge import reactionSuite
rs = reactionSuite.ReactionSuite.readXML_file("fit.azr+.gnds")
print(rs.projectile, rs.target, rs.projectileFrame)          # H1 C12 lab
print(rs.PoPs["C12"].getMass("amu"))
RM = rs.resonances.resolved.evaluated                          # RMatrix
print(RM.boundaryCondition, RM.reducedWidthAmplitudes)         # Brune False
for sg in RM.spinGroups:
    E = sg.resonanceParameters.table.getColumn("energy", "MeV")  # pole energies, projectile (lab) frame
    print(sg.spin, sg.parity, E)
```

Pole energies in a Ferdinand GNDS are **lab-frame energies relative to the
entrance threshold**, not compound-nucleus excitation energies. For example,
13N Ex = 0 appears as -2.1066 MeV = -1.94351 × 13.00728/12.

The compound-nucleus mass Ferdinand writes is a placeholder (13N = `13.0`
amu exactly). Do not read Q values from it.

## Cross-section reconstruction (TensorFlow): particle channels OK with -g, no capture

TensorFlow 2.21.0 was installed in the user site on 2026-10-03 (CPU only,
~1.9 GB). The install only added packages and did not change numpy (2.5.3).
`reconstructCrossSections.py`, `reconstructCrossSectionsBatching.py` and
`reconstructLegendre.py` can now run. Set `TF_CPP_MIN_LOG_LEVEL=3` to hide
the harmless `cuInit` / "Could not find cuda drivers" messages (no GPU on the
login nodes).

```bash
python3 .../ferdinand/reconstructCrossSections.py fit.azr+.gnds -E 3.0      # lab energies
python3 .../ferdinand/reconstructCrossSections.py fit.azr+.gnds -E 3.0 -G   # c.m. energies, files named *_cs-G*
```

It takes ~9 s for 12C+p and writes one file per channel next to the input:
`<in>_cs-ch_p-to-g`, `_cs-ch_p-to-p`, `_cs-cap_p`, `_cs-reac_p`, and so on.
Each has columns E, σ; `#`/`@`/`&` lines are Grace-style headers.

**Ferdinand does not handle radiative capture (γ) channels.** It treats the
photon pair like a particle channel (log: `photon + N13 :0.5- , 1.94351
radii 1.0`). With the 12C+p Brune-basis fit, the γ channels moved **every
level**: elastic and (p,γ) peaks came out at E_cm ≈ 0.54 MeV instead of 0.424
(the 2.365 MeV 13N level). The user expected this ("I don't think that
Ferdinand knows how to deal with radiative capture reactions"); the test
confirmed it on 2026-10-03.

For particle channels, **convert without γ channels** (`ferdinand.py fit.azr
gnds -g`; the log says "Channel photon excluded"). γ widths are eV-scale, so
this does not change elastic or particle-reaction results. On 12C+p with
`-g`, the angle-integrated 12C(p,p) agrees with the `.azr`:
- 1/2⁺: E_cm 0.4236 vs 0.4239 expected; FWHM 32.5 keV vs Γp 34.2 keV.
- 5/2⁺: 1.5993 vs 1.6013.
- 3/2⁻ (broad, interfering with the 5/2⁺): a shoulder at 1.566 vs 1.557.
  Interference moves apparent maxima, so these small offsets are expected.

Do not use Ferdinand for capture cross sections; AZURE2 is the reference.

`reconstructLegendre.py` (angular distributions, `-A ang ang`) fails on this
fit with or without γ channels: TensorFlow `MatrixInverse ... Input is not
invertible` in `LM2T_transformsTF`. The cause is not known. Until it works,
there is no point-by-point differential comparison with AZURE2's
`AZUREOut_aa=1_R=1.out`. The checks above use resonance positions and widths
only.
