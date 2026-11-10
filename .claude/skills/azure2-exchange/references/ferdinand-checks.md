# Ferdinand round trips: checks and GNDS reading

## Token-by-token comparison of two `.azr` level blocks

Run this after any `.azr → GNDS → .azr` round trip. It counts, per token position, every value that changed by more
than 10⁻⁴ relative.

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

Token names, in order (`io/include/azure2/io/LevelRows.h`): `levelJ levelPi levelE levelFix aa ir s l levelID isActive channelFix
gamma j1 pi1 j2 pi2 e2 m1 m2 z1 z2 entranceSepE sepE j3 pi3 e3 pType chRad g1 g2 ecMultMask`.

## Reading GNDS with FUDGE

```python
from fudge import reactionSuite
rs = reactionSuite.ReactionSuite.readXML_file("fit.gnds")
print(rs.projectile, rs.target, rs.projectileFrame)          # e.g. H1 C12 lab
print(rs.PoPs["C12"].getMass("amu"))
RM = rs.resonances.resolved.evaluated                          # RMatrix
print(RM.boundaryCondition, RM.reducedWidthAmplitudes)
for sg in RM.spinGroups:
    E = sg.resonanceParameters.table.getColumn("energy", "MeV")  # pole energies
    print(sg.spin, sg.parity, E)
```

- **Pole energies in a Ferdinand GNDS are lab energies relative to the entrance threshold**, not E_x. Example:
  ¹³N E_x = 0 appears as −1.94351 × 13.00728/12 = −2.1066 MeV.
- **The compound-nucleus mass Ferdinand writes is a placeholder** (¹³N = 13.0 u exactly). Do not take Q values
  from it.

## Cross-section reconstruction

`reconstructCrossSections.py`, `reconstructCrossSectionsBatching.py` and `reconstructLegendre.py` in
`$FERDINAND_HOME` need TensorFlow (CPU is enough; set `TF_CPP_MIN_LOG_LEVEL=3` to hide the CUDA messages).

```bash
python3 $FERDINAND_HOME/reconstructCrossSections.py fit.gnds -E 3.0      # lab energies up to 3 MeV
python3 $FERDINAND_HOME/reconstructCrossSections.py fit.gnds -E 3.0 -G   # c.m. energies, files *_cs-G*
```

One file per channel is written next to the input (`<in>_cs-ch_p-to-p`, `_cs-ch_p-to-g`, `_cs-cap_p`,
`_cs-reac_p`, ...), columns E and σ, with Grace-style `#`/`@`/`&` header lines. A 30-level ¹²C + p model takes about
9 s.

Reference result (¹²C + p, converted with `-g`), angle-integrated ¹²C(p,p) against the `.azr`:
1/2⁺ peak at E_c.m. 0.4236 MeV (expected 0.4239), FWHM 32.5 keV (Γ_p 34.2 keV); 5/2⁺ at 1.5993 MeV (1.6013);
the broad 3/2⁻ interfering with the 5/2⁺ shows a shoulder at 1.566 MeV (1.557). Interference moves apparent maxima,
so offsets of this size are expected; larger ones point to a mass or separation-energy error.
