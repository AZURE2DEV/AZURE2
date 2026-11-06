# light_excitation

The light particle of a pair in an excited state (2026-11-05): an optional 32nd column `e1` (MeV) at the end
of a `<levels>` line, the light particle's excitation energy. It is written only when nonzero, so a 31-column
line means `e1 = 0` and every older file reads, and is written back, byte for byte as before. The light
particle's J^pi of that state goes in `j1`/`pi1`, as always; masses stay ground-state masses.

In the engine `PPair` keeps `e1` (`GetLightExE`) and the heavy partner's `e2` (`GetHeavyExE`) apart, and
`GetExE()` is their sum -- the channel's excitation above the ground-state threshold, which is what every
threshold, Q-value and kinematics use. Two particles are identical when Z, mass, spin, parity, g-factor
**and** excitation are equal, so 12C(2+)+12C(2+) is entered as `e1 = e2 = 4.43891`, `j1 = j2 = 2`.

The models are built at run time by `make_models.py`; `analyze.py` does the numerical comparisons (each
prints a number first, so an empty or crashed analysis is a FAIL).

| check | what |
|---|---|
| (a) | An explicit `0` in column 32 on every line gives byte-identical output files and run log (progress bars aside): alpha + 12C in modes 1 and 3, and `tests/13N` (capture, external capture) in mode 1. |
| (b) | alpha' + 12C(2+): the excitation on the heavy partner (`e2 = E*`, `J2 = 2+`) and on the light one (`e1 = E*`, `J1 = 2+`, `J2 = 0+`) give the same (a,a'), (a',a), a + 12C and a' + 12C(2+) elastic cross sections -- angle-integrated, c.m. and lab differential, extrapolated (mode 3) and at lab-frame data points (mode 1) -- to 1e-10 (they come out bit-identical). Without column 32 the light version differs by O(1). |
| (c) | 12C(2+)+12C(2+) entered as `e1 = e2`: identical bosons. With `includeAllJPi=1` and l <= 4 the 16 hard-sphere J^pi groups hold exactly the 39 channels with l + s even (an independent enumeration), the elastic c.m. cross section is symmetric about 90 deg, and a level listing l = 2, s = 3 is refused. Entered the old way (`e1 = 0`, `e2 = 2E*`) the pair is not identical: all 76 channels, and an asymmetric cross section. |
| (d) | Refused with an ERROR and no output: `e1` on a gamma pair, a negative `e1`, two different `e1` for one pair key (a pair's columns are read from its first line), and a UPOS (secondary-gamma) segment whose exit pair has `e1` (the coefficient is that of the heavy partner's decay); the same UPOS segment with the excitation on the heavy partner runs. |

The symmetry in (c) is checked to 1e-6, not 1e-10: the Coulomb amplitude uses `pi = 3.141592650`
(`include/Constants.h`), which breaks sigma(theta) = sigma(180 - theta) by ~1e-7 at 10/170 deg for every
identical pair -- 1.4e-7 for ground-state 12C+12C, a model this change does not touch -- against 8.3e-8 here.

An engine before 2026-11-05 ignores column 32 silently, so (b) and (c) fail on it (the light version opens
pair 2 at S, the e1 = e2 pair is not identical) and (d) runs instead of refusing.

    ./tests/light_excitation/check.sh path/to/AZURE2      # ~30 s
