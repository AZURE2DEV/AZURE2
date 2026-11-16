# udr_low_energy — a time-of-flight resolution kernel at thermal to keV energies

dev's `tests/udr_low_energy` (dev `a68474f0`), brought here on 2026-11-15 with its repair
(`docs/upgrade/energy_grid.md` section 13).

`model.azr`: the Oliva 2024 17O + n level scheme with the Koehler & Graff 17O(n,α)14C cross section (EXFOR 13534,
`data/koehler1991.dat`, 67 points from 25 meV to 1 MeV; 1/v below a few keV), converted with `AZURE2 --convert`
(no note, no question). `udr_delta.dat`: a 0.2 ns wide delay distribution. The script adds a `detector` record with
a `udr` kernel (7 m flight path) to the data set, three times, and compares each folded cross section with the
bare one:

| | kernel | must hold | measured |
|---|---|---|---|
| (a) | the 0.2 ns delta alone | bare at all 67 points to 1e-4, no zero | 8.6e-6 |
| (b) | a channel of 7.55 bins per decade (±7.6 % of the flight time) | bare on the 1/v region (E ≤ 10 eV) to 3e-3: σ ∝ t is linear in the delay | 6.5e-5 |
| (c) | a 64 ns channel under a 125 ns burst | bare below 10 keV to 1e-3 | 5.0e-5 |

A negative channel width is |n| bins per decade of energy, as on the old `<targetInt>` line (the format writes it
with the unit ns; `HANDOFF.md`, entry of 2026-11-15).

Before the repair the three fail: zeros at 17, 7 and 3 points. A few seconds.
