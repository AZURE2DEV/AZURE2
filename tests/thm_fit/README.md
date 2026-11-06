# thm_fit — a fit through the THM experiment features

`check.sh` fits (menu 2, MIGRAD) one level energy of `tests/18O_p_a_thm`,
unfolded and started 10 keV high, with five `<thm>` experiment lines: an
experiment alone (profiled norm), `distortion=coulomb`, `vertexModel=dw`,
`lineshape=on` and a `ps=` window. Each fit must end, move the energy, not end
above its starting chi-squared, and write a chi-squared that a calculation
from its `param.sav` reproduces (rel. 1e-6).

The calculations of these features are pinned by the other `tests/thm_*`
checks; before this one, no test ran a fit through them (review item P2-35).
