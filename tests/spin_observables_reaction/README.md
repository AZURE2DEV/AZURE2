# spin_observables_reaction — the quantity spin_observable in a reaction: time reversal

A 4He toy with p + t (pair 1) and n + 3He (pair 2). The same nine spin observables are computed
for t(p,n)3He (data sets `grid1_2` ... `grid1_10`) and for 3He(n,p)t (`grid2_2` ... `grid2_10`)
on matching grids (equal excitation energy, equal c.m. angle); `grid1` and `grid2` are the cross
sections. The data values are placeholders (1 with error 0.1): the chi2 of `expected/` records the
engine's values, not an agreement.

**Truth.** Time reversal (Ohlsen, Rep. Prog. Phys. 35 (1972) 717): an analysing power of the
forward reaction is the outgoing polarization of the inverse one, and the transfer coefficients
map onto each other,

    forward b:y      = inverse e1:y      forward t:y      = inverse e2:y
    forward e1:y     = inverse b:y       forward e2:y     = inverse t:y
    forward b:y>e1:y = inverse b:y>e1:y  forward b:y>e2:y = inverse t:y>e1:y
    forward t:y>e1:y = inverse b:y>e2:y  forward t:y>e2:y = inverse t:y>e2:y

which nothing in the engine imposes: the two directions are separate data sets with their own
entrance pair, Coulomb functions and pathways. `check.sh` compares them (the engine reaches 1e-9;
the limit is 1e-6, with every observable above 0.05 somewhere).

`check.sh` also checks

- the refusals, each with its message: an operator outside the grammar (`E_SPIN_OPERATOR`, from the
  validator), a particle named twice, an exit particle on the entrance side of `>`, a Cartesian
  tensor on a spin-1/2 beam and a rank above 2j (`E_REACTION_NOT_ALLOWED`, where the pairs' spins
  decide);
- the analytic gradient (`fit --use-gradient`, `AZURE_GRAD_CHECK=exit`) against finite differences
  on this project and on `spin_observables_p3He` and `spin_observables_dd` (worst 1e-7 of the
  largest component; the limit is 1e-5).

**Origin.** dev's `tests/spin_observables_reaction`, written in format 3 by `AZURE2 --convert`
(0 notes, 0 questions).

**Run.** `./tests/spin_observables_reaction/check.sh build/src/AZURE2` (4 s);
`tests/run_tests.sh` runs it and compares `output/chiSquared.out` with `expected/`.
