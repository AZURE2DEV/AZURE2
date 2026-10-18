# spin_observables_reaction — isDiff 9 in a reaction, checked by time reversal

A 4He toy with p + t (pair 1) and n + 3He (pair 2). Nine spin observables are computed for
t(p,n)3He (segments 2-10) and for 3He(n,p)t (segments 12-20) at matching excitation energies and c.m.
angles. `check.sh` verifies the time-reversal identities (beam analyzing power = outgoing polarization
of the inverse reaction, target analyzing power = the other outgoing polarization, transfer
coefficients) to 1e-6 (the engine reaches 1e-9), and that start-up refuses tokens the pairs cannot
carry (bad operator, Cartesian tensor on spin 1/2, rank above 2j, a particle twice, an exit particle
before '>'). The data files are dummies (1.0 +- 0.1); expected/chiSquared.out only guards the numbers.
