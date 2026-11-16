# ec_unbound_final — external capture into an unbound final state (a narrow resonance)

A capture final state above the threshold of a particle pair with a channel in its J^pi group is unbound in that
channel. The engine treats it as a narrow resonance (`ECIntegral::IntegrateUnbound`): its exterior is the outgoing
wave O(k_f r)/O(k_f a), O = G + iF, its strength the level's observed width in that channel, and the N_f factor is
continued as (hbar^2/2 mu a) dS/dE at E_f. Below the final state there is no photon: the capture is 0.

**Models** (`models/`): 12C + p (13N) with a capture pair 3 into the 1/2+ level at E_x = 2.3682 MeV, 424.69 keV above
the proton threshold (S = 1.94351 MeV), the 12C(p,g1) transition of Rolfs and Azuma, NPA 227, 291 (1974); E1 lines
from the 1/2- and 3/2- levels, no internal gamma widths into it. `full` is that model; `pure_u` / `pure_b` have every
entrance-side particle width zero and the final state unbound / at the mirror bound energy S - 0.42469 MeV;
`c10b` ... `c40u` have the final state 10, 20 and 40 keV below (b) and above (u) the threshold with the same reduced
width; `grad` is `full` with pseudo-data (`models/data/pseudo.dat`, the Brune capture of `full`, 5 %) and the final
state's width moved from 34 to 30 keV.

**Truth.**

- (b) Pure direct capture into the unbound state over the same transition into the mirror bound state, with
  E_gamma^3 and g^2 taken out, is |I_O|^2 / I_W^2: computed independently with a damping factor exp(-eps r) and
  Richardson extrapolation in eps (five factors, orders 3 and 4 agree to 2e-5), where the code integrates by parts.
  The four numbers are in `check.sh` (`REF`); agreement within 0.5 % is required (0.18 % measured).
- (c) The cross section over E_gamma^3 g^2 is continuous through the threshold: unbound over bound minus 1 at 10, 20
  and 40 keV, fitted with a parabola, extrapolates to 0 at the threshold (|c0| < 0.005; it is up to 0.12 at 10 keV).
- (d) The formal parametrization (`--no-brune`) gives the Brune capture (1e-7).
- (e) The analytic gradient equals finite differences (`AZURE_GRAD_CHECK`).
- (a) No NaN anywhere, 0 below the photon threshold, and a NOTE in the run's log.

**Origin.** dev's `tests/ec_unbound_final` (2026-11-13). Its models are built there by a script from
`tests/formal_ec_final_state/models/A_orig.azr`; here they are that script's ten models written in format 3 by
`AZURE2 --convert` (0 notes each). This engine's curves equal dev's on all of them (1e-10 relative).

**What is refused** is in `tests/ec_final_state_unbound`: a final state at a threshold, or within about a keV above
a charged pair's (eta > 100).

**Run.** `./tests/ec_unbound_final/check.sh build/src/AZURE2` (30 s).
