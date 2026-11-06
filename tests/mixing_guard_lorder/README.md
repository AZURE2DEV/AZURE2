# mixing_guard_lorder

Secondary-gamma multipole mixing at every Legendre order, and a warning for channels beyond the
Legendre cap (audit 2026-10-25, checklist t3-24).

**Mixing.** A UPOS segment with a mixed L/L+1 gamma decay uses the radiation parameter
`R_k = (R_k(L,L) + 2 delta R_k(L,L') + delta^2 R_k(L',L')) / (1 + delta^2)` (Brune and deBoer 2020,
Eq. 7) at each even order k. Up to dev 7d4558b `GenMatrixFunc` and `AMatrixFunc` applied the mixed
form only when `R_k(L',L') != 0`, so at an order where that coefficient happens to vanish while
`R_k(L,L)` does not, the interference term and the `1/(1 + delta^2)` were dropped for that order. The
first such case is a 3/2 to 3/2 M1/E2 decay at k = 2: `R_2(2,2) = 0`, `R_2(1,1) = -0.400`,
`R_2(1,2) = +0.775`, so the `a_2` coefficient stayed at its pure-M1 value whatever delta (at
delta = 0.3 it should change sign and shrink to 0.15 of it). The form is now applied unconditionally;
every `R_k` is zero when its coupling is not allowed, so nothing else changes.

**Legendre cap.** `Config::maxLOrder = 20` is the highest Legendre order computed, so a particle
channel with l > 10 contributes terms (l + l' > 20) that were silently dropped from every differential
cross section it entered. `CNuc::Initialize` now prints a WARNING per such channel; angle-integrated
values were and are complete.

## Checks (`check.sh path/to/AZURE2`)

- (a), (b): `tests/secondary_gamma/models/mixing.azr` with the UPOS final state set to J = 3/2 and
  delta = 0.3 and 0: both distributions are `a0 + a2 P2` at all 13 angles, and `a2/a0` at delta = 0.3
  over delta = 0 equals `(R_2(1,1) + 0.6 R_2(1,2)) / (1.09 R_2(1,1)) = -0.1485` to 1e-5;
- (c): `tests/13N` plus a J = 21/2- level whose p + 12C channel has l = 11 logs the warning and
  completes; (d): the unmodified 13N logs none.

dev 7d4558b: (b) gives a ratio of exactly 1, (c) prints nothing. `tests/secondary_gamma` (final state
1/2, where `R_2(2,2)` is not zero) is unchanged.
