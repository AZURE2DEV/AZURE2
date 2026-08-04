# Polarization observables in R-matrix theory

What a polarized beam measures, why it carries information no cross section can,
and how the calculation adapts to the spins of the particles involved. Two worked
examples at the end.

---

## 1. The idea, before any formulas

Fire an unpolarized beam at a target and count what comes out at angle theta. You
measure a cross section, and you learn how *much* scattering happens.

Now spin-polarize the beam, so every projectile enters with its spin pointing the
same way. Count again — and count separately on the left and on the right of the
beam. In general **the two counts differ**. The asymmetry between them is the
analyzing power:

$$
A_y = \frac{N_{\mathrm{left}} - N_{\mathrm{right}}}{N_{\mathrm{left}} + N_{\mathrm{right}}}
\Big/ P_{\mathrm{beam}}
$$

where the division by the beam polarization normalizes to a perfectly polarized
beam. It runs from $-1$ to $+1$.

Why should there be an asymmetry at all? Because the nuclear force depends on
spin. A projectile whose spin points "up" relative to its orbital motion feels a
different potential than one pointing "down" — this is the spin–orbit force, the
same one that splits nuclear shells. Left and right correspond to opposite
relative orientations, so they scatter differently.

**The reason this is worth measuring** is subtler and more valuable. A cross
section is a sum of squared amplitudes. Squaring destroys phase information. An
analyzing power is a *ratio*, whose numerator is an interference between two
amplitudes and whose denominator is the cross section — so the overall size and
the overall phase both cancel, and what survives is the **relative phase**
between amplitudes. Where a cross section shows a smooth bump, an analyzing power
can swing from $+1$ to $-1$ in a few tens of keV. That is why a single
analyzing-power measurement can settle a spin-parity assignment that a great deal
of cross-section data leaves open.

---

## 2. What the calculation must keep track of

With spin, a single complex amplitude is no longer enough. The scattering must
say what it does to the spins as well as where it sends the particle, so the
amplitude becomes a **matrix in spin space**. Everything below is bookkeeping for
that matrix.

### The symbols, and what they actually mean

| Symbol | Name | What it means, plainly | Values it takes |
|---|---|---|---|
| $\theta$ | scattering angle | where the detector sits, measured from the beam | $0$ to $180$ degrees |
| $I_1$ | projectile spin | intrinsic spin of the incoming particle | $1/2$ for a proton |
| $I_2$ | target spin | intrinsic spin of the target nucleus | $0$ for carbon-12, $1/2$ for nitrogen-15 |
| $m_1, m_2$ | spin projections | which way each spin points along the chosen axis | $-I$ to $+I$ in steps of 1 |
| $s$ | **channel spin** | the two intrinsic spins added together as vectors | $\lvert I_1-I_2\rvert$ to $I_1+I_2$ |
| $\nu$ | its projection | which way the combined spin points | $-s$ to $+s$ |
| $l$ | orbital angular momentum | how much the pair is "swinging around" each other on the way in | $0, 1, 2, \ldots$ (s-wave, p-wave, …) |
| $l'$ | the same, on the way out | | |
| $\mu = \nu - \nu'$ | orbital projection, outgoing | the angular momentum the *motion* must absorb if the spins flipped | $-l'$ to $+l'$ |
| $J$ | total angular momentum | orbital plus channel spin; conserved, so it labels resonances | $\lvert l-s\rvert$ to $l+s$ |
| $\pi$ | parity | conserved; decides which $l$ can appear | $+$ or $-$ |
| $U$ | collision matrix | what the R-matrix calculation produces: how much of each entrance channel turns into each exit channel | complex, one per $(J, \text{in}, \text{out})$ |
| $M$ | **amplitude matrix** | the amplitude to come in with one spin state and leave with another | complex, one per (in state, out state) |
| $C(\theta)$ | Coulomb amplitude | pure electrostatic (Rutherford) scattering | complex |
| $\omega_l$ | Coulomb phase | the phase the long-range electric field adds | real |
| $Y_{l'}^{\mu}$ | spherical harmonic | the angular pattern of the outgoing motion | function of $\theta$ |
| $k$ | wave number | momentum of relative motion | $1/\mathrm{fm}$ |
| $\hat{n}$ | normal to the scattering plane | the only direction a polarization can point (see §5) | unit vector |

A prime always means "on the way out". So $M_{s'\nu' s\nu}$ reads: *the amplitude
to enter with combined spin $s$ pointing $\nu$, and leave with combined spin $s'$
pointing $\nu'$.*

### Why we add the spins together first

It would seem simpler to track the projectile and target spins separately. The
reason not to is that **the nuclear force conserves total angular momentum $J$ and
parity**, and $J$ is built from the orbital motion plus the *combined* spin:

$$\mathbf{J} = \mathbf{l} + \mathbf{s}, \qquad \mathbf{s} = \mathbf{I}_1 + \mathbf{I}_2 .$$

Resonances have definite $J$ and parity. So if you organize the calculation by
$(l, s)$, each resonance couples to a small, definite set of channels, and the
whole problem block-diagonalizes. Organized by $(m_1, m_2)$ instead, every
resonance would smear across everything. The channel spin is the basis in which
the physics is simple — and, not coincidentally, the basis in which the model
parameters (level energies and reduced widths) are defined.

This choice has one consequence that dominates everything in §6 and §8: **the
channel spin is not the projectile's spin**, except in the special case where the
target has none.

---

## 3. The master formula

The amplitude matrix follows from the collision matrix as

$$
M_{s'\nu' s\nu}(\theta) = \frac{\sqrt{\pi}}{k}\left[
  -C(\theta)\,\delta_{ss'}\delta_{\nu\nu'}
  + i\sum_{J,l,l'} \sqrt{2l+1}\;
    \langle s\,\nu\,l\,0 \mid J\,\nu \rangle\;
    \langle s'\,\nu'\,l'\,\mu \mid J\,\nu \rangle\;
    e^{i(\omega_l+\omega_{l'})}\,
    \left(\delta_{ss'}\delta_{ll'}-U^J_{s'l' sl}\right)
    Y_{l'}^{\mu}(\theta)
\right]
$$

It looks forbidding, but every piece is forced by something physical:

**The Coulomb term** $C(\theta)$ carries the deltas $\delta_{ss'}\delta_{\nu\nu'}$
because an electric field does not touch spin: whatever spin state goes in comes
out unchanged. It appears only in elastic scattering, since Rutherford scattering
cannot transmute one nucleus into another.

**The first bracket** $\langle s\,\nu\,l\,0 \mid J\,\nu \rangle$ says how the
entrance channel spin and orbital motion combine into $J$. Its orbital projection
is fixed at **zero** — that is not an approximation but a choice of axis. Point
the $z$-axis along the beam; a plane wave carries no angular momentum about its
own direction of travel, so the incoming orbital projection vanishes and the
total projection is $\nu$ alone.

**The second bracket** does the same on the way out, and enforces conservation:
whatever the exit spins take as $\nu'$, the orbital motion must carry the
remainder $\mu = \nu - \nu'$.

**The collision matrix term** $\delta_{ss'}\delta_{ll'} - U$ is the part that
actually scatters: subtracting the identity removes the piece of the wave that
went straight through.

**The spherical harmonic** $Y_{l'}^{\mu}(\theta)$ is the angular pattern of the
outgoing motion.

### The one line that explains why polarization needs new machinery

Look at $\mu = \nu - \nu'$ in that spherical harmonic.

If the spins do not flip, $\nu = \nu'$, so $\mu = 0$, and $Y_l^{0}$ is just a
Legendre polynomial $P_l(\cos\theta)$ — the familiar angular distributions of
ordinary cross sections.

If the spins **do** flip, $\mu \neq 0$, and you need the *associated* Legendre
functions, which describe angular patterns that are not symmetric about the beam
axis in the same way.

So a calculation equipped only with $P_l(\cos\theta)$ cannot produce an analyzing
power at all — not for want of an option, but because the amplitudes whose
interference *is* the analyzing power cannot be written down in that angular
basis. Spin flip is inseparable from sideways angular momentum.

---

## 4. From amplitudes to what is measured

Describe the beam by a spin density matrix $\rho_{\mathrm{in}}$ — a compact way of
saying "this fraction of the beam has its spin here, that fraction there". Then
the outgoing spin state is

$$\rho_{\mathrm{out}} = M\,\rho_{\mathrm{in}}\,M^{\dagger},$$

and every observable is a trace of $\rho_{\mathrm{out}}$ against whatever operator
the apparatus is sensitive to. Two cases:

$$
\frac{d\sigma}{d\Omega} = \frac{\mathrm{Tr}\left(M M^{\dagger}\right)}{(2I_1+1)(2I_2+1)},
\qquad
A_y = \frac{\mathrm{Tr}\left(M\,\sigma_y\,M^{\dagger}\right)}{\mathrm{Tr}\left(M M^{\dagger}\right)} .
$$

The denominator in the cross section counts the spin states an unpolarized beam
populates equally. In $A_y$ that factor cancels top and bottom — which is the
formal reason **an overall normalization cannot change an analyzing power**. It is
already a ratio.

### Why the polarization can only point one way

$M$ is a matrix in spin space, so it can be written as a piece proportional to
the identity plus a piece proportional to the Pauli matrices, with coefficients
built from the only two vectors in the problem, the incoming and outgoing momenta.

Now apply parity. The Pauli matrices form a *pseudovector* (they do not change
sign under reflection), while momenta are ordinary vectors (they do). For the
whole expression to have definite parity, the coefficient multiplying the Pauli
matrices must itself be a pseudovector. Out of two ordinary vectors there is
exactly one pseudovector available: their cross product.

$$\hat{n} = \frac{\mathbf{k}_{\mathrm{in}} \times \mathbf{k}_{\mathrm{out}}}
                 {\lvert \mathbf{k}_{\mathrm{in}} \times \mathbf{k}_{\mathrm{out}}\rvert}$$

So a vector polarization can only point along the normal to the scattering plane.
The components in the plane vanish identically, which is why there is one vector
analyzing power and not three. Writing $\sigma_y$ above is shorthand for
"the Pauli matrix along $\hat n$".

### What the numerator is really doing

Write $u$ for the amplitude when the beam spin is up, $d$ for spin down, at the
same outgoing configuration. The trace works out to

$$
A_y = \frac{2\,\sum \mathrm{Im}\left[\, u\, d^{*} \,\right]}
           {\sum \left( \lvert u \rvert^{2} + \lvert d \rvert^{2} \right)} ,
$$

summed over everything not measured. Read it slowly, because it contains the
whole phenomenology:

- if spin-up and spin-down scatter **identically**, then $u = d$, the product
  $u d^{*}$ is real, and $A_y = 0$;
- if they scatter differently but **in phase**, the product is still real, and
  $A_y = 0$ again;
- a non-zero analyzing power needs **both** a spin dependence **and** a relative
  phase between the two.

That is why analyzing powers are large near resonances — a resonance sweeps its
phase through 180 degrees — and why they are such sharp probes of interference
between overlapping levels.

### Four places it must vanish

| Where | Why |
|---|---|
| exactly forward or backward | there is no scattering plane, so $\hat n$ is undefined; the sideways angular functions vanish there |
| far below any resonance | scattering is pure Coulomb, which is spin-blind, so up and down are identical |
| if there were no spin–orbit force | nothing would distinguish the two orientations |
| a single isolated resonance with no background | one common phase, which cancels in the ratio |

The second is a trap when checking a calculation: at low energy, zero is the
*correct* answer, so agreement there proves nothing.

---

## 5. Example A — protons on carbon-12

Carbon-12 has spin **zero**. This is the case where the bookkeeping of §2
collapses to nothing.

**The channel spin.** With $I_1 = 1/2$ and $I_2 = 0$, the only possible value is

$$s = 1/2 .$$

Its projection $\nu$ is then simply the proton's own spin projection: with a
spinless target there is nothing else to combine with. **Channel spin and beam
spin coincide.**

**Size of the problem.** Two entrance spin states, two exit states, so the
amplitude matrix has $2 \times 2 = 4$ entries.

**Its structure.** From the parity argument of §4, a $2\times2$ amplitude matrix
for a spin-1/2 particle on a spinless target must have the form

$$
M = g(\theta)\,\mathbf{1} + h(\theta)\,\boldsymbol{\sigma}\cdot\hat{n}
  = \begin{pmatrix} g & -ih \\ ih & g \end{pmatrix},
$$

with $g$ the non-flip amplitude and $h$ the spin-flip one. So the two diagonal
entries must be **equal**, and the two off-diagonal entries **equal and
opposite**. Calculating the four amplitudes independently, at 1.75 MeV and 80
degrees, gives

| | leave spin up | leave spin down |
|---|---|---|
| **enter spin up** | $-1951 - 1955\,i$ | $+0.0027 - 0.0080\,i$ |
| **enter spin down** | $-0.0027 + 0.0080\,i$ | $-1951 - 1955\,i$ |

exactly the required symmetry, recovered from a sum over many different $J$, $l$
and $l'$ pathways that had no reason to conspire unless the couplings and phases
are right.

Substituting into the trace gives the classical result

$$A_y = \frac{2\,\mathrm{Re}\left(g\,h^{*}\right)}{\lvert g\rvert^{2}+\lvert h\rvert^{2}} .$$

Notice how small the flip amplitude is — around $10^{-2}$ against $2\times10^{3}$,
because the non-flip amplitude is dominated by Rutherford scattering. And yet:

**Result at 1.75 MeV**, between the $3/2^-$ level at 3.503 MeV excitation and the
$5/2^+$ at 3.545 MeV:

| angle | 20 | 40 | 60 | 80 | 100 | 120 | 140 | 160 |
|---|---|---|---|---|---|---|---|---|
| $A_y$ | +0.18 | +0.21 | −0.25 | **−0.99** | −0.05 | +0.79 | +0.91 | +0.40 |

At 80 degrees the analyzing power reaches $-0.99$: essentially every scattered
proton goes to one side. A spin-flip amplitude four orders of magnitude smaller
than the non-flip one produces a nearly complete asymmetry — because $A_y$
measures an *interference*, not a magnitude, and interference is first order in
the small amplitude while the cross section is second order. Two overlapping
resonances of different $J$ and opposite parity supply the phase difference that
makes it possible.

(This point is one of four where Baumann and co-workers measured $\lvert A_y\rvert$
reaching unity in this system, and they place it at the same energy and angle.)

---

## 6. Example B — protons on nitrogen-15

Nitrogen-15 has spin **1/2**. Now the bookkeeping matters.

**The channel spin.** With $I_1 = I_2 = 1/2$,

$$s = 0 \quad \text{or} \quad s = 1, \qquad \text{but never } 1/2 .$$

Two spin-1/2 particles combine into a singlet and a triplet. **The channel spin is
no longer the proton's spin**, and a polarized proton is not a state of definite
channel spin at all.

**Size of the problem.** Entrance states: one from $s=0$, three from $s=1$, so
four — as it must be, since two spin-1/2 particles have $2\times2$ orientations.
Same on the way out, so the amplitude matrix has $4 \times 4 = 16$ entries,
against 4 for carbon.

**The translation between the two languages.** The singlet and triplet states, in
terms of the individual spins (first arrow the proton, second the nitrogen), are

$$
\lvert 1,+1\rangle = \lvert \uparrow\uparrow \rangle, \qquad
\lvert 1,-1\rangle = \lvert \downarrow\downarrow \rangle,
$$

$$
\lvert 1,0\rangle = \frac{1}{\sqrt{2}}\left( \lvert \uparrow\downarrow \rangle + \lvert \downarrow\uparrow \rangle \right),
\qquad
\lvert 0,0\rangle = \frac{1}{\sqrt{2}}\left( \lvert \uparrow\downarrow \rangle - \lvert \downarrow\uparrow \rangle \right).
$$

Turn that around. A proton with spin up, on a target nucleus with spin down, is

$$
\lvert \uparrow\downarrow \rangle = \frac{1}{\sqrt{2}}\left( \lvert 1,0\rangle + \lvert 0,0\rangle \right).
$$

**This is the physical heart of the general case.** The polarized beam prepares a
*coherent superposition of singlet and triplet channels*, and the analyzing power
is sensitive to the relative phase between them. To compute it you must first
translate the amplitudes out of the channel-spin language and back into
"which way is the proton pointing", which for each outgoing configuration reads

$$
\begin{aligned}
M_{\uparrow\uparrow} &= M_{s=1,\,\nu=+1}, &
M_{\uparrow\downarrow} &= \frac{1}{\sqrt{2}}\left(M_{1,0} + M_{0,0}\right), \\[2pt]
M_{\downarrow\uparrow} &= \frac{1}{\sqrt{2}}\left(M_{1,0} - M_{0,0}\right), &
M_{\downarrow\downarrow} &= M_{s=1,\,\nu=-1}.
\end{aligned}
$$

The two beam orientations are then compared at fixed *target* orientation, and the
target orientations are summed **incoherently** — nobody prepared or measured
them, so they are averaged over rather than added as amplitudes.

**Why a cross section could never substitute.** In an unpolarized cross section
the entrance orientations are averaged, and that average destroys precisely the
cross terms between $s=0$ and $s=1$. Their relative phase is invisible to *every*
cross section one can measure on this system. It becomes observable only with a
polarized beam. This is not a statement about a particular code — it is why the
measurement exists.

**Result at 3.0 MeV**, elastic:

| angle | 20 | 40 | 60 | 80 | 100 | 120 | 140 | 160 |
|---|---|---|---|---|---|---|---|---|
| $A_y$ | −0.00 | +0.00 | −0.01 | −0.08 | −0.16 | **−0.20** | −0.17 | −0.09 |

Smaller than the carbon example — though the comparison is not fair, since the
carbon numbers sit deliberately on top of two overlapping resonances. There is
nonetheless a real systematic effect: with a spin-carrying target the denominator
collects amplitudes from every target orientation, while the numerator only
collects the beam-spin interference at each fixed orientation. Averaging over
something you did not measure dilutes the asymmetry.

---

## 7. How the treatment changes with other particles

| Situation | What happens | Why |
|---|---|---|
| spin-1/2 beam, **target of any spin** (0, 1/2, 1, 3/2, …) | works exactly as in §6 | the translation between channel spin and individual spins is general |
| **any parity, any resonance spin** | no change at all | parity and $J$ only decide which $(l,s)$ channels exist; the spin algebra is untouched |
| **inelastic or rearrangement**, such as $(\vec p, p')$ or $(\vec p, \alpha)$ | works | entrance and exit are independent throughout; only the Coulomb term is restricted to elastic scattering |
| **spin-1 beam** (a polarized deuteron) | needs different observables | a spin-1 particle has a $3\times3$ density matrix, so besides a vector polarization it can be *aligned* — stretched along an axis without pointing. That carries tensor moments, which need rank-2 operators, not a Pauli matrix |
| **capture**, $(\vec p, \gamma)$ | needs a separate treatment | the outgoing photon is described by multipole radiation rather than by orbital angular momentum in a channel |
| **identical particles**, such as $p+p$ | needs symmetrization | you cannot tell "beam scattered by $\theta$" from "target recoiled at $180-\theta$", so the two possibilities must be added as amplitudes before squaring |

---

## 8. Two practical consequences

**A thick target destroys an analyzing power.** A cross section measured on a
thick target is an integral: the beam loses energy as it goes, and every depth
contributes. An analyzing power is a ratio, so it cannot be integrated that way.
What the experiment returns is the ratio of the two integrated yields,

$$
\langle A_y \rangle = \frac{\int A_y(E)\,\sigma(E)\,dE}{\int \sigma(E)\,dE},
$$

weighted by the cross section, so the depths where the reaction is likely count
for more.

For charged-particle scattering this is severe. Rutherford scattering makes the
cross section blow up as the beam slows down, exactly where the analyzing power
goes to zero — so the low-energy tail of the target dominates the weighting and
drowns the asymmetry. A resonant $A_y$ of $0.8$ can average down to $10^{-6}$
through a thick gas target. This is why analyzing powers are measured on thin
foils: Baumann's carbon targets were 85 nm, about 3 keV of energy loss.

**Uncertainties should be absolute, not a percentage.** An analyzing power passes
through zero at many angles and energies. A point sitting near a zero crossing
does not have a correspondingly small uncertainty — the measurement is a
difference of two counts, and its error depends on the counts, not on how close
their difference happens to be to zero. Quoting a fixed percentage of the value
gives those points an enormous and entirely artificial statistical weight, and
any fit will chase them at the expense of everything else.
