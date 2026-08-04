# The physics of the analyzing power, and why ¹⁵N broke it

A walk through what a polarization observable actually measures, why an R-matrix
code needs extra machinery to produce one, and why an implementation that worked
perfectly for $^{12}\mathrm{C}(\vec p,p)$ returned exactly zero for
$^{15}\mathrm{N}(\vec p,p)$. Formulas are given with the reason they take the
form they do, rather than quoted.

---

## 1. Why bother with polarization at all

A differential cross section is a sum of squared amplitudes,

$$\frac{d\sigma}{d\Omega} \;=\; \sum \left| \text{amplitudes} \right|^2 .$$

Squaring throws away phase. That is a real loss, because phase is where much of
the physics sits: whether two levels interfere constructively or destructively,
whether a resonance is $3/2^-$ or $5/2^+$, whether a background amplitude has
the sign a direct-capture calculation predicts. A fit to a cross section can
usually trade a phase against a width and land on the same curve.

A polarization observable is built differently. It is a **ratio** whose numerator
is an interference between amplitudes that differ in the projectile's spin
projection, and whose denominator is the cross section:

$$A_y \;=\; \frac{\text{interference between the two beam spin states}}{\text{cross section}} .$$

The normalization cancels. The overall phase cancels. What survives is precisely
the *relative* phase. This is why $A_y$ can swing from $+1$ to $-1$ across a few
tens of keV where the cross section shows nothing but a smooth bump, and why a
single analyzing-power data set can settle a spin-parity assignment that a lot of
cross-section data leaves ambiguous.

---

## 2. With spin, the amplitude becomes a matrix

For two spinless particles one complex function $f(\theta)$ says everything, and
$d\sigma/d\Omega=|f|^2$. Once the particles carry spin, the scattering amplitude
must say what it does to the spins, so it becomes a matrix in spin space:

$$M_{s'\nu'\,s\nu}(\theta),$$

the amplitude to come in with channel spin $s$, projection $\nu$, and leave with
$s'$, $\nu'$. The **channel spin** is the vector sum of the two intrinsic spins,

$$\mathbf{s} \;=\; \mathbf{I}_1 + \mathbf{I}_2 .$$

Every observable is then a trace against the operator the apparatus measures. If
the beam is described by a spin density matrix $\rho_{\rm in}$,

$$\rho_{\rm out} \;=\; M\,\rho_{\rm in}\,M^{\dagger},$$

and

$$\frac{d\sigma}{d\Omega} = \frac{1}{2s+1}\sum_{s'\nu'}\sum_{\nu}\left|M_{s'\nu' s\nu}\right|^2
\qquad\text{(trace with the identity)},$$

$$A_y = \frac{\mathrm{Tr}\!\left(M\,\sigma_y\,M^{\dagger}\right)}{\mathrm{Tr}\!\left(M M^{\dagger}\right)}
\qquad\text{(trace with a Pauli matrix)} .$$

The averaging factor $1/(2s+1)$ in the cross section is there because an
unpolarized beam populates every entrance projection equally; in $A_y$ it cancels
between numerator and denominator, which is the formal reason a normalization
cannot affect an analyzing power.

**Why $\sigma_y$ and not $\sigma_x$ or $\sigma_z$.** Parity conservation. Under a
reflection in the scattering plane, a vector polarization along $\hat x$ or
$\hat z$ would change sign while the cross section does not, so those observables
must vanish identically. Only the component along the normal to the scattering
plane survives:

$$\hat{\mathbf n} \;=\; \frac{\mathbf k_{\rm in}\times\mathbf k_{\rm out}}
                              {|\mathbf k_{\rm in}\times\mathbf k_{\rm out}|}
\qquad\text{(Madison convention)} .$$

So $A_x=A_z=0$ is not an approximation; $A_y$ *is* the vector analyzing power.

---

## 3. Where the amplitude comes from in an R-matrix code

An R-matrix code does not naturally produce $M$. It produces the collision matrix
$U^J_{s'l'sl}$. The bridge is Seyler's Eq. (4)
(*Nucl. Phys.* **A124** (1969) 253):

$$
M_{s'\nu' s\nu}(\theta) = \frac{\sqrt{\pi}}{k}\Bigg[
  -C(\theta)\,\delta_{ss'}\delta_{\nu\nu'}
  + i\sum_{J l l'} \sqrt{2l+1}\;
    \underbrace{(s\,l\,\nu\,0|J\nu)}_{\text{entrance}}\;
    \underbrace{(s'\,l'\,\nu'\,\nu-\nu'|J\nu)}_{\text{exit}}\;
    e^{i(\omega_l+\omega_{l'})}
    \big(\delta_{ss'}\delta_{ll'}-U^J_{s'l'sl}\big)\,
    Y_{l'}^{\nu-\nu'}(\theta,0)
\Bigg]
$$

Each factor is forced:

- **$C(\theta)$** is the Rutherford amplitude. A pure Coulomb field does not
  touch spin, so it is diagonal in $s$ and $\nu$ and exists only for elastic
  scattering.
- **$(s\,l\,\nu\,0|J\nu)$** has the entrance orbital projection fixed at $0$.
  That is the choice of axis: the beam defines $\hat z$, so the incoming wave has
  $\mathbf l\cdot\hat z = 0$, and the total projection is $\nu$ alone.
- **$(s'\,l'\,\nu'\,\nu-\nu'|J\nu)$** conserves that projection: whatever the
  exit channel spin takes as $\nu'$, the exit orbital motion must carry the
  remainder $\mu=\nu-\nu'$.
- **$Y_{l'}^{\nu-\nu'}$** is the angular function of that outgoing orbital
  motion.

### The structural point

Look at that last factor. For an unpolarized cross section only $\mu=0$ ever
appears, and $Y_{l}^{0}\propto P_l(\cos\theta)$ — plain Legendre polynomials.
That is why AZURE2, like most R-matrix codes, had `GetLegendreP(L)` and nothing
else.

Spin flip means $\nu\neq\nu'$, hence $\mu\neq 0$, hence **associated** Legendre
functions. So:

> A code equipped only with $P_l(\cos\theta)$ cannot produce a non-zero analyzing
> power — not because of a missing option, but because the amplitudes that would
> interfere are not representable in its angular basis.

Adding $Y_l^m$ was therefore the one genuinely new ingredient; everything else
already existed.

---

## 4. Why $A_y$ needs two things at once

Write the trace out. With $u_i \equiv M_{i,+1/2}$ and $d_i \equiv M_{i,-1/2}$,
where $i$ runs over exit configurations, and using
$(\sigma_y)_{+-}=-i$, $(\sigma_y)_{-+}=+i$:

$$
A_y = \frac{2\sum_i \mathrm{Im}\!\left[u_i\,d_i^{*}\right]}
           {\sum_i \left(|u_i|^2+|d_i|^2\right)} .
$$

Read this slowly, because it explains the phenomenology:

- if the two beam spin states scatter **identically**, $u_i=d_i$, then
  $u_i d_i^*=|u_i|^2$ is real and $A_y=0$;
- if they scatter differently but **in phase**, $u_i d_i^*$ is still real and
  $A_y=0$ again;
- $A_y\neq 0$ requires **both** a spin dependence *and* a relative phase between
  the two.

Which is why analyzing powers are large near resonances — phases move fast there
— and why they are such sharp discriminators of interference.

### Four places it must vanish

| where | why |
|---|---|
| $\theta=0^\circ,180^\circ$ | $\hat{\mathbf n}$ is undefined; formally $Y_{l}^{\mu}(0)=0$ for all $\mu\neq0$, so no spin-flip amplitude survives |
| pure Coulomb limit | Rutherford scattering is spin-independent, so both projections get the same amplitude *and* the same phase |
| no spin–orbit force | the sum over $J$ collapses and the spin-flip amplitude cancels |
| a single partial wave, no background | one common phase, which cancels in the ratio |

The second is a trap worth naming: checking an implementation only at low energy
shows a beautiful zero that proves nothing at all, because zero is the right
answer there.

And $|A_y|\le1$ follows from Cauchy–Schwarz — the cheapest sanity check
available, though a necessary condition only.

---

## 5. What went wrong with ¹⁵N

Here is the whole bug in one sentence: **$\sigma_y$ acts on the projectile, but
$M$ is stored in the channel-spin basis, where the projectile and target spins
are already coupled together.**

### The spin-0 case, where you get away with it

For $^{12}\mathrm{C}(\vec p,p)$ the target has $I_2=0$, so

$$s = \left|I_1 - I_2\right| \ldots I_1+I_2 = \tfrac12 \quad\text{(only)},$$

and the channel spin projection $\nu$ *is* the proton's own $m_p$. The
channel-spin label and the projectile spin label are the same thing, so
$A_y$ can be read straight off the channel-spin amplitudes. This is what the
first implementation did, and it was right — for a spin-0 target.

### The spin-½ case, where you do not

For $^{15}\mathrm{N}(\vec p,p)$ both particles have spin $\tfrac12$, so

$$s = 0 \ \text{or}\ 1, \qquad \textbf{never } \tfrac12 .$$

The code looked for a channel spin of $\tfrac12$, found none, and returned
exactly $0$ — for an observable that is large and perfectly well defined. The
error was not in the physics of the formula but in the assumption that the
channel spin and the projectile spin are interchangeable.

### The fix, and why it is what it is

Undo the coupling before applying $\sigma_y$:

$$
M_{\text{out};\,m_1 m_2} \;=\; \sum_{s}
  \big\langle I_1 m_1\, I_2 m_2 \,\big|\, s,\, m_1{+}m_2 \big\rangle\;
  M_{\text{out};\, s,\, m_1+m_2},
$$

with $m_1$ the **projectile** projection and $m_2$ the target's. Then apply the
trace of §4 to $M_{\text{out};m_1m_2}$, summing over $m_2$ along with the exit
configurations — because nobody measured the target's spin, so it is traced over
incoherently.

Concretely, for two spin-$\tfrac12$ particles the Condon–Shortley coupling is

$$
\begin{aligned}
|1,+1\rangle &= |{\uparrow\uparrow}\rangle \\
|1,\,0\rangle &= \tfrac{1}{\sqrt2}\big(|{\uparrow\downarrow}\rangle+|{\downarrow\uparrow}\rangle\big) \\
|0,\,0\rangle &= \tfrac{1}{\sqrt2}\big(|{\uparrow\downarrow}\rangle-|{\downarrow\uparrow}\rangle\big) \\
|1,-1\rangle &= |{\downarrow\downarrow}\rangle
\end{aligned}
$$

so a proton with $m_1=+\tfrac12$ against a target with $m_2=-\tfrac12$ is *not* a
state of definite channel spin at all — it is
$\tfrac{1}{\sqrt2}\left(|1,0\rangle+|0,0\rangle\right)$. The polarized beam
therefore probes a **coherent superposition of $s=0$ and $s=1$**, and their
relative phase is exactly what $A_y$ is sensitive to.

That last observation also explains why no existing check could have caught the
bug: in an *unpolarized* cross section the two channel spins add **incoherently**
(you average over $\nu$, which destroys the cross terms). The relative phase
between $s=0$ and $s=1$ is invisible to every cross section AZURE2 has ever
computed. It becomes observable only when the beam is polarized.

### Which ordering, and who says so

Since the coupling order is a convention — $\langle I_1 m_1 I_2 m_2|s\nu\rangle$
or $\langle I_2 m_2 I_1 m_1|s\nu\rangle$, differing by $(-1)^{I_1+I_2-s}$ — and
since nothing in AZURE2's cross-section machinery ever fixes it, it has to be
taken from the formalism the code implements. Lane and Thomas are explicit
(*Rev. Mod. Phys.* **30** (1958) 257, §III.2a): the channel spin *"is formed by
coupling $I_1$ and $I_2$ together: $s=I_1+I_2$"*, with coefficients
*"$(I_1I_2i_1i_2|s\nu)$ … from the $(I_1i_1,I_2i_2)$ scheme to the $(I_1I_2,s\nu)$
scheme … as discussed by Condon and Shortley."* Particle 1 first.

Numerically the choice matters at the $10^{-3}$–$10^{-2}$ level in $A_y$ (about
5–7% where it is largest) — real, but not enough to explain a factor of two.

---

## 6. A ratio cannot be averaged over a target

One more place the ratio structure bites. A cross section on a thick target is an
integral over the energy loss,

$$Y=\int_{E_{\rm back}}^{E_{\rm surf}} \frac{\sigma(E)}{\varepsilon(E)}\,dE .$$

An analyzing power is not. Averaging $A_y(E)$ with equal weight along the target
has no meaning; what the experiment returns is the ratio of the polarized and
unpolarized **yields**,

$$\big\langle A_y\big\rangle=\frac{\displaystyle\int A_y(E)\,\sigma(E)\,dE}{\displaystyle\int\sigma(E)\,dE},$$

so the depths where the reaction is likely count for more.

The consequence is drastic for charged-particle elastic scattering. Rutherford
scattering makes $\sigma\propto E^{-2}$ diverge at low energy, exactly where
$A_y\to0$. A thick target therefore drowns the analyzing power: for the gas
target in `tests/13N`, a resonant $A_y\approx0.8$ averages down to $\sim10^{-6}$.
That is physics, not a bug — and it is why analyzing powers are measured on thin
targets (Baumann's were 85 nm of $^{12}$C, about 3.1 keV of energy loss).

---

## 7. So what works now?

| case | status | why |
|---|---|---|
| **Target of any spin** (0, ½, 1, 3/2, …) | **works** | the decomposition of §5 is general; verified norm-preserving to $10^{-16}$ for $I_2=0,\tfrac12,1,\tfrac32,2,\tfrac52$ |
| **Any parity** | **works** | parity never enters the formulas; it only decides which $(l,s)$ channels exist, which AZURE2 already handles |
| **Any $J^\pi$ of the compound levels** | **works** | the sum over $J$, $l$, $l'$ is unrestricted |
| **Inelastic and rearrangement**, $(\vec p,p')$, $(\vec p,\alpha)$ | **works** | entrance and exit pairs are independent; Coulomb is added only when they coincide |
| **Spin-1 projectile** (deuteron) | **returns 0** | correct refusal, not a bug: a spin-1 beam has $iT_{11}$ and the tensor moments $T_{2q}$, which are different observables needing rank-2 operators |
| **Capture**, $(\vec p,\gamma)$ | **refused** | needs the photon analogue, Seyler & Weller, *PRC* **20** (1979) 453 |
| **Identical particles**, $p+p$, $\alpha+\alpha$ | **not handled — see below** | |

### The identical-particle gap

For identical particles the amplitude must be symmetrised: the detector cannot
tell "beam scattered by $\theta$" from "target recoiled at $\pi-\theta$". AZURE2
handles this in the cross section by restricting which $(l,s)$ channels exist and
then multiplying the resonant term by 4 and the interference term by 2
(`GenMatrixFunc.cpp:285`) — which is exactly equivalent to doubling the nuclear
amplitude, since $|2M|^2=4|M|^2$ and $2\,\mathrm{Re}(C^*\!\cdot\!2M)=2\times2\,\mathrm{Re}(C^*M)$.

The Coulomb side is already correct: `EPoint::GetCoulombAmplitude` returns the
Mott-symmetrised amplitude $f_C(\theta)+\epsilon f_C(\pi-\theta)$
(`EPoint.cpp:1267`), and the polarization code consumes it as-is.

The **nuclear** side is not. `PolarizationFunc.cpp` contains no reference to
`IsIdentical`, so for an identical pair the nuclear amplitude is a factor of two
too small relative to the Coulomb one, and $A_y$ would be wrong — not zero,
which is worse, because it looks plausible.

For $\alpha+\alpha$ this is harmless in practice ($I_1=0$, so $A_y=0$ anyway),
but $p+p$ elastic scattering with a polarized beam is a real measurement and
would be wrong today. The fix looks like one line — double the nuclear pathway
contribution when `aa == ir && IsIdentical()` — but it is **not implemented,
because it is not tested**. The test that would validate it already exists: the
angle-independence gate (the ratio of the Blatt–Biedenharn cross section to the
amplitude-matrix spin sum must be constant in angle) would currently *fail* for
an identical pair, and should pass once the factor is right.

---

## 8. The lesson worth keeping

Every check in place before the ¹⁵N report passed, and the answer was still
identically zero for a whole class of reactions. The checks were good ones — they
verified the amplitudes to nine digits against production code, and fixed the
sign against measurement — but they all lived on the spin-0 side of an assumption
nobody had written down.

The assumption was visible in the source the whole time. The comment said *"for a
spin-1/2 projectile on a spin-0 target"* while the documentation said only
*"spin-1/2 projectile"*. The gap between those two sentences was the bug.
