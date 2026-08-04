# How AZURE2 computes polarization observables

An R-matrix code produces a collision matrix. An experiment with a polarized beam
measures a spin asymmetry. This note explains what sits between the two, why each
step has the form it does, and how the treatment adapts to the spins of the
particles involved — worked through for two real cases.

---

## 1. The object in the middle

Everything hinges on one quantity, the **amplitude matrix**

$$M_{s'\nu'\,s\nu}(\theta),$$

the amplitude to enter with channel spin $s$ and projection $\nu$ and leave with
$s'$, $\nu'$. Why this object and not something else:

- **Below it**, the R-matrix machinery ends at the collision matrix
  $U^J_{s'l'sl}$, which knows about total angular momentum $J$ and orbital
  motion $l$ but is not an amplitude for anything you can point a detector at.
- **Above it**, every observable — cross section, analyzing power, spin
  correlation, tensor moment — is a trace of $M\rho M^{\dagger}$ against some
  operator. No further angular-momentum algebra is needed once $M$ exists.

So $M$ is the natural waist of the calculation: build it once, and the observables
are bookkeeping. This is why the implementation is organized around
`Polarization::AmplitudeMatrix` rather than around $A_y$ specifically.

### The chain, and where it lives

| step | quantity | code |
|---|---|---|
| level energies, reduced widths | $E_\lambda,\ \gamma_{\lambda c}$ | the `.azr` file |
| level matrix inversion | $A_{\lambda\mu}$ | `AMatrixFunc::InvertMatrices` |
| collision matrix | $U^J_{s'l'sl}$ | `AMatrixFunc::CalculateTMatrix` |
| transition matrix | $\delta-U$, with Coulomb phases | the same, as `tmatrix` |
| **amplitude matrix** | $M_{s'\nu' s\nu}(\theta)$ | `Polarization::AmplitudeMatrix` |
| observable | $A_y$ | `AnalyzingPowerAy` |

One design decision is worth stating because it is load-bearing. AZURE2 already
forms

$$\texttt{tmatrix} \;=\; e^{i(\omega_l+\omega_{l'})}\big(\delta_{ss'}\delta_{ll'}-U^J_{s'l'sl}\big),$$

which is exactly the bracket the amplitude needs, Coulomb phases included. The
polarization code consumes that rather than rebuilding $U$, so it is guaranteed to
see the same collision matrix as the cross section — including every
boundary-condition and Brune-parametrization subtlety. Two independent
constructions of $U$ could disagree with neither being obviously wrong.

---

## 2. Channel spin: the basis, and why the code uses it

A reaction channel in R-matrix theory is labelled by the pair of nuclei, their
relative orbital angular momentum $l$, and the **channel spin**

$$\mathbf s \;=\; \mathbf I_1 + \mathbf I_2, \qquad
  |I_1-I_2| \le s \le I_1+I_2 ,$$

the vector sum of the two intrinsic spins. Channels are then grouped by total
angular momentum and parity,

$$\mathbf J = \mathbf l + \mathbf s, \qquad \pi = \pi_1\pi_2(-1)^l .$$

The reason this is the useful basis is that the nuclear Hamiltonian conserves
$J$ and $\pi$, so the R matrix is block-diagonal in $J^\pi$, and within a block
the channels are exactly the allowed $(l,s)$ combinations. Every level in the
`.azr` file carries one reduced width amplitude per such channel. So the
channel-spin basis is not a convenience — it is the basis in which the model
parameters are defined.

The coupling order follows Lane and Thomas, whose formalism AZURE2 implements
(*Rev. Mod. Phys.* **30** (1958) 257, §III.2a): the channel spin is *"formed by
coupling $I_1$ and $I_2$ together"*, with vector-addition coefficients
$(I_1I_2i_1i_2|s\nu)$ transforming *"from the $(I_1i_1, I_2i_2)$ scheme to the
$(I_1I_2, s\nu)$ scheme … as discussed by Condon and Shortley"*. Particle 1 —
which AZURE2 defines as the lighter of the pair — comes first.

---

## 3. Building the amplitude matrix

The bridge from $U$ to $M$ is Seyler's Eq. (4)
(*Nucl. Phys.* **A124** (1969) 253):

$$
M_{s'\nu' s\nu}(\theta) = \frac{\sqrt{\pi}}{k}\Bigg[
  -C(\theta)\,\delta_{ss'}\delta_{\nu\nu'}
  + i\sum_{J l l'} \sqrt{2l+1}\;
    (s\,l\,\nu\,0|J\nu)\,
    (s'\,l'\,\nu'\,\nu{-}\nu'|J\nu)\,
    e^{i(\omega_l+\omega_{l'})}
    \big(\delta_{ss'}\delta_{ll'}-U^J_{s'l'sl}\big)\,
    Y_{l'}^{\nu-\nu'}(\theta,0)
\Bigg]
$$

Every factor is forced by something:

**$C(\theta)\,\delta_{ss'}\delta_{\nu\nu'}$** — the Coulomb amplitude. A pure
Coulomb field is spin-independent, so it cannot change $s$ or $\nu$; and it
exists only when the entrance and exit pairs are the same. In the code,
`AddCoulomb` is called only for elastic scattering.

**$(s\,l\,\nu\,0|J\nu)$** — the entrance coupling, with the orbital projection
fixed at zero. This is the choice of quantization axis: take $\hat z$ along the
beam, and the incoming plane wave has no angular momentum about its own
direction, so $m_l=0$. The total projection is therefore $\nu$ alone.

**$(s'\,l'\,\nu'\,\nu-\nu'|J\nu)$** — the exit coupling. $J$ and its projection
are conserved, so whatever the exit channel spin takes as $\nu'$, the outgoing
orbital motion must carry the remainder $\mu=\nu-\nu'$.

**$Y_{l'}^{\nu-\nu'}(\theta,0)$** — the angular function of that outgoing orbital
motion, evaluated at azimuth zero because we chose the scattering plane.

### Why this needs machinery a cross-section code does not have

Look at the spherical harmonic. For an unpolarized cross section only $\mu=0$
survives the spin averaging, and $Y_l^0 \propto P_l(\cos\theta)$ — plain Legendre
polynomials, which is all AZURE2 had. A spin-flip amplitude has $\nu\neq\nu'$,
hence $\mu\neq0$, hence **associated** Legendre functions.

That is a structural statement, not a missing feature: a code whose angular basis
is $\{P_l\}$ cannot represent the amplitudes whose interference *is* the analyzing
power. Adding $Y_l^m$ (`AngCoeff::SphericalHarmonic`, on GSL's
`gsl_sf_legendre_sphPlm`) was the one genuinely new ingredient.

---

## 4. From the amplitude matrix to an observable

Given a beam described by a spin density matrix $\rho_{\rm in}$,

$$\rho_{\rm out} = M\,\rho_{\rm in}\,M^{\dagger},$$

and any observable is a trace of $\rho_{\rm out}$ against the operator the
apparatus measures. Two examples:

$$\frac{d\sigma}{d\Omega} = \frac{1}{(2I_1+1)(2I_2+1)}\,\mathrm{Tr}\big(MM^{\dagger}\big),
\qquad
A_y = \frac{\mathrm{Tr}\big(M\,\sigma_y\,M^{\dagger}\big)}{\mathrm{Tr}\big(MM^{\dagger}\big)} .$$

The averaging factor in the cross section counts the entrance spin states an
unpolarized beam populates equally. In $A_y$ it cancels between numerator and
denominator — which is the formal reason **a normalization factor cannot affect
an analyzing power**, and why the GUI disables *Vary Norm?* for it.

**Why $\sigma_y$ and nothing else.** $M$ is a matrix in spin space, so it can be
expanded in the identity and the Pauli matrices, with coefficients built from the
only vectors available, $\mathbf k_{\rm in}$ and $\mathbf k_{\rm out}$. Under
parity, $\boldsymbol\sigma$ is a pseudovector while $\mathbf k$ is a vector, so
the coefficient of $\boldsymbol\sigma$ must be a pseudovector — and the only one
available is $\mathbf k_{\rm in}\times\mathbf k_{\rm out}$. Hence

$$\hat{\mathbf n}=\frac{\mathbf k_{\rm in}\times\mathbf k_{\rm out}}{|\mathbf k_{\rm in}\times\mathbf k_{\rm out}|},
\qquad A_x=A_z=0 \ \text{identically}.$$

$A_y$ *is* the vector analyzing power; the Madison convention just fixes the sign
of $\hat{\mathbf n}$.

### The one place the particles' spins enter the recipe

$\sigma_y$ acts on the **beam**. The matrix $M$ is indexed by **channel** spin.
These are the same thing only when the target has no spin. In general the
entrance index must be un-coupled before the Pauli matrix can be applied:

$$
M_{\text{out};\,m_1 m_2} \;=\; \sum_{s}
  \big\langle I_1 m_1\, I_2 m_2 \,\big|\, s,\, m_1{+}m_2 \big\rangle\;
  M_{\text{out};\,s,\,m_1+m_2},
$$

with $m_1$ the projectile projection and $m_2$ the target's. The observable is
then

$$
A_y = \frac{2\sum_{\text{out}}\sum_{m_2}
   \mathrm{Im}\!\big[\,M_{\text{out};+\frac12,m_2}\;M^{*}_{\text{out};-\frac12,m_2}\big]}
  {\sum_{\text{out}}\sum_{m_1 m_2}\big|M_{\text{out};m_1m_2}\big|^{2}} .
$$

The sum over $m_2$ is **incoherent**: nobody prepared or measured the target
spin, so its projections are averaged over, not added as amplitudes. The sum over
exit configurations is incoherent for the same reason. Only the two beam
projections are kept coherent, because that is the pair the polarization
distinguishes.

Two structural readings of that formula, both useful in the laboratory:

- The numerator is an **interference** between the two beam spin states. If they
  scatter identically it vanishes; if they scatter differently but *in phase*,
  the imaginary part vanishes and it is still zero. A non-zero analyzing power
  requires **both** spin dependence and a relative phase — which is why $A_y$ is
  large near resonances, where phases move fast, and why it is such a sharp
  discriminator of interference that a cross section cannot see.
- $A_y$ must vanish at $\theta=0^\circ$ and $180^\circ$, where $\hat{\mathbf n}$ is
  undefined — formally because $Y_l^{\mu}(0)=0$ for every $\mu\neq0$ — and in the
  pure Coulomb limit, where the scattering is spin-independent.

---

## 5. Example A — $^{12}\mathrm{C}(\vec p,p)$: spin-½ on spin-0

The simplest case, and the one where the channel-spin basis needs no unpacking.

**Channel spins.** $I_1=\tfrac12$ (proton), $I_2=0$ ($^{12}$C ground state), so

$$s = \left|\tfrac12-0\right|\ldots\tfrac12+0 = \tfrac12 \quad\text{only}.$$

There is one channel spin, its projections are $\nu=\pm\tfrac12$, and because
$I_2=0$ the Clebsch–Gordan coefficient in §4 is $\langle\tfrac12 m_1\,0\,0|\tfrac12 m_1\rangle=1$.
**The channel spin projection is the proton's own spin projection.** No
decomposition is needed.

**Size of the problem.** Entrance states: 2. Exit states: 2. So $M$ has
$2\times2 = 4$ elements — confirmed by the code, which reports `n=4` amplitudes.

**Structure.** Being a $2\times2$ matrix in the proton spin, $M$ must take the
Wolfenstein form allowed by parity (§4):

$$M = g(\theta)\,\mathbb{1} + h(\theta)\,\boldsymbol\sigma\!\cdot\!\hat{\mathbf n}
    = \begin{pmatrix} g & -ih \\ ih & g \end{pmatrix},$$

so the two **non-flip** elements must be equal, and the two **flip** elements
equal and opposite. The code's diagnostic dump confirms exactly this:

```
nonflip(++) = (-1.9508e+03, -1.9546e+03)
nonflip(--) = (-1.9508e+03, -1.9546e+03)      <- equal
flip(+-)    = (-2.7472e-03,  7.9561e-03)
flip(-+)    = (+2.7472e-03, -7.9561e-03)      <- equal and opposite
```

This is a strong statement about the whole construction: the Clebsch–Gordan
ordering, the spherical-harmonic indices and the phases all have to be right for
that symmetry to come out of a sum over many $(J,l,l')$ pathways.

Substituting the matrix into the trace gives the classical result

$$A_y = \frac{2\,\mathrm{Re}\!\left(g\,h^{*}\right)}{|g|^2+|h|^2}$$

(texts that write $M=g+ih\,\boldsymbol\sigma\!\cdot\!\hat{\mathbf n}$ get
$\mathrm{Im}$ instead; the difference is only where the $i$ is put). Note that
$|h|\ll|g|$ in the numbers above — the spin-flip amplitude is tiny compared with
the Coulomb-dominated non-flip one — yet $A_y$ can still be large, because it is
a *ratio of an interference to a magnitude*, not a ratio of magnitudes.

**Result.** At $E_p = 1.75$ MeV, between the $3/2^-$ resonance at
$E_x = 3.503$ MeV and the $5/2^+$ at 3.545 MeV:

| $\theta_{\rm c.m.}$ | 20° | 40° | 60° | 80° | 100° | 120° | 140° | 160° |
|---|---|---|---|---|---|---|---|---|
| $A_y$ | +0.178 | +0.207 | −0.254 | **−0.992** | −0.054 | +0.793 | +0.909 | +0.401 |

The value at 80° is one of the four points Baumann *et al.* (1992) mark where
$|A_y|$ reaches unity — they quote $(1.750\ \mathrm{MeV},\,80^\circ,\,\text{negative})$.
Reaching $-0.99$ means the two proton spin states are scattering almost perfectly
*out of phase* at that angle: nearly all of the cross section there is spin
asymmetry. That is only possible because two resonances of opposite parity and
different $J$ overlap, giving the interference somewhere to come from.

---

## 6. Example B — $^{15}\mathrm{N}(\vec p,p)$: spin-½ on spin-½

Now the target carries spin, and the machinery of §4 is needed in full.

**Channel spins.** $I_1=\tfrac12$, $I_2=\tfrac12$ ($^{15}$N is $\tfrac12^-$), so

$$s = 0 \ \text{or}\ 1 \qquad\text{— and \emph{never} } \tfrac12 .$$

This is the crucial difference. The channel spin is no longer the proton's spin;
it is the total spin of the proton–nucleus system, and a polarized proton is not
a state of definite channel spin at all.

**Size of the problem.** Entrance states: $s=0$ contributes 1 ($\nu=0$), $s=1$
contributes 3 ($\nu=0,\pm1$) — four in total, as expected for
$2\times2$ spin states. Same on the exit side, so $M$ has $4\times4=16$ elements.
The code reports `n=16`.

**The decomposition, explicitly.** For two spin-½ particles the Condon–Shortley
coupling is

$$
\begin{aligned}
|1,+1\rangle &= |{\uparrow\uparrow}\rangle, &
|1,0\rangle &= \tfrac{1}{\sqrt2}\big(|{\uparrow\downarrow}\rangle+|{\downarrow\uparrow}\rangle\big), \\
|1,-1\rangle &= |{\downarrow\downarrow}\rangle, &
|0,0\rangle &= \tfrac{1}{\sqrt2}\big(|{\uparrow\downarrow}\rangle-|{\downarrow\uparrow}\rangle\big),
\end{aligned}
$$

where the first arrow is the proton. Inverting, a proton with $m_1=+\tfrac12$ on a
target with $m_2=-\tfrac12$ is

$$\big|{\uparrow\downarrow}\big\rangle = \tfrac{1}{\sqrt2}\Big(|1,0\rangle+|0,0\rangle\Big).$$

So the polarized beam probes a **coherent superposition of the singlet and
triplet channels**, and $A_y$ is sensitive to their relative phase. This is worth
dwelling on, because it explains why a cross section can never substitute: in an
unpolarized cross section the entrance projections are *averaged*, which destroys
exactly those cross terms. The relative phase between $s=0$ and $s=1$ is
invisible to every cross section AZURE2 computes and becomes observable only with
a polarized beam.

Working through §4 for this case, the four $(m_1,m_2)$ amplitudes are built from
the channel-spin ones as

$$
\begin{aligned}
M_{\uparrow\uparrow} &= M_{s=1,\nu=+1}, &
M_{\uparrow\downarrow} &= \tfrac{1}{\sqrt2}\big(M_{1,0}+M_{0,0}\big), \\
M_{\downarrow\downarrow} &= M_{s=1,\nu=-1}, &
M_{\downarrow\uparrow} &= \tfrac{1}{\sqrt2}\big(M_{1,0}-M_{0,0}\big),
\end{aligned}
$$

for each exit configuration, and then $A_y$ pairs $m_1=+\tfrac12$ against
$m_1=-\tfrac12$ at fixed $m_2$, summing $m_2$ incoherently.

**Result.** At $E_p = 3.0$ MeV, elastic:

| $\theta_{\rm c.m.}$ | 20° | 40° | 60° | 80° | 100° | 120° | 140° | 160° |
|---|---|---|---|---|---|---|---|---|
| $A_y$ | −0.002 | +0.002 | −0.008 | −0.081 | −0.158 | **−0.199** | −0.175 | −0.091 |

Smaller than the $^{12}$C example, but that comparison is not a like-for-like
one: the $^{12}$C numbers sit deliberately on top of two overlapping resonances,
while 3.0 MeV in $^{16}$O is not a comparably resonant place. There is,
nonetheless, a genuine systematic effect — with a spin-carrying target the
denominator collects $(2I_2+1)$ times as many amplitudes while the numerator only
collects the beam-spin interferences at fixed $m_2$, so a target with spin tends
to dilute $A_y$ relative to an otherwise identical spin-0 case.

---

## 7. What changes with other particles

The recipe above is general in the *target*. What varies is the projectile and
the channel type.

| entrance | how it is handled | status |
|---|---|---|
| **spin-½ projectile, target of any spin** (0, ½, 1, 3/2, …) | §4 in full: decompose the entrance index, trace over $m_2$ | works; the decomposition was checked to be norm-preserving to $10^{-16}$ for $I_2 = 0,\tfrac12,1,\tfrac32,2,\tfrac52$ |
| **any parity, any $J^\pi$** | parity never enters the formulas; it decides only which $(l,s)$ channels exist, which the channel enumeration already handles | works |
| **inelastic and rearrangement**, $(\vec p,p')$, $(\vec p,\alpha)$ | entrance and exit pairs are independent throughout; the Coulomb term is added only when they coincide | works |
| **spin-1 projectile** (deuteron) | a spin-1 beam is described by a $3\times3$ density matrix, so it carries a *vector* moment $iT_{11}$ and three *tensor* moments $T_{2q}$. These need rank-2 spherical operators, not $\sigma_y$ | returns 0 — a deliberate refusal, not a silent gap |
| **capture**, $(\vec p,\gamma)$ | the photon channel needs its own multipole expansion; the reference is Seyler & Weller, *PRC* **20** (1979) 453 | refused |
| **identical particles**, $p+p$ | the amplitude must be symmetrised. The Coulomb side already is — `GetCoulombAmplitude` returns the Mott amplitude $f_C(\theta)+\epsilon f_C(\pi-\theta)$ — but the nuclear side is not | **not handled; see below** |

### The identical-particle case

When the two particles are identical the detector cannot distinguish "projectile
scattered through $\theta$" from "target recoiled at $\pi-\theta$", so the
amplitude must be symmetrised. AZURE2 does this for the cross section by
restricting which $(l,s)$ channels exist and then multiplying the resonant term by
4 and the interference term by 2 — which is algebraically the same as doubling the
nuclear amplitude, since $|2M|^2=4|M|^2$ and
$2\,\mathrm{Re}(C^*\!\cdot\!2M)=2\times 2\,\mathrm{Re}(C^*M)$.

The amplitude matrix consumes the already-symmetrised Coulomb amplitude but never
applies that doubling to the nuclear part. For an identical pair the two are
therefore mismatched by a factor of two and $A_y$ would be wrong — not zero,
which is the more dangerous kind of wrong. It is harmless for $\alpha+\alpha$,
where the spin-0 projectile makes $A_y$ vanish anyway, but polarized $p+p$
elastic scattering is a real measurement and is not currently supported. The
correction looks like a single factor; it is not applied because it is untested,
and the test that would settle it already exists (the angle-independence gate of
the amplitude matrix should fail for an identical pair today and pass once the
factor is right).

---

## 8. Two consequences for setting up a calculation

**A ratio cannot be averaged over a target.** A cross section on a thick target is
an integral over the energy loss; an analyzing power is not, because it is a
ratio. What the experiment returns is the ratio of the polarized and unpolarized
*yields*,

$$\langle A_y\rangle=\frac{\int A_y(E)\,\sigma(E)\,dE}{\int\sigma(E)\,dE},$$

so the depths where the reaction is likely count for more. AZURE2 performs this
cross-section-weighted average automatically, using the same yield integrator,
quadrature and straggling treatment as the cross section.

The consequence is severe for charged-particle elastic scattering: Rutherford
scattering makes $\sigma\propto E^{-2}$ diverge at low energy, exactly where
$A_y\to0$, so a thick target drowns the asymmetry. For the gas target in
`tests/13N` a resonant $A_y\approx0.8$ averages down to $\sim10^{-6}$. This is
physics, and it is why analyzing powers are measured on thin targets — Baumann's
were 85 nm of $^{12}$C, about 3.1 keV of energy loss. **A segment compared against
thin-target data should carry no target integration.**

**Uncertainties should be absolute, not relative.** $A_y$ passes through zero at
several angles and energies, and a point sitting near a zero crossing does not
have a correspondingly small uncertainty. Assigning a fixed percentage of the
value gives those points enormous statistical weight and the fit will chase them.
A roughly constant absolute uncertainty is the physically sensible choice.
