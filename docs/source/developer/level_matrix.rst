Level Matrix Construction and Inversion
========================================

This page describes how AZURE2 constructs and inverts the level matrix
(A-matrix) or channel matrix (R-matrix) to compute cross sections. The
two formalisms are mathematically equivalent but differ in computational
efficiency depending on the number of levels vs. channels.

Class Hierarchy
---------------

The matrix calculation is organized through an inheritance hierarchy:

``GenMatrixFunc`` (abstract base)
   Defines the interface and implements the cross section calculation from
   T-matrix elements.

``AMatrixFunc`` (level matrix formalism, the default)
   Works with the A-matrix, a level-by-level matrix.  Each J-group is solved
   in channel space when it has fewer channels than levels and in level space
   otherwise (see *Channel-Space Solution* below), for every formalism, with
   the analytic gradient.

``RMatrixFunc`` (channel matrix formalism)
   Works with the R-matrix, a channel-by-channel matrix (standard and Brune
   formalisms; no Park matrices and no gradient). It is used only when the
   A-matrix is switched off, and not even then with ignore-externals or Park
   on: ``Config::ResolveMatrixRoute`` then takes the A-matrix with a warning.

``THMMatrixFunc`` (an ``AMatrixFunc``)
   The THM (HOES) cross section, which reads elements of the level matrix
   itself (``include/THMMatrixFunc.h``).

**Key files:** ``src/GenMatrixFunc.cpp``, ``src/AMatrixFunc.cpp``,
``src/RMatrixFunc.cpp``, ``src/ChannelSpace.cpp``, ``src/LevelMatrixSolver.cpp``
(headers of the same names in ``include/``), ``src/EnergyPass.cpp``.

Computation Pipeline
--------------------

The solver objects belong to the scratch of the thread that runs an energy
task (``azure2::TaskScratch``, in the evaluation's ``Results``). For every
energy task of the point pass (:doc:`data_model`), ``SolveAtOneEnergy``
calls, through the ``GenMatrixFunc`` interface on a ``Site``:

1. **ClearMatrices()** -- allocate and zero the matrices for each :math:`J^\pi` group.
2. **FillMatrices(Site*)** -- compute the matrix elements from the R-matrix
   parameters, with the task's first unit that reads a T-matrix.
3. **InvertMatrices()** -- factor the matrix (level route; the channel route
   is solved in ``FillMatrices``).
4. **CalculateTMatrix(Site*)** -- compute the T-matrix (transition matrix)
   from the inverted matrix, once per (entrance, exit) group of the task.
5. **CalculateCrossSection(Site*, Results&)** -- compute each unit's
   observable from the T-matrix into its record.

A-Matrix Formalism
------------------

Matrix Construction: ``AMatrixFunc::FillMatrices()``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The code constructs the **inverse** A-matrix (:math:`A^{-1}`) directly, then
inverts it in the next step. The matrix is indexed by level pairs
:math:`(\lambda, \lambda')` within each :math:`J^\pi` group.

The input energy in the center-of-mass frame is:

.. math::

   E_\text{in} = E_\text{CM} + E_\text{sep} + E_\text{exc}

where :math:`E_\text{sep}` and :math:`E_\text{exc}` are the separation and
excitation energies of the entrance pair.

**Standard R-matrix parameterization:**

.. math::

   A^{-1}_{\lambda\lambda'} = (E_\lambda - E_\text{in})\,\delta_{\lambda\lambda'}
   - \sum_c \gamma_{\lambda c}\,\gamma_{\lambda' c}\,L^0_c

where:

- :math:`E_\lambda` is the resonance energy of level :math:`\lambda` (MeV)
- :math:`\gamma_{\lambda c}` is the reduced width amplitude of level :math:`\lambda` in channel :math:`c`
- :math:`L^0_c = S_c - B_c + iP_c` at the energy of the data point, with
  :math:`S_c` the shift function, :math:`P_c` the penetrability and
  :math:`B_c` the boundary condition of channel :math:`c`
  (``Site::GetLoElement``; Eq. :eq:`conv-ps` of
  :doc:`../theory/conventions`). Photon and β-decay channels have
  :math:`L^0_c = 0`.
- The sum runs over all channels in the :math:`J^\pi` group

**Brune parameterization** (``USE_BRUNE_FORMALISM``):

When the Brune parameterization is enabled, additional terms modify the
channel sum for penetrability-type channels (``radType == 'P'``):

.. math::

   A^{-1}_{\lambda\lambda'} = (E_\lambda - E_\text{in})\,\delta_{\lambda\lambda'}
   - \sum_c \gamma_{\lambda c}\,\gamma_{\lambda' c}
   \left[L^0_c + B_c - S_c(\lambda, \lambda')\right]

where:

- :math:`B_c` is the boundary condition for channel :math:`c`
- For diagonal elements (:math:`\lambda = \lambda'`):
  :math:`S_c(\lambda, \lambda) = S_c^\lambda` (the shift function at level :math:`\lambda`)
- For off-diagonal elements (:math:`\lambda \neq \lambda'`):

  .. math::

     S_c(\lambda, \lambda') = \frac{S_c^\lambda (E_\text{in} - E_{\lambda'})
     - S_c^{\lambda'} (E_\text{in} - E_\lambda)}{E_\lambda - E_{\lambda'}}

This divided difference, its slope :math:`s^1` and its value :math:`s^0` at
:math:`E_\text{in} = 0` (below), and their derivatives by the level energies
come from one function, ``azure2::physics::level_pair_shift``
(``physics/include/azure2/physics/LevelPairShift.h``), as fractions over the
one denominator :math:`E_\lambda - E_{\lambda'}`; ``PairShift``
(``include/LevelState.h``) takes the level values from the level state.

Two levels of one :math:`J^\pi` closer than ``kCoincidentLevelGap_MeV``
(:math:`10^{-6}` MeV) make the divided difference 0/0 (at the gap) or lose
:math:`2\epsilon |S| / \Delta E` to round-off (near it). There the function
takes the Taylor forms about the pair, with the mean slope
:math:`\bar S' = (S'_\lambda + S'_{\lambda'})/2` and mean curvature
:math:`\bar S'' = (S''_\lambda + S''_{\lambda'})/2`:

.. math::

   S_c(\lambda, \lambda') = S_c^\lambda + \bar S'_c\, (E_\text{in} - E_\lambda),
   \qquad s^1 = \bar S'_c, \qquad s^0 = S_c^\lambda - \bar S'_c E_\lambda,

exact to :math:`\Delta E^2 S'''/12`, and in the analytic gradients
:math:`\partial s^1/\partial E_\lambda = \partial s^1/\partial E_{\lambda'} =
\bar S''/2`, :math:`\partial S_c(\lambda,\lambda')/\partial E_\lambda =
(S'_\lambda - S'_{\lambda'})/2 + \bar S''(E_\text{in} - E_\lambda)/2`,
:math:`\partial S_c(\lambda,\lambda')/\partial E_{\lambda'} = \bar S''
(E_\text{in} - E_\lambda)/2`. :math:`S'` and :math:`S''` at the level
energies are computed for such levels in every basis
(``CompoundSystem::CalcShiftFunctions``). Across the gap :math:`\chi^2` steps
by the accuracy of the numerical :math:`S'` (about :math:`10^{-8}` relative);
test ``pyazr_coincident_levels``.

**RMC formalism** (``USE_RMC_FORMALISM``):

For radiative channels (``radType == 'M'`` or ``'E'``), an imaginary term is
subtracted from the diagonal elements (``AMatrixFunc::FillMatrices``):

.. math::

   A^{-1}_{\lambda\lambda} \mathrel{-}= i\,\gamma_{\lambda c}^2

The diagonal element then carries :math:`-i\sum_\gamma
\gamma_{\lambda\gamma}^2`, i.e. :math:`-i` times half the radiative width,
since :math:`P_\gamma = 1` with ``--use-rmc`` (Reich and Moore). The
Reich-Moore capture is computed with constant boundary conditions only:
``--use-rmc`` is given with ``--no-brune``, and the program refuses it
otherwise.

**Implementation detail:** Only active levels are stored in the matrix.
An index map (``level_active_index_``) translates original level indices to
compact, zero-based indices for dense matrix storage. Gamma values and shift
functions are pre-cached in buffers to avoid repeated lookups. Channel
contributions with :math:`|\gamma| < 10^{-12}` are skipped for efficiency.

Matrix Factorization: ``LevelMatrixSolver``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

On the level route (below) :math:`A^{-1}` is LU-factored with partial
pivoting (``LevelMatrixSolver``, reusable storage).  The T-matrix needs only
the bilinear forms :math:`g_c^T A\, g_{c'}`, one triangular solve per channel;
the full inverse is formed only on request (``GetAMatrixElement``).

A level whose widths are all below the cutoff decouples: its row and column of
:math:`A^{-1}` are zero but for the diagonal :math:`E_\lambda - E`, which
vanishes at :math:`E = E_\lambda`.  The level contributes nothing at any
energy, so ``FillLevelMatrix`` puts a unit diagonal there and the matrix stays
regular.  Any other singular matrix leaves the group out at that energy,
with a warning printed once.

Channel-Space Solution (default)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Filling and factoring :math:`A^{-1}` costs :math:`O(n^2 c)` and
:math:`O(n^3)` per J-group and energy for :math:`n` levels and :math:`c`
channels.  In every formalism the level matrix has the form

.. math::

   A^{-1}(E) = K - E\,N - G\,\Lambda(E)\,G^T

with :math:`K, N` real symmetric and independent of the energy, :math:`G`
the :math:`n \times c` reduced widths and :math:`\Lambda` diagonal in the
channels (``include/ChannelSpace.h``):

- standard: :math:`K = \mathrm{diag}(E_\lambda)`, :math:`N = 1`,
  :math:`\Lambda_c = L^0_c = S_c - B_c + iP_c`;
- Brune: :math:`K_{\lambda\mu} = \tilde E_\lambda \delta_{\lambda\mu}
  + \sum_c \gamma_{\lambda c}\gamma_{\mu c}\, s^0_{\lambda\mu c}` and
  :math:`N_{\lambda\mu} = \delta_{\lambda\mu} - \sum_c \gamma_{\lambda c}\gamma_{\mu c}\,
  s^1_{\lambda\mu c}`, where the shift term above is
  :math:`S_c(\lambda,\lambda') = s^0 + E_\text{in} s^1`, and
  :math:`\Lambda_c = S_c + iP_c` (the boundary condition cancels);
- Park: :math:`K_{\lambda\mu} = E_\lambda J_{\lambda\mu} + \sum_c
  \gamma_{\lambda c}\gamma_{\mu c} S_{\lambda c}`, :math:`N = J`.

For the particle channels of the Brune and Park matrices a reference shift
:math:`\bar S_c` (the mean of :math:`S_c(E_\lambda)` over the coupled levels)
is moved from :math:`\Lambda` into :math:`K`, which keeps
:math:`\mathrm{Re}\,\Lambda` small; the result does not depend on it.

The generalized eigenproblem :math:`K v = e N v`, :math:`V^T N V = 1`, gives a
formal (standard) R-matrix with energies :math:`e_\lambda` and widths
:math:`\Gamma = V^T G` (once per parameter set and J-group, in extended
precision: Cholesky reduction and Jacobi rotations), and by the Woodbury
identity

.. math::

   G^T A\, G = (1 - R\Lambda)^{-1} R, \qquad
   R_{cc'} = \sum_\lambda \frac{\Gamma_{\lambda c}\Gamma_{\lambda c'}}{e_\lambda - E_\text{in}},

a :math:`c \times c` solve costing :math:`O(n c^2 + c^3)`
(``AMatrixFunc::SolveChannelSpace``, ``SolveFormal``, ``FColumn``).  For the
standard level matrix :math:`V = 1`.  :math:`N` is :math:`a^T a` for the
eigenvectors :math:`a` of Brune's transformation, positive definite when the
parameters have a real standard equivalent.  When it is not, the pencil is
reduced once to Hessenberg-triangular form, :math:`Q^T K Z = H`,
:math:`Q^T N Z = T`, and at each energy :math:`R = (Z^T G)^T (H - ET)^{-1} (Q^T G)`
is an upper Hessenberg solve in real arithmetic (``SolvePencil``): still no
complex fill and factorization of :math:`A^{-1}`.

Individual elements of :math:`A`, general bilinear forms (channel capture)
and the solves :math:`A g_c` of the gradient come from

.. math::

   A = V \left[D^{-1} + D^{-1}\Gamma\, \Lambda (1 - R\Lambda)^{-1}\, \Gamma^T D^{-1}\right] V^T,
   \qquad D = \mathrm{diag}(e_\lambda - E_\text{in})

(``ApplyA``, ``ChannelBilinear``, ``GetAMatrix``).

**Route of a J-group.**  At each energy each J-group takes the channel route
or the level route (fill and LU).  The level route is taken

- when the J-group has no fewer coupled channels than levels (no gain);
- when no formal representation exists (non-finite :math:`K, N`: equal
  level energies; RMC widths with the Brune or Park matrix);
- near a formal pole, when the level there is too narrow to be deflated
  (below).  Near a pole :math:`R` is dominated by one level and the Woodbury
  form loses digits the level matrix keeps.  With
  :math:`s_\lambda = \sum_c \Gamma_{\lambda c}^2 |\Lambda_c| / |e_\lambda - E|`,
  a level with :math:`s_\lambda > 10^2` or
  :math:`s_\lambda\, \delta e / |e_\lambda - E| > 10^{-11}`
  (:math:`\delta e`, two ulps of the energies involved, is the round-off of
  :math:`e_\lambda - E`; :math:`e_\lambda` is kept to extended precision as a
  sum of two doubles) is *deflated*: left out of :math:`R` and solved with the
  channels in the bordered system

  .. math::

     \begin{pmatrix} D_p & -\Gamma_p \\ -\Lambda\Gamma_p^T & 1 - \Lambda R \end{pmatrix}
     \begin{pmatrix} y_p \\ z \end{pmatrix} = \begin{pmatrix} w_p \\ \Lambda\Gamma_r^T D_r^{-1} w_r \end{pmatrix},
     \qquad y_r = D_r^{-1}(w_r + \Gamma_r z)

  for :math:`y = A_f w` (``SolveBordered``, ``FormalApply``).  The level matrix
  is taken when more than eight levels would be deflated or when a deflated
  level is so narrow (width :math:`w_\lambda = \sum_c \Gamma_{\lambda c}^2|\Lambda_c|`)
  that the round-off of the bordered solve matters,
  :math:`2\epsilon\, w_\lambda / (|e_\lambda - E|^2 + w_\lambda^2) > 10^{-11}`.
  On the pencil route there are no levels to deflate: the round-off of column
  :math:`c` of :math:`R\Lambda`, two ulps of :math:`\max|H - ET|` times
  :math:`|(K - EN)^{-1} G_c|^2 |\Lambda_c|`, must stay below :math:`10^{-11}`
  (and :math:`|R\Lambda|` below :math:`10^3`);
- when the bordered system (:math:`1 - \Lambda R` without deflated levels)
  is ill-conditioned (pivot ratio below :math:`10^{-10}`);
- with ``AZURE_LEVEL_SPACE=1``, or for an ``AMatrixFunc`` switched to it with
  ``SetLevelSpaceOnly(true)`` (callers that want the level matrix itself).
  ``THMMatrixFunc`` is always switched to it: the THM (HOES) amplitude reads
  individual elements of the level matrix over the active levels and its
  adjoint the whole matrix (``LevelMatrixABarAdjoint``), so every THM point is
  solved in level space.

The two routes agree to about :math:`10^{-13}` of the largest bilinear form
(``AZURE_CHANNEL_VERIFY=1`` compares them, and both with an extended-precision
LU of the level matrix, at every channel-route solve).  The formal
representation is cached per thread and J-group and rebuilt when a level's
parameters change (``LevelState::Version``).

**Gradient.**  In the formal representation the adjoint of a bilinear term
:math:`\beta F_{ab}` is :math:`O(n c)`: with :math:`u_c = A_f \Gamma_c` and
:math:`Y = 1 + \Lambda F`,
:math:`\bar e_\lambda = -\mathrm{Re}(\bar\beta u_{a\lambda} u_{b\lambda})` and
:math:`\bar\Gamma_{\lambda d} = \mathrm{Re}(\bar\beta (u_{a\lambda} Y_{db} + u_{b\lambda} Y_{da}))`
(``FormalBilinearAdjoint``).  For the standard level matrix these are the
parameter cotangents.  For Brune and Park a scalar gradient pass sums them over
all points (``GradAccum::fe``/``fgamma``) and resolves them once, through the
adjoint of the eigen-decomposition and of :math:`K, N`
(``ResolveFormalCotangents``).  Where the scalar pass takes the level-matrix
adjoint (pencil route, level route, nearly degenerate formal levels), it
accumulates :math:`\sum_E \mathrm{Re}\,\bar M` and :math:`\sum_E E\,\mathrm{Re}\,\bar M`
and contracts them with :math:`dK, dN` once per pass.  The rows of the
residual Jacobian and of the model gradients take, per free parameter, the
rank-two change of :math:`K - EN` it makes (``LevelAdjoint``,
``ChannelSpace::CurrentLevelDerivatives``), in every formalism.  Park's level
matrix is Brune's off the diagonal; on it :math:`N_{\lambda\lambda} = J_\lambda
= 1 - \sum_c \gamma_{\lambda c}^2 S'_{\lambda c}` and
:math:`K_{\lambda\lambda} = E_\lambda J_\lambda + \sum_c \gamma_{\lambda c}^2 S_{\lambda c}`,
so a width :math:`\gamma_{\lambda c}` changes them by
:math:`-2\gamma_{\lambda c}S'_{\lambda c}` and
:math:`2\gamma_{\lambda c}(S_{\lambda c} - E_\lambda S'_{\lambda c})`, and the
energy by :math:`-\sum_c \gamma_{\lambda c}^2 S''_{\lambda c}` and
:math:`J_\lambda + \sum_c \gamma_{\lambda c}^2 (S'_{\lambda c} - E_\lambda S''_{\lambda c})`
(:math:`S', S''` the energy derivatives at :math:`E_\lambda`, which the levels
hold).  The solves :math:`A g_c` come from the channel solves
of the route (``GetRawChannelSolve``, ``ApplyA``).

In the Park formalism :math:`S'` and :math:`S''` are numerical derivatives of the
shift function (Richardson-extrapolated differences of the exact one: four and
five more Coulomb-function evaluations per open channel; for a charged channel
the energy derivative of the Coulomb functions has no closed form, the
:math:`\eta` derivative being an integral; :math:`S'` is
``ChannelFunc::ShiftDerivative``, the :math:`dS/dE` of every other caller,
with its own branches at and below threshold).  ``CompoundSystem::CalcShiftFunctions``
takes them after the shift functions, for all levels at once and in parallel;
the Coulomb functions are computed at every call (``physics/coulomb``, no
memo keyed by the history of the run), so that they do not depend on the
number of threads; channels of a level that differ only in the channel spin
share one, and the working copies of a fit share them through a memo by exact
argument (``ShiftDerivativeMemo``, in the run's scratch).  It computes
:math:`S''`, which only the gradient needs, only when asked: a chi-squared
alone (``AZURECalc::operator()``, ``Chi2On``, ``ResidualsOnly``, MCMC) leaves it
out, and the next gradient fills it in for the levels that moved
(``LevelState::SecondDerivativesValidAt``).  The gradient's table of
:math:`S'` (``BuildShiftDerivTable``) is read from the level state instead of being
differentiated again.  With the rows in place the analytic-gradient minimizer
takes Park's gradient from the residual Jacobian as for the other formalisms
(``AZURECalc::GradientFromJacobian``, which adds the gradient of the
:math:`J > 0` wall when a level is past it): MIGRAD asks for G2 with every
gradient, so the adjoint pass it took besides is saved.  A project with THM
(HOES) segments keeps the adjoint gradient in every formalism: their rows are
not in the residual Jacobian (``ComputeResidualJacobian`` refuses them), so
``GradientFromJacobian`` declines and ``Gradient`` takes the HOES adjoint.

The scattering lengths (:doc:`../theory/scattering_length`) take the same
routes at each pair's threshold: :math:`M_{cc}` is the channel route's
bilinear form (``GetUBilinear``), deflated where a level lies at the
threshold, and its adjoint (``AMatrixFunc::BilinearAdjoint``) is taken on
either route.  It runs after the scalar pass has resolved its formal
cotangents, so for Brune and Park on the channel route it takes the
level-matrix adjoint directly.  ``tests/channel_space/check.sh`` (d) compares
the lengths, and the gradients with the constraints, between the routes and
with finite differences, including a level 0.5 keV above the threshold.

T-Matrix Calculation: ``AMatrixFunc::CalculateTMatrix()``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The collision matrix is computed from the A-matrix in Lane and Thomas's
level-matrix form (Eq. :eq:`conv-u` of :doc:`../theory/conventions`):

.. math::

   U_{cc'} = \Omega_c\,\Omega_{c'}\Big[\delta_{cc'} + 2i\,P_c^{1/2} P_{c'}^{1/2}
   \sum_{\lambda,\lambda'} \gamma_{\lambda c}\,\gamma_{\lambda' c'}\,A_{\lambda\lambda'}\Big],
   \qquad
   T_{cc'} = e^{2i\omega_c}\,\delta_{cc'} - U_{cc'},

where :math:`P_c` is the penetrability of channel :math:`c` and

.. math::

   \Omega_c = e^{i(\omega_c - \phi_c)}, \qquad
   \omega_c = \sigma_l - \sigma_0 = \sum_{n=1}^{l} \arctan\frac{\eta}{n}, \qquad
   \phi_c = \arctan\frac{F_l(\rho)}{G_l(\rho)} .

:math:`\omega_c` is the Coulomb phase of partial wave :math:`l` relative to
that of :math:`l = 0`, so every amplitude is relative to the common factor
:math:`e^{2i\sigma_0}`. :math:`\phi_c` is the hard-sphere phase of Lane and
Thomas at :math:`\rho = k_c a_c`; the hard-sphere phase *shift* is
:math:`-\phi_c`. Written out,

- **Elastic scattering** (:math:`c = c'`):

  .. math::

     T_{cc} = e^{2i\omega_c} - \Omega_c^2\,\big(1 + 2i P_c\,\gamma_c^T A\,\gamma_c\big)

- **Reactions** (:math:`c \neq c'`):

  .. math::

     T_{cc'} = -\Omega_c\,\Omega_{c'}\; 2i P_c^{1/2} P_{c'}^{1/2}\,\gamma_c^T A\,\gamma_{c'}

with :math:`\gamma_c` the column of reduced widths of channel :math:`c`. In
the code (``AMatrixFunc::CalculateTMatrix``) ``umatrix`` is the bracket without
:math:`\delta_{cc'}`, ``uphase`` is :math:`\Omega_c\Omega_{c'}`, and
``GetExpCoulombPhase`` and ``GetExpHardSpherePhase`` return
:math:`e^{i\omega_c}` and :math:`e^{-i\phi_c}`. Photon and β-decay channels
have :math:`\Omega_c = 1`. With ``--use-rmc`` the capture element is
:math:`-i \cdot 2i P_c^{1/2} P_{c'}^{1/2}\,\gamma_c^T A\,\gamma_{c'}`, without
the phases :math:`\Omega`.

Phase factors, penetrabilities and Coulomb amplitudes are computed at set-up
and held by the ``Site`` (its channel rows: one table per compound energy,
shared by every site at that energy and by its node of the ``NodeTable``).

R-Matrix Formalism
-------------------

The R-matrix formalism (``RMatrixFunc``) works with the channel matrix instead.

**R-matrix construction:**

.. math::

   R_{cc'} = \sum_\lambda \frac{\gamma_{\lambda c}\,\gamma_{\lambda c'}}
   {E_\lambda - E_\text{in}}

in the standard formalism (``RMatrixFunc::FillMatrices``). With ``--use-rmc``
the denominator is
:math:`E_\lambda - E_\text{in} - i\sum_\gamma \gamma_{\lambda\gamma}^2`.

**The [1 - RL] matrix** is then constructed. In the standard formalism

.. math::

   [1 - RL]_{cc'} = \delta_{cc'} - R_{cc'}\,L^0_{c'}

with :math:`L^0_{c'} = S_{c'} - B_{c'} + iP_{c'}` as in the level matrix
(``RMatrixFunc::FillMatrices``). In the Brune formalism the R-matrix built
from the Q-matrix (below) already contains the energy-dependent shift, and

.. math::

   [1 - RL]_{cc'} = \delta_{cc'} - i\,R_{cc'}\,P_{c'}

with :math:`P_{c'} = 0` for channels other than particle channels
(``RMatrixFunc::FillMatrices``). This matrix is inverted using the same GSL
LU decomposition, and then
multiplied by the R-matrix to obtain :math:`[1-RL]^{-1}R`, which enters
the T-matrix calculation.

When the Brune formalism is enabled, an intermediate Q-matrix inversion is
performed at the level stage before building the R-matrix from the result.

**At a pole of R.** :math:`R` is infinite at the energy of every level with a
width -- the formal energy in the standard formalism; with Brune, the observed
energy, where that level's row of the real Q-matrix vanishes -- although
:math:`[1-RL]^{-1}R` is finite there.  When the energy is within 1 eV of such
a level, ``FillMatrices`` solves the group in the level form

.. math::

   [1-RL]^{-1}R = \gamma^T \left(M - \gamma L \gamma^T\right)^{-1} \gamma,
   \qquad [1-RL]^{-1} = 1 + [1-RL]^{-1} R L,

(push-through identity; :math:`M` is the matrix whose inverse gives
:math:`R`, :math:`L` the channels' entries of :math:`[1-RL]`), and
``InvertMatrices`` takes these in place of inverting :math:`[1-RL]`.  The
mesh of a target effect places a node exactly at a level's energy, so this
happens in ordinary runs (``tests/pyazr/rmatrix_pole_test.py``).

Cross Section Calculation
-------------------------

The cross section is computed in ``GenMatrixFunc::CalculateCrossSection()``
from the T-matrix elements. The calculation depends on the data type:

**Angle-integrated cross section:**

T-matrix elements from pathways with the same :math:`(J, l, l')` values are
coherently summed (via the ``TempTMatrix`` structure), enabling interference
between internal and external capture pathways:

.. math::

   \sigma = \frac{1}{100} \sum_{J,l,l'} g\,(2J+1)\,I_{12}\,|T_{J,l,l'}|^2

where :math:`g = \pi/k^2` (fm\ :sup:`2`) is the geometrical factor, with
:math:`k` the entrance wave number, :math:`I_{12} = [(2j_1+1)(2j_2+1)]^{-1}`
is the spin factor of the entrance pair, and the factor :math:`1/100`
converts fm\ :sup:`2` to b. An identical entrance pair multiplies the result
by 2 (Eq. :eq:`conv-sigma` of :doc:`../theory/conventions`).

**Differential cross section (with Coulomb):**

When the entrance and exit channels are the same particle pair, the full
Coulomb + nuclear + interference cross section is computed
(``GenMatrixFunc::CalculateBareCrossSection``). The code does not store the physical Coulomb
amplitude but the dimensionless

.. math::

   C_\alpha(\theta) = \frac{\eta}{2\sqrt\pi}\,\csc^2\frac\theta2\,
   \exp\Big(-2i\eta\ln\sin\frac\theta2\Big)

(``Site::CoulombAmplitudesAt``; :math:`\theta` is the c.m. angle, :math:`\eta`
the Sommerfeld parameter). The point-Coulomb amplitude relative to
:math:`e^{2i\sigma_0}` is :math:`f_C = -(\sqrt\pi/k)\,C_\alpha` (fm). With
:math:`T` as above,

.. math::

   \frac{d\sigma}{d\Omega} = g\,|C_\alpha|^2 + |f_N|^2
   + g\,\frac{1}{\sqrt{\pi}}\,
   \mathrm{Re}\Big[i\,C_\alpha(\theta) \sum_{J s l} (2J+1)\,I_{12}\,
   T_{cc}^*\,P_l(\cos\theta)\Big],
   \qquad g = \frac{\pi}{k^2},

where the sum runs over the elastic channels :math:`c = (l, s)` of each
:math:`J^\pi` and :math:`|f_N|^2` is the Blatt–Biedenharn sum
(Eq. :eq:`conv-ad` of :doc:`../theory/conventions`). The first term is :math:`|f_C|^2`. The last is
:math:`2\,\mathrm{Re}[f_C^* f_N]`: for spinless particles
:math:`f_N = (i/2k)\sum_l (2l+1)\,T_l\,P_l(\cos\theta)`, and substituting
:math:`f_C` and :math:`f_N` gives the term above. In terms of :math:`C_\alpha`
the interference is therefore :math:`+\mathrm{Re}[i\,C_\alpha T^*]`, which
equals :math:`-\mathrm{Re}[i\,C_\alpha^* T]`. ``tests/coulomb_nuclear_sign``
checks the sign and the size of this term against an independent partial-wave
sum (:doc:`../theory/conventions`, Eqs. :eq:`conv-elastic` and
:eq:`conv-coulomb`). For identical particles the Coulomb amplitude is
symmetrized per channel spin and the nuclear and interference terms are
multiplied by 4 and 2 (same page, "Identical particles").

Type Definitions
----------------

``include/Constants.h`` defines the matrix types: ``complex`` is
``std::complex<double>``, ``matrix_c`` a dense complex matrix stored as nested
vectors, ``vector_matrix_c`` one such matrix per :math:`J^\pi` group.
