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
   Works with the R-matrix, a channel-by-channel matrix: the older route
   without the Brune and Park matrices, caching or gradient, used only when the
   A-matrix is switched off.

**Key files:**

- ``include/GenMatrixFunc.h``, ``src/GenMatrixFunc.cpp``
- ``include/AMatrixFunc.h``, ``src/AMatrixFunc.cpp``
- ``include/RMatrixFunc.h``, ``src/RMatrixFunc.cpp``
- ``include/ChannelSpace.h``, ``src/ChannelSpace.cpp``
- ``include/LevelMatrixSolver.h``, ``src/LevelMatrixSolver.cpp``

Computation Pipeline
--------------------

For every data point, the following sequence is executed:

1. **ClearMatrices()** -- allocate and zero the matrices for each :math:`J^\pi` group.
2. **FillMatrices(EPoint*)** -- compute the matrix elements from the R-matrix parameters.
3. **InvertMatrices()** -- factor the matrix (level route; the channel route
   is solved in ``FillMatrices``).
4. **CalculateTMatrix(EPoint*)** -- compute the T-matrix (transition matrix) from the inverted matrix.
5. **CalculateCrossSection(EPoint*)** -- compute the cross section from the T-matrix.

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
   - \sum_c \gamma_{\lambda c}\,\gamma_{\lambda' c}\,L_c

where:

- :math:`E_\lambda` is the resonance energy of level :math:`\lambda`
- :math:`\gamma_{\lambda c}` is the reduced width amplitude of level :math:`\lambda` in channel :math:`c`
- :math:`L_c = L_0(c)` is the shift function (Blatt-Biedenharn) evaluated at the data point energy
- The sum runs over all channels in the :math:`J^\pi` group

**Brune parameterization** (``USE_BRUNE_FORMALISM``):

When the Brune parameterization is enabled, additional terms modify the
channel sum for penetrability-type channels (``radType == 'P'``):

.. math::

   A^{-1}_{\lambda\lambda'} = (E_\lambda - E_\text{in})\,\delta_{\lambda\lambda'}
   - \sum_c \gamma_{\lambda c}\,\gamma_{\lambda' c}
   \left[L_c + B_c - S_c(\lambda, \lambda')\right]

where:

- :math:`B_c` is the boundary condition for channel :math:`c`
- For diagonal elements (:math:`\lambda = \lambda'`):
  :math:`S_c(\lambda, \lambda) = S_c^\lambda` (the shift function at level :math:`\lambda`)
- For off-diagonal elements (:math:`\lambda \neq \lambda'`):

  .. math::

     S_c(\lambda, \lambda') = \frac{S_c^\lambda (E_\text{in} - E_{\lambda'})
     - S_c^{\lambda'} (E_\text{in} - E_\lambda)}{E_\lambda - E_{\lambda'}}

**RMC formalism** (``USE_RMC_FORMALISM``):

For radiative channels (``radType == 'M'`` or ``'E'``), an imaginary term is
added to the diagonal elements:

.. math::

   A^{-1}_{\lambda\lambda} \mathrel{+}= i\,\gamma_{\lambda c}^2

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
regular.  Until 2026-10-31 the factorization failed instead and the whole
J-group was left out at that energy (a zero-width level placed on a data
energy changed :math:`\chi^2`); the R-matrix route (``RMatrixFunc``) divided
0 by 0 or failed its GSL inversion.  Any other singular matrix still leaves the
group out at that energy, with a warning printed once.

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
  :math:`\Lambda_c = L_c = S_c - B_c + iP_c`;
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
or the level route (fill and LU, as before).  The level route is taken

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

The two routes agree to about :math:`10^{-13}` of the largest bilinear form
(``AZURE_CHANNEL_VERIFY=1`` compares them, and both with an extended-precision
LU of the level matrix, at every channel-route solve).  The formal
representation is cached per thread and J-group and rebuilt when a level's
parameters change (``ALevel::Stamp``).

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
:math:`\eta` derivative being an integral).  ``CNuc::CalcShiftFunctions``
takes them after the shift functions, for all levels at once and in parallel,
without the Coulomb-function memo (``CoulFuncCache``: per thread and matching
energies to :math:`10^{-12}` MeV, it made :math:`S'` depend on the run's history
by up to :math:`10^{-11}`), so that they do not depend on the number of threads;
channels of a level that differ only in the channel spin share one, and the
working copies of a fit share them through a memo by exact argument.  It computes
:math:`S''`, which only the gradient needs, only when asked: a chi-squared
alone (``AZURECalc::operator()``, ``Chi2On``, ``ResidualsOnly``, MCMC) leaves it
out, and the next gradient fills it in for the levels that moved
(``ALevel::ShiftSecondDerivativesValidAt``).  The gradient's table of
:math:`S'` (``BuildShiftDerivTable``) is read from the levels instead of being
differentiated again.  With the rows in place the analytic-gradient minimizer
takes Park's gradient from the residual Jacobian as for the other formalisms
(``AZURECalc::GradientFromJacobian``, which adds the gradient of the
:math:`J > 0` wall when a level is past it): MIGRAD asks for G2 with every
gradient, so the adjoint pass it took besides is saved.  On the 19-segment
project of :doc:`edata_and_target_effects` in Park mode (one thread) an
evaluation at new parameters went from 1.55 s to 1.09 s, the residual Jacobian
from 5.9 s to 2.2 s and 20 MIGRAD calls from 53 s to 26 s (two threads: an
evaluation 1.31 s to 0.64 s, MIGRAD 36 s to 15 s; level route,
``AZURE_LEVEL_SPACE=1``, one thread: MIGRAD 83 s to 49 s).

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

The collision matrix element :math:`U` is computed from the A-matrix:

.. math::

   U_{cc'} = 2i\,\sqrt{P_c}\,\sqrt{P_{c'}} \sum_{\lambda,\lambda'}
   \gamma_{\lambda c}\,\gamma_{\lambda' c'}\,A_{\lambda\lambda'}

where :math:`P_c` is the penetrability of channel :math:`c`.

The T-matrix is then:

- **Elastic scattering** (:math:`c = c'`):

  .. math::

     T_{cc} = e^{2i\omega_c} - e^{i\omega_c}\,e^{i\phi_c}\,(1 + U_{cc})

- **Transfer reactions** (:math:`c \neq c'`):

  .. math::

     T_{cc'} = -e^{i\omega_c}\,e^{i\phi_{c'}}\,U_{cc'}

where :math:`\omega_c = \sigma_l` is the Coulomb phase shift and
:math:`\phi_c` is the hard-sphere phase shift.

Phase factors, penetrabilities, and Coulomb amplitudes are pre-computed and
cached on the ``EPoint`` object.

R-Matrix Formalism
-------------------

The R-matrix formalism (``RMatrixFunc``) works with the channel matrix instead.

**R-matrix construction:**

.. math::

   R_{cc'} = \sum_\lambda \frac{\gamma_{\lambda c}\,\gamma_{\lambda c'}}
   {E_\lambda - E_\text{in}}

**The [1 - RL] matrix** is then constructed:

.. math::

   [1 - RL]_{cc'} = \delta_{cc'} - i\,R_{cc'}\,P_{c'}

This matrix is inverted using the same GSL LU decomposition, and then
multiplied by the R-matrix to obtain :math:`[1-RL]^{-1}R`, which enters
the T-matrix calculation.

When the Brune formalism is enabled, an intermediate Q-matrix inversion is
performed at the level stage before building the R-matrix from the result.

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

where :math:`g` is the geometrical factor and :math:`I_{12}` is the channel
spin statistical factor.

**Differential cross section (with Coulomb):**

When the entrance and exit channels are the same particle pair, the full
Coulomb + nuclear + interference cross section is computed:

.. math::

   \frac{d\sigma}{d\Omega} = |f_C|^2 + |f_N|^2
   + \frac{i}{\sqrt{\pi}}\,g\,f_C^* \cdot T \cdot P_l(\cos\theta) + \text{c.c.}

Type Definitions
----------------

The matrix types used throughout the code are defined in ``Constants.h``:

.. code-block:: cpp

   typedef std::complex<double> complex;
   typedef std::vector<std::vector<std::complex<double>>> matrix_c;
   typedef std::vector<std::vector<std::vector<std::complex<double>>>> vector_matrix_c;

- ``matrix_c`` -- a 2D complex matrix (dense storage as nested vectors)
- ``vector_matrix_c`` -- a collection of matrices, one per :math:`J^\pi` group
