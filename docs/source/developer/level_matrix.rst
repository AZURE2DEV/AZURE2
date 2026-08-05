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

``AMatrixFunc`` (level matrix formalism)
   Works with the A-matrix, a level-by-level matrix. Preferred when the
   number of channels exceeds the number of levels.

``RMatrixFunc`` (channel matrix formalism)
   Works with the R-matrix, a channel-by-channel matrix. Preferred when the
   number of levels exceeds the number of channels.

**Key files:**

- ``include/GenMatrixFunc.h``, ``src/GenMatrixFunc.cpp``
- ``include/AMatrixFunc.h``, ``src/AMatrixFunc.cpp``
- ``include/RMatrixFunc.h``, ``src/RMatrixFunc.cpp``
- ``include/MatrixInv.h``, ``src/MatrixInv.cpp``

Computation Pipeline
--------------------

For every data point, the following sequence is executed:

1. **ClearMatrices()** -- allocate and zero the matrices for each :math:`J^\pi` group.
2. **FillMatrices(EPoint*)** -- compute the matrix elements from the R-matrix parameters.
3. **InvertMatrices()** -- invert the matrix using LU decomposition.
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

Matrix Inversion: ``MatrixInv``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The matrix inversion is performed by the ``MatrixInv`` class using **GSL
complex LU decomposition**:

1. The :math:`A^{-1}` matrix (stored as ``std::vector<std::vector<complex>>``)
   is copied into a ``gsl_matrix_complex``.
2. ``gsl_linalg_complex_LU_decomp()`` performs the LU factorization with
   partial pivoting.
3. ``gsl_linalg_complex_LU_invert()`` computes the inverse from the LU factors.
4. The result is copied back into a ``matrix_c`` (``std::vector<std::vector<complex>>``).

This is done independently for each :math:`J^\pi` group. The inversion runs
**sequentially** (not parallelized) to avoid race conditions on the shared
matrix data structures. Move semantics are used to avoid copying the result
matrices.

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
