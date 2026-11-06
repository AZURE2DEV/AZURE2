# AZURE2 upgrade: documentation style

All documentation written for the upgrade (Markdown in `docs/upgrade/`, test-project READMEs, the Sphinx theory
pages, docstrings of physics functions) follows one register: that of a reaction-theory text written for
experimentalists, in the manner of I. J. Thompson's FRESCO paper (Comput. Phys. Rep. 7 (1988) 167) and of Thompson and
Nunes, *Nuclear Reactions for Astrophysics* (Cambridge, 2009). The model is the manner, not the text: nothing is
copied from these sources.

## 1. Structure

1. Each document starts with one paragraph stating what is computed or measured, from which inputs, and why the
   reader needs it. No preamble, no history unless it changes how the result is used.
2. Sections follow the physics (the system, the channels, the collision matrix, the observables, the apparatus, the
   comparison with data), not the code layout. Code locations are given after the physics, in one line.
3. The defining equations come first, numbered where they are referred to again; prose explains them, it does not
   replace them.

## 2. Symbols, units, conventions

1. Every symbol is defined where it first appears, with its unit: "the channel radius a_c (fm)", "the reduced width
   amplitude γ_λc (MeV^½)".
2. Energies state their frame every time: E_lab, E_c.m., or the excitation energy E_x of the compound nucleus.
   Angles state lab or c.m.
3. Conventions are stated, not assumed: angular-momentum coupling order, phase factors (iˡ), sign of reduced widths,
   the boundary condition B_c, the Coulomb-phase reference, the basis of the level parameters (formal, Brune,
   Park). Where a choice differs from Lane and Thomas (Rev. Mod. Phys. 30 (1958) 257), say how.
4. Units follow nuclear-physics practice: MeV and keV, fm, b and mb, keV b for S factors, cm³ mol⁻¹ s⁻¹ for
   N_A⟨σv⟩. SI only where it is the practice (e.g. densities in g cm⁻³, areal densities in atoms cm⁻²).

## 3. Approximations and validity

1. Each approximation is named, its condition of validity stated, and its size estimated where possible ("the
   thin-target approximation; valid while the energy loss is small compared with the resonance width").
2. Numerical choices carry their tolerance: grid spacing, quadrature order, ℓ_max, convergence criteria, with the
   test that shows the tolerance is met.
3. Distinguish what the data determine from what the model assumes (fixed radii, held widths, background poles).

## 4. Results and tests

1. A test or check is described by what it would detect: "an error in the Coulomb-phase reference changes the
   interference term and fails this test at the 10⁻³ level".
2. Numbers are reported with their uncertainty and the quantity it refers to; agreement is stated as a number
   ("agree to 2 × 10⁻⁹ relative"), not as "good".
3. References are to the primary literature (journal, volume, page, year) or to a file and line in the repository.

## 5. Prose

1. Short declarative sentences. Present tense for what the code does; past tense only for what was measured or run.
2. One idea per paragraph. No rhetorical questions, no emphasis by adjectives ("robust", "powerful", "seamless").
3. Plain words; a technical term is used only after it is defined.
4. Lists for parallel items; prose for reasoning.
