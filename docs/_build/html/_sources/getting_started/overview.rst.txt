Overview
========

AZURE2 is designed with two main modes of operation:

1. **Data-driven mode** -- fitting R-matrix calculations to experimental data using
   least-squares minimization (MINUIT2).
2. **Calculation mode** -- computing cross sections, S-factors, or reaction rates
   at user-specified energies and angles, without fitting to data. This mode can
   also be used to extrapolate or interpolate results from a previous fit.

GUI Organization
----------------

The graphical interface is organized into tabs, ordered from left to right in
the sequence that information should be entered:

1. **Particle Pairs** -- define the reaction participants
2. **Levels and Channels** -- define compound nucleus levels and reaction channels
3. **Segments** -- specify experimental data and calculation regions
4. **Experimental Effects** -- configure target integration, convolution, and other corrections
5. **Fitting** -- manage fit parameters and their limits
6. **Calculate** -- select and execute the calculation
7. **Plot** -- visualize results (requires QWT)
8. **MCMC** -- Bayesian parameter inference via Markov Chain Monte Carlo (if enabled)

.. tip::

   Always start by entering information in the **Particle Pairs** tab. AZURE2
   automatically calculates allowed channels and other quantities based on the
   particle pair information. Changes in any tab are automatically propagated to
   all other tabs.

Menus
-----

In addition to the tabs, three menus are available from the menu bar:

- **File** -- open, save, and manage project files
- **Configure** -- set the R-matrix formalism, check files, directories, and runtime options
- **Documentation** -- access in-app help for the current tab

Project Files
-------------

All setup information for an AZURE2 calculation is stored in a single **Input File**.
This file can have any name, but the convention is to use the ``.azr`` extension.
The file is a text-formatted file that is normally created and edited through the GUI.

The recommended project directory structure is::

   my_project/
   ├── my_project.azr       # Input file
   ├── data/                 # Experimental data files
   ├── output/               # Output files from calculations
   └── checks/               # Debug check files

.. tip::

   Use relative paths (relative to the Input file location) for data files and
   output directories. This makes it easy to share projects with collaborators.

Workflow Summary
----------------

A typical AZURE2 workflow follows these steps:

1. **Define particle pairs** -- specify the entrance and exit channels of the reaction
2. **Define levels** -- enter the compound nucleus levels with their energies, spins, and parities
3. **Set channel parameters** -- enter initial values for partial widths or ANCs
4. **Create data segments** -- link experimental data files and specify energy/angle ranges
5. **Configure experimental effects** -- add target integration or convolution corrections if needed
6. **Calculate with data** -- make an initial calculation to check starting parameters
7. **Fit with data** -- perform automated least-squares fitting
8. **Extrapolate** -- compute cross sections at energies not covered by data
9. **Calculate reaction rate** -- integrate the cross section over a Maxwell-Boltzmann distribution
