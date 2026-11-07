Graphical Interface
===================

This version of AZURE2 has no graphical setup utility. A project is an ``.azr``
file; it is written with ``pyazr`` or a text editor, and run with the
``AZURE2`` command-line program or from Python. The released versions (the
``qt5`` and ``dev`` lines of the repository) keep the Qt setup utility, and a
browser interface served by ``pyazr`` is planned to replace it.

How to do what the setup utility did
------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Task
     - Without the setup utility
   * - Create or edit a project
     - ``pyazr`` reads and writes every block of the ``.azr`` file
       (:doc:`pyazr`); the format is in :doc:`../reference/data_formats`.
       The pages of this guide describe each block through the fields the
       setup utility showed; the field names are those of the file.
   * - Calculate, fit, extrapolate, compute a reaction rate
     - ``AZURE2 project.azr`` and the menu it prints, or the flags of
       :doc:`../reference/command_line`; from Python, :doc:`pyazr`.
   * - Runtime options (Brune, Park, relativistic kinematics, THM, ...)
     - Command-line flags (:doc:`../reference/command_line`), or the
       keyword arguments of ``pyazr.azure2`` (``use_brune``, ``use_park``,
       ...). The setup utility kept them in its own settings, not in the
       ``.azr`` file.
   * - Build the Coulomb-function table
     - ``AZURE2 --build-coulomb-table`` (:doc:`configure_menu`, "Coulomb
       functions from a table").
   * - Plot results
     - The files in ``output/`` (:doc:`../reference/output_files`) with any
       plotting program, or the arrays ``pyazr`` returns.
   * - MCMC
     - Option 6 of the command-line menu, with the ``<mcmc>`` block of the
       ``.azr`` file (:doc:`mcmc`), or ``pyazr`` with an external sampler.

The planned browser interface
-----------------------------

The browser interface is started by ``pyazr`` and calls the same engine
through ``pyazr``; it holds no physics of its own. Its tabs follow the
project (nucleus, data, fit, results, physics functions, THM). Before it
replaces the setup utility it is checked, item by item, against the list of
the setup utility's capabilities kept in the repository
(``docs/upgrade/qt_parity.md``).
