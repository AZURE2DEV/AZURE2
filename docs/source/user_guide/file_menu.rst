File Menu
=========

The **File** menu is located in the upper-left corner of AZURE2 and provides
standard options for managing project files.

.. list-table::
   :widths: 30 70

   * - **About AZURE2**
     - Display version and attribution information.
   * - **New**
     - Create a new, empty Input File.
   * - **Open...**
     - Open an existing Input File (``.azr``).
   * - **Open Recent**
     - Quick access to the five most recently opened files. Includes an option
       to clear the recent files list.
   * - **Save**
     - Save the current project to the open Input File.
   * - **Save As...**
     - Save the current project to a new file.
   * - **Quit**
     - Close AZURE2.

The Input File stores all setup information for an AZURE2 calculation: particle
pairs, levels, segments, experimental effects, fitting parameters, and
configuration options.

.. tip::

   The developers recommend that each project be stored in its own directory, with
   the Input file at the top level, and ``data/``, ``output/``, and ``checks/``
   subdirectories for the associated files. All paths can be specified as relative
   to the Input file location.
