Installation
============

AZURE2 requires several external packages. The program is built using the CMake
toolkit and the graphical interface is written using the Qt framework.

Dependencies
------------

**Required:**

- **CMake** 2.8 or later -- build system (https://cmake.org)
- **GNU Scientific Library (GSL)** -- mathematical functions (https://www.gnu.org/software/gsl/)
- **Minuit2** -- minimization routines, distributed within ROOT (https://root.cern.ch)
- **OpenMP** -- parallel processing support

**Optional:**

- **Qt5** -- graphical user interface framework
- **QWT** -- plotting library for the built-in Plot tab (``libqwt-qt5-dev``, ``libqt5svg5-dev``)
- **Readline** -- enhanced command-line input

Building from Source
--------------------

1. Unpack the AZURE2 archive in your directory of choice:

   .. code-block:: bash

      tar zxvf azure2_v1.tar.gz

2. Create a build directory:

   .. code-block:: bash

      cd AZURE2
      mkdir build
      cd build

3. Run CMake to generate the Makefile:

   .. code-block:: bash

      cmake ..

4. Build and install:

   .. code-block:: bash

      make && make install

CMake Options
^^^^^^^^^^^^^

The following options can be passed to CMake using ``-D[OPTION]=[VALUE]``:

.. list-table::
   :header-rows: 1
   :widths: 30 10 60

   * - Option
     - Default
     - Description
   * - ``BUILD_GUI``
     - ON
     - Build the graphical setup utility (recommended)
   * - ``USE_QWT``
     - OFF
     - Enable built-in plotting tab (requires QWT libraries)
   * - ``USE_STAT``
     - ON
     - Use ``stat()`` function for directory checking (disable for Windows)
   * - ``USE_READLINE``
     - ON
     - Enable readline library for CLI input
   * - ``BUILD_LIBRARY``
     - OFF
     - Build as a library instead of an executable
   * - ``CODE_COVERAGE``
     - ON
     - Enable coverage reporting
   * - ``MINUIT_PATH``
     - (auto)
     - Custom path for Minuit2 libraries
   * - ``GSL_PATH``
     - (auto)
     - Custom path for GSL libraries

Example with QWT plotting enabled:

.. code-block:: bash

   cmake .. -DUSE_QWT=ON

Docker
------

AZURE2 supports Docker containerization for easy deployment:

.. code-block:: bash

   # Build Docker container
   source scripts/build.sh

   # Run GUI container
   source scripts/run_gui.sh

For HPC environments, the Docker image can be converted to Singularity/Apptainer:

.. code-block:: bash

   sudo apptainer build AZURE2.sif docker-daemon://azure2:latest
