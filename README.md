# AZURE2

[![Build Status](https://travis-ci.com/phScholz/AZURE2.svg?token=Pqu1U2LEwBHgCJpiVM1f&branch=master)](https://travis-ci.com/phScholz/AZURE2) 
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/8150ba411b2445fbbf9cffb9b61909cd)](https://www.codacy.com?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=phScholz/AZURE2&amp;utm_campaign=Badge_Grade)
[![codecov](https://codecov.io/gh/phScholz/AZURE2/branch/master/graph/badge.svg?token=WfVDllGHjP)](https://codecov.io/gh/phScholz/AZURE2)


Thank you for choosing AZURE2 to for all your R-Matrix needs.

This file contains a brief description of how to compile the AZURE2 package.  

## Building a Docker container

Docker provides a way to create a portable ``container" for your software.  This is particularly useful if you would like to avoid installing the software from scratch or need to run your software in a High-Performance Computing environment.

Scripts for building the docker container for AZURE2 and running the container for AZURE2 are in the `scripts` directory.  To build and then run the docker container, 

```
$ cd /path/to/AZURE2
$ bash scripts/build_docker.sh
$ bash scripts/run_gui_mac.sh
```

## Building a Singularity/Apptainer container

Docker containers can be run on your local machine but are not usable in HPC environments.  HPC environments require you to use Singularity (which has recently been renamed Apptainer).  To create a singularity image file (.sif) that you can use on an HPC resource, you can convert your Docker image with apptainer:

```
sudo apptainer build AZURE2.sif docker-daemon://azure2:latest
```

This will create a file, `AZURE2.sif`, that you can copy to your HPC resource.

# Building AZURE2 from source

## Dependencies

**CMake**

AZURE2 is compiled using the CMake package.  Version 2.8 or higher of CMake is (probably) required to build AZURE2.  CMake can be downloaded from http://www.cmake.org. 

**GNU Scientific Library**

Additionally, much of the mathematics in AZURE2 uses the GSL routines.  This library can be obtained from http://www.gnu.org/software/gsl/.  

**ROOT, MINUIT2, and OpenMP**

The minimization routines utilized by AZURE are from the Minuit2 package, distributed as part of the ROOT distribution.  If ROOT is compiled from source, the --enable-minuit2 flag must be set when running the configure script.  The libraries are available as a stand-alone package from http://seal.web.cern.ch/seal/snapshot/work-packages/mathlibs/minuit/release/download.html.  It is important to note that the ROOT library DOES NOT build Minuit2 with OpenMP support by default, while the stand-alone version does.  On a multi-core machine the fit process will be much slower without OpenMP.  To build ROOT with OpenMP support for Minuit2, set the enviromnet variables USE_PARALLEL_MINUIT2 and USE_OPENMP prior to building.  

**Readline Development Libraries**

If not already available on your system, you can obtain the readline development libraries via your package manager. For instance, on ubuntu you can run:

```bash
sudo apt-get install libreadline-dev
```

**Qt4.X**

The final compile time dependency is Qt.  Qt is a cross-platform interface API, and required to compile the graphical setup program.  While AZURE2 can be compiled without the graphical setup program, this is not recommended. Qt is available from http://qt.nokia.com. Qt 4.4 or greater is required to build AZURE2.

Qt4 can also be obtained via package managers of some Linux distributions. For instance, on ubuntu you can simply run:

``` bash
sudo apt-get install qt4-dev-tools
```

## Build

### Basic Building

The following steps should be performed to build the AZURE2 package:

1.  Create a subdirectory of the AZURE2 root named build [mkdir build], and change to that directory [cd build].

2.  Run CMake to generate Makefiles from AZURE2 root [cmake ..].  The reference to the parent directory tells CMake where the root of the source tree is located.  

3.  Build the package [make && make install].  The resulting binary will be created in the current directory.

Alternatively, run the *build.sh* script.

### Options

Compiling options (i.e. switching compilers, etc.) are available with flags to CMake.  See the CMake documentation for more details.

A few options are available to the user when building AZURE2.  The can be passed to cmake using the -D[OPTION]=[VALUE] syntax.  CMake also provides a utility to switch these options ON/OFF.   After an initial configuration of the build directory (step 2 above), the command [ccmake ..] can be run to view and edit the configured options. 

These are:

BUILD_GUI [ON/OFF] - Toggles whether the graphical setup utility is built and linked to AZURE2.  This is ON by default, and it is not recommended to turn it OFF.

USE_QWT [ON/OFF] - Toggles whether the built-in plotting tab is added to AZURE2.  This is OFF by default. The plotting features are recommended but require the additional QWT libraries to be installed on the build system.  

USE_STAT [ON/OFF] - Toggles whether the stat() function should be used to test for properly set directories at start.  This only applied when running AZURE2 in console (--no-gui) mode, and is ON by default.  This function has been seen to not work properly for Windows, and it is recommended to select OFF if building on/for a Windows system.

MINUIT_PATH [dir] - If Minuit2 is in a non-standard path, this will add the directory to the search path.

GSL_PATH [dir] - If GSL is in a non-standard path, this will add the directory to the search path.

After changing options, execute [make clean] and then step 3 above to build/rebuild AZURE2.  



**Qt5**

Install libqwt-qt5-dev libqt5svg5-dev

Probably only some of the below are required...

Install qt5-default qtscript5-dev
