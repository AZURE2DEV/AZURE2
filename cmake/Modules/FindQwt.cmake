#=============================================================================
# FindQwt.cmake - Modified to use pkg-config more conservatively
#
# This module will first try to find Qwt using pkg-config, but only extract
# the Qwt library itself, not all Qt dependencies (since Qt is handled separately).
#=============================================================================
 
# First, try to find Qwt using pkg-config
find_package(PkgConfig QUIET)
pkg_check_modules(PC_QWT QUIET Qt5Qwt6)
 
if (PC_QWT_FOUND)
  message(STATUS "Found Qwt via pkg-config")

  # pkg-config's include directory is not necessarily the one the headers are
  # actually in. Homebrew's qwt-qt5 reports <prefix>/include while qwt_plot.h
  # sits in <prefix>/include/qwt, so taking the reported path verbatim produced
  # a configure that "found" Qwt and a build that could not open qwt_plot.h.
  # Accept pkg-config's answer only if the header is really there, and
  # otherwise search beneath it with the usual suffixes.
  set(_qwt_pc_include "")
  foreach(_dir IN LISTS PC_QWT_INCLUDE_DIRS)
    if (EXISTS "${_dir}/qwt_plot.h")
      set(_qwt_pc_include "${_dir}")
      break()
    endif()
  endforeach()

  if (_qwt_pc_include)
    set(QWT_INCLUDE_DIR "${_qwt_pc_include}")
  else()
    find_path (QWT_INCLUDE_DIR
      NAMES qwt_plot.h
      HINTS ${PC_QWT_INCLUDE_DIRS} ${QT_INCLUDE_DIR}
      PATH_SUFFIXES qwt qwt-qt5 qwt-qt4 qwt-qt3
    )
    if (QWT_INCLUDE_DIR)
      message(STATUS "pkg-config include dir has no qwt_plot.h; using ${QWT_INCLUDE_DIR}")
    endif()
  endif()
  unset(_qwt_pc_include)
  
  # Extract only the qwt library, filtering out Qt libraries since they're handled separately
  set(QWT_LIBRARIES "")
  set(QWT_LIBRARY "")
  
  # Look for qwt library in the library directories returned by pkg-config
  find_library(QWT_LIBRARY_FOUND 
    NAMES qwt qwt-qt5 libqwt libqwt-qt5
    HINTS ${PC_QWT_LIBRARY_DIRS}
    NO_DEFAULT_PATH
  )
  
  if(QWT_LIBRARY_FOUND)
    set(QWT_LIBRARY ${QWT_LIBRARY_FOUND})
    set(QWT_LIBRARIES ${QWT_LIBRARY_FOUND})
    message(STATUS "Found Qwt library: ${QWT_LIBRARY}")
  else()
    # Fallback: just use -lqwt and let the linker find it
    set(QWT_LIBRARY qwt)
    set(QWT_LIBRARIES qwt)
    message(STATUS "Using fallback Qwt library: qwt")
  endif()
  
  set(QWT_VERSION_STRING ${PC_QWT_VERSION})
  set(_QWT_VERSION_MATCH TRUE) # Set version match to TRUE for pkg-config
else()
  message(STATUS "pkg-config failed, falling back to path search...")

  # Original path search for Qwt headers
  find_path (QWT_INCLUDE_DIR
    NAMES qwt_plot.h
    HINTS ${QT_INCLUDE_DIR} /usr/include/qt5/qwt/ /usr/include/qwt/ /usr/local/include/qwt/ /opt/local/include/qwt/
    PATH_SUFFIXES qwt qwt-qt3 qwt-qt4 qwt-qt5
  )

  # Original path search for Qwt library
  find_library (QWT_LIBRARY
    NAMES qwt qwt-qt3 qwt-qt4 qwt-qt5
  )

  # Version check (original logic)
  set(_VERSION_FILE ${QWT_INCLUDE_DIR}/qwt_global.h)
  if (EXISTS ${_VERSION_FILE})
    file(STRINGS ${_VERSION_FILE} _VERSION_LINE REGEX "define[ ]+QWT_VERSION_STR")
    if (_VERSION_LINE)
      string(REGEX REPLACE ".*define[ ]+QWT_VERSION_STR[ ]+\"(.*)\".*" "\\1" QWT_VERSION_STRING "${_VERSION_LINE}")
      string(REGEX REPLACE "([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\1" QWT_MAJOR_VERSION "${QWT_VERSION_STRING}")
      string(REGEX REPLACE "([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\2" QWT_MINOR_VERSION "${QWT_VERSION_STRING}")
      string(REGEX REPLACE "([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\3" QWT_PATCH_VERSION "${QWT_VERSION_STRING}")
    endif()
  endif()
  set(_VERSION_FILE)
  # check version
  set(_QWT_VERSION_MATCH TRUE)
  if (Qwt_FIND_VERSION AND QWT_VERSION_STRING)
    if (Qwt_FIND_VERSION_EXACT)
      if (NOT Qwt_FIND_VERSION VERSION_EQUAL QWT_VERSION_STRING)
        set(_QWT_VERSION_MATCH FALSE)
      endif()
    else()
      if (QWT_VERSION_STRING VERSION_LESS Qwt_FIND_VERSION)
        set(_QWT_VERSION_MATCH FALSE)
      endif()
    endif()
  endif()

  set(QWT_LIBRARIES ${QWT_LIBRARY})
endif()
 
# handle the QUIETLY and REQUIRED arguments
include ( FindPackageHandleStandardArgs )
find_package_handle_standard_args( Qwt DEFAULT_MSG QWT_LIBRARY QWT_INCLUDE_DIR _QWT_VERSION_MATCH )

# Handle the QUIETLY and REQUIRED arguments
# This will set Qwt_FOUND based on the package name passed
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Qwt DEFAULT_MSG QWT_LIBRARY QWT_INCLUDE_DIR _QWT_VERSION_MATCH)
 
mark_as_advanced(
   QWT_LIBRARY
   QWT_LIBRARIES
   QWT_INCLUDE_DIR
)