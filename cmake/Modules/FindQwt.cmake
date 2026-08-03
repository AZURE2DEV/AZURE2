#=============================================================================
# FindQwt.cmake - Modified to use pkg-config more conservatively
#
# This module will first try to find Qwt using pkg-config, but only extract
# the Qwt library itself, not all Qt dependencies (since Qt is handled separately).
#=============================================================================
 
# First, try to find Qwt using pkg-config
find_package(PkgConfig QUIET)
pkg_check_modules(PC_QWT QUIET Qt5Qwt6)
 
# --- headers ---------------------------------------------------------------
# One search for the header, whatever pkg-config did or did not say. Qwt's
# layout varies more than pkg-config admits: the headers sit directly in
# include/ under conda, in include/qwt on most Linux distributions, and inside
# lib/qwt.framework/Headers when Homebrew builds it as a macOS framework.
# Taking pkg-config's include directory on faith produced a configure that
# reported success and a build that could not open qwt_plot.h.
#
# find_path is a no-op when QWT_INCLUDE_DIR is already set, so a caller can
# always override this with -DQWT_INCLUDE_DIR=...
find_path (QWT_INCLUDE_DIR
  NAMES qwt_plot.h
  HINTS
    ${PC_QWT_INCLUDE_DIRS}
    ${PC_QWT_PREFIX}
    ${QT_INCLUDE_DIR}
    /usr/include /usr/local/include /opt/homebrew/include /opt/local/include
    /usr/local/opt/qwt-qt5 /opt/homebrew/opt/qwt-qt5
  PATH_SUFFIXES
    qwt qwt-qt3 qwt-qt4 qwt-qt5
    include include/qwt include/qwt-qt5
    lib/qwt.framework/Headers qwt.framework/Headers
)

# --- library ---------------------------------------------------------------
if (PC_QWT_FOUND)
  message(STATUS "Found Qwt via pkg-config")

  # Only the Qwt library itself; the Qt libraries pkg-config also reports are
  # handled separately by find_package(Qt5).
  find_library (QWT_LIBRARY_FOUND
    NAMES qwt qwt-qt5 libqwt libqwt-qt5
    HINTS ${PC_QWT_LIBRARY_DIRS}
    NO_DEFAULT_PATH
  )

  if(QWT_LIBRARY_FOUND)
    set(QWT_LIBRARY ${QWT_LIBRARY_FOUND})
    message(STATUS "Found Qwt library: ${QWT_LIBRARY}")
  else()
    # Fallback: let the linker find it.
    set(QWT_LIBRARY qwt)
    message(STATUS "Using fallback Qwt library: qwt")
  endif()

  set(QWT_VERSION_STRING ${PC_QWT_VERSION})
else()
  message(STATUS "pkg-config failed, falling back to path search...")
  find_library (QWT_LIBRARY NAMES qwt qwt-qt3 qwt-qt4 qwt-qt5)
endif()

set(QWT_LIBRARIES ${QWT_LIBRARY})

# --- version ---------------------------------------------------------------
if (NOT QWT_VERSION_STRING AND QWT_INCLUDE_DIR)
  set(_VERSION_FILE ${QWT_INCLUDE_DIR}/qwt_global.h)
  if (EXISTS ${_VERSION_FILE})
    file(STRINGS ${_VERSION_FILE} _VERSION_LINE REGEX "define[ ]+QWT_VERSION_STR")
    if (_VERSION_LINE)
      string(REGEX REPLACE ".*define[ ]+QWT_VERSION_STR[ ]+\"(.*)\".*" "\\1" QWT_VERSION_STRING "${_VERSION_LINE}")
    endif()
  endif()
  unset(_VERSION_FILE)
endif()

set(_QWT_VERSION_MATCH TRUE)
if (Qwt_FIND_VERSION AND QWT_VERSION_STRING)
  if (Qwt_FIND_VERSION_EXACT)
    if (NOT Qwt_FIND_VERSION VERSION_EQUAL QWT_VERSION_STRING)
      set(_QWT_VERSION_MATCH FALSE)
    endif()
  elseif (QWT_VERSION_STRING VERSION_LESS Qwt_FIND_VERSION)
    set(_QWT_VERSION_MATCH FALSE)
  endif()
endif()

# handle the QUIETLY and REQUIRED arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Qwt DEFAULT_MSG QWT_LIBRARY QWT_INCLUDE_DIR _QWT_VERSION_MATCH)
 
mark_as_advanced(
   QWT_LIBRARY
   QWT_LIBRARIES
   QWT_INCLUDE_DIR
)