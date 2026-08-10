include_directories(../include ${AZURE2_SOURCE_DIR}/include)

set (AZUREAPI_SOURCES AZUREAPI.cpp)

set (AZUREAPI_HEADERS ../include/AZUREAPI.h)

add_library(AZUREAPI ${AZUREAPI_SOURCES} ${AZUREAPI_HEADERS_MOC})

# AZUREAPI implements the engine-facing API (AZUREAPI.cpp); the core engine
# itself lives in AZURE2Core, so the two link together in both the executable
# and the Python module below.
target_link_libraries(AZUREAPI PUBLIC AZURE2Core)

set (TEMP_VAR AZUREAPI )
set (AZUREAPI_LIBRARIES ${TEMP_VAR} PARENT_SCOPE)

# ---------------------------------------------------------------------------
# In-process Python bindings (_azure2).
#
# This replaces the old TCP/IP socket server.  The module is built against
# pybind11 and links the same AZUREAPI + AZURE2Core objects as the executable,
# so a ``pyazr.azure2`` session is a real AZUREAPI object in the interpreter.
# The module is placed in the repository's ``pyazr/`` package directory so
# ``import pyazr`` picks it up (and pip install -e . packages it).
# ---------------------------------------------------------------------------
find_package(pybind11 CONFIG QUIET)
if(NOT pybind11_FOUND)
    # pyazr cannot run without the module, but the CLI/GUI still can; build
    # those and tell the user what is missing rather than failing outright.
    message(WARNING "pybind11 not found -- skipping the _azure2 Python module. "
                    "Install pybind11 (e.g. pip install pybind11) and reconfigure "
                    "so `pyazr` can drive AZURE2 in-process.")
else()
    pybind11_add_module(_azure2 pybind_api.cpp)
    target_link_libraries(_azure2 PRIVATE
        AZUREAPI
        AZURE2Core
        ${GSL_LIBRARIES}
        Minuit2::Minuit2
        ${EXTRA_LIBS}
        ${BLAS_LIBRARIES}
        ${OpenMP_LIBRARIES})
    if(USE_NLOPT)
        target_link_libraries(_azure2 PRIVATE nlopt)
    endif()

    # Build directly into the pyazr package so ``import pyazr._azure2`` resolves
    # without any path fiddling.
    set_target_properties(_azure2 PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/pyazr)
endif()
