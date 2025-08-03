CMAKE_MINIMUM_REQUIRED(VERSION 4.0)

PROJECT(stragg VERSION 1.0.0 DESCRIPTION "Straggling Utilities")

cmake_policy(SET CMP0078 OLD)
cmake_policy(SET CMP0086 OLD)

FIND_PACKAGE(SWIG REQUIRED)
FIND_PACKAGE(Python3 COMPONENTS Development)

INCLUDE(${SWIG_USE_FILE})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${PYTHON_INCLUDE_PATH})

ADD_SUBDIRECTORY(pugixml)

SET_SOURCE_FILES_PROPERTIES(python/stragg.i PROPERTIES CPLUSPLUS ON)

SWIG_ADD_LIBRARY(stragg LANGUAGE python SOURCES python/stragg.i 
    src/ChebyshevPolynomial.cc
    src/HermitePolynomial.cc
    src/GaussFunction.cc
    src/DiracFunction.cc
    src/LandauFunction.cc
    src/VavilovAiryFunction.cc
    src/VavilovEdgeworthFunction.cc
    src/VavilovMoyalFunction.cc
    src/PhysicsDistribution.cc
    src/ResonanceFunction.cc
    src/Layer.cc
    src/SRIM.cc
    src/SRIMTable.cc)
SWIG_LINK_LIBRARIES(stragg Python3::Module pugixml)

EXECUTE_PROCESS(COMMAND python3 -c "from distutils.sysconfig import get_python_lib; print( get_python_lib() )" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
INSTALL(TARGETS _stragg DESTINATION ${PYTHON_SITE_PACKAGES})
INSTALL(FILES ${CMAKE_BINARY_DIR}/stragg.py DESTINATION ${PYTHON_SITE_PACKAGES})
INSTALL(FILES data/SRIM2013.xml DESTINATION /opt/)