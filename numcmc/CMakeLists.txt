cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_COMPILER_ARCHITECTURE_ID "x64")
##############################################
# Declare dependencies

find_package(OpenMP REQUIRED)
##############################################
# Create target and set properties

add_library(numcmc SHARED src/ensemble.cpp src/mcmc.cpp src/walker.cpp)

#Add an alias so that library can be used inside the build tree, e.g. when testing
add_library(numcmc::numcmc ALIAS numcmc)

#Set target properties
target_include_directories(numcmc
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(numcmc
    PUBLIC
        OpenMP::OpenMP_CXX
)

target_compile_features(numcmc PRIVATE cxx_auto_type)
# target_compile_options(hitdic PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wall>)
# target_compile_options(hitdic PRIVATE -Wall -Wextra -pedantic -Werror)
target_compile_options(numcmc PRIVATE -Wall -Wno-sign-compare -Wno-unknown-pragmas)

#This is required so that the exported target has the name Hitdic and not hitdic
set_target_properties(numcmc PROPERTIES EXPORT_NAME numcmc)
