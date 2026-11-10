# physics/: the physical objects and stages of the R-matrix calculation (docs/upgrade/ARCHITECTURE.md sections 2-3).
# It uses nothing else of AZURE2: no experiment/, inference/, io/, engine (src/), Qt, GSL or OpenMP.  The rule is
# enforced three times: the target links nothing and sees only its own include directory; the ctest
# physics_standalone configures and builds this directory as a project of its own (tests/physics_standalone.cmake),
# where no AZURE2 header or library exists to be found; and repo_consistency (tests/consistency) refuses any
# include but azure2/physics/ headers and the C++ standard library, which also catches a GSL, Qt or OpenMP header
# that a standalone build would find on the system include path.
#
# Layout: public headers in include/azure2/physics/ (included as "azure2/physics/X.h"), sources in src/, tests that
# need only physics/ in tests/.
#
# It works as a subdirectory of AZURE2 and as a top-level project (cmake -S physics -B build-physics).
cmake_minimum_required(VERSION 3.16...4.0)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(azure2_physics LANGUAGES C CXX)
  set(AZURE2_PHYSICS_STANDALONE ON)
  enable_testing()
endif()

add_library(azure2_physics STATIC
  src/AngularCoupling.cpp
  src/ChannelFunctions.cpp
  src/Derivatives.cpp
  src/Elements.cpp
  src/ExternalCapture.cpp
  src/HybridPotential.cpp
  src/JpiBlock.cpp
  src/LevelParameters.cpp
  src/ObservableCoupling.cpp
  src/Pathways.cpp
  src/Recoupling.cpp
  src/ResonanceQuantities.cpp
  src/SymmetricEigen.cpp
  src/System.cpp)
add_library(azure2::physics ALIAS azure2_physics)
target_include_directories(azure2_physics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(azure2_physics PUBLIC cxx_std_17)
set_target_properties(azure2_physics PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_EXTENSIONS OFF)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(azure2_physics PRIVATE -Wall -Wextra -Wpedantic)
endif()

# The Coulomb functions (docs/upgrade/coulomb.md): sources from pycoulomb, one directory per algorithm.
include(coulomb/coulomb.cmake)

# Unit test of the types themselves (construction, validation, index conventions); links physics only.
add_executable(physics_system_unit_test tests/system_unit_test.cpp)
target_link_libraries(physics_system_unit_test PRIVATE azure2_physics)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(physics_system_unit_test PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME physics_system_unit COMMAND physics_system_unit_test)

# Coupling coefficients (3j, 6j, Clebsch-Gordan) and the Convention record's recoupling coefficients, against exact
# values, orthogonality, the four-3j sum and m-scheme overlaps (docs/upgrade/level_parameters.md section 7).
# The eigensolvers against their defining identities.  level_parameters on a synthetic System: the Brune <-> formal
# level-matrix identity, round trips, convention invariance, reported gaps.  external_capture: N_f, the selection
# rules, gamma^ext and the formal-basis mask (docs/upgrade/external_capture.md section 5).  The 9j and the observable
# coupling coefficients against the six-3j sum and m-scheme angular distributions (docs/upgrade/cnuc_split.md section 9).
# The resonance strength, branching ratio and total width against numbers computed by hand and their derivatives
# against central differences (docs/upgrade/data_kinds.md section 4.7).  The hybrid potential as a value and its
# inward radial integration against the Coulomb functions with V_N = 0 (docs/upgrade/spin_orbit_potential.md).
# (Sources named in full: repo_consistency checks that every tests/*_test.cpp is named here.)
foreach(_src angular_coupling_test.cpp symmetric_eigen_test.cpp level_parameters_unit_test.cpp
             external_capture_unit_test.cpp observable_coupling_test.cpp pathways_unit_test.cpp
             resonance_quantities_test.cpp hybrid_potential_test.cpp)
  string(REPLACE "_test.cpp" "" _t ${_src})
  add_executable(physics_${_t}_test tests/${_src})
  target_link_libraries(physics_${_t}_test PRIVATE azure2_physics)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(physics_${_t}_test PRIVATE -Wall -Wextra -Wpedantic)
  endif()
  add_test(NAME physics_${_t} COMMAND physics_${_t}_test)
endforeach()

# The Coulomb functions against mpmath (tests/data/coulomb.csv, pycoulomb's reference set), the Wronskian, the
# algorithms against each other, sigma_l, argument checks and thread independence (docs/upgrade/coulomb.md section 5).
find_package(Threads REQUIRED)
add_executable(physics_coulomb_test tests/coulomb_test.cpp)
target_link_libraries(physics_coulomb_test PRIVATE azure2_physics Threads::Threads)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(physics_coulomb_test PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME physics_coulomb COMMAND physics_coulomb_test ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/coulomb.csv)

# channel_functions: neutral closed forms, P, S, phi against mpmath's Coulomb functions, S' against another
# algorithm's difference, the kinematics (docs/upgrade/coulomb.md section 5).
add_executable(physics_channel_functions_test tests/channel_functions_test.cpp)
target_link_libraries(physics_channel_functions_test PRIVATE azure2_physics)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(physics_channel_functions_test PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME physics_channel_functions COMMAND physics_channel_functions_test ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/coulomb.csv)

# Bound channel functions and the derivatives: e^x K_n against mpmath (tests/data/bessel_k.csv, make_bessel_k.py),
# S(0) against the threshold limit of rho W'/W, the neutral closed forms, the band, exact rules on polynomials
# (docs/upgrade/coulomb.md section 1.4, "Bound channel functions in physics/").
add_executable(physics_bound_functions_test tests/bound_functions_test.cpp)
target_link_libraries(physics_bound_functions_test PRIVATE azure2_physics)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(physics_bound_functions_test PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME physics_bound_functions COMMAND physics_bound_functions_test ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/bessel_k.csv)

# The Whittaker functions of bound channels against mpmath (tests/data/whittaker.csv, make_whittaker_reference.py):
# W and rho W'/W, the double range, argument checks, thread independence (docs/upgrade/coulomb.md section 1.4).
add_executable(physics_whittaker_test tests/whittaker_test.cpp)
target_link_libraries(physics_whittaker_test PRIVATE azure2_physics Threads::Threads)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(physics_whittaker_test PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME physics_whittaker COMMAND physics_whittaker_test ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/whittaker.csv)

# Element symbols by Z against the AME2020 table the masses come from (pyazr/data/mass_1.mas20).
add_executable(physics_elements_test tests/elements_test.cpp)
target_link_libraries(physics_elements_test PRIVATE azure2_physics)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(physics_elements_test PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME physics_elements COMMAND physics_elements_test ${CMAKE_CURRENT_SOURCE_DIR}/../pyazr/data/mass_1.mas20)

if(NOT AZURE2_PHYSICS_STANDALONE)
  # physics/ configured, built and tested on its own, in a fresh build tree under this one.
  add_test(NAME physics_standalone
           COMMAND ${CMAKE_COMMAND}
                   -DSRC=${CMAKE_CURRENT_SOURCE_DIR}
                   -DBIN=${CMAKE_CURRENT_BINARY_DIR}/standalone
                   -DGENERATOR=${CMAKE_GENERATOR}
                   -DCXX=${CMAKE_CXX_COMPILER}
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/physics_standalone.cmake)
  set_tests_properties(physics_standalone PROPERTIES TIMEOUT 300)
endif()
