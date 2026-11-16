# inference/: the fit's parameter vector as plain data, its priors and the rule of the objective's term list
# (docs/upgrade/inference.md): what a parameter is (the kinds' table and the record of one parameter), where the
# level block and every other parameter sit in the full and the packed vector, the support of the posterior, one
# Gaussian prior, the list of the objective's diagonal priors, the prior covariance (grammar, set-up, arithmetic) and the
# sums of a list of terms over a context.  It uses nothing else of AZURE2: no engine header (src/, include/,
# api/), no other library, no Minuit2, GSL or OpenMP.  The rule is enforced three times, as for physics/: the target
# links nothing; the ctest inference_standalone configures and builds this directory as a project of its own; and
# repo_consistency refuses any include but azure2/inference/ and the C++ standard library.  The engine keeps the
# walks that fill these types (BuildParamIndexMap, PosteriorSupportOfFit), the values of the vector (AZUREParams), the
# terms themselves (src/Objective.cpp) and the opening of the prior covariance file (src/PriorCovarianceFile.cpp).
#
# Layout: public headers in include/azure2/inference/ (included as "azure2/inference/X.h"), sources in src/, tests in
# tests/.
cmake_minimum_required(VERSION 3.16...4.0)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(azure2_inference LANGUAGES CXX)
  set(AZURE2_INFERENCE_STANDALONE ON)
  enable_testing()
endif()

add_library(azure2_inference STATIC
  src/LevelParameterLayout.cpp
  src/LinearProfile.cpp
  src/ParameterIndexMap.cpp
  src/PosteriorSupport.cpp
  src/PriorCovariance.cpp)
add_library(azure2::inference ALIAS azure2_inference)
target_include_directories(azure2_inference PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(azure2_inference PUBLIC cxx_std_17)
set_target_properties(azure2_inference PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_EXTENSIONS OFF)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(azure2_inference PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Each test builds its objects by hand and needs no project file.
#   parameter_kinds          one row per kind in the order of the codes, keys unique, the record's factories
#   level_parameter_layout   three J groups of unequal size: every slot's index, name and level number
#   parameter_index_map      a vector with a parameter of every kind: each lookup, packed <-> full with fixed entries
#   posterior_support        inside, on and outside a limit; a vector too short; a component that is not a number
#   gaussian_prior           value by hand; gradient and curvature against central differences
#   prior_covariance         a 2 x 2 and a 3 x 3 covariance by hand: chi2, gradient, rows, marginalisation, the
#                            normalisation, each refusal of the grammar
#   term_list                terms that add known numbers: the order of the sum, the breakdown, the sampler's sum
#   diagonal_priors          a list of five priors by hand: the running sum and the subtotals bit for bit, the
#                            gradient at p and the curvature against central differences, an empty list
#   abundance_priors         two samples of two nuclei by hand (a replaced prior, a sample without a sum term): the
#                            value bit for bit the objective's loop, the gradient and the sum rows against central
#                            differences, the rows' squares = the sums' part
#   park_wall                one level of two particle channels with a linear S'(E): J by hand, the value, the
#                            gradient's and the row's forms against central differences, J >= 0 and NaN off the wall
#   linear_profile           the line s m + a_0 + a_1 E + a_2 E^2 by hand (0-3 terms), chi2 against a scan, a singular
#                            normal matrix, s <= 0, the profiled norm against a scan and by hand,
#                            ThmProfileDerivative against central differences
# (Sources named in full: repo_consistency checks that every tests/*_test.cpp is named here.)
foreach(_src parameter_kinds_test.cpp level_parameter_layout_test.cpp parameter_index_map_test.cpp
             posterior_support_test.cpp gaussian_prior_test.cpp prior_covariance_test.cpp term_list_test.cpp
             diagonal_priors_test.cpp abundance_priors_test.cpp park_wall_test.cpp linear_profile_test.cpp)
  string(REPLACE "_test.cpp" "" _t ${_src})
  add_executable(inference_${_t}_test tests/${_src})
  target_link_libraries(inference_${_t}_test PRIVATE azure2_inference)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(inference_${_t}_test PRIVATE -Wall -Wextra -Wpedantic)
  endif()
  add_test(NAME inference_${_t} COMMAND inference_${_t}_test)
endforeach()

if(NOT AZURE2_INFERENCE_STANDALONE)
  # inference/ configured, built and tested on its own, in a fresh build tree under this one.
  add_test(NAME inference_standalone
           COMMAND ${CMAKE_COMMAND}
                   -DSRC=${CMAKE_CURRENT_SOURCE_DIR}
                   -DBIN=${CMAKE_CURRENT_BINARY_DIR}/standalone
                   -DGENERATOR=${CMAKE_GENERATOR}
                   -DCXX=${CMAKE_CXX_COMPILER}
                   -DNAME=inference
                   -DTEST_REGEX=^inference_
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/../physics/tests/physics_standalone.cmake)
  set_tests_properties(inference_standalone PROPERTIES TIMEOUT 600)
endif()
