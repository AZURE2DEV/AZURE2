# Compiled tests of the upgrade (docs/upgrade/).  The test projects T1-T3 and their golden outputs are run by
# tests/upgrade/check.sh (tests/run_tests.sh); the identity checks are pyazr tests (tests/pyazr/CMakeLists.txt).

# The physics/ System built from today's CNuc (src/CNucSystem.cpp), on T1, T2 and T3: every number against the
# CNuc accessor it came from (exact), the J^pi block structure against the .azr files (docs/upgrade/physics_system.md).
# Reads the projects in place and writes nothing.
add_executable(system_from_cnuc_test system_from_cnuc_test.cpp)
target_link_libraries(system_from_cnuc_test PRIVATE AZURE2Core)
add_test(NAME physics_system_from_cnuc
         COMMAND system_from_cnuc_test
                 ${CMAKE_CURRENT_SOURCE_DIR}/T1/T1.azr
                 ${CMAKE_CURRENT_SOURCE_DIR}/T2/T2.azr
                 ${CMAKE_CURRENT_SOURCE_DIR}/T3/T3.azr
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(physics_system_from_cnuc PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 300)

# physics/ level_parameters against the engine's own transformation (CNuc::TransformIn, CalcBoundaryConditions,
# TransformOut) on T1, T2 and T3 in the Brune, Park and formal bases; round trips; convention invariance; a missing
# channel reported (docs/upgrade/level_parameters.md section 7).  Reads the projects in place and writes nothing.
add_executable(level_parameters_from_cnuc_test level_parameters_from_cnuc_test.cpp)
target_link_libraries(level_parameters_from_cnuc_test PRIVATE AZURE2Core)
add_test(NAME physics_level_parameters_from_cnuc
         COMMAND level_parameters_from_cnuc_test
                 ${CMAKE_CURRENT_SOURCE_DIR}/T1/T1.azr
                 ${CMAKE_CURRENT_SOURCE_DIR}/T2/T2.azr
                 ${CMAKE_CURRENT_SOURCE_DIR}/T3/T3.azr
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(physics_level_parameters_from_cnuc PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 600)

# physics/coulomb's default algorithm against GSL's Coulomb functions, an independent algorithm physics/ does not
# link (docs/upgrade/coulomb.md section 5).
add_executable(coulomb_gsl_crosscheck_test coulomb_gsl_crosscheck_test.cpp)
target_link_libraries(coulomb_gsl_crosscheck_test PRIVATE azure2_physics ${GSL_LIBRARIES})
target_include_directories(coulomb_gsl_crosscheck_test PRIVATE ${GSL_INCLUDE_DIRS})
add_test(NAME physics_coulomb_gsl_crosscheck COMMAND coulomb_gsl_crosscheck_test)

# The per-stage channel table (docs/upgrade/channel_table.md): rows depend on their exact key only, not on the order,
# the duplicates or the thread count; a previous table is reused only for the same system.
add_executable(channel_table_test channel_table_test.cpp)
target_link_libraries(channel_table_test PRIVATE AZURE2Core)
add_test(NAME channel_table_unit
         COMMAND channel_table_test
                 ${CMAKE_CURRENT_SOURCE_DIR}/T1/T1.azr
                 ${CMAKE_CURRENT_SOURCE_DIR}/T2/T2.azr
                 ${CMAKE_CURRENT_SOURCE_DIR}/T3/T3.azr
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(channel_table_unit PROPERTIES TIMEOUT 300)

# The observable coupling coefficients of physics/ against the engine's expressions before step 5, written with
# AngCoeff and GSL's 9j (docs/upgrade/cnuc_split.md section 9 (b)).
add_executable(observable_coupling_gsl_test observable_coupling_gsl_test.cpp)
target_link_libraries(observable_coupling_gsl_test PRIVATE AZURE2Core)
add_test(NAME physics_observable_coupling_gsl COMMAND observable_coupling_gsl_test)

# The pathway trees built from physics/ reaction_pathways against a verbatim copy of the SortPathways they replaced,
# on every .azr under tests/, Brune and formal basis (docs/upgrade/cnuc_split.md section 10 (b)).
add_executable(pathways_from_cnuc_test pathways_from_cnuc_test.cpp)
target_link_libraries(pathways_from_cnuc_test PRIVATE AZURE2Core)
add_test(NAME physics_pathways_from_cnuc COMMAND pathways_from_cnuc_test ${CMAKE_SOURCE_DIR}/tests
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(physics_pathways_from_cnuc PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 900)

# physics/ central_derivative against GSL's gsl_deriv_central, bit for bit (docs/upgrade/coulomb.md section 1.4).
add_executable(derivatives_gsl_test derivatives_gsl_test.cpp)
target_link_libraries(derivatives_gsl_test PRIVATE azure2_physics ${GSL_LIBRARIES})
target_include_directories(derivatives_gsl_test PRIVATE ${GSL_INCLUDE_DIRS})
add_test(NAME physics_derivatives_gsl COMMAND derivatives_gsl_test)

# The results of an evaluation belong to its caller (include/Results.h, docs/upgrade/data_model.md S3): working copies
# of T2, one copied from another, evaluated at different parameters one after the other and at the same time, each
# buffer bit for bit equal to a fresh copy's.  Runs on a copy of T2's input files.
add_executable(results_owner_test results_owner_test.cpp)
target_link_libraries(results_owner_test PRIVATE AZURE2Core)
set(_results_owner_dir ${CMAKE_CURRENT_BINARY_DIR}/results_owner)
file(MAKE_DIRECTORY ${_results_owner_dir}/output ${_results_owner_dir}/checks)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/T2/data ${CMAKE_CURRENT_SOURCE_DIR}/T2/cov_diff.txt
          ${CMAKE_CURRENT_SOURCE_DIR}/T2/prior_cov.txt DESTINATION ${_results_owner_dir})
add_test(NAME upgrade_results_owner COMMAND results_owner_test ${CMAKE_CURRENT_SOURCE_DIR}/T2/T2.azr
         WORKING_DIRECTORY ${_results_owner_dir})
set_tests_properties(upgrade_results_owner PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" TIMEOUT 600)

# A working copy that evaluated other parameters before gives a fresh copy's bits (docs/upgrade/data_model.md S1):
# per-segment chi^2 and model values, on T2, natural_neon, shift_sqrt_subpoints (sqrt, shiftfn), 18O_p_a_thm and a
# 7Li(p,n) point an energy shift moves out of the kinematic range, each on a copy of its input files.
add_executable(history_test history_test.cpp)
target_link_libraries(history_test PRIVATE AZURE2Core)
set(_history_dir ${CMAKE_CURRENT_BINARY_DIR}/history)
foreach(_h T2:${CMAKE_CURRENT_SOURCE_DIR}/T2 natural_neon:${CMAKE_SOURCE_DIR}/tests/natural_neon
        shift_sqrt:${CMAKE_SOURCE_DIR}/tests/shift_sqrt_subpoints/models thm:${CMAKE_SOURCE_DIR}/tests/18O_p_a_thm)
  string(REPLACE ":" ";" _hp ${_h})
  list(GET _hp 0 _name)
  list(GET _hp 1 _src)
  file(MAKE_DIRECTORY ${_history_dir}/${_name}/output ${_history_dir}/${_name}/checks)
  file(COPY ${_src}/ DESTINATION ${_history_dir}/${_name} PATTERN "output" EXCLUDE PATTERN "checks" EXCLUDE
       PATTERN "expected" EXCLUDE)
endforeach()
# pn: tests/7Li_p_ay with a lab-differential 7Li(p,n) point at 1.90 MeV, 13 deg and a fitted shift of -15 keV: at A
# (1.885 MeV, theta_max 11.7 deg with the engine's masses) the shift moves it beyond the kinematic range, at B (+1 keV,
# 15.0 deg) back inside, so a copy's angle and Jacobian at A must not be those B left (data_model.md H4).
file(MAKE_DIRECTORY ${_history_dir}/pn/output ${_history_dir}/pn/checks ${_history_dir}/pn/data)
file(READ ${CMAKE_SOURCE_DIR}/tests/7Li_p_ay/7Li_p_ay.azr _pn)
string(REPLACE "<segmentsData>\n" "<segmentsData>\n1 5 2 1.8 2.0 0 180 1 1 0 0 -0.015 0 1 data/pn13.dat 0 0\n" _pn "${_pn}")
file(WRITE ${_history_dir}/pn/pn.azr "${_pn}")
file(WRITE ${_history_dir}/pn/data/pn13.dat "1.90 13 1.0 0.1\n")
file(COPY ${CMAKE_SOURCE_DIR}/tests/7Li_p_ay/data/ DESTINATION ${_history_dir}/pn/data)
add_test(NAME upgrade_history
         COMMAND history_test ${_history_dir}/T2/T2.azr ${_history_dir}/natural_neon/natural_neon.azr
                 ${_history_dir}/shift_sqrt/sqrt.azr ${_history_dir}/shift_sqrt/shiftfn.azr
                 ${_history_dir}/thm/18O_p_a_thm.azr ${_history_dir}/pn/pn.azr
         WORKING_DIRECTORY ${_history_dir})
set_tests_properties(upgrade_history PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" TIMEOUT 900)

# A fit does not depend on the number of threads: pooled working copies give a fresh copy's bits, and fits at 1 and
# 4 threads evaluate the same vectors to the same values (tests/fit_threads/check.sh, docs/upgrade/fit_threads.md).
add_test(NAME fit_threads COMMAND bash ${CMAKE_SOURCE_DIR}/tests/fit_threads/check.sh $<TARGET_FILE:AZURE2>)
set_tests_properties(fit_threads PROPERTIES TIMEOUT 600)
# Results read before they are written (docs/upgrade/data_model.md, S3 and check 16): T1-T3 calculated, extrapolated
# and (T1) fitted with and without AZURE_POISON_RESULTS=1, outputs and the API's arrays equal; known cases listed.
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  add_test(NAME upgrade_poison_results
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/poison_results_test.py --binary $<TARGET_FILE:AZURE2>)
  set_tests_properties(upgrade_poison_results PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" TIMEOUT 1800)
endif()
