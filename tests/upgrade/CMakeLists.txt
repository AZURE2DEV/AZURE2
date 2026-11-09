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
