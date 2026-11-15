# Compiled tests of the upgrade (docs/upgrade/).  The test projects T1-T3 and their golden outputs are run by
# tests/upgrade/check.sh (tests/run_tests.sh); the identity checks are pyazr tests (tests/pyazr/CMakeLists.txt).

# The physics/ System built from today's CompoundSystem (src/SystemOfCompound.cpp), on T1, T2 and T3: every number against the
# CompoundSystem accessor it came from (exact), the J^pi block structure against the .azr files (docs/upgrade/physics_system.md).
# Reads the projects in place and writes nothing.
add_executable(system_of_compound_test system_of_compound_test.cpp)
target_link_libraries(system_of_compound_test PRIVATE AZURE2Core)
add_test(NAME physics_system_of_compound
         COMMAND system_of_compound_test
                 ${CMAKE_CURRENT_SOURCE_DIR}/T1/T1.azr
                 ${CMAKE_CURRENT_SOURCE_DIR}/T2/T2.azr
                 ${CMAKE_CURRENT_SOURCE_DIR}/T3/T3.azr
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(physics_system_of_compound PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 300)

# physics/ level_parameters against the engine's own transformation (CompoundSystem::TransformIn, CalcBoundaryConditions,
# TransformOut) on T1, T2 and T3 in the Brune, Park and formal bases; round trips; convention invariance; a missing
# channel reported (docs/upgrade/level_parameters.md section 7).  Reads the projects in place and writes nothing.
add_executable(level_parameters_of_compound_test level_parameters_of_compound_test.cpp)
target_link_libraries(level_parameters_of_compound_test PRIVATE AZURE2Core)
add_test(NAME physics_level_parameters_of_compound
         COMMAND level_parameters_of_compound_test
                 ${CMAKE_CURRENT_SOURCE_DIR}/T1/T1.azr
                 ${CMAKE_CURRENT_SOURCE_DIR}/T2/T2.azr
                 ${CMAKE_CURRENT_SOURCE_DIR}/T3/T3.azr
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(physics_level_parameters_of_compound PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 600)

# physics/coulomb's default algorithm against GSL's Coulomb functions, an independent algorithm physics/ does not
# link (docs/upgrade/coulomb.md section 5).
add_executable(coulomb_gsl_crosscheck_test coulomb_gsl_crosscheck_test.cpp)
target_link_libraries(coulomb_gsl_crosscheck_test PRIVATE azure2_physics ${GSL_LIBRARIES})
target_include_directories(coulomb_gsl_crosscheck_test PRIVATE ${GSL_INCLUDE_DIRS})
add_test(NAME physics_coulomb_gsl_crosscheck COMMAND coulomb_gsl_crosscheck_test)

# The options of the calculation (include/ProjectOptions.h, docs/upgrade/project_format3.md section 3.12): every
# option row of io::ConfigKeys() has its place in Config and Config::Reset's default; flag/keyword > key > default;
# contradictions with a key refused.
add_executable(project_options_test project_options_test.cpp)
target_link_libraries(project_options_test PRIVATE AZURE2Core)
add_test(NAME upgrade_project_options COMMAND project_options_test)

# The planner of the shared mesh (include/MeshPlanner.h) against today's private grids: never coarser, shared,
# independent of order, local, deterministic (docs/upgrade/energy_grid.md §2.8).
add_executable(mesh_planner_test mesh_planner_test.cpp)
target_link_libraries(mesh_planner_test PRIVATE AZURE2Core)
add_test(NAME upgrade_mesh_planner COMMAND mesh_planner_test)
# A measuring tool (mesh_fold_check.cpp): a Breit-Wigner folded on today's grids and on the planned mesh against a
# converged reference; no test runs it.
add_executable(mesh_fold_check mesh_fold_check.cpp)
target_link_libraries(mesh_fold_check PRIVATE AZURE2Core)

# The per-stage channel table (docs/upgrade/channel_table.md): rows depend on their exact key only, not on the order,
# the duplicates or the thread count; a previous table is reused only for the same system.
add_executable(channel_table_test channel_table_test.cpp)
target_link_libraries(channel_table_test PRIVATE AZURE2Core)
add_test(NAME channel_table_unit
         COMMAND channel_table_test
                 ${CMAKE_CURRENT_SOURCE_DIR}/T1/T1.azr
                 ${CMAKE_CURRENT_SOURCE_DIR}/T2/T2.azr
                 ${CMAKE_CURRENT_SOURCE_DIR}/T3/T3.azr
                 ${CMAKE_CURRENT_SOURCE_DIR}/../spin_orbit_potential/spin_orbit_potential.azr
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(channel_table_unit PROPERTIES TIMEOUT 300)

# The observable coupling coefficients of physics/ against the engine's expressions before step 5, written with
# AngCoeff and GSL's 9j (docs/upgrade/cnuc_split.md section 9 (b)).
add_executable(observable_coupling_gsl_test observable_coupling_gsl_test.cpp)
target_link_libraries(observable_coupling_gsl_test PRIVATE AZURE2Core)
add_test(NAME physics_observable_coupling_gsl COMMAND observable_coupling_gsl_test)

# Properties of the pathway table checked against the loaded CompoundSystem without the builder's loops, on every .azr under
# tests/ (docs/upgrade/cnuc_split.md section 11.6 step 1).
add_executable(pathway_table_test pathway_table_test.cpp)
target_link_libraries(pathway_table_test PRIVATE AZURE2Core)
add_test(NAME physics_pathway_table COMMAND pathway_table_test ${CMAKE_SOURCE_DIR}/tests
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(physics_pathway_table PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 900)

# The external-capture exit is found by the photon pair's number, and an exit without a pathway is refused naming the
# data segment (include/ExternalCaptureExit.h; docs/upgrade/cnuc_split.md section 11.9, defect 1).
add_executable(external_capture_exit_test external_capture_exit_test.cpp)
target_link_libraries(external_capture_exit_test PRIVATE AZURE2Core)
file(COPY ${CMAKE_SOURCE_DIR}/tests/ec_exit_position/models/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ec_exit_position)
add_test(NAME upgrade_external_capture_exit
         COMMAND external_capture_exit_test ${CMAKE_CURRENT_BINARY_DIR}/ec_exit_position/neutron_pair_first.azr)
set_tests_properties(upgrade_external_capture_exit PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 300)

# The one set-up of the compound and the data: what the last part writes refuses to be read in the middle, the
# requests, the refusals by step, the structure alone (include/CompoundSystemSetUp.h; cnuc_split.md section 12.3).
add_executable(compound_setup_test compound_setup_test.cpp)
target_link_libraries(compound_setup_test PRIVATE AZURE2Core)
add_test(NAME upgrade_compound_setup COMMAND compound_setup_test ${CMAKE_SOURCE_DIR}/tests)
set_tests_properties(upgrade_compound_setup PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 300)
# A session rebuilt at another channel radius (AZUREAPI::SetRadius) equals a session opened at that radius: a rate
# session, and a project of model data alone (tests/upgrade/session_rebuild_test.cpp).
if(TARGET AZUREAPI)
  add_executable(session_rebuild_test session_rebuild_test.cpp)
  target_include_directories(session_rebuild_test PRIVATE ${CMAKE_SOURCE_DIR}/api/include)
  target_link_libraries(session_rebuild_test PRIVATE AZUREAPI)
  add_test(NAME upgrade_session_rebuild COMMAND session_rebuild_test ${CMAKE_SOURCE_DIR}/tests)
  set_tests_properties(upgrade_session_rebuild PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 300)
endif()

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

# The residual in the frame of the datum (docs/upgrade/energy_grid.md section 4): on T2, every pull equals the c.m.
# comparison to rounding, bit for bit where J = 1, at the start and with the shifts moved.  On a copy of T2's inputs.
add_executable(lab_frame_residual_test lab_frame_residual_test.cpp)
target_link_libraries(lab_frame_residual_test PRIVATE AZURE2Core)
set(_lab_frame_dir ${CMAKE_CURRENT_BINARY_DIR}/lab_frame_residual)
file(MAKE_DIRECTORY ${_lab_frame_dir}/output ${_lab_frame_dir}/checks)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/T2/data ${CMAKE_CURRENT_SOURCE_DIR}/T2/cov_diff.txt
          ${CMAKE_CURRENT_SOURCE_DIR}/T2/prior_cov.txt DESTINATION ${_lab_frame_dir})
add_test(NAME upgrade_lab_frame_residual COMMAND lab_frame_residual_test ${CMAKE_CURRENT_SOURCE_DIR}/T2/T2.azr
         WORKING_DIRECTORY ${_lab_frame_dir})
set_tests_properties(upgrade_lab_frame_residual PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" TIMEOUT 600)

# The gradient over the energy tasks: the same bits at 1 and 3 threads, and the residual pass with no free parameter.
add_executable(energy_task_gradient_test energy_task_gradient_test.cpp)
target_link_libraries(energy_task_gradient_test PRIVATE AZURE2Core)
set(_energy_grad_dir ${CMAKE_CURRENT_BINARY_DIR}/energy_task_gradient)
file(MAKE_DIRECTORY ${_energy_grad_dir}/output ${_energy_grad_dir}/checks)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/T2/data ${CMAKE_CURRENT_SOURCE_DIR}/T2/cov_diff.txt
          ${CMAKE_CURRENT_SOURCE_DIR}/T2/prior_cov.txt ${CMAKE_CURRENT_SOURCE_DIR}/T2/prior_cov.f3.txt
          DESTINATION ${_energy_grad_dir})
add_test(NAME upgrade_energy_task_gradient COMMAND energy_task_gradient_test ${CMAKE_CURRENT_SOURCE_DIR}/T2/T2.azr
         WORKING_DIRECTORY ${_energy_grad_dir})
set_tests_properties(upgrade_energy_task_gradient PROPERTIES TIMEOUT 600)

# A working copy that evaluated other parameters before gives a fresh copy's bits (docs/upgrade/data_model.md S1):
# per-segment chi^2 and model values, on T2, natural_neon, shift_sqrt_subpoints (sqrt, shiftfn), 18O_p_a_thm and a
# 7Li(p,n) point an energy shift moves out of the kinematic range, each on a copy of its input files.
add_executable(history_test history_test.cpp)
target_link_libraries(history_test PRIVATE AZURE2Core)
set(_history_dir ${CMAKE_CURRENT_BINARY_DIR}/history)
# name|directory: '|' and not ':', which would split a Windows drive letter (D:/...)
foreach(_h "T2|${CMAKE_CURRENT_SOURCE_DIR}/T2" "natural_neon|${CMAKE_SOURCE_DIR}/tests/natural_neon"
        "shift_sqrt|${CMAKE_SOURCE_DIR}/tests/shift_sqrt_subpoints/models" "thm|${CMAKE_SOURCE_DIR}/tests/18O_p_a_thm")
  string(REPLACE "|" ";" _hp "${_h}")
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
# ec: tests/ec_exit_position/models/photon_pair_first.azr (external capture into the 13N ground state, a component)
# with a fitted energy shift of 5 keV on the LUNA capture data set and the ground state's energy free.  At A and B the
# shifted points take their capture amplitudes from the memo of the shifted amplitude (Unit::GetECAmplitudeWithShift);
# B moves the final state's energy by 2 keV.  A copy that evaluated B before A must give fresh A's bits: the final
# state enters external capture at the energy of the set-up (ALevel::GetECEnergy), whatever the copy evaluated first.
file(MAKE_DIRECTORY ${_history_dir}/ec/output ${_history_dir}/ec/checks)
file(COPY ${CMAKE_SOURCE_DIR}/tests/ec_exit_position/models/data DESTINATION ${_history_dir}/ec)
file(READ ${CMAKE_SOURCE_DIR}/tests/ec_exit_position/models/photon_pair_first.azr _ec)
set(_ec_edits
    "dataset luna_hpge\n    reaction       12C(p,g)13N\n    quantity       sigma\n    unit           b\n    file           data/skowronski_luna_hpge.dat\n    columns        E_lab[MeV] _ value error\n    energy_window  0 5 MeV\n    norm           1 fixed\n"
    "dataset luna_hpge\n    reaction       12C(p,g)13N\n    quantity       sigma\n    unit           b\n    file           data/skowronski_luna_hpge.dat\n    columns        E_lab[MeV] _ value error\n    energy_window  0 5 MeV\n    norm           1 fixed\n    shift          offset 0.005 MeV fit\n"
    "\n 0.5   -1            0    1    1    1    1    2    1 " "\n 0.5   -1            0    0    1    1    1    2    1 "
    "\n 0.5   -1            0    1    1    2    1    2    1 " "\n 0.5   -1            0    0    1    2    1    2    1 ")
list(LENGTH _ec_edits _ec_n)
math(EXPR _ec_last "${_ec_n} - 1")
foreach(_i RANGE 0 ${_ec_last} 2)
  math(EXPR _j "${_i} + 1")
  list(GET _ec_edits ${_i} _from)
  list(GET _ec_edits ${_j} _to)
  string(FIND "${_ec}" "${_from}" _at)
  if(_at EQUAL -1)
    message(FATAL_ERROR "history ec model: tests/ec_exit_position/models/photon_pair_first.azr no longer holds '${_from}'")
  endif()
  string(REPLACE "${_from}" "${_to}" _ec "${_ec}")
endforeach()
file(WRITE ${_history_dir}/ec/ec.azr "${_ec}")
add_test(NAME upgrade_history
         COMMAND history_test ${_history_dir}/T2/T2.azr ${_history_dir}/natural_neon/natural_neon.azr
                 ${_history_dir}/shift_sqrt/sqrt.azr ${_history_dir}/shift_sqrt/shiftfn.azr
                 ${_history_dir}/thm/18O_p_a_thm.azr ${_history_dir}/pn/pn.azr
                 --external-capture ${_history_dir}/ec/ec.azr
         WORKING_DIRECTORY ${_history_dir})
set_tests_properties(upgrade_history PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" TIMEOUT 900)

# The calculation plan (include/CalculationPlan.h, docs/upgrade/data_model.md S5): every point of every project of
# tests/ is calculated by exactly one task, checked when the plan is built at set-up, again here and on a copy.
add_executable(plan_test plan_test.cpp)
target_link_libraries(plan_test PRIVATE AZURE2Core)
# A measuring tool for the shared mesh (node_count.cpp: units, distinct compound energies, the mesh estimate per
# folded project); no test runs it.
add_executable(node_count node_count.cpp)
target_link_libraries(node_count PRIVATE AZURE2Core)
add_test(NAME upgrade_plan COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/plan_test.sh $<TARGET_FILE:plan_test> ${CMAKE_SOURCE_DIR}/tests)
set_tests_properties(upgrade_plan PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" TIMEOUT 1800)
# ... and at run time (AZURE_PLAN_CHECK=1: writes, and every unit's shared solve against its own) on T1-T3, 13N,
# natural_neon, a THM project, a transmission and a yield, thick_target_resonance and beam_spread_ay.
add_test(NAME upgrade_plan_runtime COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/plan_runtime.sh $<TARGET_FILE:AZURE2>)
set_tests_properties(upgrade_plan_runtime PROPERTIES TIMEOUT 1800)
# The node tables at run time (AZURE_NODE_CHECK=1): every unit of every plan task at the node its slot names, on the
# stage table or the working copy's overflow table (re-placed grids, shifted points).
add_test(NAME upgrade_node_check COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/node_check.sh $<TARGET_FILE:AZURE2>)
set_tests_properties(upgrade_node_check PROPERTIES TIMEOUT 1800)

# A fit does not depend on the number of threads: pooled working copies give a fresh copy's bits, and fits at 1 and
# 4 threads evaluate the same vectors to the same values (tests/fit_threads/check.sh, docs/upgrade/fit_threads.md).
add_test(NAME fit_threads COMMAND bash ${CMAKE_SOURCE_DIR}/tests/fit_threads/check.sh $<TARGET_FILE:AZURE2>)
set_tests_properties(fit_threads PROPERTIES TIMEOUT 600)
# The four data check files of T1 are written, with their line counts (tests/check_files_data/check.sh): no other
# project sets these checks (src/DataCheckFiles.cpp, docs/upgrade/edata_split.md N2).
add_test(NAME check_files_data COMMAND bash ${CMAKE_SOURCE_DIR}/tests/check_files_data/check.sh $<TARGET_FILE:AZURE2>)
set_tests_properties(check_files_data PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" TIMEOUT 300)
# A capture final state above a particle threshold is refused by --validate and the run, with external capture and
# without; one just below the threshold runs with the start commit's chi^2 (tests/ec_final_state_unbound/check.sh;
# experiment/src/FinalStates.cpp UnboundFinalStates); about 3 s.
add_test(NAME ec_final_state_unbound COMMAND bash ${CMAKE_SOURCE_DIR}/tests/ec_final_state_unbound/check.sh $<TARGET_FILE:AZURE2>)
set_tests_properties(ec_final_state_unbound PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" TIMEOUT 300)
# The sampler sets up the compound and the data as the run does: the same check files, the same text on a refused
# last part (tests/sampler_setup/check.sh, SetUpProgramCompound in src/RunSetUp.cpp).
add_test(NAME sampler_setup_as_run COMMAND bash ${CMAKE_SOURCE_DIR}/tests/sampler_setup/check.sh $<TARGET_FILE:AZURE2>)
set_tests_properties(sampler_setup_as_run PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 300)
# The program's commands give the menu's outputs and refuse a wrong command line (tests/cli_commands/check.sh,
# src/RunRequest.cpp); about 30 s at one thread.
add_test(NAME cli_commands COMMAND bash ${CMAKE_SOURCE_DIR}/tests/cli_commands/check.sh $<TARGET_FILE:AZURE2>)
set_tests_properties(cli_commands PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 600)
# A sum in a project with two pairs of the same particles: the record's pair keys reach every term whose label names
# both, the old composite line converts to the same file, a line whose terms take both pairs is refused
# (tests/composite_pair_keys/check.sh, project_format3.md section 3.5); about 2 s.
add_test(NAME composite_pair_keys COMMAND bash ${CMAKE_SOURCE_DIR}/tests/composite_pair_keys/check.sh $<TARGET_FILE:AZURE2>)
set_tests_properties(composite_pair_keys PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 300)
# Results read before they are written (docs/upgrade/data_model.md, S3 and check 16): T1-T3 calculated, extrapolated
# and (T1) fitted with and without AZURE_POISON_RESULTS=1, outputs and the API's arrays equal; known cases listed.
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  add_test(NAME upgrade_poison_results
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/poison_results_test.py --binary $<TARGET_FILE:AZURE2>)
  set_tests_properties(upgrade_poison_results PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" TIMEOUT 1800)
endif()
# The memory of eight folded excitation functions on one mesh (docs/upgrade/data_model.md section 9): the peak resident
# memory of an extrapolation made from examples/c12c12_tumino2018 above that of the same project unfolded, bounded;
# about 2 s with a Coulomb table, 25 s without; skipped where getrusage is missing.
if(Python3_Interpreter_FOUND)
  add_test(NAME fold_memory
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/fold_memory_test.py --binary $<TARGET_FILE:AZURE2>)
  set_tests_properties(fold_memory PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 600 SKIP_RETURN_CODE 77)
endif()
# What every data set and extrapolation of tests/ and examples/ predicts (docs/upgrade/data_model.md section 3, step
# S10): `AZURE2 --print-predictions` against the committed table tests/upgrade/predictions.txt (regenerate with
# predictions_test.py --write).
if(Python3_Interpreter_FOUND)
  add_test(NAME upgrade_predictions
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/predictions_test.py --binary $<TARGET_FILE:AZURE2>)
  set_tests_properties(upgrade_predictions PROPERTIES TIMEOUT 1800)
endif()

# The capture-A_y terms built at set-up for the reactions the data sets name (docs/upgrade/cnuc_split.md section 11.5):
# every capture-A_y data set of tests/ has its terms, none for a beta entrance, the table independent of the order of
# the requests, and a capture-A_y point without its request refused.
add_executable(capture_ay_requests_test capture_ay_requests_test.cpp)
target_link_libraries(capture_ay_requests_test PRIVATE AZURE2Core)
add_test(NAME upgrade_capture_ay_requests
         COMMAND capture_ay_requests_test ${CMAKE_SOURCE_DIR}/tests ${CMAKE_CURRENT_BINARY_DIR}/capture_ay_requests)
set_tests_properties(upgrade_capture_ay_requests PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" TIMEOUT 900)

# The level block of the parameter vector (include/LevelParameterLayout.h): a layout by hand, and on projects with
# external capture, a gamma_decay record and several nuclei the layout against the groups' sizes and the names
# CompoundSystem::FillMnParams gives the block.  Reads the projects in place and writes nothing.
add_executable(level_parameter_layout_test level_parameter_layout_test.cpp)
target_link_libraries(level_parameter_layout_test PRIVATE AZURE2Core)
add_test(NAME upgrade_level_parameter_layout
         COMMAND level_parameter_layout_test
                 ${CMAKE_CURRENT_SOURCE_DIR}/T1/T1.azr
                 ${CMAKE_CURRENT_SOURCE_DIR}/T2/T2.azr
                 ${CMAKE_CURRENT_SOURCE_DIR}/T3/T3.azr
                 ${CMAKE_SOURCE_DIR}/tests/13N/13N.azr
                 ${CMAKE_SOURCE_DIR}/tests/gamma_decay/models/record.azr
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(upgrade_level_parameter_layout PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 300)

# The levels' input values (LevelInputs.h): after Fill the <levels> rows' (a gamma_decay record's partial widths at its
# photon channels), after Initialize transformed without writing a level, a physical vector at the layout's slots; in
# the Brune and the formal basis.  Reads the projects in place and writes nothing.
add_executable(level_inputs_test level_inputs_test.cpp)
target_link_libraries(level_inputs_test PRIVATE AZURE2Core)
add_test(NAME upgrade_level_inputs
         COMMAND level_inputs_test
                 ${CMAKE_CURRENT_SOURCE_DIR}/T1/T1.azr
                 ${CMAKE_CURRENT_SOURCE_DIR}/T2/T2.azr
                 ${CMAKE_SOURCE_DIR}/tests/13N/13N.azr
                 ${CMAKE_SOURCE_DIR}/tests/gamma_decay/models/record.azr
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(upgrade_level_inputs PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 300)
