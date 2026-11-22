# rmfit's own tests: one ctest per script of this directory, named rmfit_<script without test_ and .py>.
#
# Each script is run as it is by hand (python3 rmfit/tests/test_X.py); it prints its checks and exits non-zero when
# one fails, 77 when it cannot run here (no AZURE2 binary).  AZURE2_ROOT names this source tree, whose pyazr/ is the
# package the tests import, with the engine module the build made (in pyazr/, or in PYAZR_OUTPUT_DIR, which
# PYAZR_ENGINE_DIR then names), so the tests exercise this build and not a pyazr elsewhere on the machine.
# rmfit's engine pin (enginepin.py) looks for the module the same way: $AZURE2_ROOT/pyazr first, then
# $PYAZR_ENGINE_DIR, which it then records and its job-script guards check (2026-11-19).
#
# rmfit runs on pyazr 2 and pyazr 3 (rmfit/backend.py).  A script that is the same on both is test_<name>.py; where the
# two differ, test_<name>.py is the pyazr 2 version and test_p3_<name>.py the pyazr 3 one, and each exits 77 (skipped)
# on the other backend, so this directory's ctest list holds both and the build's pyazr picks.
#
# The scripts that start engine sessions on the 13N project (campaign rounds, polishes, level scans, the sampler)
# take minutes to tens of minutes each: label slow.  The scripts listed in RMFIT_FAST take seconds and run with
# `ctest -LE slow`.

find_package(Python3 COMPONENTS Interpreter)
if(NOT Python3_Interpreter_FOUND)
    message(STATUS "python3 not found; skipping the rmfit tests")
    return()
endif()

set(RMFIT_FAST
    classify_broad classify_gate classify_width cluster_rank cluster_width dataset_table_split effective_k enginepin
    flip_floor keys ledger_signs new_level_cap nudge_zero_widths pair_family patterns_signs_memory point_weights
    sanity_same_jpi scoring_segments separate_close_levels sign_distance site_mail site_shell site_slurm
    split_evidence template_channels total_width_floor uq_rerun width_bound width_bound_pins width_floor window_delta
    window_patterns threshold_guard azure2_root spec_decisions seed_omitted audit_subagent py_compat drive_ensdf_flags)

file(GLOB RMFIT_TEST_SCRIPTS CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/test_*.py)
foreach(script ${RMFIT_TEST_SCRIPTS})
    get_filename_component(stem ${script} NAME_WE)
    string(REGEX REPLACE "^test_" "" name ${stem})
    string(REGEX REPLACE "^p3_" "" base ${name})
    add_test(NAME rmfit_${name} COMMAND ${Python3_EXECUTABLE} ${script})
    set_tests_properties(rmfit_${name} PROPERTIES
                         ENVIRONMENT "AZURE2_ROOT=${CMAKE_SOURCE_DIR};OMP_NUM_THREADS=1;AZURE2_BIN=$<TARGET_FILE:AZURE2>"
                         SKIP_RETURN_CODE 77)
    # A build with -DPYAZR_OUTPUT_DIR (work/build.sh, sanitizer and scratch builds) writes _azure2 there, not into
    # pyazr/: PYAZR_ENGINE_DIR points the package at it, as tests/pyazr/CMakeLists.txt does for the pyazr tests.
    if(PYAZR_OUTPUT_DIR)
        set_property(TEST rmfit_${name} APPEND PROPERTY ENVIRONMENT "PYAZR_ENGINE_DIR=${PYAZR_OUTPUT_DIR}")
    endif()
    list(FIND RMFIT_FAST ${base} fast)
    if(fast EQUAL -1)
        set_tests_properties(rmfit_${name} PROPERTIES LABELS "rmfit;slow" TIMEOUT 3600)
    else()
        set_tests_properties(rmfit_${name} PROPERTIES LABELS "rmfit" TIMEOUT 300)
    endif()
endforeach()
