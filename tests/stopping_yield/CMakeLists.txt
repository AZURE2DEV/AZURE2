# tests/stopping_yield (README.md): the convention of the stopping per active nucleus, against numbers put together by
# hand from the library's per-atom values.  check.sh (tests/run_tests.sh) runs the projects.
add_executable(stopping_convention_test stopping_convention_test.cpp)
target_link_libraries(stopping_convention_test PRIVATE AZURE2Core)
add_test(NAME stopping_convention COMMAND stopping_convention_test)
# AZUREParams::Add refuses a second parameter of a name (duplicate_parameter_test.cpp).
add_executable(duplicate_parameter_test duplicate_parameter_test.cpp)
target_link_libraries(duplicate_parameter_test PRIVATE AZURE2Core)
add_test(NAME duplicate_parameter COMMAND duplicate_parameter_test)

# The names of the stopping scales (scale_names_test.py): a target with `phase solid` beside one without a phase has a
# scale of its own name, and two scales of one name are refused (the 789 parameters / 790 records of 2026-11-18).
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  add_test(NAME stopping_scale_names
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scale_names_test.py)
  set_tests_properties(stopping_scale_names PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>"
                       SKIP_RETURN_CODE 77)
  if(PYAZR_OUTPUT_DIR)
    set_property(TEST stopping_scale_names APPEND PROPERTY ENVIRONMENT "PYAZR_ENGINE_DIR=${PYAZR_OUTPUT_DIR}")
  endif()
endif()
