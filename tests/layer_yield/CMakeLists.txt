# The yield of a layer of finite thickness (README.md): the kernel alone on hand-chosen sub-point values (to rounding),
# and pyazr building the records and running them.  The project's limits through the CLI are check.sh
# (tests/run_tests.sh runs it).
add_executable(layer_kernel_test layer_kernel_test.cpp)
target_link_libraries(layer_kernel_test PRIVATE AZURE2Core)
add_test(NAME layer_kernel COMMAND layer_kernel_test)

find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  add_test(NAME pyazr_layer_yield
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pyazr_layer_test.py)
  set_tests_properties(pyazr_layer_yield PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;AZURE2_BIN=$<TARGET_FILE:AZURE2>"
                       SKIP_RETURN_CODE 77)
  # as tests/pyazr/CMakeLists.txt: a scratch build keeps the engine module out of the checkout; without this line the
  # test was skipped there (found November 2026)
  if(PYAZR_OUTPUT_DIR)
    set_property(TEST pyazr_layer_yield APPEND PROPERTY ENVIRONMENT "PYAZR_ENGINE_DIR=${PYAZR_OUTPUT_DIR}")
  endif()
endif()
