# Configure, build and test physics/ as a project of its own (cmake -P; called by the ctest physics_standalone).
# In a fresh tree nothing else of AZURE2 exists: an include of an engine header (CompoundSystem.h, Config.h) or a call into
# the engine or a library it links fails to compile or link here, whatever the main build puts on its include path.
# A header on the system include path (GSL, say) is still found; repo_consistency refuses those includes.
#   cmake -DSRC=<physics dir> -DBIN=<scratch build dir> [-DGENERATOR=...] [-DCXX=...] -P physics_standalone.cmake
# experiment/ uses the same script for itself (-DNAME=experiment, -DTEST_REGEX=^experiment_: physics/ is built there
# only as the library it links, its tests are not).
if(NOT SRC OR NOT BIN)
  message(FATAL_ERROR "physics_standalone.cmake: SRC and BIN are required")
endif()
if(NOT NAME)
  set(NAME physics)
endif()
set(_regex)
if(TEST_REGEX)
  set(_regex -R "${TEST_REGEX}")
endif()
file(REMOVE_RECURSE "${BIN}")
set(_gen)
if(GENERATOR)
  set(_gen -G "${GENERATOR}")
endif()
set(_cxx)
if(CXX)
  set(_cxx -DCMAKE_CXX_COMPILER=${CXX})
endif()
execute_process(COMMAND ${CMAKE_COMMAND} -S "${SRC}" -B "${BIN}" ${_gen} ${_cxx} -DCMAKE_BUILD_TYPE=Release
                RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "${NAME}/ alone: configure failed (${rc})")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build "${BIN}" RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "${NAME}/ alone: build failed (${rc})")
endif()
execute_process(COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure ${_regex} WORKING_DIRECTORY "${BIN}" RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "${NAME}/ alone: tests failed (${rc})")
endif()
file(REMOVE_RECURSE "${BIN}")
message(STATUS "${NAME}/ alone: configured, built and tested")
