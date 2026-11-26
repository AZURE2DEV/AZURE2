# The front ends' facade (docs/upgrade/qt_gui.md): azure2::front, the one header a front end in its own process (the Qt
# setup program, gui/) includes.  It links the program's run (AZURE2Program: RunProgram, --validate's and --convert's
# functions), so that a run from a front end is the program's run.
add_library(azure2_front STATIC src/FrontProject.cpp src/FrontRun.cpp src/FrontSchema.cpp src/FrontOutputs.cpp
            src/FrontParameters.cpp src/FrontThm.cpp src/FrontLevels.cpp src/FrontPhysics.cpp src/FrontForms.cpp
            src/FrontPhysicsPlots.cpp src/FrontPolish.cpp)
# The Potential panel of the Pairs tab and the Levels tab's check of the partial waves (docs/upgrade/qt_gui.md section 20).
target_sources(azure2_front PRIVATE src/FrontPotential.cpp)
# The Overview page's calls (docs/upgrade/qt_gui.md section 17).
target_sources(azure2_front PRIVATE src/FrontOverview.cpp)
# The Coulomb table of this computer: its state and the command line's build (qt_gui.md section 19.12).
target_sources(azure2_front PRIVATE src/FrontCoulombTable.cpp)
target_include_directories(azure2_front PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
# PRIVATE: the engine's include directories do not reach a front end; only Engine.h does.
# AZUREAPI: the parameter vector of a session (ReadParameters, as pyazr's Session sets it up).
# azure2_physics_table: the rows of the physics table (PhysicsFunctions, EvaluatePhysics), the binding's too.
target_link_libraries(azure2_front PRIVATE AZURE2Program AZUREAPI azure2_physics_table)
target_include_directories(azure2_front PRIVATE ${CMAKE_SOURCE_DIR}/api/include)
set_target_properties(azure2_front PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(BUILD_TESTING)
  # Headless: no Qt.  The facade against the program on the test projects (front_test.cpp).
  add_executable(front_test tests/front_test.cpp)
  target_link_libraries(front_test PRIVATE azure2_front)
  add_test(NAME api_front COMMAND front_test ${CMAKE_SOURCE_DIR} $<TARGET_FILE:AZURE2>)
  # The facade under a decimal-comma numeric locale, without a window: read, check, edit, save, convert give the bytes
  # of the "C" locale (docs/upgrade/qt_gui.md section 19.10; skipped where the machine has no such locale).
  add_executable(front_locale_test tests/locale_test.cpp)
  target_link_libraries(front_locale_test PRIVATE azure2_front)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/locale_scratch)
  add_test(NAME api_front_locale COMMAND front_locale_test ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/locale_scratch)
  set_tests_properties(api_front_locale PROPERTIES SKIP_RETURN_CODE 77)
  # The limit is a guard against a hang, twice what the slow runner needs and more: the macOS runner took 165 to
  # 320 s in five runs of 2026-11-24.
  set_tests_properties(api_front PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" TIMEOUT 1200)
endif()
