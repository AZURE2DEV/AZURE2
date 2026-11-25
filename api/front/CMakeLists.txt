# The front ends' facade (docs/upgrade/qt_gui.md): azure2::front, the one header a front end in its own process (the Qt
# setup program, gui/) includes.  It links the program's run (AZURE2Program: RunProgram, --validate's and --convert's
# functions), so that a run from a front end is the program's run.
add_library(azure2_front STATIC src/FrontProject.cpp src/FrontRun.cpp src/FrontSchema.cpp src/FrontOutputs.cpp
            src/FrontParameters.cpp src/FrontThm.cpp src/FrontLevels.cpp src/FrontPhysics.cpp src/FrontForms.cpp
            src/FrontPhysicsPlots.cpp src/FrontPolish.cpp)
# The Potential panel of the Pairs tab and the Levels tab's check of the partial waves (docs/upgrade/qt_gui.md section 19).
target_sources(azure2_front PRIVATE src/FrontPotential.cpp)
# The Overview page's calls (docs/upgrade/qt_gui.md section 17).
target_sources(azure2_front PRIVATE src/FrontOverview.cpp)
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
  set_tests_properties(api_front PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" TIMEOUT 600)
endif()
