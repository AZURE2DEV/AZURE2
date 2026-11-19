# The Qt setup program, headless (docs/upgrade/qt_gui.md section 10): registered only in a build with -DBUILD_GUI=ON,
# label gui, QT_QPA_PLATFORM=offscreen.  AZURE2_QT_SHOTS: where the screenshots of every tab go.
add_executable(gui_window_test gui_window_test.cpp)
target_link_libraries(gui_window_test PRIVATE AZURESetupLib)
set(AZURE2_QT_SHOTS "${CMAKE_BINARY_DIR}/qt-shots" CACHE PATH "Where gui_window_test writes the tabs' screenshots")
add_test(NAME gui_window COMMAND gui_window_test ${CMAKE_SOURCE_DIR} $<TARGET_FILE:AZURE2> ${AZURE2_QT_SHOTS})
# The Tools menu's check runs python -m pyazr with this build's module (PYAZR_ENGINE_DIR) and this checkout's pyazr.
find_package(Python3 COMPONENTS Interpreter)
set_tests_properties(gui_window PROPERTIES LABELS gui TIMEOUT 1500
                     ENVIRONMENT "QT_QPA_PLATFORM=offscreen;OMP_NUM_THREADS=2;AZURE2_TEST_PYTHON=${Python3_EXECUTABLE};PYTHONPATH=${CMAKE_SOURCE_DIR};PYAZR_ENGINE_DIR=${PYAZR_OUTPUT_DIR}")
