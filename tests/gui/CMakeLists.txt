# The Qt setup program, headless (docs/upgrade/qt_gui.md section 5): registered only in a build with -DBUILD_GUI=ON,
# label gui, QT_QPA_PLATFORM=offscreen.  AZURE2_QT_SHOTS: where the screenshots of every tab go.
add_executable(gui_window_test gui_window_test.cpp)
target_link_libraries(gui_window_test PRIVATE AZURESetupLib)
set(AZURE2_QT_SHOTS "${CMAKE_BINARY_DIR}/qt-shots" CACHE PATH "Where gui_window_test writes the tabs' screenshots")
add_test(NAME gui_window COMMAND gui_window_test ${CMAKE_SOURCE_DIR} $<TARGET_FILE:AZURE2> ${AZURE2_QT_SHOTS})
set_tests_properties(gui_window PROPERTIES LABELS gui TIMEOUT 900
                     ENVIRONMENT "QT_QPA_PLATFORM=offscreen;OMP_NUM_THREADS=2")
