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
# The window's shell (docs/upgrade/qt_gui.md section 11): the look, the start page, the tool and status bars, the
# Diagnostics dock, the Calculate tab's editors and summary, the plot's zoom, pan and CSV.  Its own work directory
# (gui_shell_work) and settings; its screenshots shell_*.png beside gui_window's.
add_executable(gui_shell_test gui_shell_test.cpp)
target_link_libraries(gui_shell_test PRIVATE AZURESetupLib)
add_test(NAME gui_shell COMMAND gui_shell_test ${CMAKE_SOURCE_DIR} $<TARGET_FILE:AZURE2> ${AZURE2_QT_SHOTS})
set_tests_properties(gui_shell PROPERTIES LABELS gui TIMEOUT 900 ENVIRONMENT "QT_QPA_PLATFORM=offscreen;OMP_NUM_THREADS=2")


# The Levels and Particle Pairs tabs on the go (docs/upgrade/qt_gui.md section 11): keys typed into the widgets
# (QTest's key events), the project text read back.  Needs Qt 6's Test module.
find_package(Qt6 COMPONENTS Test QUIET)
if(Qt6Test_FOUND)
  add_executable(gui_levels_test gui_levels_test.cpp)
  target_link_libraries(gui_levels_test PRIVATE AZURESetupLib Qt6::Test)
  add_test(NAME gui_levels COMMAND gui_levels_test ${CMAKE_SOURCE_DIR})
  set_tests_properties(gui_levels PROPERTIES LABELS gui TIMEOUT 600
                       ENVIRONMENT "QT_QPA_PLATFORM=offscreen;OMP_NUM_THREADS=2")
else()
  message(STATUS "Qt 6 Test not found: gui_levels is not built")
endif()

# The forms of the record tabs (docs/upgrade/qt_gui.md section 11): every record through its widgets byte for byte,
# values set through the widgets, the table, the engine's words as the user types.  AZURE2_QT_SHOTS/u7-forms: grabs.
add_executable(gui_forms_test gui_forms_test.cpp)
target_link_libraries(gui_forms_test PRIVATE AZURESetupLib)
add_test(NAME gui_forms COMMAND gui_forms_test ${CMAKE_SOURCE_DIR} ${AZURE2_QT_SHOTS}/forms)
set_tests_properties(gui_forms PROPERTIES LABELS gui TIMEOUT 600 ENVIRONMENT "QT_QPA_PLATFORM=offscreen;OMP_NUM_THREADS=2")

# The Physics tab (docs/upgrade/qt_gui.md section 15): each page through its widgets, the plotted arrays against the
# facade's calls, the pole map's cells and read-out, the pole table; screenshots AZURE2_QT_SHOTS/u8-physqt.
add_executable(gui_physqt_test gui_physqt_test.cpp)
target_link_libraries(gui_physqt_test PRIVATE AZURESetupLib)
add_test(NAME gui_physqt COMMAND gui_physqt_test ${CMAKE_SOURCE_DIR} ${AZURE2_QT_SHOTS}/u8-physqt)
set_tests_properties(gui_physqt PROPERTIES LABELS gui TIMEOUT 900 ENVIRONMENT "QT_QPA_PLATFORM=offscreen;OMP_NUM_THREADS=2")
