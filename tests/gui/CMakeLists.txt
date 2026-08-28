# Headless GUI tests.  Only meaningful when the GUI is built, since they link
# against the AZURESetup library.
#
# The Nuclear Potential tab reads its particle pairs from a model that is still
# empty when the tab is constructed, so its selector has to follow that model
# rather than take a snapshot.  Nothing in the physics regression suite can see
# that, which is how it shipped broken once.
#
# Skipped when cross-compiling for Windows: the result would not run on the
# build machine anyway.  The CI Windows job builds natively under MSYS2, so it
# still gets covered there.

if(NOT BUILD_GUI OR CROSS_COMPILE_WINDOWS)
    return()
endif()

# Qt's imported targets are found in gui/src and are scoped to that directory,
# so ask for them here as well rather than reaching into another scope.
find_package(Qt5 COMPONENTS Core Gui Widgets Network REQUIRED)

add_executable(nuclear_potential_tab_test nuclear_potential_tab_test.cpp)

# The test declares no QObject of its own; the tab's moc comes from AZURESetup.
set_target_properties(nuclear_potential_tab_test PROPERTIES AUTOMOC OFF)

target_include_directories(nuclear_potential_tab_test PRIVATE
    ${AZURE2_SOURCE_DIR}/gui/include
    ${AZURE2_SOURCE_DIR}/include
)

# Same order as the AZURE2 executable: the GUI calls into the engine, and
# AZUREAPI pulls AZURE2Core in as a PUBLIC dependency of its own.
target_link_libraries(nuclear_potential_tab_test PRIVATE
    ${AZURESETUP_LIBRARIES}
    ${AZUREAPI_LIBRARIES}
    AZURE2Core
    Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network
)

if(USE_QWT)
    target_link_libraries(nuclear_potential_tab_test PRIVATE ${QWT_LIBRARIES})
endif()

if(USE_ERYA)
    target_link_libraries(nuclear_potential_tab_test PRIVATE pugixml srim_utilities)
endif()

add_test(NAME nuclear_potential_tab COMMAND nuclear_potential_tab_test)
# No display on a CI runner, and none needed: the tab is driven through its
# widgets rather than through the window system.
set_tests_properties(nuclear_potential_tab PROPERTIES
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen")

# Round-trip of the Experimental Effects tab's optional energy-range tokens:
# read, model storage, conditional write-back, and legacy-line stability.
add_executable(target_int_tab_test target_int_tab_test.cpp)
set_target_properties(target_int_tab_test PROPERTIES AUTOMOC OFF)
target_include_directories(target_int_tab_test PRIVATE
    ${AZURE2_SOURCE_DIR}/gui/include
    ${AZURE2_SOURCE_DIR}/include
)
target_link_libraries(target_int_tab_test PRIVATE
    ${AZURESETUP_LIBRARIES}
    ${AZUREAPI_LIBRARIES}
    AZURE2Core
    Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network
)
if(USE_QWT)
    target_link_libraries(target_int_tab_test PRIVATE ${QWT_LIBRARIES})
endif()
if(USE_ERYA)
    target_link_libraries(target_int_tab_test PRIVATE pugixml srim_utilities)
endif()
add_test(NAME target_int_tab COMMAND target_int_tab_test)
set_tests_properties(target_int_tab PROPERTIES
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen")

# Moving a data segment must carry the whole row: the per-column rebuild the
# up/down buttons used dropped the UPOS fields (James deBoer, Aug 2026).
add_executable(segments_move_test segments_move_test.cpp)
set_target_properties(segments_move_test PROPERTIES AUTOMOC OFF)
target_include_directories(segments_move_test PRIVATE
    ${AZURE2_SOURCE_DIR}/gui/include
    ${AZURE2_SOURCE_DIR}/include
)
target_link_libraries(segments_move_test PRIVATE
    ${AZURESETUP_LIBRARIES}
    ${AZUREAPI_LIBRARIES}
    AZURE2Core
    Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network
)
if(USE_QWT)
    target_link_libraries(segments_move_test PRIVATE ${QWT_LIBRARIES})
endif()
if(USE_ERYA)
    target_link_libraries(segments_move_test PRIVATE pugixml srim_utilities)
endif()
add_test(NAME segments_move COMMAND segments_move_test)
set_tests_properties(segments_move PROPERTIES
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen")
