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

# The <thm> block has no tab; open + save must carry it verbatim and must not
# add one to a project that had none.
add_executable(thm_block_test thm_block_test.cpp)
set_target_properties(thm_block_test PROPERTIES AUTOMOC OFF)
target_include_directories(thm_block_test PRIVATE
    ${AZURE2_SOURCE_DIR}/gui/include
    ${AZURE2_SOURCE_DIR}/include
)
target_compile_definitions(thm_block_test PRIVATE
    AZURE2_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
target_link_libraries(thm_block_test PRIVATE
    ${AZURESETUP_LIBRARIES}
    ${AZUREAPI_LIBRARIES}
    AZURE2Core
    Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network
)
if(USE_QWT)
    target_link_libraries(thm_block_test PRIVATE ${QWT_LIBRARIES})
endif()
if(USE_ERYA)
    target_link_libraries(thm_block_test PRIVATE pugixml srim_utilities)
endif()
add_test(NAME thm_block COMMAND thm_block_test)
# A modal dialog on a failed open would wait forever offscreen; fail instead.
set_tests_properties(thm_block PROPERTIES
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen"
    TIMEOUT 120)

# <parameterSettings> rows name levels and channels the way AZURE2 numbers the
# saved file, and open + save of the GUI's own output changes nothing.
add_executable(parameter_settings_test parameter_settings_test.cpp)
set_target_properties(parameter_settings_test PROPERTIES AUTOMOC OFF)
target_include_directories(parameter_settings_test PRIVATE
    ${AZURE2_SOURCE_DIR}/gui/include
    ${AZURE2_SOURCE_DIR}/include
)
target_compile_definitions(parameter_settings_test PRIVATE
    AZURE2_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
target_link_libraries(parameter_settings_test PRIVATE
    ${AZURESETUP_LIBRARIES}
    ${AZUREAPI_LIBRARIES}
    AZURE2Core
    Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network
)
if(USE_QWT)
    target_link_libraries(parameter_settings_test PRIVATE ${QWT_LIBRARIES})
endif()
if(USE_ERYA)
    target_link_libraries(parameter_settings_test PRIVATE pugixml srim_utilities)
endif()
add_test(NAME parameter_settings COMMAND parameter_settings_test)
set_tests_properties(parameter_settings PROPERTIES
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen"
    TIMEOUT 300)

# Open + save must not round the model: every <levels> number comes back
# exactly, and the engine's chi-squared of the saved file is the original's.
add_executable(level_precision_test level_precision_test.cpp)
set_target_properties(level_precision_test PROPERTIES AUTOMOC OFF)
target_include_directories(level_precision_test PRIVATE
    ${AZURE2_SOURCE_DIR}/gui/include
    ${AZURE2_SOURCE_DIR}/include
)
target_compile_definitions(level_precision_test PRIVATE
    AZURE2_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
    AZURE2_BINARY="$<TARGET_FILE:AZURE2>")
target_link_libraries(level_precision_test PRIVATE
    ${AZURESETUP_LIBRARIES}
    ${AZUREAPI_LIBRARIES}
    AZURE2Core
    Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network
)
if(USE_QWT)
    target_link_libraries(level_precision_test PRIVATE ${QWT_LIBRARIES})
endif()
if(USE_ERYA)
    target_link_libraries(level_precision_test PRIVATE pugixml srim_utilities)
endif()
add_dependencies(level_precision_test AZURE2)
add_test(NAME level_precision COMMAND level_precision_test)
set_tests_properties(level_precision PROPERTIES
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen"
    TIMEOUT 300)

# The THM options editor: options set on the dialog are written as the <thm>
# block, come back on reopening, an untouched dialog keeps the block byte for
# byte, and the engine runs the GUI's block (weight file included).
add_executable(thm_options_dialog_test thm_options_dialog_test.cpp)
set_target_properties(thm_options_dialog_test PROPERTIES AUTOMOC OFF)
target_include_directories(thm_options_dialog_test PRIVATE
    ${AZURE2_SOURCE_DIR}/gui/include
    ${AZURE2_SOURCE_DIR}/include
)
target_compile_definitions(thm_options_dialog_test PRIVATE
    AZURE2_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
    AZURE2_BINARY="$<TARGET_FILE:AZURE2>")
target_link_libraries(thm_options_dialog_test PRIVATE
    ${AZURESETUP_LIBRARIES}
    ${AZUREAPI_LIBRARIES}
    AZURE2Core
    Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network
)
if(USE_QWT)
    target_link_libraries(thm_options_dialog_test PRIVATE ${QWT_LIBRARIES})
endif()
if(USE_ERYA)
    target_link_libraries(thm_options_dialog_test PRIVATE pugixml srim_utilities)
endif()
add_dependencies(thm_options_dialog_test AZURE2)
add_test(NAME thm_options_dialog COMMAND thm_options_dialog_test)
set_tests_properties(thm_options_dialog PROPERTIES
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen"
    TIMEOUT 300)

# The THM workspace (Configure > THM Workspace...): its Model, Experiments and
# Channels pages write the <thm> block and the THM columns of <levels>, come
# back on reopening, refuse what the engine refuses, leave an untouched
# project byte for byte, and show the energies the engine prints.
add_executable(thm_workspace_test thm_workspace_test.cpp)
set_target_properties(thm_workspace_test PROPERTIES AUTOMOC OFF)
target_include_directories(thm_workspace_test PRIVATE
    ${AZURE2_SOURCE_DIR}/gui/include
    ${AZURE2_SOURCE_DIR}/include
)
target_compile_definitions(thm_workspace_test PRIVATE
    AZURE2_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
    AZURE2_BINARY="$<TARGET_FILE:AZURE2>")
target_link_libraries(thm_workspace_test PRIVATE
    ${AZURESETUP_LIBRARIES}
    ${AZUREAPI_LIBRARIES}
    AZURE2Core
    Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network
)
if(USE_QWT)
    target_link_libraries(thm_workspace_test PRIVATE ${QWT_LIBRARIES})
endif()
if(USE_ERYA)
    target_link_libraries(thm_workspace_test PRIVATE pugixml srim_utilities)
endif()
add_dependencies(thm_workspace_test AZURE2)
add_test(NAME thm_workspace COMMAND thm_workspace_test)
set_tests_properties(thm_workspace PROPERTIES
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen"
    TIMEOUT 300)

# A pair is its key, as in the engine: examples/c12c12_tumino2018's pairs 1
# and 6 (12C+12C, with and without the THM binding energy) both open, every
# channel keeps its pair, open + save changes no number, the GUI's own output
# saves byte for byte, an edit round-trips, the engine's chi-squared is the
# original's, and the THM workspace leaves the project alone.
add_executable(pair_identity_test pair_identity_test.cpp)
set_target_properties(pair_identity_test PROPERTIES AUTOMOC OFF)
target_include_directories(pair_identity_test PRIVATE
    ${AZURE2_SOURCE_DIR}/gui/include
    ${AZURE2_SOURCE_DIR}/include
)
target_compile_definitions(pair_identity_test PRIVATE
    AZURE2_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
    AZURE2_BINARY="$<TARGET_FILE:AZURE2>")
target_link_libraries(pair_identity_test PRIVATE
    ${AZURESETUP_LIBRARIES}
    ${AZUREAPI_LIBRARIES}
    AZURE2Core
    Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network
)
if(USE_QWT)
    target_link_libraries(pair_identity_test PRIVATE ${QWT_LIBRARIES})
endif()
if(USE_ERYA)
    target_link_libraries(pair_identity_test PRIVATE pugixml srim_utilities)
endif()
add_dependencies(pair_identity_test AZURE2)
add_test(NAME pair_identity COMMAND pair_identity_test)
# A modal dialog on a failed open would wait forever offscreen; fail instead.
set_tests_properties(pair_identity PROPERTIES
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen"
    TIMEOUT 300)
