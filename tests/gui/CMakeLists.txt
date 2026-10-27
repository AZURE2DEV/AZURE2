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

# The Add/Edit Level dialog lists only the J^pi values that give channels under
# the current pairs and limits, with a fallback when no pairs are defined.
add_executable(level_dialog_test level_dialog_test.cpp)
set_target_properties(level_dialog_test PROPERTIES AUTOMOC OFF)
target_include_directories(level_dialog_test PRIVATE
    ${AZURE2_SOURCE_DIR}/gui/include
    ${AZURE2_SOURCE_DIR}/include
)
target_link_libraries(level_dialog_test PRIVATE
    ${AZURESETUP_LIBRARIES}
    ${AZUREAPI_LIBRARIES}
    AZURE2Core
    Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network
)
if(USE_QWT)
    target_link_libraries(level_dialog_test PRIVATE ${QWT_LIBRARIES})
endif()
if(USE_ERYA)
    target_link_libraries(level_dialog_test PRIVATE pugixml srim_utilities)
endif()
add_test(NAME level_dialog COMMAND level_dialog_test)
set_tests_properties(level_dialog PROPERTIES
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen")

# Not a ctest: a manual backward-compatibility check of the Levels tab against
# real project files (`levels_compat_check a.azr b.azr ...`).
add_executable(levels_compat_check levels_compat_check.cpp)
set_target_properties(levels_compat_check PROPERTIES AUTOMOC OFF)
target_include_directories(levels_compat_check PRIVATE
    ${AZURE2_SOURCE_DIR}/gui/include
    ${AZURE2_SOURCE_DIR}/include
)
target_link_libraries(levels_compat_check PRIVATE
    ${AZURESETUP_LIBRARIES}
    ${AZUREAPI_LIBRARIES}
    AZURE2Core
    Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network
)
if(USE_QWT)
    target_link_libraries(levels_compat_check PRIVATE ${QWT_LIBRARIES})
endif()
if(USE_ERYA)
    target_link_libraries(levels_compat_check PRIVATE pugixml srim_utilities)
endif()

# The particle-pair drop-downs that replaced typed pair keys: labels, key
# round trip, exit items relabelled as reactions with the entrance pair.
add_executable(pair_key_combo_test pair_key_combo_test.cpp)
set_target_properties(pair_key_combo_test PROPERTIES AUTOMOC OFF)
target_include_directories(pair_key_combo_test PRIVATE
    ${AZURE2_SOURCE_DIR}/gui/include
    ${AZURE2_SOURCE_DIR}/include
)
target_link_libraries(pair_key_combo_test PRIVATE
    ${AZURESETUP_LIBRARIES}
    ${AZUREAPI_LIBRARIES}
    AZURE2Core
    Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network
)
if(USE_QWT)
    target_link_libraries(pair_key_combo_test PRIVATE ${QWT_LIBRARIES})
endif()
if(USE_ERYA)
    target_link_libraries(pair_key_combo_test PRIVATE pugixml srim_utilities)
endif()
add_test(NAME pair_key_combo COMMAND pair_key_combo_test)
set_tests_properties(pair_key_combo PROPERTIES
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen")

# Segment lines round-trip through the Segments tab with their trailing
# keyword blocks (sqrtshift, natural) intact; the sqrtshift block used to be
# dropped on load.
add_executable(segments_roundtrip_test segments_roundtrip_test.cpp)
set_target_properties(segments_roundtrip_test PROPERTIES AUTOMOC OFF)
target_include_directories(segments_roundtrip_test PRIVATE
    ${AZURE2_SOURCE_DIR}/gui/include
    ${AZURE2_SOURCE_DIR}/include
)
target_link_libraries(segments_roundtrip_test PRIVATE
    ${AZURESETUP_LIBRARIES}
    ${AZUREAPI_LIBRARIES}
    AZURE2Core
    Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network
)
if(USE_QWT)
    target_link_libraries(segments_roundtrip_test PRIVATE ${QWT_LIBRARIES})
endif()
if(USE_ERYA)
    target_link_libraries(segments_roundtrip_test PRIVATE pugixml srim_utilities)
endif()
add_test(NAME segments_roundtrip COMMAND segments_roundtrip_test)
set_tests_properties(segments_roundtrip PROPERTIES
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen")

# Several compound nuclei in one project: nuclei of the pairs, the abundance
# table and its <abundances> block, levels per nucleus, the "natural 1"
# segment keyword in the dialogs and models, abundance fitting parameters.
add_executable(natural_target_test natural_target_test.cpp)
set_target_properties(natural_target_test PROPERTIES AUTOMOC OFF)
target_include_directories(natural_target_test PRIVATE
    ${AZURE2_SOURCE_DIR}/gui/include
    ${AZURE2_SOURCE_DIR}/include
)
target_link_libraries(natural_target_test PRIVATE
    ${AZURESETUP_LIBRARIES}
    ${AZUREAPI_LIBRARIES}
    AZURE2Core
    Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network
)
if(USE_QWT)
    target_link_libraries(natural_target_test PRIVATE ${QWT_LIBRARIES})
endif()
if(USE_ERYA)
    target_link_libraries(natural_target_test PRIVATE pugixml srim_utilities)
endif()
add_test(NAME natural_target COMMAND natural_target_test)
set_tests_properties(natural_target PROPERTIES
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen")

# The run thread's Config copy must carry every option the GUI sets on its
# own Config (the hard-sphere option was missing).
add_executable(main_thread_config_test main_thread_config_test.cpp)
set_target_properties(main_thread_config_test PROPERTIES AUTOMOC OFF)
target_include_directories(main_thread_config_test PRIVATE
    ${AZURE2_SOURCE_DIR}/gui/include
    ${AZURE2_SOURCE_DIR}/include
)
target_link_libraries(main_thread_config_test PRIVATE
    ${AZURESETUP_LIBRARIES}
    ${AZUREAPI_LIBRARIES}
    AZURE2Core
    Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network
)
if(USE_QWT)
    target_link_libraries(main_thread_config_test PRIVATE ${QWT_LIBRARIES})
endif()
if(USE_ERYA)
    target_link_libraries(main_thread_config_test PRIVATE pugixml srim_utilities)
endif()
add_test(NAME main_thread_config COMMAND main_thread_config_test)
set_tests_properties(main_thread_config PROPERTIES
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen")

# <parameterSettings> "Level N ..." entries in the engine's numbering of the file they are written into: a
# project whose levels are not in the written order, with an inactive level first, through a GUI read/write,
# checked against CNuc::Fill of the written file; the older row-numbered block; a second read/write.
add_executable(parameter_settings_numbering_test parameter_settings_numbering_test.cpp)
set_target_properties(parameter_settings_numbering_test PROPERTIES AUTOMOC OFF)
target_include_directories(parameter_settings_numbering_test PRIVATE
    ${AZURE2_SOURCE_DIR}/gui/include
    ${AZURE2_SOURCE_DIR}/include
)
target_link_libraries(parameter_settings_numbering_test PRIVATE
    ${AZURESETUP_LIBRARIES}
    ${AZUREAPI_LIBRARIES}
    AZURE2Core
    Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network
)
if(USE_QWT)
    target_link_libraries(parameter_settings_numbering_test PRIVATE ${QWT_LIBRARIES})
endif()
if(USE_ERYA)
    target_link_libraries(parameter_settings_numbering_test PRIVATE pugixml srim_utilities)
endif()
add_test(NAME parameter_settings_numbering COMMAND parameter_settings_numbering_test)
set_tests_properties(parameter_settings_numbering PROPERTIES
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen")

# The Segments tab's "Measured Scattering Lengths" table and its optional
# <scatteringLengths> block: allowed channel spins, round trip, flags for rows
# whose pair is gone, reset, and a whole project file through AZURESetup.
add_executable(scattering_lengths_test scattering_lengths_test.cpp)
set_target_properties(scattering_lengths_test PROPERTIES AUTOMOC OFF)
target_include_directories(scattering_lengths_test PRIVATE
    ${AZURE2_SOURCE_DIR}/gui/include
    ${AZURE2_SOURCE_DIR}/include
)
target_compile_definitions(scattering_lengths_test PRIVATE AZURE2_TEST_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
target_link_libraries(scattering_lengths_test PRIVATE
    ${AZURESETUP_LIBRARIES}
    ${AZUREAPI_LIBRARIES}
    AZURE2Core
    Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network
)
if(USE_QWT)
    target_link_libraries(scattering_lengths_test PRIVATE ${QWT_LIBRARIES})
endif()
if(USE_ERYA)
    target_link_libraries(scattering_lengths_test PRIVATE pugixml srim_utilities)
endif()
add_test(NAME scattering_lengths COMMAND scattering_lengths_test)
set_tests_properties(scattering_lengths PROPERTIES
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen")

# The segment dialogs' observable controls: a short observable list, the
# Lab / C.M. frame of a differential cross section, the E1 / E2 multipolarity
# of angle-integrated capture and the "Total capture" exit item, all mapping
# onto the unchanged isDiff codes, with edit round trips for every code.
add_executable(segment_observable_test segment_observable_test.cpp)
set_target_properties(segment_observable_test PROPERTIES AUTOMOC OFF)
target_include_directories(segment_observable_test PRIVATE
    ${AZURE2_SOURCE_DIR}/gui/include
    ${AZURE2_SOURCE_DIR}/include
)
target_link_libraries(segment_observable_test PRIVATE
    ${AZURESETUP_LIBRARIES}
    ${AZUREAPI_LIBRARIES}
    AZURE2Core
    Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network
)
if(USE_QWT)
    target_link_libraries(segment_observable_test PRIVATE ${QWT_LIBRARIES})
endif()
if(USE_ERYA)
    target_link_libraries(segment_observable_test PRIVATE pugixml srim_utilities)
endif()
add_test(NAME segment_observable COMMAND segment_observable_test)
set_tests_properties(segment_observable PROPERTIES
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen")

# The EXFOR client's ratio-to-Rutherford conversion: a lab ANG column is taken to the centre of mass and
# divided by the lab-to-c.m. factor; an ANG-CM column is used as is (audit item 8, GUI sibling of pyazr.nds).
add_executable(exfor_rutherford_test exfor_rutherford_test.cpp)
set_target_properties(exfor_rutherford_test PROPERTIES AUTOMOC OFF)
target_include_directories(exfor_rutherford_test PRIVATE
    ${AZURE2_SOURCE_DIR}/gui/include
    ${AZURE2_SOURCE_DIR}/include
)
target_link_libraries(exfor_rutherford_test PRIVATE
    ${AZURESETUP_LIBRARIES}
    ${AZUREAPI_LIBRARIES}
    AZURE2Core
    Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network
)
if(USE_QWT)
    target_link_libraries(exfor_rutherford_test PRIVATE ${QWT_LIBRARIES})
endif()
if(USE_ERYA)
    target_link_libraries(exfor_rutherford_test PRIVATE pugixml srim_utilities)
endif()
add_test(NAME exfor_rutherford COMMAND exfor_rutherford_test)
set_tests_properties(exfor_rutherford PROPERTIES
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen")

# The Fitting tab's Screening rows: screening_<n> per active Experimental Effects line with electron
# screening, synced both ways with that tab and written back to the project file.
add_executable(screening_fitting_test screening_fitting_test.cpp)
set_target_properties(screening_fitting_test PROPERTIES AUTOMOC OFF)
target_include_directories(screening_fitting_test PRIVATE
    ${AZURE2_SOURCE_DIR}/gui/include
    ${AZURE2_SOURCE_DIR}/include
)
target_compile_definitions(screening_fitting_test PRIVATE AZURE2_TEST_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
target_link_libraries(screening_fitting_test PRIVATE
    ${AZURESETUP_LIBRARIES}
    ${AZUREAPI_LIBRARIES}
    AZURE2Core
    Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network
)
if(USE_QWT)
    target_link_libraries(screening_fitting_test PRIVATE ${QWT_LIBRARIES})
endif()
if(USE_ERYA)
    target_link_libraries(screening_fitting_test PRIVATE pugixml srim_utilities)
endif()
add_test(NAME screening_fitting COMMAND screening_fitting_test)
set_tests_properties(screening_fitting PROPERTIES
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen")
