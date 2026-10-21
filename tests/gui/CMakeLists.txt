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

# THM is opt-in (Configure > Runtime Options), as the nuclear potential is:
# off for a classic project, on for one with THM content; switching it hides
# or shows the THM controls and changes no byte of a save.
add_executable(thm_opt_in_test thm_opt_in_test.cpp)
set_target_properties(thm_opt_in_test PROPERTIES AUTOMOC OFF)
target_include_directories(thm_opt_in_test PRIVATE
    ${AZURE2_SOURCE_DIR}/gui/include
    ${AZURE2_SOURCE_DIR}/include
)
target_compile_definitions(thm_opt_in_test PRIVATE
    AZURE2_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
target_link_libraries(thm_opt_in_test PRIVATE
    ${AZURESETUP_LIBRARIES}
    ${AZUREAPI_LIBRARIES}
    AZURE2Core
    Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network
)
if(USE_QWT)
    target_link_libraries(thm_opt_in_test PRIVATE ${QWT_LIBRARIES})
endif()
if(USE_ERYA)
    target_link_libraries(thm_opt_in_test PRIVATE pugixml srim_utilities)
endif()
add_test(NAME thm_opt_in COMMAND thm_opt_in_test)
# A modal dialog left open would wait forever offscreen; fail instead.
set_tests_properties(thm_opt_in PROPERTIES
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen"
    TIMEOUT 120)

# Moving, adding or deleting a data segment line: the <thm> experiments and
# weights, the Experimental Effects lists and the segment_N rows of
# <parameterSettings> follow the line; a save round-trips.
add_executable(segments_renumber_test segments_renumber_test.cpp)
set_target_properties(segments_renumber_test PROPERTIES AUTOMOC OFF)
target_include_directories(segments_renumber_test PRIVATE
    ${AZURE2_SOURCE_DIR}/gui/include
    ${AZURE2_SOURCE_DIR}/include
)
target_compile_definitions(segments_renumber_test PRIVATE
    AZURE2_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
target_link_libraries(segments_renumber_test PRIVATE
    ${AZURESETUP_LIBRARIES}
    ${AZUREAPI_LIBRARIES}
    AZURE2Core
    Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network
)
if(USE_QWT)
    target_link_libraries(segments_renumber_test PRIVATE ${QWT_LIBRARIES})
endif()
if(USE_ERYA)
    target_link_libraries(segments_renumber_test PRIVATE pugixml srim_utilities)
endif()
add_test(NAME segments_renumber COMMAND segments_renumber_test)
set_tests_properties(segments_renumber PROPERTIES
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen"
    TIMEOUT 300)

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

# The MCMC tab's parameter list: each varying parameter's kind and automatic
# prior from its name, as the sampler takes them (the sqrt(E) shift
# coefficient has its own prior, not the energy shift's).
if(USE_MCMC)
    add_executable(mcmc_priors_test mcmc_priors_test.cpp)
    set_target_properties(mcmc_priors_test PROPERTIES AUTOMOC OFF)
    target_include_directories(mcmc_priors_test PRIVATE
        ${AZURE2_SOURCE_DIR}/gui/include
        ${AZURE2_SOURCE_DIR}/include
    )
    target_compile_definitions(mcmc_priors_test PRIVATE
        AZURE2_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        AZURE2_BINARY="$<TARGET_FILE:AZURE2>")
    target_link_libraries(mcmc_priors_test PRIVATE
        ${AZURESETUP_LIBRARIES}
        ${AZUREAPI_LIBRARIES}
        AZURE2Core
        Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network
    )
    if(USE_QWT)
        target_link_libraries(mcmc_priors_test PRIVATE ${QWT_LIBRARIES})
    endif()
    if(USE_ERYA)
        target_link_libraries(mcmc_priors_test PRIVATE pugixml srim_utilities)
    endif()
    add_dependencies(mcmc_priors_test AZURE2)
    add_test(NAME mcmc_priors COMMAND mcmc_priors_test)
    # A modal dialog on a failed open would wait forever offscreen; fail instead.
    set_tests_properties(mcmc_priors PROPERTIES
        ENVIRONMENT "QT_QPA_PLATFORM=offscreen"
        TIMEOUT 300)
endif()

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
