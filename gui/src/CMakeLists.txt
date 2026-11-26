# The Qt setup program AZURESetup (docs/upgrade/qt_gui.md): Qt 6 Widgets on the engine's facade azure2::front
# (api/front/), the one engine header it includes.  Built with -DBUILD_GUI=ON (off by default); Qt 6 Core, Gui,
# Widgets, Network and Svg (the plot's SVG export).
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets Network Svg)

# The headers with Q_OBJECT, so that AUTOMOC finds them (they sit in ../include, not beside the sources).
file(GLOB AZURESETUP_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/../include/*.h)

set(AZURESETUP_SOURCES
    AZURESetup.cpp AZURESetupFile.cpp main.cpp
    ProjectHolder.cpp RecordsPanel.cpp RecordDialog.cpp ValueEditor.cpp RecordTableModel.cpp SectionTextPanel.cpp RecordTabs.cpp FittingTab.cpp RunTab.cpp RunWorker.cpp
    PlotTab.cpp PlotWidget.cpp MCMCTab.cpp ThmTab.cpp PythonTools.cpp
    ConfigDialogs.cpp DiagnosticsPanel.cpp
    PairsModel.cpp PairsTab.cpp AddPairDialog.cpp PairKeyCombo.cpp MassTable.cpp MassCheckDialog.cpp
    LevelsModel.cpp LevelsTab.cpp LevelsFile.cpp ChannelsModel.cpp ChannelDetails.cpp AddLevelDialog.cpp
    IAEALevelData.cpp IAEALevelDialog.cpp
    LevelsTabEdits.cpp LevelsTabNavigation.cpp LevelsUndo.cpp LevelsViews.cpp LevelsPartialWaves.cpp
    PotentialPanel.cpp
    ChooseFileButton.cpp RichTextDelegate.cpp InfoDialog.cpp InTabDocs.cpp AboutAZURE2Dialog.cpp
    AzureTheme.cpp FileDialogPolicy.cpp StartPage.cpp NotesBar.cpp CoulombTableUi.cpp ConvertRepairDialog.cpp RecordExtras.cpp
    PhysicsTab.cpp PhysicsSessionPages.cpp PhysicsForms.cpp PhysicsWorker.cpp HeatmapWidget.cpp
    AZURESetup.qrc ShellIcons.qrc)

# The Overview page and the watch of the project's files (docs/upgrade/qt_gui.md sections 15 and 16).
list(APPEND AZURESETUP_SOURCES OverviewPage.cpp OverviewCharts.cpp ProjectWatcher.cpp)

# The Plot tab is drawn with QPainter (PlotWidget.cpp) and exported with QtSvg and QPdfWriter: no plotting library.
# dev's Qwt plot (USE_QWT) is not kept: it read the output files itself (docs/upgrade/qt_gui.md section 5).

# The program's library and the window, without main(), so that the headless tests (tests/gui) link the same code.
list(REMOVE_ITEM AZURESETUP_SOURCES main.cpp)
add_library(AZURESetupLib STATIC ${AZURESETUP_SOURCES} ${AZURESETUP_HEADERS})
target_include_directories(AZURESetupLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(AZURESetupLib PUBLIC Qt6::Core Qt6::Gui Qt6::Widgets Qt6::Network Qt6::Svg azure2_front)
# The start page lists the example projects of the checkout the program was built from (StartPage::examplePlaces).
target_compile_definitions(AZURESetupLib PRIVATE AZURE2_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
# A Qt built by a newer compiler than the one building AZURE2 (conda's Qt with an older system g++) needs its own
# libstdc++ and libgcc_s at link and run time: e.g. -DGUI_LINK_OPTIONS="-L$CONDA_PREFIX/lib;-Wl,-rpath,$CONDA_PREFIX/lib".
set(GUI_LINK_OPTIONS "" CACHE STRING "Extra link options of the GUI's programs (a list)")
if(GUI_LINK_OPTIONS)
  target_link_options(AZURESetupLib INTERFACE ${GUI_LINK_OPTIONS})
endif()
if(APPLE)
  target_compile_definitions(AZURESetupLib PUBLIC MACX_SPACING)
elseif(WIN32)
  target_compile_definitions(AZURESetupLib PUBLIC WIN_SPACING)
endif()

add_executable(AZURESetup main.cpp)
target_link_libraries(AZURESetup PRIVATE AZURESetupLib)

# Installed, the start page finds share/azure2/examples/index.json above the program (StartPage::examplePlaces) and
# lists what the manifest names, at the same paths relative to share/azure2 as in the checkout (qt_gui.md section
# 19.9).  The program, the manifest, examples/ and every other project folder the manifest names (the tutorials of
# tests/upgrade) are installed; the run products of a checkout (output/, checks/ contents, run.log) are not.
include(GNUInstallDirs)
install(TARGETS AZURESetup RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT gui)
set(AZURE2_EXAMPLE_DIRS examples)
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.19)
  file(READ ${CMAKE_SOURCE_DIR}/examples/index.json _manifest)
  string(JSON _groups LENGTH "${_manifest}" groups)
  math(EXPR _glast "${_groups} - 1")
  foreach(_g RANGE ${_glast})
    string(JSON _entries LENGTH "${_manifest}" groups ${_g} entries)
    math(EXPR _elast "${_entries} - 1")
    foreach(_e RANGE ${_elast})
      string(JSON _project GET "${_manifest}" groups ${_g} entries ${_e} project)
      get_filename_component(_dir "${_project}" DIRECTORY)
      list(APPEND AZURE2_EXAMPLE_DIRS ${_dir})
    endforeach()
  endforeach()
  list(REMOVE_DUPLICATES AZURE2_EXAMPLE_DIRS)
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/examples/index.json)
endif()
foreach(_dir ${AZURE2_EXAMPLE_DIRS})
  string(FIND "${_dir}" "examples/" _inExamples)
  if(_inExamples EQUAL 0)
    continue()  # inside examples/, installed with it
  endif()
  get_filename_component(_parent "${_dir}" DIRECTORY)
  install(DIRECTORY ${CMAKE_SOURCE_DIR}/${_dir} DESTINATION ${CMAKE_INSTALL_DATADIR}/azure2/${_parent} COMPONENT examples
          REGEX "/(output|checks)/.+" EXCLUDE REGEX "/run\\.log$" EXCLUDE REGEX "/check\\.log$" EXCLUDE
          REGEX "/__pycache__" EXCLUDE)
endforeach()
