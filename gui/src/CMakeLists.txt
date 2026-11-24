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
    LevelsTabEdits.cpp LevelsTabNavigation.cpp LevelsUndo.cpp LevelsViews.cpp
    ChooseFileButton.cpp RichTextDelegate.cpp InfoDialog.cpp InTabDocs.cpp AboutAZURE2Dialog.cpp
    AzureTheme.cpp FileDialogPolicy.cpp StartPage.cpp NotesBar.cpp ConvertRepairDialog.cpp
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
