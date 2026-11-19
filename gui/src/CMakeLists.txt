# The Qt setup program AZURESetup (docs/upgrade/qt_gui.md): Qt 6 Widgets on the engine's facade azure2::front
# (api/front/), the one engine header it includes.  Built with -DBUILD_GUI=ON (off by default); the Plot tab needs Qwt
# (-DUSE_QWT=ON) and is left out without it.
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets Network)

# The headers with Q_OBJECT, so that AUTOMOC finds them (they sit in ../include, not beside the sources).
file(GLOB AZURESETUP_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/../include/*.h)

set(AZURESETUP_SOURCES
    AZURESetup.cpp AZURESetupFile.cpp main.cpp
    ProjectHolder.cpp RecordsPanel.cpp SectionTextPanel.cpp RecordTabs.cpp FittingTab.cpp RunTab.cpp RunWorker.cpp
    ConfigDialogs.cpp DiagnosticsPanel.cpp
    PairsModel.cpp PairsTab.cpp AddPairDialog.cpp PairKeyCombo.cpp MassTable.cpp MassCheckDialog.cpp
    LevelsModel.cpp LevelsTab.cpp LevelsFile.cpp ChannelsModel.cpp ChannelDetails.cpp AddLevelDialog.cpp
    IAEALevelData.cpp IAEALevelDialog.cpp
    ChooseFileButton.cpp RichTextDelegate.cpp InfoDialog.cpp InTabDocs.cpp AboutAZURE2Dialog.cpp
    AzureTheme.cpp FileDialogPolicy.cpp
    AZURESetup.qrc)

if(USE_QWT)
  # The Plot tab of dev reads the output files with Qwt; not ported to the facade (docs/upgrade/qt_gui.md section 4).
  message(WARNING "USE_QWT: the Plot tab is not ported to this branch; the GUI is built without it.")
endif()

# The program's library and the window, without main(), so that the headless tests (tests/gui) link the same code.
list(REMOVE_ITEM AZURESETUP_SOURCES main.cpp)
add_library(AZURESetupLib STATIC ${AZURESETUP_SOURCES} ${AZURESETUP_HEADERS})
target_include_directories(AZURESetupLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(AZURESetupLib PUBLIC Qt6::Core Qt6::Gui Qt6::Widgets Qt6::Network azure2_front)
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
