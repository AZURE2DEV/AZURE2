
IF(CROSS_COMPILE_WINDOWS)
    # Hardcode Qt and Qwt include/lib paths for MinGW
    set(MINGW_PREFIX /usr/x86_64-w64-mingw32/sys-root/mingw)

    # Qt includes and libraries
    include_directories(${MINGW_PREFIX}/include)
    include_directories(${MINGW_PREFIX}/include/qt5)
    include_directories(${MINGW_PREFIX}/include/qt5/QtCore)
    include_directories(${MINGW_PREFIX}/include/qt5/QtGui)
    include_directories(${MINGW_PREFIX}/include/qt5/QtWidgets)
    include_directories(${MINGW_PREFIX}/include/qt5/QtPrintSupport)
    include_directories(${MINGW_PREFIX}/include/qt5/QtNetwork)

    # Qt-specific macros and MOC support
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTORCC ON)
    set(CMAKE_AUTOUIC ON)

    # Link to Qt5 .dll.a import libs directly
    set(QT_LIBRARIES
        ${MINGW_PREFIX}/lib/libQt5Core.dll.a
        ${MINGW_PREFIX}/lib/libQt5Gui.dll.a
        ${MINGW_PREFIX}/lib/libQt5Widgets.dll.a
        ${MINGW_PREFIX}/lib/libQt5PrintSupport.dll.a
        ${MINGW_PREFIX}/lib/libQt5Network.dll.a
    )

    if(USE_QWT)
        include_directories(${MINGW_PREFIX}/include/qt5/qwt)
        set(QWT_LIBRARIES ${MINGW_PREFIX}/lib/libqwt-qt5.dll.a)
        add_definitions(-DUSE_QWT)
    endif()
else()
    # Qt-specific macros and MOC support for native builds
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTORCC ON)
    set(CMAKE_AUTOUIC ON)
     
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")
    find_package (Qt5 COMPONENTS Core PrintSupport Gui Widgets Network REQUIRED)
    set(QT_LIBRARIES Qt5::Core Qt5::Gui Qt5::PrintSupport Qt5::Widgets Qt5::Network)
    
    # Only try to find QWT if USE_QWT is enabled
    if(USE_QWT)
        # Try to find Qwt using standard CMake find_package
        find_package(Qwt QUIET)
        if(NOT Qwt_FOUND) 
            message(FATAL_ERROR "Could not find QWT library. Please install libqwt-qt5-dev or set USE_QWT=OFF") 
        else()
            message(STATUS "Found QWT using find_package") 
            # Since the custom FindQwt.cmake module finds the library and include directories,
            # we need to use those variables to set up our build.
            # The find_package call will set QWT_LIBRARY and QWT_INCLUDE_DIR.
            set(QWT_LIBRARIES ${QWT_LIBRARY}) 
            include_directories(${QWT_INCLUDE_DIR}) 
            add_definitions(-DUSE_QWT) 
        endif()
    endif()
endif()

# Source and header files
set(AZURESETUP_SOURCES
    ChooseFileButton.cpp PairsModel.cpp AddLevelDialog.cpp EditChecksDialog.cpp
    PairsTab.cpp AddPairDialog.cpp EditDirsDialog.cpp RichTextDelegate.cpp
    AddSegDataDialog.cpp EditOptionsDialog.cpp RunTab.cpp AddSegTestDialog.cpp
    SegmentsDataModel.cpp AddTargetIntDialog.cpp SegmentsTab.cpp AZURESetup.cpp
    LevelsModel.cpp SegmentsTestModel.cpp ChannelDetails.cpp LevelsTab.cpp
    TargetIntModel.cpp ChannelsModel.cpp main.cpp TargetIntTab.cpp FittingTab.cpp
    AboutAZURE2Dialog.cpp InfoDialog.cpp InTabDocs.cpp
    IAEALevelData.cpp IAEALevelDialog.cpp NuclearPotentialTab.cpp
    ExforData.cpp ExforDialog.cpp
)

set(AZURESETUP_HEADERS
    ../include/AddLevelDialog.h ../include/AddPairDialog.h ../include/AddSegDataDialog.h
    ../include/AddSegTestDialog.h ../include/AddTargetIntDialog.h ../include/AZUREMainThread.h
    ../include/AZURESetup.h ../include/ChannelDetails.h ../include/ChannelsModel.h
    ../include/ChooseFileButton.h ../include/EditChecksDialog.h ../include/EditDirsDialog.h
    ../include/EditOptionsDialog.h ../include/ElementMap.h ../include/FilteredTextEdit.h
    ../include/LevelsModel.h ../include/LevelsTab.h ../include/PairsModel.h
    ../include/PairsTab.h ../include/RichTextDelegate.h ../include/RunTab.h
    ../include/SegmentsDataModel.h ../include/SegmentsTab.h ../include/SegmentsTestModel.h
    ../include/TargetIntModel.h ../include/TargetIntTab.h ../include/FittingTab.h
    ../include/TextEditBuffer.h ../include/AboutAZURE2Dialog.h ../include/InfoDialog.h
    ../include/IAEALevelData.h ../include/IAEALevelDialog.h ../include/NuclearPotentialTab.h
    ../include/ExforData.h ../include/ExforDialog.h
)

set(AZURESETUP_RCCS AZURESetup.qrc)

if(USE_QWT)
    list(APPEND AZURESETUP_SOURCES PlotTab.cpp AZUREPlot.cpp)
    list(APPEND AZURESETUP_HEADERS ../include/PlotTab.h ../include/AZUREPlot.h)
endif()

if(USE_MCMC)
    list(APPEND AZURESETUP_SOURCES MCMCTab.cpp AZUREMCMCThread.cpp)
    list(APPEND AZURESETUP_HEADERS ../include/MCMCTab.h ../include/AZUREMCMCThread.h)
endif()

# Add GUI library target
add_library(AZURESetup
    ${AZURESETUP_SOURCES}
    ${AZURESETUP_HEADERS}
    ${AZURESETUP_RCCS}
)

# Ensure AUTOMOC can find all necessary headers
if(USE_MCMC)
    set_target_properties(AZURESetup PROPERTIES AUTOMOC_MOC_OPTIONS "-I${CMAKE_CURRENT_SOURCE_DIR}/../include")
endif()

# Include your own project headers
target_include_directories(AZURESetup PRIVATE
    ../include
    ${AZURE2_SOURCE_DIR}/include
)

# Platform definitions
if(APPLE)
    target_compile_definitions(AZURESetup PRIVATE -DMACX_SPACING)
elseif(WIN32)
    target_compile_definitions(AZURESetup PRIVATE -DWIN_SPACING)
endif()

# Link Qt, and optionally Qwt and ERYA utilities
target_link_libraries(AZURESetup PRIVATE ${QT_LIBRARIES})

if(USE_QWT)
    target_link_libraries(AZURESetup PRIVATE ${QWT_LIBRARIES})
endif()

if(USE_ERYA)
    target_link_libraries(AZURESetup PRIVATE pugixml srim_utilities)
endif()

# Export to parent
set(AZURESETUP_LIBRARIES AZURESetup PARENT_SCOPE)
# The Qt (and Qwt) link lists are resolved here -- imported targets on a
# native build, import libraries when cross-compiling.  tests/gui needs the
# same ones, so hand them up rather than repeat the detection.
set(QT_LIBRARIES ${QT_LIBRARIES} PARENT_SCOPE)
set(QWT_LIBRARIES ${QWT_LIBRARIES} PARENT_SCOPE)
