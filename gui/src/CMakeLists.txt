
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
    
    find_package (Qt5 COMPONENTS Core PrintSupport Gui Widgets REQUIRED)
    set(QT_LIBRARIES Qt5::Core Qt5::Gui Qt5::PrintSupport Qt5::Widgets)
    
    if(USE_QWT)
        # Try to find Qwt using standard CMake find_package
        find_package(Qwt QUIET)
        
        if(NOT Qwt_FOUND)
            message(STATUS "Standard Qwt package not found, searching system locations...")
            
            # Define common search paths
            set(QWT_SEARCH_PATHS
                # Conda environments
                "$ENV{CONDA_PREFIX}"
                "$ENV{CONDA_PREFIX}/lib"
                "$ENV{HOME}/miniconda3/envs/main"
                "$ENV{HOME}/miniforge3/envs/main"
                # Standard system paths
                "/usr/local"
                "/opt/homebrew"
                "/usr"
                # Additional common paths
                "/opt/local"
                "$ENV{HOME}/.local"
            )
            
            # Search for QWT headers
            find_path(QWT_INCLUDE_DIR
                NAMES qwt.h qwt_plot.h
                PATHS ${QWT_SEARCH_PATHS}
                PATH_SUFFIXES include include/qwt qwt
                NO_DEFAULT_PATH
            )
            
            # Search for QWT library
            find_library(QWT_LIBRARY
                NAMES qwt qwt6 libqwt libqwt6
                PATHS ${QWT_SEARCH_PATHS}
                PATH_SUFFIXES lib lib64
                NO_DEFAULT_PATH
            )
            
            # Check if we found both header and library
            if(QWT_INCLUDE_DIR AND QWT_LIBRARY)
                message(STATUS "Found QWT manually:")
                message(STATUS "  Headers: ${QWT_INCLUDE_DIR}")
                message(STATUS "  Library: ${QWT_LIBRARY}")
                
                # Set up the QWT target manually
                add_library(Qwt::Qwt UNKNOWN IMPORTED)
                set_target_properties(Qwt::Qwt PROPERTIES
                    IMPORTED_LOCATION "${QWT_LIBRARY}"
                    INTERFACE_INCLUDE_DIRECTORIES "${QWT_INCLUDE_DIR}"
                )
                
                set(QWT_LIBRARIES Qwt::Qwt)
                set(Qwt_FOUND TRUE)
                
                # Add include directory
                include_directories(${QWT_INCLUDE_DIR})
            else()
                message(FATAL_ERROR "Could not find QWT library. Please install QWT or set USE_QWT=OFF")
            endif()
        else()
            message(STATUS "Found QWT using find_package")
            set(QWT_LIBRARIES ${QWT_LIBRARIES})
        endif()
        
        add_definitions(-DUSE_QWT)
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
)

set(AZURESETUP_RCCS AZURESetup.qrc)

if(USE_QWT)
    list(APPEND AZURESETUP_SOURCES PlotTab.cpp AZUREPlot.cpp)
    list(APPEND AZURESETUP_HEADERS ../include/PlotTab.h ../include/AZUREPlot.h)
endif()

# Add GUI library target
add_library(AZURESetup
    ${AZURESETUP_SOURCES}
    ${AZURESETUP_HEADERS}
    ${AZURESETUP_RCCS}
)

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

# Link Qt and Qwt directly
target_link_libraries(AZURESetup PRIVATE ${QT_LIBRARIES})

if(USE_QWT)
    target_link_libraries(AZURESetup PRIVATE ${QWT_LIBRARIES})
endif()

# Export to parent
set(AZURESETUP_LIBRARIES AZURESetup PARENT_SCOPE)
