find_package (Qt4 4.4 REQUIRED)
IF( USE_QWT )
find_package (Qwt REQUIRED)
ENDIF( USE_QWT )

IF( USE_QWT )
include_directories(../include ${AZURE2_SOURCE_DIR}/include ${Qwt6_INCLUDE_DIR})
ELSE( USE_QWT )
include_directories(../include ${AZURE2_SOURCE_DIR}/include)
ENDIF( USE_QWT )

set (AZURESETUP_SOURCES ChooseFileButton.cpp PairsModel.cpp
AddLevelDialog.cpp EditChecksDialog.cpp PairsTab.cpp AddPairDialog.cpp
EditDirsDialog.cpp RichTextDelegate.cpp AddSegDataDialog.cpp
EditOptionsDialog.cpp RunTab.cpp AddSegTestDialog.cpp
SegmentsDataModel.cpp AddTargetIntDialog.cpp SegmentsTab.cpp
AZURESetup.cpp LevelsModel.cpp SegmentsTestModel.cpp
ChannelDetails.cpp LevelsTab.cpp TargetIntModel.cpp ChannelsModel.cpp
main.cpp TargetIntTab.cpp AboutAZURE2Dialog.cpp InfoDialog.cpp InTabDocs.cpp)

set (AZURESETUP_HEADERS ../include/AddLevelDialog.h
../include/AddPairDialog.h ../include/AddSegDataDialog.h
../include/AddSegTestDialog.h ../include/AddTargetIntDialog.h
../include/AZUREMainThread.h ../include/AZURESetup.h
../include/ChannelDetails.h ../include/ChannelsModel.h
../include/ChooseFileButton.h ../include/EditChecksDialog.h
../include/EditDirsDialog.h ../include/EditOptionsDialog.h
../include/ElementMap.h ../include/FilteredTextEdit.h
../include/LevelsModel.h ../include/LevelsTab.h
../include/PairsModel.h ../include/PairsTab.h
../include/RichTextDelegate.h ../include/RunTab.h
../include/SegmentsDataModel.h ../include/SegmentsTab.h
../include/SegmentsTestModel.h ../include/TargetIntModel.h
../include/TargetIntTab.h ../include/TextEditBuffer.h
../include/AboutAZURE2Dialog.h ../include/InfoDialog.h)

set (AZURESETUP_RCCS AZURESetup.qrc)

IF( USE_QWT )
set (AZURESETUP_SOURCES ${AZURESETUP_SOURCES} PlotTab.cpp AZUREPlot.cpp)
set (AZURESETUP_HEADERS ${AZURESETUP_HEADERS} ../include/PlotTab.h ../include/AZUREPlot.h)
add_definitions(-DUSE_QWT)
ENDIF( USE_QWT )

include (${QT_USE_FILE})
add_definitions(${QT_DEFINITIONS})

if (APPLE)
add_definitions(-DMACX_SPACING)
endif (APPLE)
if (WIN32)
add_definitions(-DWIN_SPACING)
endif (WIN32)

QT4_ADD_RESOURCES(AZURESETUP_SOURCES ${AZURESETUP_RCCS})
QT4_WRAP_CPP(AZURESETUP_HEADERS_MOC ${AZURESETUP_HEADERS})

add_library(AZURESetup ${AZURESETUP_SOURCES} ${AZURESETUP_HEADERS_MOC})

set (TEMP_VAR AZURESetup )
IF( USE_QWT )
FIND_LIBRARY( QtSvg_LIBRARY NAMES QtSvg PATHS ${QT_LIBRARY_DIR} )
set (TEMP_VAR ${TEMP_VAR} ${Qwt6_Qt4_LIBRARY} ${QtSvg_LIBRARY} )
ENDIF( USE_QWT )
set (TEMP_VAR ${TEMP_VAR} ${QT_LIBRARIES})

if(NOT QT_GOBJECT_LIBRARY)
find_library(QT_GOBJECT_LIBRARY NAMES gobject-2.0)
endif (NOT QT_GOBJECT_LIBRARY)
if (QT_GOBJECT_LIBRARY)
set (TEMP_VAR ${TEMP_VAR} ${QT_GOBJECT_LIBRARY})
mark_as_advanced(QT_GOBJECT_LIBRARY)
endif (QT_GOBJECT_LIBRARY)
set (AZURESETUP_LIBRARIES ${TEMP_VAR} PARENT_SCOPE)
