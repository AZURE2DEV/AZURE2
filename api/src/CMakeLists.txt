include_directories(../include ${AZURE2_SOURCE_DIR}/include)

set (AZUREAPI_SOURCES main.cpp AZUREAPI.cpp AZURESocket.cpp)

set (AZUREAPI_HEADERS ../include/AZUREAPI.h
../include/AZURESocket.h)

add_library(AZUREAPI ${AZUREAPI_SOURCES} ${AZUREAPI_HEADERS_MOC})

# Winsock. AZURESocket.h carries #pragma comment(lib, "ws2_32.lib"), which only
# MSVC honours -- MinGW silently ignores it, so every socket call came out as an
# undefined __imp_* reference at link time. PUBLIC so anything linking AZUREAPI
# inherits it.
if(WIN32)
    target_link_libraries(AZUREAPI PUBLIC ws2_32)
endif()

# SuperSubatomic - added this
#target_link_libraries (AZURESetup Qt5::Core Qt5::PrintSupport)
#add_executable(AZURE2 AZURESetup.qrc)

set (TEMP_VAR AZUREAPI )
set (AZUREAPI_LIBRARIES ${TEMP_VAR} PARENT_SCOPE)