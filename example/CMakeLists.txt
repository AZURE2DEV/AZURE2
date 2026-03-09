cmake_minimum_required (VERSION 2.8)
project (AZURE2)

if (MINUIT_PATH) 
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${MINUIT_PATH})
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${MINUIT_PATH})
endif (MINUIT_PATH) 

if (GSL_PATH)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${GSL_PATH})
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${GSL_PATH})
endif (GSL_PATH) 


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package (gsl REQUIRED)
include_directories(${GSL_INCLUDE_DIR})
find_package (Minuit2 REQUIRED)
include_directories(${MINUIT2_INCLUDE_DIR})
find_package (OpenMP REQUIRED)
set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS})

option (BUILD_GUI 
       "Build AZURE2 with graphical setup program" ON)
option (USE_STAT
          "Use stat() function to check for directories.  Not
          recommended for windows builds" YES)
option (USE_QWT
          "Use Qwt to build Plot tab." ON)

if (BUILD_GUI)
add_subdirectory (gui)
add_definitions(-DGUI_BUILD)
endif (BUILD_GUI)

include_directories(${AZURE2_SOURCE_DIR}/coul/include)
add_subdirectory (coul)
set (EXTRA_LIBS coul)

if (NOT USE_STAT) 
add_definitions (-DNO_STAT) 
endif (NOT USE_STAT)

add_subdirectory (src)

