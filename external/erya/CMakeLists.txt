# ERYA SRIM Utilities CMakeLists.txt
# This is a subdirectory of the main AZURE2 project

# Add subdirectory for pugixml
add_subdirectory(pugixml)

# Create SRIM utilities library
add_library(srim_utilities STATIC
    src/SRIMUtilities.cpp
)

# Set include directories for the library
target_include_directories(srim_utilities PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/pugixml/src
)

# Link with pugixml
target_link_libraries(srim_utilities PRIVATE pugixml)

# Make the library target available to parent scope
set(SRIM_UTILITIES_LIBRARIES srim_utilities PARENT_SCOPE)

# Export include directories to parent scope
set(SRIM_UTILITIES_INCLUDE_DIRS 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/pugixml/src
    PARENT_SCOPE)

# Install targets (only if this is the main project, not a subdirectory)
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    install(FILES include/SRIMUtilities.h DESTINATION include/erya)
    install(TARGETS srim_utilities DESTINATION lib)
    install(FILES data/SRIM2013.xml DESTINATION share/azure2/)
endif()