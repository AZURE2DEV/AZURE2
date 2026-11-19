# The table of physics functions in plain C++ (docs/upgrade/physics_tab.md section 2): the rows the binding
# (api/src/pybind_physics.cpp: pyazr.physics) and the front ends' facade (api/front) read.  No pybind11, no Qt, no
# engine header: it uses physics/, physics/transfer and experiment/ (Frames.h) only, the modules below it.
add_library(azure2_physics_table STATIC src/PhysicsTable.cpp)
target_include_directories(azure2_physics_table PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(azure2_physics_table PRIVATE azure2_physics azure2_physics_transfer azure2_experiment)
set_target_properties(azure2_physics_table PROPERTIES POSITION_INDEPENDENT_CODE ON)
