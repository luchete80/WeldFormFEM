cmake_minimum_required(VERSION 3.7.0...${CMAKE_VERSION})

project(Omega_h VERSION 9.34.13 LANGUAGES CXX)

include_directories(../../lib/omega_h-9.34.13/src) #FOR: Omega_h_config.h

add_executable(remesh_test main.cpp)
target_link_libraries(remesh_test PRIVATE omega_h)
add_test(NAME remesh_test COMMAND $<TARGET_FILE:remesh_test>)

#add_custom_command(OUTPUT square.msh
#	COMMAND ${GMSH_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/square.geo
#	-2 -o square.msh
#	MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/square.geo
#)
#add_custom_target(osh_field_on_square_mesh ALL DEPENDS square.msh)
