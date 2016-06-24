cmake_minimum_required(VERSION 2.8)
SET(CMAKE_BUILD_TYPE Release)

PROJECT(sphere_pack)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(${PROJECT_NAME} ${PROJECT_NAME} jsoncpp.cpp parse_json)
add_executable(sphere_resize sphere_resize)

if(VTK_LIBRARIES)
   target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})
   target_link_libraries(sphere_resize ${VTK_LIBRARIES})
endif()
