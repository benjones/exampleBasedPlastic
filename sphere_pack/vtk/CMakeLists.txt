cmake_minimum_required(VERSION 2.8)
SET(CMAKE_BUILD_TYPE Release)

PROJECT(signed_dist)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(signed_dist signed_dist.cpp)

if(VTK_LIBRARIES)
   target_link_libraries(signed_dist ${VTK_LIBRARIES})
endif()
