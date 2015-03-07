cmake_minimum_required(VERSION 2.8.6)

## We need C++11
set(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++")

## We need Eigen
find_package(Eigen3 REQUIRED)
include_directories( "${EIGEN3_INCLUDE_DIR}" )

## We need OpenSubdiv
set(OpenSubdiv_DIR "${CMAKE_SOURCE_DIR}/OpenSubdiv")

find_library(OPENSUBDIV_LIBRARY libosdCPU.a "${CMAKE_SOURCE_DIR}/OpenSubdiv/build/lib")
include_directories("${OpenSubdiv_DIR}/opensubdiv")

## To prevent duplicate compilation for the static and dynamic libraries
# http://www.cmake.org/Wiki/CMake/Tutorials/Object_Library
add_library(subdivision_skinning_common_objects
    OBJECT
    subdivision_matrices.cpp
    subdivision_matrices_eigen.cpp
    subdivision_matrices_osd.cpp
    subdivision_matrices_osd_eigen.cpp
    subdivision_skinning.cpp
    subdivision_engine.cpp
    subdivision_limit_mesh.cpp
    save_obj.cpp
    Timing.cpp
    )

## To prevent duplicate compilation for the static and dynamic libraries
# http://www.cmake.org/Wiki/CMake/Tutorials/Object_Library
add_library(subdivision_skinning_library_sources
    OBJECT
    subdivision_skinning_wrapper.cpp
    )
add_library(subdivision_skinning_library_static STATIC $<TARGET_OBJECTS:subdivision_skinning_library_sources> $<TARGET_OBJECTS:subdivision_skinning_common_objects>)
add_library(subdivision_skinning_library_shared SHARED $<TARGET_OBJECTS:subdivision_skinning_library_sources> $<TARGET_OBJECTS:subdivision_skinning_common_objects>)

set_target_properties(subdivision_skinning_library_static PROPERTIES OUTPUT_NAME subdivision_skinning CLEAN_DIRECT_OUTPUT 1)
set_target_properties(subdivision_skinning_library_shared PROPERTIES OUTPUT_NAME subdivision_skinning CLEAN_DIRECT_OUTPUT 1)
target_link_libraries(
    subdivision_skinning_library_shared
    "${OPENSUBDIV_LIBRARY}"
)
