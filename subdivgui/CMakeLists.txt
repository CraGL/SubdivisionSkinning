cmake_minimum_required(VERSION 2.8.6)

## CGAL
find_package(CGAL REQUIRED)
## This conflicts with many of my directives if it isn't first.
include(${CGAL_USE_FILE})
# CGAL's monkeying with all of the flags. Rather than change the CGAL_USE_FILE
# just get ride of this flag.
# http://stackoverflow.com/a/18234926/148668
macro(remove_cxx_flag flag)
  string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endmacro()
remove_cxx_flag("-stdlib=libc++")

# This is absolutely necessary for Exact Construction
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -frounding-math -fsignaling-nans")
#set(SUBDIVGUI_LIBS ${SUBDIVGUI_LIBS} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

## We need Eigen
find_package(Eigen3 REQUIRED)
include_directories( "${EIGEN3_INCLUDE_DIR}" )

## We need libigl

if(DEFINED ENV{LIBIGL})
  set(LIBIGL $ENV{LIBIGL})
else(DEFINED ENV{LIBIGL})
  set(LIBIGL "${CMAKE_SOURCE_DIR}/libigl")
endif(DEFINED ENV{LIBIGL})
include_directories("${LIBIGL}/include")

## We don't have/want MOSEK
add_definitions(-DIGL_NO_MOSEK)

## AntTweakBar
include_directories("${LIBIGL}/external/AntTweakBar/include")
find_library(ANTTWEAKBAR_LIB AntTweakBar "${LIBIGL}/external/AntTweakBar/lib")
set(SUBDIVGUI_LIBS ${SUBDIVGUI_LIBS} ${ANTTWEAKBAR_LIB})

## Embree
include_directories("${LIBIGL}/external/embree")
include_directories("${LIBIGL}/external/embree/include")
find_library(EMBREE_MAIN_LIB NAMES embree PATHS "${LIBIGL}/external/embree/build")
find_library(EMBREE_SYS_LIB NAMES sys PATHS "${LIBIGL}/external/embree/build")
set(SUBDIVGUI_LIBS ${SUBDIVGUI_LIBS} ${EMBREE_MAIN_LIB} ${EMBREE_SYS_LIB})

## Tetgen
include_directories("${LIBIGL}/external/tetgen")
find_library(TETGEN_LIB tet "${LIBIGL}/external/tetgen")
set(SUBDIVGUI_LIBS ${SUBDIVGUI_LIBS} ${TETGEN_LIB})

## YImage
include_directories("${LIBIGL}/external/yimg")
find_library(YIMG_LIB yimg "${LIBIGL}/external/yimg")
find_package(PNG REQUIRED)
set(SUBDIVGUI_LIBS ${SUBDIVGUI_LIBS} ${YIMG_LIB} ${PNG_LIBRARIES})

## OpenGL/GLU/GLUT
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
#include_directories("${OPENGL_INCLUDE_DIRS}")
#include_directories("${GLUT_INCLUDE_DIRS}")
## I don't have Alec's patched GLUT with GLUT_ACTIVE_COMMAND
add_definitions(-DGLUT_ACTIVE_COMMAND=8)
set(SUBDIVGUI_LIBS ${SUBDIVGUI_LIBS} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

## Our library
include_directories("${CMAKE_SOURCE_DIR}/lib")
## I'm not good enough at CMake to use find_library() on a previous subdirectories output.
# find_library(SUBDIVISION_SKINNING_LIB subdivision_skinning "${CMAKE_BINARY_DIR}/lib")
set(SUBDIVISION_SKINNING_LIB "${CMAKE_BINARY_DIR}/lib/libsubdivision_skinning.a")
## OpenSubdiv library is defined by the CMakeList from our lib directory, which
## is processed first.
set(SUBDIVGUI_LIBS ${SUBDIVGUI_LIBS} ${SUBDIVISION_SKINNING_LIB} ${OPENSUBDIV_LIBRARY})

## We need C++11. Put this directive after CGAL's include.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(subdivgui
    subdivgui.cpp
    clean.cpp
    robust_bbw.cpp
    subdiv_weights.cpp
    Mesh.cpp
    )
target_link_libraries(subdivgui
    ${SUBDIVGUI_LIBS}
    )
