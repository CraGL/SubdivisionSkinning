cmake_minimum_required(VERSION 2.8.6)

## The path to the FindEigen file.
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

## Build the library.
add_subdirectory(lib)

## Build the GUI
add_subdirectory(subdivgui)
