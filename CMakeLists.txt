cmake_minimum_required(VERSION 2.8)
project(bnet)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(gsl REQUIRED)

if (UNIX)
  find_library(MATH_LIBRARY m)
endif (UNIX)

add_subdirectory(src)


