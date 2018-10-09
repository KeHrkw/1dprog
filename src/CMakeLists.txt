cmake_minimum_required(VERSION 2.8)


message(STATUS "CMAKE_SOURCE_DIR = ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR = ${CMAKE_CURRENT_BINARY_DIR}")
### Project settings
project(SALMON Fortran C)

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
message(STATUS "MAKE_DIRECTORY : ${CMAKE_CURRENT_BINARY_DIR}/output")

set(CMAKE_VERBOSE_MAKEFILE None)
message(STATUS "CMAKE_VERBOSE_MAKEFILE = ${CMAKE_VERBOSE_MAKEFILE}")

set(CMAKE_Fortran_COMPILER ifort)

option(IF_LOW OFF)
message(STATUS "IF_LOW = ${IF_LOW}")
if(IF_LOW)
  set(CMAKE_Fortran_FLAGS_RELEASE "-fpp -O3 -parallel -openmp ")
else()
  set(CMAKE_Fortran_FLAGS_RELEASE "-fpp -O3 -parallel -qopenmp ")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "-fpp -CB -traceback -g")

message(STATUS "Fortran_Flag_Debug = ${CMAKE_Fortran_FLAGS_DEBUG}")
message(STATUS "Fortran_Flag_Release = ${CMAKE_Fortran_FLAGS_RELEASE}")

option(Debug OFF)
message(STATUS "Debug = ${Debug}")
if(Debug)
set(CMAKE_BUILD_TYPE Debug)
else()
set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

set(sources main.f90 FourieT.f90 Initialize.f90 Module.f90 gnuplot.f90 input.f90)
#set(sources_Length ${CMAKE_SOURCE_DIR}/Length/Current_Length_operation0.f90 ${CMAKE_SOURCE_DIR}/Length/Gauge_Length0.f90 ${CMAKE_SOURCE_DIR}/Length/zh_Length_operation.f90)
set(sources_Length  ${CMAKE_SOURCE_DIR}/Length/Gauge_Length.f90)
#set(sources_Coef  ${CMAKE_SOURCE_DIR}/Coef/Coef.f90)
set(sources_Velocity ${CMAKE_SOURCE_DIR}/Velocity/Current_Velocity_operation.f90 ${CMAKE_SOURCE_DIR}/Velocity/Gauge_Velocity.f90 ${CMAKE_SOURCE_DIR}/Velocity/zh_Velocity_operation.f90)
set(sources_base ${CMAKE_SOURCE_DIR}/base/base_plot.f90 ${CMAKE_SOURCE_DIR}/base/CG.f90 ${CMAKE_SOURCE_DIR}/base/Phase_matching.f90)

add_executable(myapp ${sources} ${sources_Length} ${sources_Velocity} ${sources_Coef} ${sources_base})
