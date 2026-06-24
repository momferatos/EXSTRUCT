# CMake build for the exstruct project.
#
# Target:
#   exstruct.exe   structure extraction from an HDF5 field (needs HDF5 Fortran)
#
# Quick start:
#   cmake -B build -DHDF5_ROOT=/path/to/hdf5
#   cmake --build build -j
#
# Common options:
#   -DCMAKE_Fortran_COMPILER=nvfortran   (or gfortran)
#   -DEXSTRUCT_OPENMP=ON                  enable OpenMP threading   [OFF]
#   -DCMAKE_BUILD_TYPE=Release|Debug|...                           [Release]
#
# Note: the HDF5 Fortran modules must have been built with the SAME Fortran
# compiler used here, otherwise hdf5.mod cannot be read. Point CMake at the
# right install with -DHDF5_ROOT=<prefix>.

cmake_minimum_required(VERSION 3.18)   # for the Fortran_PREPROCESS property
project(exstruct LANGUAGES Fortran)

# Default to an optimised build.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY
               STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

option(EXSTRUCT_OPENMP "Enable OpenMP threading" OFF)

if(EXSTRUCT_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

find_package(HDF5 REQUIRED COMPONENTS Fortran)

add_executable(exstruct exstruct.f90)
set_target_properties(exstruct PROPERTIES
  OUTPUT_NAME      exstruct.exe
  Fortran_PREPROCESS ON)               # exstruct.f90 uses #ifdef
target_include_directories(exstruct PRIVATE ${HDF5_Fortran_INCLUDE_DIRS})
target_link_libraries(exstruct PRIVATE ${HDF5_Fortran_LIBRARIES})
if(EXSTRUCT_OPENMP)
  target_link_libraries(exstruct PRIVATE OpenMP::OpenMP_Fortran)
  target_compile_definitions(exstruct PRIVATE _OPENMP_)
endif()

message(STATUS "exstruct: OpenMP=${EXSTRUCT_OPENMP} build=${CMAKE_BUILD_TYPE}")
