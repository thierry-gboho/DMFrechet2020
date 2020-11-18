cmake_minimum_required (VERSION 2.6)

project(frechet)

set(CMAKE_CXX_STANDARD 11)

message(STATUS "Checking OpenMP")
find_package(OpenMP)
# add flags for OpenMP
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Wall")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${OpenMP_SHARED_LINKER_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

set(CMAKE_CONFIGURATION_TYPES "Debug" CACHE STRING "" FORCE)


include_directories(
        src/includes
)

add_executable(
        frechet
        src/std/frechet.cpp
)
add_executable(
        frechetOmp
        src/omp/frechet.cpp
)