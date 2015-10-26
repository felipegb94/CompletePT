cmake_minimum_required (VERSION 2.8)

project(PermTesting_GPU)

# Armadillo needs to be compiled for C++11 in order to use large matrices
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
# Find ArrayFire and Armadillo
#find_package(OpenCL)
find_package(ArrayFire)
find_package(Armadillo REQUIRED)

#include_directories(${OpenCL_INCLUDE_DIRS})
include_directories(${ArrayFire_INCLUDE_DIRS})
include_directories(${Armadillo_INCLUDE_DIRS})

if(${Armadillo_FOUND})
    message("Armadillo was found!")
endif()
if(${ArrayFire_FOUND})
    message("ArrayFire was found!")
    if(${ArrayFire_CPU_FOUND})
        message("ArrayFire CPU libraries were found!")
    endif()
    if(${ArrayFire_OpenCL_FOUND})
        message("ArrayFire OpenCL libraries were found!")
    endif()
endif()


add_executable(PermTesting_AFMatrix PermTesting_AFMatrix.cpp)
target_link_libraries(PermTesting_AFMatrix ${ARMADILLO_LIBRARIES})
target_link_libraries(PermTesting_AFMatrix ${ArrayFire_LIBRARIES})



