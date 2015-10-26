cmake_minimum_required (VERSION 2.8)
project(PermTestingToolbox CXX)


# Armadillo needs to be compiled for C++11 in order to use large matrices
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
find_package(Armadillo REQUIRED)


include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(Utils)
add_subdirectory(ArmadilloUtils)

set(permTesting_source_files
    PermTestingCPU.cpp
    )
add_library(permTesting ${permTesting_source_files})
target_link_libraries(permTesting ${ARMADILLO_LIBRARIES})
target_link_libraries(permTesting utils)
target_link_libraries(permTesting armadilloUtils)


# Finda CUDA and OPENMP
#if(OPENMP_ENABLED)
#    find_package("OpenMP")
#    if (OPENMP_FOUND)
#        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    endif()
#endif()

set(PermTestingMatrix_source_files PermTestingMatrix.cpp)
add_executable(PermTestingMatrix ${PermTestingMatrix_source_files})
target_link_libraries(PermTestingMatrix ${ARMADILLO_LIBRARIES})
target_link_libraries(PermTestingMatrix permTesting)




