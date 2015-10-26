cmake_minimum_required (VERSION 2.8)
project(PermTestingToolbox CXX)


# Armadillo needs to be compiled for C++11 in order to use large matrices
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

add_subdirectory(Utils)
add_subdirectory(ArmadilloUtils)
add_subdirectory(PermTestingCPU)

set(PermTestingMatrix_source_files PermTestingMatrix.cpp)
# Finda CUDA and OPENMP
#if(OPENMP_ENABLED)
#    find_package("OpenMP")
#    if (OPENMP_FOUND)
#        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    endif()
#endif()

add_executable(PermTestingMatrix ${PermTestingMatrix_source_files})
target_link_libraries(PermTestingMatrix ${ARMADILLO_LIBRARIES})








