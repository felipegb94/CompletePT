cmake_minimum_required (VERSION 2.8)
project(PermTestingToolbox CXX)

set(GPU_ENABLED 1)
set(OPENMP_ENABLED 0)

# Armadillo needs to be compiled for C++11 in order to use large matrices
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ")
find_package(Armadillo REQUIRED)
find_package(ArrayFire)
find_package(CUDA)
find_package(Boost)

if(Boost_FOUND)
    message("Boost was found")
endif()
if(CUDA_FOUND)
    message("CUDA FOUND!")
endif()

include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
if(GPU_ENABLED)
    include_directories(${ArrayFire_INCLUDE_DIRS})
endif()

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
    if(${ArrayFire_CUDA_FOUND})
        message("ArrayFire CUDA found")
    endif()
    if(${ArrayFire_CUDA_LIBRARIES})
        message("ArrayFire CUDA libraries found!")
    endif()    
endif()

add_subdirectory(Utils)
add_subdirectory(ArmadilloUtils)

if(GPU_ENABLED)
    set(PermTesting_source_files PermTestingShared.cpp PermTestingCPU.cpp PermTestingGPU.cpp)
    add_library(PermTesting ${PermTesting_source_files})
    target_link_libraries(PermTesting ${ArrayFire_LIBRARIES})
else()
    set(PermTesting_source_files PermTestingShared.cpp PermTestingCPU.cpp)
    add_library(PermTesting ${PermTesting_source_files})
endif()
target_link_libraries(PermTesting ${ARMADILLO_LIBRARIES})
target_link_libraries(PermTesting ${Boost_LIBRARIES})
target_link_libraries(PermTesting utils)
target_link_libraries(PermTesting armadilloUtils)


# Finda CUDA and OPENMP
#if(OPENMP_ENABLED)
#    find_package("OpenMP")
#    if (OPENMP_FOUND)
#        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    endif()
#endif()



#add_executable(demo_TwoSamplePermTestingCPU demo_TwoSamplePermTestingCPU.cpp)
#target_link_libraries(demo_TwoSamplePermTestingCPU ${ARMADILLO_LIBRARIES})
#target_link_libraries(demo_TwoSamplePermTestingCPU PermTesting)
#if(GPU_ENABLED)
#    add_executable(demo_TwoSamplePermTestingGPU demo_TwoSamplePermTestingGPU.cpp)
#    target_link_libraries(demo_TwoSamplePermTestingGPU ${ARMADILLO_LIBRARIES})
#    target_link_libraries(demo_TwoSamplePermTestingGPU PermTesting)
#endif()

add_executable(demo_OneSamplePermTestingCPU demo_OneSamplePermTestingCPU.cpp)
target_link_libraries(demo_OneSamplePermTestingCPU ${ARMADILLO_LIBRARIES})
target_link_libraries(demo_OneSamplePermTestingCPU PermTesting)
if(GPU_ENABLED)
    add_executable(demo_OneSamplePermTestingGPU demo_OneSamplePermTestingGPU.cpp)
    target_link_libraries(demo_OneSamplePermTestingGPU ${ARMADILLO_LIBRARIES})
    target_link_libraries(demo_OneSamplePermTestingGPU PermTesting)
endif()




