
# Armadillo needs to be compiled for C++11 in order to use large matrices
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

set(permTestingCPU_source_files
    PermTestingCPU.cpp
    PermTestingCPU.h
    )

add_library(PermTestingCPU ${permTestingCPU_source_files})
target_link_libraries(PermTestingCPU ${ARMADILLO_LIBRARIES})
