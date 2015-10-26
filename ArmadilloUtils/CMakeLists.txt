
# Armadillo needs to be compiled for C++11 in order to use large matrices
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

set(armadilloUtils_source_files
    PermTestingArmadilloUtils.cpp
    )

add_library(armadilloUtils ${armadilloUtils_source_files})
target_link_libraries(armadilloUtils ${ARMADILLO_LIBRARIES})


