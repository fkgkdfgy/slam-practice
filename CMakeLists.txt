cmake_minimum_required(VERSION 2.8)
project(LidarSLAM2d)

set( CMAKE_CXX_COMPILER "g++" )
set(CMAKE_CXX_FLAGS "-std=c++11 -march=native -o3")
set(CMAKE_BUILD_TYPE "Release")

include_directories("/usr/include/eigen3")

find_package(Sophus REQUIRED)
include_directories(${Sophus_INCLUDE_DIRS})

find_package(OpenCV 3 REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

list( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules )
find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})
set(the_third_party_libs
        ${OpenCV_LIBS}
        ${Sophus_LIBRARIES}
        ${CERES_LIBRARIES}
        )

include_directories(include)
add_subdirectory(${PROJECT_SOURCE_DIR}/src)

