cmake_minimum_required(VERSION 3.15)
project(SmoothLineSegment)

set(CMAKE_CXX_STANDARD 14)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

find_package(Glog REQUIRED)

include_directories(${Glog_INCLUDE_DIRS})

add_executable(SmoothLineSegment main.cpp)
target_link_libraries(SmoothLineSegment ${GLOG_LIBRARY})