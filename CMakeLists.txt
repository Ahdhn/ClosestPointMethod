cmake_minimum_required(VERSION 3.11)

project(cpm-base VERSION 0.1.0 LANGUAGES CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
 
find_package(Eigen3 3.3 REQUIRED NO_MODULE) 

add_subdirectory(deps)
add_subdirectory(src)
add_subdirectory(examples)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)   
        