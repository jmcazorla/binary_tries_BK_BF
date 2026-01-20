cmake_minimum_required(VERSION 3.5)

project(COMP_BINARY_TRIES_BK_BF VERSION 1.0)

include_directories(~/include)
# include_directories(${COMP_BINARY_TRIES_BK_BF_SOURCE_DIR}/include)
link_directories(~/lib)

set(CMAKE_CXX_STANDARD 14)

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    add_compile_definitions("LINUX")
else()
    add_compile_definitions("MACOS")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

add_subdirectory(test)

