cmake_minimum_required(VERSION 3.20)
project(CCR)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
IF(CMAKE_SYSTEM_NAME MATCHES "Linux")
    set(CMAKE_CXX_COMPILER "/usr/bin/g++")
    set(CMAKE_C_COMPILER "/usr/bin/gcc")
ELSEIF(CMAKE_SYSTEM_NAME MATCHES "Windows")
    set(CMAKE_CXX_COMPILER "D:/MinGW/mingw64/bin/g++.exe")
    set(CMAKE_C_COMPILER "D:/MinGW/mingw64/bin/gcc.exe")
ENDIF()
find_package(gflags)

add_definitions()

find_package(Threads REQUIRED)
SET(ROOT_DIR ${CMAKE_SOURCE_DIR})
SET(THIRD_DIR ${ROOT_DIR}/third)
SET(SPD_DIR ${THIRD_DIR}/spdlog/include)
SET(OREACH ${ROOT_DIR}/third/O-Reach)

include_directories(include ${SPD_DIR})
set(DCMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "-O3 -fopenmp")

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILE Computation.cpp)
file(GLOB_RECURSE SRC ${CMAKE_SOURCE_DIR}/src/*.cpp)


add_executable(Computation ${SOURCE_FILE} ${SRC})
target_link_libraries (Computation gflags Threads::Threads)
add_subdirectory(tool)



file(GLOB_RECURSE ALL_SOURCE_FILES 
     ${CMAKE_SOURCE_DIR}/*.c 
     ${CMAKE_SOURCE_DIR}/*.cpp 
     ${CMAKE_SOURCE_DIR}/*.h)

# 排除 third 目录
list(FILTER ALL_SOURCE_FILES EXCLUDE REGEX "${CMAKE_SOURCE_DIR}/third/.*")

# 自定义命令以运行 clang-format
add_custom_command(OUTPUT clang_format
                   COMMAND clang-format
                   ARGS -i -style=file ${ALL_SOURCE_FILES}
                   COMMENT "Running clang-format on all source files except those under 'third/' directory"
                   VERBATIM)

# 创建伪目标
add_custom_target(format DEPENDS clang_format)
