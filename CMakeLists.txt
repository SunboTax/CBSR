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
      ${CMAKE_SOURCE_DIR}/*.cc
     ${CMAKE_SOURCE_DIR}/*.h)

# 排除 third 目录
list(FILTER ALL_SOURCE_FILES EXCLUDE REGEX "${CMAKE_SOURCE_DIR}/third/.*")

if(CLANG_FORMAT)
    file(GLOB_RECURSE ALL_SOURCE_FILES
        "*.cpp" "*.h" "*.cc" 
    )
    list(FILTER ALL_SOURCE_FILES EXCLUDE REGEX ".*third.*")
    
    # 将命令和参数作为列表传递
    add_custom_target(clang-format
        COMMAND "${CLANG_FORMAT}" "-i" "-style=LLVM" ${ALL_SOURCE_FILES}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Running clang-format on all source files"
    )
endif()
