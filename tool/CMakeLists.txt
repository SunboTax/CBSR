aux_source_directory(. code)
include_directories(include)
find_package(gflags)
set(SRC 
${CMAKE_SOURCE_DIR}/src/graph.cpp
${CMAKE_SOURCE_DIR}/src/utils.cpp
)

# foreach (_variableName ${SRC})
#     message(STATUS "${_variableName}=${${_variableName}}")
# endforeach()
add_executable(tool ${code} ${SRC})
target_link_libraries (tool gflags)