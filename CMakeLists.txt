cmake_minimum_required (VERSION 3.4)

project (lrcaller CXX)

## DEPENDENCIES

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/submodules/seqan/include/seqan/version.h")
    set (CMAKE_INCLUDE_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/submodules/seqan/include
        ${CMAKE_INCLUDE_PATH})
    set (CMAKE_PREFIX_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/submodules/seqan/util/cmake
        ${CMAKE_PREFIX_PATH})
    set (CMAKE_MODULE_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/submodules/seqan/util/cmake
        ${CMAKE_MODULE_PATH})
    message (STATUS "Found a local SeqAn library provided with the Lambda source code.")
    message (    "   This will be preferred over system global headers.")
endif ()

find_package (ZLIB)
find_package (OpenMP)
find_package (SeqAn REQUIRED)

include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

## GENERAL SETTINGS

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable          (lrcaller lrcaller.cpp)
target_link_libraries   (lrcaller ${SEQAN_LIBRARIES})
target_compile_features (lrcaller PUBLIC cxx_std_20)
target_compile_options  (lrcaller PRIVATE -Wall -Wextra -pedantic -DSEQAN_DISABLE_VERSION_CHECK=1)

# Add tests
add_subdirectory(test)
enable_testing()
