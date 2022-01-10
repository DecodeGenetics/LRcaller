cmake_minimum_required (VERSION 3.4)

project (lrcaller CXX)

set (LRCALLER_VERSION "1.0")
add_definitions(-DLRCALLER_VERSION="${LRCALLER_VERSION}")

#============================================= DEPENDENCIES =================================================

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

find_package (ZLIB REQUIRED)
find_package (OpenMP REQUIRED)
find_package (SeqAn REQUIRED)

include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

#============================================= GENERAL SETTINGS =============================================

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable          (lrcaller lrcaller.cpp)
target_link_libraries   (lrcaller ${SEQAN_LIBRARIES})
target_compile_features (lrcaller PUBLIC cxx_std_20)
target_compile_options  (lrcaller PRIVATE -Wall -Wextra -pedantic -DSEQAN_DISABLE_VERSION_CHECK=1)

#============================================= TESTS ========================================================

## Add Unit tests
add_subdirectory(test)
enable_testing()

## Add clang-format
find_program(CLANG_FORMAT "clang-format-12")

if (CLANG_FORMAT STREQUAL "CLANG_FORMAT-NOTFOUND")
    find_program(CLANG_FORMAT "clang-format")
endif()

if (NOT CLANG_FORMAT STREQUAL "CLANG_FORMAT-NOTFOUND")
    add_custom_target (check_format "${CLANG_FORMAT}" "-style=file" "-n" "-Werror" "${CMAKE_CURRENT_SOURCE_DIR}/*.?pp"
                       COMMENT "Checking all files' format with clang-format.")
    add_custom_target (do_format "${CLANG_FORMAT}" "-style=file" "-i" "${CMAKE_CURRENT_SOURCE_DIR}/*.?pp"
                       COMMENT "Enforce all files' format with clang-format.")
endif()
