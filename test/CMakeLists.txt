cmake_minimum_required (VERSION 3.4)

enable_testing ()
include (CTest)

## GITHUB UNIT TESTS
add_test (NAME small_test
          COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/small_test.sh" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

## DECODE INTERNAL UNIT TESTS
if (EXISTS "/nfs/odinn/users/hannesha/data/DONOTMOVEME/lrcaller-test")
    add_test (NAME test1
            COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test1.sh "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" "1")

    add_test (NAME test2a
            COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test1.sh "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" "2a")

    add_test (NAME test2b
            COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test1.sh "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" "2b")

    add_test (NAME test2c
            COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test1.sh "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" "2c")
endif ()
