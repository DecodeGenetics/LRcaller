cmake_minimum_required (VERSION 3.4)

enable_testing ()
include (CTest)

add_test (NAME small_test
          COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/small_test.sh" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

if(DEFINED ENV{DECODE_INTERNAL_TESTS})
    set(TESTDIR "$ENV{DECODE_INTERNAL_TESTS}/lrcaller-test")
    if (EXISTS ${TESTDIR})
        add_test (NAME test1
                  COMMAND ${TESTDIR}/test.sh "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" "1")

        add_test (NAME test2a
                  COMMAND ${TESTDIR}/test.sh "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" "2a")

        add_test (NAME test2b
                  COMMAND ${TESTDIR}/test.sh "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" "2b")

        add_test (NAME test2c
                  COMMAND ${TESTDIR}/test.sh "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" "2c")
    endif ()
endif ()



