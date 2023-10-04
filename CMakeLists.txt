cmake_minimum_required(VERSION 3.2)

project (CMaths)

include(CTest)
enable_testing()


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
include(CPack)

add_subdirectory(src)

OPTION(TESTS "Enable testing of the project" ON)
if (TESTS)
    add_compile_definitions(COMPILE_WITH_TESTS)
    add_subdirectory(tests)
    message("TESTING ENABLED")
else()
    MESSAGE("TESTING DISABLED")
endif()

