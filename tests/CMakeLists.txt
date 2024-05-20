# Add Google Test
add_subdirectory(googletest-main)

# Link to Google Test Library
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# Add test executables
add_executable(tests
#task1.cpp
#task2.cpp
#test_matrix.cpp
2ex_2task.cpp
#daywork.cpp
)


# Link test executable against Google Test and our library

#target_link_libraries(run_tests gtest gtest_main SLAE_methods)
target_link_libraries(tests gtest gtest_main comp_maths)
