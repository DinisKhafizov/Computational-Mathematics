# Add Google Test
add_subdirectory(googletest-main)

# Link to Google Test Library
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# Add test executables
add_executable(tests
task2.cpp
task3.cpp
task4.cpp
task5.cpp
interpolation.cpp) 


# Link test executable against Google Test and our library
#target_link_libraries(run_tests gtest gtest_main SLAE_methods)
target_link_libraries(tests gtest gtest_main comp_maths)
