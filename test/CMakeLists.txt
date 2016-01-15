# ---------------------------------
# https://cmake.org/cmake/help/v3.0/module/FindGTest.html
# https://github.com/snikulov/google-test-examples/blob/master/CMakeLists.txt

enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/slib)


SET(ALL_LIBS ${GTEST_LIBRARY} ${EXTERNAL_LIBS} icebin)


foreach(TEST grid)
	add_executable(test_${TEST} test_${TEST}.cpp)
	target_link_libraries(test_${TEST} ${ALL_LIBS})
	add_test(AllTests test_${TEST})
endforeach()
