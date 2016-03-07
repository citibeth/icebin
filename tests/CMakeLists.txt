# ---------------------------------
# https://cmake.org/cmake/help/v3.0/module/FindGTest.html
# https://github.com/snikulov/google-test-examples/blob/master/CMakeLists.txt

enable_testing()
find_package(Gtest REQUIRED)
include_directories(${GTEST_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/slib)

message("GTEST_LIBRARY " ${GTEST_LIBRARY})

SET(ALL_LIBS icebin ${EXTERNAL_LIBS} ${GTEST_LIBRARY})


foreach(TEST grid)
	add_executable(test_${TEST} test_${TEST}.cpp)
	target_link_libraries(test_${TEST} ${ALL_LIBS})
	add_test(AllTests test_${TEST})
endforeach()
