# ---------------------------------
# https://cmake.org/cmake/help/v3.0/module/FindGTest.html
# https://github.com/snikulov/google-test-examples/blob/master/CMakeLists.txt

enable_testing()
find_package(Gtest REQUIRED)
include_directories(${GTEST_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/slib)

message("GTEST_LIBRARY " ${GTEST_LIBRARY})

SET(ALL_LIBS icebin ${EXTERNAL_LIBS} ${GTEST_LIBRARY})


foreach(TEST grid z1qx1n_bs1)
    add_executable(test_${TEST} test_${TEST}.cpp)
    target_link_libraries(test_${TEST} ${ALL_LIBS})
    add_test(AllTests test_${TEST})
endforeach()


# This test has a second Fortran file in it
add_executable(test_hntr test_hntr.cpp Z1QX1N.BS1.F help_hntr.F90)
target_link_libraries(test_hntr ${ALL_LIBS} gfortran)   # Linking in gfortran is GCC-only
add_test(AllTests test_hntr)
