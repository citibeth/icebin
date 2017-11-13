find_package(Python)
find_package(Cython)
find_package(Numpy)


# ----------Stuff to help with running Python before this project is installed.

## --- Convenience link to Python source
#execute_process(COMMAND ln -s ${CMAKE_CURRENT_SOURCE_DIR} src
#    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# -------------------------------------------------------------

include_directories(${PYTHON_INCLUDES} ${NUMPY_INCLUDE_DIRS})
list(APPEND EXTERNAL_LIBS ${PYTHON_LIBRARY})
# Compile everything (including Cython extension) as a single
# "compilation unit" for Numpy.  See:
# http://stackoverflow.com/questions/12253389/numpy-c-api-link-several-object-files
add_definitions(-DPY_ARRAY_UNIQUE_SYMBOL=icebin_array_unique_symbol)

# Only ONE compilation unit should NOT have -DNO_IMPORT_ARRAY
foreach (SOURCE_FILE icebin_cython.cpp ibmisc_cython.cpp )
    list(APPEND PYLIB_SOURCES ${SOURCE_FILE})
    set_property(SOURCE ${SOURCE_FILE}
        PROPERTY COMPILE_DEFINITIONS NO_IMPORT_ARRAY )
endforeach()

# Construct Cython -I flags
get_directory_property(CYTHON_INCLUDE_DIRECTORIES INCLUDE_DIRECTORIES)
foreach (DIR ${CYTHON_INCLUDE_DIRECTORIES})
    list(APPEND CYTHON_INCLUDE_FLAGS -I${DIR})
endforeach()

foreach (CYTHON_MODULE _icebin)

    add_custom_command(OUTPUT ${CYTHON_MODULE}.cpp
        MAIN_DEPENDENCY ${PROJECT_SOURCE_DIR}/pylib/${CYTHON_MODULE}.pyx
        DEPENDS cicebin.pxd
        COMMAND ${CYTHON_EXECUTABLE} --cplus -3
            --include-dir ${PROJECT_SOURCE_DIR}/pylib
            ${CYTHON_INCLUDE_FLAGS}
            ${PROJECT_SOURCE_DIR}/pylib/${CYTHON_MODULE}.pyx
            -o ${CYTHON_MODULE}.cpp)

    if (BUILD_COUPLER)
      set(EVERYTRACE_REFADDR ${EVERYTRACE_cf_mpi_REFADDR})
    else()
      set(EVERYTRACE_REFADDR ${EVERYTRACE_cf_REFADDR})
    endif()



    ADD_LIBRARY(${CYTHON_MODULE}_so MODULE ${CYTHON_MODULE}.cpp ${PYLIB_SOURCES} ${EVERYTRACE_REFADDR})

    target_link_libraries(${CYTHON_MODULE}_so icebin ${PYTHON_LIBRARY})

    set_target_properties(${CYTHON_MODULE}_so
        PROPERTIES
        OUTPUT_NAME ${CYTHON_MODULE}
        PREFIX ""
        SUFFIX ".so"
        OUTPUT ${CYTHON_MODULE}.so)



    # Set RPATH in the installed library
    # http://www.cmake.org/Wiki/CMake_RPATH_handling
    # http://www.kitware.com/blog/home/post/510
    set_target_properties(${CYTHON_MODULE}_so PROPERTIES
        INSTALL_RPATH_USE_LINK_PATH TRUE)

    install(TARGETS ${CYTHON_MODULE}_so DESTINATION lib/python${PYTHON_VERSION}/site-packages)

endforeach()

install(DIRECTORY ${PROJECT_SOURCE_DIR}/pylib/ DESTINATION lib/python${PYTHON_VERSION}/site-packages
          FILES_MATCHING PATTERN "*.py")
install(DIRECTORY ${PROJECT_SOURCE_DIR}/pylib/ DESTINATION pylib
          FILES_MATCHING PATTERN "*.pxd")

include_directories(${PROJECT_SOURCE_DIR}/slib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

