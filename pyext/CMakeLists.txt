
find_package(python)
find_package(cython)
find_package(numpy)

add_custom_command(OUTPUT _icebin.cpp
	MAIN_DEPENDENCY ${PROJECT_SOURCE_DIR}/pyext/_icebin.pyx
	DEPENDS cpp.pxd
	COMMAND ${CYTHON_EXECUTABLE} --cplus
		--include-dir ${PROJECT_SOURCE_DIR}/pyext
		${PROJECT_SOURCE_DIR}/pyext/_icebin.pyx
		-o _icebin.cpp)


include_directories(${PROJECT_SOURCE_DIR}/slib)
include_directories(${PYTHON_INCLUDES} ${NUMPY_INCLUDE_DIRS})

ADD_LIBRARY(_icebin MODULE _icebin.cpp
	${EVERYTRACE_REFADDR}
)


#set(CFLAGS ${CFLAGS} ${PYTHON_INCLUES} ${NUMPY_INCLUDE_DIRS})

target_link_libraries(_icebin icebin ${PYTHON_LIBRARY})

set_target_properties(_icebin
	PROPERTIES
	PREFIX ""
	SUFFIX ".so"
	OUTPUT _icebin.so)


# Set RPATH in the installed library
# http://www.cmake.org/Wiki/CMake_RPATH_handling
# http://www.kitware.com/blog/home/post/510
set_target_properties(_icebin PROPERTIES
	INSTALL_RPATH_USE_LINK_PATH TRUE)



# NOTE: Some pythons want .../python/site-packages instead of .../python2.7/site-packages
# In this case, create a symbolic link!

message(STATUS "PY_SYSCONFIG_PLATLIB: ${PY_SYSCONFIG_PLATLIB}")
install(TARGETS _icebin DESTINATION ${PY_SYSCONFIG_PLATLIB})


#install(FILES
#	${CMAKE_CURRENT_SOURCE_DIR}/icebin/__init__.py
#	${CMAKE_CURRENT_SOURCE_DIR}/icebin/cext.py
#	${CMAKE_CURRENT_SOURCE_DIR}/icebin/Grid.py
#	${CMAKE_CURRENT_SOURCE_DIR}/icebin/io.py
#		DESTINATION ${PY_SYSCONFIG_PURELIB}/icebin)



