

include_directories(${PROJECT_SOURCE_DIR}/slib)
include_directories(${PYTHON_INCLUDES})


ADD_LIBRARY(_glint2 MODULE
	pyutil.cpp
	NcFile_py.cpp
	Grid_py.cpp
	# IceSheet_py.cpp
	MatrixMaker_py.cpp
	matrix_ops_py.cpp
	_glint2_module_py.cpp
)

include_directories(${PYTHON_INCLUES} ${NUMPY_INCLUDE_DIRS})

#set(CFLAGS ${CFLAGS} ${PYTHON_INCLUES} ${NUMPY_INCLUDE_DIRS})

target_link_libraries(_glint2 glint2 ${PYTHON_LIBRARY})

set_target_properties(_glint2
	PROPERTIES
	PREFIX ""
	SUFFIX ".so"
	OUTPUT _glint2.so)


# Set RPATH in the installed library
# http://www.cmake.org/Wiki/CMake_RPATH_handling
# http://www.kitware.com/blog/home/post/510
set_target_properties(_glint2 PROPERTIES
	INSTALL_RPATH_USE_LINK_PATH TRUE)



# NOTE: Some pythons want .../python/site-packages instead of .../python2.7/site-packages
# In this case, create a symbolic link!

message(STATUS "PY_SYSCONFIG_PLATLIB: ${PY_SYSCONFIG_PLATLIB}")
install(TARGETS _glint2 DESTINATION ${PY_SYSCONFIG_PLATLIB})

install(FILES
	${CMAKE_CURRENT_SOURCE_DIR}/glint2/__init__.py
	${CMAKE_CURRENT_SOURCE_DIR}/glint2/cext.py
	${CMAKE_CURRENT_SOURCE_DIR}/glint2/Grid.py
	${CMAKE_CURRENT_SOURCE_DIR}/glint2/io.py
		DESTINATION ${PY_SYSCONFIG_PURELIB}/glint2)



