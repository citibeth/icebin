
find_package(python)
find_package(cython)
find_package(numpy)

foreach (CYTHON_MODULE icebin ibmisc)
	add_custom_command(OUTPUT ${CYTHON_MODULE}.cpp
		MAIN_DEPENDENCY ${PROJECT_SOURCE_DIR}/pyext/${CYTHON_MODULE}.pyx
		DEPENDS cpp.pxd
		COMMAND ${CYTHON_EXECUTABLE} --cplus -3
			--include-dir ${PROJECT_SOURCE_DIR}/pyext
			${PROJECT_SOURCE_DIR}/pyext/${CYTHON_MODULE}.pyx
			-o ${CYTHON_MODULE}.cpp)

	ADD_LIBRARY(${CYTHON_MODULE}_so MODULE ${CYTHON_MODULE}.cpp ${EVERYTRACE_REFADDR})

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

	install(TARGETS ${CYTHON_MODULE}_so DESTINATION ${PY_SYSCONFIG_PLATLIB})

endforeach()




include_directories(${PROJECT_SOURCE_DIR}/slib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${PYTHON_INCLUDES} ${NUMPY_INCLUDE_DIRS})



#set(CFLAGS ${CFLAGS} ${PYTHON_INCLUES} ${NUMPY_INCLUDE_DIRS})






# NOTE: Some pythons want .../python/site-packages instead of .../python2.7/site-packages
# In this case, create a symbolic link!

message(STATUS "PY_SYSCONFIG_PLATLIB: ${PY_SYSCONFIG_PLATLIB}")


#install(FILES
#	${CMAKE_CURRENT_SOURCE_DIR}/icebin/__init__.py
#	${CMAKE_CURRENT_SOURCE_DIR}/icebin/cext.py
#	${CMAKE_CURRENT_SOURCE_DIR}/icebin/Grid.py
#	${CMAKE_CURRENT_SOURCE_DIR}/icebin/io.py
#		DESTINATION ${PY_SYSCONFIG_PURELIB}/icebin)



