
pism_set_dependencies()
set(Glint2_EXTERNAL_LIBS ${Pism_EXTERNAL_LIBS})
glint2_set_dependencies()



# Make sure the linker can find the Glint2 library once it is built.
link_directories (${PROJECT_BINARY_DIR}/slib)


###############################################
# Find Python libraries.
#
set(PYTHON_LIBRARY        "C:/Python27/libs/python27.lib")
set(PYTHON_INCLUDE_DIR    "C:/Python27/include")
FIND_PACKAGE(PythonLibs)
IF(PYTHONLIBS_FOUND)
  INCLUDE_DIRECTORIES("${PYTHON_INCLUDE_DIRS}")
ELSE()
  MESSAGE(FATAL_ERROR "Unable to find PythonLibs.")
ENDIF()

###############################################

ADD_LIBRARY(_glint2 MODULE
	Grid_py.cpp
	IceSheet_py.cpp
	MatrixMaker_py.cpp
	NcFile_py.cpp
	_glint2_module_py.cpp
	matrix_ops_py.cpp
	pyutil.cpp
)

set_target_properties(_glint2
	PROPERTIES
	SUFFIX ".pyd")

TARGET_LINK_LIBRARIES(_glint2 glint2 ${PYTHON_LIBRARIES})

