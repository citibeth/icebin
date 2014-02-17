
pism_set_dependencies()
set(Glint2_EXTERNAL_LIBS ${Pism_EXTERNAL_LIBS})
glint2_set_dependencies()



# Make sure the linker can find the Glint2 library once it is built.
link_directories (${PROJECT_BINARY_DIR}/slib)

# Add executable called "helloDemo" that is built from the source files
# "demo.cxx" and "demo_b.cxx". The extensions are automatically found.
add_executable (overlap overlap.cpp)

# Link the executable to the Glint2 library.
target_link_libraries (overlap glint2 ${Glint2_EXTERNAL_LIBS}) 

add_executable (desm desm.cpp)
target_link_libraries (desm glint2 ${Glint2_EXTERNAL_LIBS}) 



# ================================================

install(TARGETS overlap desm
	DESTINATION bin)

# Set RPATH in the installed executable
# http://www.cmake.org/pipermail/cmake/2010-February/035157.html
set_target_properties(overlap 		# more targets here...
	PROPERTIES
	INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib
	INSTALL_RPATH_USE_LINK_PATH TRUE)
