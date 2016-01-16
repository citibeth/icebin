# WARNING: This setup.py is not working
# It is being kept for posterity.  We use the native CMake-based build

from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

include_path = [numpy.get_include()]

# http://docs.cython.org/src/userguide/wrapping_CPlusPlus.html
setup(ext_modules = cythonize(Extension(
	# The extension name
	"icebin",

	# The Cython source and additional C++ source files
	sources=["icebin.pyx"],

	language="c++",                        # generate and compile C++ code

	# We're wrapping the IceBin C++ library
	# http://docs.cython.org/src/tutorial/clibraries.html
#	libraries="icebin",
	libraries="/Users/rpfische/git/icebin/build/slib/libicebin.dylib",

	library_dirs = [],				# -L options
	runtime_library_dirs = [],		# -rpath options

	include_dirs = [],
)))


# Compiling from command line:
# CFLAGS="-I/usr/local/otherdir/calg/include"  \
# LDFLAGS="-L/usr/local/otherdir/calg/lib"     \
#    python setup.py build_ext -i


# For CMake see:
# https://github.com/thewtex/cython-cmake-example
