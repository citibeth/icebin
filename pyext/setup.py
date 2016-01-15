from distutils.core import setup
from Cython.Build import cythonize

# http://docs.cython.org/src/userguide/wrapping_CPlusPlus.html
setup(ext_modules = cythonize(Extension(
	# The extension name
	"icebin",

	# The Cython source and additional C++ source files
	sources=["icebin.pyx"],

	language="c++",                        # generate and compile C++ code

	# We're wrapping the IceBin C++ library
	# http://docs.cython.org/src/tutorial/clibraries.html
	libraries="icebin"
)))


# Compiling from command line:
# CFLAGS="-I/usr/local/otherdir/calg/include"  \
# LDFLAGS="-L/usr/local/otherdir/calg/lib"     \
#    python setup.py build_ext -i


# For CMake see:
# https://github.com/thewtex/cython-cmake-example