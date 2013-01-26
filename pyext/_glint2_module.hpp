#pragma once

#include <Python.h>

/* See http://dsnra.jpl.nasa.gov/software/Python/numpydoc/numpy-13.html

If your extension does not reside in a single file, there is an
additional step that is necessary. Be sure to define the symbol
PY_ARRAY_UNIQUE_SYMBOL to some name (the same name in all the files
comprising the extension), upstream from the include of
arrayobject.h. Typically this would be in some header file that is
included before arrayobject.h. The import_array statement goes into
the init function for the module as before, and not in any of the
other files. Of course, it is ok to define PY_ARRAY_UNIQUE_SYMBOL
symbol even if you only use one file for the extension.

See also: http://www.scipy.org/Cookbook/C_Extensions
http://www.scribd.com/doc/51820831/189/NO-IMPORT-ARRAY

http://projects.scipy.org/numpy/browser/trunk/doc/source/reference/c-api.array.rst?rev=7066

*/

#define PY_ARRAY_UNIQUE_SYMBOL glint2_ARRAY_API

#include <numpy/arrayobject.h>
