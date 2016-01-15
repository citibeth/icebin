/*
 * GLINT2: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013 by Robert Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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
