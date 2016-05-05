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

#include <cstdio>
#include <glint2/util.hpp>
//#include <everytrace.h>
#include <cstdio>

#ifdef __GNUC__
#include <execinfo.h>
// http://stackoverflow.com/questions/9053658/correct-format-specifier-to-print-pointer-address
#include <inttypes.h>		// C-99
#include <cstdlib>
#endif

namespace giss {

static void my_exit(int err)
{
    fprintf(stderr, "Exiting now...\n");
    int *p = 0;
    *p=17;
}

std::function<void(int)> exit(&my_exit);

}
