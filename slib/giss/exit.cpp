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

#ifdef __GNUC__
#include <execinfo.h>
// http://stackoverflow.com/questions/9053658/correct-format-specifier-to-print-pointer-address
#include <inttypes.h>		// C-99
#include <cstdlib>
#endif

namespace giss {

void exit_exception(int errcode)
{
	throw std::exception();
}

void exit_segfault(int errcode)
{
	int *ptr = 0;
	*ptr = 17;
}

#ifdef __GNUC__
// http://stackoverflow.com/questions/77005/how-to-generate-a-stacktrace-when-my-gcc-c-app-crashes

void exit_stacktrace(int errcode)
{
	const int MAX_TRACE = 200;
	void *trace[MAX_TRACE];

	fprintf(stderr, "User stacktrace:\n");

	size_t ntrace = backtrace(trace, MAX_TRACE);

	for (size_t i=0; i<ntrace; ++i) {
		fprintf(stderr, "#%d 0x%lx\n", i, (uintptr_t)(trace[i]));
	}

	throw std::exception();
}
	std::function<void(int)> exit(&exit_stacktrace);
#else
	std::function<void(int)> exit(&exit_exception);
#endif


}
