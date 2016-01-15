#include <cstdio>
#include <cstdarg>
#include <exception>
#include <icebin/error.hpp>
#ifdef USE_EVERYTRACE
#include <everytrace.h>
#endif

namespace icebin {

void default_error(int retcode, const char *format, ...)
{
	va_list arglist;

	va_start(arglist, format);
	vfprintf(stderr, format, arglist);
	va_end(arglist);
	fprintf(stderr, "\n");

#ifdef USE_EVERYTRACE
	everytrace_exit(-1);
#endif
	throw icebin::Exception();
//	exit(-1);
}

error_ptr icebin_error = &default_error;

} 	// Namespace
