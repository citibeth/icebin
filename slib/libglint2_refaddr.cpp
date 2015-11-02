#include <cstdio>
#include <netcdfcpp.h>
#include <everytrace.h>

extern "C"
void libglint2_ncerror_segfault(void)
{
	// Permanently set NetCDF's legacy C++ interface to segfault on error
	NcError *err = new NcError(NcError::verbose_fatal, &everytrace_exit);
}
