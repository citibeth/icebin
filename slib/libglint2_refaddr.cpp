#include <cstdio>
#include <netcdfcpp.h>

// Produce a reference address that can be used to match symbols to stacktrace
extern "C"
void libglint2_refaddr(void)
{
    fprintf(stderr, "REFERENCE_ADDRESS libglint2_refaddr %p\n", libglint2_refaddr);
    fflush(stderr);
}


static void glint2_exit(int code)
{
    printf("Exiting with code %d\n", code);
    int *ptr = 0;
	*ptr = 17;			// Cause a segfault
}

extern "C"
void libglint2_ncerror_segfault(void)
{
	// Permanently set NetCDF's legacy C++ interface to segfault on error
	NcError *err = new NcError(NcError::verbose_fatal, &glint2_exit);
}
