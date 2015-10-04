#include <cstdio>

// Produce a reference address that can be used to match symbols to stacktrace
extern "C"
void libglint2_refaddr(void)
{
    fprintf(stderr, "REFERENCE_ADDRESS libglint2_refaddr %p\n", libglint2_refaddr);
    fflush(stderr);
}
