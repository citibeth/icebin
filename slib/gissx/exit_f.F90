
#ifdef __GNUC__
subroutine print_stacktrace() bind(c)
	use, intrinsic :: iso_c_binding
	call backtrace
end subroutine print_stacktrace
#endif
