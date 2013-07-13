! Extensions to c_loc() functionality

module c_loc_x

implicit none

CONTAINS

	function c_loc_array_double(A)
        use, intrinsic :: iso_c_binding
        real*8, target :: A(*)
        type(c_ptr) :: c_loc_array_double

        c_loc_array_double = c_loc(A)
	end function c_loc_array_double

	function c_loc_array_int(A)
        use, intrinsic :: iso_c_binding
        integer, target :: A(*)
        type(c_ptr) :: c_loc_array_int

        c_loc_array_int = c_loc(A)
	end function c_loc_array_int

end module c_loc_x
