! ======== DO NOT EDIT!!!!  Macine Generated!!!!
module f90blitz

use, intrinsic :: iso_c_binding
implicit none
	
	type, bind(c) :: arr_spec_3
		type(c_ptr) :: base
		type(c_ptr) :: deltas(3)
		integer(c_int) :: lbounds(3)
		integer(c_int) :: ubounds(3)
	end type
contains

	function c_loc_float(x)
		use, intrinsic :: iso_c_binding
		real*4, target :: x
		type(c_ptr) :: c_loc_float

		c_loc_float = c_loc(x)
	end function
	! ================ Type (real*4, float), Rank 3
	subroutine get_spec_float_3(arr, low1,low2,low3, spec)
	implicit none
	real*4, dimension(:,:,:), target :: arr
	integer :: low1,low2,low3
	type(arr_spec_3) :: spec

		spec%base = c_loc_float( arr(lbound(arr,1),lbound(arr,2),lbound(arr,3)) )
		

		! ------- Dimension 1
		spec%lbounds(1) = low1
		spec%ubounds(1) = low1 + ubound(arr,1) - 1
		if (spec%lbounds(1) < spec%ubounds(1)) then
			spec%deltas(1) = c_loc_float(arr(lbound(arr,1)+1,lbound(arr,2),lbound(arr,3)))
		else
			spec%deltas(1) = spec%base
		end if

		! ------- Dimension 2
		spec%lbounds(2) = low2
		spec%ubounds(2) = low2 + ubound(arr,2) - 1
		if (spec%lbounds(2) < spec%ubounds(2)) then
			spec%deltas(2) = c_loc_float(arr(lbound(arr,1),lbound(arr,2)+1,lbound(arr,3)))
		else
			spec%deltas(2) = spec%base
		end if

		! ------- Dimension 3
		spec%lbounds(3) = low3
		spec%ubounds(3) = low3 + ubound(arr,3) - 1
		if (spec%lbounds(3) < spec%ubounds(3)) then
			spec%deltas(3) = c_loc_float(arr(lbound(arr,1),lbound(arr,2),lbound(arr,3)+1))
		else
			spec%deltas(3) = spec%base
		end if
	end subroutine


end module f90blitz
