! ======== DO NOT EDIT!!!!  Machine Generated!!!!
module f90blitz

use, intrinsic :: iso_c_binding
implicit none
	
	type, bind(c) :: arr_spec_1
		type(c_ptr) :: base
		type(c_ptr) :: deltas(1)
		integer(c_int) :: lbounds(1)
		integer(c_int) :: ubounds(1)
	end type
	type, bind(c) :: arr_spec_2
		type(c_ptr) :: base
		type(c_ptr) :: deltas(2)
		integer(c_int) :: lbounds(2)
		integer(c_int) :: ubounds(2)
	end type
	type, bind(c) :: arr_spec_3
		type(c_ptr) :: base
		type(c_ptr) :: deltas(3)
		integer(c_int) :: lbounds(3)
		integer(c_int) :: ubounds(3)
	end type
contains

	function c_loc_double(x)
		use, intrinsic :: iso_c_binding
		real*8, target :: x
		type(c_ptr) :: c_loc_double

		c_loc_double = c_loc(x)
	end function
	! ================ Type (real*8, double), Rank 1
	subroutine get_spec_double_1(arr, low1, spec)
	implicit none
	real*8, dimension(:), target :: arr
	integer :: low1
	type(arr_spec_1) :: spec

		spec%base = c_loc_double( arr(lbound(arr,1)) )
		

		! ------- Dimension 1
		spec%lbounds(1) = low1
		spec%ubounds(1) = low1 + ubound(arr,1) - 1
		if (spec%lbounds(1) < spec%ubounds(1)) then
			spec%deltas(1) = c_loc_double(arr(lbound(arr,1)+1))
		else
			spec%deltas(1) = spec%base
		end if
	end subroutine

	! ================ Type (real*8, double), Rank 2
	subroutine get_spec_double_2(arr, low1,low2, spec)
	implicit none
	real*8, dimension(:,:), target :: arr
	integer :: low1,low2
	type(arr_spec_2) :: spec

		spec%base = c_loc_double( arr(lbound(arr,1),lbound(arr,2)) )
		

		! ------- Dimension 1
		spec%lbounds(1) = low1
		spec%ubounds(1) = low1 + ubound(arr,1) - 1
		if (spec%lbounds(1) < spec%ubounds(1)) then
			spec%deltas(1) = c_loc_double(arr(lbound(arr,1)+1,lbound(arr,2)))
		else
			spec%deltas(1) = spec%base
		end if

		! ------- Dimension 2
		spec%lbounds(2) = low2
		spec%ubounds(2) = low2 + ubound(arr,2) - 1
		if (spec%lbounds(2) < spec%ubounds(2)) then
			spec%deltas(2) = c_loc_double(arr(lbound(arr,1),lbound(arr,2)+1))
		else
			spec%deltas(2) = spec%base
		end if
	end subroutine

	! ================ Type (real*8, double), Rank 3
	subroutine get_spec_double_3(arr, low1,low2,low3, spec)
	implicit none
	real*8, dimension(:,:,:), target :: arr
	integer :: low1,low2,low3
	type(arr_spec_3) :: spec

		spec%base = c_loc_double( arr(lbound(arr,1),lbound(arr,2),lbound(arr,3)) )
		

		! ------- Dimension 1
		spec%lbounds(1) = low1
		spec%ubounds(1) = low1 + ubound(arr,1) - 1
		if (spec%lbounds(1) < spec%ubounds(1)) then
			spec%deltas(1) = c_loc_double(arr(lbound(arr,1)+1,lbound(arr,2),lbound(arr,3)))
		else
			spec%deltas(1) = spec%base
		end if

		! ------- Dimension 2
		spec%lbounds(2) = low2
		spec%ubounds(2) = low2 + ubound(arr,2) - 1
		if (spec%lbounds(2) < spec%ubounds(2)) then
			spec%deltas(2) = c_loc_double(arr(lbound(arr,1),lbound(arr,2)+1,lbound(arr,3)))
		else
			spec%deltas(2) = spec%base
		end if

		! ------- Dimension 3
		spec%lbounds(3) = low3
		spec%ubounds(3) = low3 + ubound(arr,3) - 1
		if (spec%lbounds(3) < spec%ubounds(3)) then
			spec%deltas(3) = c_loc_double(arr(lbound(arr,1),lbound(arr,2),lbound(arr,3)+1))
		else
			spec%deltas(3) = spec%base
		end if
	end subroutine


	function c_loc_int(x)
		use, intrinsic :: iso_c_binding
		integer, target :: x
		type(c_ptr) :: c_loc_int

		c_loc_int = c_loc(x)
	end function
	! ================ Type (integer, int), Rank 1
	subroutine get_spec_int_1(arr, low1, spec)
	implicit none
	integer, dimension(:), target :: arr
	integer :: low1
	type(arr_spec_1) :: spec

		spec%base = c_loc_int( arr(lbound(arr,1)) )
		

		! ------- Dimension 1
		spec%lbounds(1) = low1
		spec%ubounds(1) = low1 + ubound(arr,1) - 1
		if (spec%lbounds(1) < spec%ubounds(1)) then
			spec%deltas(1) = c_loc_int(arr(lbound(arr,1)+1))
		else
			spec%deltas(1) = spec%base
		end if
	end subroutine

	! ================ Type (integer, int), Rank 2
	subroutine get_spec_int_2(arr, low1,low2, spec)
	implicit none
	integer, dimension(:,:), target :: arr
	integer :: low1,low2
	type(arr_spec_2) :: spec

		spec%base = c_loc_int( arr(lbound(arr,1),lbound(arr,2)) )
		

		! ------- Dimension 1
		spec%lbounds(1) = low1
		spec%ubounds(1) = low1 + ubound(arr,1) - 1
		if (spec%lbounds(1) < spec%ubounds(1)) then
			spec%deltas(1) = c_loc_int(arr(lbound(arr,1)+1,lbound(arr,2)))
		else
			spec%deltas(1) = spec%base
		end if

		! ------- Dimension 2
		spec%lbounds(2) = low2
		spec%ubounds(2) = low2 + ubound(arr,2) - 1
		if (spec%lbounds(2) < spec%ubounds(2)) then
			spec%deltas(2) = c_loc_int(arr(lbound(arr,1),lbound(arr,2)+1))
		else
			spec%deltas(2) = spec%base
		end if
	end subroutine

	! ================ Type (integer, int), Rank 3
	subroutine get_spec_int_3(arr, low1,low2,low3, spec)
	implicit none
	integer, dimension(:,:,:), target :: arr
	integer :: low1,low2,low3
	type(arr_spec_3) :: spec

		spec%base = c_loc_int( arr(lbound(arr,1),lbound(arr,2),lbound(arr,3)) )
		

		! ------- Dimension 1
		spec%lbounds(1) = low1
		spec%ubounds(1) = low1 + ubound(arr,1) - 1
		if (spec%lbounds(1) < spec%ubounds(1)) then
			spec%deltas(1) = c_loc_int(arr(lbound(arr,1)+1,lbound(arr,2),lbound(arr,3)))
		else
			spec%deltas(1) = spec%base
		end if

		! ------- Dimension 2
		spec%lbounds(2) = low2
		spec%ubounds(2) = low2 + ubound(arr,2) - 1
		if (spec%lbounds(2) < spec%ubounds(2)) then
			spec%deltas(2) = c_loc_int(arr(lbound(arr,1),lbound(arr,2)+1,lbound(arr,3)))
		else
			spec%deltas(2) = spec%base
		end if

		! ------- Dimension 3
		spec%lbounds(3) = low3
		spec%ubounds(3) = low3 + ubound(arr,3) - 1
		if (spec%lbounds(3) < spec%ubounds(3)) then
			spec%deltas(3) = c_loc_int(arr(lbound(arr,1),lbound(arr,2),lbound(arr,3)+1))
		else
			spec%deltas(3) = spec%base
		end if
	end subroutine


end module f90blitz
