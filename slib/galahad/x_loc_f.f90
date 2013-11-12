! GLINT2: A Coupling Library for Ice Models and GCMs
! Copyright (c) 2013 by Robert Fischer
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

!> Work-around of a bug in GNU gfortran, in which c_loc() does not work
!! when applied directly to elements of arrays.  The bug is fixed in GCC 4.9

!> Container modulefor c_loc() workaround.
module c_loc_x

implicit none

CONTAINS

	!> Used to get the C address of a real*8 array
	!! @param A The array in question
	function c_loc_array_double(A)
        use, intrinsic :: iso_c_binding
        real*8, target :: A(*)
        type(c_ptr) :: c_loc_array_double

        c_loc_array_double = c_loc(A)
	end function c_loc_array_double

	!> Used to get the C address of a integer array
	!! @param A The array in question
	function c_loc_array_int(A)
        use, intrinsic :: iso_c_binding
        integer, target :: A(*)
        type(c_ptr) :: c_loc_array_int

        c_loc_array_int = c_loc(A)
	end function c_loc_array_int

end module c_loc_x
