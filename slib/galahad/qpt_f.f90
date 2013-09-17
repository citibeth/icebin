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

MODULE qpt_x
	use, intrinsic :: iso_c_binding
	use hsl_zd11_double
	use zd11_x

IMPLICIT NONE

	type QPT_problem_c
		type(c_ptr) :: this_fc		! QTP_problem_type * (peer)

	    type(c_ptr) :: m       ! int &; number of constraints
	    type(c_ptr) :: n       ! int &; number of variables
	    type(c_ptr) :: f   ! double & constant term
		type(c_ptr) :: Hessian_kind

	    type(c_ptr) :: G
	    type(c_ptr) :: X_l
	    type(c_ptr) :: X_u
	    type(c_ptr) :: C_l
	    type(c_ptr) :: C_u
	    type(c_ptr) :: C
	    type(c_ptr) :: X
	    type(c_ptr) :: Y
	    type(c_ptr) :: Z
!		type(c_ptr) :: WEIGHT
		type(c_ptr) :: X0

	    type(ZD11_c) :: A, H
	end type QPT_problem_c

END MODULE qpt_x

! ======================================================
! Functions to be called from C, so they're outside of a module

! ----------------------------------------------------------
! @param m Number of constraints
! @param n Number of variables
subroutine QPT_problem_init_c(this_cc, this_fc, m, n, eqp, Hessian_kind) bind(c)
use GALAHAD_QPT_double
USE qpt_x
use zd11_x
use c_loc_x
use, intrinsic :: iso_c_binding
IMPLICIT NONE
type(c_ptr), value :: this_cc			! qpt_problem_c *
	type(qpt_problem_c), pointer :: this_c
type(c_ptr), value :: this_fc			! qpt_problem_type *
	type(QPT_problem_type), pointer :: this_f
integer(c_int), value :: m, n
logical(kind=c_bool), value :: eqp		! Are we preparing for a problem w/ equality constraints?
integer, value :: Hessian_kind

	call c_f_pointer(this_fc, this_f)
	call c_f_pointer(this_cc, this_c)


		this_f%m = m
		this_f%n = n

		this_c%this_fc = c_loc(this_f)

		this_c%m = c_loc(this_f%m)
		this_c%n = c_loc(this_f%n)
		this_c%f = c_loc(this_f%f)
		this_c%Hessian_kind = c_loc(this_f%Hessian_kind)
		this_f%Hessian_kind = Hessian_kind

		if (Hessian_kind > 0) then
			allocate(this_f%X0(n))
			this_c%X0 = c_loc_array_double(this_f%X0)
		end if
!		! NOTE: Other Hessian_kinds are not implemented

		allocate(this_f%G(n), this_f%X_l(n), this_f%X_u(n))
		this_c%G = c_loc_array_double(this_f%G)
		this_c%X_l = c_loc_array_double(this_f%X_l)
		this_c%X_u = c_loc_array_double(this_f%X_u)

	! write(6,*) 'QPT_problem-c_init() eqp=', eqp
		if (eqp) then
			allocate(this_f%C(m))
			this_c%C = c_loc_array_double(this_f%C)
	! write(6,*) 'this_c%C=', this_c%C
		else
			allocate(this_f%C_l(m), this_f%C_u(m))
			this_c%C_l = c_loc_array_double(this_f%C_l)
			this_c%C_u = c_loc_array_double(this_f%C_u)
	! write(6,*) 'this_c%C_u=', this_c%C_u
		end if

		allocate(this_f%X(n), this_f%Y(m), this_f%Z(n))
		this_c%X = c_loc_array_double(this_f%X)
		this_c%Y = c_loc_array_double(this_f%Y)
		this_c%Z = c_loc_array_double(this_f%Z)

!print *, c_loc(this_f%A)
!write(6,*) 'ZD11_init(this_c%A)', A_ne
!		ptr => this_f%A
!		call ZD11_init_f(this_c%A, ptr, m, n, A_ne)
!write(6,*) 'ZD11_init(this_c%H)'
!		call ZD11_init_f(this_c%H, this_f%H, n, n, H_ne)

end subroutine QPT_problem_init_c
! ----------------------------------------------------------
subroutine QPT_problem_alloc_H(this_cc, H_ne) bind(c)
use GALAHAD_QPT_double
USE qpt_x
use zd11_x
use c_loc_x
use, intrinsic :: iso_c_binding
IMPLICIT NONE
type(c_ptr), value :: this_cc			! qpt_problem_c *
	type(qpt_problem_c), pointer :: this_c
type(QPT_problem_type), pointer :: this_f
integer(c_int), value :: H_ne		! # of elements in H (Hessian) matrix

	call c_f_pointer(this_cc, this_c)
	call c_f_pointer(this_c%this_fc, this_f)
	write(6,*) 'ZD11_init(this_c%H)',H_ne
	call ZD11_init_f(this_c%H, this_f%H, this_f%n, this_f%n, H_ne)

end subroutine QPT_problem_alloc_H
! ----------------------------------------------------------
subroutine QPT_problem_alloc_A(this_cc, A_ne) bind(c)
use GALAHAD_QPT_double
USE qpt_x
use zd11_x
use c_loc_x
use, intrinsic :: iso_c_binding
IMPLICIT NONE
type(c_ptr), value :: this_cc			! qpt_problem_c *
	type(qpt_problem_c), pointer :: this_c
type(QPT_problem_type), pointer :: this_f
integer(c_int), value :: A_ne		! # of elements in A (constraints) matrix

	call c_f_pointer(this_cc, this_c)
	call c_f_pointer(this_c%this_fc, this_f)
	write(6,*) 'ZD11_init(this_c%A)'
	call ZD11_init_f(this_c%A, this_f%A, this_f%m, this_f%n, A_ne)

end subroutine QPT_problem_alloc_A
! ----------------------------------------------------------
function QPT_problem_new_c() bind(c)
use GALAHAD_QPT_double
USE qpt_x
use, intrinsic :: iso_c_binding
IMPLICIT NONE
type(c_ptr) :: QPT_problem_new_c

	type(c_ptr) :: self

	type(QPT_problem_type), pointer :: main
	allocate(main)
	QPT_problem_new_c = c_loc(main)
end function QPT_problem_new_c
! ---------------------------------------------------------------------
subroutine QPT_problem_delete_c(this_fc) bind(c)
use GALAHAD_QPT_double
USE qpt_x
use, intrinsic :: iso_c_binding
IMPLICIT NONE
type(c_ptr), value :: this_fc
	type(QPT_problem_type), pointer :: this_f

	call c_f_pointer(this_fc, this_f)
	deallocate(this_f)
end subroutine QPT_problem_delete_c
