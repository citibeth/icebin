MODULE qpt_x
	use, intrinsic :: iso_c_binding
	use zd11_x

IMPLICIT NONE

	type QPT_problem_c
		type(c_ptr) :: main			! QTP_problem_type

	    type(c_ptr) :: m       ! int &; number of constraints
	    type(c_ptr) :: n       ! int &; number of variables
	    type(c_ptr) :: f   ! double & constant term

	    type(c_ptr) :: G
	    type(c_ptr) :: X_l
	    type(c_ptr) :: X_u
	    type(c_ptr) :: C_l
	    type(c_ptr) :: C_u
	    type(c_ptr) :: C
	    type(c_ptr) :: X
	    type(c_ptr) :: Y
	    type(c_ptr) :: Z

	    type(ZD11_c) :: A, H
	end type QPT_problem_c

END MODULE qpt_x

! ======================================================
! Functions to be called from C, so they're outside of a module


! @param m Number of constraints
! @param n Number of variables
subroutine QPT_problem_c_init(self, main, m, n, A_ne, H_ne, eqp_bool) bind(c)
use GALAHAD_QPT_double
USE qpt_x
use c_loc_x
use, intrinsic :: iso_c_binding
IMPLICIT NONE
type(QPT_problem_c) :: self
type(QPT_problem_type), target :: main
integer(c_int), value :: m, n
integer(c_int), value :: A_ne		! # of elements in A (constraint) matrix
integer(c_int), value :: H_ne		! # of elements in H (Hessian) matrix
integer(c_int), value :: eqp_bool		! Are we preparing for a problem w/ equality constraints? (1=true, 0=false)

	main%m = m
	main%n = n

	self%main = c_loc(main)

	self%m = c_loc(main%m)
	self%n = c_loc(main%n)
	self%f = c_loc(main%f)

	allocate(main%G(n), main%X_l(n), main%X_u(n))
	self%G = c_loc_array_double(main%G)
	self%X_l = c_loc_array_double(main%X_l)
	self%X_u = c_loc_array_double(main%X_u)

! write(6,*) 'QPT_problem-c_init() eqp_bool=', eqp_bool
	if (eqp_bool /= 0) then
		allocate(main%C(m))
		self%C = c_loc_array_double(main%C)
! write(6,*) 'self%C=', self%C
	else
		allocate(main%C_l(m), main%C_u(m))
		self%C_l = c_loc_array_double(main%C_l)
		self%C_u = c_loc_array_double(main%C_u)
! write(6,*) 'self%C_u=', self%C_u
	end if

	allocate(main%X(n), main%Y(m), main%Z(n))
	self%X = c_loc_array_double(main%X)
	self%Y = c_loc_array_double(main%Y)
	self%Z = c_loc_array_double(main%Z)

write(6,*) 'ZD11_c_init(self%A)', A_ne
	call ZD11_c_init(self%A, main%A, m, n, A_ne)
write(6,*) 'ZD11_c_init(self%H)'
	call ZD11_c_init(self%H, main%H, n, n, H_ne)
end subroutine QPT_problem_c_init

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

subroutine QPT_problem_delete_c(main_c) bind(c)
use GALAHAD_QPT_double
USE qpt_x
use, intrinsic :: iso_c_binding
IMPLICIT NONE
type(c_ptr), value :: main_c
	type(QPT_problem_type), pointer :: main

	call c_f_pointer(main_c, main)
	deallocate(main)
end subroutine QPT_problem_delete_c
