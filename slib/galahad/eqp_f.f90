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

! -------------------------------------------------------------------
! Calls the EQP subroutine
! This function is meant to be called from C/C++
function eqp_solve_simple(p_c, infinity) bind(C)
use iso_c_binding
use qpt_x
! USE GALAHAD_QP_double
! USE GALAHAD_QPT_double		! Debugging
USE GALAHAD_EQP_double
IMPLICIT NONE
type(c_ptr), value :: p_c		! QPT_problem_type
real(c_double), value :: infinity
logical(kind=c_bool) :: eqp_solve_simple

	type(QPT_problem_type), pointer :: p
	integer :: errcode

	integer :: time0_ms, time1_ms
	real*8 :: delta_time

	! --- Galahad Stuff
	INTEGER, PARAMETER :: wp = KIND( 1.0D+0 ) ! set precision
!	REAL ( KIND = wp ), PARAMETER :: infinity = 10.0_wp ** 20
!	INTEGER, ALLOCATABLE, DIMENSION( : ) :: C_stat, B_stat
	TYPE ( EQP_data_type ) :: data	! Used as tmp storage, we don't touch it.
	TYPE ( EQP_control_type ) :: control
	TYPE ( EQP_inform_type ) :: inform

	! --------------------------------
	call c_f_pointer(p_c, p)

!	ALLOCATE( B_stat( p%n ), C_stat( p%m ) )
	p%new_problem_structure = .TRUE.

	! ------------ problem data complete, set up control
	CALL EQP_initialize( data, control, inform ) ! Initialize control parameters
	! control%infinity = infinity	! Only for qp_control_type

	! Set infinity
!	control%quadratic_programming_solver = 'qpa' ! use QPA.  (This is important in getting it to work at all).
!!	control%scale = 7		! Sinkhorn-Knopp scaling: Breaks things!
!!	control%scale = 1		! Fast and accurate for regridding
!!	control%scale = 0		! No scaling: slow with proper derivative weights
!	control%scale = 1

!	control%generate_sif_file = .TRUE.
	control%print_level = 5

print *,'Hessian_kind',p%Hessian_kind
print *,'m,n,A%m,A%n,A%ne',p%m,p%n,p%A%m,p%A%n,p%A%ne
!print *,'p%A%row',p%A%row
!print *,'p%A%col',p%A%col
!print *,'p%A%val',p%A%val
!call QPT_A_from_C_to_S(p, errcode)
!print *,'errcode',errcode

print *
print *,'H%m,H%n,H%ne',p%H%m,p%H%n,p%H%ne
!print *,'p%H%row',p%H%row
!print *,'p%H%col',p%H%col
!print *,'p%H%val',p%H%val


!print *,'p%X_l',p%X_l(1:20)
!print *,'p%X_u',p%X_u(1:20)

	! Causes error on samples with lat/lon and cartesian grid
	!!control%presolve = .TRUE.
	!control%presolve = .FALSE.


	! From Nick Gould:
	! When I run your examples, I find that it is best to use the
	!   control%SBLS_control%preconditioner = 1
	! This builds a "constraint preconditioner" in which the Hessian
	! of your quadratic is replaced by the identity matrix, and then
	! runs the projected CG algorithm. The number of nonzeros in the
	! factors appears to be O(n), e.g. for the 25km problem there are
	! 1412, 3528 entries in the preconditioning matrix and its factors
	! while if instead I use the full matrix
	!   control%SBLS_control%preconditioner = 2
	! there are 30563, 4565576 entries in the matrix and its factors
	! which looks like O(n^2) - this is slightly misleading as when
	! the Hessian or its approximation is diagonal, a more powerful
	! Schur-complement factorization is appropriate (and indeed used
	! by EQP). The preconditioner = 1 requires 13 conjugate-gradient
	! iterations that suggests the preconditioner is effective; of course
	! with preconditioner = 2 only 1 iteration is required, but the
	! factorization cost is far higher. The better preconditioned method
	! seems to take 0.06 seconds on my desktop machine for the 25km problem.

!    0 automatic
!    1 explicit with G = I
!    2 explicit with G = H
!    3 explicit with G = diag(max(H,min_diag))
!    4 explicit with G = band(H)
!    5 explicit with G = (optional, diagonal) D
!   11 explicit with G_11 = 0, G_21 = 0, G_22 = H_22
!   12 explicit with G_11 = 0, G_21 = H_21, G_22 = H_22
!   -1 implicit with G_11 = 0, G_21 = 0, G_22 = I
!   -2 implicit with G_11 = 0, G_21 = 0, G_22 = H_22
	control%SBLS_control%preconditioner = 1	! Random low-elevation points break this
!	control%SBLS_control%preconditioner = 3	! Works even with low-elevation points (but not when the rest of points are added back)
!	control%SBLS_control%preconditioner = 4	! Super-slow
!	control%SBLS_control%preconditioner = 5	! Breaks
!	control%SBLS_control%preconditioner = 11	! Segfault in sbls_find_basis():uls_factorize():gls_analyse()
!	control%SBLS_control%preconditioner = 12	! Segfault
!	control%SBLS_control%preconditioner = 0



	call system_clock(time0_ms)
	CALL EQP_solve( p, data, control, inform) ! Solve
	call system_clock(time1_ms)
	delta_time = time1_ms - time0_ms
	delta_time = delta_time * 1d-3
	write(6, "( 'EQP_Solve took ', F6.3, ' seconds')") delta_time

!	inform%status = 0
	write(6, "( ' condition_number_1= ', F6.3, ' condition_number_2= ', F6.3)") &
		inform%SBLS_inform%SLS_inform%condition_number_1, &
		inform%SBLS_inform%SLS_inform%condition_number_2
	IF ( inform%status /= 0 ) THEN	! Error
		eqp_solve_simple = .false.
		WRITE( 6, "( ' EQP_solve exit status = ', I6 ) " ) inform%status
		select case(inform%status)
!			case(-10)
!				write(6, "( ' inform%sils_factorize_status = ', I6 ) " ) inform%sils_factorize_status
			case(-10)
				write(6, "( ' inform%SBLS_inform%status = ', I6 ) " ) inform%SBLS_inform%status
				write(6, "( ' inform%SBLS_inform%SLS_inform%status = ', I6 ) " ) inform%SBLS_inform%SLS_inform%status

!			case(-32)
!				write(6, "( ' inform%PRESOLVE_inform%status = ', I6 ) " ) inform%PRESOLVE_inform%status
!			case(-35)
!				write(6, "( ' inform%QPC_inform%status = ', I6 ) " ) inform%QPC_inform%status
		end select
	else
		eqp_solve_simple = .true.
		! ----------- Good return
		print *,'Successful QP Solve!'

		WRITE( 6, "( ' QP: ', I0, ' QPA iterations ', /, &
			' Optimal objective value =', &
			ES12.4  )" ) &
			inform%cg_iter, inform%obj

!		WRITE( 6, "( ' Optimal solution = ', ( 5ES12.4 ) )" ) p%X

	end if
	CALL EQP_terminate( data, control, inform) ! delete internal workspace
	call flush(6)

end function eqp_solve_simple
