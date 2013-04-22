module glint2_modele_types
use f90blitz

type glint2_modele_matrix
	! k = height point/class dimension
	integer, dimension(:), allocatable :: rows_i, rows_j, rows_k
	integer, dimension(:), allocatable :: cols_i, cols_j, cols_k
	real*8, dimension(:), allocatable :: vals
end type glint2_modele_matrix

type, bind(c) :: glint2_modele_matrix_f
	type(arr_spec_1) :: rows_i_f, rows_j_f, rows_k_f	! int
	type(arr_spec_1) :: cols_i_f, cols_j_f, cols_k_f	! int
	type(arr_spec_1) :: vals_f							! double
end type

end module

module glint2_modele
use f90blitz
use iso_c_binding
use glint2_modele_types
!use MpiSupport_mod
implicit none

! ================================================
! Stuff from glint2_modele.cpp

INTERFACE


	function glint2_modele_new( &
		maker_fname_f, maker_fname_len, &
		maker_vname_f, maker_vname_len, &
		im,jm, &
		i0h,i1h,j0h,j1h, &
		i0,i1,j0,j1, &
		j0s,j1s, &
		comm_f, root) bind(c)
	use iso_c_binding
		character(c_char) :: maker_fname_f(*)
		integer(c_int), value :: maker_fname_len
		character(c_char) :: maker_vname_f(*)
		integer(c_int), value :: maker_vname_len
		integer(c_int), value :: im, jm
		integer(c_int), value :: i0h,i1h,j0h,j1h
		integer(c_int), value :: i0,i1,j0,j1
		integer(c_int), value :: j0s,j1s
		integer(c_int), value :: comm_f, root
		type(c_ptr) :: glint2_modele_new
	end function

	subroutine glint2_modele_delete(api) bind(c)
		use iso_c_binding
		use f90blitz
		type(c_ptr) :: api		! NOT VALUE here.
	end subroutine

	function glint2_modele_nhc(api) bind(c)
		use iso_c_binding
		use f90blitz
		type(c_ptr), value :: api
		integer(c_int) :: glint2_modele_nhc
	end function

	subroutine glint2_modele_compute_fgice_c(api, &
		replace_fgice_b, fgice1_glint2_f, &
		fgice1_f, fgrnd1_f, focean1_f, flake1_f) bind(c)
	use iso_c_binding
	use f90blitz
		type(c_ptr), value :: api
		integer(c_int), value :: replace_fgice_b
		type(arr_spec_2) :: fgice1_glint2_f
		type(arr_spec_2) :: fgice1_f
		type(arr_spec_2) :: fgrnd1_f
		type(arr_spec_2) :: focean1_f
		type(arr_spec_2) :: flake1_f
	end subroutine

	function glint2_modele_init_landice_com_part1(api) bind(c)
		use iso_c_binding
		type(c_ptr), value :: api
		integer(c_int) :: glint2_modele_init_landice_com_part1
	end function

	subroutine glint2_modele_init_landice_com_part2(api, &
		zatmo1_f, BYGRAV, fgice1_glint2_f, fgice1_f, &
		used1h_f, fhc1h_f, elevhp_f, hp_to_hc_f, fhp_approx1h_f) bind(c)
	use iso_c_binding
	use f90blitz
	use glint2_modele_types
		type(c_ptr), value :: api
		type(arr_spec_2) :: zatmo1_f
		real*8, value :: BYGRAV
		type(arr_spec_2) :: fgice1_glint2_f
		type(arr_spec_2) :: fgice1_f
		type(arr_spec_3) :: used1h_f
		type(arr_spec_3) :: fhc1h_f
		type(arr_spec_3) :: elevhp_f
		type(glint2_modele_matrix_f) :: hp_to_hc_f
		type(arr_spec_3) :: fhp_approx1h_f
	end subroutine

	subroutine glint2_modele_couple_to_ice(api, smb1hp_f, seb1hp_f) bind(c)
	use iso_c_binding
	use f90blitz
		type(c_ptr), value :: api
		type(arr_spec_3) :: smb1hp_f, seb1hp_f
	end subroutine

END INTERFACE

!include 'mpif.h'

! ================================================


contains



subroutine glint2_modele_compute_fgice(api, replace_fgice, &
	fgice1_glint2, &
	fgice1, fgrnd1, focean1, flake1, i0h, j0h)
type(c_ptr), value :: api
integer :: i0h, j0h
logical :: replace_fgice
real*8, dimension(i0h:,j0h:) :: fgice1_glint2			! OUT
real*8, dimension(i0h:,j0h:) :: fgice1					! OUT
real*8, dimension(i0h:,j0h:) :: fgrnd1, focean1, flake1	! INOUT

	! ----------
	integer :: replace_fgice_b
	type(arr_spec_2) :: fgice1_f, fgrnd1_f, focean1_f, flake1_f
	type(arr_spec_2) :: fgice1_glint2_f

	replace_fgice_b = 0
	if (replace_fgice) replace_fgice_b = 1

	! Grab array descriptors
	call get_spec_double_2(fgice1_glint2, i0h, j0h, fgice1_glint2_f)
	call get_spec_double_2(fgice1, i0h, j0h, fgice1_f)
	call get_spec_double_2(fgrnd1, i0h, j0h, fgrnd1_f)
	call get_spec_double_2(focean1, i0h, j0h, focean1_f)
	call get_spec_double_2(flake1, i0h, j0h, flake1_f)

	! Call the C-side of the interface
	call glint2_modele_compute_fgice_c(api, replace_fgice_b, &
		fgice1_glint2_f, fgice1_f, fgrnd1_f, focean1_f, flake1_f)
end subroutine

subroutine glint2_modele_init_landice_com(api, &
	zatmo1, BYGRAV, fgice1_glint2, fgice1, &
	used1h, fhc1h, elevhp, hp_to_hc, fhp_approx1h, &
	i0h, j0h)
type(c_ptr), value :: api
integer :: i0h, j0h
real*8, dimension(i0h:,j0h:) :: zatmo1
real*8, INTENT(IN) :: BYGRAV
real*8, dimension(i0h:,j0h:) :: fgice1_glint2, fgice1
integer, dimension(i0h:,j0h:,:) :: used1h					! OUT
real*8, dimension(i0h:,j0h:,:) :: fhc1h					! OUT
real*8, dimension(i0h:,j0h:,:) :: elevhp
type(glint2_modele_matrix) :: hp_to_hc
real*8, dimension(i0h:,j0h:,:) :: fhp_approx1h


	integer :: n

	! ------------------- local vars
	type(glint2_modele_matrix_f) :: mat_f
	type(arr_spec_2) :: zatmo1_f
	type(arr_spec_2) :: fgice1_glint2_f, fgice1_f
	type(arr_spec_3) :: used1h_f
	type(arr_spec_3) :: fhc1h_f
	type(arr_spec_3) :: elevhp_f
	type(glint2_modele_matrix_f) :: hp_to_hc_f
	type(arr_spec_3) :: fhp_approx1h_f

	! ------------------- subroutine body

	! -------- Part 1: Figure out how big we must make the arrays
	n = glint2_modele_init_landice_com_part1(api)

	! -------- Allocate those arrays
	allocate(hp_to_hc%rows_i(n))
	allocate(hp_to_hc%rows_j(n))
	allocate(hp_to_hc%rows_k(n))
	allocate(hp_to_hc%cols_i(n))
	allocate(hp_to_hc%cols_j(n))
	allocate(hp_to_hc%cols_k(n))
	allocate(hp_to_hc%vals(n))

	! -------- Grab array descriptors from arrays we just allocated
	call get_spec_int_1(hp_to_hc%rows_i, 1, hp_to_hc_f%rows_i_f)
	call get_spec_int_1(hp_to_hc%rows_j, 1, hp_to_hc_f%rows_j_f)
	call get_spec_int_1(hp_to_hc%rows_k, 1, hp_to_hc_f%rows_k_f)
	call get_spec_int_1(hp_to_hc%cols_i, 1, hp_to_hc_f%cols_i_f)
	call get_spec_int_1(hp_to_hc%cols_j, 1, hp_to_hc_f%cols_j_f)
	call get_spec_int_1(hp_to_hc%cols_k, 1, hp_to_hc_f%cols_k_f)
	call get_spec_double_1(hp_to_hc%vals, 1, hp_to_hc_f%vals_f)


	! Grab array descriptors
	call get_spec_double_2(zatmo1, i0h, j0h, zatmo1_f)
	call get_spec_double_2(fgice1_glint2, i0h, j0h, fgice1_glint2_f)
	call get_spec_double_2(fgice1, i0h, j0h, fgice1_f)
	call get_spec_int_3(used1h, i0h, j0h, 1, used1h_f)
	call get_spec_double_3(fhc1h, i0h, j0h, 1, fhc1h_f)
	call get_spec_double_3(elevhp, i0h, j0h, 1, elevhp_f)
	call get_spec_double_3(fhp_approx1h, i0h, j0h, 1, fhp_approx1h_f)

	! Call the C-side of the interface
	call glint2_modele_init_landice_com_part2(api, &
		zatmo1_f, BYGRAV, fgice1_glint2_f, fgice1_f, &
		used1h_f, fhc1h_f, elevhp_f, hp_to_hc_f, fhp_approx1h_f)

print *,'END glint2_modele_init_landice_com()'
end subroutine


end module
