module glint2_modele
use f90blitz
use iso_c_binding
!use MpiSupport_mod
implicit none

! ================================================
! Stuff from glint2_modele.cpp

type glint2_modele_matrix
	! k = height point/class dimension
	integer, dimension(:), allocatable :: rows_i, rows_j, rows_k
	integer, dimension(:), allocatable :: cols_i, cols_j, cols_k
	real*8, dimension(:), allocatable :: vals
end type glint2_modele_matrix

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

	subroutine glint2_modele_get_elevhc_c(api, elevhc)
		use iso_c_binding
		use f90blitz
		type(c_ptr), value :: api
		type(arr_spec_3) :: elevhc
	end subroutine


	subroutine glint2_modele_compute_fgice_c(api, fgice1_f, fgrnd1_f, focean1_f, flake1_f) bind(c)
		use iso_c_binding
		use f90blitz
		type(c_ptr), value :: api
		type(arr_spec_2) :: fgice1_f
		type(arr_spec_2) :: fgrnd1_f
		type(arr_spec_2) :: focean1_f
		type(arr_spec_2) :: flake1_f
	end subroutine

	subroutine glint2_modele_compute_fhc_c(api, fhc1h_f) bind(c)
		use iso_c_binding
		use f90blitz
		type(c_ptr), value :: api
		type(arr_spec_3) :: fhc1h_f
	end subroutine

	function glint2_modele_hp_to_hc_part1(api) bind(c)
		use iso_c_binding
		type(c_ptr), value :: api
		integer(c_int) :: glint2_modele_hp_to_hc_part1
	end function

	subroutine glint2_modele_hp_to_hc_part2(api, &
		rows_i_f, rows_j_f, rows_k_f, &
		cols_i_f, cols_j_f, cols_k_f, &
		vals_f) bind(c)
	use iso_c_binding
	use f90blitz
		type(c_ptr), value :: api
		type(arr_spec_1) :: rows_i_f, rows_j_f, rows_k_f
		type(arr_spec_1) :: cols_i_f, cols_j_f, cols_k_f
		type(arr_spec_1) :: vals_f
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

subroutine glint2_modele_get_elevhc(api, elevhc, i0h, j0h)
type(c_ptr), value :: api
integer :: i0h, j0h
real*8, dimension(i0h:,j0h:,:) :: elevhc

	! ----------

	type(arr_spec_3) :: elevhc_f

	! Grab array descriptors
	call get_spec_double_3(elevhc, i0h, j0h, 1, elevhc_f)

	! Call the C-side of the interface
	call glint2_modele_get_elevhc_c(api, elevhc_f)
end subroutine


subroutine glint2_modele_compute_fgice(api, fgice1, fgrnd1, focean1, flake1, i0h, j0h)
type(c_ptr), value :: api
integer :: i0h, j0h
real*8, dimension(i0h:,j0h:) :: fgice1					! OUT
real*8, dimension(i0h:,j0h:) :: fgrnd1, focean1, flake1	! INOUT

	! ----------

	type(arr_spec_2) :: fgice1_f, fgrnd1_f, focean1_f, flake1_f

	! Grab array descriptors
	call get_spec_double_2(fgice1, i0h, j0h, fgice1_f)
	call get_spec_double_2(fgrnd1, i0h, j0h, fgrnd1_f)
	call get_spec_double_2(focean1, i0h, j0h, focean1_f)
	call get_spec_double_2(flake1, i0h, j0h, flake1_f)

	! Call the C-side of the interface
	call glint2_modele_compute_fice_c(api, fgice1_f, fgrnd1_f, focean1_f, flake1_f)
end subroutine

subroutine glint2_modele_compute_fhc(api, fhc1h, i0h, j0h)
type(c_ptr), value :: api
integer :: i0h, j0h
real*8, dimension(i0h:,j0h:,:) :: fhc1h					! OUT

	! ----------

	type(arr_spec_3) :: fhc1h_f

	! Grab array descriptors
	call get_spec_double_3(fhc1h, i0h, j0h, 1, fhc1h_f)

	! Call the C-side of the interface
	call glint2_modele_compute_fhc_c(api, fhc1h_f)
end subroutine

subroutine glint2_modele_hp_to_hc(api, mat)
implicit none
	type(c_ptr), value :: api
	type(glint2_modele_matrix) :: mat
	integer :: n

	! ------------------- local vars
	type(arr_spec_1) :: rows_i_f, rows_j_f, rows_k_f
	type(arr_spec_1) :: cols_i_f, cols_j_f, cols_k_f
	type(arr_spec_1) :: vals_f

	! ------------------- subroutine body

	! -------- Part 1: Figure out how big we must make the arrays
	n = glint2_modele_hp_to_hc_part1(api)

	! -------- Allocate those arrays
	allocate(mat%rows_i(n))
	allocate(mat%rows_j(n))
	allocate(mat%rows_k(n))
	allocate(mat%cols_i(n))
	allocate(mat%cols_j(n))
	allocate(mat%cols_k(n))
	allocate(mat%vals(n))

	! -------- Grab array descriptors from arrays we just allocated
	call get_spec_int_1(mat%rows_i, 1, rows_i_f)
	call get_spec_int_1(mat%rows_j, 1, rows_j_f)
	call get_spec_int_1(mat%rows_j, 1, rows_k_f)
	call get_spec_int_1(mat%cols_i, 1, cols_i_f)
	call get_spec_int_1(mat%cols_j, 1, cols_j_f)
	call get_spec_int_1(mat%cols_j, 1, cols_k_f)
	call get_spec_double_1(mat%vals, 1, vals_f)

	! -------- Part 2: Fill the arrays we allocated
	call glint2_modele_hp_to_hc_part2(api, &
		rows_i_f, rows_j_f, rows_k_f, &
		cols_i_f, cols_j_f, cols_k_f, &
		vals_f)

end subroutine

end module
