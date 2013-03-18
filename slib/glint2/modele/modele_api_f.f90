module modele_api
use f90blitz
use iso_c_binding
implicit none

! ================================================
! Stuff from modele_api.cpp

INTERFACE

	function modele_api_new( &
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
		type(c_ptr) :: modele_api_new
	end function

	subroutine modele_api_delete(api) bind(c)
		use iso_c_binding
		use f90blitz
		type(c_ptr) :: api		! NOT VALUE here.
	end subroutine

	subroutine modele_api_compute_fhc_c(api, fhc1h_f, fgice1_f) bind(c)
		use iso_c_binding
		use f90blitz
		type(c_ptr), value :: api
		type(arr_spec_3) :: fhc1h_f
		type(arr_spec_2) :: fgice1_f
	end subroutine

	function modele_api_hp_to_hc_part1(api) bind(c)
		use iso_c_binding
		type(c_ptr), value :: api
		integer(c_int) :: modele_api_hp_to_hc_part1
	end function

	subroutine modele_api_hp_to_hc_part2(api, &
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

	subroutine modele_api_couple_to_ice(api, smb1hp_f, seb1hp_f) bind(c)
	use iso_c_binding
	use f90blitz
		type(c_ptr), value :: api
		type(arr_spec_3) :: smb1hp_f, seb1hp_f
	end subroutine

END INTERFACE

! ================================================


contains

subroutine modele_api_compute_fhc(api, fhc1h, fgice1, i0h, j0h)
type(c_ptr), value :: api
integer :: i0h, j0h
real*8, dimension(i0h:,j0h:,:) :: fhc1h
real*8, dimension(i0h:,j0h:) :: fgice1

	! ----------

	type(arr_spec_3) :: fhc1h_f
	type(arr_spec_2) :: fgice1_f

	! Grab array descriptors
	call get_spec_double_3(fhc1h, i0h, j0h, 1, fhc1h_f)
	call get_spec_double_2(fgice1, i0h, j0h, fgice1_f)

	! Call the C-side of the interface
	call modele_api_compute_fhc_c(api, fhc1h_f, fgice1_f)
end subroutine

subroutine modele_api_hp_to_hc(api, &
	rows_i, rows_j, rows_k, &		! k = height point/class dimension
	cols_i, cols_j, cols_k, vals)
implicit none
	type(c_ptr), value :: api
	integer, dimension(:), allocatable :: rows_i, rows_j, rows_k
	integer, dimension(:), allocatable :: cols_i, cols_j, cols_k
	real*8, dimension(:), allocatable :: vals
	integer :: n

	! ------------------- local vars
	type(arr_spec_1) :: rows_i_f, rows_j_f, rows_k_f
	type(arr_spec_1) :: cols_i_f, cols_j_f, cols_k_f
	type(arr_spec_1) :: vals_f

	! ------------------- subroutine body

	! -------- Part 1: Figure out how big we must make the arrays
	n = modele_api_hp_to_hc_part1(api)

	! -------- Allocate those arrays
	allocate(rows_i(n))
	allocate(rows_j(n))
	allocate(rows_k(n))
	allocate(cols_i(n))
	allocate(cols_j(n))
	allocate(cols_k(n))
	allocate(vals(n))

	! -------- Grab array descriptors from arrays we just allocated
	call get_spec_int_1(rows_i, 1, rows_i_f)
	call get_spec_int_1(rows_j, 1, rows_j_f)
	call get_spec_int_1(rows_j, 1, rows_k_f)
	call get_spec_int_1(cols_i, 1, cols_i_f)
	call get_spec_int_1(cols_j, 1, cols_j_f)
	call get_spec_int_1(cols_j, 1, cols_k_f)
	call get_spec_double_1(vals, 1, vals_f)

	! -------- Part 2: Fill the arrays we allocated
	call modele_api_hp_to_hc_part2(api, &
		rows_i_f, rows_j_f, rows_k_f, &
		cols_i_f, cols_j_f, cols_k_f, &
		vals_f)

end subroutine

end module
