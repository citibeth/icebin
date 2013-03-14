module modele_api
use f90blitz
use iso_c_binding
implicit none

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
	rows_i, rows_j, cols_i, cols_j, vals)
implicit none
	type(c_ptr), value :: api
	integer, dimension(:), allocatable :: rows_i, rows_j
	integer, dimension(:), allocatable :: cols_i, cols_j
	real*8, dimension(:), allocatable :: vals
	integer :: n

	! ------------------- local vars
	type(arr_spec_1) :: rows_i_f, rows_j_f
	type(arr_spec_1) :: cols_i_f, cols_j_f
	type(arr_spec_1) :: vals_f

	! ------------------- subroutine body

	! -------- Part 1: Figure out how big we must make the arrays
	call modele_api_hp_to_hc_part1(api, n)

	! -------- Allocate those arrays
	allocate(rows_i(n))
	allocate(rows_j(n))
	allocate(cols_i(n))
	allocate(cols_j(n))
	allocate(vals(n))

	! -------- Grab array descriptors from arrays we just allocated
	call get_spec_int_1(rows_i, 1, rows_i_f)
	call get_spec_int_1(rows_j, 1, rows_j_f)
	call get_spec_int_1(cols_i, 1, cols_i_f)
	call get_spec_int_1(cols_j, 1, cols_j_f)
	call get_spec_double_1(vals, 1, vals_f)

	! -------- Part 2: Fill the arrays we allocated
	call modele_api_hp_to_hc_part2(api, &
		rows_i_f, rows_j_f, cols_i_f, cols_j_f, vals_f)

end subroutine

end module
