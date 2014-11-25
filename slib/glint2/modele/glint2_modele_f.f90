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

module glint2_modele
use f90blitz
use iso_c_binding
!use MpiSupport_mod
implicit none

! ================================================
! Stuff from glint2_modele.cpp

INTERFACE

	function new_glint2_modele() bind(c)
	use iso_c_binding
		type(c_ptr) :: new_glint2_modele
	end function

	subroutine glint2_modele_set_const(api, &
		name_f, name_len, &
		val, &
		units_f, units_len, &
		description_f, description_len) bind(c)
	use iso_c_binding
		type(c_ptr), value :: api
		character(c_char) :: name_f(*)
		integer(c_int), value :: name_len
		real(c_double), value :: val
		character(c_char) :: units_f(*)
		integer(c_int), value :: units_len
		character(c_char) :: description_f(*)
		integer(c_int), value :: description_len
	end subroutine

	function glint2_modele_add_gcm_input(api, &
		field_f, field_len, &
		units_f, units_len, &
        grid_f, grid_len, &
		long_name_f, long_name_len) bind(c)
	use iso_c_binding
		type(c_ptr), value :: api
		character(c_char) :: field_f(*)
		integer(c_int), value :: field_len
		character(c_char) :: units_f(*)
		integer(c_int), value :: units_len
		character(c_char) :: grid_f(*)
		integer(c_int), value :: grid_len
		character(c_char) :: long_name_f(*)
		integer(c_int), value :: long_name_len
		integer(c_int) :: glint2_modele_add_gcm_input
	end function

	subroutine glint2_modele_init0(api, &
		maker_fname_f, maker_fname_len, &
		maker_vname_f, maker_vname_len, &
		im,jm, &
		i0h,i1h,j0h,j1h, &
		i0,i1,j0,j1, &
		j0s,j1s, &
		comm_f, root, &
        write_constants) bind(c)
	use iso_c_binding
		type(c_ptr), value :: api
		character(c_char) :: maker_fname_f(*)
		integer(c_int), value :: maker_fname_len
		character(c_char) :: maker_vname_f(*)
		integer(c_int), value :: maker_vname_len
		integer(c_int), value :: im, jm
		integer(c_int), value :: i0h,i1h,j0h,j1h
		integer(c_int), value :: i0,i1,j0,j1
		integer(c_int), value :: j0s,j1s
		integer(c_int), value :: comm_f, root
		integer(c_int), value :: write_constants
	end subroutine

	subroutine glint2_modele_delete(api) bind(c)
		use iso_c_binding
		use f90blitz
		type(c_ptr) :: api		! NOT VALUE here.
	end subroutine

	function glint2_modele_nhp(api) bind(c)
		use iso_c_binding
		use f90blitz
		type(c_ptr), value :: api
		integer(c_int) :: glint2_modele_nhp
	end function

	subroutine glint2_modele_set_start_time(api, iyear1, itimei, dtsrc) bind(c)
		use iso_c_binding
		type(c_ptr), value :: api
		integer(c_int), value :: iyear1
		integer(c_int), value :: itimei
		real(c_double), value :: dtsrc
	end subroutine

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

	subroutine glint2_modele_init_landice_com_c(api, &
		zatmo1_f, BYGRAV, fgice1_glint2_f, fgice1_f, &
		used1h_f, fhc1h_f, elev1h_f, &
		i0, j0, i1, j1) bind(c)
	use iso_c_binding
	use f90blitz
		type(c_ptr), value :: api
		type(arr_spec_2) :: zatmo1_f
		real(c_double), value :: BYGRAV
		type(arr_spec_2) :: fgice1_glint2_f
		type(arr_spec_2) :: fgice1_f
		type(arr_spec_3) :: used1h_f
		type(arr_spec_3) :: fhc1h_f
		type(arr_spec_3) :: elev1h_f
		integer(kind=c_int) :: i0, j0, i1, j1
	end subroutine

	subroutine glint2_modele_init_hp_to_ices(api) bind(c)
	use iso_c_binding
	use f90blitz
		type(c_ptr), value :: api
	end subroutine

	subroutine glint2_modele_couple_to_ice_c(api, itime, smb1hp_f, seb1hp_f, tg21hp_f) bind(c)
	use iso_c_binding
	use f90blitz
		type(c_ptr), value :: api
		integer(c_int), value :: itime
		type(arr_spec_3) :: smb1hp_f, seb1hp_f, tg21hp_f
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
	used1h, fhc1h, elev1h, &
	i0h, j0h, &
	i0, i1, j0, j1)
type(c_ptr), value :: api
integer, value :: i0h, j0h
integer, value :: i0, i1, j0, j1
real*8, dimension(i0h:,j0h:) :: zatmo1
real*8, INTENT(IN) :: BYGRAV
real*8, dimension(i0h:,j0h:) :: fgice1_glint2, fgice1
integer, dimension(i0h:,j0h:,:) :: used1h					! OUT
real*8, dimension(i0h:,j0h:,:) :: fhc1h					! OUT
real*8, dimension(i0h:,j0h:,:) :: elev1h

	integer :: n

	! ------------------- local vars
	type(arr_spec_2) :: zatmo1_f
	type(arr_spec_2) :: fgice1_glint2_f, fgice1_f
	type(arr_spec_3) :: used1h_f
	type(arr_spec_3) :: fhc1h_f
	type(arr_spec_3) :: elev1h_f

	! ------------------- subroutine body


	! Grab array descriptors
	call get_spec_double_2(zatmo1, i0h, j0h, zatmo1_f)
	call get_spec_double_2(fgice1_glint2, i0h, j0h, fgice1_glint2_f)
	call get_spec_double_2(fgice1, i0h, j0h, fgice1_f)
	call get_spec_int_3(used1h, i0h, j0h, 1, used1h_f)
	call get_spec_double_3(fhc1h, i0h, j0h, 1, fhc1h_f)
	call get_spec_double_3(elev1h, i0h, j0h, 1, elev1h_f)

	! Call the C-side of the interface
	call glint2_modele_init_landice_com_c(api, &
		zatmo1_f, BYGRAV, fgice1_glint2_f, fgice1_f, &
		used1h_f, fhc1h_f, elev1h_f, &
		i0, i1, j0, j1)

print *,'END glint2_modele_init_landice_com()'
end subroutine

subroutine glint2_modele_couple_to_ice(api, &
	itime, smb1h, seb1h, tg21h, &
	i0h, j0h)
type(c_ptr), value :: api
integer :: i0h, j0h
real*8, dimension(i0h:,j0h:,:) :: smb1h, seb1h, tg21h
integer, value :: itime

	integer :: n

	! ------------------- local vars
	type(arr_spec_3) :: smb1h_f, seb1h_f, tg21h_f, gcm_inputs_d_f

	! ------------------- subroutine body
print *,'BEGIN glint2_modele_couple_to_ice()'

	! Grab array descriptors
	call get_spec_double_3(smb1h, i0h, j0h, 1, smb1h_f)		! kg/m^2
	call get_spec_double_3(seb1h, i0h, j0h, 1, seb1h_f)		! J/m^2: Latent Heat
	call get_spec_double_3(tg21h, i0h, j0h, 1, tg21h_f)		! C
	call get_spec_double_3(gcm_inputs_d, 1,1,1, gcm_inputs_d_f)

	! Call the C-side of the interface
	call glint2_modele_couple_to_ice_c(api, itime, smb1h_f, seb1h_f, tg21h_f, gcm_inputs_d_f)

print *,'END glint2_modele_couple_to_ice()'
end subroutine


! ! Go from a VectorSparseVector<int, double> output from GCMCoupler, to
! ! a fully scattered array
! subroutine VectorSparseVector_to_Scattered2(grid, vsv_indices, vsv_vals, nvals, scat) bind(c)
! type(dist_grid), intent(in) :: grid
! integer(c_int) :: vsv_indices(*)
! real(c_double) :: vsv_vals(*)
! integer(c_int), value :: nvals
! real*8, dimension(:,:) :: scat
!
! end subroutine VectorSparseVector_to_Scattered2


end module
