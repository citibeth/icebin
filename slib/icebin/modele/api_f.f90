! IceBin: A Coupling Library for Ice Models and GCMs
! Copyright (c) 2013-2016 by Elizabeth Fischer
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published
! by the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module icebin_modele
use icebin_f90blitz
use iso_c_binding
!use MpiSupport_mod
implicit none

! Parameters read out of the ModelE rundeck and sent to IceBin
! These are later incorporated in gcmce_new()
type ModelEParams

    ! Segment specs: to be further parsed.
    character(c_char)*(MAX_CHAR_LEN) :: icebin_segments = 'legacy,sealand,ec'
    real(c_double) :: dtsrc
    integer(c_int) :: yeari
end type ModelEParams



! ================================================
! Stuff from icebin_modele.cpp

INTERFACE

    ! ------------------------------------------------------
    ! ------------- LISheetIceBin::allocate()


    ! Called from lisheeticebin%allocate()
    function gcmce_new(
        rdparams, &
        im,jm, &
        i0,i1,j0,j1, &
        comm_f, root) bind(c)
    use iso_c_binding
        type(c_ptr) :: gcmce_new
        type(ModelEParams) :: rdparams
        integer(c_int), value :: im, jm
        integer(c_int), value :: i0h,i1h,j0h,j1h
        integer(c_int), value :: i0,i1,j0,j1
        integer(c_int), value :: j0s,j1s
        integer(c_int), value :: comm_f, root
    end subroutine


    ! Called from lisheeticebin%allocate() (via setup_gcm_inputs)
    function gcmce_read_nhc_gcm(api) bind(c)
    use iso_c_binding
        type(c_ptr), value :: api
        integer(c_int) :: gcmce_read_nhc_gcm
    end function gcmce_gcm_inputs_nhp

    ! Called from lisheeticebin%allocate()
    subroutine gcmce_set_const(api, &
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




    function new_gcmce_c() result(ret) bind(c)
    use iso_c_binding
        type(c_ptr) :: ret
    end function

    function gcmce_add_gcm_outputE(api, &
        field_f, field_len, &
        units_f, units_len, &
        grid_f, grid_len, &
        long_name_f, long_name_len) bind(c)
    use iso_c_binding
    use icebin_f90blitz
        integer(c_int) :: gcmce_add_gcm_outputE
        type(c_ptr), value :: api
        type(arr_spec_3) :: var_f
        character(c_char) :: field_f(*)
        integer(c_int), value :: field_len
        character(c_char) :: units_f(*)
        integer(c_int), value :: units_len
        character(c_char) :: grid_f(*)
        integer(c_int), value :: grid_len
        character(c_char) :: long_name_f(*)
        integer(c_int), value :: long_name_len
    end function


    ! Called from lisheeticebin%allocate() (via setup_gcm_inputs())
    function gcmce_add_gcm_inputA(api, &
        field_f, field_len, &
        units_f, units_len, &
        grid_f, grid_len, &
        initial, &
        long_name_f, long_name_len) bind(c)
    use iso_c_binding
    use icebin_f90blitz
        integer(c_int) :: gcmce_add_gcm_inputA
        type(c_ptr), value :: api
        type(arr_spec_2) :: var_f
        character(c_char) :: field_f(*)
        integer(c_int), value :: field_len
        character(c_char) :: units_f(*)
        integer(c_int), value :: units_len
        character(c_char) :: grid_f(*)
        integer(c_int), value :: grid_len
        character(c_char) :: long_name_f(*)
        integer(c_int), value :: initial
        integer(c_int), value :: long_name_len
    end function

    ! Called from lisheeticebin%allocate() (via setup_gcm_inputs())
    function gcmce_add_gcm_inputE(api, &
        field_f, field_len, &
        units_f, units_len, &
        grid_f, grid_len, &
        initial, &
        long_name_f, long_name_len) bind(c)
    use iso_c_binding
    use icebin_f90blitz
        integer(c_int) :: gcmce_add_gcm_inputE
        type(c_ptr), value :: api
        type(arr_spec_3) :: var_f
        character(c_char) :: field_f(*)
        integer(c_int), value :: field_len
        character(c_char) :: units_f(*)
        integer(c_int), value :: units_len
        character(c_char) :: grid_f(*)
        integer(c_int), value :: grid_len
        character(c_char) :: long_name_f(*)
        integer(c_int), value :: initial
        integer(c_int), value :: long_name_len
    end function

    subroutine gcmce_reference_globals(api, &
        fhc, elevE,
        focean, flake, fgrnd, fgice, zatmo)
    use iso_c_binding
    use icebin_f90blitz
        type(c_ptr), value :: api
        type(arr_spec_3) :: fhc, elevE
        type(arr_spec_2) :: focean, flake, fgrnd, fgice, zatmo
    end subroutine gcmce_reference_globals

    subroutine gcmce_io_rsf(api, &
        fname_f, fname_len)
    use iso_c_binding
        type(c_ptr), value :: api
        character(c_char) :: fname_f(*)
        integer(c_int), value :: fname_len
    end subroutine gcmce_io_rsf

    subroutine gcmce_cold_start(api, itimei, dtsrc) bind(c)
        use iso_c_binding
        type(c_ptr), value :: api
        integer(c_int), value :: itimei
        real(c_double), value :: dtsrc
    end subroutine

    subroutine gcmce_couple_native(api, itime, &
        gcm_inputs_d_f) bind(c)
    use iso_c_binding
    use icebin_f90blitz
        type(c_ptr), value :: api
        integer(c_int), value :: itime
        type(arr_spec_3) :: massxfer_f, enthxfer_f, deltah_f
        type(arr_spec_3) :: gcm_inputs_d_f
    end subroutine




    ! ------------------------------------------------------







    subroutine gcmce_delete(api) bind(c)
        use iso_c_binding
        use icebin_f90blitz
        type(c_ptr) :: api      ! NOT VALUE here.
    end subroutine


    ! -------------------------------------------
    subroutine gcmce_init_hp_to_ices(api) bind(c)
    use iso_c_binding
    use icebin_f90blitz
        type(c_ptr), value :: api
    end subroutine

    subroutine gcmce_get_initial_state_c(api, itime, gcm_inputs_d_f) bind(c)
    use iso_c_binding
    use icebin_f90blitz
        type(c_ptr), value :: api
        integer(c_int), value :: itime
        type(arr_spec_3) :: gcm_inputs_d_f
    end subroutine

    subroutine gcmce_set_gcm_output_c(api, &
        field_name_f, field_name_len, arr_f) bind(c)
    use iso_c_binding
    use icebin_f90blitz
        type(c_ptr), value :: api
        character(c_char) :: field_name_f(*)
        integer(c_int), value :: field_name_len
        type(arr_spec_3) :: arr_f
    end subroutine

END INTERFACE

!include 'mpif.h'

! ================================================


contains

! ---------------------------------------------------
subroutine gcmce_set_gcm_output(api, field_name, arr, i0, j0, k0)
    type(c_ptr), value :: api
    character(*) :: field_name
    real*8, dimension(:,:,:), target :: arr
    integer :: i0,j0,k0
    ! --------- Locals
    type(arr_spec_3) :: arr_f

    call get_spec_double_3(arr, i0,j0,k0, arr_f)
    call gcmce_set_gcm_output_c(api, &
        field_name, len(field_name), arr_f)

end subroutine gcmce_set_gcm_output


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
