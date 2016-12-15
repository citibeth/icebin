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
implicit none





! Parameters read out of the ModelE rundeck and sent to IceBin
! These are later incorporated in gcmce_new()
integer, parameter :: MAX_CHAR_LEN = 128   ! From ModelE's Dictionary_mod.F90
type, bind(c) :: ModelEParams

    ! Segment specs: to be further parsed.
    character(MAX_CHAR_LEN, kind=c_char) :: icebin_segments = 'legacy,sealand,ec'
    real(c_double) :: dtsrc
    integer(c_int) :: yeari
end type ModelEParams



! ================================================
! Stuff from icebin_modele.cpp

INTERFACE

    ! ------------------------------------------------------
    ! ------------- LISheetIceBin::allocate()


    ! Called from lisheeticebin%allocate()
    function gcmce_new( &
        rdparams, &
        im,jm, &
        i0,i1,j0,j1, &
        comm_f, root) bind(c)
    use iso_c_binding
    import ModelEParams
        type(c_ptr) :: gcmce_new
        type(ModelEParams) :: rdparams
        integer(c_int), value :: im, jm
        integer(c_int), value :: i0,i1,j0,j1
        integer(c_int), value :: comm_f, root
    end function gcmce_new


    ! Called from lisheeticebin%allocate() (via setup_gcm_inputs)
    function gcmce_read_nhc_gcm(api) bind(c)
    use iso_c_binding
        type(c_ptr), value :: api
        integer(c_int) :: gcmce_read_nhc_gcm
    end function gcmce_read_nhc_gcm

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
    end subroutine gcmce_set_const




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
    end function gcmce_add_gcm_outputE


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
    end function gcmce_add_gcm_inputA

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
    end function gcmce_add_gcm_inputE

    subroutine gcmce_reference_globals(api, &
        fhc, elevE, &
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

END INTERFACE

end module
