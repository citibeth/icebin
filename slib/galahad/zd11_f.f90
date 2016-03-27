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

! COPYRIGHT (c) 2006 Council for the Central Laboratory
!                    of the Research Councils
! Original date 21 February 2006. Version 1.0.0.
! 6 March 2007 Version 1.1.0. Argument stat made non-optional

! ===============================================================
! ===============================================================
! ===============================================================
! !> Sparse Matrix Representation Types
! MODULE HSL_ZD11D
! 
! !  ==========================
! !  Sparse matrix derived type
! !  ==========================
! 
!   !> Sparse matrix representation
!   !! Originally derived from HSL library (via GALAHAD; same as SMT_TYPE in GALAHAD).
!   TYPE, PUBLIC :: zd11_type
!     INTEGER :: m  !< Number of rows in this matrix.
!   INTEGER :: n    !< Number of columns in this matrix.
!   INTEGER :: ne   !< Number of non-zero elements in this matrix.
! 
!   !> "Name" of this matrix, not really needed.
!     CHARACTER, ALLOCATABLE, DIMENSION(:) :: id
! 
!   !> Indicates the storage scheme used.
!   !! It is set equal to either
!   !! 'DENSE', 'COORDINATE', 'SPARSE BY ROWS' or 'DIAGONAL'.  Use ZD11_put()
!   !! to set this variable.
!     CHARACTER, ALLOCATABLE, DIMENSION(:) :: type
! 
!   !> If type == 'COORDINATE', then this holds the row indices of
!   !! the non-zero elements (index base = 1).  Otherwise, it is not needed.
!     INTEGER, ALLOCATABLE, DIMENSION(:) :: row
! 
!   !> If type == 'COORDINATE' or 'SPARSE BY ROWS', this holds the column indices
!   !! of the non-zero elements (index base = 1).  Otherwise, it is not needed.
!     INTEGER, ALLOCATABLE, DIMENSION(:) :: col
! 
!   !> If type == 'SPARSE BY ROWS', this holds the starting position
!   !! of each row in the arrays col and val (index base = 1)
!     INTEGER, ALLOCATABLE, DIMENSION(:) :: ptr
! 
!   !> The value of each non-zero element in the matrix.
!     REAL ( KIND( 1.0D+0 ) ), ALLOCATABLE, DIMENSION(:) :: val
!   END TYPE
! 
! CONTAINS
! 
!    SUBROUTINE ZD11_put(array,string,stat)
!      CHARACTER, allocatable :: array(:)
!      CHARACTER(*), intent(in) ::  string
!      INTEGER, intent(OUT) ::  stat
! 
!      INTEGER :: i,l
! 
!      l = len_trim(string)
!      if (allocated(array)) then
!         deallocate(array,stat=stat)
!         if (stat/=0) return
!      end if
!      allocate(array(l),stat=stat)
!      if (stat/=0) return
!      do i = 1, l
!        array(i) = string(i:i)
!      end do
! 
!    END SUBROUTINE ZD11_put
! 
!    FUNCTION ZD11_get(array)
!      CHARACTER, intent(in):: array(:)
!      CHARACTER(size(array)) ::  ZD11_get
! ! Give the value of array to string.
! 
!      integer :: i
!      do i = 1, size(array)
!         ZD11_get(i:i) = array(i)
!      end do
! 
!    END FUNCTION ZD11_get
! 
! !-*-*-*-*-  G A L A H A D -  S T R I N G _ g e t   F U N C T I O N  -*-*-*-*-
! 
!      FUNCTION STRING_get( array )
! 
! !  obtain the elements of a character array as a character variable
! 
! !  Dummy arguments
! 
! !  array - character array whose components hold the string
! !  string_get - equivalent character string
! 
!      CHARACTER, INTENT( IN ), DIMENSION( : ) :: array
!      CHARACTER( SIZE( array ) ) :: STRING_get
! 
! !  Local variables
! 
!      INTEGER :: i
! 
!      DO i = 1, SIZE( array )
!         STRING_get( i : i ) = array( i )
!      END DO
! 
!      RETURN
! 
! !  End of function STRING_get
! 
!      END FUNCTION STRING_get
! 
! END MODULE HSL_ZD11D
! ! ===============================================================
! ===============================================================
! ===============================================================
! ---------------------------------------------------------
! Code below by Bob Fischer
! This stuff is for interfacing

!> @file
!! Code used to bind GALAHAD's data structure ZD11_type (sparse matrix
!! representation) to a C++ peer class.
!! @see galahad,galahd::qpt_problem_c
!! @see galahad,galahd::zd11_c

!> Container module for the Fortran peer class.
module zd11_x

use, intrinsic :: iso_c_binding

    !> Fortran side of the C++ peer class to zd11_type
    !! @see giss::ZD11
    type ZD11_c
        ! Pointer to the main struct (but we don't own this pointer)
        type(c_ptr) :: this_fc      ! zd11_type *

        ! Make portions of main available to C
        type(c_ptr) :: m, n, ne     ! int * (scalar)
        type(c_ptr) :: row          ! int[m]
        type(c_ptr) :: col          ! int[n]
        type(c_ptr) :: val          ! double[ne]
    end type

CONTAINS


!> Fill in references and pointers for the C++ peer class giss::ZD11
!! from an instance of the Fortran derived type zd11_type.
!! Also initializes the Fortran derived type with type of 'COORDINATE'.
!! To be called from Fortran.
!! @param self Instance of the C++ peer class giss::ZD11 (represented here by the Fortran class zd11_c).
!! @param main Instance of the Fortran derived type
!! @param m Number of rows this matrix should have
!! @param n Number of columns this matrix should have
!! @param ne Number of non-zero elements this matrix should have
!! @see hsl_zd11d::zd11_type
subroutine ZD11_init_f(this_c, this_f, m, n, ne)
use HSL_ZD11_double
use c_loc_x
use, intrinsic :: iso_c_binding
IMPLICIT NONE
!type(c_ptr), value :: this_cc
type(ZD11_c) :: this_c
type(zd11_type), target :: this_f
integer, value :: m, n, ne

integer :: stat

    this_f%m = m
    this_f%n = n
    this_f%ne = ne

    this_c%this_fc = c_loc(this_f)
    this_c%m = c_loc(this_f%m)
    this_c%n = c_loc(this_f%n)
    this_c%ne = c_loc(this_f%ne)
    allocate(this_f%row(ne))
    this_c%row = c_loc_array_int(this_f%row)
    allocate(this_f%col(ne))
    this_c%col = c_loc_array_int(this_f%col)
    allocate(this_f%val(ne))
    this_c%val = c_loc_array_double(this_f%val)

    call ZD11_put(this_f%type, 'COORDINATE', stat)
end subroutine ZD11_init_f


end module zd11_x

! ======================================================
! Functions to be called from C, so they're outside of a module





!subroutine ZD11_c_destroy(self)
!type(ZD11_c) :: self
!   type(zd11_type), pointer :: main
!   call c_f_pointer(self%main, main)
!   deallocate(main)
!subroutine ZD11_c_destroy(self)

!> Helper function for giss::ZD11::put_type()
function ZD11_put_type_c(this_fc, string, l) bind(c)
    use HSL_ZD11_double
    use, intrinsic :: iso_c_binding
implicit none
type(c_ptr), value :: this_fc               ! zd11_f *
    type(zd11_type), pointer :: this_f
character(kind=c_char), dimension(*) :: string  ! char *
integer(c_int) :: l                     ! strlen(str)
integer(c_int) :: ZD11_put_type_c


    call c_f_pointer(this_fc, this_f)
!   call c_f_pointer(string_c, string)

     if (allocated(this_f%type)) then
        deallocate(this_f%type,stat=ZD11_put_type_c)
        if (ZD11_put_type_c/=0) return
     end if
     allocate(this_f%type(l),stat=ZD11_put_type_c)
     if (ZD11_put_type_c/=0) return
     this_f%type(1:l) = string(1:l)
end function ZD11_put_type_c

