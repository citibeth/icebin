! COPYRIGHT (c) 2006 Council for the Central Laboratory
!                    of the Research Councils
! Original date 21 February 2006. Version 1.0.0.
! 6 March 2007 Version 1.1.0. Argument stat made non-optional

! ===============================================================
! ===============================================================
! ===============================================================
!> Sparse Matrix Representation Types
MODULE HSL_ZD11D

!  ==========================
!  Sparse matrix derived type
!  ==========================

  !> Sparse matrix representation
  !! Originally derived from HSL library (via GALAHAD; same as SMT_TYPE in GALAHAD).
  TYPE, PUBLIC :: zd11_type
    INTEGER :: m	!< Number of rows in this matrix.
	INTEGER :: n	!< Number of columns in this matrix.
	INTEGER :: ne	!< Number of non-zero elements in this matrix.

	!> "Name" of this matrix, not really needed.
    CHARACTER, ALLOCATABLE, DIMENSION(:) :: id

	!> Indicates the storage scheme used.
	!! It is set equal to either
	!! 'DENSE', 'COORDINATE', 'SPARSE BY ROWS' or 'DIAGONAL'.  Use ZD11_put()
	!! to set this variable.
    CHARACTER, ALLOCATABLE, DIMENSION(:) :: type

	!> If type == 'COORDINATE', then this holds the row indices of
	!! the non-zero elements (index base = 1).  Otherwise, it is not needed.
    INTEGER, ALLOCATABLE, DIMENSION(:) :: row

	!> If type == 'COORDINATE' or 'SPARSE BY ROWS', this holds the column indices
	!! of the non-zero elements (index base = 1).  Otherwise, it is not needed.
    INTEGER, ALLOCATABLE, DIMENSION(:) :: col

	!> If type == 'SPARSE BY ROWS', this holds the starting position
	!! of each row in the arrays col and val (index base = 1)
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ptr

	!> The value of each non-zero element in the matrix.
    REAL ( KIND( 1.0D+0 ) ), ALLOCATABLE, DIMENSION(:) :: val
  END TYPE

CONTAINS

   SUBROUTINE ZD11_put(array,string,stat)
     CHARACTER, allocatable :: array(:)
     CHARACTER(*), intent(in) ::  string
     INTEGER, intent(OUT) ::  stat

     INTEGER :: i,l

     l = len_trim(string)
     if (allocated(array)) then
        deallocate(array,stat=stat)
        if (stat/=0) return
     end if
     allocate(array(l),stat=stat)
     if (stat/=0) return
     do i = 1, l
       array(i) = string(i:i)
     end do

   END SUBROUTINE ZD11_put

   FUNCTION ZD11_get(array)
     CHARACTER, intent(in):: array(:)
     CHARACTER(size(array)) ::  ZD11_get
! Give the value of array to string.

     integer :: i
     do i = 1, size(array)
        ZD11_get(i:i) = array(i)
     end do

   END FUNCTION ZD11_get

!-*-*-*-*-  G A L A H A D -  S T R I N G _ g e t   F U N C T I O N  -*-*-*-*-

     FUNCTION STRING_get( array )

!  obtain the elements of a character array as a character variable

!  Dummy arguments

!  array - character array whose components hold the string
!  string_get - equivalent character string

     CHARACTER, INTENT( IN ), DIMENSION( : ) :: array
     CHARACTER( SIZE( array ) ) :: STRING_get

!  Local variables

     INTEGER :: i

     DO i = 1, SIZE( array )
        STRING_get( i : i ) = array( i )
     END DO

     RETURN

!  End of function STRING_get

     END FUNCTION STRING_get

END MODULE HSL_ZD11D
! ===============================================================
! ===============================================================
! ===============================================================
! ---------------------------------------------------------
! Code below by Bob Fischer
! This stuff is for interfacing

module zd11_x

use, intrinsic :: iso_c_binding

	!> Fortran side of the C++ peer class to zd11_type
	!! @see giss::ZD11
	type ZD11_c
		! Pointer to the main struct (but we don't own this pointer)
		type(c_ptr) :: main			! zd11_type *

		! Make portions of main available to C
	    type(c_ptr) :: m, n, ne		! int * (scalar)
		type(c_ptr) :: row			! int[m]
		type(c_ptr) :: col			! int[n]
		type(c_ptr) :: val			! double[ne]
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
subroutine ZD11_c_init(self, main, m, n, ne)
use HSL_ZD11D
use c_loc_x
use, intrinsic :: iso_c_binding
IMPLICIT NONE
type(ZD11_c) :: self
type(zd11_type), target :: main
integer, value :: m, n, ne

integer :: stat

	main%m = m
	main%n = n
	main%ne = ne

	self%main = c_loc(main)
	self%m = c_loc(main%m)
	self%n = c_loc(main%n)
	self%ne = c_loc(main%ne)
	allocate(main%row(ne))
	self%row = c_loc_array_int(main%row)
	allocate(main%col(ne))
	self%col = c_loc_array_int(main%col)
	allocate(main%val(ne))
	self%val = c_loc_array_double(main%val)

	call ZD11_put(main%type, 'COORDINATE', stat)
end subroutine ZD11_c_init


end module zd11_x

! ======================================================
! Functions to be called from C, so they're outside of a module





!subroutine ZD11_c_destroy(self)
!type(ZD11_c) :: self
!	type(zd11_type), pointer :: main
!	call c_f_pointer(self%main, main)
!	deallocate(main)
!subroutine ZD11_c_destroy(self)

!> Helper function for giss::ZD11::put_type()
function ZD11_put_type_c(self_c, string, l) bind(c)
	use HSL_ZD11_double
	use, intrinsic :: iso_c_binding
implicit none
type(c_ptr) :: self_c				! zd11_f *
character, dimension(*) :: string	! char *
integer :: l						! strlen(str)
integer :: ZD11_put_type_c

	type(zd11_type), pointer :: self

	call c_f_pointer(self_c, self)


     if (allocated(self%type)) then
        deallocate(self%type,stat=ZD11_put_type_c)
        if (ZD11_put_type_c/=0) return
     end if
     allocate(self%type(l),stat=ZD11_put_type_c)
     if (ZD11_put_type_c/=0) return
     self%type(1:l) = string(1:l)
end function ZD11_put_type_c

