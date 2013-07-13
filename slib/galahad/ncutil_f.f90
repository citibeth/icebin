MODULE ncutil_mod

use netcdf

implicit none

!include 'netcdf.inc'

CONTAINS
! -------------------------------------------------------------
subroutine check(status)
implicit none
  integer, intent ( in) :: status

  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop 2
  end if
end subroutine check

! -------------------------------------------------------------
subroutine nc_read_dim(nc, dimname, dim)
integer, intent(in) :: nc
character*(*), intent(in) :: dimname
integer, intent(out) :: dim

integer :: dimid

	call check(nf90_inq_dimid(nc, dimname, dimid))
	call check(nf90_inquire_dimension(nc, dimid, len=dim))
end subroutine nc_read_dim
! -------------------------------------------------------------
subroutine nc_read_dims(nc, ndims, dimids, dims)
integer, intent(in) :: nc
integer, intent(in) :: ndims
integer, dimension(ndims), intent(in) :: dimids
integer, dimension(ndims) :: dims
integer :: i

	do i=1,ndims
		call check(nf90_inquire_dimension(nc, dimids(i), len=dims(i)))
	end do
end subroutine nc_read_dims

! -------------------------------------------------------------
subroutine nc_read_1d_array_double(nc, name, A)
implicit none
integer, intent(in) :: nc
character(len=*), intent(in) :: name
!character(len=len(name)) :: name2
real*8,dimension(:),allocatable, intent(out) :: A

integer, dimension(NF90_MAX_VAR_DIMS) :: dimids
integer :: ndims
integer, dimension(NF90_MAX_VAR_DIMS) :: dims

	integer :: varid

!	write (6,*) 'Reading variable',name
	call check(nf90_inq_varid(nc, name, varid))
	call check(nf90_inquire_variable(nc, varid, &
		ndims=ndims, dimids=dimids))

	if (ndims /= 1) then
		WRITE( 6, "('Expected 1 dimension, got ', I4)") ndims
    	stop 2
	end if

	call nc_read_dims(nc, ndims, dimids, dims)
	allocate(A(dims(1)))

	call check(nf90_get_var(nc, varid, A))

end subroutine nc_read_1d_array_double
! -------------------------------------------------------------
subroutine nc_read_1d_array_int(nc, name, A)
implicit none
integer, intent(in) :: nc
character(len=*), intent(in) :: name
integer,dimension(:),allocatable, intent(out) :: A

integer, dimension(NF90_MAX_VAR_DIMS) :: dimids
integer :: ndims
integer, dimension(NF90_MAX_VAR_DIMS) :: dims

	integer :: varid

!	write (6,*) 'Reading variable',name
	call check(nf90_inq_varid(nc, name, varid))
	call check(nf90_inquire_variable(nc, varid, &
		ndims=ndims, dimids=dimids))

	if (ndims /= 1) then
		WRITE( 6, "('Expected 1 dimension, got ', I4)") ndims
    	stop 2
	end if

	call nc_read_dims(nc, ndims, dimids, dims)
	allocate(A(dims(1)))

	call check(nf90_get_var(nc, varid, A))

end subroutine nc_read_1d_array_int
! -------------------------------------------------------------
subroutine nc_read_2d_array_int(nc, name, A)
implicit none
integer, intent(in) :: nc
character(len=*), intent(in) :: name
integer,dimension(:,:),allocatable, intent(out) :: A

integer, dimension(NF90_MAX_VAR_DIMS) :: dimids
integer :: ndims
integer, dimension(NF90_MAX_VAR_DIMS) :: dims

	integer :: varid

!	write (6,*) 'Reading variable',name
	call check(nf90_inq_varid(nc, name, varid))
	call check(nf90_inquire_variable(nc, varid, &
		ndims=ndims, dimids=dimids))

	if (ndims /= 2) then
		WRITE( 6, "('Expected 2 dimensions, got ', I4)") ndims
    	stop 2
	end if

	call nc_read_dims(nc, ndims, dimids, dims)
	allocate(A(dims(1), dims(2)))
	call check(nf90_get_var(nc, varid, A))

end subroutine nc_read_2d_array_int
! -------------------------------------------------------------
subroutine nc_read_attribute_str(nc, var_name, att_name, att_val)
integer, intent(in) :: nc
character(*), intent(in) :: var_name
character(*), intent(in) :: att_name
!character, dimension(:), allocatable, intent(out) :: att_val
character(*), intent(out) :: att_val

	integer :: varid, strlen

	call check(nf90_inq_varid(nc, trim(var_name), varid))
	call check(nf90_inquire_attribute(nc, varid, trim(att_name), len=strlen))
!	allocate(att_val(strlen))
	call check(nf90_get_att(nc, varid, trim(att_name), att_val))
end subroutine nc_read_attribute_str
! -------------------------------------------------------------
! subroutine nc_read_attribute_int(nc, var_name, att_name, att_val)
! integer, intent(in) :: nc
! character(*), intent(in) :: var_name
! character(*), intent(in) :: att_name
! !character, dimension(:), allocatable, intent(out) :: att_val
! integer, intent(out) :: att_val
! 
! 	integer :: varid, strlen
! 
! 	call check(nf90_inq_varid(nc, trim(var_name), varid))
! 	call check(nf90_get_att(nc, varid, trim(att_name), att_val))
! end subroutine nc_read_attribute
! -------------------------------------------------------------

! subroutine nc_read_sparse_matrix(nc, name, A)
! USE GALAHAD_QP_double
! implicit none
! 
! integer, intent(in) :: nc
! character, intent(in) :: name*(*)
! type(smt_type), intent(out) :: A
! 
! integer :: nrows, ncols
! integer :: dotpos
! integer :: err
! 
! character :: varname*100
! 	varname=name
! 	dotpos = len(name)+1
! 
! 	varname(dotpos:) = '.nrows'
! 	print *, 'AA1  ',trim(varname)
! 	call nc_read_dim(nc, varname, A%m)
! 
! 	varname(dotpos:) = '.ncols'
! 	print *, 'AA2  ',trim(varname)
! 	call nc_read_dim(nc, varname, A%n)
! 
! 	CALL SMT_put( A%type, 'SPARSE_BY_ROWS', err )
! 
! 	varname(dotpos:) = '.ptr'
! 	print *, 'AA1  ',trim(varname)
! 	call nc_read_1d_array_int(nc, varname, A%ptr)
! 	varname(dotpos:) = '.col'
! 	print *, 'AA1  ',trim(varname)
! 	call nc_read_1d_array_int(nc, varname, A%col)
! 	varname(dotpos:) = '.val'
! 	print *, 'AA1  ',trim(varname)
! 	call nc_read_1d_array_double(nc, varname, A%val)
! 	print *,'done'
! 
! 	A%ne = size(A%val,1)
! 
! end subroutine nc_read_sparse_matrix
! -------------------------------------------------------------

END MODULE ncutil_mod

