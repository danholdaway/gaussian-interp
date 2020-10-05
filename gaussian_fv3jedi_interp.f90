program test_interp

! standard libaries
use mpi
use netcdf

! fckit
use fckit_mpi_module, only: fckit_mpi_comm

! oops
use unstructured_interpolation_mod, only: unstrc_interp

! fv3jedi
use fv3jedi_constants_mod,   only: deg2rad
use fv3jedi_kinds_mod,       only: kind_real
use fv3jedi_bump_interp_mod, only: fv3jedi_bump_interp

implicit none

integer :: mpierr, rank, size
integer :: ncid, dimid, varid
integer :: window_length = 0, time_step = 3
integer :: n, number_files
integer :: grid_xt, grid_yt, grid_zt

real(kind=kind_real), allocatable :: field(:,:,:)
real(kind=kind_real), allocatable :: lon(:,:), lat(:,:)
real(kind=kind_real) :: minloc, maxloc, minglo, maxglo

character(len=2055), allocatable :: filename(:)
character(len=2055) :: field_name = 'tmp'

integer :: i, j, k
integer :: npes_x = 1, npes_y = 1
integer :: grid_xt_each, grid_yt_each
integer :: isc, iec, jsc, jec
integer :: starts(4), counts(4)

integer :: nobs
real(kind=kind_real), allocatable :: lonobs(:), latobs(:), fieldobs(:,:)

type(fckit_mpi_comm) :: comm
type(fv3jedi_bump_interp) :: bump_interp

real(kind=kind_real), allocatable :: lats_in(:), lons_in(:), field_us(:)
integer :: ngrid_in, jnode
type(unstrc_interp) :: unst_interp

! MPI
! ---
call mpi_init(mpierr)
call mpi_comm_size(mpi_comm_world, size, mpierr)
call mpi_comm_rank(mpi_comm_world, rank, mpierr)

if (size .ne. npes_x*npes_y) then
  print*, "Wrong number of tasks"
  call mpi_abort(mpi_comm_world, 1, mpierr)
endif

comm = fckit_mpi_comm("world")

! Number of files needed
! ----------------------
number_files = window_length/time_step + 1
if (rank == 0) print*, 'Number of files: ', number_files

! Set the filenames
! -----------------
allocate(filename(number_files))
filename(1) = '/gpfsm/dnb31/drholdaw/JediWork/GSIValidation/HorzInterp/Backgrounds/gdas.t18z.atmf006.nc'

!filename(2) = '/gpfsm/dnb31/drholdaw/HeraInbox/gdas.t12z.atmf001.nc'
!filename(3) = '/gpfsm/dnb31/drholdaw/HeraInbox/gdas.t12z.atmf002.nc'

!filename(1) = '/gpfsm/dnb31/drholdaw/HeraInbox/gdas.t12z.atmf006_c192.nc'

! Loop over time levels and read the data
! ---------------------------------------
do n = 1, number_files

  ! Print file name
  if (rank == 0) print*, ' Reading the file: ', trim(filename(n))

  ! Open the file
  call nccheck ( nf90_open( filename(1), NF90_NOWRITE, ncid ), "nf90_open" )

  ! Setup done after opening first file
  if (n == 1) then

    ! Get dimensions of the fields and allocate field
    call nccheck ( nf90_inq_dimid(ncid, 'grid_xt', dimid), "nf90_inq_dimid x")
    call nccheck ( nf90_inquire_dimension(ncid, dimid, len = grid_xt), "nf90_inquire_dimension x")

    call nccheck ( nf90_inq_dimid(ncid, 'grid_yt', dimid), "nf90_inq_dimid y")
    call nccheck ( nf90_inquire_dimension(ncid, dimid, len = grid_yt), "nf90_inquire_dimension y")

    call nccheck ( nf90_inq_dimid(ncid, 'pfull', dimid), "nf90_inq_dimid z")
    call nccheck ( nf90_inquire_dimension(ncid, dimid, len = grid_zt), "nf90_inquire_dimension z")

    ! Index ranges
    if (size > 1) then
      ! Domain decomposition
      if (mod(grid_xt,npes_x) .ne. 0 .or. mod(grid_yt,npes_y) .ne. 0) then
        print*, "mod(grid_xt,npes_x) .ne. 0 .or. mod(grid_yt,npes_y) .ne. 0"
        call mpi_abort(mpi_comm_world, 1, mpierr)
      endif
      grid_xt_each = floor(real(grid_xt,kind_real)/real(npes_x,kind_real))
      grid_yt_each = floor(real(grid_yt,kind_real)/real(npes_y,kind_real))

      isc = mod(rank,npes_x)*grid_xt_each + 1
      iec = mod(rank,npes_x)*grid_xt_each + grid_xt_each
      jsc = mod(rank/npes_x,npes_y)*grid_yt_each + 1
      jec = mod(rank/npes_x,npes_y)*grid_yt_each + grid_yt_each
    else
      ! Read custom range with single PE
      isc = 480
      iec = 550
      jsc = 540
      jec = 590
    endif

    ! Starts and counts for reading
    starts(1) = isc
    starts(2) = jsc
    starts(3) = 1
    starts(4) = 1
    counts(1) = iec-isc+1
    counts(2) = jec-jsc+1
    counts(3) = grid_zt
    counts(4) = 1

    ! Allocate field
    allocate(lon(isc:iec,jsc:jec))
    allocate(lat(isc:iec,jsc:jec))
    allocate(field(isc:iec, jsc:jec, grid_zt))

    ! Get lon/lat
    if (rank == 0) print*, 'Reading lon/lat'
    call nccheck ( nf90_inq_varid (ncid, 'lon', varid), "nf90_inq_varid lon")
    call nccheck ( nf90_get_var( ncid, varid, lon(isc:iec,jsc:jec), starts(1:2), counts(1:2)), "nf90_get_var lon")
    call nccheck ( nf90_inq_varid (ncid, 'lat', varid), "nf90_inq_varid lat")
    call nccheck ( nf90_get_var( ncid, varid, lat(isc:iec,jsc:jec), starts(1:2), counts(1:2)), "nf90_get_var lat")

  endif

  ! Read the variable to be interpolated
  ! ------------------------------------
  if (rank == 0) print*, 'Reading ', trim(field_name)
  call nccheck ( nf90_inq_varid (ncid, trim(field_name), varid), "nf90_inq_varid")
  call nccheck ( nf90_get_var( ncid, varid, field(isc:iec,jsc:jec,:), starts, counts), "nf90_get_var")

  ! Print field information
  ! -----------------------
  minloc = minval(field(isc:iec,jsc:jec,:))
  maxloc = maxval(field(isc:iec,jsc:jec,:))
  call mpi_reduce (minloc, minglo, 1, mpi_double, mpi_min, 0, mpi_comm_world, mpierr)
  call mpi_reduce (maxloc, maxglo, 1, mpi_double, mpi_max, 0, mpi_comm_world, mpierr)

  if (rank == 0) print*, 'Name, Max, Min: ', trim(field_name), minglo, maxglo

  minloc = minval(lon(isc:iec,jsc:jec))
  maxloc = maxval(lon(isc:iec,jsc:jec))
  call mpi_reduce (minloc, minglo, 1, mpi_double, mpi_min, 0, mpi_comm_world, mpierr)
  call mpi_reduce (maxloc, maxglo, 1, mpi_double, mpi_max, 0, mpi_comm_world, mpierr)

  if (rank == 0) print*, 'Name, Max, Min: ', 'lon', minloc, maxloc

  minloc = minval(lat(isc:iec,jsc:jec))
  maxloc = maxval(lat(isc:iec,jsc:jec))
  call mpi_reduce (minloc, minglo, 1, mpi_double, mpi_min, 0, mpi_comm_world, mpierr)
  call mpi_reduce (maxloc, maxglo, 1, mpi_double, mpi_max, 0, mpi_comm_world, mpierr)

  if (rank == 0) print*, 'Name, Max, Min: ', 'lat', minglo, maxglo

  ! Close the file
  call nccheck ( nf90_close(ncid), "nf90_close")

enddo


! Observations
! ------------
nobs = 1
if (rank == 0) then
  nobs = 1
endif

allocate(lonobs(nobs))
allocate(latobs(nobs))
allocate(fieldobs(nobs,grid_zt))

if (rank == 0) then
  lonobs(1) = 62.2_kind_real
  latobs(1) = 23.22_kind_real
else
  lonobs(1) = 0.0_kind_real + 0.05*real(rank,kind_real)
  latobs(1) = 0.0_kind_real
endif

fieldobs = 0.0_kind_real


! Interpolation
! -------------
!if (rank == 0) then
  do j = jsc, jec
    do i = isc, iec
      if (lon(i,j) < 0.0_kind_real) lon(i,j) = lon(i,j) + 360.0_kind_real
      !print*, i, j, lon(i,j), lat(i,j)
    enddo
  enddo
!endif


! Interpolation BUMP
! ------------------
call bump_interp%setup (comm, isc, iec, jsc, jec, grid_zt, deg2rad*lon, deg2rad*lat, nobs, lonobs, latobs)
call bump_interp%apply (grid_zt, field(isc:iec,jsc:jec,:), nobs, fieldobs)
call bump_interp%delete()


! Print output
! ------------
if (rank == 0) then
  do k = 1, grid_zt
    print*, 'atmp_jedi[', grid_zt-k,'] = ', fieldobs(1,k)
  enddo
endif

! Finalize
! --------
call mpi_finalize(mpierr)

contains

! --------------------------------------------------------------------------------------------------

subroutine nccheck(status,iam)

  implicit none
  integer, intent ( in) :: status
  character(len=*), optional :: iam

  character(len=1024) :: error_descr

  if(status /= nf90_noerr) then

    error_descr = "NetCDF error, aborting ... "

    if (present(iam)) then
      error_descr = trim(error_descr)//", "//trim(iam)
    endif

    error_descr = trim(error_descr)//". Error code: "//trim(nf90_strerror(status))

    print*, trim(error_descr)

    call mpi_abort(mpi_comm_world, 1, mpierr)

  end if

end subroutine nccheck

! --------------------------------------------------------------------------------------------------

end program
