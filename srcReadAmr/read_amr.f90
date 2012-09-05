program read_amr_test

  ! reconstruct AMR grid and populate data onto that grid
  use BATL_lib, ONLY: iProc, nProc, iComm, nI, nJ, nK, nDim, &
       init_mpi, clean_mpi, coord_to_xyz
  use ModReadAmr, ONLY: nVar, CoordMin_D, CoordMax_D, &
       readamr_read, readamr_get, readamr_clean
  use ModConst, ONLY: cPi
  use ModMpi, ONLY: MPI_REAL, MPI_SUM, MPI_allreduce

  implicit none
  
  ! File name
  character(len=100) :: NameFile

  integer, parameter:: n = 10 ! number of test points in all directions

  real :: Xyz_D(3), Coord_D(3), State_D(nDim), Cos2_D(nDim), Tolerance = 1e-6
  real, allocatable:: State_V(:), StateLocal_V(:)
  real, allocatable:: State_VIII(:,:,:,:), StateLocal_VIII(:,:,:,:)
  integer:: i, j, k, iError
  logical:: IsFound
  character(len=*), parameter:: NameCode='READAMRTEST'
  !----------------------------------------------------------------------------
  ! Initialize MPI (for parallel execution)
  call init_mpi

  if(iProc==0)write(*,*) NameCode,' is starting with nI, nJ, nK=',nI, nJ, nK

  ! The tests are distinguished based on the configuration
  if(nI==4.and.nJ==4.and.nK==1)then
     NameFIle = "data/z=0_mhd_1_t00000010_n0000042.out"
     Tolerance=0.05
  elseif(nI==4.and.nJ==4.and.nK==4)then
     NameFIle = "data/3d__all_3_t00000010_n0000059.out"
     Tolerance=0.05
  elseif(nI==6.and.nJ==4.and.nK==4)then
     NameFIle = "data/3d__var_4_t00000000_n0000010.out"
     Tolerance=0.05
  else
     call CON_stop(' there is no test input for this block size')
  end if

  if(iProc==0) write(*,*) NameCode, &
       ' get variables 1,3 from ', trim(NameFile)
  call readamr_read(NameFile, iVarIn_I=(/1,3/), IsVerboseIn = iProc==0)

  ! Select a point inside the domain
  Coord_D = CoordMin_D + 0.2*(CoordMax_D - CoordMin_D)
  call coord_to_xyz(Coord_D, Xyz_D)

  ! Allocate state array. Zero index is needed for interpolation weight!
  allocate(State_V(0:nVar))

  ! Get data at this point
  call readamr_get(Xyz_D, State_V, IsFound)

  if(nProc > 0)then
     ! When running in parallel the contributions need to be collected
     allocate(StateLocal_V(0:nVar))
     StateLocal_V = State_V
     call MPI_allreduce(StateLocal_V, State_V, size(State_V), &
          MPI_REAL, MPI_SUM, iComm, iError)
     deallocate(StateLocal_V)
  end if
  
  if(iProc==0)then
     write(*,*) NameCode,' at Xyz_D=', Xyz_D
     write(*,*) NameCode,'  State_V=', State_V
  end if

  deallocate(State_V)

  ! Verification test

  if(iProc==0) write(*,*) NameCode,' reread data and do verification test'
  call readamr_read(NameFile, IsVerboseIn = iProc==0, &
       IsNewGridIn = .false., UseCoordTest=.true.)

  ! Allocate state array on an (n+1)^3 grid (uniform in generalized coords)
  allocate(State_VIII(0:nVar,0:n,0:n,0:n))
  State_VIII = 0.0
  if(iProc==0)write(*,*) NameCode,' number of stored variables=', nVar

  ! obtain state values on the test grid
  do k = 0, n; do j=0, n; do i = 0, n
     ! Make a uniform grid in the generalized coordinate space
     Coord_D = CoordMin_D + (CoordMax_D - CoordMin_D) &
          *(1e-6 + (/i,j,k/))/(2e-6 + n)
     call coord_to_xyz(Coord_D, Xyz_D)

     ! Get the point value
     call readamr_get(Xyz_D, State_VIII(0:nVar,i,j,k), IsFound)

     if(.not.IsFound)then
        write(*,*) NameCode,' ERROR: could not find point at i,j,k=',i,j,k
        write(*,*) NameCode,' CoordMin_D  = ', CoordMin_D
        write(*,*) NameCode,' CoordMax_D  = ', CoordMax_D
        write(*,*) NameCode,' Coord_D     = ', Coord_D   
        write(*,*) NameCode,' Xyz_D       = ', Xyz_D
     end if
  end do; end do; end do

  if(nProc > 0)then
     ! When running in parallel the contributions need to be collected
     allocate(StateLocal_VIII(0:nVar,0:n,0:n,0:n))
     StateLocal_VIII = State_VIII
     call MPI_allreduce(StateLocal_VIII, State_VIII, size(State_VIII), &
          MPI_REAL, MPI_SUM, iComm, iError)
     deallocate(StateLocal_VIII)
  end if

  if(iProc==0)write(*,*) NameCode,' check state on (n+1)^3 grid. n=', n
  
  do k = 0, n; do j=0, n; do i = 0, n
     ! Note the division by the total weight. Usually it is 1.0 but
     ! it can be different at the outer boundaries
     ! or at the edges and corners of resolution changes.
     State_D = State_VIII(1:nDim,i,j,k)/State_VIII(0,i,j,k)

     ! The stored function is cos^2 of normalized coordinate
     Coord_D = (1e-6 + (/i,j,k/))/(2e-6 + n)
     Cos2_D = cos(cPi*Coord_D(1:nDim))**2

     if(any(abs(Cos2_D - State_D) > Tolerance))then
        write(*,*) NameCode,' Test failed for i,j,k,n=',i,j,k,n
        write(*,*) NameCode,' CoordMin_D  = ', CoordMin_D
        write(*,*) NameCode,' CoordMax_D  = ', CoordMax_D
        write(*,*) NameCode,' Coord_D     = ', Coord_D
        write(*,*) NameCode,' Cos2_D      = ', Cos2_D
        write(*,*) NameCode,' State_D     = ', State_D
        write(*,*) NameCode,' Weight      = ', State_VIII(0,i,j,k)
     end if
  end do; end do; end do

  deallocate(State_VIII)

  ! Clean READAMR storage
  call readamr_clean

  if(iProc==0) write(*,*) NameCode,' finished'

  ! Finish MPI (for parallel execution)
  call clean_mpi

end program read_amr_test
!=============================================================================
subroutine CON_stop(StringError)
  ! This subroutine has to be provided to stop cleanly

  use ModMPI
  use BATL_lib, ONLY: iProc, iComm

  implicit none
  integer:: nError, iError
  character (len=*), intent(in) :: StringError
  !--------------------------------------------------------------------------
  write(*,*)'ERROR in READAMRTEST on processor ', iProc
  write(*,*) StringError
  call MPI_abort(iComm, nError, iError)
  stop

end subroutine CON_stop
