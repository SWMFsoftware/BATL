program read_amr

  ! reconstruct AMR grid and populate data onto that grid
  use BATL_lib, ONLY: iProc, nProc, iComm, init_mpi, clean_mpi, coord_to_xyz
  use ModReadAmr, ONLY: nVar, CoordMin_D, CoordMax_D, &
       readamr_read, readamr_get, readamr_clean
  use ModMpi, ONLY: MPI_CHARACTER, MPI_BCAST

  implicit none
  
  ! File name
  character(len=100) :: NameFile

  real :: Xyz_D(3), Coord_D(3), Weight
  real, allocatable:: State_V(:)
  integer:: iError
  logical:: IsFound
  !----------------------------------------------------------------------------
  ! Initialize MPI (for parallel execution)
  call init_mpi

  if(iProc==0)then
     write(*,*)'Provide data file name:'
     read(*,'(a)') NameFile
  end if
  if(nProc > 0) call MPI_bcast(NameFile, len(NameFile), MPI_CHARACTER, 0, &
       iComm, iError)

  if(iProc==0) write(*,*) 'Reading data'
  call readamr_read(NameFile, IsVerboseIn = iProc==0, UseXyzTest=.true.)

  ! Allocate state array
  allocate(State_V(nVar))
  if(iProc==0)write(*,*)'Number of stored variables=', nVar

  ! test obtaining a point value
  Coord_D = (CoordMin_D + CoordMax_D)/2
  call coord_to_xyz(Coord_D, Xyz_D)
  if(iProc==0)write(*,*)'Testing at Coord_D, Xyz_D=', Coord_D, Xyz_D
  call readamr_get(Xyz_D, State_V, Weight, IsFound)
  if(Weight>0.0) write(*,*)'iProc,Weight,State_V=', iProc, Weight, State_V

  deallocate(State_V)

  ! Clean READAMR storage
  call readamr_clean

  ! Finish MPI (for parallel execution)
  call clean_mpi

end program read_amr
!=============================================================================
subroutine CON_stop(StringError)
  ! This subroutine has to be provided to stop cleanly

  use ModMPI
  use BATL_lib, ONLY: iProc, iComm

  implicit none
  integer:: nError, iError
  character (len=*), intent(in) :: StringError
  !--------------------------------------------------------------------------
  write(*,*)'ERROR in READAMR on processor ', iProc
  write(*,*) StringError
  call MPI_abort(iComm, nError, iError)
  stop

end subroutine CON_stop
