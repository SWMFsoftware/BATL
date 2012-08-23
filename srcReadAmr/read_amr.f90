program read_amr

  ! reconstruct AMR grid and populate data onto that grid
  use BATL_lib, ONLY: iProc, init_mpi, clean_mpi, coord_to_xyz
  use ModReadAmr, ONLY: nVar, CoordMin_D, CoordMax_D, &
       readamr_init, readamr_read, readamr_get, readamr_clean

  implicit none
  
  ! File name
  character(len=100) :: NameFile

  real :: Xyz_D(3), Coord_D(3), Weight
  real, allocatable:: State_V(:)
  logical:: IsFound
  !----------------------------------------------------------------------------
  ! Initialize MPI (for parallel execution)
  call init_mpi


  ! NameFile='input/3d__all_2_n0000010'
  NameFile='input/3d__all_1_t25.60000_n0000386'

  if(iProc==0) write(*,*) 'Initializing READAMR'
  call readamr_init(NameFile, IsVerboseIn = iProc==0)

  if(iProc==0) write(*,*) 'Reading data'
  call readamr_read(trim(NameFile)//'.out', UseXyzTest=.true.)

  ! Allocate state array
  allocate(State_V(nVar))
  if(iProc==0)write(*,*)'Number of stored variables=', nVar

  ! test obtaining a point value
  Coord_D = (CoordMin_D + CoordMax_D)/2
  call coord_to_xyz(Coord_D, Xyz_D)
  if(iProc==0)write(*,*)'Testing at Coord_D, Xyz_D=', Coord_D, Xyz_D
  call readamr_get(Xyz_D, nVar, State_V, Weight, IsFound)
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
