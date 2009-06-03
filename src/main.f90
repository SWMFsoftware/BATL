program BATL_test

  use BATL_mpi,  ONLY: init_mpi
  use BATL_tree, ONLY: test_tree
  use BATL_grid, ONLY: test_grid
  use ModMpi

  implicit none

  integer:: iError
  !--------------------------------------------------------------------------
  call MPI_init(iError)
  call init_mpi(MPI_COMM_WORLD)
  call test_tree
  call test_grid
  call MPI_finalize(iError)
  
end program BATL_test
!=============================================================================
subroutine CON_stop(String)
  use BATL_mpi
  implicit none
  character (len=*), intent(in) :: String
  integer:: iError, nError
  !--------------------------------------------------------------------------
  write(*,*)'CON_stop called with String='
  write(*,*) String
  call MPI_abort(iComm, nError, iError)
  stop
end subroutine CON_stop
