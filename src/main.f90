program BATL_test

  use BATL_tree, ONLY: test_tree
  use BATL_grid, ONLY: test_grid

  implicit none

  call test_tree
  call test_grid
  
end program BATL_test
!=============================================================================
subroutine CON_stop(String)
  character (len=*), intent(in) :: String
  write(*,*)'CON_stop called with String='
  write(*,*) String
  stop
end subroutine CON_stop
!=============================================================================
subroutine stop_mpi(String)
  character (len=*), intent(in) :: String
  write(*,*)'stop_mpi called with String='
  write(*,*) String
  stop
end subroutine stop_mpi
