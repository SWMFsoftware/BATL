module BATL_mpi

  use ModMpi

  implicit none

  SAVE

  integer:: nProc, iProc, iComm

contains

  !==========================================================================
  subroutine init_mpi(iCommIn)

    integer, intent(in):: iCommIn
    integer :: iError
    !------------------------------------------------------------------------
    iComm = iCommIn
    call MPI_COMM_RANK (iComm, iProc, iError)
    call MPI_COMM_SIZE (iComm, nProc, iError)

  end subroutine init_mpi

end module BATL_mpi
