!^CFG COPYRIGHT UM

module BATL_amr

  use BATL_size, ONLY: MaxBlock, &
       nI, nJ, nK, nIjk_D, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
       MaxDim, nDim, nDimAmr

  use BATL_mpi, ONLY: iComm, nProc, iProc, barrier_mpi

  use BATL_tree, ONLY: nNodeUsed, iNodePeano_I, &
       iTree_IA, Proc_, Block_, ProcNew_, BlockNew_, &
       Coord1_, Coord2_, Coord3_, Status_, Refine_, Coarsen_

  use ModMpi

  implicit none

  SAVE

  private ! except

  public do_amr

  integer :: nWidth

contains

  !===========================================================================
  subroutine do_amr(nVar, State_VGB)

    ! Arguments
    integer, intent(in) :: nVar
    real, intent(inout) :: &
         State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

    integer, parameter:: Send_=1, Recv_=2
    integer :: iSendRecv, iMorton, iStatus
    integer :: iNodeSend, iNodeRecv, iProcSend, iProcRecv
    !-------------------------------------------------------------------------

    do iSendRecv = Send_, Recv_
       do iMorton = 1, nNodeUsed
          iNodeSend = iNodePeano_I(iMorton)

          iStatus = iTree_IA(Status_,iNodeSend)
          if(iStatus /= Coarsen_ .and. iStatus /= Refine_) CYCLE

          iProcSend = iTree_IA(Proc_,iNodeSend)
          iProcRecv = iTree_IA(ProcNew_,iNodeSend)

          if(iSendRecv == Send_ .and. iProc == iProcSend)then
             ! put restricted and prolonged state into buffer
             if(iStatus == Coarsen_)then
                call coarsen_block
             else
                call refine_block
             end if
          elseif(iSendRecv == Recv_ .and. iProc == iProcRecv)then
             ! read buffers into state
             call buffer_to_state
          end if

       end do
       if(iSendRecv == Recv_) EXIT
       ! message pass buffers
    end do

  contains
    !=========================================================================
    subroutine coarsen_block
    end subroutine coarsen_block
    !=========================================================================
    subroutine refine_block
    end subroutine refine_block
    !=========================================================================
    subroutine buffer_to_state
    end subroutine buffer_to_state

  end subroutine do_amr

end module BATL_amr
