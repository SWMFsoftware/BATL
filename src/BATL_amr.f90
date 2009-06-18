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

    do iMorton = 1, nNodeUsed
       iNodeSend = iNodePeano_I(iMorton)

       iStatus = iTree_IA(Status_,iNodeSend)
       if(iStatus /= Coarsen_ .and. iStatus /= Refine_) CYCLE

       ! put restricted and prolonged state into buffer
       if(iStatus == Coarsen_)then
          call coarsen_block_to_buffer
       else
          call refine_block_to_buffer
       end if
    end do

    ! message pass buffers

    do iMorton = 1, nNodeUsed
       iNodeSend = iNodePeano_I(iMorton)

       iStatus = iTree_IA(Status_,iNodeSend)
       if(iStatus /= Coarsen_ .and. iStatus /= Refine_) CYCLE

       ! put restricted and prolonged state into buffer
       if(iStatus == Coarsen_)then
          call buffer_to_coarse_block
       else
          call buffer_to_fine_block
       end if
    end do

  contains

    !=========================================================================
    subroutine coarsen_block_to_buffer

      integer :: i, j, k, iVar
      !-----------------------------------------------------------------------
      iProcRecv = iTree_IA(ProcNew_,iNodeSend)

      iBufferS = iBufferS_P(iProcRecv)
      do k = 1, nK, kRatio; do j = 1, nJ, jRatio; do i=1, nI, iRatio
         do iVar = 1, nVar
            BufferS_IP(iBufferS+iVar,iProcRecv) = &
                 InvIjkRatio * &
                 sum(State_VGB(iVar,i:i+iRatio-1,j:j+jRatio-1,k:k+kRatio-1,&
                 iBlockSend))
         end do
         iBufferS = iBufferS + nVar
      end do; end do; end do

      iBufferR_P(iProcSend) = iBufferR_P(iProcSend) &
           + iBufferS - iBufferS_P(iProcRecv)
      iBufferS_P(iProcRecv) = iBufferS
      
    end subroutine coarsen_block_to_buffer
    !=========================================================================
    subroutine refine_block

      iNodeParent = iTree_IA(Parent_,iNodeSend)
      iNodeChild1 = iTree_IA(Child1_,iNodeSend)
      iProcRecv   = iTree_IA(Proc_,iNodeChild1)


    end subroutine refine_block
    !=========================================================================
    subroutine buffer_to_state
    end subroutine buffer_to_state

  end subroutine do_amr

end module BATL_amr
