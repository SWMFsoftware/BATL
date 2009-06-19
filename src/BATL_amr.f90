!^CFG COPYRIGHT UM

module BATL_amr

  use BATL_size, ONLY: MaxBlock, &
       nI, nJ, nK, nIjk_D, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
       MaxDim, nDim, nDimAmr, iRatio, jRatio, kRatio, InvIjkRatio

  use BATL_mpi, ONLY: iComm, nProc, iProc, barrier_mpi

  use BATL_tree, ONLY: nNodeUsed, iNodePeano_I, &
       iTree_IA, Proc_, Block_, ProcNew_, BlockNew_, &
       Coord1_, Coord2_, Coord3_, Status_, Refine_, RefineNew_, Coarsen_

  use ModMpi

  implicit none

  SAVE

  private ! except

  public do_amr

  integer :: iProcSend, iProcRecv
  integer :: MaxBuffer
  integer, allocatable:: iBufferR_P(:), iBufferS_P(:)
  real, allocatable:: BufferR_IP(:,:), BufferS_IP(:,:)

  integer:: iRequestR, iRequestS, iError
  integer, allocatable:: iRequestR_I(:), iRequestS_I(:), iStatus_II(:,:)

contains

  !===========================================================================
  subroutine do_amr(nVar, State_VGB)

    ! Arguments
    integer, intent(in) :: nVar
    real, intent(inout) :: &
         State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

    integer, parameter:: Send_=1, Recv_=2
    integer :: iSendRecv, iMorton, iStatus
    integer :: iNodeSend, iNodeRecv, iProcSend, iProcRecv, iBlockSend
    !-------------------------------------------------------------------------

    ! Small arrays are allocated once 
    if(.not.allocated(iBufferR_P))then
       allocate(iBufferR_P(0:nProc-1), iBufferS_P(0:nProc-1))
       allocate(iRequestR_I(nProc-1), iRequestS_I(nProc-1))
       allocate(iStatus_II(MPI_STATUS_SIZE,nProc-1))
    end if

    ! Big arrays are allocated and deallocated every time (for now)
    MaxBuffer = nVar*MaxBlock*nI*nJ*nK
    call timing_start('msg_allocate')
    allocate(BufferR_IP(MaxBuffer,0:nProc-1))
    allocate(BufferS_IP(MaxBuffer,0:nProc-1))
    call timing_stop('msg_allocate')

    do iMorton = 1, nNodeUsed
       iNodeSend = iNodePeano_I(iMorton)

       iStatus = iTree_IA(Status_,iNodeSend)

       ! put restricted and prolonged state into buffer
       if(iStatus /= Coarsen_ .and. iStatus /= RefineNew_) CYCLE

       iProcSend = iTree_IA(Proc_,iNodeSend)
       iProcRecv = iTree_IA(ProcNew_,iNodeSend)
       if(iProc /= iProcSend .and. iProc /= iProcRecv) CYCLE

       iBlockSend = iTree_IA(Block_,iNodeSend)

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
       if(iStatus /= Coarsen_ .and. iStatus /= RefineNew_) CYCLE

       iProcSend = iTree_IA(Proc_,iNodeSend)
       iProcRecv = iTree_IA(ProcNew_,iNodeSend)
       if(iProc /= iProcSend .and. iProc /= iProcRecv) CYCLE

       ! put restricted and prolonged state into buffer
       if(iStatus == Coarsen_)then
          call buffer_to_coarse_block
       else
          call buffer_to_fine_block
       end if
    end do

    deallocate(BufferR_IP, BufferS_IP)

  contains

    !=========================================================================
    subroutine coarsen_block_to_buffer

      integer :: i, j, k, iVar, iBufferS
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
    subroutine refine_block_to_buffer

    end subroutine refine_block_to_buffer
    !=========================================================================
    subroutine buffer_to_coarse_block
    end subroutine buffer_to_coarse_block
    !=========================================================================
    subroutine buffer_to_fine_block
    end subroutine buffer_to_fine_block

  end subroutine do_amr

end module BATL_amr
