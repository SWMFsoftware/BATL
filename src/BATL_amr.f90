!^CFG COPYRIGHT UM

module BATL_amr

  implicit none

  SAVE

  private ! except

  public do_amr
  public test_amr

contains

  !===========================================================================
  subroutine do_amr(nVar, State_VGB)

    use BATL_size, ONLY: MaxBlock, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
         nI, nJ, nK, nIJK
    use BATL_mpi,  ONLY: iComm, nProc, iProc

    use BATL_tree, ONLY: nNodeUsed, Unused_B, &
         iTree_IA, Proc_, Block_, ProcNew_, &
         Status_, Child1_, ChildLast_, Used_, Refine_, CoarsenNew_

    use ModMpi

    ! Arguments
    integer, intent(in) :: nVar
    real, intent(inout) :: &
         State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

    real,    allocatable :: BufferS_I(:), BufferR_I(:)
    integer, allocatable :: iBlockAvailable_P(:)
    logical, allocatable :: Unused_BP(:,:)

    integer :: iNodeSend, iNodeRecv
    integer :: iProcSend, iProcRecv, iBlockSend, iBlockRecv
    integer :: iChild, iError
    !-------------------------------------------------------------------------

    ! Small arrays are allocated once 
    if(.not.allocated(iBlockAvailable_P))then
       allocate(iBlockAvailable_P(0:nProc-1))
       allocate(Unused_BP(MaxBlock,0:nProc-1))
    end if

    allocate(BufferS_I(nVar*nIJK), BufferR_I(nVar*nIJK))

    call MPI_allgather(Unused_B, MaxBlock, MPI_LOGICAL, &
         Unused_BP, MaxBlock, MPI_LOGICAL, iComm, iError)

    ! Coarsen and move blocks
    do iNodeRecv = 1, nNodeUsed
       if(iTree_IA(Status_,iNodeSend) /= CoarsenNew_) CYCLE

       iProcRecv  = iTree_IA(Proc_,iNodeRecv)
       iBlockRecv = i_block_available(iProcRecv)

       do iChild = Child1_, ChildLast_
          iNodeSend = iTree_IA(iChild,iNodeRecv)

          iProcSend  = iTree_IA(Proc_,iNodeSend)
          iBlockSend = iTree_IA(Block_,iNodeSend)

          if(iProc == iProcSend) call send_coarsened_block

          call make_block_available(iBlockSend, iProcSend)
       end do
       if(iProc == iProcRecv) call recv_coarsened_block
    end do

    ! Move blocks
    do iNodeSend = 1, nNodeUsed
       if(iTree_IA(Status_,iNodeSend) /= Used_) CYCLE

       iProcSend = iTree_IA(Status_,Proc_)
       iProcRecv = iTree_IA(Status_,ProcNew_)

       if(iProcRecv == iProcSend) CYCLE

       iBlockSend = iTree_IA(Block_,iNodeSend)
       iBlockRecv = i_block_available(iProcRecv)

       if(iProc == iProcSend) call send_block
       if(iProc == iProcRecv) call recv_block

       call make_block_available(iBlockSend, iProcSend)
    end do

    ! Prolong and move blocks
    do iNodeSend = 1, nNodeUsed
       if(iTree_IA(Status_,iNodeSend) /= Refine_) CYCLE

       iProcSend  = iTree_IA(Proc_,iNodeSend)
       iBlockSend = iTree_IA(Block_,iNodeSend)

       if(iProc == iProcSend) call send_refined_block

       do iChild = Child1_, ChildLast_
          iNodeRecv = iTree_IA(iChild,iNodeSend)

          iProcRecv  = iTree_IA(Proc_,iNodeRecv)
          iBlockRecv = i_block_available(iProcRecv)

          if(iProc == iProcRecv) call recv_refined_block
       end do

       call make_block_available(iBlockSend, iProcSend)
    end do

    deallocate(BufferR_I, BufferS_I)

  contains

    !==========================================================================
    integer function i_block_available(iProcRecv)

      integer, intent(in):: iProcRecv
      integer :: iBlock
      !-----------------------------------------------------------------------
      iBlock = iBlockAvailable_P(iProcRecv)

      i_block_available = iBlock
      Unused_BP(iBlock,iProcRecv) = .false.
      do
         iBlock = iBlock + 1
         if(Unused_BP(iBlock,iProcRecv))EXIT
      end do

      iBlockAvailable_P(iProcRecv) = iBlock

    end function i_block_available
    !==========================================================================
    subroutine make_block_available(iBlockSend, iProcSend)

      integer, intent(in):: iBlockSend, iProcSend
      !-----------------------------------------------------------------------
      Unused_BP(iBlockSend, iProcSend) = .true.
      iBlockAvailable_P(iProcSend) = &
           min(iBlockAvailable_P(iProcSend), iBlockSend)

    end subroutine make_block_available
    !==========================================================================
    subroutine send_block

      ! Copy buffer into recv block of State_VGB

      integer:: iBufferS, i, j, k, iError
      !------------------------------------------------------------------------

      iBufferS = 0
      do k = 1, nK; do j = 1, nJ; do i = 1, nI
         BufferS_I(iBufferS+1:iBufferS+nVar) = State_VGB(:,i,j,k,iBlockSend)
         iBufferS = iBufferS + nVar
      end do; end do; end do

      call MPI_send(BufferS_I, iBufferS, MPI_REAL, iProcRecv, 1, iComm, iError)

    end subroutine send_block
    !==========================================================================
    subroutine recv_block

      ! Copy buffer into recv block of State_VGB

      integer:: iBufferR, i, j, k, iError
      integer:: Status_I(MPI_STATUS_SIZE)
      !------------------------------------------------------------------------

      iBufferR = nIJK*nVar
      call MPI_recv(BufferR_I, iBufferR, MPI_REAL, iProcRecv, 1, iComm, &
           Status_I, iError)

      iBufferR = 0
      do k = 1, nK; do j = 1, nJ; do i = 1, nI
         State_VGB(:,i,j,k,iBlockRecv) = BufferR_I(iBufferR+1:iBufferR+nVar)
         iBufferR = iBufferR + nVar
      end do; end do; end do

    end subroutine recv_block

    !=========================================================================
    subroutine send_coarsened_block

!      use BATL_size, ONLY: iRatio, jRatio, kRatio
!
!      integer :: i, j, k, iVar, iBufferS
!      !-----------------------------------------------------------------------
!      iProcRecv = iTree_IA(ProcNew_,iNodeSend)
!
!      iBufferS = iBufferS_P(iProcRecv)
!      do k = 1, nK, kRatio; do j = 1, nJ, jRatio; do i=1, nI, iRatio
!         do iVar = 1, nVar
!            BufferS_IP(iBufferS+iVar,iProcRecv) = &
!                 InvIjkRatio * &
!                 sum(State_VGB(iVar,i:i+iRatio-1,j:j+jRatio-1,k:k+kRatio-1,&
!                 iBlockSend))
!         end do
!         iBufferS = iBufferS + nVar
!      end do; end do; end do
!
!      iBufferR_P(iProcSend) = iBufferR_P(iProcSend) &
!           + iBufferS - iBufferS_P(iProcRecv)
!      iBufferS_P(iProcRecv) = iBufferS
!
!      call MPI_send(BufferS_I, iBufferS, MPI_REAL, iProcRecv, 1, iComm, iError)

    end subroutine send_coarsened_block
    !==========================================================================
    subroutine send_refined_block
    end subroutine send_refined_block
    !==========================================================================
    subroutine recv_refined_block
    end subroutine recv_refined_block
    !==========================================================================
    subroutine recv_coarsened_block
    end subroutine recv_coarsened_block

  end subroutine do_amr

  !============================================================================

  subroutine test_amr

    use BATL_mpi,  ONLY: iProc, nProc
    use BATL_size, ONLY: MaxDim, nDim, &
         MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nI, nJ, nK, nBlock
    use BATL_tree, ONLY: init_tree, set_tree_root, refine_tree_node, &
         distribute_tree, move_tree, show_tree, clean_tree, iTree_IA, &
         Unused_B, ProcNew_
    use BATL_grid, ONLY: init_grid, create_grid, clean_grid, Xyz_DGB
    use BATL_geometry, ONLY: init_geometry

    integer, parameter:: MaxBlockTest            = 15
    integer, parameter:: nRootTest_D(MaxDim)     = (/2,2,2/)
    logical, parameter:: IsPeriodicTest_D(MaxDim)= .true.
    real, parameter:: DomainMin_D(MaxDim) = (/ 1.0, 10.0, 100.0 /)
    real, parameter:: DomainMax_D(MaxDim) = (/ 2.0, 20.0, 200.0 /)
    real, parameter:: DomainSize_D(MaxDim) = DomainMax_D - DomainMin_D

    integer, parameter:: nVar = nDim
    real, allocatable:: State_VGB(:,:,:,:,:)

    integer:: iBlock

    logical:: DoTestMe
    character(len=*), parameter :: NameSub = 'test_amr'
    !-----------------------------------------------------------------------
    DoTestMe = iProc == 0

    if(DoTestMe) write(*,*) 'Starting ',NameSub

    call init_tree(MaxBlockTest)
    call init_grid( DomainMin_D(1:nDim), DomainMax_D(1:nDim) )
    call init_geometry( IsPeriodicIn_D = IsPeriodicTest_D(1:nDim) )
    call set_tree_root( nRootTest_D(1:nDim))
    call distribute_tree(.true.)
    call create_grid

    if(DoTestMe) call show_tree(NameSub,.true.)

    allocate(State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlockTest))

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       State_VGB(:,:,:,:,iBlock) = &
               Xyz_DGB(1:nDim,:,:,:,iBlock)
    end do

    ! Artificially move node 1
    iTree_IA(ProcNew_,1) = nProc-1

    call do_amr(nVar, State_VGB)

    call move_tree

    if(DoTestMe) call show_tree(NameSub,.true.)

    deallocate(State_VGB)

    call clean_grid
    call clean_tree

  end subroutine test_amr

end module BATL_amr
