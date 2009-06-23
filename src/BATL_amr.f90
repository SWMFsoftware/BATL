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
         nI, nJ, nK, nIJK, iRatio, jRatio, kRatio
    use BATL_mpi,  ONLY: iComm, nProc, iProc

    use BATL_tree, ONLY: nNodeUsed, Unused_BP, &
         iTree_IA, iProcNew_A, Proc_, Block_, Coord1_, Coord2_, Coord3_, &
         Status_, Child1_, ChildLast_, Used_, Refine_, CoarsenNew_

    use ModMpi

    ! Arguments
    integer, intent(in) :: nVar
    real, intent(inout) :: &
         State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

    real,    allocatable :: BufferS_I(:), BufferR_I(:), Coarse_VC(:,:,:,:)
    integer, allocatable :: iBlockAvailable_P(:)

    integer :: iNodeSend, iNodeRecv
    integer :: iProcSend, iProcRecv, iBlockSend, iBlockRecv
    integer :: iChild

    integer:: Status_I(MPI_STATUS_SIZE), iError

    character(len=*), parameter:: NameSub = 'BATL_AMR::do_amr'
    logical, parameter:: DoTestMe = .false.
    !-------------------------------------------------------------------------

    ! DoTestMe = iProc == 0

    if(DoTestMe)write(*,*) NameSub,' starting'

    ! Small arrays are allocated once 
    if(.not.allocated(iBlockAvailable_P))then
       allocate(iBlockAvailable_P(0:nProc-1))
    end if

    allocate(BufferS_I(nVar*nIJK), BufferR_I(nVar*nIJK), &
         Coarse_VC(nVar,nI,nJ,nK))

    ! Set iBlockAvailable_P to first available block
    iBlockAvailable_P = -1
    do iProcRecv = 0, nProc-1
       do iBlockRecv = 1, MaxBlock
          if(Unused_BP(iBlockRecv, iProcRecv))then
             iBlockAvailable_P(iProcRecv) = iBlockRecv
             EXIT
          end if
       end do
    end do

    ! Coarsen and move blocks
    do iNodeRecv = 1, nNodeUsed
       if(iTree_IA(Status_,iNodeRecv) /= CoarsenNew_) CYCLE

       if(DoTestMe) write(*,*)NameSub,' CoarsenNew iNode=',iNodeRecv

       iProcRecv  = iProcNew_A(iNodeRecv)
       iBlockRecv = i_block_available(iProcRecv, iNodeRecv)

       do iChild = Child1_, ChildLast_
          iNodeSend = iTree_IA(iChild,iNodeRecv)

          iProcSend  = iTree_IA(Proc_,iNodeSend)
          iBlockSend = iTree_IA(Block_,iNodeSend)

          if(iProc == iProcSend) call send_coarsened_block
          if(iProc == iProcRecv) call recv_coarsened_block

          call make_block_available(iBlockSend, iProcSend)
       end do
    end do

    ! Move blocks
    do iNodeSend = 1, nNodeUsed

       if(iTree_IA(Status_,iNodeSend) /= Used_) CYCLE

       iProcSend = iTree_IA(Proc_,iNodeSend)
       iProcRecv = iProcNew_A(iNodeSend)

       if(iProcRecv == iProcSend) CYCLE

       iBlockSend = iTree_IA(Block_,iNodeSend)
       iBlockRecv = i_block_available(iProcRecv, iNodeSend)

       if(DoTestMe) write(*,*)NameSub, &
            ' node to move iNode,iProcS/R,iBlockS/R=',&
            iNodeSend, iProcSend, iProcRecv, iBlockSend, iBlockRecv

       if(iProc == iProcSend) call send_block
       if(iProc == iProcRecv) call recv_block

       call make_block_available(iBlockSend, iProcSend)
    end do

    ! Prolong and move blocks
    do iNodeSend = 1, nNodeUsed
       if(iTree_IA(Status_,iNodeSend) /= Refine_) CYCLE

       iProcSend  = iTree_IA(Proc_,iNodeSend)
       iBlockSend = iTree_IA(Block_,iNodeSend)

       do iChild = Child1_, ChildLast_
          iNodeRecv = iTree_IA(iChild,iNodeSend)

          iProcRecv  = iProcNew_A(iNodeRecv)
          iBlockRecv = i_block_available(iProcRecv, iNodeRecv)

          if(iProc == iProcSend) call send_refined_block
          if(iProc == iProcRecv) call recv_refined_block
       end do

       call make_block_available(iBlockSend, iProcSend)
    end do

    deallocate(BufferR_I, BufferS_I, Coarse_VC)

  contains

    !==========================================================================
    integer function i_block_available(iProcRecv, iNodeRecv)

      integer, intent(in):: iProcRecv, iNodeRecv
      integer :: iBlock
      !-----------------------------------------------------------------------
      ! Assign the processor index
      iTree_IA(Proc_,iNodeRecv) = iProcRecv

      ! Find and assign the block index
      iBlock = iBlockAvailable_P(iProcRecv)
      iTree_IA(Block_,iNodeRecv) = iBlock

      ! Return the block index
      i_block_available = iBlock

      ! Make the new block used
      Unused_BP(iBlock,iProcRecv) = .false.

      ! Find next available block
      do
         iBlock = iBlock + 1
         if(iBlock > MaxBlock)EXIT
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

      integer:: iBufferS, i, j, k
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

      integer:: iBufferR, i, j, k
      !------------------------------------------------------------------------

      iBufferR = nIJK*nVar
      call MPI_recv(BufferR_I, iBufferR, MPI_REAL, iProcSend, 1, iComm, &
           Status_I, iError)

      iBufferR = 0
      do k = 1, nK; do j = 1, nJ; do i = 1, nI
         State_VGB(:,i,j,k,iBlockRecv) = BufferR_I(iBufferR+1:iBufferR+nVar)
         iBufferR = iBufferR + nVar
      end do; end do; end do

    end subroutine recv_block

    !=========================================================================
    subroutine send_coarsened_block

      use BATL_size, ONLY: InvIjkRatio

      integer :: i, j, k, iVar, iBufferS
      !-----------------------------------------------------------------------
      iBufferS = 0
      do k = 1, nK, kRatio; do j = 1, nJ, jRatio; do i=1, nI, iRatio
         do iVar = 1, nVar
            BufferS_I(iBufferS+iVar) = InvIjkRatio * &
                 sum(State_VGB(iVar,i:i+iRatio-1,j:j+jRatio-1,k:k+kRatio-1,&
                 iBlockSend))
         end do
         iBufferS = iBufferS + nVar
      end do; end do; end do

      call MPI_send(BufferS_I, iBufferS, MPI_REAL, iProcRecv, 1, iComm, iError)

    end subroutine send_coarsened_block
    !==========================================================================
    subroutine recv_coarsened_block

      use BATL_size, ONLY: IjkRatio

      integer:: iBufferR, iSide, jSide, kSide
      integer:: iMin, jMin, kMin, iMax, jMax, kMax
      integer:: i, j, k
      !----------------------------------------------------------------------

      iBufferR = nIJK*nVar/IjkRatio
      call MPI_recv(BufferR_I, iBufferR, MPI_REAL, iProcSend, 1, iComm, &
           Status_I, iError)

      ! Find the part of the block to be written into
      iSide = modulo(iTree_IA(Coord1_,iNodeSend)-1, iRatio)
      jSide = modulo(iTree_IA(Coord2_,iNodeSend)-1, jRatio)
      kSide = modulo(iTree_IA(Coord3_,iNodeSend)-1, kRatio)

      iMin = 1 + iSide*nI/2; iMax = iMin + nI/iRatio - 1
      jMin = 1 + jSide*nJ/2; jMax = jMin + nJ/jRatio - 1
      kMin = 1 + kSide*nK/2; kMax = kMin + nK/kRatio - 1

      iBufferR = 0
      do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
         State_VGB(:,i,j,k,iBlockRecv) = BufferR_I(iBufferR+1:iBufferR+nVar)
         iBufferR = iBufferR + nVar
      end do; end do; end do

    end subroutine recv_coarsened_block
    !==========================================================================
    subroutine send_refined_block

      integer:: iSide, jSide, kSide
      integer:: iMin, jMin, kMin, iMax, jMax, kMax
      integer:: i, j, k, iBufferS
      !----------------------------------------------------------------------

      ! Find the part of the block to be prolonged
      iSide = modulo(iTree_IA(Coord1_,iNodeRecv)-1, iRatio)
      jSide = modulo(iTree_IA(Coord2_,iNodeRecv)-1, jRatio)
      kSide = modulo(iTree_IA(Coord3_,iNodeRecv)-1, kRatio)

      ! Send parent part of the block with one ghost cell
      if(iRatio == 2)then
         iMin = iSide*nI/2; iMax = iMin + nI/2 + 1
      else
         iMin = 1; iMax = nI
      endif
      if(jRatio == 2)then
         jMin = jSide*nJ/2; jMax = jMin + nJ/2 + 1
      else
         jMin = 1; jMax = nJ
      end if
      if(kRatio == 2)then
         kMin = kSide*nK/2; kMax = kMin + nK/2 + 1
      else
         kMin = 1; kMax = nK
      end if

      iBufferS = 0
      do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
         BufferS_I(iBufferS+1:iBufferS+nVar) = State_VGB(:,i,j,k,iBlockSend)
         iBufferS = iBufferS + nVar
      enddo; end do; end do
      call MPI_send(BufferS_I, iBufferS, MPI_REAL, iProcRecv, 1, iComm, iError)

    end subroutine send_refined_block
    !==========================================================================
    subroutine recv_refined_block

      ! Copy buffer into recv block of State_VGB

      integer:: iSize = 1, jSize = 1, kSize = 1
      integer:: iBufferR, i, j, k

      integer:: iS, jS, kS, iR, kR, jR, Di, Dj, Dk
      !------------------------------------------------------------------------
      iSize = nI/iRatio; if(iRatio == 2) iSize = iSize + 2
      jSize = nJ/jRatio; if(jRatio == 2) jSize = jSize + 2
      kSize = nK/kRatio; if(kRatio == 2) kSize = kSize + 2

      iBufferR = iSize*jSize*kSize*nVar
      call MPI_recv(BufferR_I, iBufferR, MPI_REAL, iProcSend, 1, iComm, &
           Status_I, iError)

      iBufferR = 0
      do kS = 1, kSize; do jS = 1, jSize; do iS = 1, iSize
         Coarse_VC(:,iS,jS,kS) = BufferR_I(iBufferR+1:iBufferR+nVar)
         iBufferR = iBufferR + nVar
      end do; end do; end do

      Dk = 0; if(kRatio == 2) Dk = 3
      Dj = 0; if(jRatio == 2) Dj = 3
      Di = 0; if(iRatio == 2) Di = 3
      
      do kR = 1, nK
         kS = (kR + Dk)/kRatio
         do jR = 1, nJ
            jS = (jR + Dj)/jRatio
            do iR = 1, nI
               iS = (iR + Di)/iRatio
               State_VGB(:,iR,jR,kR,iBlockRecv) = Coarse_VC(:,iS,jS,kS)
            end do
         end do
      end do

    end subroutine recv_refined_block

  end subroutine do_amr

  !============================================================================

  subroutine test_amr

    use BATL_mpi,  ONLY: iProc, nProc, barrier_mpi
    use BATL_size, ONLY: MaxDim, nDim, &
         MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nI, nJ, nK, nBlock
    use BATL_tree, ONLY: init_tree, set_tree_root, refine_tree_node, &
         coarsen_tree_node,  distribute_tree, move_tree, show_tree, &
         clean_tree, iTree_IA, &
         iProcNew_A, Unused_B, nNode
    use BATL_grid, ONLY: init_grid, create_grid, clean_grid, Xyz_DGB
    use BATL_geometry, ONLY: init_geometry

    use ModIoUnit,    ONLY: STDOUT_
    use ModUtilities, ONLY: flush_unit

    integer, parameter:: MaxBlockTest            = 100
    integer, parameter:: nRootTest_D(MaxDim)     = (/3,3,3/)
    logical, parameter:: IsPeriodicTest_D(MaxDim)= .true.
    real, parameter:: DomainMin_D(MaxDim) = (/ 1.0, 10.0, 100.0 /)
    real, parameter:: DomainMax_D(MaxDim) = (/10.0,100.0,1000.0 /)
    real, parameter:: DomainSize_D(MaxDim) = DomainMax_D - DomainMin_D

    integer, parameter:: nVar = nDim
    real, allocatable:: State_VGB(:,:,:,:,:), StateOld_VCB(:,:,:,:,:)

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

    if(DoTestMe) call show_tree('after create_grid')

    allocate(State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlockTest), &
         StateOld_VCB(nVar,nI,nJ,nK,MaxBlockTest))

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       State_VGB(:,:,:,:,iBlock)    = Xyz_DGB(1:nDim,:,:,:,iBlock)
       StateOld_VCB(:,:,:,:,iBlock) = Xyz_DGB(1:nDim,1:nI,1:nJ,1:nK,iBlock)
    end do

    write(*,*)'test prolong and balance'
    call refine_tree_node(1)
    if(DoTestMe) call show_tree('after refine_tree_node')
    call distribute_tree(.false.)
    if(DoTestMe) call show_tree('after distribute_tree(.false.)')

    write(*,*)'iProc, iProcNew_A=',iProc, iProcNew_A(1:nNode)

    call do_amr(nVar, State_VGB)
    if(DoTestMe) call show_tree('after do_amr')
    call move_tree
    if(DoTestMe) call show_tree('after move_tree')
    call show_state


    write(*,*)'test restrict and balance'
    call coarsen_tree_node(1)
    if(DoTestMe) call show_tree('after coarsen_tree_node')
    call distribute_tree(.false.)
    if(DoTestMe) call show_tree('after distribute_tree(.false.)')

    write(*,*)'iProc, iProcNew_A=',iProc, iProcNew_A(1:nNode)

    call do_amr(nVar, State_VGB)
    if(DoTestMe) call show_tree('after do_amr')
    call move_tree
    if(DoTestMe) call show_tree('after move_tree')
    call show_state

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE

       if(any(abs(State_VGB(:,1:nI,1:nJ,1:nK,iBlock) &
            -     StateOld_VCB(:,:,:,:,iBlock)) > 1e-6))then
          write(*,*)NameSub,' error for iProc,iBlock,maxloc=',iProc,iBlock,&
               maxloc(abs(State_VGB(:,1:nI,1:nJ,1:nK,iBlock) &
            -     StateOld_VCB(:,:,:,:,iBlock)))
       end if
    end do


    ! Artificially move node 1 to the last processor
    !iProcNew_A(1) = nProc-1
    !
    !call do_amr(nVar, State_VGB)
    !
    !if(DoTestMe) call show_tree('after do_amr')
    !
    !call move_tree
    !
    !if(DoTestMe) call show_tree('after move_tree')
    !
    !call show_state

    deallocate(State_VGB, StateOld_VCB)

    call clean_grid
    call clean_tree

  contains

    subroutine show_state

      integer :: iBlock, iDim, iProcShow

      do iProcShow = 0, nProc - 1
         if(iProc == iProcShow)then
            do iBlock = 1, nBlock
               if(Unused_B(iBlock)) CYCLE
               write(*,'(a,2i4,100f8.4)') &
                    'iProc, iBlock, State(1,1:nI,1,1,iBlock)=', &
                    iProc, iBlock, State_VGB(1,1:nI,1,1,iBlock)

               if(nDim > 1) write(*,'(a,2i4,100f8.3)') &
                    'iProc, iBlock, State(2,1,1:nJ,1,iBlock)=', &
                    iProc, iBlock, State_VGB(2,1,1:nJ,1,iBlock)

               if(nDim > 2) write(*,'(a,2i4,100f8.2)') &
                    'iProc, iBlock, State(3,1,1,1:nK,iBlock)=', &
                    iProc, iBlock, State_VGB(3,1,1,1:nK,iBlock)
            end do
            call flush_unit(STDOUT_)
         end if
         call barrier_mpi
      end do

    end subroutine show_state

  end subroutine test_amr

end module BATL_amr
