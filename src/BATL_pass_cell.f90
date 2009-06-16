!^CFG COPYRIGHT UM

module BATL_pass_cell

  use BATL_size, ONLY: MaxBlock, &
       nI, nJ, nK, nIjk_D, nG, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
       nDim, nDimAmr

  use BATL_mpi, ONLY: iComm, nProc, iProc, barrier_mpi

  use BATL_tree, ONLY: nNodeUsed, iNodePeano_I, &
       iNodeNei_IIIA, DiLevelNei_IIIA, &
       iTree_IA, Proc_, Block_, Status_, Unset_, Unused_

  use ModMpi

  implicit none

  SAVE

  private ! except

  public message_pass_cell

  integer, parameter:: RecvMin_=1, RecvMax_=2, SendMin_=3, SendMax_=4
  integer:: iEqual_DII(3,-1:1,RecvMin_:SendMax_)

  integer :: nWidth

contains

  subroutine message_pass_cell(nVar, State_VGB, &
       nWidthIn, nProlongOrderIn, nCoarseLayerIn, DoSendCornerIn, &
       DoRestrictFaceIn)

    ! Arguments
    integer, intent(in) :: nVar
    real, intent(inout) :: &
         State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

    ! Optional arguments
    integer, optional, intent(in) :: nWidthIn
    integer, optional, intent(in) :: nCoarseLayerIn
    integer, optional, intent(in) :: nProlongOrderIn
    logical, optional, intent(in) :: DoSendCornerIn
    logical, optional, intent(in) :: DoRestrictFaceIn

    ! Local variables

    logical :: DoSendCorner
    integer :: nProlongOrder
    integer :: nCoarseLayer
    logical :: DoRestrictFace
    character(len=*), parameter:: NameSub = 'BATL_pass_cell::message_pass_cell'

    ! Named indexes for the send and receive stages 
    integer, parameter :: Send_=1, Recv_=2

    integer :: iProlongOrder, iSendRecv, iPeano
    integer :: iSend, jSend, kSend
    integer :: iDir, jDir, kDir
    integer :: iNodeRecv, iNodeSend
    integer :: iBlockRecv, iProcRecv, iBlockSend, iProcSend, DiLevel

    ! Variables related to recv and send buffers
    integer, parameter:: nGhostCell = &
         (MaxI-MinI+1)*(MaxJ-MinJ+1)*(MaxK-MinK+1) - nI*nJ*nK
    integer :: MaxBuffer
    integer, allocatable:: iBufferR_P(:), iBufferS_P(:)
    real, allocatable:: BufferR_IP(:,:), BufferS_IP(:,:)

    integer:: iRequestR, iRequestS, iError
    integer, allocatable:: iRequestR_I(:), iRequestS_I(:), iStatus_II(:,:)

    logical :: DoTest = .false., DoTestMe = .false.
    !--------------------------------------------------------------------------
    if(DoTestMe)write(*,*)NameSub,' starting with nVar=',nVar

    ! Set values or defaults for optional arguments
    nWidth = nG
    if(present(nWidthIn)) nWidth = nWidthIn

!!!    nProlongOrder = 2
    nProlongOrder = 1
    if(present(nProlongOrderIn)) nProlongOrder = nProlongOrderIn

    nCoarseLayer = 1
    if(present(nCoarseLayerIn)) nCoarseLayer = nCoarseLayerIn

!!!    DoSendCorner = .true.
    DoSendCorner = .false.
    if(present(DoSendCornerIn)) DoSendCorner = DoSendCornerIn

    DoRestrictFace = .false.
    if(present(DoRestrictFaceIn)) DoRestrictFace = DoRestrictFaceIn

    ! if(nProc>1)call show_tree('BATL_pass_cell',.true.)

    call set_range

    if(nProc > 1)then
       ! Small arrays are allocated once 
       if(.not.allocated(iBufferR_P))then
          allocate(iBufferR_P(0:nProc-1), iBufferS_P(0:nProc-1))
          allocate(iRequestR_I(nProc-1), iRequestS_I(nProc-1))
          allocate(iStatus_II(MPI_STATUS_SIZE,nProc-1))
       end if
       ! Big arrays are allocated and deallocated every time (for now)
       MaxBuffer = nVar*MaxBlock*nGhostCell
       call timing_start('msg_allocate')
       allocate(BufferR_IP(MaxBuffer,0:nProc-1))
       allocate(BufferS_IP(MaxBuffer,0:nProc-1))
       call timing_stop('msg_allocate')
    end if

    do iProlongOrder = 1, nProlongOrder
       ! Second order prolongation needs two stages: 
       ! first stage fills in equal and coarser ghost cells
       ! second stage uses these to prolong and fill in finer ghost cells

       do iSendRecv = Send_, Recv_
          ! For serial run there is no need for Recv_ stage
          if(nProc==1 .and. iSendRecv == Recv_) CYCLE

          if(nProc>1)then
             ! initialize buffer indexes
             iBufferR_P = 0
             iBufferS_P = 0
          end if

          ! Loop through all nodes
          do iPeano = 1, nNodeUsed

             iNodeSend = iNodePeano_I(iPeano)

             iBlockSend = iTree_IA(Block_,iNodeSend)
             iProcSend  = iTree_IA(Proc_, iNodeSend) 

             do kSend = 0, 3
                ! Skip ignored dimension
                if(nDim < 3 .and. kSend /= 1) CYCLE

                ! Skip second subface if no AMR in this dimension
                if(nDimAmr < 3 .and. kSend == 2) CYCLE

                kDir = (2*kSend - 3)/3

                do jSend = 0, 3
                   if(nDim < 2 .and. jSend /= 1) CYCLE
                   if(nDimAmr < 2 .and. jSend == 2) CYCLE

                   ! Skip edges
                   if(.not.DoSendCorner .and. &
                        (jSend == 0 .or. jSend == 3) .and. &
                        (kSend == 0 .or. kSend == 3)) CYCLE

                   jDir = (2*jSend - 3)/3

                   do iSend = 0, 3

                      ! Ignore inner points
                      if(  0 < iSend .and. iSend < 3 .and. &
                           0 < jSend .and. jSend < 3 .and. &
                           0 < kSend .and. kSend < 3) CYCLE

                      ! Exclude corners where i and j or k is at the edge
                      if(.not.DoSendCorner .and. &
                           (iSend == 0 .or. iSend == 3) .and. &
                           (jSend == 0 .or. jSend == 3 .or. &
                           kSend == 0 .or. kSend ==3)) CYCLE

                      iDir = (2*iSend - 3)/3

                      iNodeRecv  = iNodeNei_IIIA(iSend,jSend,kSend,iNodeSend)
                      iProcRecv  = iTree_IA(Proc_,iNodeRecv)

                      if(iProc /= iProcSend .and. iProc /= iProcRecv) CYCLE

                      iBlockRecv = iTree_IA(Block_,iNodeRecv)

                      DiLevel = DiLevelNei_IIIA(iDir,jDir,kDir,iNodeSend)
                      ! Only 1 "subface" if the level is equal or coarser
                      if(DiLevel >= 0 .and. &
                           (iSend==2 .or. jSend==2 .or. kSend==2)) CYCLE

                      if(nProlongOrder == 2 .and. &
                           (iProlongOrder == 1 .and. DiLevel == -1) .or. &
                           (iProlongOrder == 2 .and. DiLevel >=  0)) CYCLE
                      
                      if(DiLevel == 0)then
                         call timing_start('do_equal')
                         call do_equal
                         call timing_stop('do_equal')
                      elseif(DiLevel == 1)then
                         call do_coarse
                      else
                         call do_fine
                      end if
                   enddo ! iSend
                enddo ! jSend
             enddo ! kSend
          end do ! iNodeSend

          ! Messages are sent at the end of the Send_ stage only
          if(iSendRecv == Recv_ .or. nProc == 1) EXIT

          call timing_start('send_recv')
          ! post requests
          iRequestR = 0
          do iProcSend = 0, nProc - 1
             if(iBufferR_P(iProcSend) == 0) CYCLE
             iRequestR = iRequestR + 1

             !write(*,*)'!!! MPI_irecv:iProcRecv,iProcSend,iBufferR=',&
             !     iProc,iProcSend,iBufferR_P(iProcSend)

             call MPI_irecv(BufferR_IP(1,iProcSend), iBufferR_P(iProcSend), &
                  MPI_REAL, iProcSend, 1, iComm, iRequestR_I(iRequestR), &
                  iError)
          end do

          ! post sends
          iRequestS = 0
          do iProcRecv = 0, nProc-1
             if(iBufferS_P(iProcRecv) == 0) CYCLE
             iRequestS = iRequestS + 1

             !write(*,*)'!!! MPI_isend:iProcRecv,iProcSend,iBufferS=',&
             !     iProcRecv,iProc,iBufferS_P(iProcRecv)

             call MPI_isend(BufferS_IP(1,iProcRecv), iBufferS_P(iProcRecv), &
                  MPI_REAL, iProcRecv, 1, iComm, iRequestS_I(iRequestS), &
                  iError)
          end do

          ! wait for all requests to be completed
          if(iRequestR > 0) &
               call MPI_waitall(iRequestR, iRequestR_I, iStatus_II, iError)

          ! wait for all sends to be completed
          if(iRequestS > 0) &
               call MPI_waitall(iRequestS, iRequestS_I, iStatus_II, iError)

          call timing_stop('send_recv')
          
       end do ! iSendRecv
    end do ! iProlongOrder

    if(nProc>1)deallocate(BufferR_IP, BufferS_IP)

  contains

    !==========================================================================
    subroutine do_equal

      integer :: iRMin, iRMax, jRMin, jRMax, kRMin, kRMax
      integer :: iSMin, iSMax, jSMin, jSMax, kSMin, kSMax
      integer :: iBufferS, iBufferR, i, j, k, nSize
      !------------------------------------------------------------------------

      iRMin = iEqual_DII(1,iDir,RecvMin_)
      iRMax = iEqual_DII(1,iDir,RecvMax_)
      jRMin = iEqual_DII(2,jDir,RecvMin_)
      jRMax = iEqual_DII(2,jDir,RecvMax_)
      kRMin = iEqual_DII(3,kDir,RecvMin_)
      kRMax = iEqual_DII(3,kDir,RecvMax_)

      iSMin = iEqual_DII(1,iDir,SendMin_)
      iSMax = iEqual_DII(1,iDir,SendMax_)
      jSMin = iEqual_DII(2,jDir,SendMin_)
      jSMax = iEqual_DII(2,jDir,SendMax_)
      kSMin = iEqual_DII(3,kDir,SendMin_)
      kSMax = iEqual_DII(3,kDir,SendMax_)

      nSize = nVar*(iSMax-iSMin+1)*(jSMax-jSMin+1)*(kSMax-kSMin+1)

      !if(nProc>1)then
      !   write(*,*)'!!! iProc,iSendRecv =',iProc,iSendRecv
      !   write(*,*)'!!! iProc,i,j,kSend, i,j,kDir=',&
      !                        iProc,iSend,jSend,kSend,iDir,jDir,kDir
      !
      !   write(*,*)'!!! do_equal: iNodeS/R, iProcS/R, iBlockS/R, nSize=',&
      !     iNodeSend, iNodeRecv, iProcSend, iProcRecv, iBlockSend, iBlockRecv,&
      !     nSize
      !end if

      if(iSendRecv == Send_)then
         if(iProcSend == iProcRecv)then
            State_VGB(:,iRMin:iRMax,jRMin:jRMax,kRMin:kRMax,iBlockRecv)= &
                 State_VGB(:,iSMin:iSMax,jSMin:jSMax,kSMin:kSMax,iBlockSend)
         else
            iBufferS = iBufferS_P(iProcRecv)
            if(iProc == iProcSend)then
               do k=kSMin,kSmax; do j=jSMin,jSMax; do i=iSMin,iSmax
                  BufferS_IP(iBufferS+1:iBufferS+nVar,iProcRecv) = &
                       State_VGB(:,i,j,k,iBlockSend)
                  iBufferS = iBufferS + nVar
               end do; end do; end do
               iBufferS_P(iProcRecv) = iBufferS
            else
               iBufferR_P(iProcSend) = iBufferR_P(iProcSend) + nSize
            end if
         end if
      elseif(iProc == iProcRecv .and. iProc /= iProcSend)then ! Receive stage
         iBufferR = iBufferR_P(iProcSend)

         do k=kRMin,kRmax; do j=jRMin,jRMax; do i=iRMin,iRmax
            State_VGB(:,i,j,k,iBlockRecv) = &
                 BufferR_IP(iBufferR+1:iBufferR+nVar,iProcSend)
            iBufferR = iBufferR + nVar
         end do; end do; end do
         iBufferR_P(iProcSend) = iBufferR

      end if

    end subroutine do_equal

    !==========================================================================

    subroutine do_coarse




    end subroutine do_coarse

    !==========================================================================

    subroutine do_fine

    end subroutine do_fine

  end subroutine message_pass_cell

  !============================================================================
  subroutine set_range

    iEqual_DII(:,-1,SendMin_) = 1
    iEqual_DII(:,-1,SendMax_) = nWidth
    iEqual_DII(:, 0,SendMin_) = 1
    iEqual_DII(:, 0,SendMax_) = nIjk_D
    iEqual_DII(:, 1,SendMin_) = nIjk_D + 1 - nWidth
    iEqual_DII(:, 1,SendMax_) = nIjk_D

    iEqual_DII(:,-1,RecvMin_) = nIjk_D + 1
    iEqual_DII(:,-1,RecvMax_) = nIjk_D + nWidth
    iEqual_DII(:, 0,RecvMin_) = 1
    iEqual_DII(:, 0,RecvMax_) = nIjk_D
    iEqual_DII(:, 1,RecvMin_) = 1 - nWidth
    iEqual_DII(:, 1,RecvMax_) = 0

  end subroutine set_range

end module BATL_pass_cell
