!^CFG COPYRIGHT UM

module BATL_pass_cell

  use BATL_size, ONLY: MaxBlock, &
       nI, nJ, nK, nIjk_D, nG, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
       MaxDim, nDim, nDimAmr

  use BATL_mpi, ONLY: iComm, nProc, iProc, barrier_mpi

  use BATL_tree, ONLY: nNodeUsed, iNodePeano_I, &
       iNodeNei_IIIA, DiLevelNei_IIIA, &
       iTree_IA, Proc_, Block_, Coord1_, Coord2_, Coord3_, Status_, &
       Unset_, Unused_

  use ModMpi

  implicit none

  SAVE

  private ! except

  public message_pass_cell

  integer, parameter:: Min_=1, Max_=2
  integer:: iEqualS_DII(MaxDim,-1:1,Min_:Max_)
  integer:: iEqualR_DII(MaxDim,-1:1,Min_:Max_)
  integer:: iRestrictS_DII(MaxDim,-1:1,Min_:Max_)
  integer:: iRestrictR_DII(MaxDim,0:3,Min_:Max_)
  integer:: iProlongS_DII(MaxDim,0:3,Min_:Max_)
  integer:: iProlongR_DII(MaxDim,-1:1,Min_:Max_)

  integer :: nWidth

  ! Refinement ratios in the 3 dimensions (depends on nDimAmr!)
  integer, parameter:: &
       iRatio = 2, jRatio = min(2,nDimAmr), kRatio = max(1,nDimAmr-1)

  integer, parameter:: iRatio_D(MaxDim) = (/ iRatio, jRatio, kRatio /)
  
  ! Inverse volume ratio for Cartesian case !!!
  real, parameter:: InvIjkRatio = 1.0/(iRatio*jRatio*kRatio)

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
    integer :: iSend, jSend, kSend, iRecv, jRecv, kRecv, DiRecv, DjRecv, DkRecv
    integer :: iDir, jDir, kDir
    integer :: iNodeRecv, iNodeSend
    integer :: iBlockRecv, iProcRecv, iBlockSend, iProcSend, DiLevel

    ! Index range for recv and send segments of the blocks
    integer :: iRMin, iRMax, jRMin, jRMax, kRMin, kRMax
    integer :: iSMin, iSMax, jSMin, jSMax, kSMin, kSMax

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

          !write(*,*)'!!! iSendRecv =',iSendRecv

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

             ! For sideways communication from a fine to a coarser block
             ! the coordinate parity of the sender block tells 
             ! if the receiver block fills into the 
             ! lower (D*Recv = 0) or uppper (D*Rev=1) half of the block
             DiRecv = modulo(iTree_IA(Coord1_,iNodeSend)-1, 2)
             DjRecv = modulo(iTree_IA(Coord2_,iNodeSend)-1, 2)
             DkRecv = modulo(iTree_IA(Coord3_,iNodeSend)-1, 2)

             !write(*,*)'!!! DiRecv, DjRecv, DkRecv=', DiRecv, DjRecv, DkRecv 

             do kSend = 0, 3
                ! Skip ignored dimension
                if(nDim < 3 .and. kSend /= 1) CYCLE

                ! Skip second subface if no AMR in this dimension
                if(nDimAmr < 3 .and. kSend == 2) CYCLE

                ! Magic formulas
                ! Direction -1,0,1 can be obtained for part index 0,1,2,3
                kDir = (2*kSend - 3)/3

                ! Receiving part is the opposite of the sending part
                kRecv = 3 - kSend
                if(kRecv == 2) kRecv = 1 + DkRecv

                do jSend = 0, 3
                   if(nDim < 2 .and. jSend /= 1) CYCLE
                   if(nDimAmr < 2 .and. jSend == 2) CYCLE

                   ! Skip edges
                   if(.not.DoSendCorner .and. &
                        (jSend == 0 .or. jSend == 3) .and. &
                        (kSend == 0 .or. kSend == 3)) CYCLE

                   ! Magic formulas
                   jDir  = (2*jSend - 3)/3
                   jRecv = 3 - jSend
                   if(jRecv == 2) jRecv = 1 + DjRecv

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

                      ! Magic formulas
                      iDir  = (2*iSend - 3)/3
                      iRecv = 3 - iSend
                      if(iRecv == 2) iRecv = 1 + DiRecv

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

                      !write(*,*)'!!! iNodeS/R, iProcS/R, iBlockS/R=',&
                      !     iNodeSend, iNodeRecv, iProcSend, iProcRecv, &
                      !     iBlockSend, iBlockRecv
                      !
                      !write(*,*)'!!! iProc,i,j,kSend, i,j,kDir=',&
                      !     iProc,iSend,jSend,kSend,iDir,jDir,kDir

                      if(DiLevel == 0)then
                         call timing_start('do_equal')
                         call do_equal
                         call timing_stop('do_equal')
                      elseif(DiLevel == 1)then
                         ! Do not restrict diagonally in the direction
                         ! of the sibling.
                         if(iDir == -1 .and. DiRecv==1) CYCLE
                         if(iDir == +1 .and. DiRecv==0) CYCLE
                         if(jDir == -1 .and. DjRecv==1) CYCLE
                         if(jDir == +1 .and. DjRecv==0) CYCLE
                         if(kDir == -1 .and. DkRecv==1) CYCLE
                         if(kDir == +1 .and. DkRecv==0) CYCLE

                         call do_restrict
                      else
                         call do_prolong
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

      integer :: iBufferS, iBufferR, i, j, k, nSize
      !------------------------------------------------------------------------

      iRMin = iEqualR_DII(1,iDir,Min_)
      iRMax = iEqualR_DII(1,iDir,Max_)
      jRMin = iEqualR_DII(2,jDir,Min_)
      jRMax = iEqualR_DII(2,jDir,Max_)
      kRMin = iEqualR_DII(3,kDir,Min_)
      kRMax = iEqualR_DII(3,kDir,Max_)

      iSMin = iEqualS_DII(1,iDir,Min_)
      iSMax = iEqualS_DII(1,iDir,Max_)
      jSMin = iEqualS_DII(2,jDir,Min_)
      jSMax = iEqualS_DII(2,jDir,Max_)
      kSMin = iEqualS_DII(3,kDir,Min_)
      kSMax = iEqualS_DII(3,kDir,Max_)

      nSize = nVar*(iRMax-iRMin+1)*(jRMax-jRMin+1)*(kRMax-kRMin+1)

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

    subroutine do_restrict

      integer :: iR, jR, kR, iS1, jS1, kS1, iS2, jS2, kS2, iVar
      integer :: iBufferS, iBufferR, nSize
      !-----------------------------------------------------------------------
      ! Sending range depends on iDir,jDir,kDir only
      iSMin = iRestrictS_DII(1,iDir,Min_)
      iSMax = iRestrictS_DII(1,iDir,Max_)
      jSMin = iRestrictS_DII(2,jDir,Min_)
      jSMax = iRestrictS_DII(2,jDir,Max_)
      kSMin = iRestrictS_DII(3,kDir,Min_)
      kSMax = iRestrictS_DII(3,kDir,Max_)

      ! Receiving range depends on iRecv,kRecv,jRecv = 0..3
      iRMin = iRestrictR_DII(1,iRecv,Min_)
      iRMax = iRestrictR_DII(1,iRecv,Max_)
      jRMin = iRestrictR_DII(2,jRecv,Min_)
      jRMax = iRestrictR_DII(2,jRecv,Max_)
      kRMin = iRestrictR_DII(3,kRecv,Min_)
      kRMax = iRestrictR_DII(3,kRecv,Max_)

      nSize = nVar*(iRMax-iRMin+1)*(jRMax-jRMin+1)*(kRMax-kRMin+1)

      !write(*,*)'iSMin,iSmax,jSMin,jSMax,kSMin,kSmax=',&
      !     iSMin,iSmax,jSMin,jSMax,kSMin,kSmax
      !
      !write(*,*)'iRMin,iRmax,jRMin,jRMax,kRMin,kRmax=',&
      !     iRMin,iRmax,jRMin,jRMax,kRMin,kRmax
      
      if(iSendRecv == Send_)then
         if(iProcSend == iProcRecv)then
            do kR=kRMin,kRMax
               kS1 = kSMin + kRatio*(kR-kRMin)
               kS2 = kS1 + kRatio - 1
               do jR=jRMin,jRMax
                  jS1 = jSMin + jRatio*(jR-jRMin)
                  jS2 = jS1 + jRatio - 1
                  do iR=iRMin,iRMax
                     iS1 = iSMin + iRatio*(iR-iRMin)
                     iS2 = iS1 + iRatio - 1
                     do iVar = 1, nVar
                        State_VGB(iVar,iR,jR,kR,iBlockRecv) = &
                             InvIjkRatio * sum &
                             (State_VGB(iVar,iS1:iS2,jS1:jS2,kS1:kS2,iBlockSend))
                     end do
                  end do
               end do
            end do
         else
            iBufferS = iBufferS_P(iProcRecv)
            if(iProc == iProcSend)then
               do kR=kRMin,kRMax
                  kS1 = kSMin + kRatio*(kR-kRMin)
                  kS2 = kS1 + kRatio - 1
                  do jR=jRMin,jRMax
                     jS1 = jSMin + jRatio*(jR-jRMin)
                     jS2 = jS1 + jRatio - 1
                     do iR=iRMin,iRMax
                        iS1 = iSMin + iRatio*(iR-iRMin)
                        iS2 = iS1 + iRatio - 1
                        do iVar = 1, nVar
                           BufferS_IP(iBufferS+iVar,iProcRecv) = &
                                InvIjkRatio * &
                                sum(State_VGB(iVar,iS1:iS2,jS1:jS2,kS1:kS2,&
                                iBlockSend))
                        end do
                        iBufferS = iBufferS + nVar
                     end do
                  end do
               end do
               iBufferS_P(iProcRecv) = iBufferS
            else
               iBufferR_P(iProcSend) = iBufferR_P(iProcSend) + nSize
            end if
         end if 
      elseif(iProc == iProcRecv .and. iProc /= iProcSend)then ! Receive stage
         iBufferR = iBufferR_P(iProcSend)
         do kR=kRMin,kRmax; do jR=jRMin,jRMax; do iR=iRMin,iRmax
            State_VGB(:,iR,jR,kR,iBlockRecv) = &
                 BufferR_IP(iBufferR+1:iBufferR+nVar,iProcSend)
            iBufferR = iBufferR + nVar
         end do; end do; end do
         iBufferR_P(iProcSend) = iBufferR
      end if

    end subroutine do_restrict

    !==========================================================================

    subroutine do_prolong

      integer :: iR, jR, kR, iS, jS, kS
      integer :: iBufferS, iBufferR, nSize
      !-----------------------------------------------------------------------
      ! Sending range depends on iDir,jDir,kDir only
      iSMin = iProlongS_DII(1,iSend,Min_)
      iSMax = iProlongS_DII(1,iSend,Max_)
      jSMin = iProlongS_DII(2,jSend,Min_)
      jSMax = iProlongS_DII(2,jSend,Max_)
      kSMin = iProlongS_DII(3,kSend,Min_)
      kSMax = iProlongS_DII(3,kSend,Max_)

      ! Receiving range depends on iRecv,kRecv,jRecv = 0..3
      iRMin = iProlongR_DII(1,iDir,Min_)
      iRMax = iProlongR_DII(1,iDir,Max_)
      jRMin = iProlongR_DII(2,jDir,Min_)
      jRMax = iProlongR_DII(2,jDir,Max_)
      kRMin = iProlongR_DII(3,kDir,Min_)
      kRMax = iProlongR_DII(3,kDir,Max_)

      nSize = nVar*(iRMax-iRMin+1)*(jRMax-jRMin+1)*(kRMax-kRMin+1)

      !write(*,*)'iSMin,iSmax,jSMin,jSMax,kSMin,kSmax=',&
      !     iSMin,iSmax,jSMin,jSMax,kSMin,kSmax
      !
      !write(*,*)'iRMin,iRmax,jRMin,jRMax,kRMin,kRmax=',&
      !     iRMin,iRmax,jRMin,jRMax,kRMin,kRmax

      if(iSendRecv == Send_)then
         if(iProcSend == iProcRecv)then
            do kR=kRMin,kRMax
               kS = kSMin + (kR-kRMin)/kRatio
               do jR=jRMin,jRMax
                  jS = jSMin + (jR-jRMin)/jRatio
                  do iR=iRMin,iRMax
                     iS = iSMin + (iR-iRMin)/iRatio
                     State_VGB(:,iR,jR,kR,iBlockRecv) = &
                          State_VGB(:,iS,jS,kS,iBlockSend)
                  end do
               end do
            end do
         else
            iBufferS = iBufferS_P(iProcRecv)
            if(iProc == iProcSend)then
               do kR=kRMin,kRMax
                  kS = kSMin + (kR-kRMin)/kRatio
                  do jR=jRMin,jRMax
                     jS = jSMin + (jR-jRMin)/jRatio
                     do iR=iRMin,iRMax
                        iS = iSMin + (iR-iRMin)/iRatio
                        BufferS_IP(iBufferS+1:iBufferS+nVar,iProcRecv) = &
                             State_VGB(:,iS,jS,kS,iBlockSend)
                        iBufferS = iBufferS + nVar
                     end do
                  end do
               end do
               iBufferS_P(iProcRecv) = iBufferS
            else
               iBufferR_P(iProcSend) = iBufferR_P(iProcSend) + nSize
            end if
         end if
      elseif(iProc == iProcRecv .and. iProc /= iProcSend)then ! Receive stage
         iBufferR = iBufferR_P(iProcSend)

         do kR=kRMin,kRmax; do jR=jRMin,jRMax; do iR=iRMin,iRmax
            State_VGB(:,iR,jR,kR,iBlockRecv) = &
                 BufferR_IP(iBufferR+1:iBufferR+nVar,iProcSend)
            iBufferR = iBufferR + nVar
         end do; end do; end do
         iBufferR_P(iProcSend) = iBufferR
      end if

    end subroutine do_prolong

  end subroutine message_pass_cell

  !============================================================================
  subroutine set_range

    ! Indexed by iDir/jDir/kDir for sender = -1,0,1
    iEqualS_DII(:,-1,Min_) = 1
    iEqualS_DII(:,-1,Max_) = nWidth
    iEqualS_DII(:, 0,Min_) = 1
    iEqualS_DII(:, 0,Max_) = nIjk_D
    iEqualS_DII(:, 1,Min_) = nIjk_D + 1 - nWidth
    iEqualS_DII(:, 1,Max_) = nIjk_D

    ! Indexed by iDir/jDir/kDir for sender = -1,0,1
    iEqualR_DII(:,-1,Min_) = nIjk_D + 1
    iEqualR_DII(:,-1,Max_) = nIjk_D + nWidth
    iEqualR_DII(:, 0,Min_) = 1
    iEqualR_DII(:, 0,Max_) = nIjk_D
    iEqualR_DII(:, 1,Min_) = 1 - nWidth
    iEqualR_DII(:, 1,Max_) = 0

    ! Indexed by iDir/jDir/kDir for sender = -1,0,1
    iRestrictS_DII(:,-1,Min_) = 1
    iRestrictS_DII(:,-1,Max_) = iRatio_D*nWidth
    iRestrictS_DII(:, 0,Min_) = 1
    iRestrictS_DII(:, 0,Max_) = nIjk_D
    iRestrictS_DII(:, 1,Min_) = nIjk_D + 1 - iRatio_D*nWidth
    iRestrictS_DII(:, 1,Max_) = nIjk_D

    ! Indexed by iRecv/jRecv/kRecv = 0..3
    iRestrictR_DII(:,0,Min_) = 1 - nWidth
    iRestrictR_DII(:,0,Max_) = 0
    iRestrictR_DII(:,1,Min_) = 1
    iRestrictR_DII(:,1,Max_) = nIjk_D/iRatio_D
    iRestrictR_DII(:,2,Min_) = nIjk_D/iRatio_D + 1
    iRestrictR_DII(:,2,Max_) = nIjk_D
    iRestrictR_DII(:,3,Min_) = nIjk_D + 1
    iRestrictR_DII(:,3,Max_) = nIjk_D + nWidth

    ! Indexed by iSend/jSend,kSend = 0..3
    iProlongS_DII(:,0,Min_) = 1
    iProlongS_DII(:,0,Max_) = 1 + (nWidth-1)/iRatio_D ! rounding up
    iProlongS_DII(:,1,Min_) = 1
    iProlongS_DII(:,1,Max_) = nIjk_D/iRatio_D
    iProlongS_DII(:,2,Min_) = nIjk_D/iRatio_D + 1
    iProlongS_DII(:,2,Max_) = nIjk_D
    iProlongS_DII(:,3,Min_) = nIjk_D + 1 - (1+(nWidth-1)/iRatio_D)
    iProlongS_DII(:,3,Max_) = nIjk_D
    
    ! Indexed by iDir/jDir/kDir for sender = -1,0,1
    iProlongR_DII(:,-1,Min_) = nIjk_D + 1
    iProlongR_DII(:,-1,Max_) = nIjk_D + nWidth
    iProlongR_DII(:, 0,Min_) = 1
    iProlongR_DII(:, 0,Max_) = nIjk_D
    iProlongR_DII(:, 1,Min_) = 1 - nWidth
    iProlongR_DII(:, 1,Max_) = 0

  end subroutine set_range

end module BATL_pass_cell
