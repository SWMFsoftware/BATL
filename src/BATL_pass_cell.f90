!^CFG COPYRIGHT UM
module BATL_pass_cell

  ! Possible improvements:
  ! (1) Instead of sending the receiving block number
  !     and the 2**nDim range limits, we can send only the tag which
  !     we would use in a block to block communication:
  !        iTag = 100*iBlockRecv + iRecv + 4*(jRecv + 4*kRecv)
  !     There are 2 advantages:
  !     (a) The amount of info reduces: 1+2**nDim numbers --> 1 number (iTag)
  !     (b) This procedure allows to send for 1st order prolongation only
  !         one copy per 2**nDim cells
  ! (2) Instead of waiting for receiving buffers from ALL processors, we
  !     we can wait for ANY receiving and already start unpacking
  ! (3) Instead of determining the receive (and send) buffer size during
  !     the message_pass_cell, we can determining the sizes a priori:
  !     (a) We can then allocate a small known buffer size
  !     (b) we do at least two times message_pass_cell per time iteration,
  !         each time determining the buffer size. This would be reduced to
  !         only once (there is a small complication with operatore split
  !         schemes)

  implicit none

  SAVE

  private ! except

  public message_pass_cell
  public test_pass_cell

contains

  subroutine message_pass_cell(nVar, State_VGB, &
       nWidthIn, nProlongOrderIn, nCoarseLayerIn, DoSendCornerIn, &
       DoRestrictFaceIn)

    use BATL_size, ONLY: MaxBlock, &
         nBlock, nI, nJ, nK, nIjk_D, nG, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
         MaxDim, nDim, nDimAmr, iRatio, jRatio, kRatio, iRatio_D, InvIjkRatio

    use BATL_mpi, ONLY: iComm, nProc, iProc

    use BATL_tree, ONLY: nNodeUsed, iNodePeano_I, &
         iNodeNei_IIIB, DiLevelNei_IIIB, Unused_B, iNode_B, &
         iTree_IA, Proc_, Block_, Coord1_, Coord2_, Coord3_

    use ModMpi

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

    integer :: nWidth
    integer :: nCoarseLayer
    integer :: nProlongOrder
    logical :: DoSendCorner
    logical :: DoRestrictFace
    character(len=*), parameter:: NameSub = 'BATL_pass_cell::message_pass_cell'

    integer :: iProlongOrder, iSendRecv, iPeano
    integer :: iSend, jSend, kSend, iRecv, jRecv, kRecv, DiRecv, DjRecv, DkRecv
    integer :: iDir, jDir, kDir
    integer :: iNodeRecv, iNodeSend
    integer :: iBlockRecv, iProcRecv, iBlockSend, iProcSend, DiLevel

    ! Fast lookup tables for index ranges per dimension
    integer, parameter:: Min_=1, Max_=2
    integer:: iEqualS_DII(MaxDim,-1:1,Min_:Max_)
    integer:: iEqualR_DII(MaxDim,-1:1,Min_:Max_)
    integer:: iRestrictS_DII(MaxDim,-1:1,Min_:Max_)
    integer:: iRestrictR_DII(MaxDim,0:3,Min_:Max_)
    integer:: iProlongS_DII(MaxDim,0:3,Min_:Max_)
    integer:: iProlongR_DII(MaxDim,0:3,Min_:Max_)

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

    DoSendCorner = .true.
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

       if(nProc>1)then
          ! initialize buffer indexes
          iBufferR_P = 0
          iBufferS_P = 0
       end if

       ! Loop through all nodes
       do iBlockSend = 1, nBlock
          if(Unused_B(iBlockSend)) CYCLE

          iNodeSend = iNode_B(iBlockSend)

          do kDir = -1, 1
             ! Do not message pass in ignored dimensions
             if(nDim < 3 .and. kDir /= 0) CYCLE

             do jDir = -1, 1
                if(nDim < 2 .and. jDir /= 0) CYCLE
                ! Skip edges
                if(.not.DoSendCorner .and. jDir /= 0 .and. kDir /= 0) CYCLE

                do iDir = -1,1
                   ! Ignore inner parts of the sending block
                   if(iDir == 0 .and. jDir == 0 .and. kDir == 0) CYCLE

                   ! Exclude corners where i and j or k is at the edge
                   if(.not.DoSendCorner .and. iDir /= 0 .and. &
                        (jDir /= 0 .or.  kDir /= 0)) CYCLE

                   DiLevel = DiLevelNei_IIIB(iDir,jDir,kDir,iBlockSend)
                   
                   ! Do prolongation in the second stage if nProlongOrder is 2
                   if(nProlongOrder == 2 .and. &
                        (iProlongOrder == 1 .and. DiLevel == -1) .or. &
                        (iProlongOrder == 2 .and. DiLevel >=  0)) CYCLE

                   if(DiLevel == 0)then
                      call do_equal
                   elseif(DiLevel == 1)then
                      call do_restrict
                   elseif(DiLevel == -1)then
                      call do_prolong
                   endif
                end do
             end do
          end do
       end do ! iBlockSend

       ! Done for serial run
       if(nProc == 1) CYCLE

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

       call timing_start('buffer_to_state')
       call buffer_to_state
       call timing_stop('buffer_to_state')

    end do ! iProlongOrder

    if(nProc>1)deallocate(BufferR_IP, BufferS_IP)

  contains

    !==========================================================================
    subroutine buffer_to_state

      ! Copy buffer into recv block of State_VGB

      integer:: iBufferR, i, j, k
      !------------------------------------------------------------------------

      jRMin = 1; jRMax = 1
      kRMin = 1; kRMax = 1
      do iProcSend = 0, nProc - 1
         if(iBufferR_P(iProcSend) == 0) CYCLE

         iBufferR = 0
         do
            iBlockRecv = nint(BufferR_IP(iBufferR+1,iProcSend))
            iRMin      = nint(BufferR_IP(iBufferR+2,iProcSend))
            iRMax      = nint(BufferR_IP(iBufferR+3,iProcSend))
            if(nDim > 1) jRMin = nint(BufferR_IP(iBufferR+4,iProcSend))
            if(nDim > 1) jRMax = nint(BufferR_IP(iBufferR+5,iProcSend))
            if(nDim > 2) kRMin = nint(BufferR_IP(iBufferR+6,iProcSend))
            if(nDim > 2) kRMax = nint(BufferR_IP(iBufferR+7,iProcSend))

            iBufferR = iBufferR + 1 + 2*nDim

            do k = kRMin, kRmax; do j = jRMin, jRMax; do i = iRMin, iRmax
               State_VGB(:,i,j,k,iBlockRecv) = &
                    BufferR_IP(iBufferR+1:iBufferR+nVar,iProcSend)

               iBufferR = iBufferR + nVar
            end do; end do; end do

            if(iBufferR >= iBufferR_P(iProcSend)) EXIT
         end do
      end do
      
    end subroutine buffer_to_state

    !==========================================================================

    subroutine do_equal

      integer :: iBufferS, i, j, k, nSize
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

      iSend = (3*iDir + 3)/2
      jSend = (3*jDir + 3)/2
      kSend = (3*kDir + 3)/2

      iNodeRecv  = iNodeNei_IIIB(iSend,jSend,kSend,iBlockSend)
      iProcRecv  = iTree_IA(Proc_,iNodeRecv)
      iBlockRecv = iTree_IA(Block_,iNodeRecv)

      if(iProc == iProcRecv)then
         State_VGB(:,iRMin:iRMax,jRMin:jRMax,kRMin:kRMax,iBlockRecv)= &
              State_VGB(:,iSMin:iSMax,jSMin:jSMax,kSMin:kSMax,iBlockSend)
      else
         iBufferS = iBufferS_P(iProcRecv)

         BufferS_IP(            iBufferS+1,iProcRecv) = iBlockRecv
         BufferS_IP(            iBufferS+2,iProcRecv) = iRMin
         BufferS_IP(            iBufferS+3,iProcRecv) = iRMax
         if(nDim > 1)BufferS_IP(iBufferS+4,iProcRecv) = jRMin
         if(nDim > 1)BufferS_IP(iBufferS+5,iProcRecv) = jRMax
         if(nDim > 2)BufferS_IP(iBufferS+6,iProcRecv) = kRMin
         if(nDim > 2)BufferS_IP(iBufferS+7,iProcRecv) = kRMax

         iBufferS = iBufferS + 1 + 2*nDim
         
         do k = kSMin,kSmax; do j = jSMin,jSMax; do i = iSMin,iSmax
            BufferS_IP(iBufferS+1:iBufferS+nVar,iProcRecv) = &
                 State_VGB(:,i,j,k,iBlockSend)
            iBufferS = iBufferS + nVar
         end do; end do; end do

         iBufferS_P(iProcRecv) = iBufferS

         ! Number of reals to be received from the other processor
         ! is the same as sent for equal levels
         nSize = nVar*(iRMax-iRMin+1)*(jRMax-jRMin+1)*(kRMax-kRMin+1)
         iBufferR_P(iProcRecv) = iBufferR_P(iProcRecv) + 1 + 2*nDim + nSize
      end if

    end subroutine do_equal

    !==========================================================================

    subroutine do_restrict

      integer :: iR, jR, kR, iS1, jS1, kS1, iS2, jS2, kS2, iVar
      integer :: iBufferS, nSize
      !------------------------------------------------------------------------

      ! For sideways communication from a fine to a coarser block
      ! the coordinate parity of the sender block tells 
      ! if the receiver block fills into the 
      ! lower (D*Recv = 0) or upper (D*Rev=1) half of the block
      DiRecv = modulo(iTree_IA(Coord1_,iNodeSend)-1, 2)
      DjRecv = modulo(iTree_IA(Coord2_,iNodeSend)-1, 2)
      DkRecv = modulo(iTree_IA(Coord3_,iNodeSend)-1, 2)

      ! Do not restrict diagonally in the direction of the sibling.
      if(iDir == -1 .and. DiRecv==1) RETURN
      if(iDir == +1 .and. DiRecv==0) RETURN
      if(jDir == -1 .and. DjRecv==1 .and. nDimAmr>1) RETURN
      if(jDir == +1 .and. DjRecv==0 .and. nDimAmr>1) RETURN
      if(kDir == -1 .and. DkRecv==1 .and. nDimAmr>2) RETURN
      if(kDir == +1 .and. DkRecv==0 .and. nDimAmr>2) RETURN

      iSend = (3*iDir + 3 + DiRecv)/2
      jSend = (3*jDir + 3 + DjRecv)/2
      kSend = (3*kDir + 3 + DkRecv)/2

      iRecv = iSend - 3*iDir
      jRecv = jSend - 3*jDir
      kRecv = kSend - 3*kDir

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

      !write(*,*)'iDir, jDir, kDir =',iDir, jDir, kDir
      !write(*,*)'iRecv,jRecv,kRecv=',iRecv,jRecv,kRecv
      !
      !write(*,*)'iSMin,iSmax,jSMin,jSMax,kSMin,kSmax=',&
      !     iSMin,iSmax,jSMin,jSMax,kSMin,kSmax
      !
      !write(*,*)'iRMin,iRmax,jRMin,jRMax,kRMin,kRmax=',&
      !     iRMin,iRmax,jRMin,jRMax,kRMin,kRmax

      iNodeRecv  = iNodeNei_IIIB(iSend,jSend,kSend,iBlockSend)
      iProcRecv  = iTree_IA(Proc_,iNodeRecv)
      iBlockRecv = iTree_IA(Block_,iNodeRecv)

      if(iProc == iProcRecv)then
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

         BufferS_IP(            iBufferS+1,iProcRecv) = iBlockRecv
         BufferS_IP(            iBufferS+2,iProcRecv) = iRMin
         BufferS_IP(            iBufferS+3,iProcRecv) = iRMax
         if(nDim > 1)BufferS_IP(iBufferS+4,iProcRecv) = jRMin
         if(nDim > 1)BufferS_IP(iBufferS+5,iProcRecv) = jRMax
         if(nDim > 2)BufferS_IP(iBufferS+6,iProcRecv) = kRMin
         if(nDim > 2)BufferS_IP(iBufferS+7,iProcRecv) = kRMax

         iBufferS = iBufferS + 1 + 2*nDim

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


         ! This processor will receive a prolonged buffer from
         ! the other processor and the "recv" direction of the prolongations
         ! will be the same as the "send" direction for this restriction:
         ! iSend,kSend,jSend = 0..3
         iRMin = iProlongR_DII(1,iSend,Min_)
         iRMax = iProlongR_DII(1,iSend,Max_)
         jRMin = iProlongR_DII(2,jSend,Min_)
         jRMax = iProlongR_DII(2,jSend,Max_)
         kRMin = iProlongR_DII(3,kSend,Min_)
         kRMax = iProlongR_DII(3,kSend,Max_)

         nSize = nVar*(iRMax-iRMin+1)*(jRMax-jRMin+1)*(kRMax-kRMin+1)
         iBufferR_P(iProcRecv) = iBufferR_P(iProcRecv) + 1 + 2*nDim + nSize
      end if

    end subroutine do_restrict

    !==========================================================================

    subroutine do_prolong

      integer :: iR, jR, kR, iS, jS, kS
      integer :: iBufferS, nSize
      !------------------------------------------------------------------------

      ! Loop through the subfaces or subedges
      do DkRecv = (1-kDir)/2, 1-(1+kDir)/2, 3-kRatio
         kSend = (3*kDir + 3 + DkRecv)/2
         kRecv = kSend - 3*kDir
         do DjRecv = (1-jDir)/2, 1-(1+jDir)/2, 3-jRatio
            jSend = (3*jDir + 3 + DjRecv)/2
            jRecv = jSend - 3*jDir
            do DiRecv = (1-iDir)/2, 1-(1+iDir)/2
               iSend = (3*iDir + 3 + DiRecv)/2
               iRecv = iSend - 3*iDir

               ! Sending range depends on iSend,jSend,kSend = 0..3
               iSMin = iProlongS_DII(1,iSend,Min_)
               iSMax = iProlongS_DII(1,iSend,Max_)
               jSMin = iProlongS_DII(2,jSend,Min_)
               jSMax = iProlongS_DII(2,jSend,Max_)
               kSMin = iProlongS_DII(3,kSend,Min_)
               kSMax = iProlongS_DII(3,kSend,Max_)

               ! Receiving range depends on iRecv,kRecv,jRecv = 0..3
               iRMin = iProlongR_DII(1,iRecv,Min_)
               iRMax = iProlongR_DII(1,iRecv,Max_)
               jRMin = iProlongR_DII(2,jRecv,Min_)
               jRMax = iProlongR_DII(2,jRecv,Max_)
               kRMin = iProlongR_DII(3,kRecv,Min_)
               kRMax = iProlongR_DII(3,kRecv,Max_)

               iNodeRecv  = iNodeNei_IIIB(iSend,jSend,kSend,iBlockSend)
               iProcRecv  = iTree_IA(Proc_,iNodeRecv)
               iBlockRecv = iTree_IA(Block_,iNodeRecv)

               !write(*,*)'iRecv, jRecv, kRecv=',iRecv, jRecv, kRecv
               !
               !write(*,*)'iSMin,iSmax,jSMin,jSMax,kSMin,kSmax=',&
               !     iSMin,iSmax,jSMin,jSMax,kSMin,kSmax
               !
               !write(*,*)'iRMin,iRmax,jRMin,jRMax,kRMin,kRmax=',&
               !     iRMin,iRmax,jRMin,jRMax,kRMin,kRmax

               if(iProc == iProcRecv)then
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

                  BufferS_IP(            iBufferS+1,iProcRecv) = iBlockRecv
                  BufferS_IP(            iBufferS+2,iProcRecv) = iRMin
                  BufferS_IP(            iBufferS+3,iProcRecv) = iRMax
                  if(nDim > 1)BufferS_IP(iBufferS+4,iProcRecv) = jRMin
                  if(nDim > 1)BufferS_IP(iBufferS+5,iProcRecv) = jRMax
                  if(nDim > 2)BufferS_IP(iBufferS+6,iProcRecv) = kRMin
                  if(nDim > 2)BufferS_IP(iBufferS+7,iProcRecv) = kRMax

                  iBufferS = iBufferS + 1 + 2*nDim

                  do kR=kRMin,kRMax
                     kS = kSMin + (kR-kRMin)/kRatio
                     do jR=jRMin,jRMax
                        jS = jSMin + (jR-jRMin)/jRatio
                        do iR=iRMin,iRMax
                           iS = iSMin + (iR-iRMin)/iRatio
                           BufferS_IP(iBufferS+1:iBufferS+nVar,iProcRecv)= &
                                State_VGB(:,iS,jS,kS,iBlockSend)
                           iBufferS = iBufferS + nVar
                        end do
                     end do
                  end do
                  iBufferS_P(iProcRecv) = iBufferS

                  ! This processor will receive a restricted buffer from
                  ! the other processor and the "recv" direction of the
                  ! restriction will be the same as the "send" direction for
                  ! this restriction: iSend,kSend,jSend = 0..3
                  iRMin = iRestrictR_DII(1,iSend,Min_)
                  iRMax = iRestrictR_DII(1,iSend,Max_)
                  jRMin = iRestrictR_DII(2,jSend,Min_)
                  jRMax = iRestrictR_DII(2,jSend,Max_)
                  kRMin = iRestrictR_DII(3,kSend,Min_)
                  kRMax = iRestrictR_DII(3,kSend,Max_)

                  nSize = nVar*(iRMax-iRMin+1)*(jRMax-jRMin+1)*(kRMax-kRMin+1)
                  iBufferR_P(iProcRecv) = iBufferR_P(iProcRecv) &
                       + 1 + 2*nDim + nSize
               end if
            end do
         end do
      end do

    end subroutine do_prolong

    !==========================================================================

    subroutine set_range

      integer:: nWidthProlongS_D(MaxDim)
      !------------------------------------------------------------------------

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

      ! number of cells sent from fine block. Rounded up.
      nWidthProlongS_D         = 0
      nWidthProlongS_D(1:nDim) = 1 + (nWidth-1)/iRatio_D(1:nDim)

      ! Indexed by iSend/jSend,kSend = 0..3
      iProlongS_DII(:,0,Min_) = 1
      iProlongS_DII(:,0,Max_) = nWidthProlongS_D
      iProlongS_DII(:,1,Min_) = 1
      iProlongS_DII(:,1,Max_) = nIjk_D/iRatio_D
      iProlongS_DII(:,2,Min_) = nIjk_D/iRatio_D + 1
      iProlongS_DII(:,2,Max_) = nIjk_D
      iProlongS_DII(:,3,Min_) = nIjk_D + 1 - nWidthProlongS_D
      iProlongS_DII(:,3,Max_) = nIjk_D

      if(DoSendCorner)then
         ! Face + two edges + corner or edge + one corner sent together 
         ! from fine to coarse block
         iProlongS_DII(:,1,Max_) = iProlongS_DII(:,1,Max_) + nWidthProlongS_D
         iProlongS_DII(:,2,Min_) = iProlongS_DII(:,2,Min_) - nWidthProlongS_D
      end if

      ! Indexed by iRecv/jRecv/kRecv = 0,1,2,3
      iProlongR_DII(:, 0,Min_) = 1 - nWidth
      iProlongR_DII(:, 0,Max_) = 0
      iProlongR_DII(:, 1,Min_) = 1
      iProlongR_DII(:, 1,Max_) = nIjk_D
      iProlongR_DII(:, 2,Min_) = 1
      iProlongR_DII(:, 2,Max_) = nIjk_D
      iProlongR_DII(:, 3,Min_) = nIjk_D + 1
      iProlongR_DII(:, 3,Max_) = nIjk_D + nWidth

      if(DoSendCorner)then
         ! Face + one edge or edge + one corner received together
         ! from fine to coarse block
         iProlongR_DII(1:nDim,1,Max_) = iProlongR_DII(1:nDim,1,Max_) + nWidth
         iProlongR_DII(1:nDim,2,Min_) = iProlongR_DII(1:nDim,2,Min_) - nWidth
      end if

    end subroutine set_range

  end subroutine message_pass_cell

  !============================================================================

  subroutine test_pass_cell

    use BATL_mpi,  ONLY: iProc
    use BATL_size, ONLY: MaxDim, nDim, nDimAmr, &
         MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nI, nJ, nK, nBlock
    use BATL_tree, ONLY: init_tree, set_tree_root, refine_tree_node, &
         distribute_tree, show_tree, clean_tree, &
         Unused_B, iNode_B, DiLevelNei_IIIB
    use BATL_grid, ONLY: init_grid, create_grid, clean_grid, &
         Xyz_DGB, CellSize_DB
    use BATL_geometry, ONLY: init_geometry

    integer, parameter:: MaxBlockTest            = 15
    integer, parameter:: nRootTest_D(MaxDim)     = (/2,2,2/)
    logical, parameter:: IsPeriodicTest_D(MaxDim)= .true.
    real, parameter:: DomainMin_D(MaxDim) = (/ 1.0, 10.0, 100.0 /)
    real, parameter:: DomainMax_D(MaxDim) = (/ 2.0, 20.0, 200.0 /)
    real, parameter:: DomainSize_D(MaxDim) = DomainMax_D - DomainMin_D

    integer, parameter:: nVar = nDim
    real, allocatable:: State_VGB(:,:,:,:,:)

    integer:: nProlongOrder
    integer:: iSendCorner
    logical:: DoSendCorner

    real:: Xyz_D(MaxDim), Tolerance_D(nDim)
    integer:: iNode, iBlock, i, j, k, jMin, jMax, kMin, kMax, iDim

    logical:: DoTestMe
    character(len=*), parameter :: NameSub = 'test_pass_cell'
    !-----------------------------------------------------------------------
    DoTestMe = iProc == 0

    if(DoTestMe) write(*,*) 'Starting ',NameSub

    call init_tree(MaxBlockTest)
    call init_grid( DomainMin_D(1:nDim), DomainMax_D(1:nDim) )
    call init_geometry( IsPeriodicIn_D = IsPeriodicTest_D(1:nDim) )
    call set_tree_root( nRootTest_D(1:nDim))

    call refine_tree_node(1)
    call distribute_tree(.true.)
    call create_grid

    if(DoTestMe) call show_tree(NameSub,.true.)

    allocate(State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlockTest))

    jMin =  1; if(nDim > 1) jMin = MinJ
    kMin =  1; if(nDim > 2) kMin = MinK
    jMax = nJ; if(nDim > 1) jMax = MaxJ
    kMax = nK; if(nDim > 2) kMax = MaxK


    do nProlongOrder = 1,1; do iSendCorner = 1,2

       DoSendCorner = iSendCorner == 2

       if(DoTestMe)write(*,*) 'testing message_pass_cell with',&
            ' nProlongOrder=',nProlongOrder, &
            ' DoSendCorner=',DoSendCorner


       State_VGB = 0.0

       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          State_VGB(:,1:nI,1:nJ,1:nK,iBlock) = &
               Xyz_DGB(1:nDim,1:nI,1:nJ,1:nK,iBlock)
       end do

       call message_pass_cell(nVar, State_VGB, &
            nProlongOrderIn=nProlongOrder, &
            DoSendCornerIn=DoSendCorner)

       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE

          Tolerance_D = 1e-6
          if(nProlongOrder == 1) then
             !iNode = iNode_B(iBlock)
             !write(*,*)'!!! iBlock, iNode=',iBlock,iNode
             !write(*,*)'!!! DiLevelNei_IIIB(:,:,:,iBlock)=',&
             !     DiLevelNei_IIIB(:,:,:,iBlock)
             if(any(DiLevelNei_IIIB(:,:,:,iBlock) == 1)) &
                  Tolerance_D(1:nDimAmr) = Tolerance_D(1:nDimAmr) &
                  + 0.5*CellSize_DB(1:nDimAmr,iBlock)
          end if

          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI

             Xyz_D = Xyz_DGB(:,i,j,k,iBlock)

             ! Check that no info is send in the non-used dimensions,
             ! i.e. for all iDim: nDim+1 < iDim < MaxDim
             if(j<jMin .or. j>jMax .or. k<kMin .or. k>kMax)then
                do iDim = 1, nDim
                   if(abs(State_VGB(iDim,i,j,k,iBlock)) > 1e-6)then
                      write(*,*)'Face should not be set: ', &
                           'iProc,iBlock,i,j,k,iDim,State,Xyz,Tol=', &
                           iProc,iBlock,i,j,k,iDim, &
                           State_VGB(iDim,i,j,k,iBlock), &
                           Xyz_D(iDim),Tolerance_D(iDim)
                   end if
                end do

                CYCLE
             end if

             ! if we do not need to send corners and edges, check that the
             ! State_VGB in these cells is still the unset value
             if(.not.DoSendCorner.and. ( &
                  (i<1.or.i>nI).and.(j<1.or.j>nJ) .or. &
                  (i<1.or.i>nI).and.(k<1.or.k>nK) .or. &
                  (k<1.or.k>nK).and.(j<1.or.j>nJ) ))then

                do iDim = 1, nDim
                   if(abs(State_VGB(iDim,i,j,k,iBlock)) > 1e-6)then
                      write(*,*)'corner/edge should not be set: ', &
                           'iProc,iBlock,i,j,k,iDim,State,Xyz,Tol=', &
                           iProc,iBlock,i,j,k,iDim, &
                           State_VGB(iDim,i,j,k,iBlock), &
                           Xyz_D(iDim),Tolerance_D(iDim)
                   end if
                end do

                CYCLE
             end if


             ! Shift ghost cell coordinate into periodic domain
             Xyz_D = DomainMin_D + modulo(Xyz_D - DomainMin_D, DomainSize_D)

             do iDim = 1, nDim
                if(abs(State_VGB(iDim,i,j,k,iBlock) - Xyz_D(iDim)) &
                     > Tolerance_D(iDim))then
                   write(*,*)'iProc,iBlock,i,j,k,iDim,State,Xyz,Tol=', &
                        iProc,iBlock,i,j,k,iDim, &
                        State_VGB(iDim,i,j,k,iBlock), &
                        Xyz_D(iDim),Tolerance_D(iDim)
                end if
             end do
          end do; end do; end do
       end do

    end do; end do ! test parameters
    deallocate(State_VGB)

    call clean_grid
    call clean_tree

  end subroutine test_pass_cell

end module BATL_pass_cell
