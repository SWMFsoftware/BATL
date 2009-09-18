!^CFG COPYRIGHT UM
module BATL_pass_face

  ! Possible improvements:
  ! (1) Overlapping communication and calculation

  implicit none

  SAVE

  private ! except

  public message_pass_face
  public test_pass_face

contains

  subroutine message_pass_face(nVar, Flux_VXB, Flux_VYB, Flux_VZB, &
       DoResChangeOnlyIn, TimeOld_B, Time_B, DoTestIn)

    use BATL_size, ONLY: MaxBlock, &
         nBlock, nI, nJ, nK, nIjk_D, &
         MaxDim, nDim, iRatio, jRatio, kRatio

    use BATL_mpi, ONLY: iComm, nProc, iProc

    use BATL_tree, ONLY: &
         iNodeNei_IIIB, DiLevelNei_IIIB, Unused_B, iNode_B, &
         iTree_IA, Proc_, Block_, Coord1_, Coord2_, Coord3_

    use ModNumConst, ONLY: i_DD
    use ModMpi

    ! Arguments
    integer, intent(in) :: nVar
    real, intent(inout) :: &
         Flux_VXB(nVar,2,nJ,nK,MaxBlock), &
         Flux_VYB(nVar,2,nI,nK,MaxBlock), &
         Flux_VZB(nVar,2,nI,nJ,MaxBlock)

    ! Optional arguments
    logical, optional, intent(in) :: DoResChangeOnlyIn
    logical, optional, intent(in) :: DoTestIn
    real,    optional, intent(in) :: TimeOld_B(MaxBlock)
    real,    optional, intent(in) :: Time_B(MaxBlock)

    ! Send sum of fine fluxes to coarse neighbors and subtract the coarse flux
    ! (if any).
    !
    ! DoResChangeOnlyIn determines if the flux correction is applied at
    !     resolution changes only. True is the default.
    !
    ! TimeOld_B and Time_B are the simulation times associated with the
    !    fluxes on the two sides for local time stepping scheme. 
    !    To be implemented !!!

    ! Local variables

    logical :: DoReschangeOnly

    integer :: iDim, iDimSide, iRecvSide, iSign
    integer :: iSend, jSend, kSend, iSide, jSide, kSide
    integer :: iDir, jDir, kDir
    integer :: iNodeRecv, iNodeSend
    integer :: iBlockRecv, iProcRecv, iBlockSend, iProcSend, DiLevel

    ! Fast lookup tables for index ranges per dimension
    integer, parameter:: Min_=1, Max_=2
    integer:: iReceive_DII(MaxDim,0:2,Min_:Max_)

    ! Index range for recv and send segments of the blocks
    integer :: iRMin, iRMax, jRMin, jRMax, kRMin, kRMax
    integer :: iSMin, iSMax, jSMin, jSMax, kSMin, kSMax

    ! Variables related to recv and send buffers
    integer, parameter:: nFaceCell = max(nI*nJ,nI*nK,nJ*nK)

    integer, allocatable, save:: iBufferR_P(:), iBufferS_P(:)

    integer :: MaxBufferS = -1, MaxBufferR = -1, DnBuffer
    real, pointer, save:: BufferR_IP(:,:), BufferS_IP(:,:)

    integer:: iRequestR, iRequestS, iError
    integer, allocatable, save:: iRequestR_I(:), iRequestS_I(:), &
         iStatus_II(:,:)

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'BATL_pass_face::message_pass_face'
    !--------------------------------------------------------------------------
    DoTest = .false.; if(present(DoTestIn)) DoTest = DoTestIn
    if(DoTest)write(*,*)NameSub,' starting with nVar=',nVar

    ! Set values or defaults for optional arguments
    DoResChangeOnly = .true.
    if(present(DoResChangeOnlyIn)) DoResChangeOnly = DoResChangeOnlyIn

    ! Set index ranges based on arguments
    call set_range

    if(nProc > 1)then
       ! Small arrays are allocated once 
       if(.not.allocated(iBufferR_P))then
          allocate(iBufferR_P(0:nProc-1), iBufferS_P(0:nProc-1))
          allocate(iRequestR_I(nProc-1), iRequestS_I(nProc-1))
          allocate(iStatus_II(MPI_STATUS_SIZE,nProc-1))
       end if
       ! Upper estimate of the number of variables sent 
       ! from one block to another. Used for dynamic pointer buffers.
       DnBuffer = nVar*nFaceCell

    end if

    if(nProc>1)then
       ! initialize buffer indexes
       iBufferR_P = 0
       iBufferS_P = 0
    end if

    ! Loop through all nodes
    do iBlockSend = 1, nBlock

       if(Unused_B(iBlockSend)) CYCLE

       iNodeSend = iNode_B(iBlockSend)

       do iDim = 1, nDim
          do iDimSide = 1, 2
             ! Opposite side will receive the fluxes
             iRecvSide = 3 - iDimSide

             ! Set direction indexes
             iSign = 2*iDimSide - 3
             iDir = iSign*i_DD(iDim,1)
             jDir = iSign*i_DD(iDim,2)
             kDir = iSign*i_DD(iDim,3)

             ! Check for resolution change
             DiLevel = DiLevelNei_IIIB(iDir,jDir,kDir,iBlockSend)

             ! For res. change send flux from fine side to coarse side
             if(DoResChangeOnly .and. DiLevel == 0) CYCLE

             if(DiLevel == 0)then
                ! call do_equal !!!
             elseif(DiLevel == 1)then
                call do_restrict
             elseif(DiLevel == -1)then
                call do_prolong
             endif
          end do ! iDimSide
       end do ! iDim
    end do ! iBlockSend

    ! Done for serial run
    if(nProc == 1) RETURN

    call timing_start('send_recv')

    ! Make sure the receive buffer is large enough
    if(maxval(iBufferR_P) > MaxBufferR) call extend_buffer(&
         .false., MaxBufferR, 2*maxval(iBufferR_P), BufferR_IP)

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
    call buffer_to_flux
    call timing_stop('buffer_to_state')

  contains

    subroutine extend_buffer(DoCopy, MaxBufferOld, MaxBufferNew, Buffer_IP)

      logical, intent(in)   :: DoCopy
      integer, intent(in)   :: MaxBufferNew
      integer, intent(inout):: MaxBufferOld
      real, pointer:: Buffer_IP(:,:)

      real, pointer :: OldBuffer_IP(:,:)
      !------------------------------------------------------------------------
      !write(*,*)'extend_buffer called with ',&
      !     'DoCopy, MaxBufferOld, MaxBufferNew=', &
           !     DoCopy, MaxBufferOld, MaxBufferNew

      if(MaxBufferOld < 0 .or. .not.DoCopy)then
         if(MaxBufferOld > 0) deallocate(Buffer_IP)
         allocate(Buffer_IP(MaxBufferNew,0:nProc-1))
      else
         ! store old values
         OldBuffer_IP => Buffer_IP
         ! allocate extended buffer
         allocate(Buffer_IP(MaxBufferNew,0:nProc-1))
         ! copy old values
         Buffer_IP(1:MaxBufferOld,:) = OldBuffer_IP(1:MaxBufferOld,:)  
         ! free old storage
         deallocate(OldBuffer_IP)
      end if
      ! Set new buffer size
      MaxBufferOld = MaxBufferNew

    end subroutine extend_buffer

    !==========================================================================
    subroutine buffer_to_flux

      ! Copy buffer into recv block of State_VGB

      integer:: iBufferR, iTag, iDim, iDimSide, iSubFace1, iSubFace2, i, j, k
      real :: TimeSend, WeightOld, WeightNew
      !------------------------------------------------------------------------

      do iProcSend = 0, nProc - 1
         if(iBufferR_P(iProcSend) == 0) CYCLE

         iBufferR = 0
         do
            ! Read the tag from the buffer
            iTag = nint(BufferR_IP(iBufferR+1,iProcSend))
            iBufferR = iBufferR + 1

            ! Decode iTag = 100*iBlockRecv + 20*iDim + 10*(iDimSide-1) 
            !               + 3*iSubFace1 + iSubFace2
            iBlockRecv = iTag/100; iTag = iTag - 100*iBlockRecv
            iDim       = iTag/20;  iTag = iTag - 20*iDim
            iDimSide   = iTag/10;  iTag = iTag - 10*iDimSide
            iDimSide   = iDimSide + 1
            iSubFace1  = iTag/3;   iTag = iTag - 3*iSubFace1
            iSubFace2  = iTag

            if(present(Time_B))then
               ! Get time of neighbor and interpolate/extrapolate ghost cells
               iBufferR = iBufferR + 1
               TimeSend  = BufferR_IP(iBufferR,iProcSend)
               WeightOld = (TimeSend - Time_B(iBlockRecv)) &
                    / max(TimeSend - TimeOld_B(iBlockRecv), 1e-30)
               WeightNew = 1 - WeightOld
            end if

            select case(iDim)
            case(1)
               ! Get transverse index ranges for the (sub)face
               jRMin = iReceive_DII(2,iSubFace1,Min_)
               jRMax = iReceive_DII(2,iSubFace1,Max_)
               kRMin = iReceive_DII(3,iSubFace2,Min_)
               kRMax = iReceive_DII(3,iSubFace2,Max_)

               ! Add sent flux to the original
               do k = kRMin, kRmax; do j = jRMin, jRMax
                  Flux_VXB(:,iDimSide,j,k,iBlockRecv) = &
                       Flux_VXB(:,iDimSide,j,k,iBlockRecv) &
                       + BufferR_IP(iBufferR+1:iBufferR+nVar,iProcSend)
                  iBufferR = iBufferR + nVar
               end do; end do

            case(2)
               iRMin = iReceive_DII(1,iSubFace1,Min_)
               iRMax = iReceive_DII(1,iSubFace1,Max_)
               kRMin = iReceive_DII(3,iSubFace2,Min_)
               kRMax = iReceive_DII(3,iSubFace2,Max_)

               do k = kRMin, kRmax; do i = iRMin, iRMax
                  Flux_VYB(:,iDimSide,i,k,iBlockRecv) = &
                       Flux_VYB(:,iDimSide,i,k,iBlockRecv) &
                       + BufferR_IP(iBufferR+1:iBufferR+nVar,iProcSend)
                  iBufferR = iBufferR + nVar
               end do; end do

            case(3)
               iRMin = iReceive_DII(1,iSubFace1,Min_)
               iRMax = iReceive_DII(1,iSubFace1,Max_)
               jRMin = iReceive_DII(2,iSubFace2,Min_)
               jRMax = iReceive_DII(2,iSubFace2,Max_)

               do j = jRMin, jRmax; do i = iRMin, iRMax
                  Flux_VZB(:,iDimSide,i,j,iBlockRecv) = &
                       Flux_VZB(:,iDimSide,i,j,iBlockRecv) &
                       + BufferR_IP(iBufferR+1:iBufferR+nVar,iProcSend)
                  iBufferR = iBufferR + nVar
               end do; end do
            end select

            if(iBufferR >= iBufferR_P(iProcSend)) EXIT

         end do
      end do

    end subroutine buffer_to_flux

    !==========================================================================

    subroutine do_restrict

      !------------------------------------------------------------------------

      ! The coordinate parity of the sender block tells 
      ! if the receiver block fills into the 
      ! lower or upper part of the face

      iSide = 0; if(iRatio==2) iSide = modulo(iTree_IA(Coord1_,iNodeSend)-1, 2)
      jSide = 0; if(jRatio==2) jSide = modulo(iTree_IA(Coord2_,iNodeSend)-1, 2)
      kSide = 0; if(kRatio==2) kSide = modulo(iTree_IA(Coord3_,iNodeSend)-1, 2)

      iSend = (3*iDir + 3 + iSide)/2
      jSend = (3*jDir + 3 + jSide)/2
      kSend = (3*kDir + 3 + kSide)/2

      iNodeRecv  = iNodeNei_IIIB(iSend,jSend,kSend,iBlockSend)
      iProcRecv  = iTree_IA(Proc_,iNodeRecv)
      iBlockRecv = iTree_IA(Block_,iNodeRecv)

      select case(iDim)
      case(1)
         call do_flux(2, 3, nJ, nK, jRatio, kRatio, jSide, kSide, Flux_VXB)
      case(2)
         call do_flux(1, 3, nI, nK, iRatio, kRatio, iSide, kSide, Flux_VYB)
      case(3)
         call do_flux(1, 2, nI, nJ, iRatio, jRatio, iSide, jSide, Flux_VZB)
      end select

    end subroutine do_restrict
    !=======================================================================
    subroutine do_flux(iDim1, iDim2, n1, n2, Dn1, Dn2, iSide1, iSide2, &
         Flux_VFB)

      integer, intent(in):: iDim1, iDim2, n1, n2, Dn1, Dn2, iSide1, iSide2
      real, intent(inout):: Flux_VFB(nVar,2,n1,n2,MaxBlock)

      integer:: iSubFace1, iSubFace2
      integer:: iR1, iR1Min, iR1Max, iR2, iR2Min, iR2Max
      integer:: iS1Min, iS1Max, iS2Min, iS2Max
      integer:: iVar, iBufferS
      !---------------------------------------------------------------------

      ! For Dn1=1 set iSubFace1=0 (full face),
      ! for Dn1=2 set iSubFace1=1 (2) for lower (upper) half
      iSubFace1 = (1 + iSide1)*(Dn1 - 1)
      iSubFace2 = (1 + iSide2)*(Dn2 - 1)

      ! Receiving range depends on subface indexes
      iR1Min = iReceive_DII(iDim1,iSubFace1,Min_)
      iR1Max = iReceive_DII(iDim1,iSubFace1,Max_)
      iR2Min = iReceive_DII(iDim2,iSubFace2,Min_)
      iR2Max = iReceive_DII(iDim2,iSubFace2,Max_)

      if(iProc == iProcRecv)then

         !if(present(Time_B))then
         !   ! Get time of neighbor and interpolate/extrapolate ghost cells
         !   WeightOld = (Time_B(iBlockSend) - Time_B(iBlockRecv)) &
              !        /   max(Time_B(iBlockSend) - TimeOld_B(iBlockRecv), 1e-30)
         !   WeightNew = 1 - WeightOld
         !else
         !   WeightNew = 1.0
         !   WeightOld = 0.0
         !end if

         do iR2 = iR2Min, iR2Max
            iS2Min = 1 + Dn2*(iR2-iR2Min)
            iS2Max = iS2Min + Dn2 - 1
            do iR1 = iR1Min, iR1Max
               iS1Min = 1 + Dn1*(iR1-iR1Min)
               iS1Max = iS1Min + Dn1 - 1
               do iVar = 1, nVar
                  Flux_VFB(iVar,iRecvSide,iR1,iR2,iBlockRecv) = &
                       Flux_VFB(iVar,iRecvSide,iR1,iR2,iBlockRecv) + &
                       sum(Flux_VFB(iVar,iDimSide,iS1Min:iS1Max,iS2Min:iS2Max, &
                       iBlockSend))
               end do
            end do
         end do
      else
         iBufferS = iBufferS_P(iProcRecv)

         if(iBufferS + DnBuffer > MaxBufferS) call extend_buffer( &
              .true., MaxBufferS, 2*(iBufferS+DnBuffer), BufferS_IP)

         ! Encode all necessary info into a single "tag"
         BufferS_IP(iBufferS+1,iProcRecv) = &
              100*iBlockRecv + 20*iDim + 10*(iDimSide-1) &
              + 3*iSubFace1 + iSubFace2

         iBufferS = iBufferS + 1

         !if(present(Time_B))then
         !   iBufferS = iBufferS + 1
         !   BufferS_IP(iBufferS,iProcRecv) = Time_B(iBlockSend)
         !end if

         do iR2 = iR2Min, iR2Max
            iS2Min = 1 + Dn2*(iR2-iR2Min)
            iS2Max = iS2Min + Dn2 - 1
            do iR1 = iR1Min, iR1Max
               iS1Min = 1 + Dn1*(iR1-iR1Min)
               iS1Max = iS1Min + Dn1 - 1
               do iVar = 1, nVar
                  BufferS_IP(iBufferS+iVar,iProcRecv) = &
                       sum(Flux_VFB(iVar,iDimSide,iS1Min:iS1Max,iS2Min:iS2Max, &
                       iBlockSend))
               end do
               iBufferS = iBufferS + nVar
            end do
         end do
         iBufferS_P(iProcRecv) = iBufferS

      end if

    end subroutine do_flux

    !==========================================================================

    subroutine do_prolong

      integer :: nSize
      !------------------------------------------------------------------------

      ! Loop through the subfaces
      do kSide = (1-kDir)/2, 1-(1+kDir)/2, 3-kRatio
         kSend = (3*kDir + 3 + kSide)/2
         do jSide = (1-jDir)/2, 1-(1+jDir)/2, 3-jRatio
            jSend = (3*jDir + 3 + jSide)/2
            do iSide = (1-iDir)/2, 1-(1+iDir)/2
               iSend = (3*iDir + 3 + iSide)/2

               iNodeRecv  = iNodeNei_IIIB(iSend,jSend,kSend,iBlockSend)
               iProcRecv  = iTree_IA(Proc_,iNodeRecv)

               if(iProc == iProcRecv) CYCLE

               select case(iDim)
               case(1)
                  nSize = nVar*nJ*nK/(jRatio*kRatio)
               case(2)
                  nSize = nVar*nI*nK/(iRatio*kRatio)
               case(3)
                  nSize = nVar*nI*nJ/(iRatio*jRatio)
               end select
               if(present(Time_B)) nSize = nSize + 1
               iBufferR_P(iProcRecv) = iBufferR_P(iProcRecv) + 1 + nSize

            end do
         end do
      end do

    end subroutine do_prolong

    !==========================================================================

    subroutine set_range

      ! Full face
      iReceive_DII(:,0,Min_) = 1
      iReceive_DII(:,0,Max_) = nIjk_D

      ! Lower subface
      iReceive_DII(:,1,Min_) = 1
      iReceive_DII(:,1,Max_) = nIjk_D/2

      ! Upper subface
      iReceive_DII(:,2,Min_) = nIjk_D/2 + 1
      iReceive_DII(:,2,Max_) = nIjk_D

    end subroutine set_range

  end subroutine message_pass_face

  !============================================================================

  subroutine test_pass_face

!    use BATL_mpi,  ONLY: iProc
!    use BATL_size, ONLY: MaxDim, nDim, iRatio, jRatio, kRatio, &
!         MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nG, nI, nJ, nK, nBlock
!    use BATL_tree, ONLY: init_tree, set_tree_root, find_tree_node, &
!         refine_tree_node, distribute_tree, show_tree, clean_tree, &
!         Unused_B, DiLevelNei_IIIB
!    use BATL_grid, ONLY: init_grid, create_grid, clean_grid, &
!         Xyz_DGB, CellSize_DB
!    use BATL_geometry, ONLY: init_geometry
!
!    integer, parameter:: MaxBlockTest            = 50
!    integer, parameter:: nRootTest_D(MaxDim)     = (/3,3,3/)
!    logical, parameter:: IsPeriodicTest_D(MaxDim)= .true.
!    real, parameter:: DomainMin_D(MaxDim) = (/ 1.0, 10.0, 100.0 /)
!    real, parameter:: DomainMax_D(MaxDim) = (/ 4.0, 40.0, 400.0 /)
!    real, parameter:: DomainSize_D(MaxDim) = DomainMax_D - DomainMin_D
!
!    real, parameter:: Tolerance = 1e-6
!
!    integer, parameter:: nVar = nDim
!    real, allocatable:: State_VGB(:,:,:,:,:)
!
!    integer:: nWidth
!    integer:: nProlongOrder
!    integer:: nCoarseLayer
!    integer:: iSendCorner,  iRestrictFace
!    logical:: DoSendCorner, DoRestrictFace
!
!    real:: Xyz_D(MaxDim)
!    integer:: iNode, iBlock, i, j, k, iMin, iMax, jMin, jMax, kMin, kMax, iDim
!    integer:: iDir, jDir, kDir, Di, Dj, Dk
!
!    logical:: DoTestMe
!    character(len=*), parameter :: NameSub = 'test_pass_face'
!    !-----------------------------------------------------------------------
!    DoTestMe = iProc == 0
!
!    if(DoTestMe) write(*,*) 'Starting ',NameSub
!
!    call init_tree(MaxBlockTest)
!    call init_grid( DomainMin_D(1:nDim), DomainMax_D(1:nDim) )
!    call init_geometry( IsPeriodicIn_D = IsPeriodicTest_D(1:nDim) )
!    call set_tree_root( nRootTest_D(1:nDim))
!
!    call find_tree_node( (/0.5,0.5,0.5/), iNode)
!    if(DoTestMe)write(*,*) NameSub,' middle node=',iNode
!    call refine_tree_node(iNode)
!    call distribute_tree(.true.)
!    call create_grid
!
!    if(DoTestMe) call show_tree(NameSub,.true.)
!
!    allocate(State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlockTest))
!
!    do nProlongOrder = 1, 2; do nCoarseLayer = 1, 2; do nWidth = 1, nG
!
!       ! Second order prolongation does not work with sending multiple coarse 
!       ! cell layers into the fine cells with their original values. 
!       if(nProlongOrder == 2 .and. nCoarseLayer == 2) CYCLE
!
!       ! Cannot send more coarse layers than the number of ghost cell layers
!       if(nCoarseLayer > nWidth) CYCLE
!
!       if(DoTestMe)write(*,*) 'testing message_pass_face with', &
!            ' nProlongOrder=',  nProlongOrder, &
!            ' nCoarseLayer=',   nCoarseLayer,  &
!            ' nWidth=',         nWidth
!
!       ! Set the range of ghost cells that should be set
!       iMin =  1 - nWidth
!       jMin =  1; if(nDim > 1) jMin = 1 - nWidth
!       kMin =  1; if(nDim > 2) kMin = 1 - nWidth
!       iMax = nI + nWidth
!       jMax = nJ; if(nDim > 1) jMax = nJ + nWidth
!       kMax = nK; if(nDim > 2) kMax = nK + nWidth
!
!       do iSendCorner = 1, 2; do iRestrictFace = 1, 2
!
!          DoSendCorner   = iSendCorner   == 2
!          DoRestrictFace = iRestrictFace == 2
!
!          ! Second order prolongation does not work with restricting face:
!          ! the first order restricted cell cannot be used in the prolongation.
!          if(DoRestrictFace .and. nProlongOrder == 2) CYCLE
!
!          ! Sending multiple coarse cell layers is not meaningful for edges
!          ! and corners. The current algorithm does not support this !!!
!          if(nCoarseLayer > 1 .and. DoSendCorner) CYCLE
!          
!          if(DoTestMe)write(*,*) 'testing message_pass_face with', &
!               ' DoSendCorner=',   DoSendCorner, &
!               ' DoRestrictFace=', DoRestrictFace
!
!          State_VGB = 0.0
!
!          do iBlock = 1, nBlock
!             if(Unused_B(iBlock)) CYCLE
!             State_VGB(:,1:nI,1:nJ,1:nK,iBlock) = &
!                  Xyz_DGB(1:nDim,1:nI,1:nJ,1:nK,iBlock)
!          end do
!
!          call message_pass_face(nVar, State_VGB, &
!               nProlongOrderIn =nProlongOrder,    &
!               nCoarseLayerIn  =nCoarseLayer,     &
!               nWidthIn        =nWidth,           &
!               DoSendCornerIn  =DoSendCorner,     &
!               DoRestrictFaceIn=DoRestrictFace)
!
!          do iBlock = 1, nBlock
!             if(Unused_B(iBlock)) CYCLE
!
!             ! Loop through all cells including ghost cells
!             do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
!
!                ! The filled in second order accurate ghost cell value 
!                ! should be the same as the coordinates of the cell center
!                Xyz_D = Xyz_DGB(:,i,j,k,iBlock)
!
!                ! Check that no info is sent in the non-used dimensions,
!                ! i.e. for all iDim: nDim+1 < iDim < MaxDim
!                if(  i < iMin .or. i > iMax .or. &
!                     j < jMin .or. j > jMax .or. &
!                     k < kMin .or. k > kMax) then
!                   do iDim = 1, nDim
!                      if(abs(State_VGB(iDim,i,j,k,iBlock)) > 1e-6)then
!                         write(*,*)'Face should not be set: ', &
!                              'iProc,iBlock,i,j,k,iDim,State,Xyz=', &
!                              iProc,iBlock,i,j,k,iDim, &
!                              State_VGB(iDim,i,j,k,iBlock), &
!                              Xyz_D(iDim)
!                      end if
!                   end do
!
!                   CYCLE
!                end if
!
!                ! Get the direction vector
!                iDir = 0; if(i<1) iDir = -1; if(i>nI) iDir = 1
!                jDir = 0; if(j<1) jDir = -1; if(j>nJ) jDir = 1
!                kDir = 0; if(k<1) kDir = -1; if(k>nK) kDir = 1
!
!                ! if we do not send corners and edges, check that the
!                ! State_VGB in these cells is still the unset value
!                if(.not.DoSendCorner .and. ( &
!                     iDir /= 0 .and. jDir /= 0 .or. &
!                     iDir /= 0 .and. kDir /= 0 .or. &
!                     jDir /= 0 .and. kDir /= 0 ))then
!
!                   do iDim = 1, nDim
!                      if(abs(State_VGB(iDim,i,j,k,iBlock)) > 1e-6)then
!                         write(*,*)'corner/edge should not be set: ', &
!                              'iProc,iBlock,i,j,k,iDim,State,Xyz=', &
!                              iProc,iBlock,i,j,k,iDim, &
!                              State_VGB(iDim,i,j,k,iBlock), &
!                              Xyz_D
!                      end if
!                   end do
!
!                   CYCLE
!                end if
!
!                ! Shift ghost cell coordinate into periodic domain
!                Xyz_D = DomainMin_D + modulo(Xyz_D - DomainMin_D, DomainSize_D)
!
!                ! Calculate distance of ghost cell layer
!                Di = 0; Dj = 0; Dk = 0
!                if(i <  1 .and. iRatio == 2) Di = 2*i-1
!                if(i > nI .and. iRatio == 2) Di = 2*(i-nI)-1
!                if(j <  1 .and. jRatio == 2) Dj = 2*j-1
!                if(j > nJ .and. jRatio == 2) Dj = 2*(j-nJ)-1
!                if(k <  1 .and. kRatio == 2) Dk = 2*k-1
!                if(k > nK .and. kRatio == 2) Dk = 2*(k-nK)-1
!
!                if(DoRestrictFace .and. &
!                     DiLevelNei_IIIB(iDir,jDir,kDir,iBlock) == -1)then
!
!                   ! Shift coordinates if only 1 layer of fine cells
!                   ! is averaged in the orthogonal direction
!                   Xyz_D(1) = Xyz_D(1) - 0.25*Di*CellSize_DB(1,iBlock)
!                   Xyz_D(2) = Xyz_D(2) - 0.25*Dj*CellSize_DB(2,iBlock)
!                   Xyz_D(3) = Xyz_D(3) - 0.25*Dk*CellSize_DB(3,iBlock)
!
!                end if
!
!                if(nProlongOrder == 1 .and. &
!                     DiLevelNei_IIIB(iDir,jDir,kDir,iBlock) == 1)then
!
!                   ! Shift depends on the parity of the fine ghost cell
!                   ! except when there is no AMR or multiple coarse cell 
!                   ! layers are sent in that direction
!                   if(iRatio == 2 .and. (nCoarseLayer == 1 .or. iDir == 0)) &
!                        Di = 2*modulo(i,2) - 1
!                   if(jRatio == 2 .and. (nCoarseLayer == 1 .or. jDir == 0)) &
!                        Dj = 2*modulo(j,2) - 1
!                   if(kRatio == 2 .and. (nCoarseLayer == 1 .or. kDir == 0)) &
!                        Dk = 2*modulo(k,2) - 1
!
!                   Xyz_D(1) = Xyz_D(1) + 0.5*Di*CellSize_DB(1,iBlock)
!                   Xyz_D(2) = Xyz_D(2) + 0.5*Dj*CellSize_DB(2,iBlock)
!                   Xyz_D(3) = Xyz_D(3) + 0.5*Dk*CellSize_DB(3,iBlock)
!
!                end if
!
!                do iDim = 1, nDim
!                   if(abs(State_VGB(iDim,i,j,k,iBlock) - Xyz_D(iDim)) &
!                        > Tolerance)then
!                      write(*,*)'iProc,iBlock,i,j,k,iDim,State,Xyz=', &
!                           iProc,iBlock,i,j,k,iDim, &
!                           State_VGB(iDim,i,j,k,iBlock), &
!                           Xyz_D(iDim)
!                   end if
!                end do
!             end do; end do; end do
!          end do
!
!       end do; end do; end do
!    end do; end do ! test parameters
!    deallocate(State_VGB)
!
!    call clean_grid
!    call clean_tree

  end subroutine test_pass_face

end module BATL_pass_face
