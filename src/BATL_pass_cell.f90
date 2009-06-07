!^CFG COPYRIGHT UM

module BATL_pass_cell

  use BATL_size, ONLY: MaxBlock, nBlock, &
       nI, nJ, nK, nG, MinI, MaxI, MinJ, MaxJ,MinK, MaxK, nDim

  use BATL_mpi, ONLY: iComm, nProc, iProc, barrier_mpi

  use BATL_tree, ONLY: Unused_B, iNodeNei_IIIB, DiLevelNei_IIIB, &
       iTree_IA, Proc_, Block_, Status_, Unset_, Unused_

  use ModMpi

  implicit none

  SAVE

  private ! except

  public message_pass_cell

contains

  subroutine message_pass_cell(nVar, State_VGB, &
       nWidthIn, nProlongOrderIn, nCoarseLayerIn, DoSendCornerIn, DoRestrictFaceIn)

    ! Set ghost cells in State_VGB. 
    ! Number of variables is nVar.
    ! The thickness of ghost cell layers is nWidthIn. 
    ! If sendcorners=.true., we send ghostcells already updated,
    ! in effect, the edge and corner information gets propagated
    ! For prolongorder=1 the prolongation operator for coarse to fine messages
    ! is only first order accurate. For prolongorder=2 there are various 2nd
    ! order prolongation operators implemented depending on the global variable
    !
    !   prolong_type='central'  (gradient based on central differencing)
    !                'central2' (fully 2nd order using one sided if needed)
    !                'lr'       (left and right gradients are calculated)
    !                'lr2'      (fully 2nd order using opposite side if needed) 
    !                'lr3'      (2nd order assuming that equal & restr. are done
    !                            for all dir., so sendcorners must be false.)
    !                'minmod'   (minmod limited gradient is calculated)
    !                'mindod2'  (fully 2nd order using unlimited if needed)
    !
    ! If restrictface=.true., average only the 4 fine cells next to the face
    ! If DoTwoCoarseLayers=.true., prolong the 2nd coarse cell into the 
    !                              2nd layer of fine cells for TVD res change.
    !                              Only for State_VGB and prolongorder=1.

    ! Notation: _o original   (for equal blocks) 
    !           _g ghost      (for equal blocks)
    !           R  originals to be restricted
    !           _r restricted (to be sent to a coarser block)
    !           P  originals to be prolonged
    !           _p prolonged  (to be sent to a refined block)
    !           _s subface    (one quarter of a face)


    character(len=10), parameter:: prolong_type = 'lr'


    ! Arguments
    integer, intent(in) :: nVar
    real, intent(inout) :: State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

    ! Optional arguments
    integer, optional, intent(in) :: nWidthIn
    integer, optional, intent(in) :: nCoarseLayerIn
    integer, optional, intent(in) :: nProlongOrderIn
    logical, optional, intent(in) :: DoSendCornerIn
    logical, optional, intent(in) :: DoRestrictFaceIn

    ! Local variables

    logical :: DoSendCorner
    integer :: nWidth
    integer :: nProlongOrder
    logical :: DoRestrictFace

    integer :: iError

    integer :: iDim, iFaceSend, iFaceRecv, iside, isubface

    ! Components of direction vector: values are -1,0,1
    integer :: iDir, jDir, kDir

    ! Indexes for neighbor block: values are 0, 1, 2, 3
    ! where 0 is to the left, 3 is to the right and
    ! 1 and 2 are for refined neighbors. 1 is used for
    ! same level and coarser
    integer :: iSend, jSend, kSend, iRecv, jRecv, kRecv

    integer :: ibuf, nWidth1, nWidth2, iBlock

    integer :: imin_g,imax_g,jmin_g,jmax_g,kmin_g,kmax_g
    integer :: imin_o,imax_o,jmin_o,jmax_o,kmin_o,kmax_o
    integer :: imin_r,imax_r,jmin_r,jmax_r,kmin_r,kmax_r
    integer :: iminR,imaxR,jminR,jmaxR,kminR,kmaxR
    integer :: imin_p,imax_p,jmin_p,jmax_p,kmin_p,kmax_p
    integer :: iminP,imaxP,jminP,jmaxP,kminP,kmaxP
    integer :: imin_s,imax_s,jmin_s,jmax_s,kmin_s,kmax_s

    integer       :: iNodeNei, iProcNei, iBlockNei
    integer       :: neiL,isize,isize_r,ichild,itag,request
    integer       :: number_receive_requests=0
    integer       :: receive_requests(MaxBlock*6)
    integer       :: status(MPI_STATUS_SIZE, MaxBlock*6)  

    ! Buffer to hold incoming layers of ghost cell values.
    ! For all blocks and all faces. Second index is side.
    real, allocatable:: Buffer(:,:)
    integer :: MaxBuffer               ! size of buffer

    ! Equal, restricted and prolonged values are stored in these arrays
    real, dimension(:,:,:,:), allocatable :: eq_buf, re_buf, pr_buf, &
         pr_buf_s, avrg_state

    ! Small allocatable arrays for prolongation
    real, dimension(:,:,:), allocatable :: qsol,avrg,avrg1,gradx,grady,gradz,&
         gradxl,gradxr,gradyl,gradyr,gradzl,gradzr,pr_buf1

    ! Number of coarse layers to be prolonged: 
    ! 1 for 1st order, 2 for TVD res change.
    ! nCoarse1 = nCoarseLayer-1
    integer :: nCoarseLayer, nCoarse1

    character(len=*), parameter:: NameSub = 'BATL_pass_cell::message_pass_cell'

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
    
    if(DoTestMe)write(*,*)'prolong_type=',prolong_type

    if(nVar < 1)&
         call CON_stop('message_pass_dir can only be used with nVar > 0')

    if(nWidth < 1 .or. nWidth > nG)call CON_stop(NameSub // &
         ': nWidthIn must be 1 ... nG')

    if(nCoarseLayer == 2 .and. nProlongOrder /= 1) &
         call CON_stop(NameSub // &
         ': nCoarseLayerIn=2 requires nProlongOrderIn=1')

    ! These will be useful for index limits
    nCoarse1 = nCoarseLayer - 1
    nWidth1 = nWidth - 1; nWidth2 = 2*nWidth - 1

    ! Allocate buffer if necessary

    ! Maximum size of the total ghost cell layers to be sent
    MaxBuffer = 2*nDim*nVar*nWidth*max(nI*nJ,nI*nK,nJ*nK)*MaxBlock
    if(allocated(Buffer))then
       if(size(Buffer) < MaxBuffer) deallocate(Buffer)
    end if
    if(.not.allocated(Buffer)) allocate(Buffer(MaxBuffer/2,2))
       
    if(nProlongOrder == 1)then
       ! Send messages direction by direction
       do iDim = 1, nDim
          call msg_pass_faces(2*iDim-1,2*iDim,.true.,.true.,.true.)
       end do
    else
       ! nProlongOrder is 2, prolonged must be done after equal and coarse
       ! Send messages for all directions but do prolonged later
       do iDim = 1, nDim
          ! Do equal and restricted
          call msg_pass_faces(2*iDim-1,2*iDim,.true. ,.true. ,.false.)
       end do
       do iDim = 1, nDim
          ! Do prolonged
          call msg_pass_faces(2*iDim-1,2*iDim,.false.,.false.,.true. )
       end do
    end if ! nProlongOrder

    if(DoTestMe)write(*,*)NameSub,' finished'

  contains

    subroutine msg_pass_faces(&
         ifacemin,ifacemax,do_equal,do_restricted,do_prolonged)

      integer:: iNodeNei, iProcNei
      integer, intent(in):: ifacemin,ifacemax
      logical:: do_equal,do_restricted,do_prolonged
      !------------------------------------------------------------------------

      ! Debug
      !if(okdebug)buffer=0.00

      number_receive_requests = 0
      receive_requests = MPI_REQUEST_NULL

      do iFaceSend = iFacemin, iFacemax

         ! Set index ranges for the face
         call setranges

         ! Size of ghost cell layer including all the variables
         isize=(imax_g-imin_g+1)*(jmax_g-jmin_g+1)*(kmax_g-kmin_g+1)*nvar

         ! Size of restricted ghost cell layer
         isize_r=isize/4

         do iBlock = 1,nBlock
            if(Unused_B(iBlock)) CYCLE
            ! Post non-blocking receive for opposite face of neighbor block

!!!            select case(neiLEV(iFaceRecv,iBlock))
            neiL = DiLevelNei_IIIB(-iDir,-jDir,-kDir,iBlock)
            select case(neiL)
            case(0)
               if(do_equal)call remote_receive(1,isize)
            case(-1)
               if(do_restricted)call remote_receive(4,isize_r)
            case(1)
               if(do_prolonged)call remote_receive(1,isize)
            case(Unset_)
               ! Do nothing
            case default
               write(*,*)'me,iBlock,iFaceRecv,neiLEV=',&
                    iProc,iBlock,iFaceRecv,neiL
               call CON_stop(&
                    'Error in message pass: Invalid value for neiLEV')
            end select
         end do ! iBlock
      end do ! iface

      !\
      ! Wait for all receive commands to be posted for all processors
      !/
      call barrier_mpi

      !\
      ! Send blocking messages with Rsend (ready to receive)
      !/
      do iFaceSend = iFaceMin,iFaceMax

         ! Set index ranges for the face
         call setranges

         if(do_equal)&
              allocate(eq_buf(nVar,imin_o:imax_o,jmin_o:jmax_o,kmin_o:kmax_o))
         if(do_restricted)&
              allocate(re_buf(nVar,imin_r:imax_r,jmin_r:jmax_r,kmin_r:kmax_r))
         if(do_prolonged)&
              allocate(&
              pr_buf(nVar,imin_p:imax_p,jmin_p:jmax_p,kmin_p:kmax_p),&
              pr_buf_s(nVar,imin_o:imax_o,jmin_o:jmax_o,kmin_o:kmax_o),&
              avrg_state(nVar,iminP:imaxP,jminP:jmaxP,kminP:kmaxP))

         if(do_prolonged.and.nProlongOrder>1)then
            allocate(&
                 qsol(iminP-1:imaxP+1,jminP-1:jmaxP+1,kminP-1:kmaxP+1),&
                 avrg(iminP:imaxP,jminP:jmaxP,kminP:kmaxP),&
                 pr_buf1(imin_p:imax_p,jmin_p:jmax_p,kmin_p:kmax_p))
            if(prolong_type(1:7)=='central'.or.prolong_type(1:6)=='minmod')&
                 allocate(&
                 gradx(iminP:imaxP,jminP:jmaxP,kminP:kmaxP),&
                 grady(iminP:imaxP,jminP:jmaxP,kminP:kmaxP),&
                 gradz(iminP:imaxP,jminP:jmaxP,kminP:kmaxP))
            if(prolong_type(1:2)=='lr'.or.prolong_type(1:6)=='minmod')&
                 allocate(&
                 gradxl(iminP:imaxP,jminP:jmaxP,kminP:kmaxP),&
                 gradyl(iminP:imaxP,jminP:jmaxP,kminP:kmaxP),&
                 gradzl(iminP:imaxP,jminP:jmaxP,kminP:kmaxP),&
                 gradxr(iminP:imaxP,jminP:jmaxP,kminP:kmaxP),&
                 gradyr(iminP:imaxP,jminP:jmaxP,kminP:kmaxP),&
                 gradzr(iminP:imaxP,jminP:jmaxP,kminP:kmaxP))
         end if

         do iBlock = 1, nBlock
            if(Unused_B(iBlock)) CYCLE

            ! neiL=neiLEV(iface,iBlock)
            neiL = DiLevelNei_IIIB(iDir,jDir,kDir,iBlock)
            select case(neiL)
            case(0)
               if(do_equal)call send_equal
            case(1)
               if(do_restricted)call send_restricted
            case(-1)
               if(do_prolonged)call send_prolonged
            case(Unset_)
               ! There is no neighbor, do nothing
            case default
               write(*,*)'me,iBlock,iFaceSend,neiLEV=',&
                    iProc,iBlock,iFaceSend,neiL
               call CON_stop(&
                    'Error in message pass: Invalid value for neiLEV')
            end select
         end do ! iBlock

         if(do_equal)deallocate(eq_buf)
         if(do_restricted)deallocate(re_buf)
         if(do_prolonged)then
            deallocate(pr_buf,pr_buf_s)
            deallocate(avrg_state)
            if(nProlongOrder>1)then
               deallocate(qsol,avrg,pr_buf1)
               if(prolong_type(1:7)=='central'.or.prolong_type(1:6)=='minmod')&
                    deallocate(gradx,grady,gradz)
               if(prolong_type(1:2)=='lr'.or.prolong_type(1:6)=='minmod')&
                    deallocate(gradxl,gradyl,gradzl,gradxr,gradyr,gradzr)
            end if
         end if

         if(DoTestMe)write(*,*)'messages sent, me, iFaceSend=',iProc,iFaceSend
      end do ! iface

      !\
      ! WAIT FOR ALL MESSAGES TO BE RECEIVED
      !/
      if (number_receive_requests > 0) &
           call MPI_waitall(number_receive_requests,receive_requests,status,iError)

      if(DoTestMe)write(*,*)'messages received, me, iDim=',iProc, iDim


      ! Copy ghost cells received from non-local neigbors
      ! and stored in the buffer into State_VGB

      do iFaceSend = ifacemin, ifacemax

         ! Set index ranges for the face
         call setranges

         do iBlock = 1,nBlock
            if(Unused_B(iBlock))CYCLE
            ibuf=(iBlock-1)*isize+1
            neiL = DiLevelNei_IIIB(-iDir,-jDir,-kDir,iBlock)
            select case(neiL)
            case(0)

!               neiP=neiPE(1,iFaceRecv,iBlock)
!               if(do_equal .and. neiP/=iProc .and. &
!                    .not.unusedBlock_BP(neiBlock(1,iFaceRecv,iBlock),neiP))then

               iNodeNei = iNodeNei_IIIB(iRecv,jRecv,kRecv,iBlock)
               iProcNei = iTree_IA(Proc_,iNodeNei)
               if(do_equal .and. iProcNei /= iProc .and. &
                    iTree_IA(Status_,iNodeNei) /= Unused_)then

                  call buffer2face_state(buffer(ibuf,iside),&
                       imin_g,imax_g,jmin_g,jmax_g,kmin_g,kmax_g,nvar)
               end if
            case(-1)
               if(do_restricted)then
                  call buffer2subfaces_state(buffer(ibuf,iside),&
                       imin_r,imax_r,jmin_r,jmax_r,kmin_r,kmax_r,nvar)
               end if
            case(1)
!!!               neiP=neiPE(1,iFaceRecv,iBlock)
!               if(do_prolonged .and. neiP/=iProc .and. &
!                    .not.unusedBlock_BP(neiBlock(1,iFaceRecv,iBlock),neiP))then

               iNodeNei = iNodeNei_IIIB(iRecv,jRecv,kRecv,iBlock)
               iProcNei = iTree_IA(Proc_, iNodeNei)
               if(do_prolonged .and. iProcNei /= iProc .and. &
                    iTree_IA(Status_,iNodeNei) /= Unused_)then

                  call buffer2face_state(buffer(ibuf,iside),&
                       imin_g,imax_g,jmin_g,jmax_g,kmin_g,kmax_g,nvar)
               end if
            end select
         end do ! iBlock
      end do ! iface

      if(DoTest)write(*,*)'msg_pass_faces finished: me, ifacemin, ifacemax=',&
           iProc,ifacemin,ifacemax

    end subroutine msg_pass_faces

    !==========================================================================

    subroutine setranges

      ! Set ranges orthogonal to iDim based on the value of iDim
      ! Include ghostcells for directions<idir if DoSendCorner is true.

      if(iDim>1.and.DoSendCorner)then
         imin_g=-1; imax_g=nI+2  ; imin_o=-1; imax_o=nI+2;
         imin_r= 0; imax_r=nI/2+1; iminR =-1; imaxR =nI+2;
         imin_p=-1; imax_p=nI*2+2; iminP = 0; imaxP =nI+1;
      else if(iDim/=1)then
         imin_g=1;  imax_g=nI;     imin_o=1;  imax_o=nI
         imin_r=1;  imax_r=nI/2;   iminR =1;  imaxR =nI;
         imin_p=1;  imax_p=nI*2;   iminP =1;  imaxP =nI;
      end if

      if(iDim>2.and.DoSendCorner)then
         jmin_g=-1; jmax_g=nJ+2  ; jmin_o=-1; jmax_o=nJ+2;
         jmin_r= 0; jmax_r=nJ/2+1; jminR =-1; jmaxR =nJ+2;
         jmin_p=-1; jmax_p=nJ*2+2; jminP = 0; jmaxP =nJ+1;
      else if(iDim/=2)then
         jmin_g=1;  jmax_g=nJ;     jmin_o=1;  jmax_o=nJ;
         jmin_r=1;  jmax_r=nJ/2;   jminR =1;  jmaxR =nJ;
         jmin_p=1;  jmax_p=nJ*2;   jminP =1;  jmaxP =nJ;
      endif

      if(iDim/=3)then
         kmin_g=1;  kmax_g=nK;     kmin_o=1;  kmax_o=nK;
         kmin_r=1;  kmax_r=nK/2;   kminR =1;  kmaxR =nK;
         kmin_p=1;  kmax_p=nK*2;   kminP =1;  kmaxP =nK;
      end if

      ! Set ranges in direction of iDim based on the value of iface

      ! Set direction and neighbor indexes for transverse direction
      iDir = 0; jDir = 0; kDir = 0
      iSend = 1; jSend = 1; kSend = 1
      iRecv = 1; jRecv = 1; kRecv = 1

      select case(iFaceSend)
      case(1)
         iDir = -1; iSend = 0; iRecv = 3
         iFaceRecv = 2
         imin_o=1;    imax_o=nWidth;
         imin_g=nI+1; imax_g=nI+nWidth;
         iminR =1;    imaxR=2*nWidth;
         imin_r=1;    imax_r=nWidth;
         iminP =1;    imaxP=nCoarseLayer;
         imin_p=1;    imax_p=2;
      case(3)
         jDir = -1; jSend = 0; jRecv = 3
         iFaceRecv = 4
         jmin_o=1;    jmax_o=nWidth;
         jmin_g=nJ+1; jmax_g=nJ+nWidth; 
         jminR =1;    jmaxR=2*nWidth;
         jmin_r=1;    jmax_r=nWidth; 
         jminP =1;    jmaxP=nCoarseLayer;
         jmin_p=1;    jmax_p=2;
      case(5)
         kDir = -1; kSend = 0; kRecv = 3
         iFaceRecv = 6
         kmin_o=1;    kmax_o=nWidth;  
         kmin_g=nK+1; kmax_g=nK+nWidth;
         kminR =1;    kmaxR =2*nWidth; 
         kmin_r=1;    kmax_r=nWidth;
         kminP =1;    kmaxP=nCoarseLayer;
         kmin_p=1;    kmax_p=2;
      case(2)
         iDir = 1; iSend = 3; iRecv = 0
         iFaceRecv = 1
         imin_o=nI-nWidth1;   imax_o=nI;
         imin_g=-nWidth1;     imax_g=0; 
         iminR=nI-nWidth2;    imaxR =nI;
         imin_r=nI/2-nWidth1; imax_r=nI/2;
         iminP =nI-nCoarse1; imaxP =nI;
         imin_p=nI*2-1;      imax_p=nI*2;
      case(4)
         jDir = 1; jSend = 3; jRecv = 0
         iFaceRecv = 3
         jmin_o=nJ-nWidth1;   jmax_o=nJ;
         jmin_g=-nWidth1;     jmax_g=0; 
         jminR =nJ-nWidth2;   jmaxR =nJ;
         jmin_r=nJ/2-nWidth1; jmax_r=nJ/2;
         jminP =nJ-nCoarse1; jmaxP =nJ;
         jmin_p=nJ*2-1;      jmax_p=nJ*2;
      case(6)
         kDir = 1; kSend = 3; kRecv = 0
         iFaceRecv = 5
         kmin_o=nK-nWidth1;   kmax_o=nK;
         kmin_g=-nWidth1;     kmax_g=0; 
         kminR=nK-nWidth2;    kmaxR=nK;
         kmin_r=nK/2-nWidth1; kmax_r=nK/2;
         kminP=nK-nCoarse1;  kmaxP=nK;
         kmin_p=nK*2-1;      kmax_p=nK*2;
      end select

      iSide = (iFaceSend-iFaceRecv+3)/2

    end subroutine setranges

    !==========================================================================

    subroutine setsubrange_g

      ! Select appropriate quarter of ghost cell layer

      select case(iFaceSend)
      case(1,2)
         imin_s=imin_g; imax_s=imax_g
         select case(isubface)
            ! Beware, case(2) and case(3) are swapped
         case(1)
            jmin_s=jmin_r; jmax_s=jmax_r; 
            kmin_s=kmin_r; kmax_s=kmax_r
         case(3)
            jmin_s=jmin_r+nJ/2; jmax_s=jmax_r+nJ/2; 
            kmin_s=kmin_r; kmax_s=kmax_r
         case(2)
            jmin_s=jmin_r; jmax_s=jmax_r; 
            kmin_s=kmin_r+nK/2; kmax_s=kmax_r+nK/2;
         case(4)
            jmin_s=jmin_r+nJ/2; jmax_s=jmax_r+nJ/2; 
            kmin_s=kmin_r+nK/2; kmax_s=kmax_r+nK/2;
         end select
      case(3,4)
         jmin_s=jmin_g; jmax_s=jmax_g
         select case(isubface)
            ! Beware, case(2) and case(3) are swapped
         case(1)
            imin_s=imin_r; imax_s=imax_r; 
            kmin_s=kmin_r; kmax_s=kmax_r
         case(3)
            imin_s=imin_r+nI/2; imax_s=imax_r+nI/2; 
            kmin_s=kmin_r;      kmax_s=kmax_r
         case(2)
            imin_s=imin_r;      imax_s=imax_r; 
            kmin_s=kmin_r+nK/2; kmax_s=kmax_r+nK/2;
         case(4)
            imin_s=imin_r+nI/2; imax_s=imax_r+nI/2; 
            kmin_s=kmin_r+nK/2; kmax_s=kmax_r+nK/2;
         end select
      case(5,6)
         kmin_s=kmin_g; kmax_s=kmax_g
         select case(isubface)
            ! Beware, case(2) and case(3) are not swapped
         case(1)
            imin_s=imin_r; imax_s=imax_r; 
            jmin_s=jmin_r; jmax_s=jmax_r
         case(2)
            imin_s=imin_r+nI/2; imax_s=imax_r+nI/2; 
            jmin_s=jmin_r;      jmax_s=jmax_r
         case(3)
            imin_s=imin_r;      imax_s=imax_r; 
            jmin_s=jmin_r+nJ/2; jmax_s=jmax_r+nJ/2;
         case(4)
            imin_s=imin_r+nI/2; imax_s=imax_r+nI/2; 
            jmin_s=jmin_r+nJ/2; jmax_s=jmax_r+nJ/2;
         end select
      end select

    end subroutine setsubrange_g

    !==========================================================================

    subroutine setsubrange_p

      ! Select appropriate quarter of prolonged originals
      ! Select the appropriate layer in the orthogonal direction for nWidth==1

      select case(iFaceSend)
      case(1,2)
         if(iFaceSend==1)then
            imin_s=1;           imax_s=nWidth;
         else
            imin_s=nI*2-nWidth1; imax_s=nI*2;
         end if
         select case(isubface)
            ! Beware, case(2) and case(3) are not swapped
         case(1)
            jmin_s=jmin_o; jmax_s=jmax_o; 
            kmin_s=kmin_o; kmax_s=kmax_o
         case(3)
            jmin_s=jmin_o+nJ; jmax_s=jmax_o+nJ; 
            kmin_s=kmin_o; kmax_s=kmax_o
         case(2)
            jmin_s=jmin_o; jmax_s=jmax_o; 
            kmin_s=kmin_o+nK; kmax_s=kmax_o+nK;
         case(4)
            jmin_s=jmin_o+nJ; jmax_s=jmax_o+nJ; 
            kmin_s=kmin_o+nK; kmax_s=kmax_o+nK;
         end select
      case(3,4)
         if(iFaceSend==3)then
            jmin_s=1;           jmax_s=nWidth;
         else
            jmin_s=nJ*2-nWidth1; jmax_s=nJ*2;
         end if
         select case(isubface)
            ! Beware, case(2) and case(3) are swapped
         case(1)
            imin_s=imin_o; imax_s=imax_o; 
            kmin_s=kmin_o; kmax_s=kmax_o
         case(3)
            imin_s=imin_o+nI; imax_s=imax_o+nI; 
            kmin_s=kmin_o;    kmax_s=kmax_o
         case(2)
            imin_s=imin_o;    imax_s=imax_o; 
            kmin_s=kmin_o+nK; kmax_s=kmax_o+nK;
         case(4)
            imin_s=imin_o+nI; imax_s=imax_o+nI; 
            kmin_s=kmin_o+nK; kmax_s=kmax_o+nK;
         end select
      case(5,6)
         if(iFaceSend==5)then
            kmin_s=1;           kmax_s=nWidth;
         else
            kmin_s=nK*2-nWidth1; kmax_s=nK*2;
         end if
         select case(isubface)
            ! Beware, case(2) and case(3) are not swapped
         case(1)
            imin_s=imin_o; imax_s=imax_o; 
            jmin_s=jmin_o; jmax_s=jmax_o
         case(2)
            imin_s=imin_o+nI; imax_s=imax_o+nI; 
            jmin_s=jmin_o;    jmax_s=jmax_o
         case(3)
            imin_s=imin_o;    imax_s=imax_o; 
            jmin_s=jmin_o+nJ; jmax_s=jmax_o+nJ;
         case(4)
            imin_s=imin_o+nI; imax_s=imax_o+nI; 
            jmin_s=jmin_o+nJ; jmax_s=jmax_o+nJ;
         end select
      end select

    end subroutine setsubrange_p

    !==========================================================================

    subroutine restrict_state

      ! Take average of 8 small cells in State_VGB and put it into re_buf

      re_buf(:,:,:,:)=0.5**nDim *(&
           State_VGB(:,iminR  :imaxR:2,jminR  :jmaxR:2,kminR  :kmaxR:2,iBlock)+&
           State_VGB(:,iminR+1:imaxR:2,jminR  :jmaxR:2,kminR  :kmaxR:2,iBlock)+&
           State_VGB(:,iminR  :imaxR:2,jminR+1:jmaxR:2,kminR  :kmaxR:2,iBlock)+&
           State_VGB(:,iminR  :imaxR:2,jminR  :jmaxR:2,kminR+1:kmaxR:2,iBlock)+&
           State_VGB(:,iminR+1:imaxR:2,jminR+1:jmaxR:2,kminR  :kmaxR:2,iBlock)+&
           State_VGB(:,iminR+1:imaxR:2,jminR  :jmaxR:2,kminR+1:kmaxR:2,iBlock)+&
           State_VGB(:,iminR  :imaxR:2,jminR+1:jmaxR:2,kminR+1:kmaxR:2,iBlock)+&
           State_VGB(:,iminR+1:imaxR:2,jminR+1:jmaxR:2,kminR+1:kmaxR:2,iBlock))

      ! redo the cells next to resolution changes if required by restrictface
      if(DoRestrictFace)call restrict_face_state

    end subroutine restrict_state

    !==========================================================================

    subroutine restrict_face_state

      ! Take average of 4 fine cells in State_VGB next to the resolution change
      ! and put the result into re_buf

      select case(iFaceSend)
      case(1)
         re_buf(:,imin_r,:,:)=0.25*(&
              State_VGB(:,iminR,jminR  :jmaxR:2,kminR  :kmaxR:2,iBlock)+&
              State_VGB(:,iminR,jminR+1:jmaxR:2,kminR  :kmaxR:2,iBlock)+&
              State_VGB(:,iminR,jminR  :jmaxR:2,kminR+1:kmaxR:2,iBlock)+&
              State_VGB(:,iminR,jminR+1:jmaxR:2,kminR+1:kmaxR:2,iBlock))
      case(2)
         re_buf(:,imax_r,:,:)=0.25*(&
              State_VGB(:,imaxR,jminR  :jmaxR:2,kminR  :kmaxR:2,iBlock)+&
              State_VGB(:,imaxR,jminR+1:jmaxR:2,kminR  :kmaxR:2,iBlock)+&
              State_VGB(:,imaxR,jminR  :jmaxR:2,kminR+1:kmaxR:2,iBlock)+&
              State_VGB(:,imaxR,jminR+1:jmaxR:2,kminR+1:kmaxR:2,iBlock))
      case(3)
         re_buf(:,:,jmin_r,:)=0.25*(&
              State_VGB(:,iminR  :imaxR:2,jminR,kminR  :kmaxR:2,iBlock)+&
              State_VGB(:,iminR+1:imaxR:2,jminR,kminR  :kmaxR:2,iBlock)+&
              State_VGB(:,iminR  :imaxR:2,jminR,kminR+1:kmaxR:2,iBlock)+&
              State_VGB(:,iminR+1:imaxR:2,jminR,kminR+1:kmaxR:2,iBlock))
      case(4)
         re_buf(:,:,jmax_r,:)=0.25*(&
              State_VGB(:,iminR  :imaxR:2,jmaxR,kminR  :kmaxR:2,iBlock)+&
              State_VGB(:,iminR+1:imaxR:2,jmaxR,kminR  :kmaxR:2,iBlock)+&
              State_VGB(:,iminR  :imaxR:2,jmaxR,kminR+1:kmaxR:2,iBlock)+&
              State_VGB(:,iminR+1:imaxR:2,jmaxR,kminR+1:kmaxR:2,iBlock))
      case(5)
         re_buf(:,:,:,kmin_r)=0.25*(&
              State_VGB(:,iminR  :imaxR:2,jminR  :jmaxR:2,kminR,iBlock)+&
              State_VGB(:,iminR+1:imaxR:2,jminR  :jmaxR:2,kminR,iBlock)+&
              State_VGB(:,iminR  :imaxR:2,jminR+1:jmaxR:2,kminR,iBlock)+&
              State_VGB(:,iminR+1:imaxR:2,jminR+1:jmaxR:2,kminR,iBlock))
      case(6)
         re_buf(:,:,:,kmax_r)=0.25*(&
              State_VGB(:,iminR  :imaxR:2,jminR  :jmaxR:2,kmaxR,iBlock)+&
              State_VGB(:,iminR+1:imaxR:2,jminR  :jmaxR:2,kmaxR,iBlock)+&
              State_VGB(:,iminR  :imaxR:2,jminR+1:jmaxR:2,kmaxR,iBlock)+&
              State_VGB(:,iminR+1:imaxR:2,jminR+1:jmaxR:2,kmaxR,iBlock))
      end select
    end subroutine restrict_face_state

    !==========================================================================

    subroutine prolong2(iVar)

      ! Second order prolongation State_VGB into pr_buf

      integer, intent(in) :: ivar

      !------------------------------------------------------------------------
      ! Obtain qsol with an extra layer of ghost cells for calculating gradients
      qsol(iminP-1:imaxP+1,jminP-1:jmaxP+1,kminP-1:kmaxP+1)=&
           State_VGB(iVar,iminP-1:imaxP+1,jminP-1:jmaxP+1,kminP-1:kmaxP+1,iBlock)

      ! The part without the extra ghost cells
      avrg=qsol(iminP:imaxP,jminP:jmaxP,kminP:kmaxP)

      ! Although we may have ghost cell values in directions/=iDim, we don't 
      ! use those so that one directional gradient can be used without sending
      ! ghost cells in any other direction.

      select case(prolong_type)
      case('central','lr','minmod')
         ! Use first order extrapolation in the other directions
         if(iDim/=1)then
            qsol(iminP-1,jminP:jmaxP,kminP:kmaxP)=&
                 qsol(iminP,jminP:jmaxP,kminP:kmaxP)
            qsol(imaxP+1,jminP:jmaxP,kminP:kmaxP)=&
                 qsol(imaxP,jminP:jmaxP,kminP:kmaxP)
         end if
         if(iDim/=2)then
            qsol(iminP:imaxP,jminP-1,kminP:kmaxP)=&
                 qsol(iminP:imaxP,jminP,kminP:kmaxP)
            qsol(iminP:imaxP,jmaxP+1,kminP:kmaxP)=&
                 qsol(iminP:imaxP,jmaxP,kminP:kmaxP)
         end if
         if(iDim/=3)then
            qsol(iminP:imaxP,jminP:jmaxP,kminP-1)=&
                 qsol(iminP:imaxP,jminP:jmaxP,kminP)
            qsol(iminP:imaxP,jminP:jmaxP,kmaxP+1)=&
                 qsol(iminP:imaxP,jminP:jmaxP,kmaxP)
         endif
      case('central2','lr2','minmod2','lr3')
         if(prolong_type/='lr3' .or. DoSendCorner)then
            ! Use second order extrapolation in the other directions
            if(iDim/=1)then
               qsol(iminP-1,jminP:jmaxP,kminP:kmaxP)=&
                    2*qsol(iminP  ,jminP:jmaxP,kminP:kmaxP)&
                    -qsol(iminP+1,jminP:jmaxP,kminP:kmaxP)

               qsol(imaxP+1,jminP:jmaxP,kminP:kmaxP)=&
                    2*qsol(imaxP  ,jminP:jmaxP,kminP:kmaxP)&
                    -qsol(imaxP-1,jminP:jmaxP,kminP:kmaxP)
            end if
            if(iDim/=2)then
               qsol(iminP:imaxP,jminP-1,kminP:kmaxP)=&
                    2*qsol(iminP:imaxP,jminP  ,kminP:kmaxP)&
                    -qsol(iminP:imaxP,jminP+1,kminP:kmaxP)

               qsol(iminP:imaxP,jmaxP+1,kminP:kmaxP)=&
                    2*qsol(iminP:imaxP,jmaxP  ,kminP:kmaxP)&
                    -qsol(iminP:imaxP,jmaxP-1,kminP:kmaxP)
            end if
            if(iDim/=3)then
               qsol(iminP:imaxP,jminP:jmaxP,kminP-1)=&
                    2*qsol(iminP:imaxP,jminP:jmaxP,kminP  )&
                    -qsol(iminP:imaxP,jminP:jmaxP,kminP+1)

               qsol(iminP:imaxP,jminP:jmaxP,kmaxP+1)=&
                    2*qsol(iminP:imaxP,jminP:jmaxP,kmaxP  )&
                    -qsol(iminP:imaxP,jminP:jmaxP,kmaxP-1)
            endif
         end if
      case default
         call CON_stop(&
              'Error 0 in prolong, unknown prolong_type='//prolong_type)
      end select ! prolong_type

      ! Calculate the gradients
      select case(prolong_type)
      case('central','central2')
         gradx=0.125*(&
              qsol(iminP+1:imaxP+1,jminP:jmaxP,kminP:kmaxP)-&
              qsol(iminP-1:imaxP-1,jminP:jmaxP,kminP:kmaxP))
         grady=0.125*(&
              qsol(iminP:imaxP,jminP+1:jmaxP+1,kminP:kmaxP)-&
              qsol(iminP:imaxP,jminP-1:jmaxP-1,kminP:kmaxP))
         gradz=0.125*(&
              qsol(iminP:imaxP,jminP:jmaxP,kminP+1:kmaxP+1)-&
              qsol(iminP:imaxP,jminP:jmaxP,kminP-1:kmaxP-1))
      case('lr','minmod','lr2','minmod2','lr3')
         gradxl=0.25*(&
              qsol(iminP  :imaxP  ,jminP:jmaxP,kminP:kmaxP)-&
              qsol(iminP-1:imaxP-1,jminP:jmaxP,kminP:kmaxP))
         gradxr=0.25*(&
              qsol(iminP+1:imaxP+1,jminP:jmaxP,kminP:kmaxP)-&
              qsol(iminP  :imaxP  ,jminP:jmaxP,kminP:kmaxP))
         gradyl=0.25*(&
              qsol(iminP:imaxP,jminP  :jmaxP  ,kminP:kmaxP)-&
              qsol(iminP:imaxP,jminP-1:jmaxP-1,kminP:kmaxP))
         gradyr=0.25*(&
              qsol(iminP:imaxP,jminP+1:jmaxP+1,kminP:kmaxP)-&
              qsol(iminP:imaxP,jminP  :jmaxP  ,kminP:kmaxP))
         gradzl=0.25*(&
              qsol(iminP:imaxP,jminP:jmaxP,kminP  :kmaxP  )-&
              qsol(iminP:imaxP,jminP:jmaxP,kminP-1:kmaxP-1))
         gradzr=0.25*(&
              qsol(iminP:imaxP,jminP:jmaxP,kminP+1:kmaxP+1)-&
              qsol(iminP:imaxP,jminP:jmaxP,kminP  :kmaxP  ))
         if(prolong_type=='minmod'.or.prolong_type=='minmod2')then
            gradx=sign(1.,gradxl)*max(0.,min(abs(gradxl),sign(1.,gradxl)*gradxr))
            grady=sign(1.,gradyl)*max(0.,min(abs(gradyl),sign(1.,gradyl)*gradyr))
            gradz=sign(1.,gradzl)*max(0.,min(abs(gradzl),sign(1.,gradzl)*gradzr))
         end if
      case default
         call CON_stop('Error 1 in prolong, unknown prolong_type='//prolong_type)
      end select

      ! Apply gradients
      select case(prolong_type)
      case('central','central2','minmod','minmod2')
         pr_buf1(imin_p  :imax_p:2,jmin_p  :jmax_p:2,kmin_p  :kmax_p:2)=&
              avrg-gradx-grady-gradz
         pr_buf1(imin_p+1:imax_p:2,jmin_p  :jmax_p:2,kmin_p  :kmax_p:2)=&
              avrg+gradx-grady-gradz
         pr_buf1(imin_p  :imax_p:2,jmin_p+1:jmax_p:2,kmin_p  :kmax_p:2)=&
              avrg-gradx+grady-gradz
         pr_buf1(imin_p  :imax_p:2,jmin_p  :jmax_p:2,kmin_p+1:kmax_p:2)=&
              avrg-gradx-grady+gradz
         pr_buf1(imin_p+1:imax_p:2,jmin_p+1:jmax_p:2,kmin_p  :kmax_p:2)=&
              avrg+gradx+grady-gradz
         pr_buf1(imin_p+1:imax_p:2,jmin_p  :jmax_p:2,kmin_p+1:kmax_p:2)=&
              avrg+gradx-grady+gradz
         pr_buf1(imin_p  :imax_p:2,jmin_p+1:jmax_p:2,kmin_p+1:kmax_p:2)=&
              avrg-gradx+grady+gradz
         pr_buf1(imin_p+1:imax_p:2,jmin_p+1:jmax_p:2,kmin_p+1:kmax_p:2)=&
              avrg+gradx+grady+gradz
      case('lr','lr2','lr3')
         pr_buf1(imin_p  :imax_p:2,jmin_p  :jmax_p:2,kmin_p  :kmax_p:2)=&
              avrg-gradxl-gradyl-gradzl
         pr_buf1(imin_p+1:imax_p:2,jmin_p  :jmax_p:2,kmin_p  :kmax_p:2)=&
              avrg+gradxr-gradyl-gradzl
         pr_buf1(imin_p  :imax_p:2,jmin_p+1:jmax_p:2,kmin_p  :kmax_p:2)=&
              avrg-gradxl+gradyr-gradzl
         pr_buf1(imin_p  :imax_p:2,jmin_p  :jmax_p:2,kmin_p+1:kmax_p:2)=&
              avrg-gradxl-gradyl+gradzr
         pr_buf1(imin_p+1:imax_p:2,jmin_p+1:jmax_p:2,kmin_p  :kmax_p:2)=&
              avrg+gradxr+gradyr-gradzl
         pr_buf1(imin_p+1:imax_p:2,jmin_p  :jmax_p:2,kmin_p+1:kmax_p:2)=&
              avrg+gradxr-gradyl+gradzr
         pr_buf1(imin_p  :imax_p:2,jmin_p+1:jmax_p:2,kmin_p+1:kmax_p:2)=&
              avrg-gradxl+gradyr+gradzr
         pr_buf1(imin_p+1:imax_p:2,jmin_p+1:jmax_p:2,kmin_p+1:kmax_p:2)=&
              avrg+gradxr+gradyr+gradzr
      case default
         call CON_stop('Error 2 in prolong, unknown prolong_type='//prolong_type)
      end select

      pr_buf(iVar,:,:,:)=pr_buf1

    end subroutine prolong2

    !==========================================================================

    subroutine prolong1_state

      ! First order prolongation of State_VGB into pr_buf

      !------------------------------------------------------------------------
      ! First order prolongation
      avrg_state=State_VGB(:,iminP:imaxP,jminP:jmaxP,kminP:kmaxP,iBlock)

      pr_buf(:,imin_p  :imax_p:2,jmin_p  :jmax_p:2,kmin_p  :kmax_p:2)&
           = avrg_state
      pr_buf(:,imin_p+1:imax_p:2,jmin_p  :jmax_p:2,kmin_p  :kmax_p:2)&
           = avrg_state
      pr_buf(:,imin_p  :imax_p:2,jmin_p+1:jmax_p:2,kmin_p  :kmax_p:2)&
           = avrg_state
      pr_buf(:,imin_p  :imax_p:2,jmin_p  :jmax_p:2,kmin_p+1:kmax_p:2)&
           = avrg_state
      pr_buf(:,imin_p+1:imax_p:2,jmin_p+1:jmax_p:2,kmin_p  :kmax_p:2)&
           = avrg_state
      pr_buf(:,imin_p+1:imax_p:2,jmin_p  :jmax_p:2,kmin_p+1:kmax_p:2)&
           = avrg_state
      pr_buf(:,imin_p  :imax_p:2,jmin_p+1:jmax_p:2,kmin_p+1:kmax_p:2)&
           = avrg_state
      pr_buf(:,imin_p+1:imax_p:2,jmin_p+1:jmax_p:2,kmin_p+1:kmax_p:2)&
           = avrg_state

    end subroutine prolong1_state

    !==========================================================================

    subroutine prolong2_state

      ! Prolongation for 2nd order TVD resolution change.
      ! Put 2 coarse layers of State_VGB into pr_buf

      !------------------------------------------------------------------------
      ! First order prolongation
      avrg_state=State_VGB(:,iminP:imaxP,jminP:jmaxP,kminP:kmaxP,iBlock)

      select case(iFaceSend)
      case(1, 2)
         pr_buf(:,imin_p:imax_p,jmin_p  :jmax_p:2,kmin_p  :kmax_p:2) = avrg_state
         pr_buf(:,imin_p:imax_p,jmin_p+1:jmax_p:2,kmin_p  :kmax_p:2) = avrg_state
         pr_buf(:,imin_p:imax_p,jmin_p  :jmax_p:2,kmin_p+1:kmax_p:2) = avrg_state
         pr_buf(:,imin_p:imax_p,jmin_p+1:jmax_p:2,kmin_p+1:kmax_p:2) = avrg_state
      case(3, 4)
         pr_buf(:,imin_p  :imax_p:2,jmin_p:jmax_p,kmin_p  :kmax_p:2) = avrg_state
         pr_buf(:,imin_p+1:imax_p:2,jmin_p:jmax_p,kmin_p  :kmax_p:2) = avrg_state
         pr_buf(:,imin_p  :imax_p:2,jmin_p:jmax_p,kmin_p+1:kmax_p:2) = avrg_state
         pr_buf(:,imin_p+1:imax_p:2,jmin_p:jmax_p,kmin_p+1:kmax_p:2) = avrg_state
      case(5, 6)
         pr_buf(:,imin_p  :imax_p:2,jmin_p  :jmax_p:2,kmin_p:kmax_p) = avrg_state
         pr_buf(:,imin_p+1:imax_p:2,jmin_p  :jmax_p:2,kmin_p:kmax_p) = avrg_state
         pr_buf(:,imin_p  :imax_p:2,jmin_p+1:jmax_p:2,kmin_p:kmax_p) = avrg_state
         pr_buf(:,imin_p+1:imax_p:2,jmin_p+1:jmax_p:2,kmin_p:kmax_p) = avrg_state
      end select

    end subroutine prolong2_state

    !==========================================================================
    subroutine send_equal

      ! Same level

      iNodeNei = iNodeNei_IIIB(iSend,jSend,kSend,iBlock)
      iProcNei = iTree_IA(Proc_, iNodeNei)
      iBlockNei= iTree_IA(Block_,iNodeNei)

!!!      neiB=neiBlock(1,iface,iBlock)
!!!     neiP=neiPE(1,iface,iBlock)
!!!      if(unusedBlock_BP(neiB,neiP)) RETURN

      if(iProcNei==iProc)then
         ! Local copy
         State_VGB(:,imin_g:imax_g,jmin_g:jmax_g,kmin_g:kmax_g,iBlockNei)= &
              State_VGB(:,imin_o:imax_o,jmin_o:jmax_o,kmin_o:kmax_o,iBlock)
      else
         ! Remote send
         itag = 100*iBlockNei+10*iFaceSend+1
         if(DoTest)write(*,*)'Remote equal send, me,iFaceSend,itag=',&
              iProc,iFaceSend,itag

         eq_buf = State_VGB(:,imin_o:imax_o,jmin_o:jmax_o,kmin_o:kmax_o,iBlock)

         call MPI_Rsend(eq_buf,isize, MPI_REAL, iProcNei, itag, iComm, iError)
      end if
    end subroutine send_equal

    !==========================================================================
    subroutine send_restricted

      ! Neighbor is coarser

      iNodeNei = iNodeNei_IIIB(iSend,jSend,kSend,iBlock)
      iProcNei = iTree_IA(Proc_, iNodeNei)
      iBlockNei= iTree_IA(Block_,iNodeNei)

!!!   neiB=neiBlock(1,iface,iBlock)
!!!   neiP=neiPE(1,iface,iBlock)
!!!      if (unusedBlock_BP(neiB,neiP)) RETURN

      ! Restrict the R range of cells into _r

      call restrict_state

      ! Subface index =1,2,3, or 4 with respect to the coarse neighbor
!!!      ichild=BLKneighborCHILD(0,0,0,1,iBlock)
      if(ichild<1.or.ichild>8)then
         write(*,*)'me,iBlock,iFaceSend,neiL,iProcNei,iBlockNei,ichild:',&
              iProc,iBlock,iFaceSend,neiL,iProcNei,iBlockNei,ichild
         call CON_stop('error in message_pass: ichild incorrect')
      end if
!!!      isubface=child2subface(ichild,iface)
      if(isubface==0)then
         write(*,*)'me,iBlock,iFaceSend,neiL,iProcNei,iBlockNei,ichild:',&
              iProc,iBlock,iFaceSend,neiL,iProcNei,iBlockNei,ichild
         call CON_stop('error in message_pass: isubface=0')
      end if
      if(iProcNei==iProc)then
         ! Local copy into appropriate subface
         call setsubrange_g
         State_VGB(:,imin_s:imax_s,jmin_s:jmax_s,kmin_s:kmax_s,iBlockNei)=re_buf
      else
         ! Remote send
         itag = 100*iBlockNei+10*iFaceSend+isubface
         if(DoTest)write(*,*)&
              'Remote restricted send, me,iFaceSend,itag=',&
              iProc,iFaceSend,itag
         call MPI_Rsend(re_buf, isize_r, MPI_REAL, iProcNei, itag,iComm,iError)
      end if
    end subroutine send_restricted

    !==========================================================================
    subroutine send_prolonged

      integer :: iVar

      ! Neighbor is finer, so prolong the _o range of cells into _p
      if(nProlongOrder > 1)then
         ! Unoptimized 2nd order prolongation
         do iVar=1,nVar
            call prolong2(iVar)
         end do
      else 
         if(nCoarseLayer == 1)then
            ! Optimized 1st order prolongation
            call prolong1_state
         else
            ! Prolong 2 coarse layers for 2nd order TVD resolution change
            call prolong2_state
         end if
      end if

      do isubface=1,4

!         iNodeNei = iNodeNei_IIIB(iSend,jSend,kSend,iBlock)
!         iProcNei = iTree_IA(Proc_, iNodeNei)
!         iBlockNei= iTree_IA(Block_,iNodeNei)

!!!      neiB=neiBlock(isubface,iface,iBlock)
!!!      neiP=neiPE(isubface,iface,iBlock)
!!!         if(unusedBlock_BP(neiB,neiP)) CYCLE

         call setsubrange_p
         if(iProcNei==iProc)then
            ! Local copy of appropriate subface
            State_VGB(:,imin_g:imax_g,jmin_g:jmax_g,kmin_g:kmax_g,iBlockNei)=&
                 pr_buf(:,imin_s:imax_s,jmin_s:jmax_s,kmin_s:kmax_s)
         else
            ! Remote send
            itag = 100*iBlockNei+10*iFaceSend+1
            if(DoTest)write(*,*)&
                 'Remote prolong send, me,iFaceSend,itag=',&
                 iProc,iFaceSend,itag

            pr_buf_s=pr_buf(1:nvar,imin_s:imax_s,jmin_s:jmax_s,kmin_s:kmax_s)

            call MPI_Rsend(pr_buf_s,isize,MPI_REAL,iProcNei,itag,iComm,iError)

         end if
      end do

    end subroutine send_prolonged

    !==========================================================================
    subroutine remote_receive(nsubface, isize1)

      ! Receive nsubface subfaces of size isize1

      integer, intent(in) :: nsubface, isize1
      !------------------------------------------------------------------------

      do isubface=1, nsubface

!!!         neiP=neiPE(isubface,iFaceRecv,iBlock)
!!!         neiB=neiBlock(isubface,iFaceRecv,iBlock)
!!!         if(unusedBlock_BP(neiB,neiP)) CYCLE

         iNodeNei = iNodeNei_IIIB(iRecv,jRecv,kRecv,iBlock)
         if(iTree_IA(Status_,iNodeNei)==Unused_) CYCLE

         iProcNei = iTree_IA(Proc_, iNodeNei)
         if(iProcNei==iProc) CYCLE

         ! Remote receive
         itag = 100*iBlock+10*iFaceSend+isubface
         if(DoTest)write(*,*)&
              'Remote recieve, me,iFaceSend,itag,neiL=',&
              iProc,iFaceSend,itag

         ibuf = (iBlock-1)*isize+(isubface-1)*isize_r+1
         call MPI_irecv(buffer(ibuf,iside), isize1, MPI_REAL, &
              iProcNei, itag, iComm, request, iError)
         number_receive_requests = number_receive_requests + 1
         receive_requests(number_receive_requests) = request
      end do

    end subroutine remote_receive

    !==========================================================================
    subroutine buffer2face_state(&
         buf,imin_g,imax_g,jmin_g,jmax_g,kmin_g,kmax_g,nvar)

      integer, intent(in) :: imin_g,imax_g,jmin_g,jmax_g,kmin_g,kmax_g,nvar
      real, dimension(nVar,imin_g:imax_g,jmin_g:jmax_g,kmin_g:kmax_g),&
           intent(inout) :: buf
      !------------------------------------------------------------------------

      ! Read buffer into State_VGB

      State_VGB(:,imin_g:imax_g,jmin_g:jmax_g,kmin_g:kmax_g,iBlock)=buf

    end subroutine buffer2face_state

    !==========================================================================
    subroutine buffer2subfaces_state(&
         buf,imin_r,imax_r,jmin_r,jmax_r,kmin_r,kmax_r,nvar)

      integer, intent(in) :: imin_r,imax_r,jmin_r,jmax_r,kmin_r,kmax_r,nvar
      real, dimension(nVar,imin_r:imax_r,jmin_r:jmax_r,kmin_r:kmax_r,4),&
           intent(inout) :: buf
      !------------------------------------------------------------------------

      ! Loop over 4 subfaces to read restricted values
      do isubface=1,4
!!!         neiP = neiPE(isubface,iFaceRecv,iBlock)
         if(iProcNei==iProc) CYCLE
!!!         if (unusedBlock_BP(neiBlock(isubface,iFaceRecv,iBlock),neiP)) CYCLE
         call setsubrange_g

         State_VGB(:,imin_s:imax_s,jmin_s:jmax_s,kmin_s:kmax_s,iBlock)=&
              buf(:,:,:,:,isubface)
      end do ! isubface

    end subroutine buffer2subfaces_state
    !==========================================================================

  end subroutine message_pass_cell

end module BATL_pass_cell
