module BATL_pass_cell

  use BATL_size, ONLY: MaxBlock, nBlock, &
       nI, nJ, nK, MinI, MaxI, MinJ, MaxJ,MinK, MaxK, nDim
  use BATL_tree, ONLY: Unused_B, iNodeNei_IIIB, DiLevelNei_IIIB, &
       iTree_IA, Proc_, Block_
  use ModNumConst, ONLY: i_DD

  use BATL_grid, ONLY: CoordMin_DB !!!

  implicit none
  SAVE

  private ! except
  
  public message_pass_cell

contains

  subroutine message_pass_cell(nVar, State_VGB, nWidthIn)

    integer, intent(in) :: nVar
    real, intent(inout) :: &
         State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

    integer, optional, intent(in):: nWidthIn

    integer:: nWidth

    integer :: iFaceRecv, iBlock, iNodeNei, iBlockNei

    integer :: &
         iRMin, iRmax, jRMin, jRmax, kRMin, kRMax, &
         iSMin, iSmax, jSMin, jSmax, kSMin, kSMax

    integer :: Di, Dj, Dk
    !--------------------------------------------------------------------------
    nWidth = 2
    if(present(nWidthIn))nWidth = nWidthIn

    do iFaceRecv = 1, 2*nDim
       call set_ranges

       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE

          iNodeNei  = iNodeNei_IIIB(Di,Dj,Dk,iBlock)
          iBlockNei = iTree_IA(Block_,iNodeNei)

          !write(*,*)'!!! iBlock,    CoordMin_DB=', &
          !     iBlock, CoordMin_DB(:,iBlock)
          !write(*,*)'!!! iBlockNei, CoordMin_DB=', &
          !     iBlockNei, CoordMin_DB(:,iBlockNei)
          
          State_VGB(:,iRMin:iRmax,jRMin:jRmax,kRMin:kRMax,iBlock) = &
               State_VGB(:,iSMin:iSMax,jSMin:jSMax,kSMin:kSMax,iBlockNei)
          
       end do
    end do

  contains

    subroutine set_ranges

      ! Set indexes for transverse direction
      Di = 1; Dj  = 1; Dk = 1
      iSMin = 1 ; jSMin =  1; kSMin =  1
      iSMax = nI; jSMax = nJ; kSMax = nK
      iRMin = 1 ; jRMin =  1; kRMin =  1
      iRMax = nI; jRMax = nJ; kRMax = nK

      ! Set indexes in normal direction
      select case(iFaceRecv)
      case(1)
         Di = 0
         iRMin =    1 - nWidth; iRMax = 0
         iSMin = nI+1 - nWidth; iSMax = nI
      case(2)
         Di = 3
         iRMin = nI+1; iRMax = nI + nWidth
         iSMin =    1; iSMax = nWidth
      case(3)
         Dj = 0
         jRMin =    1 - nWidth; jRMax = 0
         jSMin = nJ+1 - nWidth; jSMax = nJ
      case(4)
         Dj = 3
         jRMin = nJ+1; jRMax = nJ + nWidth
         jSMin =    1; jSMax = nWidth
      case(5)
         Dk = 0
         kRMin =    1 - nWidth; kRMax = 0
         kSMin = nK+1 - nWidth; kSMax = nK
      case(6)
         Dk = 3
         kRMin = nK+1; kRMax = nK + nWidth
         kSMin =    1; kSMax = nWidth
      end select

    end subroutine set_ranges

  end subroutine message_pass_cell

end module BATL_pass_cell
