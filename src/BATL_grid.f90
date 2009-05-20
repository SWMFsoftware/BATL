module BATL_grid

  use BATL_size
  use BATL_tree
  use BATL_geometry, ONLY: IsCartesian

  implicit none

  private ! except

  public :: init_grid
  public :: create_grid_block

  real :: CoordMin_D(MaxDim)                  ! Min gen. coordinates of domain
  real :: CoordMax_D(MaxDim)                  ! Max gen. coordinates of domain
  real, allocatable :: CoordMin_DB(:,:)       ! Min gen. coordinates of a block
  real, allocatable :: CoordMax_DB(:,:)       ! Max gen. coordinates of a block
  real, allocatable :: CellSize_DB(:,:)       ! Cell size in gen. coordinates
  real, allocatable :: CellFace_DB(:,:)       ! Cell faces for Cartesian grids
  real, allocatable :: CellFace_DFB(:,:,:,:,:)! Cell faces for general grids
  real, allocatable :: CellVolume_B(:)        ! Cell volume for Cartesian grids
  real, allocatable :: CellVolume_GB(:,:,:,:) ! Cell volume for general grids
  
  real, allocatable :: Xyz_DGB(:,:,:,:,:)     ! Cartesian cell centers coords

contains
  !============================================================================
  subroutine init_grid(CoordMinIn_D, CoordMaxIn_D)

    real, intent(in):: CoordMinIn_D(MaxDim), CoordMaxIn_D(MaxDim)
    !-------------------------------------------------------------------------
    allocate(CoordMin_DB(MaxDim,MaxBlock))
    allocate(CoordMax_DB(MaxDim,MaxBlock))
    allocate(CellSize_DB(MaxDim,MaxBlock))

    allocate(CellFace_DB(MaxDim,MaxBlock))
    if(.not.IsCartesian) &
         allocate(CellFace_DFB(MaxDim,1:nI+1,1:nJ+1,1:nK+1,MaxBlock))

    allocate(CellVolume_B(MaxBlock))
    if(.not.IsCartesian) &
         allocate(CellVolume_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    allocate(Xyz_DGB(MaxDim,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))

    CoordMin_D = CoordMinIn_D
    CoordMax_D = CoordMaxIn_D

  end subroutine init_grid
  !===========================================================================
  subroutine create_grid_block(iBlock)

    integer, intent(in):: iBlock

    real :: PositionMin_D(MaxDim), PositionMax_D(MaxDim)
    integer :: i, j, k
    !----------------------------------------------------------------------
    call get_block_position(iBlock, PositionMin_D, PositionMax_D)

    CoordMin_DB(:,iBlock)= CoordMin_D + (CoordMax_D - CoordMin_D)*PositionMin_D
    CoordMax_DB(:,iBlock)= CoordMin_D + (CoordMax_D - CoordMin_D)*PositionMax_D

    CellSize_DB(:,iBlock) = (CoordMax_DB(:,iBlock) - CoordMin_DB(:,iBlock)) &
         / nIJK_D

    if(IsCartesian)then
       CellVolume_B(iBlock) = product(CellSize_DB(:,iBlock))

       if(allocated(CellVolume_GB)) &
            CellVolume_GB(:,:,:,iBlock) = CellVolume_B(iBlock)

       CellFace_DB(:,iBlock) = CellVolume_B(iBlock) / CellSize_DB(:,iBlock)

       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
          Xyz_DGB(:,i,j,k,iBlock) = CoordMin_DB(:,iBlock) + &
               ( (/i, j, k/) - 0.5 ) * CellSize_DB(:,iBlock)
       end do; end do; end do
    else
       stop
    end if

  end subroutine create_grid_block
  !===========================================================================

  subroutine test_grid

    integer :: iBlock, nBlockAll, Int_D(MaxDim)
    real:: CoordTest_D(MaxDim)
 
    character(len=*), parameter :: NameSub = 'test_tree'
    !-----------------------------------------------------------------------

    write(*,*)'Testing init_mod_tree'
    call init_mod_tree(50, 100)
    call set_root_block( (/1,2,3/), (/.true., .true., .false./) )

    
  end subroutine test_grid

end module BATL_grid
