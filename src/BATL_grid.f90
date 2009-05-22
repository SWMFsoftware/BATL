module BATL_grid

  use BATL_size
  use BATL_tree
  use BATL_geometry, ONLY: IsCartesian

  implicit none

  private ! except

  public :: init_grid
  public :: create_grid_block
  public :: test_grid

  logical :: DoInitializeGrid = .true.

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
    if(.not. DoInitializeGrid) RETURN

    DoInitializeGrid = .false.

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
  subroutine clean_grid

    if(DoInitializeGrid) RETURN

    DoInitializeGrid = .true.

    deallocate(CoordMin_DB, CoordMax_DB, CellSize_DB, CellFace_DB, &
         CellVolume_B, Xyz_DGB)
    if(allocated(CellFace_DFB)) deallocate(CellFace_DFB)
    if(allocated(CellVolume_GB))deallocate(CellVolume_GB)

    CoordMin_D =  0.0
    CoordMax_D = -1.0

  end subroutine clean_grid

  !===========================================================================

  subroutine create_grid_block(iBlock)

    integer, intent(in):: iBlock

    real :: PositionMin_D(MaxDim), PositionMax_D(MaxDim)
    integer :: i, j, k
    !----------------------------------------------------------------------
    call get_tree_position(iBlock, PositionMin_D, PositionMax_D)

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

  subroutine show_grid_block(iBlock)

    integer, intent(in):: iBlock
    !------------------------------------------------------------------------
    write(*,*)'show_grid_block for iBlock=',iBlock
    write(*,*)'CoordMin  =', CoordMin_DB(:,iBlock)
    write(*,*)'CoordMax  =', CoordMax_DB(:,iBlock)
    write(*,*)'CellSize  =', CellSize_DB(:,iBlock)
    write(*,*)'CellFace  =', CellFace_DB(:,iBlock)
    write(*,*)'CellVolume=', CellVolume_B(iBlock)
    write(*,*)'Xyz( 1, 1, 1)=', Xyz_DGB(:, 1, 1, 1,iBlock)
    write(*,*)'Xyz(nI, 1, 1)=', Xyz_DGB(:,nI, 1, 1,iBlock)
    write(*,*)'Xyz( 1,nJ, 1)=', Xyz_DGB(:, 1,nJ, 1,iBlock)
    write(*,*)'Xyz( 1, 1,nK)=', Xyz_DGB(:, 1, 1,nK,iBlock)
    write(*,*)'Xyz(nI,nJ,nK)=', Xyz_DGB(:,nI,nJ,nK,iBlock)

  end subroutine show_grid_block

  !===========================================================================

  subroutine test_grid

    integer :: iBlock, nBlockAll, Int_D(MaxDim)
    real:: DomainMin_D(MaxDim) = (/ 1.0, 2.0, 3.0 /)
    real:: DomainMax_D(MaxDim) = (/ 4.0, 6.0, 9.0 /)

    character(len=*), parameter :: NameSub = 'test_grid'
    !-----------------------------------------------------------------------

    write(*,*)'Testing init_grid'
    write(*,*)'nDimTree, nIJK_D=', nDimTree, nIJK_D
    call init_tree(50, 100)
    call init_grid( DomainMin_D, DomainMax_D )
    call set_tree_root( (/1,2,3/), (/.true., .true., .false./) )

    call refine_tree_block(6)

    write(*,*)'Testing create_grid_block'
    do iBlock = 1, 6 + 2**nDimTree
       call create_grid_block(iBlock)
       call show_grid_block(iBlock)
    end do

    write(*,*)'Testing clean_grid'
    call clean_grid
    
  end subroutine test_grid

end module BATL_grid
