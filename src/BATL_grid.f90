module BATL_grid

  use BATL_size
  use BATL_tree

  private ! except

  public :: init_grid
  public :: create_grid_block

  real :: CoordMin_D(MaxDim), CoordMax_D(MaxDim)
  real, allocatable :: CoordMin_DB(:,:), CoordMax_DB(:,:)
  real, allocatable :: CellSize_DB(:,:)
  real, allocatable :: CellVolume_CB(:,:,:,:)
  real, allocatable :: Xyz_DGB(:,:,:,:,:)

contains
  !============================================================================
  subroutine grid_init(CoordMinIn_D, CoordMaxIn_D)

    intent(in):: CoordMinIn_D, CoordMaxIn_D

    allocate(CoordMin_DB(MaxDim,MaxBlock))
    allocate(CoordMax_DB(MaxDim,MaxBlock))
    allocate(CellSize_DB(MaxDim,MaxBlock))
    allocate(CellVolume_CB(nI,nJ,nK,MaxBlock))
    allocate(Xyz_DGB(MaxDim,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))

    CoordMin_D = CoordMinIn_D
    CoordMax_D = CoordMaxIn_D

  end subroutine grid_init
  !===========================================================================
  subroutine create_grid_block(iBlock)

    integer, intent(in):: iBlock

    real :: PositionMin_D(MaxDim), PositionMax_D(MaxDim)
    !----------------------------------------------------------------------
    call get_block_position(iBlock, PositionMin_D, PositionMax_D)

    CoordMin_DB(:,iBlock)= CoordMin_D + (CoordMax_D - CoordMin_D)*PositionMin_D
    CoordMax_DB(:,iBlock)= CoordMin_D + (CoordMax_D - CoordMin_D)*PositionMax_D


  end subroutine create_grid_block

end module BATL_grid
