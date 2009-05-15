module BAMR_grid

  use BAMR_size
  use BAMR_tree

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

    integer :: iLevel, iDim
    !----------------------------------------------------------------------

    iLevel = iTree_IA(Level_, iBlock)

    CoordMin_DB(1:nDim,iBlock) = CoordMin_D(1:nDim) + &
         (CoordMax_D(1:nDim) - CoordMin_D(1:nDim))* &
         (iTree_IA(Coord1_:CoordLast_,iBlock)-1.0) &
         /MaxCoord_I(iLevel)/nRoot_D(1:nDim)

    CoordMax_DB(1:nDim,iBlock) = CoordMin_D(1:nDim) + &
         (CoordMax_D(1:nDim) - CoordMin_D(1:nDim))* &
         (iTree_IA(Coord1_:CoordLast_,iBlock)+0.0) &
         /MaxCoord_I(iLevel)/nRoot_D(1:nDim)

    ! Only one root block in the ignored direction
    do iDim = nDim+1, MaxDim
       CoordMin_DB(iDim,iBlock) =  CoordMin_D(iDim)
       CoordMax_DB(iDim,iBlock) =  CoordMax_D(iDim)
    end do

  end subroutine create_grid_block

end module BAMR_grid
