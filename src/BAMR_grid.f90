module BAMR_grid

  use BAMR_size
  use BAMR_tree

  real :: CoordMin_D(MaxDim), CoordMax_D(MaxDim)
  real, allocatable :: CoordMin_DB(:,:), CoordMax_DB(:,:)
  real, allocatable :: CellSize_DB(:,:)
  real, allocatable :: CellVolume_CB(:,:,:,:)
  real, allocatable :: Xyz_DGB(:,:,:,:,:)

contains

  subroutine grid_init

    allocate(CoordMin_DB(MaxDim,MaxBlock))
    allocate(CoordMax_DB(MaxDim,MaxBlock))
    allocate(CellSize_DB(MaxDim,MaxBlock))
    allocate(CellVolume_CB(nI,nJ,nK,MaxBlock))
    allocate(Xyz_DGB(MaxDim,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))

  end subroutine grid_init
  

end module BAMR_grid
