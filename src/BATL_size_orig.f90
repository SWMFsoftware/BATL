module BATL_size

  implicit none

  SAVE

  ! Dimensionality of grid and AMR
  integer, parameter :: MaxDim = 3    ! This has to be 3 all the time
  integer, parameter :: nDim   = 3    ! Number of not ignored dimensions
  integer, parameter :: nDimAmr = 3  ! Number of refined dimensions

  ! Note: 3 = MaxDim >= nDim >= nDimAmr >= 1

  ! Maximum number of blocks per processor
  integer :: MaxBlock = 0

  ! Largest block index
  integer :: nBlock = 0

  ! Number of cells per block in each direction
  integer, parameter :: nI = 8, nJ = 4, nK = 2

  ! Array for block size
  integer, parameter:: &
       nIJK_D(MaxDim) = (/ nI, nJ, nK /)

  ! number of ghost cells
  integer, parameter :: nG = 2  
  
  integer, parameter :: &
       MinI = 1 - nG, MaxI = nI + nG, &
       MinJ = 1 - nG, MaxJ = nJ + nG, &
       MinK = 1 - nG, MaxK = nK + nG

  ! Refinement ratios in the 3 dimensions (depends on nDimAmr)
  integer, parameter:: &
       iRatio = 2, jRatio = min(2,nDimAmr), kRatio = max(1,nDimAmr-1)

  ! Array of refinement ratios
  integer, parameter:: iRatio_D(MaxDim) = (/ iRatio, jRatio, kRatio /)
  
  ! Inverse volume ratio for Cartesian case
  real, parameter:: InvIjkRatio = 1.0/(iRatio*jRatio*kRatio)

end module BATL_size

