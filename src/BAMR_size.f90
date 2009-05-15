module BAMR_size

  implicit none

  SAVE

  ! Maximum number of blocks per processor
  integer :: MaxBlock = 0

  ! Number of cells per block in each direction
  integer, parameter :: nI = 8, nJ = 8, nK = 1

  ! number of ghost cells
  integer, parameter :: nG = 2  
  
  integer, parameter :: &
       MinI = 1 - nG, MaxI = nI + nG, &
       MinJ = 1 - nG, MaxJ = nJ + nG, &
       MinK = 1 - nG, MaxK = nK + nG

end module BAMR_size

