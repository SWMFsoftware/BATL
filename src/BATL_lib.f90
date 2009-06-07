module BATL_lib

  ! Collection of all public methods and data that an application can access

  use BATL_size
  use BATL_mpi
  use BATL_tree
  use BATL_geometry
  use BATL_grid
  use BATL_pass_cell

  implicit none

  private ! except

  ! Inherited from BATL_size
  public:: MaxDim, nDim
  public:: MaxBlock, nBlock
  public:: nI, nJ, nK, nIJK_D
  public:: MinI, MaxI, MinJ, MaxJ, MinK, MaxK

  ! Inherited from BATL_mpi
  public:: init_mpi, clean_mpi, barrier_mpi
  public:: iComm, nProc, iProc          

  ! Inherited from BATL_tree
  public:: Unused_B
  public:: DiLevelNei_IIIB, iNodeNei_IIIB

  ! Inherited from BATL_geometry
  public:: TypeGeometry, IsPeriodic_D

  ! Inherited from BATL_grid
  public:: CoordMin_D, CoordMax_D, CoordMin_DB, CoordMax_DB, CellSize_DB
  public:: Xyz_DGB
  public:: CellFace_DB, CellFace_DFB
  public:: CellVolume_B, CellVolume_GB

  ! Inherited from BATL_pass_cell
  public:: message_pass_cell

  ! Public methods of this module
  public:: init_batl
  public:: clean_batl

contains
  !============================================================================
  subroutine init_batl(&
       CoordMinIn_D, CoordMaxIn_D, MaxBlockIn, &
       TypeGeometryIn, IsPeriodicIn_D, nRootIn_D)

    integer, intent(in):: MaxBlockIn
    real,    intent(in):: CoordMinIn_D(nDim), CoordMaxIn_D(nDim)

    integer,          optional, intent(in):: nRootIn_D(nDim)
    character(len=*), optional, intent(in):: TypeGeometryIn
    logical,          optional, intent(in):: IsPeriodicIn_D(nDim)
    !-------------------------------------------------------------------------
    call init_tree(MaxBlockIn)
    call init_geometry(TypeGeometryIn, IsPeriodicIn_D)
    call init_grid(CoordMinIn_D, CoordMaxIn_D)
    call set_tree_root(nRootIn_D)
    call distribute_tree(.true.)
    call create_grid

  end subroutine init_batl
  !============================================================================
  subroutine clean_batl
    call clean_grid
    call clean_tree
  end subroutine clean_batl
  !============================================================================

end module BATL_lib
