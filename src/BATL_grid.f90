module BATL_grid

  use BATL_size
  use BATL_tree
  use BATL_geometry, ONLY: IsCartesian

  implicit none

  private ! except

  public :: init_grid
  public :: clean_grid
  public :: create_grid
  public :: create_grid_block
  public :: test_grid

  real, public:: &
       CoordMin_D(MaxDim),      & ! Min gen. coords of domain
       CoordMax_D(MaxDim)         ! Max gen. coordinates of domain

  real, public, allocatable::   &
       CoordMin_DB(:,:),        & ! Min gen. coordinates of a block
       CoordMax_DB(:,:),        & ! Max gen. coordinates of a block
       CellSize_DB(:,:),        & ! Cell size in gen. coordinates
       CellFace_DB(:,:),        & ! Cell faces for Cartesian grids
       CellFace_DFB(:,:,:,:,:), & ! Cell faces for general grids
       CellVolume_B(:),         & ! Cell volume for Cartesian grids
       CellVolume_GB(:,:,:,:),  & ! Cell volume for general grids
       Xyz_DGB(:,:,:,:,:)         ! Cartesian cell centers coords

  ! Local variables

  logical :: DoInitializeGrid = .true.
  

contains
  !============================================================================
  subroutine init_grid(CoordMinIn_D, CoordMaxIn_D)

    real, intent(in):: CoordMinIn_D(nDim), CoordMaxIn_D(nDim)
    !-------------------------------------------------------------------------
    if(.not. DoInitializeGrid) RETURN

    DoInitializeGrid = .false.

    ! Make sure that the thickness is unity in the ignored dimensions
    CoordMin_D = -0.5
    CoordMax_D = +0.5
    CoordMin_D(1:nDim) = CoordMinIn_D
    CoordMax_D(1:nDim) = CoordMaxIn_D

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

    ! Create geometrical information for block iBlock on the local PE

    integer, intent(in):: iBlock

    character(len=*), parameter:: NameSub = 'create_grid_block'

    real :: PositionMin_D(MaxDim), PositionMax_D(MaxDim)
    integer :: iNode, i, j, k
    !----------------------------------------------------------------------
    iNode = iNode_B(iBlock)
    call get_tree_position(iNode, PositionMin_D, PositionMax_D)

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
       call CON_stop(NameSub//' non-Cartesian is not yet implemented')
    end if

  end subroutine create_grid_block

  !===========================================================================
  subroutine create_grid

    integer:: iBlock
    !------------------------------------------------------------------------
    do iBlock = 1, nBlock
       if(Unused_B(iBlock))CYCLE
       call create_grid_block(iBlock)
    end do

  end subroutine create_grid
  !===========================================================================

  subroutine show_grid_block(iBlock)

    use BATL_mpi, ONLY: iProc

    integer, intent(in):: iBlock

    ! Show grid information for block iBlock

    character(len=*), parameter:: NameSub = 'show_grid_block'
    !------------------------------------------------------------------------
    if(Unused_B(iBlock))then
       write(*,*) NameSub//' WARNING unused block ',iBlock,' on proc',iProc
       RETURN
    end if
    write(*,*)'show_grid_block for iProc, iBlock=',iProc, iBlock
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

  subroutine show_grid

    use ModUtilities, ONLY: flush_unit
    use ModIoUnit,    ONLY: STDOUT_
    use BATL_mpi, ONLY: iProc, nProc, barrier_mpi

    ! Show all blocks sequentially on all processors, ie. show_grid 
    ! must be called from all processors of the MPI communicator iComm!

    integer:: iBlock, iPe
    !------------------------------------------------------------------------

    do iPe = 0, nProc - 1
       if(iPe == iProc) then
          do iBlock = 1, nBlock
             if(Unused_B(iBlock)) CYCLE
             call show_grid_block(iBlock)
          end do
       end if
       call flush_unit(STDOUT_)
       call barrier_mpi
    end do
  end subroutine show_grid

  !===========================================================================

  subroutine test_grid

    use BATL_mpi, ONLY: iProc
    use BATL_geometry, ONLY: init_geometry

    integer :: iBlock, nBlockAll, Int_D(MaxDim)

    integer, parameter:: MaxBlockTest            = 50
    integer, parameter:: nRootTest_D(MaxDim)     = (/3,2,1/)
    logical, parameter:: IsPeriodicTest_D(MaxDim)= (/.true., .true., .false./)
    real:: DomainMin_D(MaxDim) = (/ 3.0, 2.0, 1.0 /)
    real:: DomainMax_D(MaxDim) = (/ 9.0, 6.0, 4.0 /)

    logical:: DoTestMe
    character(len=*), parameter :: NameSub = 'test_grid'
    !-----------------------------------------------------------------------
    DoTestMe = iProc == 0

    if(DoTestMe) write(*,*)'Testing init_grid'
    if(DoTestMe) write(*,*)'nDimAmr, nIJK_D=', nDimAmr, nIJK_D
    call init_tree(MaxBlockTest)
    call init_grid( DomainMin_D(1:nDim), DomainMax_D(1:nDim) )
    call init_geometry( IsPeriodicIn_D = IsPeriodicTest_D(1:nDim) )
    call set_tree_root( nRootTest_D(1:nDim))

    call refine_tree_node(3)
    call distribute_tree(.true.)
    if(DoTestMe) call show_tree('After distribute_tree')

    if(DoTestMe) write(*,*)'Testing create_grid'
    call create_grid

    call show_grid

    if(DoTestMe) write(*,*)'Testing clean_grid'
    call clean_grid
    
  end subroutine test_grid

end module BATL_grid
