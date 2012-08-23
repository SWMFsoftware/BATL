module ModReadAmr

  ! reconstruct AMR grid and read data on this grid
  ! interpolate data

  use BATL_lib, ONLY: MaxDim

  implicit none
  save

  private !except

  ! Public methods 
  public:: readamr_init   ! read AMR tree information
  public:: readamr_read   ! read AMR data
  public:: readamr_get    ! get data at some point
  public:: readamr_clean  ! clean all variables

  ! Public data (optional)
  real, public, allocatable:: State_VGB(:,:,:,:,:)

  integer, public:: nVar=0                 ! number of variables in State_VGB

  character(len=20), public:: TypeGeometry = '???'
  real, public:: CoordMin_D(MaxDim) = -0.5 ! lower domain limits in gen coords
  real, public:: CoordMax_D(MaxDim) = +0.5 ! upper domain limits in gen coords

  real,    public:: TimeData = -1.0        ! simulation time of data
  integer, public:: nStepData = 0          ! time step of data
  integer, public:: nVarData = 0           ! number of variables in data file 
  integer, public:: nBlockData = 0         ! number of blocks in data file

  integer, public:: nParamData = 0         ! number of parameters in data file
  real, allocatable, public:: ParamData_I(:) ! paramters in data file

  character(len=500), public:: NameVarData  = '???' ! all variable names
  character(len=500), public:: NameUnitData = '???' ! unit names

  ! Local variables
  character(len=20):: TypeGeometryBatl = '???'
  integer:: nCellData = 0 ! number of cells in data file
  integer:: nProcData ! number of processors that wrote data file
  logical:: IsBinary  ! if the unprocessed IDL files are in binary format
  integer:: nByteReal ! number of bytes for reals in unprocessed IDL file

  character(len=5):: TypeDataFile ! type of processed file (real4/real8/ascii)

contains
  !============================================================================
  subroutine readamr_init(NameFile, IsVerboseIn)

    use ModIoUnit, ONLY: UnitTmp_
    use BATL_lib,  ONLY: MaxDim, nDim, nIJK, nProc, init_batl
    use BATL_grid, ONLY: create_grid
    use BATL_tree, ONLY: read_tree_file, distribute_tree

    character(len=*), intent(in):: NameFile
    logical, optional, intent(in):: IsVerboseIn ! provide verbose output

    logical:: IsVerbose

    integer:: i, iDim, iError

    character(len=500):: NameFileOrig

    integer:: MaxBlock
    integer:: nRgen=0
    real, allocatable:: Rgen_I(:)

    real:: CellSizePlot_D(MaxDim), CellSizeMin_D(MaxDim)

    ! Variables read from tree file
    integer:: nDimIn, nInfoIn, nNodeIn, iRatioIn_D(nDim),  nRoot_D(nDim)

    character(len=*), parameter:: NameSub = 'readamr_init'
    !-------------------------------------------------------------------------
    IsVerbose = .false.
    if(present(IsVerboseIn)) IsVerbose = IsVerboseIn

    open(UnitTmp_, file=trim(NameFile)//'.info', status='old', iostat=iError)
    if(iError /= 0) call CON_stop(NameSub// &
         ' ERROR: could not open '//trim(NameFile)//'.info')

    ! Read information from the .info file
    read(UnitTmp_,'(a)') NameFileOrig
    read(UnitTmp_,*) nProcData
    read(UnitTmp_,*) nStepData
    read(UnitTmp_,*) TimeData

    if(IsVerbose) write(*,*)'nStepData=', nStepData, ' TimeData=', TimeData

    read(UnitTmp_,*) (CoordMin_D(iDim), CoordMax_D(iDim), iDim=1,MaxDim)
    if(IsVerbose) write(*,*)'CoordMin_D=', CoordMin_D
    if(IsVerbose) write(*,*)'CoordMax_D=', CoordMax_D

    read(UnitTmp_,*) CellSizePlot_D, CellSizeMin_D, nCellData
    if(CellSizePlot_D(1) >= 0.0) call CON_stop(NameSub// &
         ': the resolution should be set to -1 for file'//trim(NameFile))
    if(IsVerbose) write(*,*)'nCellData=', nCellData

    ! Total number of blocks in the data file
    nBlockData = nCellData / nIJK

    ! Number of blocks per processor !!! this is wrong for box cut from sphere
    MaxBlock  = (nBlockData + nProc - 1) / nProc

    read(UnitTmp_,*) nVarData
    read(UnitTmp_,*) nParamData
    allocate(ParamData_I(nParamData))
    read(UnitTmp_,*) ParamData_I
    if(IsVerbose) write(*,*)'nVarData=',nVarData,' nParamData=', nParamData, &
         ' ParamData_I=', ParamData_I

    read(UnitTmp_,'(a)') NameVarData
    if(IsVerbose) write(*,*)'NameVarData =',trim(NameVarData)

    read(UnitTmp_,'(a)') NameUnitData
    if(IsVerbose)write(*,*) 'NameUnitData=', trim(NameUnitData)

    read(UnitTmp_,*) IsBinary
    if(IsBinary) read(UnitTmp_,*) nByteReal

    read(UnitTmp_,'(a)') TypeGeometry
    if(index(TypeGeometry,'genr') > 0)then
       read(UnitTmp_,*) nRgen
       allocate(Rgen_I(nRgen))
       do i = 1, nRgen
          read(UnitTmp_,*) Rgen_I(i)
       end do
       ! For some reason we store the logarithm of the radius, so take exp
       Rgen_I = exp(Rgen_I)
    end if

    read(UnitTmp_,'(a)') TypeDataFile  ! TypeFile for IDL data

    close(UnitTmp_)

    ! Read root block information from the tree file
    open(UnitTmp_, file=trim(NameFile)//'.tree', status='old', &
         form='unformatted', iostat=iError)
    if(iError /= 0) call CON_stop(NameSub// &
         ': could not open tree file'//trim(NameFile)//'.tree')
    read(UnitTmp_) nDimIn, nInfoIn, nNodeIn
    read(UnitTmp_) iRatioIn_D
    read(UnitTmp_) nRoot_D
    close(UnitTmp_)

    ! Initialize BATL (using generalized coordinates and radians)
    if(TypeGeometry(1:9)=='spherical') then
       TypeGeometryBatl = 'rlonlat'//TypeGeometry(10:20)
    else
       TypeGeometryBatl = TypeGeometry
    end if

    call init_batl(CoordMin_D, CoordMax_D, MaxBlock, &
         TypeGeometryBatl, rGenIn_I=rGen_I, nRootIn_D=nRoot_D, &
         UseRadiusIn=.false., UseDegreeIn=.false.)

    ! Read the full tree information and create grid
    call read_tree_file(trim(NameFile)//'.tree')
    call distribute_tree(.true.)
    call create_grid

  end subroutine readamr_init
  !============================================================================
  subroutine readamr_read(NameFile, iVarIn_I, UseXyzTest, UseSinTest)

    use BATL_lib, ONLY: nDim, &
         MinI, MaxI, MinJ, MaxJ, MinK, MaxK, MaxBlock, nG, iProc, &
         find_grid_block, message_pass_cell, xyz_to_coord
    use ModPlotFile, ONLY: read_plot_file
    use ModConst,    ONLY: cTwoPi

    character(len=*), intent(in):: NameFile     ! data file name
    integer, optional, intent(in):: iVarIn_I(:) ! index of variables to store
    logical, optional, intent(in):: UseXyzTest  ! store x,y,z into State
    logical, optional, intent(in):: UseSinTest  ! store sin(coord) into State

    ! Allocatable arrays for holding linear file data
    real, allocatable  :: State_VI(:,:), Xyz_DI(:,:)

    real:: Xyz_D(MaxDim) = 0.0, Coord_D(MaxDim)
    integer:: iCell, iCell_D(MaxDim), i, j, k, iBlock, iProcFound

    character(len=*), parameter:: NameSub = 'readamr_read'
    !--------------------------------------------------------------------------
    nVar = nVarData
    if(present(iVarIn_I))then
       nVar = size(iVarIn_I)
       if(nVar>nVarData .or. any(iVarIn_I<1) .or. any(iVarIn_I>nVarData))then
          write(*,*)'ERROR: nVarData, iVarIn_I=', nVarData, iVarIn_I
          call CON_stop(NameSub//': invalid iVarIn_I array')
       end if
    end if

    ! The tests need at least nDim variables to be stored
    if(present(UseXyzTest) .or. present(UseSinTest)) nVar = max(nVar, nDim)

    allocate(State_VI(nVarData,nCellData), Xyz_DI(nDim,nCellData))

    call read_plot_file(NameFile, TypeFileIn=TypeDataFile, &
         CoordOut_DI=Xyz_DI, VarOut_VI = State_VI)

    ! allocate block-AMR data structure
    allocate(&
         State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))

    State_VGB = 0.0

    ! find node containing point 
    do iCell = 1, nCellData
       ! find cell on the grid
       Xyz_D(1:nDim) = Xyz_DI(:,iCell)
       call find_grid_block(Xyz_D, iProcFound, iBlock, iCell_D)

       if(iBlock < 0)then
          write(*,*)'ERROR for iCell, Xyz_D=', iCell, Xyz_D
          call CON_stop(NameSub//': could not find cell on the grid')
       end if

       !check if point belongs on this processor
       if (iProcFound /= iProc) CYCLE

       i = iCell_D(1); j = iCell_D(2); k = iCell_D(3)
       if(present(iVarIn_I))then
          State_VGB(:,i,j,k,iBlock) = State_VI(iVarIn_I,iCell)
       else
          State_VGB(:,i,j,k,iBlock) = State_VI(:,iCell)
       end if

       ! For verification tests
       if(present(UseXyzTest))then
          ! Store x,y,z coordinates into first nDim elements
          State_VGB(1:nDim,i,j,k,iBlock) = Xyz_DI(:,iCell)
       elseif(present(UseSinTest))then
          call xyz_to_coord(Xyz_D, Coord_D)
          State_VGB(1:nDim,i,j,k,iBlock) = sin(cTwoPi*Coord_D(1:nDim) &
               /(CoordMax_D(1:nDim) - CoordMin_D(1:nDim)))
       end if
    enddo

    ! deallocate to save memory
    deallocate (State_VI, Xyz_DI) 

    ! Set ghost cells if any
    if(nG > 0) call message_pass_cell(nVar, State_VGB)

  end subroutine readamr_read

  !============================================================================
  subroutine readamr_get(Xyz_D, nVarIn, State_V, Weight, IsFound)

    use BATL_lib, ONLY: nDim, nG, iProc, &
         interpolate_grid, find_grid_block

    real,    intent(in)  :: Xyz_D(MaxDim)
    integer, intent(in)  :: nVarIn
    real,    intent(out) :: State_V(nVarIn)
    real,    intent(out) :: Weight
    logical, intent(out) :: IsFound

    ! Block and processor index for the point
    integer:: iBlock, iProcOut

    ! Variables for linear interpolation using ghost cells
    integer:: i1, j1=1, k1=1, i2, j2, k2
    real:: Dist_D(MaxDim), Dx1, Dx2, Dy1, Dy2, Dz1, Dz2

    ! Variables for AMR interpolation without ghost cells
    integer:: iCell, nCell, iCell_II(0:nDim,2**nDim), iCell_D(MaxDim), i, j, k
    real:: Weight_I(2**nDim)
    !-------------------------------------------------------------------------
    State_V = 0.0
    Weight  = 0.0
    if(nG > 0)then
       call find_grid_block(Xyz_D, iProcOut, iBlock, iCell_D, Dist_D)
       IsFound = iBlock > 0
       if(iProcOut /= iProc) RETURN
       Dx1 = Dist_D(1); Dx2 = 1 - Dx1
       i1  = iCell_D(1); i2 = i1 + 1
       if(nDim > 1)then
          Dy1 = Dist_D(2); Dy2 = 1 - Dy1
          j1 = iCell_D(2); j2 = j1 + 1
       end if
       if(nDim > 2)then
          Dz1 = Dist_D(3); Dz2 = 1 - Dz1
          k1 = iCell_D(3); k2 = k1 + 1
       end if
       if(nDim == 1)then
          State_V = Dx2*State_VGB(:,i1,j1,k1,iBlock)   &
               +    Dx1*State_VGB(:,i2,j1,k1,iBlock)
       end if
       if(nDim == 2)then
          State_V = Dy2*(Dx2*State_VGB(:,i1,j1,k1,iBlock)   &
               +         Dx1*State_VGB(:,i2,j1,k1,iBlock))  &
               +    Dy1*(Dx2*State_VGB(:,i1,j2,k1,iBlock)   &
               +         Dx1*State_VGB(:,i2,j2,k1,iBlock))
       end if
       if(nDim == 3)then
          State_V = Dz2*(Dy2*(Dx2*State_VGB(:,i1,j1,k1,iBlock)   &
               +              Dx1*State_VGB(:,i2,j1,k1,iBlock))  &
               +         Dy1*(Dx2*State_VGB(:,i1,j2,k1,iBlock)   &
               +              Dx1*State_VGB(:,i2,j2,k1,iBlock))) &
               +    Dz1*(Dy2*(Dx2*State_VGB(:,i1,j1,k2,iBlock)   &
               +              Dx1*State_VGB(:,i2,j1,k2,iBlock))  &
               +         Dy1*(Dx2*State_VGB(:,i1,j2,k2,iBlock)   &
               +              Dx1*State_VGB(:,i2,j2,k2,iBlock)))
       end if
       Weight = 1.0
    else
       ! Check if position is covered by blocks
       call find_grid_block(Xyz_D, iProcOut, iBlock)
       IsFound = iBlock > 0
       if(.not.IsFound) RETURN

       call interpolate_grid(Xyz_D, nCell, iCell_II, Weight_I)
       do iCell = 1, nCell
          iBlock  = iCell_II(0,iCell)
          iCell_D = 1
          iCell_D(1:nDim) = iCell_II(1:nDim,iCell)
          i      = iCell_D(1)
          j      = iCell_D(2)
          k      = iCell_D(3)
          Weight = Weight  + Weight_I(iCell)
          State_V= State_V + Weight_I(iCell)*State_VGB(:,i,j,k,iBlock)
       end do
    end if

  end subroutine readamr_get
  !============================================================================
  subroutine readamr_clean
    use BATL_lib, ONLY: clean_batl

    call clean_batl
    if(allocated(State_VGB)) deallocate(State_VGB)
    nVar       = 0
    nVarData   = 0
    nBlockData = 0

  end subroutine readamr_clean

end module ModReadAmr
