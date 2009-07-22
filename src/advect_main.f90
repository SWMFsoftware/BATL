program advect

  use BATL_lib, ONLY: nDim, nI, nJ, nK, &
       MinI, MaxI, MinJ, MaxJ, MinK, MaxK, iProc, barrier_mpi

  implicit none

  ! Square of the radius of the sphere
  real, parameter :: Radius2 = 25.0

  ! Final simulation time, frequency of plots
  real, parameter :: TimeMax = 10.0, DtPlot = 1.0

  ! Advection velocity. Should be positive. For now set to 2
  real :: Velocity_D(nDim) = 2.0

  ! Spatial order of accuracy and beta parameter for the TVD limiter
  integer, parameter :: nOrder = 2
  real,    parameter :: BetaLimiter = 1.5 ! 1 <= Beta <= 2 for nOrder=2

  ! Parameters for the explicit scheme
  real, parameter :: Cfl = 0.8

  ! Size of the computational domain
  real, parameter :: &
       DomainMax_D(3) = (/ 10.0, 5.0, 5.0 /), &
       DomainMin_D(3) = -DomainMax_D

  ! Time step counter and simulation time
  integer :: iStep
  real :: Time, TimePlot, Dt

  integer, parameter:: nVar = 1, Rho_ = 1

  ! Cell centered state
  real, allocatable :: State_VGB(:,:,:,:,:), StateOld_VCB(:,:,:,:,:)

  ! Face centered flux for one block
  real:: Flux_VFD(nVar,1:nI+1,1:nJ+1,1:nK+1,nDim)

  ! Total initial mass
  real:: TotalIni_I(0:nVar) = -1.0

  logical, parameter :: DoTest = .false.
  character(len=*), parameter:: NameSub = 'advect_main'
  !--------------------------------------------------------------------------
  call initialize

  if(DoTest)write(*,*)NameSub,' starting iProc=',iProc
  if(DoTest)call barrier_mpi

  call timing_start('ADVECT')
  do
     if(DoTest)write(*,*)NameSub,' advance iProc, iStep=',iProc,iStep

     ! Save plot at required frequency
     if( Time >= TimePlot - 1e-10 )then

        if(DoTest)write(*,*)NameSub,' saving plots iProc=',iProc

        call timing_start('save_plot')
        call save_plot
        call timing_stop('save_plot')

        if(DoTest)write(*,*)NameSub,' save logfile iProc=',iProc

        call timing_start('save_log')
        call save_log
        call timing_stop('save_log')

        if(DoTest)write(*,*)NameSub,' finished plots iProc=',iProc

        TimePlot = TimePlot + DtPlot
     end if

     if(DoTest)write(*,*)NameSub,' advance_expl iProc=',iProc
     if(DoTest)call barrier_mpi

     call timing_start('explicit')
     call advance_explicit
!!!     call advance_localstep !!!
     call timing_stop('explicit')

     if(DoTest)write(*,*)NameSub,' adapt grid iProc=',iProc
     if(DoTest)call barrier_mpi

     call timing_start('amr')
     call adapt_grid
     call timing_stop('amr')

     if(DoTest)write(*,*)NameSub,' update time iProc=',iProc
     if(DoTest)call barrier_mpi

     ! Update time
     iStep = iStep + 1
     Time  = Time + Dt

     call timing_step(iStep)

     if(Time >= TimeMax - 1e-10) EXIT

  end do
  call timing_stop('ADVECT')
  call timing_report_total

  if(DoTest)write(*,*)NameSub,' finalize iProc=',iProc

  call finalize

contains

  !===========================================================================
  subroutine adapt_grid

    use BATL_lib, ONLY: regrid_batl, nBlock, Unused_B, iNode_B, &
         iStatusNew_A, Refine_, Coarsen_

    integer:: iBlock
    !------------------------------------------------------------------------
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       if(maxval(State_VGB(Rho_,1:nI,1:nJ,1:nK,iBlock))>1.5)then !!! 1.1
          iStatusNew_A(iNode_B(iBlock)) = Refine_
       elseif(minval(State_VGB(Rho_,1:nI,1:nJ,1:nK,iBlock))<1.4)then !!! 1.05
          iStatusNew_A(iNode_B(iBlock)) = Coarsen_
       end if
    end do

    call regrid_batl(nVar,State_VGB)

  end subroutine adapt_grid
  !===========================================================================
  real function exact_density(Xyz_D)

    use ModNumConst, ONLY: cHalfPi

    real, intent(in):: Xyz_D(nDim)

    ! Square of the radius of the circle/sphere
    real:: XyzShift_D(nDim)

    real :: r2

    real, parameter:: DomainSize_D(nDim) = &
         DomainMax_D(1:nDim)-DomainMin_D(1:nDim)
    !-------------------------------------------------------------------------
    ! Move position back to initial point
    XyzShift_D = Xyz_D - Time*Velocity_D

    ! Take periodicity into account
    XyzShift_D = modulo(XyzShift_D - DomainMin_D(1:nDim), DomainSize_D) &
         + DomainMin_D(1:nDim)

    r2 = sum(XyzShift_D**2)
    if(r2 < Radius2)then
       exact_density = 1 + cos(cHalfPi*sqrt(r2/Radius2))**2
    else
       exact_density = 1
    end if

  end function exact_density

  !===========================================================================
  subroutine initialize

    use BATL_lib, ONLY: init_mpi, init_batl, init_grid_batl, &
         iProc, MaxBlock, nBlock, Unused_B, Xyz_DGB, iNode_B, &
         iTree_IA, iStatusNew_A, MaxLevel_, Refine_

    integer :: i, j, k, iBlock, iError
    !------------------------------------------------------------------------

    call init_mpi
    call init_batl( &
         MaxBlockIn     = 2000, &          
         CoordMinIn_D   = DomainMin_D, &
         CoordMaxIn_D   = DomainMax_D, &
         nRootIn_D      = (/4,4,2/),   & 
         IsPeriodicIn_D = (/.true.,.true.,.true./) )

    ! Allow only one level of refinement
    iTree_IA(MaxLevel_,:) = 1

    LOOPBLOCK: do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(sum(Xyz_DGB(:,i,j,k,iBlock)**2) < Radius2)then !!!
             iStatusNew_A(iNode_B(iBlock)) = Refine_
             CYCLE LOOPBLOCK
          end if
       end do; end do; end do
    end do LOOPBLOCK

    call init_grid_batl

    ! Initial time step and time
    iStep    = 0
    Time     = 0.0
    TimePlot = 0.0

    allocate( &
         State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock), &
         StateOld_VCB(nVar,nI,nJ,nK,MaxBlock) )

    State_VGB = 0.0
    Flux_VFD = 0.0

    ! Initialize the state as a sphere with a cos^2 density profile
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
          State_VGB(:,i,j,k,iBlock) &
               = exact_density(Xyz_DGB(1:nDim,i,j,k,iBlock))
       end do; end do; end do
    end do

    if(iProc==0)then
       call timing_active(.true.)
       call timing_step(0)
       call timing_depth(-1)
       call timing_report_style('tree')
    end if

  end subroutine initialize

  !===========================================================================
  subroutine save_log

    ! Calculate the totals on processor 0
    use BATL_lib, ONLY: nDimAmr, nBlock, Unused_B, CellVolume_B, Xyz_DGB, &
         iComm, iProc, nProc
    use ModMpi,    ONLY: MPI_reduce, MPI_REAL, MPI_SUM
    use ModIoUnit, ONLY: UnitTmp_

    integer:: iBlock, i, j, k, iVar, iError
    integer, parameter:: Volume_=nVar+1, Error_=nVar+2
    real:: TotalPe_I(Error_), Total_I(Error_), Error

    character(len=100):: NameFile = '???'
    logical :: DoInitialize = .true.
    !------------------------------------------------------------------------
    Total_I = 0.0
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       do iVar = 1, nVar
          Total_I(iVar) = Total_I(iVar) + CellVolume_B(iBlock) &
               *sum(State_VGB(iVar,1:nI,1:nJ,1:nK,iBlock))
       end do
       Total_I(Volume_) = Total_I(Volume_) &
            + CellVolume_B(iBlock)*nI*nJ*nK

       do k=1,nK; do j=1,nJ; do i=1,nI
          Total_I(Error_) = Total_I(Error_) &
               + CellVolume_B(iBlock) &
               *abs(State_VGB(Rho_,i,j,k,iBlock)  &
               -    exact_density(Xyz_DGB(1:nDim,i,j,k,iBlock)))
       end do; end do; end do
    end do

    if(nProc > 1)then
       TotalPe_I = Total_I
       call MPI_reduce(TotalPe_I, Total_I, Error_, MPI_REAL, MPI_SUM, 0, &
            iComm, iError)
    end if

    if(iProc /= 0) RETURN

    ! Divide total error by total mass to get relative error
    Total_I(Error_) = Total_I(Error_) / Total_I(1)

    if(DoInitialize)then
       write(NameFile,'(a,i1,i1,a)') 'advect',nDim,nDimAmr,'.log'
       open(UnitTmp_,file=NameFile,status='replace')
       write(UnitTmp_,'(a)')'Advection test for BATL'
       write(UnitTmp_,'(a)')'mass volume error'

       DoInitialize = .false.
    else
       open(UnitTmp_,file=NameFile,position='append')
    end if
    write(UnitTmp_,'(100es15.6)') Total_I
    close(UnitTmp_)

  end subroutine save_log


  !===========================================================================

  subroutine save_plot

    use BATL_lib, ONLY: MaxDim, nDimAmr, nIJK, nBlock, Unused_B, &
         iComm, nProc, iProc, iNode_B, &
         TypeGeometry, CellSize_DB, Xyz_DGB, CoordMin_DB, CoordMax_DB

    use ModMpi,    ONLY: MPI_REAL, MPI_INTEGER, MPI_MIN, MPI_SUM, MPI_reduce
    use ModIoUnit, ONLY: UnitTmp_
    use ModKind,   ONLY: nByteReal

    character(len=100):: NameSnapshot, NameFile
    real:: CellSizeMin_D(MaxDim), CellSizeMinAll_D(MaxDim), &
         CellSizePlot_D(MaxDim)
    real:: PlotMax_D(MaxDim), PlotMin_D(MaxDim)
    integer :: iDim, iBlock, i, j, k, iError
    integer :: nCell, nCellAll, nPlotDim
    !-----------------------------------------------------------------------
    write(NameSnapshot,'(a,i7.7)') 'plots/cut_var_1_n',iStep

    ! write data from all processors into separate files
    write(NameFile,'(a,i4.4,a)') trim(NameSnapshot)//'_pe',iProc,'.idl'
    open(UnitTmp_, file=NameFile, status='replace', form='unformatted')
    nCell = 0
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       if(CoordMin_DB(3,iBlock) > 1e-6) CYCLE
       if(CoordMax_DB(3,iBlock) <-1e-6) CYCLE
       do k = 1, nK; 
          if(abs(Xyz_DGB(3,1,1,k,iBlock)) > 0.51*CellSize_DB(3,iBlock)) CYCLE
          do j = 1, nJ; do i = 1, nI
             nCell = nCell + 1
             write(UnitTmp_) CellSize_DB(1,iBlock), &
                  Xyz_DGB(:,i,j,k,iBlock), State_VGB(:,i,j,k,iBlock), &
                  exact_density(Xyz_DGB(:,i,j,k,iBlock)), &
                  CellSize_DB(1,iBlock), real(iNode_B(iBlock)), &
                  real(iProc), real(iBlock)
          end do; end do; 
       end do
    end do
    close(UnitTmp_)

    do iDim = 1, MaxDim
       CellSizeMin_D(iDim) = &
            minval(CellSize_DB(iDim,1:nBlock), MASK=.not.Unused_B(1:nBlock))
    end do
    call MPI_reduce(CellSizeMin_D, CellSizeMinAll_D, MaxDim, MPI_REAL, &
         MPI_MIN, 0, iComm, iError)

    call MPI_reduce(nCell, nCellAll, 1, MPI_INTEGER, MPI_SUM, 0, iComm, iError)

    if(iProc == 0)then
       nPlotDim = min(2,nDim)
       CellSizePlot_D = CellSizeMinAll_D
       CellSizePlot_D(1:nDimAmr) = -1.0

       PlotMin_D = -1e-10; PlotMin_D(1:nPlotDim) = DomainMin_D(1:nPlotDim)
       PlotMax_D = +1e-10; PlotMax_D(1:nPlotDim) = DomainMax_D(1:nPlotDim)

       NameFile = trim(NameSnapshot)//'.h'
       open(UnitTmp_,file=NameFile,status="replace")
       write(UnitTmp_,'(a)')         NameFile
       write(UnitTmp_,'(i8,a)')      nProc,        ' nProc'
       write(UnitTmp_,'(i8,a)')      iStep,        ' n_step'
       write(UnitTmp_,'(1pe13.5,a)') Time,         ' t'
       write(UnitTmp_,'(6(1pe18.10),a)') &
            (PlotMin_D(iDim),PlotMax_D(iDim),iDim=1,MaxDim),' plot_range'
       write(UnitTmp_,'(6(1pe18.10),i10,a)') &
            CellSizePlot_D, &
            CellSizeMinAll_D, nCellAll,            ' plot_dx, dxmin, ncell'
       write(UnitTmp_,'(i8,a)')     nVar+5,        ' nplotvar'
       write(UnitTmp_,'(i8,a)')     1,             ' neqpar'
       write(UnitTmp_,'(10es13.5)') 0.0            ! eqpar
       write(UnitTmp_,'(a)')        'rho exact dx node proc block none' 
       write(UnitTmp_,'(a)')        '1 1 1'        ! units
       write(UnitTmp_,'(l8,a)')     .true.         ! save binary .idl files
       write(UnitTmp_,'(i8,a)')     nByteReal,     ' nByteReal'
       write(UnitTmp_,'(a)')        TypeGeometry
       write(UnitTmp_,'(a)')        'real4'        ! type of .out file
       close(UnitTmp_)
    end if

  end subroutine save_plot
  !===========================================================================
  subroutine save_plot_block

    use ModPlotFile, ONLY: save_plot_file
    use BATL_lib, ONLY: iProc, &
         nBlock, Unused_B, CoordMin_DB, CoordMax_DB, CellSize_DB
    
    integer :: iBlock
    character(len=100):: NameFile
    character (len=10) :: TypePosition = 'rewind'
    !---------------------------------------------------------------------

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE

       write(NameFile,'(a,i3.3,a,i5.5,a)') &
            'advect_pe',iProc,'_blk',iBlock,'.out'
       
       call save_plot_file(NameFile,     &
            TypeFileIn='real4',          &
            TypePositionIn=TypePosition, &
            nStepIn = iStep, &
            TimeIn  = Time, &
            nDimIn  = nDim, &
            CoordMinIn_D = CoordMin_DB(1:nDim,iBlock)        &
            +              0.5*CellSize_DB(1:nDim,iBlock),   &
            CoordMaxIn_D = CoordMax_DB(1:nDim,iBlock)        &
            -              0.5*CellSize_DB(1:nDim,iBlock),   &
            VarIn_VIII = State_VGB(:,1:nI,1:nJ,1:nK,iBlock))
    end do

    TypePosition = 'append'

  end subroutine save_plot_block
  !===========================================================================
  subroutine finalize

    use BATL_lib, ONLY: clean_batl, clean_mpi
    !------------------------------------------------------------------------
    call save_plot
    call save_log

    call clean_batl
    call clean_mpi

  end subroutine finalize
  !===========================================================================
  subroutine calc_face_values(iBlock)
    use ModNumConst, ONLY: i_DD

    integer, intent(in):: iBlock

    real:: StateLeft_VFD( nVar,1:nI+1,1:nJ+1,1:nK+1,nDim)
    real:: StateRight_VFD(nVar,1:nI+1,1:nJ+1,1:nK+1,nDim)
    real:: Slope_VGD(nVar,0:nI+1,0:nJ+1,0:nK+1,nDim)

    integer :: iDim, i, j, k, Di, Dj, Dk
    !------------------------------------------------------------------------

    if(nOrder==2)call limit_slope(State_VGB(:,:,:,:,iBlock), Slope_VGD)

    do iDim = 1, nDim
       Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
       if(nOrder==1)then
          do k = 1, nK+Dk; do j =1, nJ+Dj; do i = 1, nI+Di
             StateLeft_VFD( :,i,j,k,iDim) = State_VGB(:,i-Di,j-Dj,k-Dk,iBlock)
             StateRight_VFD(:,i,j,k,iDim) = State_VGB(:,i,j,k,iBlock)
          end do; end do; end do
       else
          do k = 1, nK+Dk; do j =1, nJ+Dj; do i = 1, nI+Di
             StateLeft_VFD( :,i,j,k,iDim) = State_VGB(:,i-Di,j-Dj,k-Dk,iBlock)&
                  + 0.5*Slope_VGD(:,i-Di,j-Dj,k-Dk,iDim)
             StateRight_VFD(:,i,j,k,iDim) = State_VGB(:,i,j,k,iBlock) &
                  - 0.5*Slope_VGD(:,i,j,k,iDim)
          end do; end do; end do
          
       end if
    end do

    ! Calculate upwinded fluxes (assume positive velocities)
    do iDim = 1, nDim
       Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
       do k = 1, nK+Dk; do j =1, nJ+Dj; do i = 1, nI+Di
          Flux_VFD(:,i,j,k,iDim) = Velocity_D(iDim)*StateLeft_VFD(:,i,j,k,iDim)
       end do; end do; end do
    end do

  end subroutine calc_face_values

  !===========================================================================
  subroutine limit_slope(State_VG, Slope_VGD)

    ! Calculate TVD limited slopes Slope_GD of State_VG

    use ModNumConst, ONLY: i_DD

    real, intent(in) :: State_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real, intent(out):: Slope_VGD(nVar,0:nI+1,0:nJ+1,0:nK+1,nDim)

    real    :: SlopeLeft_V(nVar), SlopeRight_V(nVar)
    integer :: iDim, i, j, k, Di, Dj, Dk
    !----------------------------------------------------------------------
    do iDim = 1, nDim
       Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
       do k = 1-Dk, nK+Dk; do j = 1-Dj, nJ+Dj; do i = 1-Di, nI+Di
          SlopeLeft_V  = State_VG(:,i,j,k)  - State_VG(:,i-Di,j-Dj,k-Dk)
          SlopeRight_V = State_VG(:,i+Di,j+Dj,k+Dk) - State_VG(:,i,j,k)
          Slope_VGD(:,i,j,k,iDim) = &
               (sign(0.5, SlopeLeft_V) + sign(0.5, SlopeRight_V))*min( &
               BetaLimiter*abs(SlopeLeft_V), &
               BetaLimiter*abs(SlopeRight_V), &
               0.5*abs(SlopeLeft_V + SlopeRight_V))
       end do; end do; end do
    end do

  end subroutine limit_slope
  !===========================================================================
  subroutine advance_explicit

    use BATL_lib, ONLY: message_pass_cell, nBlock, Unused_B, CellSize_DB, &
         iComm, nProc
    use ModNumConst, ONLY: i_DD
    use ModMpi

    integer:: iStage, iDim, iBlock, i, j, k, Di, Dj, Dk, iError
    real:: DtPe

    logical, parameter:: DoTest = .false.
    character(len=*), parameter:: NameSub = 'advance_explicit'
    !--------------------------------------------------------------------------
    Dt = 1e30
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       Dt = min( Dt, 1 / sum( Velocity_D/CellSize_DB(1:nDim,iBlock) ) ) 
    end do
    Dt = min( Cfl*Dt, TimeMax - Time)

    if(nProc > 1)then
       DtPe = Dt
       call MPI_allreduce(DtPe, Dt, 1, MPI_REAL, MPI_MIN, iComm, iError)
    end if

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       StateOld_VCB(:,:,:,:,iBlock) = State_VGB(:,1:nI,1:nJ,1:nK,iBlock)
    end do

    do iStage = 1, nOrder

       if(DoTest)write(*,*)NameSub,' advance_expl iProc, iStage=',iProc,iStage       

       call timing_start('message_pass')
       call message_pass_cell(nVar, State_VGB, &
            DoSendCornerIn=.true., nProlongOrderIn=2)
       call timing_stop('message_pass')

       if(DoTest)write(*,*)NameSub,' finished message_pass iProc=',iProc

       call timing_start('update')
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE

          call calc_face_values(iBlock)

          State_VGB(:,1:nI,1:nJ,1:nK,iBlock) = StateOld_VCB(:,:,:,:,iBlock)

          do iDim = 1, nDim
             Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                State_VGB(:,i,j,k,iBlock) = State_VGB(:,i,j,k,iBlock) &
                     - (Dt*iStage)/nOrder / CellSize_DB(iDim, iBlock) &
                     *(Flux_VFD(:,i+Di,j+Dj,k+Dk,iDim) -Flux_VFD(:,i,j,k,iDim))
             end do; end do; end do

          end do
       end do
       call timing_stop('update')

    end do

  end subroutine advance_explicit
  !===========================================================================
  subroutine advance_localstep

    use BATL_lib, ONLY: message_pass_cell, MaxBlock, nBlock, Unused_B, &
         CellSize_DB, iComm, nProc, DiLevelNei_IIIB
    use ModNumConst, ONLY: i_DD
    use ModMpi

    integer:: nStage, iStage, iStageBlock, iDim, iBlock, i, j, k, Di, Dj, Dk, iError
    real :: DtMin, DtMax, DtMinPe, DtMaxPe
    real :: TimeStage, Coef
    real, save, allocatable:: Dt_B(:), Time_B(:), TimeOld_B(:)
    integer, save, allocatable:: iStage_B(:)

    logical, parameter:: DoTest = .false.
    character(len=*), parameter:: NameSub = 'advance_localstep'
    !-----------------------------------------------------------------------
    if(.not.allocated(Dt_B))then
       allocate(Dt_B(MaxBlock), Time_B(MaxBlock), TimeOld_B(MaxBlock), &
            iStage_B(MaxBlock))
       Time_B    = Time
       TimeOld_B = Time
       iStage_B  = 1
    end if

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       Dt_B(iBlock) = Cfl / sum( Velocity_D/CellSize_DB(1:nDim,iBlock) )
       Time_B(iBlock) = Time
    end do
    DtMin = minval(Dt_B(1:nBlock), MASK = .not. Unused_B(1:nBlock))
    DtMax = maxval(Dt_B(1:nBlock), MASK = .not. Unused_B(1:nBlock))

    if(nProc > 1)then
       DtMinPe = DtMin
       call MPI_allreduce(DtMinPe, DtMin, 1, MPI_REAL, MPI_MIN, iComm, iError)
       DtMaxPe = DtMax
       call MPI_allreduce(DtMaxPe, DtMax, 1, MPI_REAL, MPI_MAX, iComm, iError)
    end if

    ! Exploit the fact that time step is proportional to cell size !!!
    Dt = DtMax
    nStage = nint(DtMax / DtMin)
    TimeStage = Time

    ! write(*,*)'!!! nStage = ',nStage

    do iStage = 1, 2*nStage

       call timing_start('message_pass')
       call message_pass_cell(nVar, State_VGB, &
            DoSendCornerIn=.true., nProlongOrderIn=2, &
            TimeOld_B=TimeOld_B, Time_B=Time_B)

       call timing_stop('message_pass')

       call timing_start('update')
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          if(Time_B(iBlock) > TimeStage*(1+1e-10)) CYCLE

          iStageBlock = iStage_B(iBlock)

          call calc_face_values(iBlock)

          if(iStageBlock == 1)then
             StateOld_VCB(:,:,:,:,iBlock) = State_VGB(:,1:nI,1:nJ,1:nK,iBlock)
          else
             State_VGB(:,1:nI,1:nJ,1:nK,iBlock) = StateOld_VCB(:,:,:,:,iBlock)
          end if

          do iDim = 1, nDim
             Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
             Coef =  Dt_B(iBlock) * iStageBlock/2.0 / CellSize_DB(iDim, iBlock)
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                State_VGB(:,i,j,k,iBlock) = State_VGB(:,i,j,k,iBlock) &
                     - Coef * &
                     (Flux_VFD(:,i+Di,j+Dj,k+Dk,iDim) - Flux_VFD(:,i,j,k,iDim))
             end do; end do; end do
          end do

          !write(*,*)'!!! iBlock, Time_B -> Time_B=', &
          !     iBlock, Time_B(iBlock), Time_B(iBlock) + Dt_B(iBlock)

          TimeOld_B(iBlock) = Time_B(iBlock)
          Time_B(iBlock)    = Time_B(iBlock) + Dt_B(iBlock)/2
          iStage_B(iBlock) = 3 - iStageBlock

       end do

       TimeStage = TimeStage + DtMin/2

       call timing_stop('update')
    end do

    if( maxval(Time_B, MASK=.not.Unused_B) > &
         minval(Time_B, MASK=.not.Unused_B) + 1e-10)then

       write(*,*)'!!! minval(Time_B)=',minval(Time_B, MASK=.not.Unused_B)
       write(*,*)'!!! maxval(Time_B)=',maxval(Time_B, MASK=.not.Unused_B)
       write(*,*)'!!! minloc(Time_B)=',minloc(Time_B, MASK=.not.Unused_B)
       write(*,*)'!!! maxloc(Time_B)=',maxloc(Time_B, MASK=.not.Unused_B)

       call CON_stop('time step problem')

    end if

  end subroutine advance_localstep

end program advect

!=============================================================================
subroutine CON_stop(String)
  use ModMpi, ONLY: MPI_abort, MPI_COMM_WORLD
  implicit none
  integer:: iError, nError
  character (len=*), intent(in) :: String
  !--------------------------------------------------------------------------
  write(*,*)'CON_stop called with String='
  write(*,*) String

  call MPI_abort(MPI_COMM_WORLD, nError, iError)
  stop
end subroutine CON_stop
