program flow2d

  use BATL_size

  implicit none

  ! Final simulation time, frequency of plots
  real, parameter :: TimeMax = 10.0, DtPlot = 1.0

  ! Advection velocity. Should be positive. For now set to 2
  real :: Velocity_D(nDim) = 2.0

  ! Name plot file
  character (len=*), parameter:: NamePlotFile='advect.out'

  ! Spatial order of accuracy and beta parameter for the TVD limiter
  integer, parameter :: nOrder = 2
  real,    parameter :: BetaLimiter = 1.5 ! 1 <= Beta <= 2 for nOrder=2

  ! Parameters for the explicit scheme
  real, parameter :: Cfl = 0.8

  ! Size of the computational domain
  real, parameter :: &
       DomainMax_D(3) = (/ 10.0, 5.0, 2.0 /), &
       DomainMin_D(3) = -DomainMax_D

  ! Time step counter and simulation time
  integer :: iStep
  real :: Time, TimePlot, Dt

  integer, parameter:: nVar = 1

  ! Cell centered state
  real, allocatable :: State_VGB(:,:,:,:,:), StateOld_VCB(:,:,:,:,:)

  ! Face centered flux for one block
  real:: Flux_VFD(nVar,1:nI+1,1:nJ+1,1:nK+1,nDim)
  !--------------------------------------------------------------------------
  call initialize
  do
     ! Save plot at required frequency
     if( Time >= TimePlot - 1e-10 )then
        call save_plot
        TimePlot = TimePlot + DtPlot
     end if

     call advance_explicit

     ! Update time
     iStep = iStep + 1
     Time  = Time + Dt

     if(Time >= TimeMax - 1e-10) EXIT

  end do
  call finalize

contains

  !===========================================================================
  subroutine initialize

    use ModMpi,    ONLY: MPI_COMM_WORLD, MPI_init
    use BATL_mpi,  ONLY: init_mpi
    use BATL_tree, ONLY: init_tree, Unused_B, set_tree_root, distribute_tree
    use BATL_grid, ONLY: init_grid, create_grid_block, Xyz_DGB
    use ModNumConst, ONLY: cHalfPi

    ! Square of the radius of the circle/sphere
    real, parameter :: Radius2 = 25.0
    real:: r2
    integer :: i, j, k, iBlock, iError
    !------------------------------------------------------------------------

    call MPI_init(iError)
    call init_mpi(MPI_COMM_WORLD)
    call init_tree( &
         MaxBlockIn = 10, &
         MaxNodeIn  = 100)
    call init_grid( &
         CoordMinIn_D = DomainMin_D, &
         CoordMaxIn_D = DomainMax_D )
    call set_tree_root( &
         nRootIn_D      = (/1,1,1/), &
         IsPeriodicIn_D = (/.true.,.true.,.true./) )

    call distribute_tree(.true.)
    do iBlock = 1, nBlock
       if(Unused_B(iBlock))CYCLE
       call create_grid_block(iBlock)
    end do

    allocate( &
         State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock), &
         StateOld_VCB(nVar,nI,nJ,nK,MaxBlock) )

    State_VGB = 0.0
    Flux_VFD = 0.0

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
          r2 = sum(Xyz_DGB(1:nDimTree,i,j,k,iBlock)**2)
          if(r2 < Radius2)then
             State_VGB(:,i,j,k,iBlock) &
                  = 1.0 + cos( cHalfPi*sqrt(r2/Radius2) )**2
          else
             State_VGB(:,i,j,k,iBlock) = 1.0
          end if
       end do; end do; end do
    end do

    write(*,*)'!!! maxval(State_VGB)=',maxval(State_VGB(:,:,:,:,1:nBlock))
    write(*,*)'!!! minval(State_VGB)=',minval(State_VGB(:,:,:,:,1:nBlock))

    ! Initial time step and time
    iStep    = 0
    Time     = 0.0
    TimePlot = 0.0

  end subroutine initialize
  !=======================================================================
  subroutine save_plot

    use ModPlotFile, ONLY: save_plot_file
    
    integer, parameter :: nPlotVar=3
    integer, parameter :: nEqpar=1
    integer, parameter :: nDim=2
    character (len=*), parameter :: &
         StringHead = "water_rho11", NameCoord = "x y z "

    character (len=10) :: TypePosition = 'rewind'
    !---------------------------------------------------------------------

    call save_plot_file(NamePlotFile, &
         TypeFileIn='real8',          &
         TypePositionIn=TypePosition, &
         nStepIn = iStep, &
         TimeIn  = Time, &
         nDimIn  = nDim, &
         CoordMinIn_D = DomainMin_D(1:nDim), &
         CoordMaxIn_D = DomainMax_D(1:nDim), &
         VarIn_VIII = State_VGB(:,1:nI,1:nJ,1:nK,1))

    TypePosition = 'append'

  end subroutine save_plot
  !===========================================================================
  subroutine finalize
    use BATL_tree, ONLY: clean_tree
    use BATL_grid, ONLY: clean_grid
    use ModMpi, ONLY: MPI_finalize
    integer :: iError
    !------------------------------------------------------------------------
    call save_plot
    call clean_grid
    call clean_tree
    call MPI_finalize(iError)
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

    ! Apply boundary conditions for h here: 1 block periodic !!!
    State_VGB(:,  -1,:,:,iBlock) = State_VGB(:,nI-1,:,:,iBlock)
    State_VGB(:,   0,:,:,iBlock) = State_VGB(:,nI  ,:,:,iBlock)
    State_VGB(:,nI+1,:,:,iBlock) = State_VGB(:,   1,:,:,iBlock)
    State_VGB(:,nI+2,:,:,iBlock) = State_VGB(:,   2,:,:,iBlock)
    State_VGB(:,:,  -1,:,iBlock) = State_VGB(:,:,nJ-1,:,iBlock)
    State_VGB(:,:,   0,:,iBlock) = State_VGB(:,:,nJ  ,:,iBlock)
    State_VGB(:,:,nJ+1,:,iBlock) = State_VGB(:,:,   1,:,iBlock)
    State_VGB(:,:,nJ+2,:,iBlock) = State_VGB(:,:,   2,:,iBlock)
    State_VGB(:,:,:,  -1,iBlock) = State_VGB(:,:,:,nK-1,iBlock)
    State_VGB(:,:,:,   0,iBlock) = State_VGB(:,:,:,nK  ,iBlock)
    State_VGB(:,:,:,nK+1,iBlock) = State_VGB(:,:,:,   1,iBlock)
    State_VGB(:,:,:,nK+2,iBlock) = State_VGB(:,:,:,   2,iBlock)

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

    use BATL_grid, ONLY: CellSize_DB
    use BATL_tree, ONLY: Unused_B
    use ModNumConst, ONLY: i_DD

    integer:: iStage, iDim, iBlock, i, j, k, Di, Dj, Dk
    !--------------------------------------------------------------------------
    Dt = 1e30
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       Dt = min( Dt, 1 / sum( Velocity_D/CellSize_DB(1:nDim,iBlock) ) ) 
    end do
    Dt = min( Cfl*Dt, TimeMax - Time)

    !!! MPI_allreduce

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       StateOld_VCB(:,:,:,:,iBlock) = State_VGB(:,1:nI,1:nJ,1:nK,iBlock)
    end do

    do iStage = 1, nOrder

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

    end do

  end subroutine advance_explicit

end program flow2d

!=============================================================================
subroutine CON_stop(String)
  use BATL_mpi
  implicit none
  character (len=*), intent(in) :: String
  integer:: iError, nError
  !--------------------------------------------------------------------------
  write(*,*)'CON_stop called with String='
  write(*,*) String
  call MPI_abort(iComm, nError, iError)
  stop
end subroutine CON_stop
