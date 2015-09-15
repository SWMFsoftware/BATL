
! Extrernal subroutines to call from other languages
!==============================================================================
subroutine wrapreadamr_set_mpi_param(iProcIn, nProcIn, iCommIn) bind(C)
  use BATL_lib,  ONLY: iProc, nProc, iComm

  ! Set the MPI parameters
  
  implicit none
  integer,intent(in):: iProcIn, nProcIn, iCommIn
  
  iProc = iProcIn
  nProc = nProcIn
  iComm = iCommIn

end subroutine wrapreadamr_set_mpi_param

!==============================================================================
subroutine wrapreadamr_get_nvar(nVarOut) bind(C)
  use ModReadAmr, ONLY: nVar

  implicit none
  integer,intent(out):: nVarOut
  
  nVarOut = nVar 
end subroutine wrapreadamr_get_nvar    
!==============================================================================
subroutine wrapreadamr_get_namevardata(l, NameVarOut) bind(C)

  ! Names of variables that are returned by the interpolation routine 

  use ModReadAmr, ONLY:NameVarData

  use iso_c_binding, ONLY: c_char
  
  implicit none
  integer,          intent(in) :: l
  character(kind=c_char), intent(out):: NameVarOut
 
  NameVarOut = NameVarData(1:l)
  write(*,*) "Fort: ", NameVarData(1:l)

end subroutine wrapreadamr_get_namevardata  
!==============================================================================
subroutine wrapreadamr_varnames(NameVarOut)

  ! Names of variables that are returned by the interpolation routine 

  use ModReadAmr, ONLY:NameVarData

  implicit none
  character(len=*), intent(out):: NameVarOut
  NameVarOut = NameVarData

end subroutine wrapreadamr_varnames
!==============================================================================
subroutine wrapreadamr_get_nameunitdata(l, NameUnitOut) bind(C)

  ! units of the coordinates followed by the units of the interpolated
  ! variables
  use ModReadAmr, ONLY: NameUnitData

  use iso_c_binding, ONLY: c_char
  
  implicit none
  integer, intent(in):: l
  character(kind=c_char), intent(out):: NameUnitOut
 
  NameUnitOut = NameUnitData(1:l)

end subroutine wrapreadamr_get_nameunitdata  
!==============================================================================
subroutine wrapreadamr_units(NameUnitOut)

  use ModReadAmr, ONLY: NameUnitData

  
  implicit none
  character(len=*), intent(out):: NameUnitOut
 
  NameUnitOut = NameUnitData

end subroutine wrapreadamr_units 
!==============================================================================
subroutine wrapreadamr_domain_limits(CoordMinOut_D, CoordMaxOut_D) bind(C)

  ! Boundary of the computational domain in the normalized units

  use ModKind,    ONLY: Real8_
  use ModReadAmr, ONLY: CoordMin_D, CoordMax_D
  
  implicit none
  real(Real8_), intent(out):: CoordMinOut_D(3), CoordMaxOut_D(3)

  CoordMinOut_D = CoordMin_D
  CoordMaxOut_D = CoordMax_D

end subroutine wrapreadamr_domain_limits  
!==============================================================================
subroutine wrapreadamr_read_file(l, FileName, lNewGrid, lVerbose) bind(C)

  ! Read data from file Filename.
  ! If lNewGrid == 1 then read all information, otherwise data only
  ! If lVerbose == 1 then write out verbose information

  use ModReadAmr, ONLY:readamr_read

  use iso_c_binding, ONLY: c_char

  implicit none

  integer,                intent(in):: l
  character(kind=c_char), intent(in):: FileName
  integer,                intent(in):: lNewGrid
  integer,                intent(in):: lVerbose

  call readamr_read(FileName, &
       IsNewGridIn = lNewGrid==1, IsVerboseIn=lVerbose==1)

end subroutine wrapreadamr_read_file
!==============================================================================

subroutine wrapreadamr_read_file_py(IsNewGrid, IsVerbose, NameFile)

  ! Read data from file Filename. Python interface.
  ! Call has to pass the length of the file name as an extra c_long argument.

  use ModReadAmr, ONLY:readamr_read

  implicit none

  logical,          intent(in):: IsNewGrid
  logical,          intent(in):: IsVerbose
  character(len=*), intent(in):: NameFile  ! base name

  call readamr_read(NameFile, IsVerboseIn=IsVerbose, IsNewGridIn=IsNewGrid)

end subroutine wrapreadamr_read_file_py

!=============================================================================
subroutine wrapreadamr_read_header(l, FileName) bind(C)

  ! Read header information from .h or .info
  ! This sets number and names of variables, domain size, etc.
  use ModReadAmr, ONLY:readamr_init

  use iso_c_binding, ONLY: c_char
  
  implicit none
  integer,intent(in):: l
  character(kind=c_char), intent(in):: FileName
  integer:: i
  !---------------------------------------------------------------------------
  ! Cut off extension from the file name (if any)
  i = index(FileName,".",BACK=.true.)
  call readamr_init(FileName(1:i-1), .true.)

end subroutine wrapreadamr_read_header
  
!==============================================================================
subroutine wrapreadamr_deallocate() bind(C)

  ! Deallocate all memory used by READAMR

  use ModReadAmr, ONLY: readamr_clean
  
  call readamr_clean
  
end subroutine wrapreadamr_deallocate

!=============================================================================
subroutine wrapreadamr_get_data(XyzIn_D, StateOut_V, iFound) bind(C)

  ! Get the interpolated values StateOut_V at the point XyzOut_V
  ! The first index of State_V is the interpolation weight, 
  ! so StateOut_V has nVar+1 elements.
  ! For parallel execution, an MPI_SUM is needed and a division by total weight.
  ! iFound is set to 0 if point is not found (outside domain), 1 otherwise.

  use ModReadAmr, ONLY: nVar, readamr_get
  use BATL_lib,   ONLY: MaxDim
  use ModKind,    ONLY: Real8_
  
  implicit none
  real(Real8_), intent(in) :: XyzIn_D(MaxDim)
  real(Real8_), intent(out):: StateOut_V(0:nVar)
  integer,intent(out):: iFound

  real   :: State_V(0:nVar)
  real   :: Xyz_D(MaxDim)
  logical:: IsFound
  !----------------------------------------------------------------------------

  ! This copy converts real precision if needed
  Xyz_D = XyzIn_D   
  call readamr_get(Xyz_D, State_V, IsFound)
  
  ! This copy converts real precision if needed
  StateOut_V = State_V

  ! Set integer found flag
  iFound = 0
  if(IsFound) iFound = 1

end subroutine wrapreadamr_get_data

!=============================================================================
subroutine wrapreadamr_get_more_data(X, S, iFound, nPoints) bind(C)
  
  ! interpolation loop over nPoints pairs of coordinates.
  ! X = array containing coordinates to interpolate at. 
  ! S = array containing the variables from interpolation.
  ! nPoints = number of coordinate pairs.


  use ModReadAmr, ONLY: nVar, readamr_get
  use BATL_lib,   ONLY: MaxDim
  use ModKind,    ONLY: Real8_
  
  implicit none
  real(Real8_), intent(in) :: X(MaxDim,nPoints)
  real(Real8_), intent(out):: S(1:nVar,nPoints)
  integer,intent(out)      :: iFound
  integer,intent(in)       :: nPoints

  integer:: i 
  real   :: State_V(0:nVar)
  real   :: Xyz_D(MaxDim)
  logical:: IsFound
  !----------------------------------------------------------------------------
  DO i=1,nPoints
      Xyz_D(1:MaxDim) = X(1:MaxDim,i)   
      call readamr_get(Xyz_D, State_V, IsFound)
      S(1:nVar,i) = State_V(1:nVar) / State_V(0)
  END DO

end subroutine wrapreadamr_get_more_data
!=============================================================================
subroutine wrapreadamr_get_data_serial(XyzIn_D, StateOut_V, iFound) bind(C)

  ! Get the interpolated values StateOut_V at the point XyzOut_V.
  ! Division by sum of weights is done internally, 
  ! so StateOut_V has nVar elements.
  ! iFound is set to 0 if point is not found (outside domain), 1 otherwise.

  use ModReadAmr, ONLY: nVar, readamr_get
  use BATL_lib,   ONLY: MaxDim
  use ModKind,    ONLY: Real8_
  
  implicit none
  real(Real8_), intent(in) :: XyzIn_D(MaxDim)
  real(Real8_), intent(out):: StateOut_V(nVar)
  integer,intent(out):: iFound

  real   :: State_V(0:nVar)
  real   :: Xyz_D(MaxDim)
  logical:: IsFound
  !----------------------------------------------------------------------------

  ! This copy converts real precision if needed
  Xyz_D = XyzIn_D   
  call readamr_get(Xyz_D, State_V, IsFound)

  ! Divide by weight.
  ! Also this converts real precision if needed.
  StateOut_V = State_V(1:nVar)/State_V(0)

  ! Set integer found flag
  iFound = 0
  if(IsFound) iFound = 1
  write(*,*) "get_data_serial called"

end subroutine wrapreadamr_get_data_serial

!=============================================================================
subroutine wrapreadamr_get_data_cell(XyzIn_D, StateOut_V, &
     CellSizeOut_D, iFound) bind(C)


  ! Get the interpolated values StateOut_V at the point XyzOut_V
  ! iFound is set to 0 if point is not found (outside domain), 1 otherwise.
  ! Set CellSizeOut_D to size of the cell containing the point.

  use ModReadAmr, ONLY: nVar, readamr_get
  use BATL_lib, ONLY: MaxDim
  use ModKind,    ONLY: Real8_
  
  implicit none
  real(Real8_), intent(in) :: XyzIn_D(MaxDim)
  real(Real8_), intent(out):: StateOut_V(0:nVar)
  real(Real8_), intent(out):: CellSizeOut_D(3)
  integer,intent(out):: iFound

  real   :: State_V(0:nVar)
  real   :: Xyz_D(MaxDim)
  real   :: CellSize_D(3)
  logical:: IsFound
  !----------------------------------------------------------------------------

  ! This copy converts real precision if needed
  Xyz_D = XyzIn_D   
  call readamr_get(Xyz_D, State_V, IsFound, CellSize_D)
 
  ! These copies convert real precision if needed
  StateOut_V = State_V
  CellSizeOut_D = CellSize_D

  ! Set integer found flag
  iFound = 0
  if(IsFound) iFound = 1

end subroutine wrapreadamr_get_data_cell

!=============================================================================
subroutine CON_stop(StringError)
  ! This subroutine has to be provided to stop cleanly

  use ModMPI
  use BATL_lib, ONLY: iProc, iComm

  implicit none
  integer:: nError, iError
  character (len=*), intent(in) :: StringError
  !--------------------------------------------------------------------------
  write(*,*)'ERROR in READAMR on processor ', iProc
  write(*,*) StringError
  call MPI_abort(iComm, nError, iError)
  stop

end subroutine CON_stop
