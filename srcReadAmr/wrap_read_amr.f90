
! Extrernal subroutines to call from other languages
!==============================================================================
subroutine wrapreadamr_set_mpi_param(iProcIn, nProcIn, iCommIn)
  use BATL_lib,  ONLY: iProc, nProc, iComm

  ! Set the MPI parameters
  
  implicit none
  integer,intent(in):: iProcIn, nProcIn, iCommIn
  
  iProc = iProcIn
  nProc = nProcIn 
  iComm = iCommIn

end subroutine wrapreadamr_set_mpi_param

!==============================================================================
subroutine wrapreadamr_get_nvar(nVarOut)
  use ModReadAmr, ONLY: nVar

  implicit none
  integer,intent(out):: nVarOut
  
  nVarOut = nVar 
end subroutine wrapreadamr_get_nvar    
!==============================================================================
subroutine wrapreadamr_get_namevardata(l, NameVarOut)

  ! Names of variables that are returned by the interpolation routine 

  use ModReadAmr, ONLY:NameVarData
  
  implicit none
  integer,          intent(in) :: l
  character(len=l), intent(out):: NameVarOut
 
  NameVarOut = NameVarData(1:l)

end subroutine wrapreadamr_get_namevardata  

!==============================================================================
subroutine wrapreadamr_get_nameunitdata(l, NameUnitOut)

  ! units of the coordinates followed by the units of the interpolated variables
  use ModReadAmr, ONLY: NameUnitData
  
  implicit none
  integer, intent(in):: l
  character(len=l), intent(out):: NameUnitOut
 
  NameUnitOut = NameUnitData(1:l)

end subroutine wrapreadamr_get_nameunitdata  
  
!===============================================================================
subroutine wrapreadamr_domain_limits(CoordMinOut_D, CoordMaxOut_D) 

  ! Boundary of the computational domain in the normalized units

  use ModKind,    ONLY: Real8_
  use ModReadAmr, ONLY: CoordMin_D, CoordMax_D
  
  implicit none
  real(Real8_), intent(out):: CoordMinOut_D(3), CoordMaxOut_D(3)

  CoordMaxOut_D = CoordMin_D
  CoordMinOut_D = CoordMax_D

end subroutine wrapreadamr_domain_limits  

!==============================================================================
subroutine wrapreadamr_read_file(l, FileName, lNewGrid, lVerbose)

  ! Read data from file Filename.
  ! If lNewGrid == 1 then read all information, otherwise data only
  ! If lVerbose == 1 then write out verbose information

  use ModReadAmr, ONLY:readamr_read 

  implicit none

  integer,intent(in):: l
  character(len=l), intent(in):: FileName
  integer,          intent(in):: lNewGrid
  integer,          intent(in):: lVerbose

  call readamr_read(FileName, &
       IsNewGridIn = lNewGrid==1, IsVerboseIn=lVerbose==1)

end subroutine wrapreadamr_read_file

!=============================================================================
subroutine wrapreadamr_read_header(l, FileName)

  ! Read header information from .h or .info
  ! This sets number and names of variables, domain size, etc.
  use ModReadAmr, ONLY:readamr_init
  
  implicit none
  integer,intent(in):: l
  character(len=l), intent(in):: FileName
  integer:: i
  !---------------------------------------------------------------------------
  ! Cut off extension from the file name (if any)
  i = index(FileName,".",BACK=.true.)
  call readamr_init(FileName(1:i-1), .true.)

end subroutine wrapreadamr_read_header
  
!==============================================================================
subroutine wrapreadamr_deallocate

  ! Deallocate all memory used by READAMR

  use ModReadAmr, ONLY: readamr_clean
  
  call readamr_clean
  
end subroutine wrapreadamr_deallocate

!=============================================================================
subroutine wrapreadamr_get_data(XyzIn_D, StateOut_V, iFound) 

  ! Get the interpolated values StateOut_V at the point XyzOut_V
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
subroutine wrapreadamr_get_data_cell(XyzIn_D, StateOut_V, &
     CellSizeOut_D, iFound)


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

  
  
