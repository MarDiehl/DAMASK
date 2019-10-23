!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief holds some global variables and gets extra information for commercial FEM
!--------------------------------------------------------------------------------------------------
module FEsolving
  use prec
  use IO
  use DAMASK_interface
   
  implicit none
  private
 
  logical, public :: &
    restartWrite      = .false., &                                                                  !< write current state to enable restart
    terminallyIll     = .false.                                                                     !< at least one material point is terminally ill

  integer, dimension(:,:), allocatable, public :: &
    FEsolving_execIP                                                                                !< for ping-pong scheme always range to max IP, otherwise one specific IP
  integer, dimension(2),                public :: &
    FEsolving_execElem                                                                              !< for ping-pong scheme always whole range, otherwise one specific element
    
#if defined(Marc4DAMASK) || defined(Abaqus)
  logical, public, protected :: & 
    symmetricSolver   = .false.                                                                     !< use a symmetric FEM solver (only Abaqus)
  logical, dimension(:,:), allocatable, public :: &
    calcMode                                                                                        !< do calculation or simply collect when using ping pong scheme

  public :: FE_init
#endif

contains


#if defined(Marc4DAMASK) || defined(Abaqus)
!--------------------------------------------------------------------------------------------------
!> @brief determine whether a symmetric solver is used
!--------------------------------------------------------------------------------------------------
subroutine FE_init
 
  integer, parameter :: &
    FILEUNIT = 222
  character(len=pStringLen) :: tag, line
  integer, allocatable, dimension(:) :: chunkPos

  write(6,'(/,a)')   ' <<<+-  FEsolving init  -+>>>'

  call IO_open_inputFile(FILEUNIT)
  rewind(FILEUNIT)
  do
    read (FILEUNIT,'(a256)',END=100) line
    chunkPos = IO_stringPos(line)
    tag = IO_lc(IO_stringValue(line,chunkPos,1))
    select case(tag)
      case ('solver')
        read (FILEUNIT,'(a256)',END=100) line                                                       ! next line
        chunkPos = IO_stringPos(line)
        symmetricSolver = (IO_intValue(line,chunkPos,2) /= 1)
    end select
  enddo
  100 close(FILEUNIT)

end subroutine FE_init
#endif

end module FEsolving
