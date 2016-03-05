!--------------------------------------------------------------------------------------------------
!> @author   Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @brief    setting precision for real and int type depending on makros "FLOAT" and "INT"
!> @details  setting precision for real and int type and for DAMASK_NaN. Definition is made 
!!           depending on makros "FLOAT" and "INT" defined during compilation
!!           for details on NaN see https://software.intel.com/en-us/forums/topic/294680
!--------------------------------------------------------------------------------------------------
module prec

#if !(defined(__GFORTRAN__) && __GNUC__ < 5)
 use, intrinsic :: &                                                                                ! unfortunately not avialable in gfortran <= 5
  IEEE_arithmetic
#endif

 implicit none
 private 
#if (FLOAT==4)
#if defined(Spectral) || defined(FEM)
 SPECTRAL SOLVER AND OWN FEM DO NOT SUPPORT SINGLE PRECISION, STOPPING COMPILATION
#endif
 integer,     parameter, public :: pReal = 4                                                        !< floating point single precition (was selected_real_kind(6,37), number with 6 significant digits, up to 1e+-37)
#ifdef __INTEL_COMPILER
 real(pReal), parameter, public :: DAMASK_NaN = Z'7F800001'                                         !< quiet NaN for single precision (from http://www.hpc.unimelb.edu.au/doc/f90lrm/dfum_035.html, copy can be found in documentation/Code/Fortran)
#endif
#ifdef __GFORTRAN__
 real(pReal), parameter, public :: DAMASK_NaN = real(Z'7F800001', pReal)                            !< quiet NaN for single precision (from http://www.hpc.unimelb.edu.au/doc/f90lrm/dfum_035.html, copy can be found in documentation/Code/Fortran)
#endif
#elif (FLOAT==8)
 integer,     parameter, public :: pReal = 8                                                        !< floating point double precision (was selected_real_kind(15,300), number with 15 significant digits, up to 1e+-300)
#ifdef __INTEL_COMPILER
 real(pReal), parameter, public :: DAMASK_NaN = Z'7FF8000000000000'                                 !< quiet NaN for double precision (from http://www.hpc.unimelb.edu.au/doc/f90lrm/dfum_035.html, copy can be found in documentation/Code/Fortran)
#endif
#ifdef __GFORTRAN__
 real(pReal), parameter, public :: DAMASK_NaN = real(Z'7FF8000000000000',pReal)                     !< quiet NaN for double precision (from http://www.hpc.unimelb.edu.au/doc/f90lrm/dfum_035.html, copy can be found in documentation/Code/Fortran)
#endif
#else
 NO SUITABLE PRECISION FOR REAL SELECTED, STOPPING COMPILATION
#endif

#if (INT==4)
 integer,     parameter, public :: pInt  = 4                                                        !< integer representation 32 bit (was selected_int_kind(9), number with at least up to +- 1e9)
#elif (INT==8)
 integer,     parameter, public :: pInt  = 8                                                        !< integer representation 64 bit (was selected_int_kind(12), number with at least up to +- 1e12)
#else
 NO SUITABLE PRECISION FOR INTEGER SELECTED, STOPPING COMPILATION
#endif

 integer,     parameter, public :: pLongInt  = 8                                                    !< integer representation 64 bit (was selected_int_kind(12), number with at least up to +- 1e12)
 real(pReal), parameter, public :: tol_math_check = 1.0e-8_pReal                                    !< tolerance for internal math self-checks (rotation)

 integer(pInt), allocatable, dimension(:) :: realloc_lhs_test
 
 type, public :: p_vec                                                                              !< variable length datatype used for storage of state
   real(pReal), dimension(:), pointer :: p
 end type p_vec

 type, public :: p_intvec
   integer(pInt), dimension(:), pointer :: p
 end type p_intvec

!http://stackoverflow.com/questions/3948210/can-i-have-a-pointer-to-an-item-in-an-allocatable-array
 type, public :: tState
   integer(pInt) :: &
     sizeState = 0_pInt , &                                                                         !< size of state
     sizeDotState = 0_pInt, &                                                                       !< size of dot state, i.e. parts of the state that are integrated
     sizeDeltaState = 0_pInt, &                                                                     !< size of delta state, i.e. parts of the state that have discontinuous rates
     sizePostResults = 0_pInt                                                                       !< size of output data
   real(pReal), pointer,     dimension(:), contiguous :: &
     atolState
   real(pReal), pointer,     dimension(:,:), contiguous :: &                                        ! a pointer is needed here because we might point to state/doState. However, they will never point to something, but are rather allocated and, hence, contiguous 
     state, &                                                                                       !< state
     dotState, &                                                                                    !< state rate
     state0
   real(pReal), allocatable, dimension(:,:) :: &
     partionedState0, &
     subState0, &
     state_backup, &
     deltaState, &
     previousDotState, &                                                                            !< state rate of previous xxxx
     previousDotState2, &                                                                           !< state rate two xxxx ago
     dotState_backup, &                                                                             !< backup of state rate
     RK4dotState
   real(pReal), allocatable, dimension(:,:,:) :: &
     RKCK45dotState
 end type

 type, extends(tState), public :: tPlasticState
   integer(pInt) :: &
     nSlip = 0_pInt , &
     nTwin = 0_pInt, &
     nTrans = 0_pInt
   logical :: & 
     nonlocal = .false.                                                                             !< absolute tolerance for state integration
   real(pReal), pointer,     dimension(:,:), contiguous :: &
     slipRate, &                                                                                    !< slip rate
     accumulatedSlip                                                                                !< accumulated plastic slip
 end type

 type, public :: tSourceState
   type(tState), dimension(:), allocatable :: p                                                     !< tState for each active source mechanism in a phase
 end type
 
 type, public :: tHomogMapping
   integer(pInt), pointer, dimension(:,:) :: p                                  
 end type 

 type, public :: tPhaseMapping
   integer(pInt), pointer, dimension(:,:,:) :: p
 end type 

#ifdef FEM
 type, public :: tOutputData
   integer(pInt) :: &
     sizeIpCells = 0_pInt , &
     sizeResults = 0_pInt
   real(pReal), allocatable, dimension(:,:) :: &
     output                                                                                         !< output data
 end type 
#endif

 public :: &
   prec_init, &
   prec_isNaN
 
contains


!--------------------------------------------------------------------------------------------------
!> @brief reporting precision and checking if DAMASK_NaN is set correctly
!--------------------------------------------------------------------------------------------------
subroutine prec_init
 use, intrinsic :: &
   iso_fortran_env                                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)

 implicit none
 integer(pInt) :: worldrank = 0_pInt
#ifdef PETSc
#include <petsc/finclude/petscsys.h>
 PetscErrorCode :: ierr
#endif
 external :: &
   quit, &
   MPI_Comm_rank, &
   MPI_Abort
   
#ifdef PETSc
 call MPI_Comm_rank(PETSC_COMM_WORLD,worldrank,ierr);CHKERRQ(ierr)
#endif

 mainProcess: if (worldrank == 0) then
   write(6,'(/,a)') ' <<<+-  prec init  -+>>>'
#include "compilation_info.f90"
   write(6,'(a,i3)')    ' Bytes for pReal:    ',pReal
   write(6,'(a,i3)')    ' Bytes for pInt:     ',pInt
   write(6,'(a,i3)')    ' Bytes for pLongInt: ',pLongInt
   write(6,'(a,e10.3)') ' NaN:           ',     DAMASK_NaN
   write(6,'(a,l3)')    ' NaN != NaN:         ',DAMASK_NaN /= DAMASK_NaN
   write(6,'(a,l3,/)')  ' NaN check passed    ',prec_isNAN(DAMASK_NaN)
 endif mainProcess

 if ((.not. prec_isNaN(DAMASK_NaN)) .or. (DAMASK_NaN == DAMASK_NaN)) call quit(9000)
 realloc_lhs_test = [1_pInt,2_pInt]
 if (realloc_lhs_test(2)/=2_pInt) call quit(9000)
 

end subroutine prec_init


!--------------------------------------------------------------------------------------------------
!> @brief figures out if a floating point number is NaN
! basically just a small wrapper, because gfortran < 5.0 does not have the IEEE module
!--------------------------------------------------------------------------------------------------
logical elemental pure function prec_isNaN(a)

 implicit none
 real(pReal), intent(in) :: a

#if (defined(__GFORTRAN__) && __GNUC__ < 5)
 intrinsic :: isNaN
 prec_isNaN = isNaN(a)
#else
 prec_isNaN = IEEE_is_NaN(a)
#endif
end function prec_isNaN



!--------------------------------------------------------------------------------------------------
!> @brief equality comparison for double precision
! replaces "==" but for certain (absolute) tolerance. Counterpart to dNeq
!--------------------------------------------------------------------------------------------------
logical elemental pure function dEq(a,b,tol)
  real(pReal), intent(in) :: a,b
  real(pReal), intent(in), optional :: tol
  real(pReal), parameter  :: eps = 1.0e-15_pReal 
  dEq = merge(.True., .False.,abs(a-b) < merge(tol,eps,present(tol)))
end function dEq



!--------------------------------------------------------------------------------------------------
!> @brief inequality comparison for double precision
! replaces "!=" but for certain (absolute) tolerance. Counterpart to dEq
!--------------------------------------------------------------------------------------------------
logical elemental pure function dNeq(a,b,tol)
  real(pReal), intent(in)           :: a,b
  real(pReal), intent(in), optional :: tol
  real(pReal), parameter  :: eps = 1.0e-15_pReal 
  dNeq = merge(.True., .False.,abs(a-b) < merge(tol,eps,present(tol)))
end function dNeq

end module prec
