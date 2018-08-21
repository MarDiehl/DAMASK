!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for purely mechanical material
!--------------------------------------------------------------------------------------------------
module currentDensity_none
 use prec, only: &
   pReal, &
   pInt, &
   p_vec

 implicit none
 private
 integer(pInt),                       dimension(:),     allocatable,          public, protected :: &
   currrentDensity_none_sizePostResults

 integer(pInt),                       dimension(:,:),   allocatable, target,  public :: &
   currentDensity_none_sizePostResult                                                                 !< size of each post result output

 public :: &
   currentDensity_none_init

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine currentDensity_none_init
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use debug, only: &
   debug_level, &
   debug_constitutive, &
   debug_levelBasic
 use IO, only: &
   IO_timeStamp
 use numerics, only: &
   numerics_integrator
 use material, only: &
   currentDensityMapping, &
   currentDensity, &
   phaseconstmemberAt, &
   phase_currentDensity, &
   phase_currentDensityInstance, &
   CURRENTDENSITY_NONE_label, &
   material_phase, &
   CURRENTDENSITY_none_ID, &
   currentDensityState
 implicit none

 integer(pInt) :: &
   instance, &
   maxNinstance, &
   phase, &
   NipcMyPhase, &
   sizeState, &
   sizeDotState, &
   sizeDeltaState
 
 write(6,'(/,a)')   ' <<<+-  currentDensity_'//CURRENTDENSITY_NONE_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
 maxNinstance = int(count(phase_currentDensity == CURRENTDENSITY_none_ID),pInt)
 if (maxNinstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance

 initializeInstances: do phase = 1_pInt, size(phase_currentDensity)
   if (phase_currentDensity(phase) == CURRENTDENSITY_none_ID) then
     NipcMyPhase = count(material_phase == phase)                                                    ! number of IPCs containing my phase  
     instance =  phase_currentDensityInstance(phase)                                                 ! which instance of my phase
!--------------------------------------------------------------------------------------------------
! allocate state arrays
     sizeState = 0_pInt
     sizeDotState = 0_pInt
     sizeDeltaState = 0_pInt
     currentDensityState(phase)%sizeState = sizeState
     currentDensityState(phase)%sizeDotState = sizeDotState
     currentDensityState(phase)%sizeDeltaState = sizeDeltaState
     currentDensityState(phase)%sizePostResults = 0.0_pReal
     allocate(currentDensityState(phase)%aTolState          (   sizeState),             source=0.0_pReal)
     allocate(currentDensityState(phase)%state0             (   sizeState,NipcMyPhase), source=0.0_pReal)
     allocate(currentDensityState(phase)%partionedState0    (   sizeState,NipcMyPhase), source=0.0_pReal)
     allocate(currentDensityState(phase)%subState0          (   sizeState,NipcMyPhase), source=0.0_pReal)
     allocate(currentDensityState(phase)%state              (   sizeState,NipcMyPhase), source=0.0_pReal)
     allocate(currentDensityState(phase)%dotState           (sizeDotState,NipcMyPhase), source=0.0_pReal)
     allocate(currentDensityState(phase)%deltaState       (sizeDeltaState,NipcMyPhase), source=0.0_pReal)
     if (any(numerics_integrator == 1_pInt)) then
       allocate(currentDensityState(phase)%previousDotState (sizeDotState,NipcMyPhase),source=0.0_pReal)
       allocate(currentDensityState(phase)%previousDotState2(sizeDotState,NipcMyPhase),source=0.0_pReal)
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(currentDensityState(phase)%RK4dotState      (sizeDotState,NipcMyPhase), source=0.0_pReal)
     if (any(numerics_integrator == 5_pInt)) &
       allocate(currentDensityState(phase)%RKCK45dotState (6,sizeDotState,NipcMyPhase), source=0.0_pReal)
       
     currentDensityMapping  (phase)%p => phaseconstmemberAt
     allocate(currentDensity(phase)%p(3,1), source=0.0_pReal)
   endif
 enddo initializeInstances
 
allocate(currrentDensity_none_sizePostResults(maxNinstance), source=0_pInt)
 
end subroutine currentDensity_none_init

end module currentDensity_none
