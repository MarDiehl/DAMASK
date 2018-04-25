!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for purely mechanical material
!--------------------------------------------------------------------------------------------------
module chemicalFE_none
 use prec, only: &
   pReal, &
   pInt, &
   p_vec

 implicit none
 private
 integer(pInt),                       dimension(:),     allocatable,          public, protected :: &
   chemicalFE_none_sizePostResults

 integer(pInt),                       dimension(:,:),   allocatable, target,  public :: &
   chemicalFE_none_sizePostResult                                                                 !< size of each post result output

 public :: &
   chemicalFE_none_init

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine chemicalFE_none_init
#ifdef __GFORTRAN__
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
   chemConcMapping, &
   chemicalConc, &
   chemicalConcRate, &
   phase_Ncomponents, &
   phaseconstmemberAt, &
   phase_chemicalFE, &
   CHEMICALFE_NONE_label, &
   material_phase, &
   chemicalState, &
   CHEMICALFE_none_ID

 implicit none

 integer(pInt) :: &
   maxNinstance, &
   phase, &
   comp, &
   NofMyPhase, &
   sizeState, &
   sizeDotState, &
   sizeDeltaState
 
 write(6,'(/,a)')   ' <<<+-  chemicalFE_'//CHEMICALFE_NONE_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
 maxNinstance = int(count(phase_chemicalFE == CHEMICALFE_none_ID),pInt)
 if (maxNinstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance

 initializeInstances: do phase = 1_pInt, size(phase_chemicalFE)
   if (phase_chemicalFE(phase) == CHEMICALFE_none_ID) then
     NofMyPhase=count(material_phase==phase)

     sizeState    = 0_pInt
     chemicalState(phase)%sizeState = sizeState
     sizeDotState = sizeState
     chemicalState(phase)%sizeDotState = sizeDotState
     sizeDeltaState = 0_pInt
     chemicalState(phase)%sizeDeltaState = sizeDeltaState
     chemicalState(phase)%sizePostResults = 0_pInt
     allocate(chemicalState(phase)%aTolState          (sizeState))
     allocate(chemicalState(phase)%state0             (sizeState,NofMyPhase))
     allocate(chemicalState(phase)%partionedState0    (sizeState,NofMyPhase))
     allocate(chemicalState(phase)%subState0          (sizeState,NofMyPhase))
     allocate(chemicalState(phase)%state              (sizeState,NofMyPhase))

     allocate(chemicalState(phase)%dotState           (sizeDotState,NofMyPhase))
     allocate(chemicalState(phase)%deltaState        (sizeDeltaState,NofMyPhase))
     if (any(numerics_integrator == 1_pInt)) then
       allocate(chemicalState(phase)%previousDotState (sizeDotState,NofMyPhase))
       allocate(chemicalState(phase)%previousDotState2(sizeDotState,NofMyPhase))
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(chemicalState(phase)%RK4dotState      (sizeDotState,NofMyPhase))
     if (any(numerics_integrator == 5_pInt)) &
       allocate(chemicalState(phase)%RKCK45dotState (6,sizeDotState,NofMyPhase))
     chemConcMapping(phase)%p => phaseconstmemberAt
     do comp = 1_pInt, phase_Ncomponents(phase)
       allocate(chemicalConc    (comp,phase)%p(1), source=0.0_pReal)
       allocate(chemicalConcRate(comp,phase)%p(1), source=0.0_pReal)
     enddo
   endif
 enddo initializeInstances

 allocate(chemicalFE_none_sizePostResults(maxNinstance), source=0_pInt)

end subroutine chemicalFE_none_init

end module chemicalFE_none
