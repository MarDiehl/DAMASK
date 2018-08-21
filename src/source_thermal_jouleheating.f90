!--------------------------------------------------------------------------------------------------
!> @author Arka Lahiri, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for thermal source due to joule heating
!> @details to be done
!--------------------------------------------------------------------------------------------------
module source_thermal_jouleheating
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   source_thermal_jouleheating_sizePostResults, &                                                        !< cumulative size of post results
   source_thermal_jouleheating_offset, &                                                                 !< which source is my current thermal joule heating mechanism?
   source_thermal_jouleheating_instance                                                                  !< instance of thermal joule heating source mechanism

 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   source_thermal_jouleheating_sizePostResult                                                            !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   source_thermal_jouleheating_output                                                                    !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public :: &
   source_thermal_jouleheating_Noutput                                                                   !< number of outputs per instance of this source 


 public :: &
   source_thermal_jouleheating_init, &
   source_thermal_jouleheating_getRateAndItsTangent

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine source_thermal_jouleheating_init(fileUnit)
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use debug, only: &
   debug_level,&
   debug_constitutive,&
   debug_levelBasic
 use IO, only: &
   IO_read, &
   IO_lc, &
   IO_getTag, &
   IO_isBlank, &
   IO_stringPos, &
   IO_stringValue, &
   IO_floatValue, &
   IO_intValue, &
   IO_warning, &
   IO_error, &
   IO_timeStamp, &
   IO_EOF
 use material, only: &
   phase_source, &
   phase_Nsources, &
   phase_Noutput, &
   SOURCE_thermal_jouleheating_label, &
   SOURCE_thermal_jouleheating_ID, &
   material_phase, &  
   sourceState
 use config, only: &
   material_Nphase, &
   MATERIAL_partPhase  
 use numerics,only: &
   numerics_integrator

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: maxNinstance,phase,instance,source,sourceOffset
 integer(pInt) :: sizeState, sizeDotState, sizeDeltaState
 integer(pInt) :: NofMyPhase   
 character(len=65536) :: &
   tag  = '', &
   line = ''

 write(6,'(/,a)')   ' <<<+-  source_'//SOURCE_thermal_jouleheating_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
 maxNinstance = int(count(phase_source == SOURCE_thermal_jouleheating_ID),pInt)
 if (maxNinstance == 0_pInt) return
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(source_thermal_jouleheating_offset(material_Nphase), source=0_pInt)
 allocate(source_thermal_jouleheating_instance(material_Nphase), source=0_pInt)
 do phase = 1, material_Nphase
   source_thermal_jouleheating_instance(phase) = count(phase_source(:,1:phase) == SOURCE_thermal_jouleheating_ID)
   do source = 1, phase_Nsources(phase)
     if (phase_source(source,phase) == SOURCE_thermal_jouleheating_ID) &
       source_thermal_jouleheating_offset(phase) = source
   enddo    
 enddo
   
 allocate(source_thermal_jouleheating_sizePostResults(maxNinstance),                     source=0_pInt)
 allocate(source_thermal_jouleheating_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(source_thermal_jouleheating_output  (maxval(phase_Noutput),maxNinstance))
          source_thermal_jouleheating_output = ''
 allocate(source_thermal_jouleheating_Noutput(maxNinstance),                             source=0_pInt) 

 
 initializeInstances: do phase = 1_pInt, material_Nphase
   if (any(phase_source(:,phase) == source_thermal_jouleheating_ID)) then
     NofMyPhase=count(material_phase==phase)
     instance = source_thermal_jouleheating_instance(phase)
     sourceOffset = source_thermal_jouleheating_offset(phase)

     sizeDotState              =   0_pInt
     sizeDeltaState            =   0_pInt
     sizeState                 =   0_pInt
     sourceState(phase)%p(sourceOffset)%sizeState =       sizeState
     sourceState(phase)%p(sourceOffset)%sizeDotState =    sizeDotState
     sourceState(phase)%p(sourceOffset)%sizeDeltaState =  sizeDeltaState
     sourceState(phase)%p(sourceOffset)%sizePostResults = source_thermal_jouleheating_sizePostResults(instance)
     allocate(sourceState(phase)%p(sourceOffset)%aTolState           (sizeState),                source=0.0_pReal)
     allocate(sourceState(phase)%p(sourceOffset)%state0              (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(sourceState(phase)%p(sourceOffset)%partionedState0     (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(sourceState(phase)%p(sourceOffset)%subState0           (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(sourceState(phase)%p(sourceOffset)%state               (sizeState,NofMyPhase),     source=0.0_pReal)

     allocate(sourceState(phase)%p(sourceOffset)%dotState            (sizeDotState,NofMyPhase),  source=0.0_pReal)
     allocate(sourceState(phase)%p(sourceOffset)%deltaState        (sizeDeltaState,NofMyPhase),  source=0.0_pReal)
     if (any(numerics_integrator == 1_pInt)) then
       allocate(sourceState(phase)%p(sourceOffset)%previousDotState  (sizeDotState,NofMyPhase),  source=0.0_pReal)
       allocate(sourceState(phase)%p(sourceOffset)%previousDotState2 (sizeDotState,NofMyPhase),  source=0.0_pReal)
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(sourceState(phase)%p(sourceOffset)%RK4dotState       (sizeDotState,NofMyPhase),  source=0.0_pReal)
     if (any(numerics_integrator == 5_pInt)) &
       allocate(sourceState(phase)%p(sourceOffset)%RKCK45dotState    (6,sizeDotState,NofMyPhase),source=0.0_pReal)      

   endif
 
 enddo initializeInstances
end subroutine source_thermal_jouleheating_init

!--------------------------------------------------------------------------------------------------
!> @brief returns local joule heating rate 
!--------------------------------------------------------------------------------------------------
subroutine source_thermal_jouleheating_getRateAndItsTangent(TDot, dTDOT_dT, ipc, ip, el)
 use material, only: &
   phaseAt, &
   currentDensityMapping, &
   currentDensity
 use math, only: &  
   math_inv33, &
   math_mul33X3, &
   math_mul3X3
 use currentDensity_ohm, only: &
   currentDensity_ohm_getFluxTangent
   
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),  intent(out) :: &
   TDot, &
   dTDOT_dT
 integer(pInt) :: &
   phase, offset
 real(pReal), dimension(3,3) :: &
   electrical_conductivity
 real(pReal), dimension(3) :: &
   electricalField, &
   currentDensityLocal

 phase                 = phaseAt(ipc,ip,el)
 offset                = currentDensityMapping(phase)%p(ipc,ip,el)
 currentDensityLocal   = currentDensity(phase)%p(1:3,offset)
 
 call currentDensity_ohm_getFluxTangent(electrical_conductivity,ipc,ip,el)
  
 electricalField = math_mul33X3(math_inv33(electrical_conductivity),currentDensityLocal)
  
 TDot = math_mul3x3(currentDensityLocal,electricalField)
 
 dTDOT_dT = 0.0_pReal       
 
end subroutine source_thermal_jouleheating_getRateAndItsTangent
 
end module source_thermal_jouleheating
