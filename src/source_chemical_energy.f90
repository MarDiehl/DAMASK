!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for driving force due to chemical energy
!> @details to be done
!--------------------------------------------------------------------------------------------------
module source_chemical_energy
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   source_chemical_energy_sizePostResults, &                                                        !< cumulative size of post results
   source_chemical_energy_offset, &                                                                 !< which source is my current thermal dissipation mechanism?
   source_chemical_energy_instance                                                                  !< instance of thermal dissipation source mechanism

 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   source_chemical_energy_sizePostResult                                                            !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   source_chemical_energy_output                                                                    !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public :: &
   source_chemical_energy_Noutput                                                                   !< number of outputs per instance of this source 

 public :: &
   source_chemical_energy_init, &
   source_chemical_energy_getRateAndItsTangent

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine source_chemical_energy_init
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
   IO_timeStamp
 use material, only: &
   phase_source, &
   phase_Nsources, &
   phase_Noutput, &
   SOURCE_chemical_energy_label, &
   SOURCE_chemical_energy_ID, &
   material_Nphase, &
   material_phase, &  
   sourceState
 use numerics,only: &
   numerics_integrator

 implicit none
 integer(pInt) :: maxNinstance,phase,instance,source,sourceOffset
 integer(pInt) :: sizeState, sizeDotState, sizeDeltaState
 integer(pInt) :: NofMyPhase   

 write(6,'(/,a)')   ' <<<+-  source_'//SOURCE_chemical_energy_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
 maxNinstance = int(count(phase_source == SOURCE_chemical_energy_ID),pInt)
 if (maxNinstance == 0_pInt) return
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(source_chemical_energy_offset(material_Nphase), source=0_pInt)
 allocate(source_chemical_energy_instance(material_Nphase), source=0_pInt)
 do phase = 1, material_Nphase
   source_chemical_energy_instance(phase) = count(phase_source(:,1:phase) == SOURCE_chemical_energy_ID)
   do source = 1, phase_Nsources(phase)
     if (phase_source(source,phase) == SOURCE_chemical_energy_ID) &
       source_chemical_energy_offset(phase) = source
   enddo    
 enddo
   
 allocate(source_chemical_energy_sizePostResults(maxNinstance),                     source=0_pInt)
 allocate(source_chemical_energy_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(source_chemical_energy_output  (maxval(phase_Noutput),maxNinstance))
          source_chemical_energy_output = ''
 allocate(source_chemical_energy_Noutput(maxNinstance),                             source=0_pInt) 

 initializeInstances: do phase = 1_pInt, material_Nphase
   if (any(phase_source(:,phase) == SOURCE_chemical_energy_ID)) then
     NofMyPhase=count(material_phase==phase)
     instance = source_chemical_energy_instance(phase)
     sourceOffset = source_chemical_energy_offset(phase)

     sizeDotState              =   0_pInt
     sizeDeltaState            =   0_pInt
     sizeState                 =   0_pInt
     sourceState(phase)%p(sourceOffset)%sizeState =       sizeState
     sourceState(phase)%p(sourceOffset)%sizeDotState =    sizeDotState
     sourceState(phase)%p(sourceOffset)%sizeDeltaState =  sizeDeltaState
     sourceState(phase)%p(sourceOffset)%sizePostResults = source_chemical_energy_sizePostResults(instance)
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
end subroutine source_chemical_energy_init

!--------------------------------------------------------------------------------------------------
!> @brief returns chemical driving force 
!--------------------------------------------------------------------------------------------------
subroutine source_chemical_energy_getRateAndItsTangent(TDot, dTDOT_dT, ipc, ip, el)
 use material, only: &
   phase_chemicalFE, &
   material_phase, &
   CHEMICALFE_none_ID, &
   CHEMICALFE_quadenergy_ID, &
   CHEMICALFE_thermodynamic_ID
 use chemicalFE_quadenergy, only: &
   chemicalFE_quadenergy_getEnergy
 use chemicalFE_thermodynamic, only: &
   chemicalFE_thermodynamic_getEnergy

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),  intent(out) :: &
   TDot, &
   dTDOT_dT

 chemicalFEType: select case(phase_chemicalFE(material_phase(ipc, ip, el)))
   case (CHEMICALFE_none_ID) chemicalFEType
     TDot = 0.0_pReal
     dTDOT_dT = 0.0_pReal

   case (CHEMICALFE_quadenergy_ID) chemicalFEType
     TDot = chemicalFE_quadenergy_getEnergy(ipc,ip,el)
     dTDOT_dT = 0.0_pReal

   case (CHEMICALFE_thermodynamic_ID) chemicalFEType
     TDot = chemicalFE_thermodynamic_getEnergy(ipc,ip,el)
     dTDOT_dT = 0.0_pReal

 end select chemicalFEType
 
end subroutine source_chemical_energy_getRateAndItsTangent
 
end module source_chemical_energy
