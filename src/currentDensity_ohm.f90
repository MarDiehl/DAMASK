!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for the current density using ohm's law
!--------------------------------------------------------------------------------------------------
module currentDensity_ohm
 use prec, only: &
   pReal,&
   pInt, &
   p_vec

 implicit none
 private
 real(pReal),                                          parameter,           private :: &
   R = 8.314459848_pReal                                                                          !< Universal gas constant in J/(mol Kelvin)

 integer(pInt),                       dimension(:,:),   allocatable, target, public :: &
   currentDensity_ohm_sizePostResult                                                        !< size of each post result output

 character(len=64),                   dimension(:,:),   allocatable, target, public :: &
   currentDensity_ohm_output                                                                !< name of each post result output

 enum, bind(c)
   enumerator :: undefined_ID, &
                 currentdensity_ID
 end enum

 type, private :: tParameters                                                                         !< container type for internal constitutive parameters
   integer(kind(undefined_ID)), dimension(:), allocatable :: & 
     outputID
   real(pReal) :: &
     conductivity = 0.0_pReal
 end type

 type(tParameters), dimension(:), allocatable, private :: param                                       !< containers of constitutive parameters (len Ninstance)

 public :: &
   currentDensity_ohm_init, &
   currentDensity_ohm_getFlux, &
   currentDensity_ohm_getFluxTangent, &
   currentDensity_ohm_getElectricField, &
   currentDensity_ohm_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine currentDensity_ohm_init
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use prec, only: &
   dEq0
 use debug, only: &
   debug_level, &
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
   currentDensityMapping, &
   currentDensity, &
   phasememberAt, &
   phase_currentDensity, &
   phase_currentDensityInstance, &
   phase_Noutput, &
   CURRENTDENSITY_OHM_label, &
   CURRENTDENSITY_OHM_ID, &
   material_phase, &
   currentDensityState
 use config, only: &
   MATERIAL_partPhase, &
   config_phase
 use numerics,only: &
   numerics_integrator

 implicit none

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: &
   maxNinstance, &
   instance,phase,o,p, &
   NipcMyPhase, &
   outputSize, &
   sizeState,sizeDotState, sizeDeltaState
 character(len=65536), dimension(:), allocatable :: outputs
 character(len=65536), dimension(0), parameter :: emptyStringArray = [character(len=65536)::]
 integer(kind(undefined_ID)) :: outputID  

 write(6,'(/,a)')   ' <<<+-  currentDensity_'//CURRENTDENSITY_OHM_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = int(count(phase_currentDensity == CURRENTDENSITY_OHM_ID),pInt)
 if (maxNinstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance


 allocate(currentDensity_ohm_sizePostResult(maxval(phase_Noutput),maxNinstance), &
                                                                                 source=0_pInt)
 allocate(currentDensity_ohm_output(maxval(phase_Noutput),maxNinstance))
          currentDensity_ohm_output               = ''
  
 allocate(param(maxNinstance))
 
 parsingFile:do p = 1_pInt, size(phase_currentDensityInstance)
   if (phase_currentDensity(p) /= CURRENTDENSITY_ohm_ID) cycle
   instance = phase_currentDensityInstance(p)
      
   param(instance)%conductivity & 
              = config_phase(p)%getFloat('electrical_conductivity',defaultVal=0.0_pReal)
    
   outputs = config_phase(p)%getStrings('(output)',defaultVal=emptyStringArray)
   allocate(param(instance)%outputID(0))   
   do o=1_pInt, size(outputs)
     outputID = undefined_ID
     select case(outputs(o))
       case ('currentdensity')
        outputID = currentdensity_ID
        outputSize = 3_pInt
       end select 
       
     if (outputID /= undefined_ID) then
      currentDensity_ohm_output(o,instance) = outputs(o)
      currentDensity_ohm_sizePostResult(o,instance) = outputSize
      param(instance)%outputID = [param(instance)%outputID , outputID]
     endif
   enddo
  enddo parsingFile

 initializeInstances: do phase = 1_pInt, size(phase_currentDensity)
   myPhase2: if (phase_currentDensity(phase) == CURRENTDENSITY_ohm_ID) then                                    ! only consider my phase
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
     currentDensityState(phase)%sizePostResults = sum(currentDensity_ohm_sizePostResult(:,instance))
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
     
     currentDensityMapping(phase)%p => phasememberAt
     allocate(currentDensity(phase)%p(3,NipcMyPhase), source = 0.0_pReal)
   endif myPhase2
 enddo initializeInstances

end subroutine currentDensity_ohm_init


!--------------------------------------------------------------------------------------------------
!> @brief returns the current density per (ipc,ip,el)
!--------------------------------------------------------------------------------------------------
subroutine currentDensity_ohm_getFlux(CurrentDensity,ElectricalField,ipc,ip,el)
 use material, only: &
   phase_currentDensityInstance, &
   material_phase
   
 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
real(pReal), dimension(3),   intent(in)  :: &
   ElectricalField
real(pReal), dimension(3),   intent(out)  :: &
   CurrentDensity   
 integer(pInt) :: &
   instance
   
 instance = phase_currentDensityInstance(material_phase(ipc,ip,el)) 
 CurrentDensity = param(instance)%conductivity*ElectricalField
 
end subroutine currentDensity_ohm_getFlux



!--------------------------------------------------------------------------------------------------
!> @brief returns the current density tangent (conductivity) per (ipc,ip,el)
!--------------------------------------------------------------------------------------------------
subroutine currentDensity_ohm_getFluxTangent(CurrentDensityTangent,ipc,ip,el)
 use material, only: &
   phase_currentDensityInstance, &
   material_phase
   
 use math, only: &
   MATH_I3
   
 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(3,3),   intent(out)  :: &
   CurrentDensityTangent  
 integer(pInt) :: &
   instance
   
 instance = phase_currentDensityInstance(material_phase(ipc,ip,el)) 
 CurrentDensityTangent = param(instance)%conductivity*MATH_I3

end subroutine currentDensity_ohm_getFluxTangent



!--------------------------------------------------------------------------------------------------
!> @brief returns the local Electric field per (ipc,ip,el)
!--------------------------------------------------------------------------------------------------   
subroutine currentDensity_ohm_getElectricField(ElectricalField,CurrentDensity,ipc,ip,el)
 use material, only: &
   material_phase
 use math, only: &  
   math_inv33, &
   math_mul33X3
    
 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
real(pReal), dimension(3),   intent(out)  :: &
   ElectricalField
real(pReal), dimension(3),   intent(in)  :: &
   CurrentDensity
real(pReal), dimension(3,3) :: &
   electrical_conductivity   
    
 call currentDensity_ohm_getFluxTangent(electrical_conductivity,ipc,ip,el)
  
 ElectricalField = math_mul33X3(math_inv33(electrical_conductivity),CurrentDensity)
  
end subroutine currentDensity_ohm_getElectricField



!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results 
!--------------------------------------------------------------------------------------------------
function currentDensity_ohm_postResults(ipc,ip,el)
 use material, only: &
   currentDensity, &
   currentDensityMapping, &
   currentDensityState, &
   phase_currentDensityInstance, &
   material_phase
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),  dimension(currentDensityState(material_phase(ipc,ip,el))%sizePostResults) :: &
   currentDensity_ohm_postResults
 
 integer(pInt) :: &
   o, c, phase, instance
   
 phase = material_phase(ipc,ip,el)
 instance = phase_currentDensityInstance(phase)

 currentDensity_ohm_postResults = 0.0_pReal  
    
 c = 0_pInt
 
 outputsLoop:do o = 1_pInt, size(param(instance)%outputID)
   select case(param(instance)%outputID(o))
     case (currentdensity_ID)
    
       currentDensity_ohm_postResults(c+1:c+3) = &
           currentDensity(phase)%p(1:3,currentDensityMapping(phase)%p(ipc,ip,el))
      
       c = c + 3

   end select
 enddo outputsLoop 

end function currentDensity_ohm_postResults


end module currentDensity_ohm
