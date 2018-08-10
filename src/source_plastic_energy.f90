!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for driving force due to plastic energy
!> @details to be done
!--------------------------------------------------------------------------------------------------
module source_plastic_energy
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   source_plastic_energy_sizePostResults, &                                                        !< cumulative size of post results
   source_plastic_energy_offset, &                                                                 !< which source is my current thermal dissipation mechanism?
   source_plastic_energy_instance                                                                  !< instance of thermal dissipation source mechanism

 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   source_plastic_energy_sizePostResult                                                            !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   source_plastic_energy_output                                                                    !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public :: &
   source_plastic_energy_Noutput                                                                   !< number of outputs per instance of this source 
 
 type, private :: tParameters                                                                      !< container type for internal constitutive parameters
   real(pReal) :: &
     source_plastic_work_storeCoeff
 end type  
 
 type(tParameters), dimension(:), allocatable, private :: param                                    !< containers of constitutive parameters (len Ninstance)
   
 public :: &
   source_plastic_energy_init, &
   source_plastic_energy_dotState, &
   source_plastic_energy_getRateAndItsTangent

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine source_plastic_energy_init(fileUnit)
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
   SOURCE_plastic_energy_label, &
   SOURCE_plastic_energy_ID, &
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

 write(6,'(/,a)')   ' <<<+-  source_'//SOURCE_plastic_energy_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
 maxNinstance = int(count(phase_source == SOURCE_plastic_energy_ID),pInt)
 if (maxNinstance == 0_pInt) return
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(source_plastic_energy_offset(material_Nphase), source=0_pInt)
 allocate(source_plastic_energy_instance(material_Nphase), source=0_pInt)
 do phase = 1, material_Nphase
   source_plastic_energy_instance(phase) = count(phase_source(:,1:phase) == SOURCE_plastic_energy_ID)
   do source = 1, phase_Nsources(phase)
     if (phase_source(source,phase) == SOURCE_plastic_energy_ID) &
       source_plastic_energy_offset(phase) = source
   enddo    
 enddo
   
 allocate(source_plastic_energy_sizePostResults(maxNinstance),                     source=0_pInt)
 allocate(source_plastic_energy_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(source_plastic_energy_output  (maxval(phase_Noutput),maxNinstance))
          source_plastic_energy_output = ''
 allocate(source_plastic_energy_Noutput(maxNinstance),                             source=0_pInt) 
  
 allocate(param(maxNinstance))
 
 rewind(fileUnit)
 phase = 0_pInt
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= MATERIAL_partPhase)         ! wind forward to <phase>
   line = IO_read(fileUnit)
 enddo
 
 parsingFile: do while (trim(line) /= IO_EOF)                                                       ! read through sections of phase part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit                                                                                           
   endif   
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next phase section
     phase = phase + 1_pInt                                                                         ! advance phase section counter
     cycle                                                                                          ! skip to next line
   endif

   if (phase > 0_pInt ) then; if (any(phase_source(:,phase) == SOURCE_plastic_energy_ID)) then      ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran

     instance = source_plastic_energy_instance(phase)                                               ! which instance of my source is present phase
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                              ! extract key
     select case(tag)
       case ('plastic_work_store_coeff')
         param(instance)%source_plastic_work_storeCoeff = IO_floatValue(line,chunkPos,2_pInt)

     end select
   endif; endif
 enddo parsingFile

 initializeInstances: do phase = 1_pInt, material_Nphase
   if (any(phase_source(:,phase) == SOURCE_plastic_energy_ID)) then
     NofMyPhase=count(material_phase==phase)
     instance = source_plastic_energy_instance(phase)
     sourceOffset = source_plastic_energy_offset(phase)

     sizeDotState              =   1_pInt
     sizeDeltaState            =   0_pInt
     sizeState                 =   1_pInt
     sourceState(phase)%p(sourceOffset)%sizeState =       sizeState
     sourceState(phase)%p(sourceOffset)%sizeDotState =    sizeDotState
     sourceState(phase)%p(sourceOffset)%sizeDeltaState =  sizeDeltaState
     sourceState(phase)%p(sourceOffset)%sizePostResults = source_plastic_energy_sizePostResults(instance)
     allocate(sourceState(phase)%p(sourceOffset)%aTolState           (sizeState),                source=0.1_pReal)
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
end subroutine source_plastic_energy_init

!--------------------------------------------------------------------------------------------------
!> @brief rate of change of state

!--------------------------------------------------------------------------------------------------
subroutine source_plastic_energy_dotState(Tstar_v, Lp, ipc, ip, el)
 use material, only: &
   phaseAt, phasememberAt, &
   sourceState
 use math, only: &  
   math_Mandel6to33, &
   math_mul33xx33
   
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)
 real(pReal),  intent(in), dimension(3,3) :: &                                                     
   Lp                                                                                               !velocity gradient 
 real(pReal),  dimension(3,3) :: &
   Tstar    
 integer(pInt) :: &
   phase, &
   constituent, &
   sourceOffset

 phase = phaseAt(ipc,ip,el)
 constituent = phasememberAt(ipc,ip,el)
 sourceOffset = source_plastic_energy_offset(phase)
 
 Tstar = math_Mandel6to33(Tstar_v)
 
 sourceState(phase)%p(sourceOffset)%dotState(1,constituent) = math_mul33xx33(Tstar,Lp)               ! state is current time
 
 end subroutine source_plastic_energy_dotState

!--------------------------------------------------------------------------------------------------
!> @brief returns plastic driving force 
!--------------------------------------------------------------------------------------------------
subroutine source_plastic_energy_getRateAndItsTangent(TDot, dTDOT_dT, ipc, ip, el)
 use material, only: &
   phaseAt, phasememberAt, &
   phase_plasticity, &
   sourceState, &
   material_phase, &
   PLASTICITY_NONE_ID, &
   PLASTICITY_ISOTROPIC_ID, &
   PLASTICITY_PHENOPOWERLAW_ID, &
   PLASTICITY_DISLOTWIN_ID, &
   PLASTICITY_DISLOUCLA_ID, &
   PLASTICITY_NONLOCAL_ID

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),  intent(out) :: &
   TDot, &
   dTDOT_dT
 integer(pInt) :: &
   phase, &
   constituent, &
   sourceOffset, &
   instance

 phase = phaseAt(ipc,ip,el)
 constituent = phasememberAt(ipc,ip,el)
 sourceOffset = source_plastic_energy_offset(phase) 
 instance = source_plastic_energy_instance(phase)

 plasticityType: select case(phase_plasticity(material_phase(ipc, ip, el)))
   case (PLASTICITY_NONE_ID) plasticityType
     TDot = 0.0_pReal
     dTDOT_dT = 0.0_pReal

   case (PLASTICITY_ISOTROPIC_ID) plasticityType
     TDot = param(instance)%source_plastic_work_storeCoeff* & 
            sourceState(phase)%p(sourceOffset)%State(1,constituent)
     dTDOT_dT = 0.0_pReal

   case (PLASTICITY_PHENOPOWERLAW_ID) plasticityType
     TDot = param(instance)%source_plastic_work_storeCoeff* & 
            sourceState(phase)%p(sourceOffset)%State(1,constituent)
     dTDOT_dT = 0.0_pReal

   case (PLASTICITY_DISLOTWIN_ID) plasticityType
     TDot = param(instance)%source_plastic_work_storeCoeff* & 
            sourceState(phase)%p(sourceOffset)%State(1,constituent)
     dTDOT_dT = 0.0_pReal

   case (PLASTICITY_DISLOUCLA_ID) plasticityType
     TDot = param(instance)%source_plastic_work_storeCoeff* & 
            sourceState(phase)%p(sourceOffset)%State(1,constituent)
     dTDOT_dT = 0.0_pReal

   case (PLASTICITY_NONLOCAL_ID) plasticityType
     TDot = param(instance)%source_plastic_work_storeCoeff* & 
            sourceState(phase)%p(sourceOffset)%State(1,constituent)
     dTDOT_dT = 0.0_pReal

 end select plasticityType
 
end subroutine source_plastic_energy_getRateAndItsTangent
 
end module source_plastic_energy
