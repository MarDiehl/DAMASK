!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for vacancy generation due to irradiation
!> @details to be done
!--------------------------------------------------------------------------------------------------
module source_stochastic_phase_nucleation
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   source_stochastic_phase_nucleation_sizePostResults, &                                                        !< cumulative size of post results
   source_stochastic_phase_nucleation_offset, &                                                                 !< which source is my current damage mechanism?
   source_stochastic_phase_nucleation_instance                                                                  !< instance of damage source mechanism

 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   source_stochastic_phase_nucleation_sizePostResult                                                            !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   source_stochastic_phase_nucleation_output                                                                    !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public :: &
   source_stochastic_phase_nucleation_Noutput                                                                   !< number of outputs per instance of this damage 

 real(pReal),                         dimension(:),           allocatable,        private :: &
   source_stochastic_phase_nucleation_probability, &
   source_stochastic_phase_nucleation_sourcestrength

 public :: &
   source_stochastic_phase_nucleation_init, &
   source_stochastic_phase_nucleation_deltaState, &
   source_stochastic_phase_nucleation_getRateAndItsTangent

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine source_stochastic_phase_nucleation_init(fileUnit)
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
   SOURCE_stochastic_phase_nucleation_label, &
   SOURCE_stochastic_phase_nucleation_ID, &
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

 write(6,'(/,a)')   ' <<<+-  source_'//SOURCE_stochastic_phase_nucleation_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
 maxNinstance = int(count(phase_source == SOURCE_stochastic_phase_nucleation_ID),pInt)
 if (maxNinstance == 0_pInt) return
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(source_stochastic_phase_nucleation_offset(material_Nphase), source=0_pInt)
 allocate(source_stochastic_phase_nucleation_instance(material_Nphase), source=0_pInt)
 do phase = 1, material_Nphase
   source_stochastic_phase_nucleation_instance(phase) = count(phase_source(:,1:phase) == source_stochastic_phase_nucleation_ID)
   do source = 1, phase_Nsources(phase)
     if (phase_source(source,phase) == source_stochastic_phase_nucleation_ID) &
       source_stochastic_phase_nucleation_offset(phase) = source
   enddo    
 enddo
   
 allocate(source_stochastic_phase_nucleation_sizePostResults(maxNinstance),                     source=0_pInt)
 allocate(source_stochastic_phase_nucleation_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(source_stochastic_phase_nucleation_output(maxval(phase_Noutput),maxNinstance))
          source_stochastic_phase_nucleation_output = ''
 allocate(source_stochastic_phase_nucleation_Noutput(maxNinstance),                                source=0_pInt) 
 allocate(source_stochastic_phase_nucleation_probability(maxNinstance),                         source=0.0_pReal) 
 allocate(source_stochastic_phase_nucleation_sourcestrength(maxNinstance),                      source=0.0_pReal) 

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

   if (phase > 0_pInt ) then; if (any(phase_source(:,phase) == SOURCE_stochastic_phase_nucleation_ID)) then ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran

     instance = source_stochastic_phase_nucleation_instance(phase)                                          ! which instance of my vacancy is present phase
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key
     select case(tag)
       case ('nucleation_probability')
         source_stochastic_phase_nucleation_probability(instance) = IO_floatValue(line,chunkPos,2_pInt)

       case ('nucleation_sourcestrength')
         source_stochastic_phase_nucleation_sourcestrength(instance) = IO_floatValue(line,chunkPos,2_pInt)

     end select
   endif; endif
 enddo parsingFile
 
 initializeInstances: do phase = 1_pInt, material_Nphase
   if (any(phase_source(:,phase) == SOURCE_stochastic_phase_nucleation_ID)) then
     NofMyPhase=count(material_phase==phase)
     instance = source_stochastic_phase_nucleation_instance(phase)
     sourceOffset = source_stochastic_phase_nucleation_offset(phase)

     sizeDotState              =   2_pInt
     sizeDeltaState            =   2_pInt
     sizeState                 =   2_pInt
     sourceState(phase)%p(sourceOffset)%sizeState =       sizeState
     sourceState(phase)%p(sourceOffset)%sizeDotState =    sizeDotState
     sourceState(phase)%p(sourceOffset)%sizeDeltaState =  sizeDeltaState
     sourceState(phase)%p(sourceOffset)%sizePostResults = source_stochastic_phase_nucleation_sizePostResults(instance)
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
end subroutine source_stochastic_phase_nucleation_init

!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine source_stochastic_phase_nucleation_deltaState(ipc, ip, el)
 use prec, only: &
   dEq0
 use material, only: &
   phaseAt, phasememberAt, &
   sourceState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 integer(pInt) :: &
   phase, constituent, sourceOffset, instance
 real(pReal) :: &
   randNo   

 phase       = phaseAt(ipc,ip,el)
 constituent = phasememberAt(ipc,ip,el)
 sourceOffset = source_stochastic_phase_nucleation_offset(phase)
 instance = source_stochastic_phase_nucleation_instance(phase)
 
 sourceState(phase)%p(sourceOffset)%deltaState(:,constituent) = 0.0_pReal
 
 call random_number(randNo)
 sourceState(phase)%p(sourceOffset)%deltaState(1,constituent) = &
   randNo - sourceState(phase)%p(sourceOffset)%state(1,constituent)
   
 if (      (   sourceState(phase)%p(sourceOffset)%state(1,constituent) &
             > source_stochastic_phase_nucleation_probability(instance)) &
     .and. (dEq0(sourceState(phase)%p(sourceOffset)%State(2,constituent)))) &
   
   sourceState(phase)%p(sourceOffset)%deltaState(2,constituent) = 1.0_pReal

end subroutine source_stochastic_phase_nucleation_deltaState

!--------------------------------------------------------------------------------------------------
!> @brief returns local nucleation source
!--------------------------------------------------------------------------------------------------
subroutine source_stochastic_phase_nucleation_getRateAndItsTangent(TDot, dTDOT_dT, Phi, ipc, ip, el)
 use prec, only: &
   dNeq0
 use material, only: &
   phaseAt, phasememberAt, &
   sourceState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), intent(in) :: &
   Phi
 real(pReal), intent(out) :: &
   TDot, dTDOT_dT
 integer(pInt) :: &
   instance, phase, constituent, sourceOffset 

 phase = phaseAt(ipc,ip,el)
 constituent = phasememberAt(ipc,ip,el)
 instance = source_stochastic_phase_nucleation_instance(phase)
 sourceOffset = source_stochastic_phase_nucleation_offset(phase)
 
 TDot = 0.0_pReal
 dTDot_dT = 0.0_pReal
 if (dNeq0(sourceState(phase)%p(sourceOffset)%State(2,constituent))) then
   TDot = source_stochastic_phase_nucleation_sourcestrength(instance)*(1.0_pReal - Phi)
   dTDot_dT = -source_stochastic_phase_nucleation_sourcestrength(instance)
 endif  
 
end subroutine source_stochastic_phase_nucleation_getRateAndItsTangent
 
end module source_stochastic_phase_nucleation
