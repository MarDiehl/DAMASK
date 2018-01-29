!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for thermodynamic chemical free energy
!--------------------------------------------------------------------------------------------------
module chemicalFE_thermodynamic
 use prec, only: &
   pReal,&
   pInt, &
   p_vec

 implicit none
 private
 real(pReal),                                          parameter,           private :: &
   kB = 1.3806488e-23_pReal                                                                          !< Boltzmann constant in J/Kelvin

 integer(pInt),                       dimension(:),     allocatable,         public, protected :: &
   chemicalFE_thermodynamic_sizePostResults                                                       !< cumulative size of post results

 integer(pInt),                       dimension(:,:),   allocatable, target, public :: &
   chemicalFE_thermodynamic_sizePostResult                                                        !< size of each post result output

 character(len=64),                   dimension(:,:),   allocatable, target, public :: &
   chemicalFE_thermodynamic_output                                                                !< name of each post result output

 integer(pInt),                       dimension(:),     allocatable, target, public :: &
   chemicalFE_thermodynamic_Noutput                                                               !< number of outputs per instance of this constitution

 integer(pInt),                       dimension(:),     allocatable,         private, protected :: &
   chemicalFE_thermodynamic_Ncomponents                                                           !< no. of chemical components

 real(pReal),                         dimension(:,:),   allocatable,         private :: &
   chemicalFE_thermodynamic_Mobility, &                                                         
   chemicalFE_thermodynamic_SolutionEnergy                                                        
 real(pReal),                         dimension(:,:,:), allocatable,         private :: &
   chemicalFE_thermodynamic_InteractionEnergy

 enum, bind(c)
   enumerator :: undefined_ID, &
                 chemicalFE_ID, &
                 chemicalPot_ID, &
                 chemicalConc_ID
 end enum
 integer(kind(undefined_ID)),         dimension(:,:),   allocatable,          private :: &
   chemicalFE_thermodynamic_outputID                                                              !< ID of each post result output

 public :: &
   chemicalFE_thermodynamic_init, &
   chemicalFE_thermodynamic_getEnergy, &
   chemicalFE_thermodynamic_getChemPot, &
   chemicalFE_thermodynamic_getConc, &
   chemicalFE_thermodynamic_putConc, &
   chemicalFE_thermodynamic_getConcTangent, &
   chemicalFE_thermodynamic_getMobility, &
   chemicalFE_thermodynamic_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine chemicalFE_thermodynamic_init(fileUnit)
#ifdef __GFORTRAN__
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
   chemicalMapping, &
   chemicalConc, &
   phasememberAt, &
   phase_chemicalFE, &
   phase_chemicalFEInstance, &
   phase_Noutput, &
   material_maxNcomponents, &
   CHEMICALFE_THERMODYNAMIC_label, &
   CHEMICALFE_THERMODYNAMIC_ID, &
   material_phase, &
   chemicalState, &
   MATERIAL_partPhase
 use numerics,only: &
   numerics_integrator

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: &
   maxNinstance, &
   instance,phase,j,k,o, &
   NipcMyPhase, &
   mySize=0_pInt,sizeState,sizeDotState, sizeDeltaState
 character(len=65536) :: &
   tag  = '', &
   line = ''

 write(6,'(/,a)')   ' <<<+-  chemicalFE_'//CHEMICALFE_THERMODYNAMIC_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = int(count(phase_chemicalFE == CHEMICALFE_THERMODYNAMIC_ID),pInt)
 if (maxNinstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance


 allocate(chemicalFE_thermodynamic_sizePostResults(maxNinstance),                   source=0_pInt)
 allocate(chemicalFE_thermodynamic_sizePostResult(maxval(phase_Noutput),maxNinstance), &
                                                                                 source=0_pInt)
 allocate(chemicalFE_thermodynamic_output(maxval(phase_Noutput),maxNinstance))
          chemicalFE_thermodynamic_output               = ''
 allocate(chemicalFE_thermodynamic_outputID(maxval(phase_Noutput),maxNinstance),source=undefined_ID)
 allocate(chemicalFE_thermodynamic_Noutput(maxNinstance),                       source=0_pInt)
 allocate(chemicalFE_thermodynamic_Ncomponents                              (maxNinstance),source=  0_pInt )
 allocate(chemicalFE_thermodynamic_Mobility         (material_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(chemicalFE_thermodynamic_SolutionEnergy   (material_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(chemicalFE_thermodynamic_InteractionEnergy(material_maxNcomponents, &
                                                     material_maxNcomponents,maxNinstance),source=0.0_pReal)

 rewind(fileUnit)
 phase = 0_pInt
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= material_partPhase)         ! wind forward to <phase>
   line = IO_read(fileUnit)
 enddo

 parsingFile: do while (trim(line) /= IO_EOF)                                                       ! read through sections of phase part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit
   endif
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next phase
     phase = phase + 1_pInt                                                                         ! advance phase section counter
     cycle                                                                                          ! skip to next line
   endif
   if (phase > 0_pInt ) then; if (phase_chemicalFE(phase) == CHEMICALFE_THERMODYNAMIC_ID) then      ! one of my phases. Do not short-circuit here (.and. between if-statements), it's not safe in Fortran
     instance = phase_chemicalFEInstance(phase)                                                     ! which instance of my chemicalFE is present phase
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case ('chemicalfe')
             chemicalFE_thermodynamic_Noutput(instance) = chemicalFE_thermodynamic_Noutput(instance) + 1_pInt
             chemicalFE_thermodynamic_outputID(chemicalFE_thermodynamic_Noutput(instance),instance) = chemicalFE_ID
             chemicalFE_thermodynamic_output(chemicalFE_thermodynamic_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('chemicalpot')
             chemicalFE_thermodynamic_Noutput(instance) = chemicalFE_thermodynamic_Noutput(instance) + 1_pInt
             chemicalFE_thermodynamic_outputID(chemicalFE_thermodynamic_Noutput(instance),instance) = chemicalPot_ID
             chemicalFE_thermodynamic_output(chemicalFE_thermodynamic_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('chemicalconc')
             chemicalFE_thermodynamic_Noutput(instance) = chemicalFE_thermodynamic_Noutput(instance) + 1_pInt
             chemicalFE_thermodynamic_outputID(chemicalFE_thermodynamic_Noutput(instance),instance) = chemicalConc_ID
             chemicalFE_thermodynamic_output(chemicalFE_thermodynamic_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))

           case default

         end select

       case ('ncomponents')
         chemicalFE_thermodynamic_Ncomponents(instance) = IO_intValue(line,chunkPos,2_pInt)

!--------------------------------------------------------------------------------------------------
! parameters depending on number of components
       case ('component_mobility')
         if (chunkPos(1) /= chemicalFE_thermodynamic_Ncomponents(instance) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_THERMODYNAMIC_label//')')
         do j = 1_pInt, chemicalFE_thermodynamic_Ncomponents(instance)
           chemicalFE_thermodynamic_Mobility(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
       case ('component_solutione')
         if (chunkPos(1) /= chemicalFE_thermodynamic_Ncomponents(instance) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_THERMODYNAMIC_label//')')
         do j = 1_pInt, chemicalFE_thermodynamic_Ncomponents(instance)
           chemicalFE_thermodynamic_SolutionEnergy(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
       case ('component_interactione')
         if (chunkPos(1) /= chemicalFE_thermodynamic_Ncomponents(instance)* &
                            (chemicalFE_thermodynamic_Ncomponents(instance) + 1_pInt) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_THERMODYNAMIC_label//')')
         o = 0_pInt
         do j = 1_pInt, chemicalFE_thermodynamic_Ncomponents(instance)
           do k = j, chemicalFE_thermodynamic_Ncomponents(instance)
             o = o + 1_pInt
             chemicalFE_thermodynamic_InteractionEnergy(j,k,instance) = IO_floatValue(line,chunkPos,1_pInt+o)
             chemicalFE_thermodynamic_InteractionEnergy(k,j,instance) = IO_floatValue(line,chunkPos,1_pInt+o)
           enddo
         enddo

       case default

     end select
   endif; endif
 enddo parsingFile

 sanityChecks: do phase = 1_pInt, size(phase_chemicalFE)
   myPhase: if (phase_chemicalFE(phase) == CHEMICALFE_thermodynamic_ID) then
     instance = phase_chemicalFEInstance(phase)
     if (any(chemicalFE_thermodynamic_Mobility(:,instance) < 0.0_pReal)) &
       call IO_error(211_pInt,el=instance,ext_msg='mobility ('//CHEMICALFE_THERMODYNAMIC_label//')')
   endif myPhase
 enddo sanityChecks

 initializeInstances: do phase = 1_pInt, size(phase_chemicalFE)                                     ! loop through all phases in material.config
   myPhase2: if (phase_chemicalFE(phase) == CHEMICALFE_thermodynamic_ID) then                       ! only consider my phase
     NipcMyPhase = count(material_phase == phase)                                                   ! number of IPCs containing my phase
     instance = phase_chemicalFEInstance(phase)                                                     ! which instance of my phase

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,chemicalFE_thermodynamic_Noutput(instance)
       select case(chemicalFE_thermodynamic_outputID(o,instance))
         case(chemicalFE_ID)
           mySize = 1_pInt
         case(chemicalPot_ID, &
              chemicalConc_ID &
              )
           mySize = chemicalFE_thermodynamic_Ncomponents(instance)
         case default
       end select

       outputFound: if (mySize > 0_pInt) then
         chemicalFE_thermodynamic_sizePostResult(o,instance) = mySize
         chemicalFE_thermodynamic_sizePostResults(instance)  = chemicalFE_thermodynamic_sizePostResults(instance) + mySize
       endif outputFound
     enddo outputsLoop
!--------------------------------------------------------------------------------------------------
! allocate state arrays
     sizeState = 0_pInt
     sizeDotState = sizeState
     sizeDeltaState = 0_pInt
     chemicalState(phase)%sizeState = sizeState
     chemicalState(phase)%sizeDotState = sizeDotState
     chemicalState(phase)%sizeDeltaState = sizeDeltaState
     chemicalState(phase)%sizePostResults = chemicalFE_thermodynamic_sizePostResults(instance)
     allocate(chemicalState(phase)%aTolState          (   sizeState),             source=0.0_pReal)
     allocate(chemicalState(phase)%state0             (   sizeState,NipcMyPhase), source=0.0_pReal)
     allocate(chemicalState(phase)%partionedState0    (   sizeState,NipcMyPhase), source=0.0_pReal)
     allocate(chemicalState(phase)%subState0          (   sizeState,NipcMyPhase), source=0.0_pReal)
     allocate(chemicalState(phase)%state              (   sizeState,NipcMyPhase), source=0.0_pReal)
     allocate(chemicalState(phase)%dotState           (sizeDotState,NipcMyPhase), source=0.0_pReal)
     allocate(chemicalState(phase)%deltaState       (sizeDeltaState,NipcMyPhase), source=0.0_pReal)
     if (any(numerics_integrator == 1_pInt)) then
       allocate(chemicalState(phase)%previousDotState (sizeDotState,NipcMyPhase),source=0.0_pReal)
       allocate(chemicalState(phase)%previousDotState2(sizeDotState,NipcMyPhase),source=0.0_pReal)
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(chemicalState(phase)%RK4dotState      (sizeDotState,NipcMyPhase), source=0.0_pReal)
     if (any(numerics_integrator == 5_pInt)) &
       allocate(chemicalState(phase)%RKCK45dotState (6,sizeDotState,NipcMyPhase), source=0.0_pReal)
     chemicalMapping(phase)%p => phasememberAt
     do j = 1_pInt, chemicalFE_thermodynamic_Ncomponents(instance)
       allocate(chemicalConc(j,phase)%p(NipcMyPhase), source=0.0_pReal)
     enddo
   endif myPhase2
 enddo initializeInstances

end subroutine chemicalFE_thermodynamic_init


!--------------------------------------------------------------------------------------------------
!> @brief returns the chemical energy for a given instance of this model
!--------------------------------------------------------------------------------------------------
function chemicalFE_thermodynamic_getEnergy(conc,T,ipc,ip,el)
 use material, only: &
   material_phase, &
   phase_chemicalFEInstance, &
   material_maxNcomponents

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(material_maxNcomponents), intent(in) :: &
   conc
 real(pReal), intent(in) :: &
   T
 real(pReal) :: &
   chemicalFE_thermodynamic_getEnergy
 integer(pInt) :: &
   cpI, cpJ, &
   phase, &
   instance    

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 chemicalFE_thermodynamic_getEnergy = &
   kB*T*(1.0_pReal - sum(conc(1:chemicalFE_thermodynamic_Ncomponents(instance))))* &
   log(1.0_pReal - sum(conc(1:chemicalFE_thermodynamic_Ncomponents(instance))))
 do cpI = 1_pInt, chemicalFE_thermodynamic_Ncomponents(instance)
   chemicalFE_thermodynamic_getEnergy = &
     chemicalFE_thermodynamic_getEnergy + &
     chemicalFE_thermodynamic_SolutionEnergy(cpI,instance)*conc(cpI) + &
     kB*T*conc(cpI)*log(conc(cpI))
   do cpJ = 1_pInt, chemicalFE_thermodynamic_Ncomponents(instance)  
     chemicalFE_thermodynamic_getEnergy = &
       chemicalFE_thermodynamic_getEnergy + &
       chemicalFE_thermodynamic_InteractionEnergy(cpI,cpJ,instance)*conc(cpI)*conc(cpJ)
   enddo
 enddo

end function chemicalFE_thermodynamic_getEnergy

!--------------------------------------------------------------------------------------------------
!> @brief returns the component chemical potential for a given instance of this model
!--------------------------------------------------------------------------------------------------
function chemicalFE_thermodynamic_getChemPot(conc,T,ipc,ip,el)
 use material, only: &
   material_phase, &
   phase_chemicalFEInstance, &
   material_maxNcomponents

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(material_maxNcomponents), intent(in) :: &
   conc
 real(pReal), intent(in) :: &
   T
 real(pReal), dimension(material_maxNcomponents) :: &
   chemicalFE_thermodynamic_getChemPot
 integer(pInt) :: &
   cpI, cpJ, &
   phase, &
   instance    

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 chemicalFE_thermodynamic_getChemPot = 0.0_pReal
 do cpI = 1_pInt, chemicalFE_thermodynamic_Ncomponents(instance)
   chemicalFE_thermodynamic_getChemPot(cpI) = &
     chemicalFE_thermodynamic_SolutionEnergy(cpI,instance) + &
     kB*T* &
     log(conc(cpI)/(1.0_pReal - sum(conc(1:chemicalFE_thermodynamic_Ncomponents(instance)))))
   do cpJ = 1_pInt, chemicalFE_thermodynamic_Ncomponents(instance)
     chemicalFE_thermodynamic_getChemPot(cpI) = &
       chemicalFE_thermodynamic_getChemPot(cpI) + &
       chemicalFE_thermodynamic_InteractionEnergy(cpI,cpJ,instance)*conc(cpJ)
   enddo
 enddo

end function chemicalFE_thermodynamic_getChemPot


!--------------------------------------------------------------------------------------------------
!> @brief returns the component concentration for a given instance of this model
!--------------------------------------------------------------------------------------------------
function chemicalFE_thermodynamic_getConc(chempot,conc0,T,ipc,ip,el)
 use material, only: &
   material_phase, &
   phase_chemicalFEInstance, &
   material_maxNcomponents

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(material_maxNcomponents), intent(in) :: &
   chempot, conc0
 real(pReal), intent(in) :: &
   T
 real(pReal), dimension(material_maxNcomponents) :: &
   chemicalFE_thermodynamic_getConc
 integer(pInt) :: &
   cpI, cpJ, &
   phase, &
   instance    

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 chemicalFE_thermodynamic_getConc = 0.0_pReal
 do cpI = 1_pInt, chemicalFE_thermodynamic_Ncomponents(instance)
   chemicalFE_thermodynamic_getConc(cpI) = &
     chempot(cpI) - chemicalFE_thermodynamic_SolutionEnergy(cpI,instance)
   do cpJ = 1_pInt, chemicalFE_thermodynamic_Ncomponents(instance)
     chemicalFE_thermodynamic_getConc(cpI) = &
       chemicalFE_thermodynamic_getConc(cpI) - &
       chemicalFE_thermodynamic_InteractionEnergy(cpI,cpJ,instance)*conc0(cpJ)
   enddo
 enddo
 chemicalFE_thermodynamic_getConc(1:chemicalFE_thermodynamic_Ncomponents(instance)) = &
   exp(chemicalFE_thermodynamic_getConc(1:chemicalFE_thermodynamic_Ncomponents(instance))/(kB*T))

 chemicalFE_thermodynamic_getConc(1:chemicalFE_thermodynamic_Ncomponents(instance)) = &
   chemicalFE_thermodynamic_getConc(1:chemicalFE_thermodynamic_Ncomponents(instance))/ &
   (1.0_pReal + &
    sum(chemicalFE_thermodynamic_getConc(1:chemicalFE_thermodynamic_Ncomponents(instance))))
 
end function chemicalFE_thermodynamic_getConc


!--------------------------------------------------------------------------------------------------
!> @brief returns the component concentration for a given instance of this model
!--------------------------------------------------------------------------------------------------
subroutine chemicalFE_thermodynamic_putConc(chempot,conc0,T,ipc,ip,el)
 use material, only: &
   material_phase, &
   phase_chemicalFEInstance, &
   material_maxNcomponents, &
   chemicalMapping, &
   chemicalConc

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(material_maxNcomponents), intent(in) :: &
   chempot, conc0
 real(pReal), intent(in) :: &
   T
 real(pReal), dimension(material_maxNcomponents) :: &
   conc
 integer(pInt) :: &
   cp, &
   phase, &
   instance    

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 conc = chemicalFE_thermodynamic_getConc(chempot,conc0,T,ipc,ip,el)
 do cp = 1_pInt, chemicalFE_thermodynamic_Ncomponents(instance)
   chemicalConc(cp,phase)%p(chemicalMapping(phase)%p(ipc,ip,el)) = &
     conc(cp)
 enddo
 
end subroutine chemicalFE_thermodynamic_putConc


!--------------------------------------------------------------------------------------------------
!> @brief returns the component concentration tangent for a given instance of this model
!--------------------------------------------------------------------------------------------------
function chemicalFE_thermodynamic_getConcTangent(chempot,conc0,T,ipc,ip,el)
 use material, only: &
   material_phase, &
   phase_chemicalFEInstance, &
   material_maxNcomponents

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(material_maxNcomponents), intent(in) :: &
   chempot, conc0
 real(pReal), intent(in) :: &
   T
 real(pReal), dimension(material_maxNcomponents,material_maxNcomponents) :: &
   chemicalFE_thermodynamic_getConcTangent
 real(pReal), dimension(material_maxNcomponents) :: &
   conc
 integer(pInt) :: &
   cpI, cpJ, &
   phase, &
   instance    

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 conc = chemicalFE_thermodynamic_getConc(chempot,conc0,T,ipc,ip,el)
 
 chemicalFE_thermodynamic_getConcTangent = 0.0_pReal
 do cpI = 1_pInt, chemicalFE_thermodynamic_Ncomponents(instance)
   do cpJ = 1_pInt, chemicalFE_thermodynamic_Ncomponents(instance)
     chemicalFE_thermodynamic_getConcTangent(cpI,cpJ) = &
       conc(cpI) - conc(cpI)*conc(cpJ)
   enddo
 enddo

end function chemicalFE_thermodynamic_getConcTangent


!--------------------------------------------------------------------------------------------------
!> @brief returns the component mobility for a given instance of this model
!--------------------------------------------------------------------------------------------------
function chemicalFE_thermodynamic_getMobility(ipc,ip,el)
 use material, only: &
   material_phase, &
   phase_chemicalFEInstance, &
   material_maxNcomponents

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(material_maxNcomponents) :: &
   chemicalFE_thermodynamic_getMobility
 integer(pInt) :: &
   cp, &
   phase, &
   instance    

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 chemicalFE_thermodynamic_getMobility = 0.0_pReal
 do cp = 1_pInt, chemicalFE_thermodynamic_Ncomponents(instance)
   chemicalFE_thermodynamic_getMobility(cp) = &
     chemicalFE_thermodynamic_Mobility(cp,instance)
 enddo

end function chemicalFE_thermodynamic_getMobility


!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function chemicalFE_thermodynamic_postResults(conc,T,ipc,ip,el)
 use material, only: &
   material_maxNcomponents, &
   material_phase, &
   phase_chemicalFEInstance

 implicit none
 real(pReal), dimension(material_maxNcomponents), intent(in) :: &
   conc
 real(pReal),                                     intent(in) :: &
   T
 integer(pInt),                                   intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element                                                                                        !< microstructure state

 real(pReal), dimension(chemicalFE_thermodynamic_sizePostResults(phase_chemicalFEInstance(material_phase(ipc,ip,el)))) :: &
   chemicalFE_thermodynamic_postResults
 real(pReal), dimension(material_maxNcomponents) :: &
   tempPerComponent

 integer(pInt) :: &
   instance, &
   o,c

 instance = phase_chemicalFEInstance(material_phase(ipc,ip,el))
 chemicalFE_thermodynamic_postResults = 0.0_pReal
 c = 0_pInt

 outputsLoop: do o = 1_pInt,chemicalFE_thermodynamic_Noutput(instance)
   select case(chemicalFE_thermodynamic_outputID(o,instance))
     case (chemicalFE_ID)
       chemicalFE_thermodynamic_postResults(c+1_pInt) = chemicalFE_thermodynamic_getEnergy(conc,T,ipc,ip,el)
       c = c + 1_pInt

     case (chemicalPot_ID)
       tempPerComponent = chemicalFE_thermodynamic_getChemPot(conc,T,ipc,ip,el)
       chemicalFE_thermodynamic_postResults(c+1_pInt:c+chemicalFE_thermodynamic_Ncomponents(instance)) = &
         tempPerComponent(1_pInt:chemicalFE_thermodynamic_Ncomponents(instance))
       c = c + chemicalFE_thermodynamic_Ncomponents(instance)

     case (chemicalConc_ID)
       chemicalFE_thermodynamic_postResults(c+1_pInt:c+chemicalFE_thermodynamic_Ncomponents(instance)) = &
         conc(1_pInt:chemicalFE_thermodynamic_Ncomponents(instance))
       c = c + chemicalFE_thermodynamic_Ncomponents(instance)

   end select
 enddo outputsLoop

end function chemicalFE_thermodynamic_postResults

end module chemicalFE_thermodynamic
