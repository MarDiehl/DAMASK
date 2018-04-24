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

 enum, bind(c)
   enumerator :: undefined_ID, &
                 chemicalFE_ID, &
                 chemicalPot_ID, &
                 chemicalConc_ID
 end enum

 type, private :: tParameters                                                                         !< container type for internal constitutive parameters
   integer(kind(undefined_ID)), dimension(:), allocatable :: & 
     outputID
   real(pReal), dimension(:,:), allocatable :: &
     InteractionEnergy
   real(pReal), dimension(:  ), allocatable :: &
     Mobility, &
     SolutionEnergy, &
     InitialConc
   real(pReal), dimension(:,:), pointer     :: &  
     chemicalConc0
 end type

 type(tParameters), dimension(:), allocatable, private :: param                                       !< containers of constitutive parameters (len Ninstance)

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
   chemConcMapping, &
   chemicalConc, &
   phasememberAt, &
   phase_chemicalFE, &
   phase_chemicalFEInstance, &
   phase_Noutput, &
   phase_maxNcomponents, &
   CHEMICALFE_THERMODYNAMIC_label, &
   CHEMICALFE_THERMODYNAMIC_ID, &
   material_phase, &
   chemicalState, &
   phase_Ncomponents, &
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
 real(pReal), dimension(:,:), allocatable :: &
   tempMobility, &
   tempSolutionEnergy, &
   tempInitialConc
 real(pReal), dimension(:,:,:), allocatable :: &
   tempInteractionEnergy
 integer(kind(undefined_ID)), dimension(:,:), allocatable :: & 
   tempOutputID

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
 allocate(chemicalFE_thermodynamic_Noutput(maxNinstance),                       source=0_pInt)
 
 allocate(param(maxNinstance))

 allocate(tempOutputID         (maxval(phase_Noutput)  ,maxNinstance))
 allocate(tempMobility         (phase_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(tempSolutionEnergy   (phase_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(tempInitialConc      (phase_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(tempInteractionEnergy(phase_maxNcomponents, &
                                phase_maxNcomponents,maxNinstance),source=0.0_pReal)

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
             chemicalFE_thermodynamic_output(chemicalFE_thermodynamic_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))
             tempOutputID(chemicalFE_thermodynamic_Noutput(instance),instance) = chemicalFE_ID
           case ('chemicalpot')
             chemicalFE_thermodynamic_Noutput(instance) = chemicalFE_thermodynamic_Noutput(instance) + 1_pInt
             chemicalFE_thermodynamic_output(chemicalFE_thermodynamic_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))
             tempOutputID(chemicalFE_thermodynamic_Noutput(instance),instance) = chemicalPot_ID
           case ('chemicalconc')
             chemicalFE_thermodynamic_Noutput(instance) = chemicalFE_thermodynamic_Noutput(instance) + 1_pInt
             chemicalFE_thermodynamic_output(chemicalFE_thermodynamic_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))
             tempOutputID(chemicalFE_thermodynamic_Noutput(instance),instance) = chemicalConc_ID

           case default

         end select

!--------------------------------------------------------------------------------------------------
! parameters depending on number of components
       case ('component_mobility')
         if (chunkPos(1) /= phase_Ncomponents(phase) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_THERMODYNAMIC_label//')')
         do j = 1_pInt, phase_Ncomponents(phase)
           tempMobility(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
       case ('component_solutione')
         if (chunkPos(1) /= phase_Ncomponents(phase) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_THERMODYNAMIC_label//')')
         do j = 1_pInt, phase_Ncomponents(phase)
           tempSolutionEnergy(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
       case ('component_initialconc')
         if (chunkPos(1) /= phase_Ncomponents(phase) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_THERMODYNAMIC_label//')')
         do j = 1_pInt, phase_Ncomponents(phase)
           tempInitialConc(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
       case ('component_interactione')
         if (chunkPos(1) /= phase_Ncomponents(phase)* &
                            (phase_Ncomponents(phase) + 1_pInt) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_THERMODYNAMIC_label//')')
         o = 0_pInt
         do j = 1_pInt, phase_Ncomponents(phase)
           do k = j, phase_Ncomponents(phase)
             o = o + 1_pInt
             tempInteractionEnergy(j,k,instance) = IO_floatValue(line,chunkPos,1_pInt+o)
             tempInteractionEnergy(k,j,instance) = IO_floatValue(line,chunkPos,1_pInt+o)
           enddo
         enddo

       case default

     end select
   endif; endif
 enddo parsingFile

 sanityChecks: do phase = 1_pInt, size(phase_chemicalFE)
   myPhase: if (phase_chemicalFE(phase) == CHEMICALFE_thermodynamic_ID) then
     instance = phase_chemicalFEInstance(phase)
     allocate(param(instance)%OutputID         (chemicalFE_thermodynamic_Noutput(instance)))
     param(instance)%OutputID = tempOutputID(1:chemicalFE_thermodynamic_Noutput(instance),instance)
     allocate(param(instance)%Mobility         (phase_Ncomponents(phase)))
     allocate(param(instance)%SolutionEnergy   (phase_Ncomponents(phase)))
     allocate(param(instance)%InitialConc      (phase_Ncomponents(phase)))
     allocate(param(instance)%InteractionEnergy(phase_Ncomponents(phase), &
                                                phase_Ncomponents(phase)))
     param(instance)%Mobility          = tempMobility         (1:phase_Ncomponents(phase),instance)
     param(instance)%SolutionEnergy    = tempSolutionEnergy   (1:phase_Ncomponents(phase),instance)
     param(instance)%InitialConc       = tempInitialConc      (1:phase_Ncomponents(phase),instance)
     param(instance)%InteractionEnergy = tempInteractionEnergy(1:phase_Ncomponents(phase), &
                                                               1:phase_Ncomponents(phase),instance)
     if (any(param(instance)%Mobility < 0.0_pReal)) &
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
       select case(param(instance)%outputID(o))
         case(chemicalFE_ID)
           mySize = 1_pInt
         case(chemicalPot_ID, &
              chemicalConc_ID &
              )
           mySize = phase_Ncomponents(phase)
         case default
       end select

       outputFound: if (mySize > 0_pInt) then
         chemicalFE_thermodynamic_sizePostResult(o,instance) = mySize
         chemicalFE_thermodynamic_sizePostResults(instance)  = chemicalFE_thermodynamic_sizePostResults(instance) + mySize
       endif outputFound
     enddo outputsLoop
!--------------------------------------------------------------------------------------------------
! allocate state arrays
     sizeState = phase_Ncomponents(phase)
     sizeDotState = 0_pInt
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
     chemConcMapping(phase)%p => phasememberAt
     do j = 1_pInt, phase_Ncomponents(phase)
       chemicalState(phase)%state0   (j,1:NipcMyPhase) = param(instance)%InitialConc(j)
       chemicalState(phase)%subState0(j,1:NipcMyPhase) = param(instance)%InitialConc(j)
       chemicalState(phase)%state    (j,1:NipcMyPhase) = param(instance)%InitialConc(j)
       chemicalConc(j,phase)%p => chemicalState(phase)%state(j,1:NipcMyPhase)
       param(instance)%chemicalConc0 => chemicalState(phase)%state0
     enddo
   endif myPhase2
 enddo initializeInstances

end subroutine chemicalFE_thermodynamic_init


!--------------------------------------------------------------------------------------------------
!> @brief returns the chemical energy for a given instance of this model
!--------------------------------------------------------------------------------------------------
function chemicalFE_thermodynamic_getEnergy(ipc,ip,el)
 use material, only: &
   chemicalConc, &
   chemConcMapping, &
   material_phase, &
   material_homog, &
   temperature, &
   thermalMapping, &
   phase_chemicalFEInstance, &
   phase_Ncomponents

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal) :: &
   T
 real(pReal) :: &
   chemicalFE_thermodynamic_getEnergy
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))) :: &
   conc
 integer(pInt) :: &
   cpI, cpJ, &
   phase, &
   instance    

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 T = temperature(material_homog(ip,el))%p(thermalMapping(material_homog(ip,el))%p(ip,el))
 conc = 0.0_pReal
 do cpI = 1_pInt, phase_Ncomponents(phase)
   conc(cpI) = chemicalConc(cpI,phase)%p(chemConcMapping(phase)%p(ipc, ip, el))
 enddo
 chemicalFE_thermodynamic_getEnergy = &
   kB*T*(1.0_pReal - sum(conc(1:phase_Ncomponents(phase))))* &
   log(1.0_pReal - sum(conc(1:phase_Ncomponents(phase))))
 do cpI = 1_pInt, phase_Ncomponents(phase)
   chemicalFE_thermodynamic_getEnergy = &
     chemicalFE_thermodynamic_getEnergy + &
     param(instance)%SolutionEnergy(cpI)*conc(cpI) + &
     kB*T*conc(cpI)*log(conc(cpI))
   do cpJ = 1_pInt, phase_Ncomponents(phase)  
     chemicalFE_thermodynamic_getEnergy = &
       chemicalFE_thermodynamic_getEnergy + &
       param(instance)%InteractionEnergy(cpI,cpJ)*conc(cpI)*conc(cpJ)
   enddo
 enddo

end function chemicalFE_thermodynamic_getEnergy

!--------------------------------------------------------------------------------------------------
!> @brief returns the component chemical potential for a given instance of this model
!--------------------------------------------------------------------------------------------------
function chemicalFE_thermodynamic_getChemPot(ipc,ip,el)
 use material, only: &
   material_phase, &
   chemicalConc, &
   chemConcMapping, &
   material_homog, &
   temperature, &
   thermalMapping, &
   phase_chemicalFEInstance, &
   phase_Ncomponents

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))) :: &
   conc
 real(pReal) :: &
   T
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))) :: &
   chemicalFE_thermodynamic_getChemPot
 integer(pInt) :: &
   cpI, cpJ, &
   phase, &
   instance    

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 conc = 0.0_pReal
 do cpI = 1_pInt, phase_Ncomponents(phase)
   conc(cpI) = chemicalConc(cpI,phase)%p(chemConcMapping(phase)%p(ipc, ip, el))
 enddo
 T = temperature(material_homog(ip,el))%p(thermalMapping(material_homog(ip,el))%p(ip,el))
 chemicalFE_thermodynamic_getChemPot = 0.0_pReal
 do cpI = 1_pInt, phase_Ncomponents(phase)
   chemicalFE_thermodynamic_getChemPot(cpI) = &
     param(instance)%SolutionEnergy(cpI) + &
     kB*T* &
     log(conc(cpI)/(1.0_pReal - sum(conc(1:phase_Ncomponents(phase)))))
   do cpJ = 1_pInt, phase_Ncomponents(phase)
     chemicalFE_thermodynamic_getChemPot(cpI) = &
       chemicalFE_thermodynamic_getChemPot(cpI) + &
       param(instance)%InteractionEnergy(cpI,cpJ)*conc(cpJ)
   enddo
 enddo

end function chemicalFE_thermodynamic_getChemPot


!--------------------------------------------------------------------------------------------------
!> @brief returns the component concentration for a given instance of this model
!--------------------------------------------------------------------------------------------------
function chemicalFE_thermodynamic_getConc(chempot,ipc,ip,el)
 use material, only: &
   material_phase, &
   chemConcMapping, &
   material_homog, &
   temperature, &
   thermalMapping, &
   phase_chemicalFEInstance, &
   phase_Ncomponents

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))), intent(in) :: &
   chempot
 real(pReal) :: &
   T
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))) :: &
   conc0, &
   chemicalFE_thermodynamic_getConc
 integer(pInt) :: &
   cpI, cpJ, &
   phase, &
   instance    

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 T = temperature(material_homog(ip,el))%p(thermalMapping(material_homog(ip,el))%p(ip,el))
 chemicalFE_thermodynamic_getConc = 0.0_pReal
 do cpI = 1_pInt, phase_Ncomponents(phase)
   conc0(cpI) = param(instance)%chemicalConc0(cpI,chemConcMapping(phase)%p(ipc, ip, el))
 enddo
 do cpI = 1_pInt, phase_Ncomponents(phase)
   chemicalFE_thermodynamic_getConc(cpI) = &
     chempot(cpI) - param(instance)%SolutionEnergy(cpI)
   do cpJ = 1_pInt, phase_Ncomponents(phase)
     chemicalFE_thermodynamic_getConc(cpI) = &
       chemicalFE_thermodynamic_getConc(cpI) - &
       param(instance)%InteractionEnergy(cpI,cpJ)*conc0(cpJ)
   enddo
 enddo
 chemicalFE_thermodynamic_getConc(1:phase_Ncomponents(phase)) = &
   exp(chemicalFE_thermodynamic_getConc(1:phase_Ncomponents(phase))/(kB*T))

 chemicalFE_thermodynamic_getConc(1:phase_Ncomponents(phase)) = &
   chemicalFE_thermodynamic_getConc(1:phase_Ncomponents(phase))/ &
   (1.0_pReal + &
    sum(chemicalFE_thermodynamic_getConc(1:phase_Ncomponents(phase))))
 
end function chemicalFE_thermodynamic_getConc


!--------------------------------------------------------------------------------------------------
!> @brief returns the component concentration for a given instance of this model
!--------------------------------------------------------------------------------------------------
subroutine chemicalFE_thermodynamic_putConc(chempot,ipc,ip,el)
 use material, only: &
   material_phase, &
   phase_chemicalFEInstance, &
   phase_Ncomponents, &
   chemConcMapping, &
   chemicalConc

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))), intent(in) :: &
   chempot
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))) :: &
   conc
 integer(pInt) :: &
   cp, &
   phase, &
   instance    

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 conc = chemicalFE_thermodynamic_getConc(chempot,ipc,ip,el)
 do cp = 1_pInt, phase_Ncomponents(phase)
   chemicalConc(cp,phase)%p(chemConcMapping(phase)%p(ipc,ip,el)) = &
     conc(cp)
 enddo
 
end subroutine chemicalFE_thermodynamic_putConc


!--------------------------------------------------------------------------------------------------
!> @brief returns the component concentration tangent for a given instance of this model
!--------------------------------------------------------------------------------------------------
subroutine chemicalFE_thermodynamic_getConcTangent(dConcdChemPot,dConcdGradC,chempot,ipc,ip,el)
 use material, only: &
   material_phase, &
   chemConcMapping, &
   material_homog, &
   temperature, &
   thermalMapping, &
   phase_chemicalFEInstance, &
   phase_Ncomponents

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))), intent(in) :: &
   chempot
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el)), &
                        phase_Ncomponents(material_phase(ipc,ip,el))), intent(out) :: &
   dConcdChemPot, &
   dConcdGradC
 real(pReal) :: &
   T
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))) :: &
   conc, conc0
 integer(pInt) :: &
   cpI, cpJ, &
   phase, &
   instance    

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 T = temperature(material_homog(ip,el))%p(thermalMapping(material_homog(ip,el))%p(ip,el))
 do cpI = 1_pInt, phase_Ncomponents(phase)
   conc0(cpI) = param(instance)%chemicalConc0(cpI,chemConcMapping(phase)%p(ipc, ip, el))
 enddo
 conc = chemicalFE_thermodynamic_getConc(chempot,ipc,ip,el)
 
 dConcdChemPot = 0.0_pReal
 dConcdGradC = 0.0_pReal
 do cpI = 1_pInt, phase_Ncomponents(phase)
   do cpJ = 1_pInt, phase_Ncomponents(phase)
     dConcdChemPot(cpI,cpJ) = conc(cpI) - conc(cpI)*conc(cpJ)
   enddo
 enddo

end subroutine chemicalFE_thermodynamic_getConcTangent


!--------------------------------------------------------------------------------------------------
!> @brief returns the component mobility for a given instance of this model
!--------------------------------------------------------------------------------------------------
function chemicalFE_thermodynamic_getMobility(ipc,ip,el)
 use material, only: &
   material_phase, &
   phase_chemicalFEInstance, &
   phase_Ncomponents

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))) :: &
   chemicalFE_thermodynamic_getMobility
 integer(pInt) :: &
   cp, &
   phase, &
   instance    

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 chemicalFE_thermodynamic_getMobility = 0.0_pReal
 do cp = 1_pInt, phase_Ncomponents(phase)
   chemicalFE_thermodynamic_getMobility(cp) = &
     param(instance)%Mobility(cp)
 enddo

end function chemicalFE_thermodynamic_getMobility


!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function chemicalFE_thermodynamic_postResults(ipc,ip,el)
 use material, only: &
   phase_Ncomponents, &
   material_phase, &
   chemicalConc, &
   chemConcMapping, &
   material_homog, &
   temperature, &
   thermalMapping, &
   phase_chemicalFEInstance

 implicit none
 integer(pInt),                                   intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element                                                                                        !< microstructure state
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))) :: &
   conc
 real(pReal) :: &
   T

 real(pReal), dimension(chemicalFE_thermodynamic_sizePostResults(phase_chemicalFEInstance(material_phase(ipc,ip,el)))) :: &
   chemicalFE_thermodynamic_postResults
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))) :: &
   tempPerComponent

 integer(pInt) :: &
   phase, &
   instance, &
   o,c,cp

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 do cp = 1_pInt, phase_Ncomponents(phase)
   conc(cp) = chemicalConc(cp,phase)%p(chemConcMapping(phase)%p(ipc,ip,el))
 enddo
 T = temperature(material_homog(ip,el))%p(thermalMapping(material_homog(ip,el))%p(ip,el))
 chemicalFE_thermodynamic_postResults = 0.0_pReal
 c = 0_pInt

 outputsLoop: do o = 1_pInt,chemicalFE_thermodynamic_Noutput(instance)
   select case(param(instance)%outputID(o))
     case (chemicalFE_ID)
       chemicalFE_thermodynamic_postResults(c+1_pInt) = chemicalFE_thermodynamic_getEnergy(ipc,ip,el)
       c = c + 1_pInt

     case (chemicalPot_ID)
       tempPerComponent = chemicalFE_thermodynamic_getChemPot(ipc,ip,el)
       chemicalFE_thermodynamic_postResults(c+1_pInt:c+phase_Ncomponents(phase)) = &
         tempPerComponent(1_pInt:phase_Ncomponents(phase))
       c = c + phase_Ncomponents(phase)

     case (chemicalConc_ID)
       chemicalFE_thermodynamic_postResults(c+1_pInt:c+phase_Ncomponents(phase)) = &
         conc(1_pInt:phase_Ncomponents(phase))
       c = c + phase_Ncomponents(phase)

   end select
 enddo outputsLoop

end function chemicalFE_thermodynamic_postResults

end module chemicalFE_thermodynamic
