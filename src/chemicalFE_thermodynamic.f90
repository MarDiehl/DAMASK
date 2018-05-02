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
   R = 8.314459848_pReal                                                                          !< Universal gas constant in J/(mol Kelvin)

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
                 chemicalConc_ID
 end enum

 type, private :: tParameters                                                                         !< container type for internal constitutive parameters
   integer(kind(undefined_ID)), dimension(:), allocatable :: & 
     outputID
   integer(pInt) :: &
     maxNIter = 20_pInt
   real(pReal) :: &
     MolarVolume, &
     aTol = 1e-12_pReal
   real(pReal), dimension(:,:), allocatable :: &
     InteractionEnergy
   real(pReal), dimension(:  ), allocatable :: &
     Mobility, &
     SolutionEnergy, &
     InitialAvgConc, &
     InitialDeltaConc, &
     GradientCoeff
 end type

 type(tParameters), dimension(:), allocatable, private :: param                                       !< containers of constitutive parameters (len Ninstance)

 public :: &
   chemicalFE_thermodynamic_init, &
   chemicalFE_thermodynamic_dotState, &
   chemicalFE_thermodynamic_getEnergy, &
   chemicalFE_thermodynamic_calConcandTangent, &
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
   chemicalConc0, &
   chemicalConcRate, &
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
   tempInitialAvgConc, &
   tempInitialDeltaConc, &
   tempGradCoeff
 real(pReal), dimension(:,:,:), allocatable :: &
   tempInteractionEnergy
 real(pReal) :: &
   randNum, &
   InitialConc
 integer(kind(undefined_ID)), dimension(:,:), allocatable :: & 
   tempOutputID

 write(6,'(/,a)')   ' <<<+-  chemicalFE_'//CHEMICALFE_THERMODYNAMIC_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = int(count(phase_chemicalFE == CHEMICALFE_THERMODYNAMIC_ID),pInt)
 if (maxNinstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance


 allocate(chemicalFE_thermodynamic_sizePostResults(maxNinstance),                source=0_pInt)
 allocate(chemicalFE_thermodynamic_sizePostResult(maxval(phase_Noutput),maxNinstance), &
                                                                                 source=0_pInt)
 allocate(chemicalFE_thermodynamic_output(maxval(phase_Noutput),maxNinstance))
          chemicalFE_thermodynamic_output               = ''
 allocate(chemicalFE_thermodynamic_Noutput(maxNinstance),                        source=0_pInt)
 
 allocate(param(maxNinstance))

 allocate(tempOutputID         (maxval(phase_Noutput)  ,maxNinstance))
 allocate(tempMobility         (phase_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(tempSolutionEnergy   (phase_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(tempInitialAvgConc   (phase_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(tempInitialDeltaConc (phase_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(tempGradCoeff        (phase_maxNcomponents,maxNinstance),source=0.0_pReal)
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
           case ('chemicalconc')
             chemicalFE_thermodynamic_Noutput(instance) = chemicalFE_thermodynamic_Noutput(instance) + 1_pInt
             chemicalFE_thermodynamic_output(chemicalFE_thermodynamic_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))
             tempOutputID(chemicalFE_thermodynamic_Noutput(instance),instance) = chemicalConc_ID

           case default

         end select

!--------------------------------------------------------------------------------------------------
! parameters depending on number of components
       case ('component_tolerance')
         param(instance)%aTol = IO_floatValue(line,chunkPos,2_pInt)
         
       case ('component_maxniter')
         param(instance)%maxNIter = IO_intValue(line,chunkPos,2_pInt)
         
       case ('molarvolume')
         param(instance)%MolarVolume = IO_floatValue(line,chunkPos,2_pInt)
         
       case ('component_mobility')
         if (chunkPos(1) /= phase_Ncomponents(phase) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_THERMODYNAMIC_label//')')
         do j = 1_pInt, phase_Ncomponents(phase)
           tempMobility(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
       case ('component_solutione','component_solutionenergy')
         if (chunkPos(1) /= phase_Ncomponents(phase) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_THERMODYNAMIC_label//')')
         do j = 1_pInt, phase_Ncomponents(phase)
           tempSolutionEnergy(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
       case ('component_initialavgconc')
         if (chunkPos(1) /= phase_Ncomponents(phase) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_THERMODYNAMIC_label//')')
         do j = 1_pInt, phase_Ncomponents(phase)
           tempInitialAvgConc(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
       case ('component_initialdeltaconc')
         if (chunkPos(1) /= phase_Ncomponents(phase) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_THERMODYNAMIC_label//')')
         do j = 1_pInt, phase_Ncomponents(phase)
           tempInitialDeltaConc(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
       case ('component_gradcoeff','component_gradientcoeff')
         if (chunkPos(1) /= phase_Ncomponents(phase) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_THERMODYNAMIC_label//')')
         do j = 1_pInt, phase_Ncomponents(phase)
           tempGradCoeff(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
       case ('component_interactione','component_interactionenergy')
         if (chunkPos(1) /= phase_Ncomponents(phase)* &
                            (phase_Ncomponents(phase) + 1_pInt)/2_pInt + 1_pInt) &
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
     allocate(param(instance)%InitialAvgConc   (phase_Ncomponents(phase)))
     allocate(param(instance)%InitialDeltaConc (phase_Ncomponents(phase)))
     allocate(param(instance)%GradientCoeff    (phase_Ncomponents(phase)))
     allocate(param(instance)%InteractionEnergy(phase_Ncomponents(phase), &
                                                phase_Ncomponents(phase)))
     param(instance)%Mobility          = tempMobility         (1:phase_Ncomponents(phase),instance)
     param(instance)%SolutionEnergy    = tempSolutionEnergy   (1:phase_Ncomponents(phase),instance)
     param(instance)%InitialAvgConc    = tempInitialAvgConc   (1:phase_Ncomponents(phase),instance)
     param(instance)%InitialDeltaConc  = tempInitialDeltaConc (1:phase_Ncomponents(phase),instance)
     param(instance)%GradientCoeff     = tempGradCoeff        (1:phase_Ncomponents(phase),instance)
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
         case(chemicalConc_ID)
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
     sizeDotState = phase_Ncomponents(phase)
     sizeDeltaState = 0_pInt
     chemicalState(phase)%sizeState = sizeState
     chemicalState(phase)%sizeDotState = sizeDotState
     chemicalState(phase)%sizeDeltaState = sizeDeltaState
     chemicalState(phase)%sizePostResults = chemicalFE_thermodynamic_sizePostResults(instance)
     allocate(chemicalState(phase)%aTolState          (   sizeState),             source=param(instance)%aTol)
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
       do o = 1_pInt, NipcMyPhase
         call random_number(randNum)
         InitialConc = param(instance)%InitialAvgConc(j) + param(instance)%InitialDeltaConc(j)*(randNum - 0.5_pReal)
         chemicalState(phase)%state0   (j,o) = InitialConc
         chemicalState(phase)%subState0(j,o) = InitialConc
         chemicalState(phase)%state    (j,o) = InitialConc
       enddo  
       chemicalConc (j,phase)%p => chemicalState(phase)%state (j,1:NipcMyPhase)
       chemicalConc0(j,phase)%p => chemicalState(phase)%state0(j,1:NipcMyPhase)
       allocate(chemicalConcRate(j,phase)%p(NipcMyPhase), source=0.0_pReal)
     enddo
   endif myPhase2
 enddo initializeInstances

end subroutine chemicalFE_thermodynamic_init


!--------------------------------------------------------------------------------------------------
!> @brief returns the chemical concentration rates
!--------------------------------------------------------------------------------------------------
subroutine chemicalFE_thermodynamic_dotState(ipc,ip,el)
 use material, only: &
   chemicalState, &
   chemicalConcRate, &
   chemConcMapping, &
   material_phase, &
   phase_Ncomponents
   
 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 integer(pInt) :: &
   cp, &
   phase

 phase = material_phase(ipc,ip,el) 
 do cp = 1_pInt, phase_Ncomponents(phase)
   chemicalState(phase)%dotState(cp,chemConcMapping(phase)%p(ipc,ip,el)) = &
     chemicalConcRate(cp,phase)%p(chemConcMapping(phase)%p(ipc,ip,el))
 enddo

end subroutine chemicalFE_thermodynamic_dotState

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
   R*T*(1.0_pReal - sum(conc(1:phase_Ncomponents(phase))))* &
   log(1.0_pReal - sum(conc(1:phase_Ncomponents(phase))))
 do cpI = 1_pInt, phase_Ncomponents(phase)
   chemicalFE_thermodynamic_getEnergy = &
     chemicalFE_thermodynamic_getEnergy + &
     param(instance)%SolutionEnergy(cpI)*conc(cpI) + &
     R*T*conc(cpI)*log(conc(cpI))
   do cpJ = 1_pInt, phase_Ncomponents(phase)  
     chemicalFE_thermodynamic_getEnergy = &
       chemicalFE_thermodynamic_getEnergy + &
       param(instance)%InteractionEnergy(cpI,cpJ)*conc(cpI)*conc(cpJ)
   enddo
 enddo

end function chemicalFE_thermodynamic_getEnergy


!--------------------------------------------------------------------------------------------------
!> @brief returns the component concentration tangent for a given instance of this model
!--------------------------------------------------------------------------------------------------
subroutine chemicalFE_thermodynamic_calConcandTangent(Conc,dConcdChemPot,dConcdGradC, &
                                                      ChemPot,GradC,MechChemPot,ipc,ip,el)
 use IO, only: &
   IO_error
 use numerics, only: &
   charLength
 use material, only: &
   material_phase, &
   chemicalState, &
   chemicalConc0, &
   chemConcMapping, &
   material_homog, &
   temperature, &
   thermalMapping, &
   phase_chemicalFEInstance, &
   phase_Ncomponents

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))), intent(in)  :: &
   ChemPot, &
   GradC, &
   MechChemPot
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))), intent(out) :: &
   Conc
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el)), &
                        phase_Ncomponents(material_phase(ipc,ip,el))), intent(out) :: &
   dConcdChemPot, &
   dConcdGradC
 real(pReal) :: &
   T
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))) :: &
   Conc0, &
   ConcLastIter, &
   TempPerComponent, &
   Residual
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el)), &
                        phase_Ncomponents(material_phase(ipc,ip,el))) :: &
   Jacobian
 integer(pInt) :: &
   iter, &
   cpI, cpJ, &
   phase, &
   instance, &
   ierr, &
   ipiv(phase_Ncomponents(material_phase(ipc,ip,el)))

 external :: &
   dgesv

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 T = temperature(material_homog(ip,el))%p(thermalMapping(material_homog(ip,el))%p(ip,el))
 do cpI = 1_pInt, phase_Ncomponents(phase)
   Conc0(cpI) = chemicalConc0(cpI,phase)%p(chemConcMapping(phase)%p(ipc,ip,el))
 enddo
 Conc = Conc0
 iter = 0_pInt
 Residual = huge(1.0_pReal)
 do while (any(Residual > chemicalState(phase)%aTolState) .and. iter < param(instance)%maxNIter)
   iter = iter + 1_pInt
   ConcLastIter = Conc
   Jacobian = 0.0_pReal
   do cpI = 1_pInt, phase_Ncomponents(phase)
     Jacobian(cpI,cpI) = &
       1.0_pReal + &
       Conc(cpI)*param(instance)%GradientCoeff(cpI)/charLength/charLength/(R*T)
     TempPerComponent(cpI) = &
       ChemPot(cpI) - &
       param(instance)%SolutionEnergy(cpI) - &
       param(instance)%GradientCoeff(cpI)*(Conc(cpI) - GradC(cpI))/charLength/charLength - &
       param(instance)%MolarVolume*MechChemPot(cpI)
     do cpJ = 1_pInt, phase_Ncomponents(phase)
       Jacobian(cpI,cpJ) = &
         Jacobian(cpI,cpJ) - & 
         Conc(cpI)*Conc(cpJ)* &
         param(instance)%GradientCoeff(cpJ)/charLength/charLength/(R*T)
       TempPerComponent(cpI) = &
         tempPerComponent(cpI) - &
         param(instance)%InteractionEnergy(cpI,cpJ)*Conc0(cpJ)
     enddo
   enddo
   Residual = Conc - &
              exp(TempPerComponent(1:phase_Ncomponents(phase))/(R*T))/ &
              (1.0_pReal + sum(exp(TempPerComponent(1:phase_Ncomponents(phase))/(R*T))))

   call dgesv(phase_Ncomponents(phase),1,Jacobian,phase_Ncomponents(phase),ipiv,Residual,phase_Ncomponents(phase),ierr) 
   if (ierr /= 0_pInt) &
     call IO_error(400_pInt,el=el,ip=ip,g=ipc,ext_msg='chemicalFE thermodynamic concentration calculation did not invert')
   Conc = Conc - Residual
 enddo  
 if (iter == param(instance)%maxNIter) &
   call IO_error(400_pInt,el=el,ip=ip,g=ipc,ext_msg='chemicalFE thermodynamic concentration calculation did not converge')
 
 do cpI = 1_pInt, phase_Ncomponents(phase)
   TempPerComponent(cpI) = Conc(cpI)/(R*T + Conc(cpI)*param(instance)%GradientCoeff(cpI)/charLength/charLength)
 enddo

 dConcdChemPot   = 0.0_pReal
 do cpI = 1_pInt, phase_Ncomponents(phase)
   dConcdChemPot(cpI,cpI) = TempPerComponent(cpI) 
   do cpJ = 1_pInt, phase_Ncomponents(phase)
     dConcdChemPot(cpI,cpJ) = &
       dConcdChemPot(cpI,cpJ) - & 
       TempPerComponent(cpI)*TempPerComponent(cpJ)*(R*T)/ &
       (1.0_pReal - sum(Conc))/ &
       (1.0_pReal + sum(TempPerComponent)*(R*T)/(1.0_pReal - sum(Conc)))
   enddo
 enddo
 dConcdGradC   = 0.0_pReal
 do cpI = 1_pInt, phase_Ncomponents(phase)
   dConcdGradC(cpI,cpI) = TempPerComponent(cpI)* &
                          param(instance)%GradientCoeff(cpI)/charLength/charLength 
   do cpJ = 1_pInt, phase_Ncomponents(phase)
     dConcdGradC(cpI,cpJ) = &
       dConcdGradC(cpI,cpJ) - & 
       TempPerComponent(cpI)*TempPerComponent(cpJ)*(R*T)/ &
       (1.0_pReal - sum(Conc))/ &
       (1.0_pReal + sum(TempPerComponent)*(R*T)/(1.0_pReal - sum(Conc)))* &
       param(instance)%GradientCoeff(cpJ)/charLength/charLength
   enddo
 enddo

end subroutine chemicalFE_thermodynamic_calConcandTangent


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
   phase_chemicalFEInstance

 implicit none
 integer(pInt),                                   intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element                                                                                        !< microstructure state
 real(pReal), dimension(chemicalFE_thermodynamic_sizePostResults(phase_chemicalFEInstance(material_phase(ipc,ip,el)))) :: &
   chemicalFE_thermodynamic_postResults
 integer(pInt) :: &
   phase, &
   instance, &
   o,c,cp

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)

 chemicalFE_thermodynamic_postResults = 0.0_pReal
 c = 0_pInt
 outputsLoop: do o = 1_pInt,chemicalFE_thermodynamic_Noutput(instance)
   select case(param(instance)%outputID(o))
     case (chemicalFE_ID)
       chemicalFE_thermodynamic_postResults(c+1_pInt) = chemicalFE_thermodynamic_getEnergy(ipc,ip,el)
       c = c + 1_pInt

     case (chemicalConc_ID)
       do cp = 1_pInt, phase_Ncomponents(phase)
         chemicalFE_thermodynamic_postResults(c+cp) = &
           chemicalConc(cp,phase)%p(chemConcMapping(phase)%p(ipc,ip,el))
       enddo
       c = c + phase_Ncomponents(phase)

   end select
 enddo outputsLoop

end function chemicalFE_thermodynamic_postResults

end module chemicalFE_thermodynamic
