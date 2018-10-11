!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for quadratic chemical free energy
!--------------------------------------------------------------------------------------------------
module chemicalFE_quadenergy
 use prec, only: &
   pReal,&
   pInt, &
   p_vec

 implicit none
 private
 integer(pInt),                       dimension(:),     allocatable,         public, protected :: &
   chemicalFE_quadenergy_sizePostResults                                                       !< cumulative size of post results

 integer(pInt),                       dimension(:,:),   allocatable, target, public :: &
   chemicalFE_quadenergy_sizePostResult                                                        !< size of each post result output

 character(len=64),                   dimension(:,:),   allocatable, target, public :: &
   chemicalFE_quadenergy_output                                                                !< name of each post result output

 integer(pInt),                       dimension(:),     allocatable, target, public :: &
   chemicalFE_quadenergy_Noutput                                                               !< number of outputs per instance of this constitution

 enum, bind(c)
   enumerator :: undefined_ID, &
                 chemicalFE_ID, &
                 chemicalConc_ID
 end enum

 type, private :: tParameters                                                                         !< container type for internal constitutive parameters
   integer(kind(undefined_ID)), dimension(:), allocatable :: & 
     outputID
   real(pReal), dimension(:,:), allocatable :: &
     QuadraticCoeff, &
     QuadraticCoeffInv
   real(pReal), dimension(:  ), allocatable :: &
     Mobility, &
     EffChargeNumber, &
     EqConc, &
     InitialAvgConc, &
     InitialDeltaConc, &
     LinearCoeff, &
     GradientCoeff
   real(pReal) :: &
     MolarVolume, &
     ConstantCoeff, &
     aTol
 end type

 type(tParameters), dimension(:), allocatable, private :: param                                       !< containers of constitutive parameters (len Ninstance)

 public :: &
   chemicalFE_quadenergy_init, &
   chemicalFE_quadenergy_dotState, &
   chemicalFE_quadenergy_getEnergy, &
   chemicalFE_quadenergy_calConcandTangent, &
   chemicalFE_quadenergy_getMobility, &
   chemicalFE_quadenergy_getEffChargeNumber, &
   chemicalFE_quadenergy_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine chemicalFE_quadenergy_init(fileUnit)
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
 use math, only: &
   math_invert
 use material, only: &
   chemConcMapping, &
   chemicalConc, &
   chemicalConc0, &
   chemicalConcRate, &
   phasememberAt, &
   phase_chemicalFE, &
   phase_chemicalFEInstance, &
   phase_Noutput, &
   CHEMICALFE_QUADENERGY_label, &
   CHEMICALFE_QUADENERGY_ID, &
   material_phase, &
   phase_Ncomponents, &
   phase_maxNcomponents, &
   chemicalState
 use config, only: &
   MATERIAL_partPhase
 use numerics,only: &
   charLength, &
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
   tempEffChargeNumber, &   
   tempEqConc, &
   tempInitialAvgConc, &
   tempInitialDeltaConc, &
   tempLinCoeff, &
   tempGradCoeff
 real(pReal), dimension(:,:,:), allocatable :: &
   tempQuadCoeff, &
   tempQuadCoeffInv
 real(pReal) :: &
   randNum, &
   InitialConc
 integer(kind(undefined_ID)), dimension(:,:), allocatable :: & 
   tempOutputID
 logical :: error  

 write(6,'(/,a)')   ' <<<+-  chemicalFE_'//CHEMICALFE_QUADENERGY_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = int(count(phase_chemicalFE == CHEMICALFE_QUADENERGY_ID),pInt)
 if (maxNinstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance


 allocate(chemicalFE_quadenergy_sizePostResults(maxNinstance),                   source=0_pInt)
 allocate(chemicalFE_quadenergy_sizePostResult(maxval(phase_Noutput),maxNinstance), &
                                                                                 source=0_pInt)
 allocate(chemicalFE_quadenergy_output(maxval(phase_Noutput),maxNinstance))
          chemicalFE_quadenergy_output               = ''
 allocate(chemicalFE_quadenergy_Noutput(maxNinstance),                       source=0_pInt)
 
 allocate(param(maxNinstance))
 
 allocate(tempOutputID        (maxval(phase_Noutput)  ,maxNinstance))
 allocate(tempMobility        (phase_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(tempEffChargeNumber (phase_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(tempEqConc          (phase_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(tempInitialAvgConc  (phase_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(tempInitialDeltaConc(phase_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(tempLinCoeff        (phase_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(tempGradCoeff       (phase_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(tempQuadCoeff       (phase_maxNcomponents, &
                               phase_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(tempQuadCoeffInv    (phase_maxNcomponents, &
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
   if (phase > 0_pInt ) then; if (phase_chemicalFE(phase) == CHEMICALFE_QUADENERGY_ID) then      ! one of my phases. Do not short-circuit here (.and. between if-statements), it's not safe in Fortran
     instance = phase_chemicalFEInstance(phase)                                                     ! which instance of my chemicalFE is present phase
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case ('chemicalfe')
             chemicalFE_quadenergy_Noutput(instance) = chemicalFE_quadenergy_Noutput(instance) + 1_pInt
             chemicalFE_quadenergy_output(chemicalFE_quadenergy_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))
             tempOutputID(chemicalFE_quadenergy_Noutput(instance),instance) = chemicalFE_ID
           case ('chemicalconc')
             chemicalFE_quadenergy_Noutput(instance) = chemicalFE_quadenergy_Noutput(instance) + 1_pInt
             chemicalFE_quadenergy_output(chemicalFE_quadenergy_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))
             tempOutputID(chemicalFE_quadenergy_Noutput(instance),instance) = chemicalConc_ID

           case default

         end select

!--------------------------------------------------------------------------------------------------
! parameters depending on number of components
       case ('molarvolume')
         param(instance)%MolarVolume = IO_floatValue(line,chunkPos,2_pInt)
         
       case ('component_tolerance')
         param(instance)%aTol = IO_floatValue(line,chunkPos,2_pInt)
         
       case ('component_mobility')
         if (chunkPos(1) /= phase_Ncomponents(phase) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_QUADENERGY_label//')')
         do j = 1_pInt, phase_Ncomponents(phase)
           tempMobility(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
        case ('effective_charge_number')
         if (chunkPos(1) /= phase_Ncomponents(phase) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_QUADENERGY_label//')')
         do j = 1_pInt, phase_Ncomponents(phase)
           tempEffChargeNumber(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
       enddo  
  
       case ('component_eqconc')
         if (chunkPos(1) /= phase_Ncomponents(phase) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_QUADENERGY_label//')')
         do j = 1_pInt, phase_Ncomponents(phase)
           tempEqConc(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
       case ('component_initialavgconc')
         if (chunkPos(1) /= phase_Ncomponents(phase) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_QUADENERGY_label//')')
         do j = 1_pInt, phase_Ncomponents(phase)
           tempInitialAvgConc(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
       case ('component_initialdeltaconc')
         if (chunkPos(1) /= phase_Ncomponents(phase) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_QUADENERGY_label//')')
         do j = 1_pInt, phase_Ncomponents(phase)
           tempInitialDeltaConc(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
       case ('component_constcoeff','component_constantcoeff')
         param(instance)%ConstantCoeff = IO_floatValue(line,chunkPos,2_pInt)
         
       case ('component_lincoeff','component_linearcoeff')
         if (chunkPos(1) /= phase_Ncomponents(phase) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_QUADENERGY_label//')')
         do j = 1_pInt, phase_Ncomponents(phase)
           tempLinCoeff(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
       case ('component_gradcoeff','component_gradientcoeff')
         if (chunkPos(1) /= phase_Ncomponents(phase) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_QUADENERGY_label//')')
         do j = 1_pInt, phase_Ncomponents(phase)
           tempGradCoeff(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
       case ('component_quadcoeff','component_quadraticcoeff')
         if (chunkPos(1) /= phase_Ncomponents(phase)* &
                            (phase_Ncomponents(phase) + 1_pInt)/2_pInt + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_QUADENERGY_label//')')
         o = 0_pInt
         do j = 1_pInt, phase_Ncomponents(phase)
           do k = j, phase_Ncomponents(phase)
             o = o + 1_pInt
             tempQuadCoeff(j,k,instance) = IO_floatValue(line,chunkPos,1_pInt+o)
             tempQuadCoeff(k,j,instance) = IO_floatValue(line,chunkPos,1_pInt+o)
           enddo
         enddo

       case default

     end select
   endif; endif
 enddo parsingFile

 sanityChecks: do phase = 1_pInt, size(phase_chemicalFE)
   myPhase: if (phase_chemicalFE(phase) == CHEMICALFE_quadenergy_ID) then
     instance = phase_chemicalFEInstance(phase)
     allocate(param(instance)%OutputID         (chemicalFE_quadenergy_Noutput(instance)))
     param(instance)%OutputID = tempOutputID(1:chemicalFE_quadenergy_Noutput(instance),instance)
     allocate(param(instance)%Mobility         (phase_Ncomponents(phase)))
     allocate(param(instance)%EffChargeNumber  (phase_Ncomponents(phase)))
     allocate(param(instance)%EqConc           (phase_Ncomponents(phase)))
     allocate(param(instance)%InitialAvgConc   (phase_Ncomponents(phase)))
     allocate(param(instance)%InitialDeltaConc (phase_Ncomponents(phase)))
     allocate(param(instance)%LinearCoeff      (phase_Ncomponents(phase)))
     allocate(param(instance)%GradientCoeff    (phase_Ncomponents(phase)))
     allocate(param(instance)%QuadraticCoeff   (phase_Ncomponents(phase), &
                                                phase_Ncomponents(phase)))
     allocate(param(instance)%QuadraticCoeffInv(phase_Ncomponents(phase), &
                                                phase_Ncomponents(phase)))
     param(instance)%Mobility        = tempMobility        (1:phase_Ncomponents(phase),instance)
     param(instance)%EffChargeNumber = tempEffChargeNumber (1:phase_Ncomponents(phase),instance)
     param(instance)%EqConc          = tempEqConc          (1:phase_Ncomponents(phase),instance)
     param(instance)%InitialAvgConc  = tempInitialAvgConc  (1:phase_Ncomponents(phase),instance)
     param(instance)%InitialDeltaConc= tempInitialDeltaConc(1:phase_Ncomponents(phase),instance)
     param(instance)%LinearCoeff     = tempLinCoeff        (1:phase_Ncomponents(phase),instance)
     param(instance)%GradientCoeff   = tempGradCoeff       (1:phase_Ncomponents(phase),instance)
     param(instance)%QuadraticCoeff  = tempQuadCoeff       (1:phase_Ncomponents(phase), &
                                                            1:phase_Ncomponents(phase),instance)
     do j = 1_pInt, phase_Ncomponents(phase)
       tempQuadCoeff(j,j,instance) = tempQuadCoeff(j,j,instance) + &
                                     param(instance)%GradientCoeff(j)/charLength/charLength 
     enddo
     call math_invert(phase_Ncomponents(phase), &
                      tempQuadCoeff(1:phase_Ncomponents(phase), &
                                    1:phase_Ncomponents(phase),instance), &
                      param(instance)%QuadraticCoeffInv, &
                      error)
     if (error) &
       call IO_error(400_pInt,el=instance,ext_msg='quad energy ('//CHEMICALFE_QUADENERGY_label//')')
     if (any(param(instance)%Mobility < 0.0_pReal)) &
       call IO_error(211_pInt,el=instance,ext_msg='mobility ('//CHEMICALFE_QUADENERGY_label//')')
     if (any(param(instance)%EqConc   < 0.0_pReal)) &
       call IO_error(211_pInt,el=instance,ext_msg='eqconc ('//CHEMICALFE_QUADENERGY_label//')')
   endif myPhase
 enddo sanityChecks

 initializeInstances: do phase = 1_pInt, size(phase_chemicalFE)                                     ! loop through all phases in material.config
   myPhase2: if (phase_chemicalFE(phase) == CHEMICALFE_quadenergy_ID) then                       ! only consider my phase
     NipcMyPhase = count(material_phase == phase)                                                   ! number of IPCs containing my phase
     instance = phase_chemicalFEInstance(phase)                                                     ! which instance of my phase

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,chemicalFE_quadenergy_Noutput(instance)
       select case(param(instance)%outputID(o))
         case(chemicalFE_ID)
           mySize = 1_pInt
         case(chemicalConc_ID)
           mySize = phase_Ncomponents(phase)
         case default
       end select

       outputFound: if (mySize > 0_pInt) then
         chemicalFE_quadenergy_sizePostResult(o,instance) = mySize
         chemicalFE_quadenergy_sizePostResults(instance)  = chemicalFE_quadenergy_sizePostResults(instance) + mySize
       endif outputFound
     enddo outputsLoop
!--------------------------------------------------------------------------------------------------
! allocate state arrays
     sizeState = phase_Ncomponents(phase)
     sizeDotState = sizeState
     sizeDeltaState = 0_pInt
     chemicalState(phase)%sizeState = sizeState
     chemicalState(phase)%sizeDotState = sizeDotState
     chemicalState(phase)%sizeDeltaState = sizeDeltaState
     chemicalState(phase)%sizePostResults = chemicalFE_quadenergy_sizePostResults(instance)
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
     do o = 1_pInt, NipcMyPhase
       do j = 1_pInt, phase_Ncomponents(phase)
         call random_number(randNum)
         InitialConc = param(instance)%InitialAvgConc(j) + param(instance)%InitialDeltaConc(j)*(randNum - 0.5_pReal)
         chemicalState(phase)%state0   (j,o) = InitialConc
         chemicalState(phase)%subState0(j,o) = InitialConc
         chemicalState(phase)%state    (j,o) = InitialConc
       enddo  
     enddo
     chemicalConc (phase)%p => chemicalState(phase)%state 
     chemicalConc0(phase)%p => chemicalState(phase)%state0
     allocate(chemicalConcRate(phase)%p(phase_Ncomponents(phase),NipcMyPhase), source=0.0_pReal)
   endif myPhase2
 enddo initializeInstances   

end subroutine chemicalFE_quadenergy_init


!--------------------------------------------------------------------------------------------------
!> @brief returns the chemical concentration rates
!--------------------------------------------------------------------------------------------------
subroutine chemicalFE_quadenergy_dotState(ipc,ip,el)
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
     chemicalConcRate(phase)%p(cp,chemConcMapping(phase)%p(ipc,ip,el))
 enddo

end subroutine chemicalFE_quadenergy_dotState

!--------------------------------------------------------------------------------------------------
!> @brief returns the chemical energy for a given instance of this model
!--------------------------------------------------------------------------------------------------
function chemicalFE_quadenergy_getEnergy(ipc,ip,el)
 use material, only: &
   chemicalConc, &
   chemConcMapping, &
   material_phase, &
   phase_Ncomponents, &
   phase_chemicalFEInstance
   
 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))) :: &
   conc, &
   chemPot
 real(pReal) :: &
   chemicalFE_quadenergy_getEnergy
 integer(pInt) :: &
   cpI, cpJ, &
   phase, &
   instance    

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 conc = 0.0_pReal
 do cpI = 1_pInt, phase_Ncomponents(phase)
   conc(cpI) = chemicalConc(phase)%p(cpI,chemConcMapping(phase)%p(ipc, ip, el))
 enddo
 chemPot = 0.0_pReal
 do cpI = 1_pInt, phase_Ncomponents(phase)
   chemPot(cpI) = chemPot(cpI) + &
                  param(instance)%LinearCoeff(cpI)
    do cpJ = 1_pInt, phase_Ncomponents(phase)
      chemPot(cpI) = chempot(cpI) + &
                     param(instance)%QuadraticCoeff(cpI,cpJ)* &
                     (conc(cpJ)-param(instance)%EqConc(cpJ))
    enddo                 
 enddo
 chemicalFE_quadenergy_getEnergy = &
   param(instance)%ConstantCoeff
 do cpI = 1_pInt, phase_Ncomponents(phase)
   chemicalFE_quadenergy_getEnergy = &
     chemicalFE_quadenergy_getEnergy + &
     param(instance)%LinearCoeff(cpI)* &
     (conc(cpI) - param(instance)%EqConc(cpI)) - &
     chemPot(cpI)*conc(cpI)
   do cpJ = 1_pInt, phase_Ncomponents(phase)  
     chemicalFE_quadenergy_getEnergy = &
       chemicalFE_quadenergy_getEnergy + &
       0.5_pReal* &
       param(instance)%QuadraticCoeff(cpI,cpJ)* &
       (conc(cpI) - param(instance)%EqConc(cpI))* &
       (conc(cpJ) - param(instance)%EqConc(cpJ))
   enddo
 enddo
 chemicalFE_quadenergy_getEnergy = &
   chemicalFE_quadenergy_getEnergy/param(instance)%MolarVolume

end function chemicalFE_quadenergy_getEnergy


!--------------------------------------------------------------------------------------------------
!> @brief returns the component concentration tangent for a given instance of this model
!--------------------------------------------------------------------------------------------------
subroutine chemicalFE_quadenergy_calConcandTangent(Conc,dConcdChemPot,dConcdGradC, & 
                                                   ChemPot,GradC,MechChemPot,IntfChemPot,ElectroChemPot,ipc,ip,el)
 use numerics, only: &
   charLength
 use material, only: &
   material_phase, &
   phase_Ncomponents, &
   phase_chemicalFEInstance

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))), intent(in)  :: &
   MechChemPot, &
   IntfChemPot, &
   ElectroChemPot, &
   ChemPot, &
   GradC
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))), intent(out) :: &
   Conc
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el)), &
                        phase_Ncomponents(material_phase(ipc,ip,el))), intent(out) :: &
   dConcdChemPot, &
   dConcdGradC
 integer(pInt) :: &
   cpI, cpJ, &
   phase, &
   instance    

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 Conc = 0.0_pReal
 dConcdChemPot = 0.0_pReal
 dConcdGradC = 0.0_pReal
 do cpI = 1_pInt, phase_Ncomponents(phase)
   Conc(cpI) = param(instance)%EqConc(cpI)
   do cpJ = 1_pInt, phase_Ncomponents(phase)
     Conc(cpI) = Conc(cpI) + &
                 param(instance)%QuadraticCoeffInv(cpI,cpJ)* &
                 (ChemPot(cpJ) - MechChemPot(cpJ) - IntfChemPot(cpJ) - &
                  ElectroChemPot(cpJ) - param(instance)%LinearCoeff(cpJ) + &
                  GradC(cpJ)*param(instance)%GradientCoeff(cpJ)/charLength/charLength)
     dConcdChemPot(cpI,cpJ) = &
       param(instance)%QuadraticCoeffInv(cpI,cpJ)
     dConcdGradC  (cpI,cpJ) = &
       param(instance)%QuadraticCoeffInv(cpI,cpJ)* &
       param(instance)%GradientCoeff(cpJ)/charLength/charLength
   enddo
 enddo

end subroutine chemicalFE_quadenergy_calConcandTangent


!--------------------------------------------------------------------------------------------------
!> @brief returns the component mobility for a given instance of this model
!--------------------------------------------------------------------------------------------------
function chemicalFE_quadenergy_getMobility(ipc,ip,el)
 use material, only: &
   material_phase, &
   phase_Ncomponents, &
   phase_chemicalFEInstance

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))) :: &
   chemicalFE_quadenergy_getMobility
 integer(pInt) :: &
   cp, &
   phase, &
   instance    

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 chemicalFE_quadenergy_getMobility = 0.0_pReal
 do cp = 1_pInt, phase_Ncomponents(phase)
   chemicalFE_quadenergy_getMobility(cp) = &
     param(instance)%Mobility(cp)
 enddo

end function chemicalFE_quadenergy_getMobility


!--------------------------------------------------------------------------------------------------
!> @brief returns the component charge number for a given instance of this model
!--------------------------------------------------------------------------------------------------
function chemicalFE_quadenergy_getEffChargeNumber(ipc,ip,el)
 use material, only: &
   material_phase, &
   phase_Ncomponents, &
   phase_chemicalFEInstance

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))) :: &
   chemicalFE_quadenergy_getEffChargeNumber
 integer(pInt) :: &
   cp, &
   phase, &
   instance    

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 chemicalFE_quadenergy_getEffChargeNumber = 0.0_pReal
 do cp = 1_pInt, phase_Ncomponents(phase)
   chemicalFE_quadenergy_getEffChargeNumber(cp) = &
     param(instance)%EffChargeNumber(cp)
 enddo

end function chemicalFE_quadenergy_getEffChargeNumber




!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function chemicalFE_quadenergy_postResults(ipc,ip,el)
 use material, only: &
   material_phase, &
   chemicalConc, &
   chemConcMapping, &
   phase_Ncomponents, &
   phase_chemicalFEInstance

 implicit none
 integer(pInt),                                   intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element                                                                                        !< microstructure state

 real(pReal), dimension(chemicalFE_quadenergy_sizePostResults(phase_chemicalFEInstance(material_phase(ipc,ip,el)))) :: &
   chemicalFE_quadenergy_postResults
 integer(pInt) :: &
   phase, &
   instance, &
   o,c,cp

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 
 chemicalFE_quadenergy_postResults = 0.0_pReal
 c = 0_pInt
 outputsLoop: do o = 1_pInt,chemicalFE_quadenergy_Noutput(instance)
   select case(param(instance)%outputID(o))
     case (chemicalFE_ID)
       chemicalFE_quadenergy_postResults(c+1_pInt) = chemicalFE_quadenergy_getEnergy(ipc,ip,el)
       c = c + 1_pInt

     case (chemicalConc_ID)
       do cp = 1_pInt, phase_Ncomponents(phase)
         chemicalFE_quadenergy_postResults(c+cp) = &
           chemicalConc(phase)%p(cp,chemConcMapping(phase)%p(ipc,ip,el))
       enddo
       c = c + phase_Ncomponents(phase)

   end select
 enddo outputsLoop

end function chemicalFE_quadenergy_postResults

end module chemicalFE_quadenergy
