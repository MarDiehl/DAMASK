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
                 chemicalPot_ID, &
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
     EqConc, &
     InitialConc, &
     LinearCoeff
   real(pReal) :: &
     ConstantCoeff
 end type

 type(tParameters), dimension(:), allocatable, private :: param                                       !< containers of constitutive parameters (len Ninstance)

 public :: &
   chemicalFE_quadenergy_init, &
   chemicalFE_quadenergy_getEnergy, &
   chemicalFE_quadenergy_getChemPot, &
   chemicalFE_quadenergy_getConc, &
   chemicalFE_quadenergy_putConc, &
   chemicalFE_quadenergy_getConcTangent, &
   chemicalFE_quadenergy_getMobility, &
   chemicalFE_quadenergy_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine chemicalFE_quadenergy_init(fileUnit)
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
 use math, only: &
   math_invert
 use material, only: &
   chemConcMapping, &
   chemicalConc, &
   phasememberAt, &
   phase_chemicalFE, &
   phase_chemicalFEInstance, &
   phase_Noutput, &
   CHEMICALFE_QUADENERGY_label, &
   CHEMICALFE_QUADENERGY_ID, &
   material_phase, &
   phase_Ncomponents, &
   phase_maxNcomponents, &
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
 real(pReal), dimension(:,:), allocatable :: &
   tempMobility, &
   tempEqConc, &
   tempInitialConc, &
   tempLinCoeff
 real(pReal), dimension(:,:,:), allocatable :: &
   tempQuadCoeff, &
   tempQuadCoeffInv
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
 
 allocate(tempOutputID   (maxval(phase_Noutput)  ,maxNinstance))
 allocate(tempMobility   (phase_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(tempEqConc     (phase_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(tempInitialConc(phase_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(tempLinCoeff   (phase_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(tempQuadCoeff  (phase_maxNcomponents, &
                          phase_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(tempQuadCoeffInv(phase_maxNcomponents, &
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
           case ('chemicalpot')
             chemicalFE_quadenergy_Noutput(instance) = chemicalFE_quadenergy_Noutput(instance) + 1_pInt
             chemicalFE_quadenergy_output(chemicalFE_quadenergy_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))
             tempOutputID(chemicalFE_quadenergy_Noutput(instance),instance) = chemicalPot_ID
           case ('chemicalconc')
             chemicalFE_quadenergy_Noutput(instance) = chemicalFE_quadenergy_Noutput(instance) + 1_pInt
             chemicalFE_quadenergy_output(chemicalFE_quadenergy_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))
             tempOutputID(chemicalFE_quadenergy_Noutput(instance),instance) = chemicalConc_ID

           case default

         end select

!--------------------------------------------------------------------------------------------------
! parameters depending on number of components
       case ('component_mobility')
         if (chunkPos(1) /= phase_Ncomponents(phase) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_QUADENERGY_label//')')
         do j = 1_pInt, phase_Ncomponents(phase)
           tempMobility(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
       case ('component_eqconc')
         if (chunkPos(1) /= phase_Ncomponents(phase) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_QUADENERGY_label//')')
         do j = 1_pInt, phase_Ncomponents(phase)
           tempEqConc(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
       case ('component_initialconc')
         if (chunkPos(1) /= phase_Ncomponents(phase) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_QUADENERGY_label//')')
         do j = 1_pInt, phase_Ncomponents(phase)
           tempInitialConc(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
       case ('component_constcoeff')
         param(instance)%ConstantCoeff = IO_floatValue(line,chunkPos,2_pInt)
         
       case ('component_lincoeff')
         if (chunkPos(1) /= phase_Ncomponents(phase) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_QUADENERGY_label//')')
         do j = 1_pInt, phase_Ncomponents(phase)
           tempLinCoeff(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
       case ('component_quadcoeff')
         if (chunkPos(1) /= phase_Ncomponents(phase)* &
                            (phase_Ncomponents(phase) + 1_pInt) + 1_pInt) &
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
     allocate(param(instance)%EqConc           (phase_Ncomponents(phase)))
     allocate(param(instance)%InitialConc      (phase_Ncomponents(phase)))
     allocate(param(instance)%LinearCoeff      (phase_Ncomponents(phase)))
     allocate(param(instance)%QuadraticCoeff   (phase_Ncomponents(phase), &
                                                phase_Ncomponents(phase)))
     allocate(param(instance)%QuadraticCoeffInv(phase_Ncomponents(phase), &
                                                phase_Ncomponents(phase)))
     param(instance)%Mobility       = tempMobility   (1:phase_Ncomponents(phase),instance)
     param(instance)%EqConc         = tempEqConc     (1:phase_Ncomponents(phase),instance)
     param(instance)%InitialConc    = tempInitialConc(1:phase_Ncomponents(phase),instance)
     param(instance)%LinearCoeff    = tempLinCoeff   (1:phase_Ncomponents(phase),instance)
     param(instance)%QuadraticCoeff = tempQuadCoeff  (1:phase_Ncomponents(phase), &
                                                      1:phase_Ncomponents(phase),instance)
     call math_invert(phase_Ncomponents(phase), &
                      param(instance)%QuadraticCoeff, &
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
         case(chemicalPot_ID, &
              chemicalConc_ID &
              )
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
     sizeState = 0_pInt
     sizeDotState = sizeState
     sizeDeltaState = 0_pInt
     chemicalState(phase)%sizeState = sizeState
     chemicalState(phase)%sizeDotState = sizeDotState
     chemicalState(phase)%sizeDeltaState = sizeDeltaState
     chemicalState(phase)%sizePostResults = chemicalFE_quadenergy_sizePostResults(instance)
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
       allocate(chemicalConc(j,phase)%p(NipcMyPhase), source=param(instance)%InitialConc(j))
     enddo
   endif myPhase2
 enddo initializeInstances   

end subroutine chemicalFE_quadenergy_init


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
   conc
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
   conc(cpI) = chemicalConc(cpI,phase)%p(chemConcMapping(phase)%p(ipc, ip, el))
 enddo
 chemicalFE_quadenergy_getEnergy = &
   param(instance)%ConstantCoeff
 do cpI = 1_pInt, phase_Ncomponents(phase)
   chemicalFE_quadenergy_getEnergy = &
     chemicalFE_quadenergy_getEnergy + &
     param(instance)%LinearCoeff(cpI)* &
     (conc(cpI) - param(instance)%EqConc(cpI))
   do cpJ = 1_pInt, phase_Ncomponents(phase)  
     chemicalFE_quadenergy_getEnergy = &
       chemicalFE_quadenergy_getEnergy + &
       0.5_pReal* &
       param(instance)%QuadraticCoeff(cpI,cpJ)* &
       (conc(cpI) - param(instance)%EqConc(cpI))* &
       (conc(cpJ) - param(instance)%EqConc(cpJ))
   enddo
 enddo

end function chemicalFE_quadenergy_getEnergy

!--------------------------------------------------------------------------------------------------
!> @brief returns the component chemical potential for a given instance of this model
!--------------------------------------------------------------------------------------------------
function chemicalFE_quadenergy_getChemPot(ipc,ip,el)
 use material, only: &
   material_phase, &
   chemicalConc, &
   chemConcMapping, &
   phase_Ncomponents, &
   phase_chemicalFEInstance

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))) :: &
   conc
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))) :: &
   chemicalFE_quadenergy_getChemPot
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
 chemicalFE_quadenergy_getChemPot = 0.0_pReal
 do cpI = 1_pInt, phase_Ncomponents(phase)
   chemicalFE_quadenergy_getChemPot(cpI) = &
     param(instance)%LinearCoeff(cpI)
   do cpJ = 1_pInt, phase_Ncomponents(phase)
     chemicalFE_quadenergy_getChemPot(cpI) = &
       chemicalFE_quadenergy_getChemPot(cpI) + &
       param(instance)%QuadraticCoeff(cpI,cpJ)* &
       (conc(cpJ) - param(instance)%EqConc(cpJ))
   enddo
 enddo

end function chemicalFE_quadenergy_getChemPot


!--------------------------------------------------------------------------------------------------
!> @brief returns the component concentration for a given instance of this model
!--------------------------------------------------------------------------------------------------
function chemicalFE_quadenergy_getConc(chempot,ipc,ip,el)
 use material, only: &
   material_phase, &
   phase_Ncomponents, &
   phase_chemicalFEInstance

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))), intent(in) :: &
   chempot
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))) :: &
   chemicalFE_quadenergy_getConc
 integer(pInt) :: &
   cpI, cpJ, &
   phase, &
   instance    

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 chemicalFE_quadenergy_getConc = 0.0_pReal
 do cpI = 1_pInt, phase_Ncomponents(phase)
   chemicalFE_quadenergy_getConc(cpI) = &
     param(instance)%EqConc(cpI)
   do cpJ = 1_pInt, phase_Ncomponents(phase)
     chemicalFE_quadenergy_getConc(cpI) = &
       chemicalFE_quadenergy_getConc(cpI) + &
       param(instance)%QuadraticCoeffInv(cpI,cpJ)* &
       (chempot(cpJ) - param(instance)%LinearCoeff(cpJ))
   enddo
 enddo

end function chemicalFE_quadenergy_getConc


!--------------------------------------------------------------------------------------------------
!> @brief stores the component concentration for a given instance of this model
!--------------------------------------------------------------------------------------------------
subroutine chemicalFE_quadenergy_putConc(chempot,ipc,ip,el)
 use material, only: &
   material_phase, &
   phase_Ncomponents, &
   phase_chemicalFEInstance, &
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
 conc = chemicalFE_quadenergy_getConc(chempot,ipc,ip,el)
 do cp = 1_pInt, phase_Ncomponents(phase)
   chemicalConc(cp,phase)%p(chemConcMapping(phase)%p(ipc,ip,el)) = &
     conc(cp)
 enddo

end subroutine chemicalFE_quadenergy_putConc


!--------------------------------------------------------------------------------------------------
!> @brief returns the component concentration tangent for a given instance of this model
!--------------------------------------------------------------------------------------------------
subroutine chemicalFE_quadenergy_getConcTangent(dConcdChemPot,dConcdGradC,ipc,ip,el)
 use material, only: &
   material_phase, &
   phase_Ncomponents, &
   phase_chemicalFEInstance

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
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
 dConcdChemPot = 0.0_pReal
 dConcdGradC = 0.0_pReal
 do cpI = 1_pInt, phase_Ncomponents(phase)
   do cpJ = 1_pInt, phase_Ncomponents(phase)
     dConcdChemPot(cpI,cpJ) = &
       param(instance)%QuadraticCoeffInv(cpI,cpJ)
   enddo
 enddo

end subroutine chemicalFE_quadenergy_getConcTangent


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
 real(pReal), dimension(phase_Ncomponents(material_phase(ipc,ip,el))) :: &
   conc, tempPerComponent

 integer(pInt) :: &
   phase, &
   instance, &
   o,c,cp

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 do cp = 1_pInt, phase_Ncomponents(phase)
   conc(cp) = chemicalConc(cp,phase)%p(chemConcMapping(phase)%p(ipc,ip,el))
 enddo
 
 chemicalFE_quadenergy_postResults = 0.0_pReal
 c = 0_pInt
 outputsLoop: do o = 1_pInt,chemicalFE_quadenergy_Noutput(instance)
   select case(param(instance)%outputID(o))
     case (chemicalFE_ID)
       chemicalFE_quadenergy_postResults(c+1_pInt) = chemicalFE_quadenergy_getEnergy(ipc,ip,el)
       c = c + 1_pInt

     case (chemicalPot_ID)
       tempPerComponent = chemicalFE_quadenergy_getChemPot(ipc,ip,el)
       chemicalFE_quadenergy_postResults(c+1_pInt:c+phase_Ncomponents(phase)) = &
         tempPerComponent(1_pInt:phase_Ncomponents(phase))
       c = c + phase_Ncomponents(phase)

     case (chemicalConc_ID)
       chemicalFE_quadenergy_postResults(c+1_pInt:c+phase_Ncomponents(phase)) = &
         conc(1_pInt:phase_Ncomponents(phase))
       c = c + phase_Ncomponents(phase)

   end select
 enddo outputsLoop

end function chemicalFE_quadenergy_postResults

end module chemicalFE_quadenergy
