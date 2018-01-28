!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for quadratic chemical free energy
!--------------------------------------------------------------------------------------------------
module chemicalFE_quadenergy
 use prec, only: &
   pReal,&
   pInt

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

 integer(pInt),                       dimension(:),     allocatable,         private, protected :: &
   chemicalFE_quadenergy_Ncomponents                                                           !< no. of chemical components

 real(pReal),                         dimension(:),     allocatable,         private :: &
   chemicalFE_quadenergy_ConstCoeff                                                          
 real(pReal),                         dimension(:,:),   allocatable,         private :: &
   chemicalFE_quadenergy_Mobility, &                                                         
   chemicalFE_quadenergy_EqConc, &                                                          
   chemicalFE_quadenergy_LinCoeff                                                        
 real(pReal),                         dimension(:,:,:), allocatable,         private :: &
   chemicalFE_quadenergy_QuadCoeff, &
   chemicalFE_quadenergy_QuadCoeffInv

 enum, bind(c)
   enumerator :: undefined_ID, &
                 chemicalFE_ID, &
                 chemicalPot_ID, &
                 chemicalConc_ID
 end enum
 integer(kind(undefined_ID)),         dimension(:,:),   allocatable,          private :: &
   chemicalFE_quadenergy_outputID                                                              !< ID of each post result output

 public :: &
   chemicalFE_quadenergy_init, &
   chemicalFE_quadenergy_getEnergy, &
   chemicalFE_quadenergy_getChemPot, &
   chemicalFE_quadenergy_getConc, &
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
   phase_chemicalFE, &
   phase_chemicalFEInstance, &
   phase_Noutput, &
   material_maxNcomponents, &
   CHEMICALFE_QUADENERGY_label, &
   CHEMICALFE_QUADENERGY_ID, &
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
 allocate(chemicalFE_quadenergy_outputID(maxval(phase_Noutput),maxNinstance),source=undefined_ID)
 allocate(chemicalFE_quadenergy_Noutput(maxNinstance),                       source=0_pInt)
 allocate(chemicalFE_quadenergy_Ncomponents                       (maxNinstance),source=  0_pInt )
 allocate(chemicalFE_quadenergy_ConstCoeff                        (maxNinstance),source=0.0_pReal)
 allocate(chemicalFE_quadenergy_Mobility  (material_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(chemicalFE_quadenergy_EqConc    (material_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(chemicalFE_quadenergy_LinCoeff  (material_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(chemicalFE_quadenergy_QuadCoeff (material_maxNcomponents, &
                                           material_maxNcomponents,maxNinstance),source=0.0_pReal)
 allocate(chemicalFE_quadenergy_QuadCoeffInv(material_maxNcomponents, &
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
   if (phase > 0_pInt ) then; if (phase_chemicalFE(phase) == CHEMICALFE_QUADENERGY_ID) then      ! one of my phases. Do not short-circuit here (.and. between if-statements), it's not safe in Fortran
     instance = phase_chemicalFEInstance(phase)                                                     ! which instance of my chemicalFE is present phase
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case ('chemicalfe')
             chemicalFE_quadenergy_Noutput(instance) = chemicalFE_quadenergy_Noutput(instance) + 1_pInt
             chemicalFE_quadenergy_outputID(chemicalFE_quadenergy_Noutput(instance),instance) = chemicalFE_ID
             chemicalFE_quadenergy_output(chemicalFE_quadenergy_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('chemicalpot')
             chemicalFE_quadenergy_Noutput(instance) = chemicalFE_quadenergy_Noutput(instance) + 1_pInt
             chemicalFE_quadenergy_outputID(chemicalFE_quadenergy_Noutput(instance),instance) = chemicalPot_ID
             chemicalFE_quadenergy_output(chemicalFE_quadenergy_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('chemicalconc')
             chemicalFE_quadenergy_Noutput(instance) = chemicalFE_quadenergy_Noutput(instance) + 1_pInt
             chemicalFE_quadenergy_outputID(chemicalFE_quadenergy_Noutput(instance),instance) = chemicalConc_ID
             chemicalFE_quadenergy_output(chemicalFE_quadenergy_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))

           case default

         end select

       case ('ncomponents')
         chemicalFE_quadenergy_Ncomponents(instance) = IO_intValue(line,chunkPos,2_pInt)

!--------------------------------------------------------------------------------------------------
! parameters depending on number of components
       case ('component_mobility')
         if (chunkPos(1) /= chemicalFE_quadenergy_Ncomponents(instance) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_QUADENERGY_label//')')
         do j = 1_pInt, chemicalFE_quadenergy_Ncomponents(instance)
           chemicalFE_quadenergy_Mobility(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
       case ('component_eqconc')
         if (chunkPos(1) /= chemicalFE_quadenergy_Ncomponents(instance) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_QUADENERGY_label//')')
         do j = 1_pInt, chemicalFE_quadenergy_Ncomponents(instance)
           chemicalFE_quadenergy_EqConc(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
       case ('component_constcoeff')
         chemicalFE_quadenergy_ConstCoeff(instance) = IO_floatValue(line,chunkPos,2_pInt)
         
       case ('component_lincoeff')
         if (chunkPos(1) /= chemicalFE_quadenergy_Ncomponents(instance) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_QUADENERGY_label//')')
         do j = 1_pInt, chemicalFE_quadenergy_Ncomponents(instance)
           chemicalFE_quadenergy_LinCoeff(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
       case ('component_quadcoeff')
         if (chunkPos(1) /= chemicalFE_quadenergy_Ncomponents(instance)* &
                            (chemicalFE_quadenergy_Ncomponents(instance) + 1_pInt) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//CHEMICALFE_QUADENERGY_label//')')
         o = 0_pInt
         do j = 1_pInt, chemicalFE_quadenergy_Ncomponents(instance)
           do k = j, chemicalFE_quadenergy_Ncomponents(instance)
             o = o + 1_pInt
             chemicalFE_quadenergy_QuadCoeff(j,k,instance) = IO_floatValue(line,chunkPos,1_pInt+o)
             chemicalFE_quadenergy_QuadCoeff(k,j,instance) = IO_floatValue(line,chunkPos,1_pInt+o)
           enddo
         enddo

       case default

     end select
   endif; endif
 enddo parsingFile

 sanityChecks: do phase = 1_pInt, size(phase_chemicalFE)
   myPhase: if (phase_chemicalFE(phase) == CHEMICALFE_quadenergy_ID) then
     instance = phase_chemicalFEInstance(phase)
     call math_invert(chemicalFE_quadenergy_Ncomponents(instance), &
                      chemicalFE_quadenergy_QuadCoeff(1:chemicalFE_quadenergy_Ncomponents(instance), &
                                                      1:chemicalFE_quadenergy_Ncomponents(instance), &
                                                      instance), &
                      chemicalFE_quadenergy_QuadCoeffInv(1:chemicalFE_quadenergy_Ncomponents(instance), &
                                                         1:chemicalFE_quadenergy_Ncomponents(instance), &
                                                         instance), &
                      error)
     if (error) &
       call IO_error(400_pInt,el=instance,ext_msg='quad energy ('//CHEMICALFE_QUADENERGY_label//')')
     if (any(chemicalFE_quadenergy_Mobility(:,instance) < 0.0_pReal)) &
       call IO_error(211_pInt,el=instance,ext_msg='mobility ('//CHEMICALFE_QUADENERGY_label//')')
     if (any(chemicalFE_quadenergy_EqConc(:,instance) < 0.0_pReal)) &
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
       select case(chemicalFE_quadenergy_outputID(o,instance))
         case(chemicalFE_ID)
           mySize = 1_pInt
         case(chemicalPot_ID, &
              chemicalConc_ID &
              )
           mySize = chemicalFE_quadenergy_Ncomponents(instance)
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
   endif myPhase2
 enddo initializeInstances

end subroutine chemicalFE_quadenergy_init


!--------------------------------------------------------------------------------------------------
!> @brief returns the chemical energy for a given instance of this model
!--------------------------------------------------------------------------------------------------
function chemicalFE_quadenergy_getEnergy(conc,ipc,ip,el)
 use material, only: &
   material_phase, &
   phase_chemicalFEInstance, &
   material_maxNcomponents

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(material_maxNcomponents), intent(in) :: &
   conc
 real(pReal) :: &
   chemicalFE_quadenergy_getEnergy
 integer(pInt) :: &
   cpI, cpJ, &
   phase, &
   instance    

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 chemicalFE_quadenergy_getEnergy = &
   chemicalFE_quadenergy_ConstCoeff(instance)
 do cpI = 1_pInt, chemicalFE_quadenergy_Ncomponents(instance)
   chemicalFE_quadenergy_getEnergy = &
     chemicalFE_quadenergy_getEnergy + &
     chemicalFE_quadenergy_LinCoeff(cpI,instance)* &
     (conc(cpI) - chemicalFE_quadenergy_EqConc(cpI,instance))
   do cpJ = 1_pInt, chemicalFE_quadenergy_Ncomponents(instance)  
     chemicalFE_quadenergy_getEnergy = &
       chemicalFE_quadenergy_getEnergy + &
       0.5_pReal* &
       chemicalFE_quadenergy_QuadCoeff(cpI,cpJ,instance)* &
       (conc(cpI) - chemicalFE_quadenergy_EqConc(cpI,instance))* &
       (conc(cpJ) - chemicalFE_quadenergy_EqConc(cpJ,instance))
   enddo
 enddo

end function chemicalFE_quadenergy_getEnergy

!--------------------------------------------------------------------------------------------------
!> @brief returns the component chemical potential for a given instance of this model
!--------------------------------------------------------------------------------------------------
function chemicalFE_quadenergy_getChemPot(conc,ipc,ip,el)
 use material, only: &
   material_phase, &
   phase_chemicalFEInstance, &
   material_maxNcomponents

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(material_maxNcomponents), intent(in) :: &
   conc
 real(pReal), dimension(material_maxNcomponents) :: &
   chemicalFE_quadenergy_getChemPot
 integer(pInt) :: &
   cpI, cpJ, &
   phase, &
   instance    

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 chemicalFE_quadenergy_getChemPot = 0.0_pReal
 do cpI = 1_pInt, chemicalFE_quadenergy_Ncomponents(instance)
   chemicalFE_quadenergy_getChemPot(cpI) = &
     chemicalFE_quadenergy_LinCoeff(cpI,instance)
   do cpJ = 1_pInt, chemicalFE_quadenergy_Ncomponents(instance)
     chemicalFE_quadenergy_getChemPot(cpI) = &
       chemicalFE_quadenergy_getChemPot(cpI) + &
       chemicalFE_quadenergy_QuadCoeff(cpI,cpJ,instance)* &
       (conc(cpJ) - chemicalFE_quadenergy_EqConc(cpJ,instance))
   enddo
 enddo

end function chemicalFE_quadenergy_getChemPot


!--------------------------------------------------------------------------------------------------
!> @brief returns the component concentration for a given instance of this model
!--------------------------------------------------------------------------------------------------
function chemicalFE_quadenergy_getConc(chempot,ipc,ip,el)
 use material, only: &
   material_phase, &
   phase_chemicalFEInstance, &
   material_maxNcomponents

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(material_maxNcomponents), intent(in) :: &
   chempot
 real(pReal), dimension(material_maxNcomponents) :: &
   chemicalFE_quadenergy_getConc
 integer(pInt) :: &
   cpI, cpJ, &
   phase, &
   instance    

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 chemicalFE_quadenergy_getConc = 0.0_pReal
 do cpI = 1_pInt, chemicalFE_quadenergy_Ncomponents(instance)
   chemicalFE_quadenergy_getConc(cpI) = &
     chemicalFE_quadenergy_EqConc(cpI,instance)
   do cpJ = 1_pInt, chemicalFE_quadenergy_Ncomponents(instance)
     chemicalFE_quadenergy_getConc(cpI) = &
       chemicalFE_quadenergy_getConc(cpI) + &
       chemicalFE_quadenergy_QuadCoeffInv(cpI,cpJ,instance)* &
       (chempot(cpJ) - chemicalFE_quadenergy_LinCoeff(cpJ,instance))
   enddo
 enddo

end function chemicalFE_quadenergy_getConc


!--------------------------------------------------------------------------------------------------
!> @brief returns the component concentration tangent for a given instance of this model
!--------------------------------------------------------------------------------------------------
function chemicalFE_quadenergy_getConcTangent(ipc,ip,el)
 use material, only: &
   material_phase, &
   phase_chemicalFEInstance, &
   material_maxNcomponents

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(material_maxNcomponents,material_maxNcomponents) :: &
   chemicalFE_quadenergy_getConcTangent
 integer(pInt) :: &
   cpI, cpJ, &
   phase, &
   instance    

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 chemicalFE_quadenergy_getConcTangent = 0.0_pReal
 do cpI = 1_pInt, chemicalFE_quadenergy_Ncomponents(instance)
   do cpJ = 1_pInt, chemicalFE_quadenergy_Ncomponents(instance)
     chemicalFE_quadenergy_getConcTangent(cpI,cpJ) = &
       chemicalFE_quadenergy_QuadCoeffInv(cpI,cpJ,instance)
   enddo
 enddo

end function chemicalFE_quadenergy_getConcTangent


!--------------------------------------------------------------------------------------------------
!> @brief returns the component mobility for a given instance of this model
!--------------------------------------------------------------------------------------------------
function chemicalFE_quadenergy_getMobility(ipc,ip,el)
 use material, only: &
   material_phase, &
   phase_chemicalFEInstance, &
   material_maxNcomponents

 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal), dimension(material_maxNcomponents) :: &
   chemicalFE_quadenergy_getMobility
 integer(pInt) :: &
   cp, &
   phase, &
   instance    

 phase = material_phase(ipc,ip,el)
 instance = phase_chemicalFEInstance(phase)
 chemicalFE_quadenergy_getMobility = 0.0_pReal
 do cp = 1_pInt, chemicalFE_quadenergy_Ncomponents(instance)
   chemicalFE_quadenergy_getMobility(cp) = &
     chemicalFE_quadenergy_Mobility(cp,instance)
 enddo

end function chemicalFE_quadenergy_getMobility


!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function chemicalFE_quadenergy_postResults(conc,ipc,ip,el)
 use material, only: &
   material_maxNcomponents, &
   material_phase, &
   phase_chemicalFEInstance

 implicit none
 real(pReal), dimension(material_maxNcomponents), intent(in) :: &
   conc
 integer(pInt),                                   intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element                                                                                        !< microstructure state

 real(pReal), dimension(chemicalFE_quadenergy_sizePostResults(phase_chemicalFEInstance(material_phase(ipc,ip,el)))) :: &
   chemicalFE_quadenergy_postResults
 real(pReal), dimension(material_maxNcomponents) :: &
   tempPerComponent

 integer(pInt) :: &
   instance, &
   o,c

 instance = phase_chemicalFEInstance(material_phase(ipc,ip,el))
 chemicalFE_quadenergy_postResults = 0.0_pReal
 c = 0_pInt

 outputsLoop: do o = 1_pInt,chemicalFE_quadenergy_Noutput(instance)
   select case(chemicalFE_quadenergy_outputID(o,instance))
     case (chemicalFE_ID)
       chemicalFE_quadenergy_postResults(c+1_pInt) = chemicalFE_quadenergy_getEnergy(conc,ipc,ip,el)
       c = c + 1_pInt

     case (chemicalPot_ID)
       tempPerComponent = chemicalFE_quadenergy_getChemPot(conc,ipc,ip,el)
       chemicalFE_quadenergy_postResults(c+1_pInt:c+chemicalFE_quadenergy_Ncomponents(instance)) = &
         tempPerComponent(1_pInt:chemicalFE_quadenergy_Ncomponents(instance))
       c = c + chemicalFE_quadenergy_Ncomponents(instance)

     case (chemicalConc_ID)
       chemicalFE_quadenergy_postResults(c+1_pInt:c+chemicalFE_quadenergy_Ncomponents(instance)) = &
         conc(1_pInt:chemicalFE_quadenergy_Ncomponents(instance))
       c = c + chemicalFE_quadenergy_Ncomponents(instance)

   end select
 enddo outputsLoop

end function chemicalFE_quadenergy_postResults

end module chemicalFE_quadenergy
