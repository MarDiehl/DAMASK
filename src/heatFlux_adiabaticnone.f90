!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for adiabatic heat flux
!--------------------------------------------------------------------------------------------------
module heatflux_adiabaticnone
 use prec, only: &
   pReal,&
   pInt, &
   p_vec

 implicit none
 private
 integer(pInt),                       dimension(:),     allocatable,         public, protected :: &
   heatflux_adiabaticnone_sizePostResults                                                       !< cumulative size of post results

 integer(pInt),                       dimension(:,:),   allocatable, target, public :: &
   heatflux_adiabaticnone_sizePostResult                                                        !< size of each post result output

 character(len=64),                   dimension(:,:),   allocatable, target, public :: &
   heatflux_adiabaticnone_output                                                                !< name of each post result output

 integer(pInt),                       dimension(:),     allocatable, target, public :: &
   heatflux_adiabaticnone_Noutput                                                               !< number of outputs per instance of this constitution

 enum, bind(c)
   enumerator :: undefined_ID, &
                 temperature_ID
 end enum

 type, private :: tParameters                                                                         !< container type for internal constitutive parameters
   integer(kind(undefined_ID)), dimension(:), allocatable :: & 
     outputID
   real(pReal) :: &
     atol =   0.01_pReal, &
     T0   = 300.00_pReal, &
     rho  =   0.00_pReal, &
     Cp   =   0.00_pReal
 end type

 type(tParameters), dimension(:), allocatable, private :: param                                       !< containers of constitutive parameters (len Ninstance)

 public :: &
   heatflux_adiabaticnone_init, &
   heatflux_adiabaticnone_dotState, &
   heatflux_adiabaticnone_getThermalViscosity, &
   heatflux_adiabaticnone_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine heatflux_adiabaticnone_init(fileUnit)
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
   thermalMapping, &
   temperature, &
   temperatureRate, &
   phasememberAt, &
   phase_heatFlux, &
   phase_heatFluxInstance, &
   phase_Noutput, &
   HEATFLUX_ADIABATICNONE_label, &
   HEATFLUX_ADIABATICNONE_ID, &
   material_phase, &
   heatfluxState
 use config, only: &
   MATERIAL_partPhase
 use numerics,only: &
   numerics_integrator

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: &
   maxNinstance, &
   instance,phase,o, &
   NipcMyPhase, &
   mySize=0_pInt,sizeState,sizeDotState, sizeDeltaState
 character(len=65536) :: &
   tag  = '', &
   line = ''
 integer(kind(undefined_ID)), dimension(:,:), allocatable :: & 
   tempOutputID

 write(6,'(/,a)')   ' <<<+-  heatflux_'//HEATFLUX_ADIABATICNONE_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = int(count(phase_heatflux == HEATFLUX_ADIABATICNONE_ID),pInt)
 if (maxNinstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance


 allocate(heatflux_adiabaticnone_sizePostResults(maxNinstance),                   source=0_pInt)
 allocate(heatflux_adiabaticnone_sizePostResult(maxval(phase_Noutput),maxNinstance), &
                                                                                 source=0_pInt)
 allocate(heatflux_adiabaticnone_output(maxval(phase_Noutput),maxNinstance))
          heatflux_adiabaticnone_output               = ''
 allocate(heatflux_adiabaticnone_Noutput(maxNinstance),                       source=0_pInt)
 
 allocate(param(maxNinstance))
 
 allocate(tempOutputID        (maxval(phase_Noutput)  ,maxNinstance))

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
   if (phase > 0_pInt ) then; if (phase_heatflux(phase) == HEATFLUX_ADIABATICNONE_ID) then      ! one of my phases. Do not short-circuit here (.and. between if-statements), it's not safe in Fortran
     instance = phase_heatfluxInstance(phase)                                                     ! which instance of my heatflux is present phase
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case ('temperature')
             heatflux_adiabaticnone_Noutput(instance) = heatflux_adiabaticnone_Noutput(instance) + 1_pInt
             heatflux_adiabaticnone_output(heatflux_adiabaticnone_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))
             tempOutputID(heatflux_adiabaticnone_Noutput(instance),instance) = temperature_ID

           case default

         end select

!--------------------------------------------------------------------------------------------------
! parameters depending on number of components
       case ('temperature_atol')
         param(instance)%aTol = IO_floatValue(line,chunkPos,2_pInt)
         
       case ('t0')
         param(instance)%T0 = IO_floatValue(line,chunkPos,2_pInt)
         
       case ('mass_density')
         param(instance)%rho = IO_floatValue(line,chunkPos,2_pInt)
         
       case ('specific_heat')
         param(instance)%Cp = IO_floatValue(line,chunkPos,2_pInt)
         
       case default

     end select
   endif; endif
 enddo parsingFile

 sanityChecks: do phase = 1_pInt, size(phase_heatflux)
   myPhase: if (phase_heatflux(phase) == HEATFLUX_adiabaticnone_ID) then
     instance = phase_heatfluxInstance(phase)
     allocate(param(instance)%OutputID         (heatflux_adiabaticnone_Noutput(instance)))
     param(instance)%OutputID = tempOutputID(1:heatflux_adiabaticnone_Noutput(instance),instance)

     if (param(instance)%aTol <= 0.0_pReal) &
       call IO_error(211_pInt,el=instance,ext_msg='tolerance ('//HEATFLUX_ADIABATICNONE_label//')')
     if (param(instance)%T0 <= 0.0_pReal) &
       call IO_error(211_pInt,el=instance,ext_msg='T0 ('//HEATFLUX_ADIABATICNONE_label//')')
     if (param(instance)%rho <= 0.0_pReal) &
       call IO_error(211_pInt,el=instance,ext_msg='mass density ('//HEATFLUX_ADIABATICNONE_label//')')
     if (param(instance)%Cp <= 0.0_pReal) &
       call IO_error(211_pInt,el=instance,ext_msg='specific heat capacity ('//HEATFLUX_ADIABATICNONE_label//')')
   endif myPhase
 enddo sanityChecks

 initializeInstances: do phase = 1_pInt, size(phase_heatflux)                                     ! loop through all phases in material.config
   myPhase2: if (phase_heatflux(phase) == HEATFLUX_adiabaticnone_ID) then                       ! only consider my phase
     NipcMyPhase = count(material_phase == phase)                                                   ! number of IPCs containing my phase
     instance = phase_heatfluxInstance(phase)                                                     ! which instance of my phase

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,heatflux_adiabaticnone_Noutput(instance)
       select case(param(instance)%outputID(o))
         case(temperature_ID)
           mySize = 1_pInt
         case default
       end select

       outputFound: if (mySize > 0_pInt) then
         heatflux_adiabaticnone_sizePostResult(o,instance) = mySize
         heatflux_adiabaticnone_sizePostResults(instance)  = heatflux_adiabaticnone_sizePostResults(instance) + mySize
       endif outputFound
     enddo outputsLoop
!--------------------------------------------------------------------------------------------------
! allocate state arrays
     sizeState = 1_pInt
     sizeDotState = sizeState
     sizeDeltaState = 0_pInt
     heatfluxState(phase)%sizeState = sizeState
     heatfluxState(phase)%sizeDotState = sizeDotState
     heatfluxState(phase)%sizeDeltaState = sizeDeltaState
     heatfluxState(phase)%sizePostResults = heatflux_adiabaticnone_sizePostResults(instance)
     allocate(heatfluxState(phase)%aTolState          (   sizeState),             source=param(instance)%aTol)
     allocate(heatfluxState(phase)%state0             (   sizeState,NipcMyPhase), source=param(instance)%T0)
     allocate(heatfluxState(phase)%partionedState0    (   sizeState,NipcMyPhase), source=param(instance)%T0)
     allocate(heatfluxState(phase)%subState0          (   sizeState,NipcMyPhase), source=param(instance)%T0)
     allocate(heatfluxState(phase)%state              (   sizeState,NipcMyPhase), source=param(instance)%T0)
     allocate(heatfluxState(phase)%dotState           (sizeDotState,NipcMyPhase), source=0.0_pReal)
     allocate(heatfluxState(phase)%deltaState       (sizeDeltaState,NipcMyPhase), source=0.0_pReal)
     if (any(numerics_integrator == 1_pInt)) then
       allocate(heatfluxState(phase)%previousDotState (sizeDotState,NipcMyPhase),source=0.0_pReal)
       allocate(heatfluxState(phase)%previousDotState2(sizeDotState,NipcMyPhase),source=0.0_pReal)
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(heatfluxState(phase)%RK4dotState      (sizeDotState,NipcMyPhase), source=0.0_pReal)
     if (any(numerics_integrator == 5_pInt)) &
       allocate(heatfluxState(phase)%RKCK45dotState (6,sizeDotState,NipcMyPhase), source=0.0_pReal)

     thermalMapping(phase)%p => phasememberAt
     temperature(phase)%p => heatfluxState(phase)%state(1,1:NipcMyPhase) 
     allocate(temperatureRate(phase)%p(NipcMyPhase), source=0.0_pReal)
   endif myPhase2
 enddo initializeInstances   

end subroutine heatflux_adiabaticnone_init


!--------------------------------------------------------------------------------------------------
!> @brief returns the temperature rates
!--------------------------------------------------------------------------------------------------
subroutine heatflux_adiabaticnone_dotState(ipc,ip,el)
 use material, only: &
   heatfluxState, &
   temperatureRate, &
   thermalMapping, &
   material_phase
   
 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 integer(pInt) :: &
   phase

 phase = material_phase(ipc,ip,el) 
 heatfluxState(phase)%dotState(1,thermalMapping(phase)%p(ipc,ip,el)) = &
     temperatureRate(phase)%p(thermalMapping(phase)%p(ipc,ip,el))

end subroutine heatflux_adiabaticnone_dotState

!--------------------------------------------------------------------------------------------------
!> @brief returns the thermal viscosiy, i.e. rho*Cp
!--------------------------------------------------------------------------------------------------
pure function heatflux_adiabaticnone_getThermalViscosity(ipc,ip,el)
 use material, only: &
   phase_heatfluxInstance, &
   material_phase
   
 implicit none
 integer(pInt), intent(in) :: &
   ipc, ip, el
 real(pReal) :: &
   heatflux_adiabaticnone_getThermalViscosity
 integer(pInt) :: &
   phase, &
   instance

 phase = material_phase(ipc,ip,el)
 instance = phase_heatfluxInstance(phase)
 heatflux_adiabaticnone_getThermalViscosity = param(instance)%rho*param(instance)%Cp

end function heatflux_adiabaticnone_getThermalViscosity

!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function heatflux_adiabaticnone_postResults(ipc,ip,el)
 use material, only: &
   material_phase, &
   temperature, &
   thermalMapping, &
   phase_heatfluxInstance

 implicit none
 integer(pInt),                                   intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element                                                                                        !< microstructure state

 real(pReal), dimension(heatflux_adiabaticnone_sizePostResults(phase_heatfluxInstance(material_phase(ipc,ip,el)))) :: &
   heatflux_adiabaticnone_postResults
 integer(pInt) :: &
   phase, &
   instance, &
   o,c

 phase = material_phase(ipc,ip,el)
 instance = phase_heatfluxInstance(phase)
 
 heatflux_adiabaticnone_postResults = 0.0_pReal
 c = 0_pInt
 outputsLoop: do o = 1_pInt,heatflux_adiabaticnone_Noutput(instance)
   select case(param(instance)%outputID(o))
     case (temperature_ID)
       heatflux_adiabaticnone_postResults(c+1_pInt) = temperature(phase)%p(thermalMapping(phase)%p(ipc,ip,el))
       c = c + 1_pInt

   end select
 enddo outputsLoop

end function heatflux_adiabaticnone_postResults

end module heatflux_adiabaticnone
