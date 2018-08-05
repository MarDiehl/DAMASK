!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for local temperature evolution
!--------------------------------------------------------------------------------------------------
module thermal_local
 use prec, only: &
   pInt, &
   pReal
 
 implicit none
 private
 integer(pInt),             dimension(:),     allocatable,         public, protected :: &
   thermal_local_sizePostResults
 integer(pInt),             dimension(:,:),   allocatable, target, public :: &
   thermal_local_sizePostResult
 
 character(len=64),         dimension(:,:),   allocatable, target, public :: &
  thermal_local_output                                                                                !< name of each post result output
 integer(pInt),             dimension(:),     allocatable, target, public :: &
   thermal_local_Noutput                                                                              !< number of outputs per homog instance

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 avgtemperature_ID
 end enum

 type, private :: tParameters                                                                       !< container type for internal constitutive parameters
   integer(kind(undefined_ID)), dimension(:), allocatable :: & 
     outputID
 end type

 type(tParameters), dimension(:), allocatable, private :: param                                     !< containers of constitutive parameters (len Ninstance)

 public :: &
   thermal_local_init, &
   thermal_local_getThermalViscosity, &  
   thermal_local_getHeatSource, &  
   thermal_local_putTemperatureRate, &
   thermal_local_postResults

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine thermal_local_init(fileUnit)
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use debug, only: &
   debug_HOMOGENIZATION, &
   debug_level, &
   debug_levelBasic
 use IO
 use material
 use config
 
 implicit none
 integer(pInt),                intent(in) :: fileUnit
 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: &
   section = 0_pInt, i, mySize = 0_pInt, o
 integer(pInt) :: sizeState
 integer :: &
   maxNinstance, &
   homog, &
   instance
 integer :: &
   NofmyThermal                                                                                       ! no pInt (stores a system dependen value from 'count'
 character(len=65536) :: &
   tag  = '', &
   line = ''
 
 write(6,'(/,a)')   ' <<<+-  thermal_'//THERMAL_local_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = count(thermal_type == THERMAL_local_ID)
 if (maxNinstance == 0) return
 
 if (iand(debug_level(debug_HOMOGENIZATION),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 allocate(thermal_local_sizePostResults(maxNinstance),                              source=0_pInt)
 allocate(thermal_local_sizePostResult(maxval(homogenization_Noutput),maxNinstance),source=0_pInt)
 allocate(thermal_local_Noutput        (maxNinstance),                              source=0_pInt)
 allocate(thermal_local_output        (maxval(homogenization_Noutput),maxNinstance))
          thermal_local_output = ''
 
 allocate(param(maxNinstance))

 rewind(fileUnit)
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= material_partHomogenization)! wind forward to <homogenization>
   line = IO_read(fileUnit)
 enddo

 parsingFile: do while (trim(line) /= IO_EOF)                                                       ! read through sections of homogenization part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit                                                                                           
   endif
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt
     if (thermal_type(section) == THERMAL_local_ID) then
       i = thermal_typeInstance(section)                                                            ! count instances of my homogenization law
       allocate(param(i)%outputID(homogenization_Noutput(section)))                                 ! allocate space for IDs of every requested output
     endif
     cycle
   endif
   if (section > 0_pInt ) then                                                                      ! do not short-circuit here (.and. with next if-statement). It's not safe in Fortran
     if (thermal_type(section) == THERMAL_local_ID) then                         ! one of my sections
       i = thermal_typeInstance(section)                                                     ! which instance of my type is present homogenization
       chunkPos = IO_stringPos(line)
       tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                            ! extract key
       select case(tag)
         case ('(output)')
           select case(IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
             case('avgtemperature')
               thermal_local_Noutput(i) = thermal_local_Noutput(i) + 1_pInt
               thermal_local_output(thermal_local_Noutput(i),i) = &
                 IO_lc(IO_stringValue(line,chunkPos,2_pInt))
               param(i)%outputID(thermal_local_Noutput(i)) = avgtemperature_ID

           end select

       end select
     endif
   endif
 enddo parsingFile

 initializeInstances: do homog = 1_pInt, material_Nhomogenization
   myThermal: if (thermal_type(homog) == THERMAL_local_ID) then
     NofmyThermal = count(material_homog == homog)
     instance = thermal_typeInstance(homog)

! *  Determine size of postResults array
     outputsLoop: do o = 1_pInt, thermal_local_Noutput(instance)
       select case(param(instance)%outputID(o))
        case(avgtemperature_ID)
          mySize = 1_pInt
        case default
          mySize = 0_pInt

       end select

       outputFound: if (mySize > 0_pInt) then
        thermal_local_sizePostResult(o,instance) = mySize
        thermal_local_sizePostResults(instance) = &
          thermal_local_sizePostResults(instance) + mySize
       endif outputFound
     enddo outputsLoop

! allocate state arrays
     sizeState = 0_pInt
     thermalState(section)%sizeState = sizeState
     thermalState(section)%sizePostResults = thermal_local_sizePostResults(instance)
     allocate(thermalState(section)%state0   (sizeState,NofmyThermal), source=0.0_pReal)
     allocate(thermalState(section)%subState0(sizeState,NofmyThermal), source=0.0_pReal)
     allocate(thermalState(section)%state    (sizeState,NofmyThermal), source=0.0_pReal)

   endif myThermal
 enddo initializeInstances
 
end subroutine thermal_local_init


!--------------------------------------------------------------------------------------------------
!> @brief return average thermal viscosity at material point
!--------------------------------------------------------------------------------------------------
function thermal_local_getThermalViscosity(ip,el)
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   material_phase, &
   phasefracMapping, &
   phasefrac, &
   phase_heatflux, &
   HEATFLUX_adiabaticnone_ID, &
   HEATFLUX_joule_ID
 use heatflux_adiabaticnone, only: &
   heatflux_adiabaticnone_getThermalViscosity  
 use heatflux_joule, only: &
   heatflux_joule_getThermalViscosity  
 
 implicit none
 integer(pInt), intent(in) :: &
   ip, el                                                                                      !< element number
 real(pReal) :: &
   thermal_local_getThermalViscosity, &
   localThermalViscosity
 integer(pInt) :: &
   gr, &
   homog

 homog = material_homog(ip,el)
 thermal_local_getThermalViscosity = 0.0_pReal
 do gr = 1_pInt, homogenization_Ngrains(homog)
   heatfluxType: select case (phase_heatflux(material_phase(gr,ip,el)))
     case (HEATFLUX_adiabaticnone_ID) heatfluxType
       localThermalViscosity = heatflux_adiabaticnone_getThermalViscosity(gr,ip,el)

     case (HEATFLUX_joule_ID) heatfluxType
       localThermalViscosity = heatflux_joule_getThermalViscosity(gr,ip,el)

     case default
       localThermalViscosity = 0.0_pReal
       
   end select heatfluxType
   
   thermal_local_getThermalViscosity = &
       thermal_local_getThermalViscosity + &
       phasefrac(homog)%p(gr,phasefracMapping(homog)%p(ip,el))*localThermalViscosity
 enddo

end function thermal_local_getThermalViscosity


!--------------------------------------------------------------------------------------------------
!> @brief return average heat source at material point
!--------------------------------------------------------------------------------------------------
function thermal_local_getHeatSource(ip,el)
 use math, only: &
   math_Mandel6to33
 use material, only: &
   homogenization_Ngrains, &
   material_homog, &
   phasefracMapping, &
   phasefrac, &
   phaseAt, &
   phase_Nsources, &
   phase_source, &
   SOURCE_thermal_dissipation_ID, &
   SOURCE_thermal_externalheat_ID
 use source_thermal_dissipation, only: &
   source_thermal_dissipation_getRateAndItsTangent
 use source_thermal_externalheat, only: &
   source_thermal_externalheat_getRateAndItsTangent
 use crystallite, only: &
   crystallite_Tstar_v, &
   crystallite_Lp  

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal) :: &
   thermal_local_getHeatSource, &
   localHeatSource, &
   localHeatSourceTangent
 integer(pInt) :: &
   phase, &
   homog, &
   grain, &
   source
   
 homog  = material_homog(ip,el)  
 thermal_local_getHeatSource = 0.0_pReal
 do grain = 1, homogenization_Ngrains(homog)
   phase = phaseAt(grain,ip,el)
   do source = 1, phase_Nsources(phase)
     select case(phase_source(source,phase))                                                   
       case (SOURCE_thermal_dissipation_ID)
        call source_thermal_dissipation_getRateAndItsTangent(localHeatSource, localHeatSourceTangent, &
                                                             crystallite_Tstar_v(1:6,grain,ip,el), &
                                                             crystallite_Lp(1:3,1:3,grain,ip,el), &
                                                             grain, ip, el)

       case (SOURCE_thermal_externalheat_ID)
        call source_thermal_externalheat_getRateAndItsTangent(localHeatSource, localHeatSourceTangent, &
                                                              grain, ip, el)

       case default
        localHeatSource = 0.0_pReal

     end select
     thermal_local_getHeatSource = &
       thermal_local_getHeatSource + &
       phasefrac(homog)%p(grain,phasefracMapping(homog)%p(ip,el))*localHeatSource
   enddo  
 enddo
 
end function thermal_local_getHeatSource


!--------------------------------------------------------------------------------------------------
!> @brief set module wide temperature rate 
!--------------------------------------------------------------------------------------------------
subroutine thermal_local_putTemperatureRate(ip,el)
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   material_phase, &
   temperatureRate, &
   thermalMapping
 implicit none
 integer(pInt),                                         intent(in) :: &
   ip, el                                                                                      !< element number
 integer(pInt) :: &
   gr, &
   phase, &
   homog

 homog = material_homog(ip,el)
 do gr = 1_pInt, homogenization_Ngrains(homog)
   phase = material_phase(gr,ip,el)
   temperatureRate(phase)%p(thermalMapping(phase)%p(gr,ip,el)) = &
     thermal_local_getHeatSource(ip,el)/thermal_local_getThermalViscosity(ip,el)
 enddo

end subroutine thermal_local_putTemperatureRate


!--------------------------------------------------------------------------------------------------
!> @brief return array of homogenization results for post file inclusion 
!--------------------------------------------------------------------------------------------------
function thermal_local_postResults(ip,el)
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   material_phase, &
   phasefrac, &
   phasefracMapping, &
   temperature, &
   thermalMapping, &
   thermal_typeInstance, &
   homogenization_Noutput
 
 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),  dimension(thermal_local_sizePostResults(thermal_typeInstance(material_homog(ip,el)))) :: &
   thermal_local_postResults
 
 integer(pInt) :: &
   homID, &
   o, c, gr
   
 c = 0_pInt
 homID = thermal_typeInstance(material_homog(ip,el))
 thermal_local_postResults = 0.0_pReal
 
 do o = 1_pInt,homogenization_Noutput(material_homog(ip,el))
   select case(param(homID)%outputID(o))
     case (avgtemperature_ID)
       do gr = 1_pInt, homogenization_Ngrains(material_homog(ip,el))
         thermal_local_postResults(c+1_pInt) = &
           thermal_local_postResults(c+1_pInt) + &
           phasefrac(material_homog(ip,el))%p(gr,phasefracMapping(material_homog(ip,el))%p(ip,el))* &
           temperature(material_phase(gr,ip,el))%p(thermalMapping(material_phase(gr,ip,el))%p(gr,ip,el))
       enddo
       c = c + 1_pInt

   end select
 enddo

end function thermal_local_postResults

end module thermal_local
