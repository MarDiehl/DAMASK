!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Mixture homogenization models
!--------------------------------------------------------------------------------------------------
module electrical_conduction
 use prec, only: &
   pInt, &
   pReal
 
 implicit none
 private
 integer(pInt),             dimension(:,:),   allocatable, target, public :: &
   electrical_conduction_sizePostResult
 
 character(len=64),         dimension(:,:),   allocatable, target, public :: &
  electrical_conduction_output                                                                                !< name of each post result output
 integer(pInt),             dimension(:),     allocatable, target, public :: &
  electrical_conduction_Noutput                                                                              !< number of outputs per homog instance

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 avgcurrentdensity_ID
 end enum

 type, private :: tParameters                                                                       !< container type for internal constitutive parameters
   integer(kind(undefined_ID)), dimension(:), allocatable :: & 
     outputID
 end type

 type(tParameters), dimension(:), allocatable, private :: param                                     !< containers of constitutive parameters (len Ninstance)

 public :: &                                                                                        !<some of these are going to be redundant            
   electrical_conduction_init, &
   electrical_conduction_getFlux, &  
   electrical_conduction_getFluxTangent, &
   electrical_conduction_calAndPutCurrentDensity, &
   electrical_conduction_getavgElectricalField_from_currentDensity, &
   electrical_conduction_postResults

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine electrical_conduction_init(fileUnit)
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
   instance, &
   outputSize
 integer :: &
   NofmyElectrical                                                                                       ! no pInt (stores a system dependen value from 'count'
 character(len=65536), dimension(:), allocatable :: outputs
 character(len=65536), dimension(0), parameter :: emptyStringArray = [character(len=65536)::]
 integer(kind(undefined_ID)) :: outputID  
 
 write(6,'(/,a)')   ' <<<+-  electrical_'//ELECTRICAL_conduction_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp() 
#include "compilation_info.f90"

 maxNinstance = count(electrical_type == ELECTRICAL_conduction_ID)
 if (maxNinstance == 0_pInt) return
 
 if (iand(debug_level(debug_HOMOGENIZATION),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 allocate(electrical_conduction_sizePostResult(maxval(homogenization_Noutput),maxNinstance),source=0_pInt)
 allocate(electrical_conduction_Noutput        (maxNinstance),                              source=0_pInt)
 allocate(electrical_conduction_output        (maxval(homogenization_Noutput),maxNinstance))
          electrical_conduction_output = ''
 
 allocate(param(maxNinstance))

 parsingFile:do homog = 1_pInt, material_Nhomogenization
   if (electrical_type(homog) /= ELECTRICAL_conduction_ID) cycle
   instance = electrical_typeInstance(homog)
      
   outputs = config_homogenization(homog)%getStrings('(output)',defaultVal=emptyStringArray)
   allocate(param(instance)%outputID(0))   
   do o=1_pInt, size(outputs)
     outputID = undefined_ID
     select case(outputs(o))
       case ('avgcurrentdensity')
        outputID = avgcurrentdensity_ID
        outputSize = 3_pInt
     end select 
       
     if (outputID /= undefined_ID) then
      electrical_conduction_output(o,instance) = outputs(o)
      electrical_conduction_sizePostResult(o,instance) = outputSize
      param(instance)%outputID = [param(instance)%outputID , outputID]
     endif
   enddo
 enddo parsingFile
 
 initializeInstances: do homog = 1_pInt, material_Nhomogenization
   myElectrical: if (electrical_type(homog) == ELECTRICAL_conduction_ID) then
     NofmyElectrical = count(material_homog == homog)
     instance = electrical_typeInstance(homog)

     ! allocate state arrays
     sizeState = 0_pInt
     electricalState(homog)%sizeState = sizeState
     electricalState(homog)%sizePostResults = sum(electrical_conduction_sizePostResult(:,instance))
     allocate(electricalState(homog)%state0   (sizeState,NofmyElectrical))
     allocate(electricalState(homog)%subState0(sizeState,NofmyElectrical))
     allocate(electricalState(homog)%state    (sizeState,NofmyElectrical))

     
   endif myElectrical
 enddo initializeInstances
 
end subroutine electrical_conduction_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates and returns current density at material point 
!--------------------------------------------------------------------------------------------------
subroutine electrical_conduction_getFlux(CurrentDensity,ElectricalField,ip,el)
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   material_phase, &
   phase_currentDensity, &
   CURRENTDENSITY_none_ID, &
   CURRENTDENSITY_ohm_ID, &
   phasefracMapping, &
   phasefrac
 use currentDensity_ohm, only: &
   currentDensity_ohm_getFlux  
 
 implicit none
 integer(pInt),             intent(in)  :: &
   ip, el                                                                                      !< element number
 real(pReal), dimension(3), intent(in)  :: &
   ElectricalField
 real(pReal), dimension(3), intent(out) :: &
   CurrentDensity
 real(pReal), dimension(3) :: &
   CurrentDensity_local
 integer(pInt) :: &
   gr, &
   homog, & 
   offset

 homog = material_homog(ip,el)
 offset = phasefracMapping(homog)%p(ip,el)
 CurrentDensity = 0.0_pReal
 do gr = 1_pInt, homogenization_Ngrains(homog)
   currentDensityType: select case (phase_currentDensity(material_phase(gr,ip,el)))
     case (CURRENTDENSITY_none_ID) currentDensityType
       CurrentDensity_local = 0.0_pReal

     case (CURRENTDENSITY_ohm_ID) currentDensityType
       call currentDensity_ohm_getFlux(CurrentDensity_local,ElectricalField,gr,ip,el)

     case default
       CurrentDensity_local = 0.0_pReal

   end select currentDensityType
   CurrentDensity = CurrentDensity + phasefrac(homog)%p(gr,offset)*CurrentDensity_local
 enddo

end subroutine electrical_conduction_getFlux


!--------------------------------------------------------------------------------------------------
!> @brief calculates and returns conductivity at material point 
!--------------------------------------------------------------------------------------------------
subroutine electrical_conduction_getFluxTangent(CurrentDensityTangent,ip,el)
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   material_phase, &
   phase_currentDensity, &
   CURRENTDENSITY_none_ID, &
   CURRENTDENSITY_ohm_ID, &
   phasefracMapping, &
   phasefrac
 use currentDensity_ohm, only: &
   currentDensity_ohm_getFluxTangent  
 
 implicit none
 integer(pInt),               intent(in)  :: &
   ip, el                                                                                      !< element number
 real(pReal), dimension(3,3), intent(out) :: &
   CurrentDensityTangent
 real(pReal), dimension(3,3) :: &
   CurrentDensityTangent_local
 integer(pInt) :: &
   gr, &
   homog, & 
   offset

 homog = material_homog(ip,el)
 offset = phasefracMapping(homog)%p(ip,el)
 CurrentDensityTangent = 0.0_pReal
 
 do gr = 1_pInt, homogenization_Ngrains(homog)
   currentDensityType: select case (phase_currentDensity(material_phase(gr,ip,el)))
     case (CURRENTDENSITY_none_ID) currentDensityType
       CurrentDensityTangent_local = 0.0_pReal

     case (CURRENTDENSITY_ohm_ID) currentDensityType
       call currentDensity_ohm_getFluxTangent(CurrentDensityTangent_local,gr,ip,el)

     case default
       CurrentDensityTangent_local = 0.0_pReal

   end select currentDensityType
   CurrentDensityTangent = CurrentDensityTangent + phasefrac(homog)%p(gr,offset)*CurrentDensityTangent_local
 enddo

end subroutine electrical_conduction_getFluxTangent





!--------------------------------------------------------------------------------------------------
!> @brief store current density at material point 
!--------------------------------------------------------------------------------------------------
subroutine electrical_conduction_calAndPutCurrentDensity(ElectricalField,ip,el)
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   material_phase, &
   phase_currentDensity, &
   CURRENTDENSITY_none_ID, &
   CURRENTDENSITY_ohm_ID, &
   currentDensity, &
   currentDensityMapping
  
use currentDensity_ohm, only: &
   currentDensity_ohm_getFlux  
 
 implicit none
 integer(pInt),             intent(in)  :: &
   ip, el                                                                                      !< element number
 real(pReal), dimension(3), intent(in)  :: &
   ElectricalField
 real(pReal), dimension(3)              :: &
   currentDensityLocal
   
  integer(pInt) :: &
    gr, &
    homog

   homog = material_homog(ip,el)
   do gr = 1_pInt, homogenization_Ngrains(homog)
    currentDensityType: select case (phase_currentDensity(material_phase(gr,ip,el)))
      case (CURRENTDENSITY_none_ID) currentDensityType
        currentDensityLocal = 0.0_pReal
      
      case (CURRENTDENSITY_ohm_ID) currentDensityType
        call currentDensity_ohm_getFlux(CurrentDensityLocal,ElectricalField,gr,ip,el)
              
      case default 
        currentdensityLocal = 0.0_pReal
       
 
     end select currentDensityType
     
     currentDensity(material_phase(gr,ip,el))% &
       p(1:3,currentDensityMapping(material_phase(gr,ip,el))%p(gr,ip,el)) = currentDensityLocal
     
   enddo

end subroutine electrical_conduction_calAndPutCurrentDensity



!--------------------------------------------------------------------------------------------------
!> @brief returns the average electrical current from stored current density
!--------------------------------------------------------------------------------------------------
subroutine electrical_conduction_getavgElectricalField_from_currentDensity(avgElfield, ip, el)
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   material_phase, &
   phase_currentDensity, &
   CURRENTDENSITY_none_ID, &
   CURRENTDENSITY_ohm_ID, &
   phasefracMapping, &
   phasefrac, & 
   currentDensity, &
   currentDensityMapping
   
 use currentDensity_ohm, only: &
   currentDensity_ohm_getElectricField
   
 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(3), intent(out) :: &
   avgElfield
 real(pReal), dimension(3)              :: &
   currentDensityLocal, & 
   Elfield_local
 integer(pInt) :: &
   gr, &
   homog, &
   offset
     
 homog      = material_homog(ip,el)
 offset     = phasefracMapping(homog)%p(ip,el)
 avgElfield = 0.0_pReal;

 do gr = 1_pInt, homogenization_Ngrains(homog)
  currentDensityLocal = currentDensity(material_phase(gr,ip,el))% &
       p(1:3,currentDensityMapping(material_phase(gr,ip,el))%p(gr,ip,el))
       
  currentDensityType: select case (phase_currentDensity(material_phase(gr,ip,el)))
     case (CURRENTDENSITY_none_ID) currentDensityType
       Elfield_local = 0.0_pReal

     case (CURRENTDENSITY_ohm_ID) currentDensityType
       call currentDensity_ohm_getElectricField(Elfield_local,currentDensityLocal,gr,ip,el)

     case default
       Elfield_local = 0.0_pReal

   end select currentDensityType
   avgElfield = avgElfield + phasefrac(homog)%p(gr,offset)*Elfield_local 
 
 enddo
 
 end subroutine electrical_conduction_getavgElectricalField_from_currentDensity
 
 
!--------------------------------------------------------------------------------------------------
!> @brief return array of homogenization results for post file inclusion 
!--------------------------------------------------------------------------------------------------
function electrical_conduction_postResults(ip,el)
 use material, only: &
   phasefrac, &
   phasefracMapping, &
   currentDensity, &
   currentDensityMapping, &
   material_homog, &
   material_phase, &
   homogenization_Ngrains, &
   electrical_typeInstance, &
   homogenization_Noutput
 
 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),  dimension(sum(electrical_conduction_sizePostResult(:,electrical_typeInstance(material_homog(ip,el))))) :: &
   electrical_conduction_postResults
 
 integer(pInt) :: &
   homID, &
   o, c, gr, offset
   
 c = 0_pInt
 homID = electrical_typeInstance(material_homog(ip,el))
 
 offset = phasefracMapping(material_homog(ip,el))%p(ip,el)
 electrical_conduction_postResults = 0.0_pReal
 
 do o = 1_pInt,size(param(homID)%outputID)
   select case(param(homID)%outputID(o))
     case (avgcurrentdensity_ID)
       electrical_conduction_postResults(c+1_pInt:c+3_pInt) = 0.0_pReal
       do gr = 1_pInt, homogenization_Ngrains(material_homog(ip,el))
         electrical_conduction_postResults(c+1_pInt:c+3_pInt) = &
           electrical_conduction_postResults(c+1_pInt:c+3_pInt) + &
           phasefrac(material_homog(ip,el))%p(gr,offset)* &
           currentDensity(material_phase(gr,ip,el))%p(1:3,currentDensityMapping(material_phase(gr,ip,el))%p(gr,ip,el))
       enddo
       c = c + 3

   end select
 enddo

end function electrical_conduction_postResults

end module electrical_conduction
