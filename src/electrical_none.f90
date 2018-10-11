!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Mixture homogenization models
!--------------------------------------------------------------------------------------------------
module electrical_none
 use prec, only: &
   pInt, &
   pReal
 
 implicit none
 private
 integer(pInt),             dimension(:),     allocatable,         public, protected :: &
   electrical_none_sizePostResults
 integer(pInt),             dimension(:,:),   allocatable, target, public :: &
   electrical_none_sizePostResult
 
 character(len=64),         dimension(:,:),   allocatable, target, public :: &
  electrical_none_output                                                                                !< name of each post result output
 integer(pInt),             dimension(:),     allocatable, target, public :: &
   electrical_none_Noutput                                                                              !< number of outputs per homog instance
   
 public :: &
   electrical_none_init

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine electrical_none_init
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
 integer(pInt) :: sizeState  
 integer :: &
   maxNinstance, &
   homog, &
   instance
 integer :: &
   NofmyElectrical                                                                                       ! no pInt (stores a system dependen value from 'count'
   
 write(6,'(/,a)')   ' <<<+-  electrical_'//ELECTRICAL_none_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = count(electrical_type == ELECTRICAL_none_ID)
 if (maxNinstance == 0) return
 
 if (iand(debug_level(debug_HOMOGENIZATION),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 allocate(electrical_none_sizePostResults(maxNinstance),                              source=0_pInt)
 allocate(electrical_none_sizePostResult(maxval(homogenization_Noutput),maxNinstance),source=0_pInt)
 allocate(electrical_none_Noutput        (maxNinstance),                              source=0_pInt)
 allocate(electrical_none_output        (maxval(homogenization_Noutput),maxNinstance))
          electrical_none_output = ''
          
  initializeInstances: do homog = 1_pInt, material_Nhomogenization
  myElectrical: if (electrical_type(homog) == ELECTRICAL_none_ID) then
     NofmyElectrical = count(material_homog == homog)
     instance = electrical_typeInstance(homog)

     ! allocate state arrays
     sizeState = 0_pInt
     electricalState(homog)%sizeState = sizeState
     electricalState(homog)%sizePostResults =  electrical_none_sizePostResults(instance)
     allocate(electricalState(homog)%state0   (sizeState,NofmyElectrical))
     allocate(electricalState(homog)%subState0(sizeState,NofmyElectrical))
     allocate(electricalState(homog)%state    (sizeState,NofmyElectrical))
     
     electricPotentialMapping(homog)%p => mappingHomogenizationConst
     allocate  (electricPotential(homog)%p(1), source = 0.0_pReal)

     
   endif myElectrical
 enddo initializeInstances
 
end subroutine electrical_none_init

end module electrical_none