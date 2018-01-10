!--------------------------------------------------------------------------------------------------
!> @author Kishan Govind, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Isostress homogenization scheme
!--------------------------------------------------------------------------------------------------
module homogenization_isostress
 use prec, only: &
   pInt
 
 implicit none
 private
 integer(pInt),               dimension(:),   allocatable,         public, protected :: &
   homogenization_isostress_sizePostResults
 integer(pInt),               dimension(:,:), allocatable, target, public :: &
   homogenization_isostress_sizePostResult
 
 character(len=64),           dimension(:,:), allocatable, target, public :: &
  homogenization_isostress_output                                                                   !< name of each post result output
 integer(pInt),               dimension(:),   allocatable, target, public :: &
   homogenization_isostress_Noutput                                                                 !< number of outputs per homog instance
 integer(pInt),               dimension(:),   allocatable,         private :: &
   homogenization_isostress_Ngrains
 enum, bind(c) 
   enumerator :: undefined_ID, &
                 nconstituents_ID, &
                 ipcoords_ID, &
                 avgdefgrad_ID, &
                 avgfirstpiola_ID
 end enum
 enum, bind(c) 
   enumerator :: parallel_ID, &
                 average_ID
 end enum
 integer(kind(undefined_ID)), dimension(:,:), allocatable,         private :: &
  homogenization_isostress_outputID                                                                 !< ID of each post result output


 public :: &
   homogenization_isostress_init, &
   homogenization_isostress_partitionDeformation, &
   homogenization_isostress_averageStressAndItsTangent, &
   homogenization_isostress_postResults

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine homogenization_isostress_init(fileUnit)
#ifdef __GFORTRAN__
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use prec, only: &
   pReal
 use debug, only: &
   debug_HOMOGENIZATION, &
   debug_level, &
   debug_levelBasic
 use IO
 use material
 
 implicit none
 integer(pInt),                                      intent(in) :: fileUnit
 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: &
   section = 0_pInt, i, mySize, o
 integer :: &
   maxNinstance, &
   homog, &
   instance
 integer :: &
   NofMyHomog                                                                                       ! no pInt (stores a system dependen value from 'count'
 character(len=65536) :: &
   tag  = '', &
   line = ''
 
 write(6,'(/,a)')   ' <<<+-  homogenization_'//HOMOGENIZATION_ISOSTRESS_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = count(homogenization_type == HOMOGENIZATION_ISOSTRESS_ID)
 if (maxNinstance == 0) return
 
 if (iand(debug_level(debug_HOMOGENIZATION),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 allocate(homogenization_isostress_sizePostResults(maxNinstance),          source=0_pInt)
 allocate(homogenization_isostress_sizePostResult(maxval(homogenization_Noutput),maxNinstance), &
                                                                           source=0_pInt)
 allocate(homogenization_isostress_Noutput(maxNinstance),                  source=0_pInt)
 allocate(homogenization_isostress_Ngrains(maxNinstance),                  source=0_pInt)
 allocate(homogenization_isostress_output(maxval(homogenization_Noutput),maxNinstance))
          homogenization_isostress_output = ''
 allocate(homogenization_isostress_outputID(maxval(homogenization_Noutput),maxNinstance), &
                                                                           source=undefined_ID)

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
     cycle
   endif
   if (section > 0_pInt ) then                                                                      ! do not short-circuit here (.and. with next if-statement). It's not safe in Fortran
     if (homogenization_type(section) == HOMOGENIZATION_ISOSTRESS_ID) then                          ! one of my sections
       i = homogenization_typeInstance(section)                                                     ! which instance of my type is present homogenization
       chunkPos = IO_stringPos(line)
       tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                           ! extract key
       select case(tag)
         case ('(output)')
           select case(IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
             case('nconstituents','ngrains')
               homogenization_isostress_Noutput(i) = homogenization_isostress_Noutput(i) + 1_pInt
               homogenization_isostress_outputID(homogenization_isostress_Noutput(i),i) = nconstituents_ID
               homogenization_isostress_output(homogenization_isostress_Noutput(i),i) = &
                 IO_lc(IO_stringValue(line,chunkPos,2_pInt))
             case('ipcoords')
               homogenization_isostress_Noutput(i) = homogenization_isostress_Noutput(i) + 1_pInt
               homogenization_isostress_outputID(homogenization_isostress_Noutput(i),i) = ipcoords_ID
               homogenization_isostress_output(homogenization_isostress_Noutput(i),i) = &
                 IO_lc(IO_stringValue(line,chunkPos,2_pInt))
             case('avgdefgrad','avgf')
               homogenization_isostress_Noutput(i) = homogenization_isostress_Noutput(i) + 1_pInt
               homogenization_isostress_outputID(homogenization_isostress_Noutput(i),i) = avgdefgrad_ID
               homogenization_isostress_output(homogenization_isostress_Noutput(i),i) = &
                 IO_lc(IO_stringValue(line,chunkPos,2_pInt))
             case('avgp','avgfirstpiola','avg1stpiola')
               homogenization_isostress_Noutput(i) = homogenization_isostress_Noutput(i) + 1_pInt
               homogenization_isostress_outputID(homogenization_isostress_Noutput(i),i) = avgfirstpiola_ID
               homogenization_isostress_output(homogenization_isostress_Noutput(i),i) = &
                 IO_lc(IO_stringValue(line,chunkPos,2_pInt))

           end select
         case ('nconstituents','ngrains')
           homogenization_isostress_Ngrains(i) = IO_intValue(line,chunkPos,2_pInt)

       end select
     endif
   endif
 enddo parsingFile

 initializeInstances: do homog = 1_pInt, material_Nhomogenization
   myHomog: if (homogenization_type(homog) == HOMOGENIZATION_ISOSTRESS_ID) then
     NofMyHomog = count(material_homog == homog)
     instance = homogenization_typeInstance(homog)

! *  Determine size of postResults array
     outputsLoop: do o = 1_pInt, homogenization_isostress_Noutput(instance)
       select case(homogenization_isostress_outputID(o,instance))
        case(nconstituents_ID)
          mySize = 1_pInt
        case(ipcoords_ID)
          mySize = 3_pInt
        case(avgdefgrad_ID, avgfirstpiola_ID)
          mySize = 9_pInt
        case default
          mySize = 0_pInt
       end select

       outputFound: if (mySize > 0_pInt) then
        homogenization_isostress_sizePostResult(o,instance) = mySize
        homogenization_isostress_sizePostResults(instance) = &
          homogenization_isostress_sizePostResults(instance) + mySize
       endif outputFound
     enddo outputsLoop

! allocate state arrays
     homogState(homog)%sizeState = 0_pInt
     homogState(homog)%sizePostResults = homogenization_isostress_sizePostResults(instance)
     allocate(homogState(homog)%state0   (0_pInt,NofMyHomog), source=0.0_pReal)
     allocate(homogState(homog)%subState0(0_pInt,NofMyHomog), source=0.0_pReal)
     allocate(homogState(homog)%state    (0_pInt,NofMyHomog), source=0.0_pReal)

   endif myHomog
 enddo initializeInstances

end subroutine homogenization_isostress_init


!--------------------------------------------------------------------------------------------------
!> @brief partitions the deformation gradient onto the constituents
!--------------------------------------------------------------------------------------------------
subroutine homogenization_isostress_partitionDeformation(F,avgF,el)
 use prec, only: &
   pReal
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_maxNgrains, &
   homogenization_Ngrains
 
 implicit none
 real(pReal),   dimension (3,3,homogenization_maxNgrains), intent(out) :: F                         !< partioned def grad per grain
 real(pReal),   dimension (3,3),                           intent(in)  :: avgF                      !< my average def grad
 integer(pInt),                                            intent(in)  :: &
   el                                                                                               !< element number
 F=0.0_pReal
 F(1:3,1:3,1:homogenization_Ngrains(mesh_element(3,el)))= &
   spread(avgF,3,homogenization_Ngrains(mesh_element(3,el)))

end subroutine homogenization_isostress_partitionDeformation


!--------------------------------------------------------------------------------------------------
!> @brief derive average stress and stiffness from constituent quantities 
!--------------------------------------------------------------------------------------------------
subroutine homogenization_isostress_averageStressAndItsTangent(avgP,dAvgPdAvgF,P,dPdF,el)
 use prec, only: &
   pReal
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_maxNgrains, &
   homogenization_Ngrains, &
   homogenization_typeInstance
 
 implicit none
 real(pReal),   dimension (3,3),                               intent(out) :: avgP                  !< average stress at material point
 real(pReal),   dimension (3,3,3,3),                           intent(out) :: dAvgPdAvgF            !< average stiffness at material point
 real(pReal),   dimension (3,3,homogenization_maxNgrains),     intent(in)  :: P                     !< array of current grain stresses
 real(pReal),   dimension (3,3,3,3,homogenization_maxNgrains), intent(in)  :: dPdF                  !< array of current grain stiffnesses
 integer(pInt),                                                intent(in)  :: el                    !< element number
 integer(pInt) :: &
   homID, & 
   Ngrains

 homID = homogenization_typeInstance(mesh_element(3,el))
 Ngrains = homogenization_Ngrains(mesh_element(3,el))
 avgP       = sum(P,3)   /real(Ngrains,pReal)
 dAvgPdAvgF = sum(dPdF,5)/real(Ngrains,pReal)

end subroutine homogenization_isostress_averageStressAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief return array of homogenization results for post file inclusion 
!--------------------------------------------------------------------------------------------------
pure function homogenization_isostress_postResults(ip,el,avgP,avgF)
 use prec, only: &
   pReal
 use mesh, only: &
   mesh_element, &
   mesh_ipCoordinates
 use material, only: &
   homogenization_typeInstance, &
   homogenization_Noutput
 
 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(3,3), intent(in) :: &
   avgP, &                                                                                          !< average stress at material point
   avgF                                                                                             !< average deformation gradient at material point
 real(pReal),  dimension(homogenization_isostress_sizePostResults &
                         (homogenization_typeInstance(mesh_element(3,el)))) :: &
   homogenization_isostress_postResults
 
 integer(pInt) :: &
   homID, &
   o, c
   
 c = 0_pInt
 homID = homogenization_typeInstance(mesh_element(3,el))
 homogenization_isostress_postResults = 0.0_pReal
 
 do o = 1_pInt,homogenization_Noutput(mesh_element(3,el))
   select case(homogenization_isostress_outputID(o,homID))
     case (nconstituents_ID)
       homogenization_isostress_postResults(c+1_pInt) = real(homogenization_isostress_Ngrains(homID),pReal)
       c = c + 1_pInt
     case (avgdefgrad_ID)
       homogenization_isostress_postResults(c+1_pInt:c+9_pInt) = reshape(avgF,[9])
       c = c + 9_pInt
     case (avgfirstpiola_ID)
       homogenization_isostress_postResults(c+1_pInt:c+9_pInt) = reshape(avgP,[9])
       c = c + 9_pInt
     case (ipcoords_ID)
       homogenization_isostress_postResults(c+1_pInt:c+3_pInt) = mesh_ipCoordinates(1:3,ip,el)                       ! current ip coordinates
       c = c + 3_pInt
    end select
 enddo

end function homogenization_isostress_postResults

end module homogenization_isostress
