!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Mixture homogenization models
!--------------------------------------------------------------------------------------------------
module solute_isoconc
 use prec, only: &
   pInt, &
   pReal
 
 implicit none
 private
 integer(pInt),             dimension(:),     allocatable,         public, protected :: &
   solute_isoconc_sizePostResults
 integer(pInt),             dimension(:,:),   allocatable, target, public :: &
   solute_isoconc_sizePostResult
 
 character(len=64),         dimension(:,:),   allocatable, target, public :: &
  solute_isoconc_output                                                                                !< name of each post result output
 integer(pInt),             dimension(:),     allocatable, target, public :: &
   solute_isoconc_Noutput                                                                              !< number of outputs per homog instance

 enum, bind(c) 
   enumerator :: undefined_ID
 end enum

 type, private :: tParameters                                                                       !< container type for internal constitutive parameters
   integer(kind(undefined_ID)), dimension(:), allocatable :: & 
     outputID
   real(pReal), dimension(:),   allocatable :: &
     initialChemPot
 end type

 type(tParameters), dimension(:), allocatable, private :: param                                     !< containers of constitutive parameters (len Ninstance)

 public :: &
   solute_isoconc_init, &
   solute_isoconc_postResults

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine solute_isoconc_init(fileUnit)
#ifdef __GFORTRAN__
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
 use mesh, only: &
   FE_Nips, &
   FE_geomtype, &
   mesh_NcpElems, &
   mesh_element
 
 implicit none
 integer(pInt),                intent(in) :: fileUnit
 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: &
   section = 0_pInt, i, mySize = 0_pInt, o, el, ip, gr
 integer(pInt) :: sizeState
 integer :: &
   maxNinstance, &
   homog, &
   phase, &
   comp, &
   instance
 integer :: &
   NofmySolute                                                                                       ! no pInt (stores a system dependen value from 'count'
 character(len=65536) :: &
   tag  = '', &
   line = ''
 
 write(6,'(/,a)')   ' <<<+-  solute_'//SOLUTE_isoconc_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = count(solute_type == SOLUTE_isoconc_ID)
 if (maxNinstance == 0) return
 
 if (iand(debug_level(debug_HOMOGENIZATION),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 allocate(solute_isoconc_sizePostResults(maxNinstance),                              source=0_pInt)
 allocate(solute_isoconc_sizePostResult(maxval(homogenization_Noutput),maxNinstance),source=0_pInt)
 allocate(solute_isoconc_Noutput        (maxNinstance),                              source=0_pInt)
 allocate(solute_isoconc_output        (maxval(homogenization_Noutput),maxNinstance))
          solute_isoconc_output = ''
 
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
     if (solute_type(section) == SOLUTE_isoconc_ID) then
       i = solute_typeInstance(section)                                                             ! count instances of my homogenization law
       allocate(param(i)%outputID(homogenization_Noutput(section)))                                 ! allocate space for IDs of every requested output
       allocate(param(i)%initialChemPot(homogenization_Ncomponents(section)), source = 0.0_pReal)
     endif
     cycle
   endif
   if (section > 0_pInt ) then                                                                      ! do not short-circuit here (.and. with next if-statement). It's not safe in Fortran
     if (solute_type(section) == SOLUTE_isoconc_ID) then                         ! one of my sections
       i = solute_typeInstance(section)                                                     ! which instance of my type is present homogenization
       chunkPos = IO_stringPos(line)
       tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                            ! extract key
       select case(tag)
         case ('(output)')
           select case(IO_lc(IO_stringValue(line,chunkPos,2_pInt)))

           end select

         case ('initial_chemicalpotential')
           if (chunkPos(1) /= homogenization_Ncomponents(section) + 1_pInt) &
             call IO_error(150_pInt,ext_msg=trim(tag)//' ('//SOLUTE_isoconc_label//')')
           do comp = 1_pInt, homogenization_Ncomponents(section)
             param(i)%initialChemPot(comp) = IO_floatValue(line,chunkPos,1_pInt+comp)
           enddo

       end select
     endif
   endif
 enddo parsingFile

 initializeInstances: do homog = 1_pInt, material_Nhomogenization
   mySolute: if (solute_type(homog) == SOLUTE_isoconc_ID) then
     NofmySolute = count(material_homog == homog)
     instance = solute_typeInstance(homog)

! *  Determine size of postResults array
     outputsLoop: do o = 1_pInt, solute_isoconc_Noutput(instance)
       select case(param(instance)%outputID(o))

       end select

       outputFound: if (mySize > 0_pInt) then
        solute_isoconc_sizePostResult(o,instance) = mySize
        solute_isoconc_sizePostResults(instance) = &
          solute_isoconc_sizePostResults(instance) + mySize
       endif outputFound
     enddo outputsLoop

! allocate state arrays
     sizeState = 0_pInt
     soluteState(homog)%sizeState = sizeState
     soluteState(homog)%sizePostResults = solute_isoconc_sizePostResults(instance)
     allocate(soluteState(homog)%state0   (sizeState,NofmySolute))
     allocate(soluteState(homog)%subState0(sizeState,NofmySolute))
     allocate(soluteState(homog)%state    (sizeState,NofmySolute))

   endif mySolute
 enddo initializeInstances
 
! state init
 elementLoop: do el = 1_pInt, mesh_NcpElems
   IpLoop: do ip = 1_pInt, FE_Nips(FE_geomtype(mesh_element(2,el)))
     homog = material_homog(ip,el)
     instance = solute_typeInstance(homog)
     mySolute2: if (solute_type(homog) == SOLUTE_isoconc_ID) then
       do gr = 1_pInt, homogenization_Ngrains(homog)
         phase = material_phase(gr,ip,el)
         if (homogenization_Ncomponents(homog) /= phase_Ncomponents(phase)) &
           call IO_error(211_pInt,el=el,ip=ip,g=gr,ext_msg='mismatch in # of homog and phase components')
       enddo
     endif mySolute2
   enddo IPLoop
 enddo elementLoop    

end subroutine solute_isoconc_init


!--------------------------------------------------------------------------------------------------
!> @brief return array of homogenization results for post file inclusion 
!--------------------------------------------------------------------------------------------------
pure function solute_isoconc_postResults(ip,el)
 use material, only: &
   material_homog, &
   solute_typeInstance, &
   homogenization_Noutput
 
 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),  dimension(solute_isoconc_sizePostResults(solute_typeInstance(material_homog(ip,el)))) :: &
   solute_isoconc_postResults
 
 integer(pInt) :: &
   homID, &
   o, c
   
 c = 0_pInt
 homID = solute_typeInstance(material_homog(ip,el))
 solute_isoconc_postResults = 0.0_pReal
 
 do o = 1_pInt,homogenization_Noutput(material_homog(ip,el))
   select case(param(homID)%outputID(o))

   end select
 enddo

end function solute_isoconc_postResults

end module solute_isoconc
