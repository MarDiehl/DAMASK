!--------------------------------------------------------------------------------------------------
!> @author Kishan Govind, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Mixture homogenization models
!--------------------------------------------------------------------------------------------------
module homogenization_multiphase
 use prec, only: &
   pInt, &
   pReal
 
 implicit none
 private
 integer(pInt),               dimension(:),   allocatable,         public, protected :: &
   homogenization_multiphase_sizePostResults
 integer(pInt),               dimension(:,:), allocatable, target, public :: &
   homogenization_multiphase_sizePostResult
 
 character(len=64),           dimension(:,:), allocatable, target, public :: &
  homogenization_multiphase_output                                                                   !< name of each post result output
 integer(pInt),               dimension(:),   allocatable, target, public :: &
   homogenization_multiphase_Noutput                                                                 !< number of outputs per homog instance
 real(pReal),                 dimension(:),   allocatable,        private :: &
   homogenization_multiphase_absTol, &
   homogenization_multiphase_relTol

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 nconstituents_ID, &
                 ipcoords_ID, &
                 avgdefgrad_ID, &
                 avgfirstpiola_ID
 end enum

 integer(kind(undefined_ID)), dimension(:,:), allocatable,         private :: &
  homogenization_multiphase_outputID                                                                 !< ID of each post result output

 enum, bind(c) 
   enumerator :: isostrain_ID, &
                 isofPKstress_ID, &
                 isoCauchystress_ID, &
                 rankone_ID, &
                 laminate_ID
 end enum

 integer(kind(isostrain_ID)), dimension(:),   allocatable,         private :: &
  homogenization_multiphase_mixtureID                                                                !< ID of mixture rule

 type, private :: tMultiphaseState                                                                     !< internal state aliases
   real(pReal), pointer, dimension(:,:) :: &                                                        ! scalars along NipcMyInstance
     newIter, &
     oldIter, &
     searchDir
   real(pReal), pointer, dimension(:)   :: &                                                        ! scalars along NipcMyInstance
     residual, &
     stepLength
 end type

 type(tMultiphaseState),      dimension(:),   allocatable,         private :: &                                       !< state aliases per instance
   state

 public :: &
   homogenization_multiphase_init, &
   homogenization_multiphase_partitionDeformation, &
   homogenization_multiphase_averageStressAndItsTangent, &
   homogenization_multiphase_updateState, &
   homogenization_multiphase_putPhaseFrac, &
   homogenization_multiphase_postResults

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine homogenization_multiphase_init(fileUnit)
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
 use math, only: &
   math_I3, &
   math_inv33, &
   math_mul33x33, &
   math_transpose33, &
   math_EulerToR, &
   math_Mandel6to33, &
   math_Mandel33to6
 use lattice, only: &
   lattice_initialPlasticStrain
 use material
 use mesh, only: &
   FE_Nips, &
   FE_geomtype, &
   mesh_NcpElems, &
   mesh_element
 use crystallite, only: &
   crystallite_F0, &
   crystallite_Fe, &
   crystallite_Fp0, &
   crystallite_Fi0  
 
 implicit none
 integer(pInt),                intent(in) :: fileUnit
 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: &
   section = 0_pInt, i, mySize = 0_pInt, o, el, ip, gr
 integer :: &
   maxNinstance, &
   homog, &
   micro, &
   instance, &
   offset, &
   xioffset
 integer :: &
   NofMyHomog                                                                                       ! no pInt (stores a system dependen value from 'count'
 character(len=65536) :: &
   tag  = '', &
   line = ''
 real(pReal) :: &
   Ldt(3,3,homogenization_maxNgrains), &
   avgLdt(3,3)  
 
 write(6,'(/,a)')   ' <<<+-  homogenization_'//HOMOGENIZATION_MULTIPHASE_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = count(homogenization_type == HOMOGENIZATION_MULTIPHASE_ID)
 if (maxNinstance == 0) return
 
 if (iand(debug_level(debug_HOMOGENIZATION),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 allocate(homogenization_multiphase_sizePostResults(maxNinstance),          source=0_pInt)
 allocate(homogenization_multiphase_sizePostResult(maxval(homogenization_Noutput),maxNinstance), &
                                                                            source=0_pInt)
 allocate(homogenization_multiphase_Noutput(maxNinstance),                  source=0_pInt)
 allocate(homogenization_multiphase_output(maxval(homogenization_Noutput),maxNinstance))
          homogenization_multiphase_output = ''
 allocate(homogenization_multiphase_outputID(maxval(homogenization_Noutput),maxNinstance), &
                                                                           source=undefined_ID)
 allocate(homogenization_multiphase_mixtureID(maxNinstance),               source=isostrain_ID)
 allocate(homogenization_multiphase_absTol(maxNinstance),                  source=0.0_pReal)
 allocate(homogenization_multiphase_relTol(maxNinstance),                  source=0.0_pReal)

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
     if (homogenization_type(section) == HOMOGENIZATION_MULTIPHASE_ID) then                          ! one of my sections
       i = homogenization_typeInstance(section)                                                     ! which instance of my type is present homogenization
       chunkPos = IO_stringPos(line)
       tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                           ! extract key
       select case(tag)
         case ('(output)')
           select case(IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
             case('nconstituents','ngrains')
               homogenization_multiphase_Noutput(i) = homogenization_multiphase_Noutput(i) + 1_pInt
               homogenization_multiphase_outputID(homogenization_multiphase_Noutput(i),i) = nconstituents_ID
               homogenization_multiphase_output(homogenization_multiphase_Noutput(i),i) = &
                 IO_lc(IO_stringValue(line,chunkPos,2_pInt))
             case('ipcoords')
               homogenization_multiphase_Noutput(i) = homogenization_multiphase_Noutput(i) + 1_pInt
               homogenization_multiphase_outputID(homogenization_multiphase_Noutput(i),i) = ipcoords_ID
               homogenization_multiphase_output(homogenization_multiphase_Noutput(i),i) = &
                 IO_lc(IO_stringValue(line,chunkPos,2_pInt))
             case('avgdefgrad','avgf')
               homogenization_multiphase_Noutput(i) = homogenization_multiphase_Noutput(i) + 1_pInt
               homogenization_multiphase_outputID(homogenization_multiphase_Noutput(i),i) = avgdefgrad_ID
               homogenization_multiphase_output(homogenization_multiphase_Noutput(i),i) = &
                 IO_lc(IO_stringValue(line,chunkPos,2_pInt))
             case('avgp','avgfirstpiola','avg1stpiola')
               homogenization_multiphase_Noutput(i) = homogenization_multiphase_Noutput(i) + 1_pInt
               homogenization_multiphase_outputID(homogenization_multiphase_Noutput(i),i) = avgfirstpiola_ID
               homogenization_multiphase_output(homogenization_multiphase_Noutput(i),i) = &
                 IO_lc(IO_stringValue(line,chunkPos,2_pInt))

           end select

         case ('mixture_rule')
           select case(IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
             case('isostrain')
               homogenization_multiphase_mixtureID(i) = isostrain_ID
             case('isofpkstress')
               homogenization_multiphase_mixtureID(i) = isofPKstress_ID
             case('isocauchystress')
               homogenization_multiphase_mixtureID(i) = isoCauchystress_ID
             case('rankone')
               homogenization_multiphase_mixtureID(i) = rankone_ID
             case('laminate')
               homogenization_multiphase_mixtureID(i) = laminate_ID
           end select

         case ('abs_tol')
           homogenization_multiphase_absTol(i) = IO_floatValue(line,chunkPos,2_pInt)

         case ('rel_tol')
           homogenization_multiphase_relTol(i) = IO_floatValue(line,chunkPos,2_pInt)

       end select
     endif
   endif
 enddo parsingFile

 allocate(state(maxNinstance))                                                                      ! internal state aliases
 initializeInstances: do homog = 1_pInt, material_Nhomogenization
   myHomog1: if (homogenization_type(homog) == HOMOGENIZATION_MULTIPHASE_ID) then
     NofMyHomog = count(material_homog == homog)
     instance = homogenization_typeInstance(homog)

! *  Determine size of postResults array
     outputsLoop: do o = 1_pInt, homogenization_multiphase_Noutput(instance)
       select case(homogenization_multiphase_outputID(o,instance))
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
        homogenization_multiphase_sizePostResult(o,instance) = mySize
        homogenization_multiphase_sizePostResults(instance) = &
          homogenization_multiphase_sizePostResults(instance) + mySize
       endif outputFound
     enddo outputsLoop

! allocate state arrays
     select case(homogenization_multiphase_mixtureID(instance))
       case(isostrain_ID)
         homogState(homog)%sizeState = 0_pInt
         
       case(isofPKstress_ID)
         homogState(homog)%sizeState = (homogenization_Ngrains(homog) - 1_pInt)*27_pInt + 2_pInt
         
       case(isoCauchystress_ID)
         homogState(homog)%sizeState = (homogenization_Ngrains(homog) - 1_pInt)*18_pInt + 2_pInt
         
       case(rankone_ID)
         homogState(homog)%sizeState = (homogenization_Ngrains(homog) - 1_pInt)* 9_pInt + 2_pInt

       case(laminate_ID)
         homogState(homog)%sizeState = (homogenization_Ngrains(homog) - 1_pInt)*18_pInt + 2_pInt

     end select
     homogState(homog)%sizePostResults = homogenization_multiphase_sizePostResults(instance)
     allocate(homogState(homog)%state0   (homogState(homog)%sizeState,NofMyHomog), source=0.0_pReal)
     allocate(homogState(homog)%subState0(homogState(homog)%sizeState,NofMyHomog), source=0.0_pReal)
     allocate(homogState(homog)%state    (homogState(homog)%sizeState,NofMyHomog), source=0.0_pReal)

     select case(homogenization_multiphase_mixtureID(instance))
       case(isostrain_ID)
         homogState(homog)%sizeState = 0_pInt
         
       case(isofPKstress_ID)
         state(instance)%newIter   => homogState(homog)% &
                                        state((homogenization_Ngrains(homog) - 1_pInt)* 0_pInt + 1_pInt: &
                                              (homogenization_Ngrains(homog) - 1_pInt)* 9_pInt,1:NofMyHomog)
         state(instance)%oldIter   => homogState(homog)% &
                                        state((homogenization_Ngrains(homog) - 1_pInt)* 9_pInt + 1_pInt: &
                                              (homogenization_Ngrains(homog) - 1_pInt)*18_pInt,1:NofMyHomog)
         state(instance)%searchDir => homogState(homog)% &
                                        state((homogenization_Ngrains(homog) - 1_pInt)*18_pInt + 1_pInt: &
                                              (homogenization_Ngrains(homog) - 1_pInt)*27_pInt,1:NofMyHomog)
         state(instance)%residual  => homogState(homog)% &
                                        state((homogenization_Ngrains(homog) - 1_pInt)*27_pInt + 1_pInt, &
                                                                                               1:NofMyHomog)
         state(instance)%stepLength=> homogState(homog)% &
                                        state((homogenization_Ngrains(homog) - 1_pInt)*27_pInt + 2_pInt, &
                                                                                               1:NofMyHomog)

       case(isoCauchystress_ID)
         state(instance)%newIter   => homogState(homog)% &
                                        state((homogenization_Ngrains(homog) - 1_pInt)* 0_pInt + 1_pInt: &
                                              (homogenization_Ngrains(homog) - 1_pInt)* 6_pInt,1:NofMyHomog)
         state(instance)%oldIter   => homogState(homog)% &
                                        state((homogenization_Ngrains(homog) - 1_pInt)* 6_pInt + 1_pInt: &
                                              (homogenization_Ngrains(homog) - 1_pInt)*12_pInt,1:NofMyHomog)
         state(instance)%searchDir => homogState(homog)% &
                                        state((homogenization_Ngrains(homog) - 1_pInt)*12_pInt + 1_pInt: &
                                              (homogenization_Ngrains(homog) - 1_pInt)*18_pInt,1:NofMyHomog)
         state(instance)%residual  => homogState(homog)% &
                                        state((homogenization_Ngrains(homog) - 1_pInt)*18_pInt + 1_pInt, &
                                                                                               1:NofMyHomog)
         state(instance)%stepLength=> homogState(homog)% &
                                        state((homogenization_Ngrains(homog) - 1_pInt)*18_pInt + 2_pInt, &
                                                                                               1:NofMyHomog)

       case(rankone_ID)
         state(instance)%newIter   => homogState(homog)% &
                                        state((homogenization_Ngrains(homog) - 1_pInt)* 0_pInt + 1_pInt: &
                                              (homogenization_Ngrains(homog) - 1_pInt)* 3_pInt,1:NofMyHomog)
         state(instance)%oldIter   => homogState(homog)% &
                                        state((homogenization_Ngrains(homog) - 1_pInt)* 3_pInt + 1_pInt: &
                                              (homogenization_Ngrains(homog) - 1_pInt)* 6_pInt,1:NofMyHomog)
         state(instance)%searchDir => homogState(homog)% &
                                        state((homogenization_Ngrains(homog) - 1_pInt)* 6_pInt + 1_pInt: &
                                              (homogenization_Ngrains(homog) - 1_pInt)* 9_pInt,1:NofMyHomog)
         state(instance)%residual  => homogState(homog)% &
                                        state((homogenization_Ngrains(homog) - 1_pInt)* 9_pInt + 1_pInt, &
                                                                                               1:NofMyHomog)
         state(instance)%stepLength=> homogState(homog)% &
                                        state((homogenization_Ngrains(homog) - 1_pInt)* 9_pInt + 2_pInt, &
                                                                                               1:NofMyHomog)

       case(laminate_ID)
         state(instance)%newIter   => homogState(homog)% &
                                        state((homogenization_Ngrains(homog) - 1_pInt)* 0_pInt + 1_pInt: &
                                              (homogenization_Ngrains(homog) - 1_pInt)* 6_pInt,1:NofMyHomog)
         state(instance)%oldIter   => homogState(homog)% &
                                        state((homogenization_Ngrains(homog) - 1_pInt)* 6_pInt + 1_pInt: &
                                              (homogenization_Ngrains(homog) - 1_pInt)*12_pInt,1:NofMyHomog)
         state(instance)%searchDir => homogState(homog)% &
                                        state((homogenization_Ngrains(homog) - 1_pInt)*12_pInt + 1_pInt: &
                                              (homogenization_Ngrains(homog) - 1_pInt)*18_pInt,1:NofMyHomog)
         state(instance)%residual  => homogState(homog)% &
                                        state((homogenization_Ngrains(homog) - 1_pInt)*18_pInt + 1_pInt, &
                                                                                               1:NofMyHomog)
         state(instance)%stepLength=> homogState(homog)% &
                                        state((homogenization_Ngrains(homog) - 1_pInt)*18_pInt + 2_pInt, &
                                                                                               1:NofMyHomog)


     end select

     nullify(phasefracMapping(homog)%p)
     phasefracMapping(homog)%p => mappingHomogenization(1,:,:)
     do gr = 1_pInt, homogenization_Ngrains(homog)
       deallocate(phasefrac(gr,homog)%p)
       allocate  (phasefrac(gr,homog)%p(NofMyHomog))
     enddo  

   endif myHomog1
 enddo initializeInstances
 
! state init
 elementLoop: do el = 1_pInt, mesh_NcpElems
   IpLoop: do ip = 1_pInt, FE_Nips(FE_geomtype(mesh_element(2,el)))
     homog = mesh_element(3,el)
     micro = mesh_element(4,el)
     instance = homogenization_typeInstance(homog)
     myHomog2: if (homogenization_type(homog) == HOMOGENIZATION_MULTIPHASE_ID) then
       do gr = 1_pInt, homogenization_Ngrains(homog)
         phasefrac(gr,homog)%p(phasefracMapping(homog)%p(ip,el)) = microstructure_fraction(gr,micro)
       enddo    
       instance = homogenization_typeInstance(homog)
       select case(homogenization_multiphase_mixtureID(instance))
         case(isostrain_ID)
           
         case(isofPKstress_ID)
           offset = mappingHomogenization(1,ip,el)
           do gr = 1_pInt, homogenization_Ngrains(homog)-1_pInt
             xioffset = 9_pInt*(gr - 1_pInt)
             state(instance)%newIter(xioffset+1:xioffset+9,offset) = &
               reshape(math_mul33x33(math_transpose33(math_EulerToR(material_EulerAngles(1:3,gr,ip,el))), &
                                     math_mul33x33(lattice_initialPlasticStrain(1:3,1:3,material_phase(gr,ip,el)), &
                                                   math_EulerToR(material_EulerAngles(1:3,gr,ip,el)))), &
                       shape = [9])
             state(instance)%oldIter(xioffset+1:xioffset+9,offset) = &
               reshape(math_mul33x33(math_transpose33(math_EulerToR(material_EulerAngles(1:3,gr,ip,el))), &
                                     math_mul33x33(lattice_initialPlasticStrain(1:3,1:3,material_phase(gr,ip,el)), &
                                                   math_EulerToR(material_EulerAngles(1:3,gr,ip,el)))), &
                       shape = [9])
           enddo    
           state(instance)%residual  (offset) = 0.0_pReal
           state(instance)%stepLength(offset) = 1.0_pReal

           homogState(homog)%state0   (1:homogState(homog)%sizeState,offset) = &
             homogState(homog)%state(1:homogState(homog)%sizeState,offset)
           homogState(homog)%subState0(1:homogState(homog)%sizeState,offset) = &
             homogState(homog)%state(1:homogState(homog)%sizeState,offset)
             
           crystallite_F0(1:3,1:3,homogenization_Ngrains(homog),ip,el) = math_I3
           do gr = 1_pInt, homogenization_Ngrains(homog)-1_pInt
             xioffset = 9_pInt*(gr - 1_pInt)
             crystallite_F0(1:3,1:3,gr,ip,el) = &
               math_I3 + &
               reshape(state(instance)%newIter(xioffset+1:xioffset+9,offset), shape = [3,3])
             crystallite_F0(1:3,1:3,homogenization_Ngrains(homog),ip,el) = &
               crystallite_F0(1:3,1:3,homogenization_Ngrains(homog),ip,el) - &
               phasefrac(gr,homog)%p(offset)* &
               reshape(state(instance)%newIter(xioffset+1:xioffset+9,offset), shape = [3,3])
             do o = 1_pInt, homogenization_Ngrains(homog)-1_pInt
               xioffset = 9_pInt*(o - 1_pInt)
               crystallite_F0(1:3,1:3,gr,ip,el) = &
                 crystallite_F0(1:3,1:3,gr,ip,el) - &
                 phasefrac(o,homog)%p(offset)* &
                 reshape(state(instance)%newIter(xioffset+1:xioffset+9,offset), shape = [3,3])
             enddo
           enddo                 
           do gr = 1_pInt, homogenization_Ngrains(homog)
             crystallite_Fe(1:3,1:3,gr,ip,el)  = &
               math_mul33x33(crystallite_F0(1:3,1:3,gr,ip,el), &
                             math_inv33(math_mul33x33(crystallite_Fi0(1:3,1:3,gr,ip,el), & 
                                                      crystallite_Fp0(1:3,1:3,gr,ip,el))))
           enddo                          
         
         case(isoCauchystress_ID)
           offset = mappingHomogenization(1,ip,el)
           state(instance)%residual  (offset) = 0.0_pReal
           state(instance)%stepLength(offset) = 1.0_pReal

           homogState(homog)%state0   (1:homogState(homog)%sizeState,offset) = &
             homogState(homog)%state(1:homogState(homog)%sizeState,offset)
           homogState(homog)%subState0(1:homogState(homog)%sizeState,offset) = &
             homogState(homog)%state(1:homogState(homog)%sizeState,offset)
             
           avgLdt = 0.0_pReal
           Ldt(1:3,1:3,homogenization_Ngrains(homog)) = avgLdt
           do gr = 1_pInt, homogenization_Ngrains(homog)-1_pInt
             xioffset = 6_pInt*(gr - 1_pInt)
             Ldt(1:3,1:3,gr) = &
               avgLdt +math_I3 - &
                       math_inv33( &
                         math_I3 + &
                         math_mul33x33( &
                           math_transpose33(math_EulerToR(material_EulerAngles(1:3,gr,ip,el))), &
                           math_mul33x33( &
                             lattice_initialPlasticStrain(1:3,1:3,material_phase(gr,ip,el)), &
                             math_EulerToR(material_EulerAngles(1:3,gr,ip,el)) &
                                        ) &
                                      ) &
                                 ) 
             Ldt(1:3,1:3,homogenization_Ngrains(homog)) = &
               Ldt(1:3,1:3,homogenization_Ngrains(homog)) - &
                 phasefrac(gr,homog)%p(offset)* ( &
                       math_I3 - &
                       math_inv33( &
                         math_I3 + &
                         math_mul33x33( &
                           math_transpose33(math_EulerToR(material_EulerAngles(1:3,gr,ip,el))), &
                           math_mul33x33( &
                             lattice_initialPlasticStrain(1:3,1:3,material_phase(gr,ip,el)), &
                             math_EulerToR(material_EulerAngles(1:3,gr,ip,el)) &
                                        ) &
                                      ) &
                                 ))
             do o = 1_pInt, homogenization_Ngrains(homog)-1_pInt
               xioffset = 6_pInt*(o - 1_pInt)
               Ldt(1:3,1:3,gr) = &
                 Ldt(1:3,1:3,gr) - &
                 phasefrac(o,homog)%p(offset)*( &
                       math_I3 - &
                       math_inv33( &
                         math_I3 + &
                         math_mul33x33( &
                           math_transpose33(math_EulerToR(material_EulerAngles(1:3,gr,ip,el))), &
                           math_mul33x33( &
                             lattice_initialPlasticStrain(1:3,1:3,material_phase(gr,ip,el)), &
                             math_EulerToR(material_EulerAngles(1:3,gr,ip,el)) &
                                        ) &
                                      ) &
                                 ))
             enddo
           enddo                 
           do gr = 1_pInt, homogenization_Ngrains(homog)
             crystallite_F0(1:3,1:3,gr,ip,el) = &
               math_inv33(math_I3 - Ldt(1:3,1:3,gr))
             crystallite_Fe(1:3,1:3,gr,ip,el)  = &
               math_mul33x33(crystallite_F0(1:3,1:3,gr,ip,el), &
                             math_inv33(math_mul33x33(crystallite_Fi0(1:3,1:3,gr,ip,el), & 
                                                      crystallite_Fp0(1:3,1:3,gr,ip,el))))
           enddo                          
           
         case(rankone_ID)

         case(laminate_ID)

       end select
     endif myHomog2
   enddo IPLoop
 enddo elementLoop    

end subroutine homogenization_multiphase_init


!--------------------------------------------------------------------------------------------------
!> @brief partitions the deformation gradient onto the constituents
!--------------------------------------------------------------------------------------------------
subroutine homogenization_multiphase_partitionDeformation(F,F0,avgF,avgF0,ip,el)
 use math, only: &
   math_I3, &
   math_inv33, &
   math_mul33x33, &
   math_Mandel6to33
 use mesh, only: &
   mesh_element
 use material, only: &
   phasefrac, &
   phasefracMapping, &
   homogenization_maxNgrains, &
   homogenization_Ngrains, &
   homogenization_typeInstance
 
 implicit none
 real(pReal),   dimension (3,3,homogenization_maxNgrains), intent(out) :: F                         !< partioned def grad per grain
 real(pReal),   dimension (3,3,homogenization_maxNgrains), intent(in)  :: F0                        !< partioned def grad per grain from previous time increment
 real(pReal),   dimension (3,3),                           intent(in)  :: avgF, &                   !< average def grad
                                                                          avgF0                     !< average def grad from previous time increment
 integer(pInt),                                            intent(in)  :: &
   el, ip                                                                                           !< element number
 integer(pInt) :: homog, instance, grI, grJ, offset, xioffset  
 real(pReal),   dimension (3,3,homogenization_maxNgrains) :: Ldt
 real(pReal),   dimension (3,3)                           :: avgLdt
 
 homog = mesh_element(3,el)
 instance = homogenization_typeInstance(homog)
 select case(homogenization_multiphase_mixtureID(instance))
   case(isostrain_ID)
     F(1:3,1:3,1:homogenization_Ngrains(mesh_element(3,el)))= &
       spread(avgF,3,homogenization_Ngrains(mesh_element(3,el)))
   
   case(isofPKstress_ID)
     offset = phasefracMapping(homog)%p(ip,el)
     F(1:3,1:3,homogenization_Ngrains(homog)) = avgF
     do grI = 1_pInt, homogenization_Ngrains(homog)-1_pInt
       xioffset = 9_pInt*(grI - 1_pInt)
       F(1:3,1:3,grI) = &
         avgF + &
         reshape(state(instance)%newIter(xioffset+1:xioffset+9,offset), shape=[3,3])
       F(1:3,1:3,homogenization_Ngrains(homog)) = &
         F(1:3,1:3,homogenization_Ngrains(homog)) - &
         phasefrac(grI,homog)%p(offset)* &
         reshape(state(instance)%newIter(xioffset+1:xioffset+9,offset), shape=[3,3])
       do grJ = 1_pInt, homogenization_Ngrains(homog)-1_pInt
         xioffset = 9_pInt*(grJ - 1_pInt)
         F(1:3,1:3,grI) = &
           F(1:3,1:3,grI) - &
           phasefrac(grJ,homog)%p(offset)* &
           reshape(state(instance)%newIter(xioffset+1:xioffset+9,offset), shape=[3,3])
       enddo
     enddo                 

   case(isoCauchystress_ID)
     offset = phasefracMapping(homog)%p(ip,el)
     avgLdt = math_mul33x33(avgF - avgF0,math_inv33(avgF))
     Ldt(1:3,1:3,homogenization_Ngrains(homog)) = avgLdt
     do grI = 1_pInt, homogenization_Ngrains(homog)-1_pInt
       xioffset = 6_pInt*(grI - 1_pInt)
       Ldt(1:3,1:3,grI) = &
         avgLdt + &
         math_Mandel6to33(state(instance)%newIter(xioffset+1:xioffset+6,offset))
       Ldt(1:3,1:3,homogenization_Ngrains(homog)) = &
         Ldt(1:3,1:3,homogenization_Ngrains(homog)) - &
         phasefrac(grI,homog)%p(offset)* &
         math_Mandel6to33(state(instance)%newIter(xioffset+1:xioffset+6,offset))
       do grJ = 1_pInt, homogenization_Ngrains(homog)-1_pInt
         xioffset = 6_pInt*(grJ - 1_pInt)
         Ldt(1:3,1:3,grI) = &
           Ldt(1:3,1:3,grI) - &
           phasefrac(grJ,homog)%p(offset)* &
           math_Mandel6to33(state(instance)%newIter(xioffset+1:xioffset+6,offset))                 
       enddo
     enddo                 
     do grI = 1_pInt, homogenization_Ngrains(homog)
       F(1:3,1:3,grI) = &
         math_mul33x33(math_inv33(math_I3 - Ldt(1:3,1:3,grI)),F0(1:3,1:3,grI))
     enddo                          

   case(rankone_ID)

   case(laminate_ID)

 end select

end subroutine homogenization_multiphase_partitionDeformation


!--------------------------------------------------------------------------------------------------
!> @brief derive average stress and stiffness from constituent quantities 
!--------------------------------------------------------------------------------------------------
subroutine homogenization_multiphase_averageStressAndItsTangent(avgP,dAvgPdAvgF, &
                                                                avgF,avgF0,P,dPdF,F,F0,ip,el)
 use IO, only: &
   IO_error
 use math, only: &
   math_I3, &
   math_Mandel3333to66, &
   math_Mandel66to3333, &
   math_tensorcomp3333, &
   math_tensorcomptransp3333, &
   math_tensorproduct3333, &
   math_mul3333xx3333, &
   math_mul33x33, &
   math_inv33, &
   math_det33, &
   math_transpose33, &
   math_invSym3333
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_Ngrains, &
   homogenization_maxNgrains, &
   homogenization_typeInstance, &
   phasefracMapping, &
   phasefrac
 
 implicit none
 real(pReal),   dimension (3,3),                               intent(out) :: avgP                  !< average stress at material point
 real(pReal),   dimension (3,3,3,3),                           intent(out) :: dAvgPdAvgF            !< average stiffness at material point
 real(pReal),   dimension (3,3),                               intent(in)  :: avgF, avgF0                  !< average stress at material point
 real(pReal),   dimension (3,3,homogenization_maxNgrains),     intent(in)  :: P, F, F0                     !< array of current grain stresses
 real(pReal),   dimension (3,3,3,3,homogenization_maxNgrains), intent(in)  :: dPdF                  !< array of current grain stiffnesses
 integer(pInt),                                                intent(in)  :: ip, el                !< element number
 integer(pInt) :: &
   gr, &
   homog, & 
   offset, &
   instance
 real(pReal),   dimension (3,3,3,3) :: &
   dAvgFdP, &
   dCauchydL, &
   dLdCauchy, &
   dAvgLdCauchy, &
   dCauchydAvgL
 real(pReal),   dimension (3,3)     :: &
   cauchy
 
 homog = mesh_element(3,el)
 offset = phasefracMapping(homog)%p(ip,el)
 instance = homogenization_typeInstance(homog)
 select case(homogenization_multiphase_mixtureID(instance))
   case(isostrain_ID)
     avgP = 0.0_pReal
     dAvgPdAvgF = 0.0_pReal
     do gr = 1_pInt, homogenization_Ngrains(homog)
       avgP = avgP + &
          phasefrac(gr,homog)%p(offset)*P(1:3,1:3,gr) 
       dAvgPdAvgF = dAvgPdAvgF + &
            phasefrac(gr,homog)%p(offset)*dPdF(1:3,1:3,1:3,1:3,gr) 
     enddo   

   case(isofPKstress_ID)                
     avgP = P(1:3,1:3,homogenization_Ngrains(homog))
     dAvgFdP = 0.0_pReal
     do gr = 1_pInt, homogenization_Ngrains(homog)
       dAvgFdP = dAvgFdP + phasefrac(gr,homog)%p(offset)* &
                           math_invSym3333(dPdF(1:3,1:3,1:3,1:3,gr))
     enddo 
     dAvgPdAvgF = math_invSym3333(dAvgFdP)

   case(isoCauchystress_ID)                        
     Cauchy = math_mul33x33(P(1:3,1:3,homogenization_Ngrains(homog)), &
                            math_transpose33(F(1:3,1:3,homogenization_Ngrains(homog))))/ &
                            math_det33(F(1:3,1:3,homogenization_Ngrains(homog)))              
     avgP = math_mul33x33(Cauchy,math_inv33(math_transpose33(avgF)))*math_det33(avgF)              
     dAvgLdCauchy = 0.0_pReal
     do gr = 1_pInt, homogenization_Ngrains(homog)
       dCauchydL = &
           math_mul3333xx3333( &
             math_mul3333xx3333( &
               math_tensorcomp3333(math_I3,math_transpose33(F(1:3,1:3,gr))/ &
                                           math_det33(F(1:3,1:3,gr))), &
               dPdF(1:3,1:3,1:3,1:3,gr) &
                               ) + &
             math_tensorcomptransp3333(P(1:3,1:3,gr),math_I3/math_det33(F(1:3,1:3,gr))) - &
             math_mul3333xx3333( &
               math_tensorproduct3333(Cauchy,math_transpose33(F(1:3,1:3,gr))), & 
               math_tensorcomp3333(math_inv33(F(1:3,1:3,gr)),math_inv33(F(1:3,1:3,gr))) &
                               ), &
             math_mul3333xx3333( &
               math_tensorcomp3333(math_I3,F0(1:3,1:3,gr)), &
               math_tensorcomp3333( &
                                   math_mul33x33(F(1:3,1:3,gr),math_inv33(F0(1:3,1:3,gr))), &
                                   math_mul33x33(F(1:3,1:3,gr),math_inv33(F0(1:3,1:3,gr)))  &
                                  ) &
                                ) &
                              )
       dLdCauchy = math_invSym3333(dCauchydL)
       dAvgLdCauchy = dAvgLdCauchy + &
          phasefrac(gr,homog)%p(offset)*dLdCauchy 
     enddo   
     dCauchydAvgL = math_invSym3333(dAvgLdCauchy)
     dAvgPdAvgF = &
       math_mul3333xx3333( &
         math_tensorcomp3333(math_I3,math_inv33(math_transpose33(avgF))*math_det33(avgF)), &
         math_mul3333xx3333( &
                            dCauchydAvgL, &
                            math_tensorcomp3333(math_I3,math_inv33(avgF0)) &
                           ) &
                       ) - &
       math_mul3333xx3333( &
         math_tensorcomptransp3333(Cauchy,math_I3*math_det33(avgF)), &
         math_tensorcomp3333(math_inv33(avgF),math_inv33(avgF)) &
                         ) + &
       math_tensorproduct3333(avgP,math_inv33(math_transpose33(avgF)))

   case(rankone_ID)
     avgP = 0.0_pReal
     dAvgPdAvgF = 0.0_pReal
     do gr = 1_pInt, homogenization_Ngrains(homog)
       avgP = avgP + &
          phasefrac(gr,homog)%p(offset)*P(1:3,1:3,gr) 
       dAvgPdAvgF = dAvgPdAvgF + &
            phasefrac(gr,homog)%p(offset)*dPdF(1:3,1:3,1:3,1:3,gr) 
     enddo   

   case(laminate_ID)
     avgP = 0.0_pReal
     dAvgPdAvgF = 0.0_pReal
     do gr = 1_pInt, homogenization_Ngrains(homog)
       avgP = avgP + &
          phasefrac(gr,homog)%p(offset)*P(1:3,1:3,gr) 
       dAvgPdAvgF = dAvgPdAvgF + &
            phasefrac(gr,homog)%p(offset)*dPdF(1:3,1:3,1:3,1:3,gr) 
     enddo   

 end select

end subroutine homogenization_multiphase_averageStressAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief state update for different mixture rules 
!--------------------------------------------------------------------------------------------------
function homogenization_multiphase_updateState(P,dPdF,F,F0,iter,ip,el)
 use IO, only: &
   IO_error
 use math, only: &
   math_I3, &
   math_Mandel33to6, &
   math_Mandel6to33, &
   math_Mandel3333to66, &
   math_tensorcomp3333, &
   math_tensorcomptransp3333, &
   math_tensorproduct3333, &
   math_mul3333xx3333, &
   math_mul33x33, &
   math_inv33, &
   math_det33, &
   math_transpose33
 use mesh, only: &
   mesh_element
 use material, only: &
   phasefrac, &
   phasefracMapping, &
   homogenization_Ngrains, &
   homogenization_maxNgrains, &
   homogenization_typeInstance
 
 implicit none
 logical,       dimension(2)                                             :: &
   homogenization_multiphase_updateState                                                            !< average stress at material point
 real(pReal),   dimension(3,3,homogenization_maxNgrains),     intent(in) :: &
   P, F, F0                                                                                         !< array of current grain stresses and deformation gradients
 real(pReal),   dimension(3,3,3,3,homogenization_maxNgrains), intent(in) :: dPdF                    !< array of current grain stiffnesses
 integer(pInt),                                               intent(in) :: &
   iter, ip, el                                                                                     !< state loop iter, ip, element number
 real(pReal),   dimension(:),                 allocatable                :: &
   residual, &
   detFiInv   
 real(pReal),   dimension(:,:),               allocatable                :: &
   jacobian, &
   cauchy 
 integer(pInt), dimension(:),                 allocatable                :: &
   ipiv  
 integer(pInt) :: &
   homog, & 
   instance, &
   offset, &
   xioffsetI, &
   xioffsetJ, &
   xisize, &
   grI, grJ, &
   ierr
 real(pReal) :: &
   stressTol, &
   dSi_dLi(6,6), &
   dSn_dLn(6,6)

 external :: &
   dgesv

 homog = mesh_element(3,el)
 instance = homogenization_typeInstance(homog)
 offset = phasefracMapping(homog)%p(ip,el)
 myMixRule: select case(homogenization_multiphase_mixtureID(instance))
   case(isostrain_ID) myMixRule
     homogenization_multiphase_updateState = [.true., .true.]

   case(isofPKstress_ID) myMixRule
     xisize = (homogenization_Ngrains(homog) - 1_pInt)*9_pInt
     allocate(residual(xisize       ), source=0.0_pReal)
     allocate(jacobian(xisize,xisize), source=0.0_pReal)
     do grI = 1_pInt, homogenization_Ngrains(homog)-1
       xioffsetI = 9_pInt*(grI - 1_pInt)
       residual(xioffsetI+1_pInt:xioffsetI+9_pInt) = &
        reshape(P(1:3,1:3,grI) - P(1:3,1:3,homogenization_Ngrains(homog)),[9])
     enddo     
     stressTol = max(homogenization_multiphase_absTol(instance), &
                     homogenization_multiphase_relTol(instance)* &
                     norm2(sum(P,dim=3)/real(homogenization_Ngrains(homog),pReal)))
     
     convergencefPK: if (norm2(residual) < stressTol) then
       homogenization_multiphase_updateState = [.true., .true.]
     elseif (     iter == 1_pInt &
             .or. norm2(residual) < state(instance)%residual(offset)) then convergencefPK              ! not converged, but improved norm of residuum (always proceed in first iteration)...
       homogenization_multiphase_updateState = [.false., .true.]
       state(instance)%residual  (offset) = norm2(residual)                                         ! ...remember old values and...
       state(instance)%stepLength(offset) = &
        min(1.0_pReal,2.0_pReal*state(instance)%stepLength(offset))                                 ! ...proceed with normal step length (calculate new search direction)
       state(instance)%oldIter(:,offset) = state(instance)%newIter(:,offset)

       do grI = 1_pInt,  homogenization_Ngrains(homog)-1
         xioffsetI = 9_pInt*(grI - 1_pInt)
         jacobian(xioffsetI+1_pInt:xioffsetI+9_pInt,xioffsetI+1_pInt:xioffsetI+9_pInt) = &
           reshape(dPdF(1:3,1:3,1:3,1:3,grI),shape=[9,9])
         do grJ = 1_pInt,  homogenization_Ngrains(homog)-1
           xioffsetJ = 9_pInt*(grJ - 1_pInt)
           jacobian(xioffsetI+1_pInt:xioffsetI+9_pInt,xioffsetJ+1_pInt:xioffsetJ+9_pInt) = &
             jacobian(xioffsetI+1_pInt:xioffsetI+9_pInt,xioffsetJ+1_pInt:xioffsetJ+9_pInt) - &
             phasefrac(grJ,homog)%p(offset)* &
             (reshape(dPdF(1:3,1:3,1:3,1:3,grI)                          ,shape=[9,9]) - &
              reshape(dPdF(1:3,1:3,1:3,1:3,homogenization_Ngrains(homog)),shape=[9,9]))
         enddo    
       enddo    
       allocate(ipiv(xisize))
       call dgesv(xisize,1,jacobian,xisize,ipiv,residual,xisize,ierr)                               !< solve Jacobian * delta state = -residual for delta state
       if (ierr == 0_pInt) then
         state(instance)%searchDir(:,offset) = residual
         state(instance)%newIter  (:,offset) = &
           state(instance)%oldIter  (:,offset) - &
           state(instance)%stepLength(offset)* &
           state(instance)%searchDir(:,offset)
       else
        call IO_error(400_pInt,el=el,ip=ip,ext_msg='homogenization multiphase')
       endif        
     else convergencefPK                                                                               ! not converged and residuum not improved...
       homogenization_multiphase_updateState = [.false., .true.]
       state(instance)%stepLength(offset) = &
         state(instance)%stepLength(offset)/2.0_pReal                                      ! ...try with smaller step length in same direction
         state(instance)%newIter  (:,offset) = &
           state(instance)%oldIter  (:,offset) - &
           state(instance)%stepLength(offset)* &
           state(instance)%searchDir(:,offset)
     endif convergencefPK    
     
   case(isoCauchystress_ID) myMixRule
     xisize = (homogenization_Ngrains(homog) - 1_pInt)*6_pInt
     allocate(residual(xisize       ), source=0.0_pReal)
     allocate(jacobian(xisize,xisize), source=0.0_pReal)
     allocate(cauchy(6,homogenization_Ngrains(homog)), source=0.0_pReal)
     allocate(detFiInv(homogenization_Ngrains(homog)), source=0.0_pReal)
     do grI = 1_pInt, homogenization_Ngrains(homog)
       detFiInv(grI) = 1.0_pReal/math_det33(F(1:3,1:3,grI))
       cauchy(1:6,grI) = &
         math_Mandel33to6(math_mul33x33(P(1:3,1:3,grI),math_transpose33(F(1:3,1:3,grI))))* &
         detFiInv(grI)
     enddo     
     do grI = 1_pInt, homogenization_Ngrains(homog)-1
       xioffsetI = 6_pInt*(grI - 1_pInt)
       residual(xioffsetI+1_pInt:xioffsetI+6_pInt) = &
        cauchy(1:6,grI) - cauchy(1:6,homogenization_Ngrains(homog))
     enddo     
     stressTol = max(homogenization_multiphase_absTol(instance), &
                     homogenization_multiphase_relTol(instance)* &
                     norm2(sum(cauchy,dim=2)/real(homogenization_Ngrains(homog),pReal)))
     
     convergenceCauchy: if (norm2(residual) < stressTol) then
       homogenization_multiphase_updateState = [.true., .true.]
     elseif (     iter == 1_pInt &
             .or. norm2(residual) < state(instance)%residual(offset)) then convergenceCauchy              ! not converged, but improved norm of residuum (always proceed in first iteration)...
       homogenization_multiphase_updateState = [.false., .true.]
       state(instance)%residual  (offset) = norm2(residual)                                         ! ...remember old values and...
       state(instance)%stepLength(offset) = &
        min(1.0_pReal,2.0_pReal*state(instance)%stepLength(offset))                                 ! ...proceed with normal step length (calculate new search direction)
       state(instance)%oldIter(:,offset) = state(instance)%newIter(:,offset)

       dSn_dLn(1:6,1:6) = &
         math_Mandel3333to66( &
           math_mul3333xx3333( &
             math_mul3333xx3333( &
               math_tensorcomp3333( &
                                   math_I3, &
                                   math_transpose33(F(1:3,1:3,homogenization_Ngrains(homog)))* &
                                   detFiInv(homogenization_Ngrains(homog)) &
                                  ), &
                                dPdF(1:3,1:3,1:3,1:3,homogenization_Ngrains(homog)) &
                               ) + &
             math_tensorcomptransp3333( &
                                       P(1:3,1:3,homogenization_Ngrains(homog)), &
                                       math_I3*detFiInv(homogenization_Ngrains(homog)) &
                                      ) - &  
             math_mul3333xx3333( &
               math_tensorproduct3333( &
                                      math_Mandel6to33(cauchy(1:6,homogenization_Ngrains(homog))), &
                                      math_transpose33(F(1:3,1:3,homogenization_Ngrains(homog))) &
                                     ), &
               math_tensorcomp3333( &
                                   math_inv33(F(1:3,1:3,homogenization_Ngrains(homog))), &
                                   math_inv33(F(1:3,1:3,homogenization_Ngrains(homog))) &
                                  ) & 
                               ), &
           math_mul3333xx3333( &
             math_tensorcomp3333(math_I3,F0(1:3,1:3,homogenization_Ngrains(homog))), &
             math_tensorcomp3333( &
                                 math_mul33x33(F(1:3,1:3,homogenization_Ngrains(homog)), &
                                               math_inv33(F0(1:3,1:3,homogenization_Ngrains(homog)))), &
                                 math_mul33x33(F(1:3,1:3,homogenization_Ngrains(homog)), &
                                               math_inv33(F0(1:3,1:3,homogenization_Ngrains(homog))))  &
                                ) &
                              ) &
                             ) &
                            )   
       do grI = 1_pInt,  homogenization_Ngrains(homog)-1
         xioffsetI = 6_pInt*(grI - 1_pInt)
         dSi_dLi(1:6,1:6) = &
           math_Mandel3333to66( &
             math_mul3333xx3333( &
               math_mul3333xx3333( &
                 math_tensorcomp3333( &
                                     math_I3, &
                                     math_transpose33(F(1:3,1:3,grI))* &
                                     detFiInv(grI) &
                                    ), &
                                  dPdF(1:3,1:3,1:3,1:3,grI) &
                                 ) + &
               math_tensorcomptransp3333( &
                                         P(1:3,1:3,grI), &
                                         math_I3*detFiInv(grI) &
                                        ) - &  
               math_mul3333xx3333( &
                 math_tensorproduct3333( &
                                        math_Mandel6to33(cauchy(1:6,grI)), &
                                        math_transpose33(F(1:3,1:3,grI)) &
                                       ), &
                 math_tensorcomp3333( &
                                     math_inv33(F(1:3,1:3,grI)), &
                                     math_inv33(F(1:3,1:3,grI)) &
                                    ) & 
                                 ), &
             math_mul3333xx3333( &
               math_tensorcomp3333(math_I3,F0(1:3,1:3,grI)), &
               math_tensorcomp3333( &
                                   math_mul33x33(F(1:3,1:3,grI), &
                                                 math_inv33(F0(1:3,1:3,grI))), &
                                   math_mul33x33(F(1:3,1:3,grI), &
                                                 math_inv33(F0(1:3,1:3,grI)))  &
                                  ) &
                                ) &
                               ) &
                              )   
         jacobian(xioffsetI+1_pInt:xioffsetI+6_pInt,xioffsetI+1_pInt:xioffsetI+6_pInt) = &
           dSi_dLi
         do grJ = 1_pInt,  homogenization_Ngrains(homog)-1
           xioffsetJ = 6_pInt*(grJ - 1_pInt)
           jacobian(xioffsetI+1_pInt:xioffsetI+6_pInt,xioffsetJ+1_pInt:xioffsetJ+6_pInt) = &
             jacobian(xioffsetI+1_pInt:xioffsetI+6_pInt,xioffsetJ+1_pInt:xioffsetJ+6_pInt) - &
             phasefrac(grJ,homog)%p(offset)*(dSi_dLi - dSn_dLn)
         enddo    
       enddo    
       allocate(ipiv(xisize))
       call dgesv(xisize,1,jacobian,xisize,ipiv,residual,xisize,ierr)                               !< solve Jacobian * delta state = -residual for delta state
       if (ierr == 0_pInt) then
         state(instance)%searchDir(:,offset) = residual
         state(instance)%newIter  (:,offset) = &
           state(instance)%oldIter  (:,offset) - &
           state(instance)%stepLength(offset)* &
           state(instance)%searchDir(:,offset)
       else
        call IO_error(400_pInt,el=el,ip=ip,ext_msg='homogenization multiphase')
       endif        
     else convergenceCauchy                                                                               ! not converged and residuum not improved...
       homogenization_multiphase_updateState = [.false., .true.]
       state(instance)%stepLength(offset) = &
         state(instance)%stepLength(offset)/2.0_pReal                                      ! ...try with smaller step length in same direction
         state(instance)%newIter  (:,offset) = &
           state(instance)%oldIter  (:,offset) - &
           state(instance)%stepLength(offset)* &
           state(instance)%searchDir(:,offset)
     endif convergenceCauchy    
     
   case(rankone_ID) myMixRule

   case(laminate_ID) myMixRule

 end select myMixRule

end function homogenization_multiphase_updateState


!--------------------------------------------------------------------------------------------------
!> @brief set module wide phase fractions 
!--------------------------------------------------------------------------------------------------
subroutine homogenization_multiphase_putPhaseFrac(frac,ip,el)
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_Ngrains, &
   homogenization_typeInstance, &
   phasefracMapping, &
   phasefrac
 
 implicit none
 integer(pInt),                                                intent(in)  :: ip, el                !< element number
 real(pReal),   dimension (homogenization_Ngrains(mesh_element(3,el))),     &
                                                               intent(in)  :: frac             !< array of current phase fractions
 integer(pInt) :: &
   gr, &
   homog, & 
   instance, &
   offset

 homog = mesh_element(3,el)
 instance = homogenization_typeInstance(homog)
 offset = phasefracMapping(homog)%p(ip,el)
 do gr = 1_pInt, homogenization_Ngrains(homog)
   phasefrac(gr,homog)%p(offset) =  frac(gr) 
 enddo

end subroutine homogenization_multiphase_putPhaseFrac


!--------------------------------------------------------------------------------------------------
!> @brief return array of homogenization results for post file inclusion 
!--------------------------------------------------------------------------------------------------
pure function homogenization_multiphase_postResults(ip,el,avgP,avgF)
 use mesh, only: &
   mesh_element, &
   mesh_ipCoordinates
 use material, only: &
   homogenization_Ngrains, &
   homogenization_typeInstance, &
   homogenization_Noutput
 
 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(3,3), intent(in) :: &
   avgP, &                                                                                          !< average stress at material point
   avgF                                                                                             !< average deformation gradient at material point
 real(pReal),  dimension(homogenization_multiphase_sizePostResults &
                         (homogenization_typeInstance(mesh_element(3,el)))) :: &
   homogenization_multiphase_postResults
 
 integer(pInt) :: &
   homID, &
   o, c
   
 c = 0_pInt
 homID = homogenization_typeInstance(mesh_element(3,el))
 homogenization_multiphase_postResults = 0.0_pReal
 
 do o = 1_pInt,homogenization_Noutput(mesh_element(3,el))
   select case(homogenization_multiphase_outputID(o,homID))
     case (nconstituents_ID)
       homogenization_multiphase_postResults(c+1_pInt) = real(homogenization_Ngrains(mesh_element(3,el)),pReal)
       c = c + 1_pInt
     case (avgdefgrad_ID)
       homogenization_multiphase_postResults(c+1_pInt:c+9_pInt) = reshape(avgF,[9])
       c = c + 9_pInt
     case (avgfirstpiola_ID)
       homogenization_multiphase_postResults(c+1_pInt:c+9_pInt) = reshape(avgP,[9])
       c = c + 9_pInt
     case (ipcoords_ID)
       homogenization_multiphase_postResults(c+1_pInt:c+3_pInt) = mesh_ipCoordinates(1:3,ip,el)                       ! current ip coordinates
       c = c + 3_pInt
    end select
 enddo

end function homogenization_multiphase_postResults

end module homogenization_multiphase
