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
 integer(pInt),             dimension(:),     allocatable,         public, protected :: &
   homogenization_multiphase_sizePostResults
 integer(pInt),             dimension(:,:),   allocatable, target, public :: &
   homogenization_multiphase_sizePostResult
 
 character(len=64),         dimension(:,:),   allocatable, target, public :: &
  homogenization_multiphase_output                                                                   !< name of each post result output
 integer(pInt),             dimension(:),     allocatable, target, public :: &
   homogenization_multiphase_Noutput                                                                 !< number of outputs per homog instance

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 nconstituents_ID, &
                 ipcoords_ID, &
                 avgdefgrad_ID, &
                 avgfirstpiola_ID, &
                 phasefrac_ID
 end enum

 enum, bind(c) 
   enumerator :: isostrain_ID, &
                 isofPKstress_ID, &
                 isoCauchystress_ID, &
                 partialrankone_ID, &
                 fullrankone_ID
 end enum

 type, private :: tParameters                                                                         !< container type for internal constitutive parameters
   integer(kind(undefined_ID)), dimension(:), allocatable :: & 
     outputID
   integer(kind(isostrain_ID)) :: &
     mixtureID = isostrain_ID
   real(pReal), dimension(:,:), allocatable :: &
     interfaceMobility, &
     interfaceEnergy
   real(pReal) :: &
     interfaceWidth = 0.0_pReal, &
     absTol= 0.0_pReal, &
     relTol= 0.0_pReal, &
     phaseFracTol = 1e-4_pReal
   real(pReal), pointer, dimension(:,:) :: &                                                        ! scalars along NipcMyInstance
     newIter, &
     oldIter, &
     searchDir, &
     interfaceNormal
   real(pReal), pointer, dimension(:)   :: &                                                        ! scalars along NipcMyInstance
     residual, &
     stepLength
 end type

 type(tParameters), dimension(:), allocatable, private :: param                                       !< containers of constitutive parameters (len Ninstance)

 public :: &
   homogenization_multiphase_init, &
   homogenization_multiphase_partitionDeformation, &
   homogenization_multiphase_averageStressAndItsTangent, &
   homogenization_multiphase_updateState, &
   homogenization_multiphase_getPhaseFlux, &
   homogenization_multiphase_getPhaseFluxTangent, &
   homogenization_multiphase_getPhaseSource, &
   homogenization_multiphase_getPhaseSourceTangent, &
   homogenization_multiphase_putPhaseFrac, &
   homogenization_multiphase_putInterfaceNormals, &
   homogenization_multiphase_getComponentConc, &
   homogenization_multiphase_getComponentConcTangent, &
   homogenization_multiphase_getComponentMobility, &
   homogenization_multiphase_putComponentConc, &
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
   math_mul33x3 , &   
   math_transpose33, &
   math_EulerToR
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
   section = 0_pInt, i, mySize = 0_pInt, o, el, ip, grI, grJ
 integer :: &
   Ninterface, &
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
     if (homogenization_type(section) == HOMOGENIZATION_MULTIPHASE_ID) then
       i = homogenization_typeInstance(section)                                                     ! count instances of my homogenization law
       allocate(param(i)%outputID(homogenization_Noutput(section)))                                 ! allocate space for IDs of every requested output
       allocate(param(i)%interfaceMobility(homogenization_Ngrains(section), &
                                           homogenization_Ngrains(section)), source = 0.0_pReal)
       allocate(param(i)%interfaceEnergy  (homogenization_Ngrains(section), &
                                           homogenization_Ngrains(section)), source = 0.0_pReal)
     endif
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
               homogenization_multiphase_output(homogenization_multiphase_Noutput(i),i) = &
                 IO_lc(IO_stringValue(line,chunkPos,2_pInt))
               param(i)%outputID(homogenization_multiphase_Noutput(i)) = nconstituents_ID
             case('ipcoords')
               homogenization_multiphase_Noutput(i) = homogenization_multiphase_Noutput(i) + 1_pInt
               homogenization_multiphase_output(homogenization_multiphase_Noutput(i),i) = &
                 IO_lc(IO_stringValue(line,chunkPos,2_pInt))
               param(i)%outputID(homogenization_multiphase_Noutput(i)) = ipcoords_ID
             case('avgdefgrad','avgf')
               homogenization_multiphase_Noutput(i) = homogenization_multiphase_Noutput(i) + 1_pInt
               homogenization_multiphase_output(homogenization_multiphase_Noutput(i),i) = &
                 IO_lc(IO_stringValue(line,chunkPos,2_pInt))
               param(i)%outputID(homogenization_multiphase_Noutput(i)) = avgdefgrad_ID
             case('avgp','avgfirstpiola','avg1stpiola')
               homogenization_multiphase_Noutput(i) = homogenization_multiphase_Noutput(i) + 1_pInt
               homogenization_multiphase_output(homogenization_multiphase_Noutput(i),i) = &
                 IO_lc(IO_stringValue(line,chunkPos,2_pInt))
               param(i)%outputID(homogenization_multiphase_Noutput(i)) = avgfirstpiola_ID
             case('phasefrac','phasefraction')
               homogenization_multiphase_Noutput(i) = homogenization_multiphase_Noutput(i) + 1_pInt
               homogenization_multiphase_output(homogenization_multiphase_Noutput(i),i) = &
                 IO_lc(IO_stringValue(line,chunkPos,2_pInt))
               param(i)%outputID(homogenization_multiphase_Noutput(i)) = phasefrac_ID

           end select

         case ('mixture_rule')
           select case(IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
             case('isostrain')
               param(i)%mixtureID = isostrain_ID
             case('isofpkstress')
               param(i)%mixtureID = isofPKstress_ID
             case('isocauchystress')
               param(i)%mixtureID = isoCauchystress_ID
             case('partialrankone')
               param(i)%mixtureID = partialrankone_ID
             case('fullrankone')
               param(i)%mixtureID = fullrankone_ID
             case default
               call IO_error(211_pInt,el=i, &
                             ext_msg='mixture_rule ('//HOMOGENIZATION_multiphase_label//')')
           end select

         case ('interface_mobility')
           Ninterface = homogenization_Ngrains(section)*(homogenization_Ngrains(section) - 1_pInt)/2_pInt
           if (chunkPos(1) /= Ninterface + 1_pInt) &
             call IO_error(150_pInt,ext_msg=trim(tag)//' ('//HOMOGENIZATION_multiphase_label//')')
           o = 0_pInt
           do grI = 1_pInt, homogenization_Ngrains(section)
             do grJ = grI+1_pInt, homogenization_Ngrains(section)
               o = o + 1_pInt
               param(i)%interfaceMobility(grI,grJ) = &
                 IO_floatValue(line,chunkPos,1_pInt+o)
               param(i)%interfaceMobility(grJ,grI) = &
                 IO_floatValue(line,chunkPos,1_pInt+o)
             enddo
           enddo

         case ('interface_energy')
           Ninterface = homogenization_Ngrains(section)*(homogenization_Ngrains(section) - 1_pInt)/2_pInt
           if (chunkPos(1) /= Ninterface + 1_pInt) &
             call IO_error(150_pInt,ext_msg=trim(tag)//' ('//HOMOGENIZATION_multiphase_label//')')
           o = 0_pInt
           do grI = 1_pInt, homogenization_Ngrains(section)
             do grJ = grI+1_pInt, homogenization_Ngrains(section)
               o = o + 1_pInt
               param(i)%interfaceEnergy(grI,grJ) = &
                 IO_floatValue(line,chunkPos,1_pInt+o)
               param(i)%interfaceEnergy(grJ,grI) = &
                 IO_floatValue(line,chunkPos,1_pInt+o)
             enddo
           enddo

         case ('interface_width')
           param(i)%InterfaceWidth = IO_floatValue(line,chunkPos,2_pInt)

         case ('abs_tol')
           param(i)%absTol = IO_floatValue(line,chunkPos,2_pInt)

         case ('rel_tol')
           param(i)%relTol = IO_floatValue(line,chunkPos,2_pInt)

         case ('phasefrac_tol')
           param(i)%phasefracTol = IO_floatValue(line,chunkPos,2_pInt)

       end select
     endif
   endif
 enddo parsingFile

 sanityChecks: do homog = 1_pInt, material_Nhomogenization
   myHomog1: if (homogenization_type(homog) == HOMOGENIZATION_MULTIPHASE_ID) then
     instance = homogenization_typeInstance(homog)
     if (any(param(instance)%InterfaceMobility < 0.0_pReal)) &
       call IO_error(211_pInt,el=instance, &
                     ext_msg='interface_mobility ('//HOMOGENIZATION_multiphase_label//')')
     if (any(param(instance)%InterfaceEnergy < 0.0_pReal)) &
       call IO_error(211_pInt,el=instance, &
                     ext_msg='interface_energy ('//HOMOGENIZATION_multiphase_label//')')
     if (param(instance)%InterfaceWidth < 0.0_pReal) &
       call IO_error(211_pInt,el=instance, &
                     ext_msg='interface_width ('//HOMOGENIZATION_multiphase_label//')')
   endif myHomog1
 enddo sanityChecks

 initializeInstances: do homog = 1_pInt, material_Nhomogenization
   myHomog2: if (homogenization_type(homog) == HOMOGENIZATION_MULTIPHASE_ID) then
     NofMyHomog = count(material_homog == homog)
     instance = homogenization_typeInstance(homog)

! *  Determine size of postResults array
     outputsLoop: do o = 1_pInt, homogenization_multiphase_Noutput(instance)
       select case(param(instance)%outputID(o))
        case(nconstituents_ID)
          mySize = 1_pInt
        case(ipcoords_ID)
          mySize = 3_pInt
        case(avgdefgrad_ID, avgfirstpiola_ID)
          mySize = 9_pInt
        case(phasefrac_ID)
          mySize = homogenization_Ngrains(homog)
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
     select case(param(instance)%mixtureID)
       case(isostrain_ID)
         homogState(homog)%sizeState = 0_pInt
         
       case(isofPKstress_ID)
         homogState(homog)%sizeState = &
           27*homogenization_Ngrains(homog)*(homogenization_Ngrains(homog) - 1_pInt)/2_pInt + 2_pInt
         
       case(isoCauchystress_ID)
         homogState(homog)%sizeState = &
           18*homogenization_Ngrains(homog)*(homogenization_Ngrains(homog) - 1_pInt)/2_pInt + 2_pInt
         
       case(partialrankone_ID)
         homogState(homog)%sizeState = &    
           9 *homogenization_Ngrains(homog)*(homogenization_Ngrains(homog) - 1_pInt)/2_pInt + 2_pInt

       case(fullrankone_ID)
         homogState(homog)%sizeState = 0_pInt

     end select
     homogState(homog)%sizePostResults = homogenization_multiphase_sizePostResults(instance)
     allocate(homogState(homog)%state0   (homogState(homog)%sizeState,NofMyHomog), source=0.0_pReal)
     allocate(homogState(homog)%subState0(homogState(homog)%sizeState,NofMyHomog), source=0.0_pReal)
     allocate(homogState(homog)%state    (homogState(homog)%sizeState,NofMyHomog), source=0.0_pReal)

     select case(param(instance)%mixtureID)
       case(isostrain_ID)
         
       case(isofPKstress_ID)
         xioffset = homogenization_Ngrains(homog)*(homogenization_Ngrains(homog) - 1_pInt)/2_pInt
         param(instance)%newIter   => homogState(homog)% &
                                        state(xioffset* 0_pInt + 1_pInt:xioffset* 9_pInt,1:NofMyHomog)
         param(instance)%oldIter   => homogState(homog)% &
                                        state(xioffset* 9_pInt + 1_pInt:xioffset*18_pInt,1:NofMyHomog)
         param(instance)%searchDir => homogState(homog)% &
                                        state(xioffset*18_pInt + 1_pInt:xioffset*27_pInt,1:NofMyHomog)
         param(instance)%residual  => homogState(homog)%state(xioffset*27_pInt + 1_pInt,1:NofMyHomog)
         param(instance)%stepLength=> homogState(homog)%state(xioffset*27_pInt + 2_pInt,1:NofMyHomog)

       case(isoCauchystress_ID)
         xioffset = homogenization_Ngrains(homog)*(homogenization_Ngrains(homog) - 1_pInt)/2_pInt
         param(instance)%newIter   => homogState(homog)% &
                                        state(xioffset* 0_pInt + 1_pInt:xioffset* 6_pInt,1:NofMyHomog)
         param(instance)%oldIter   => homogState(homog)% &
                                        state(xioffset* 6_pInt + 1_pInt:xioffset*12_pInt,1:NofMyHomog)
         param(instance)%searchDir => homogState(homog)% &
                                        state(xioffset*12_pInt + 1_pInt:xioffset*18_pInt,1:NofMyHomog)
         param(instance)%residual  => homogState(homog)%state(xioffset*18_pInt + 1_pInt,1:NofMyHomog)
         param(instance)%stepLength=> homogState(homog)%state(xioffset*18_pInt + 2_pInt,1:NofMyHomog)

       case(partialrankone_ID)
         xioffset = homogenization_Ngrains(homog)*(homogenization_Ngrains(homog) - 1_pInt)/2_pInt
         param(instance)%newIter   => homogState(homog)% &
                                        state(xioffset* 0_pInt + 1_pInt:xioffset* 3_pInt,1:NofMyHomog)
         param(instance)%oldIter   => homogState(homog)% &
                                        state(xioffset* 3_pInt + 1_pInt:xioffset* 6_pInt,1:NofMyHomog)
         param(instance)%searchDir => homogState(homog)% &
                                        state(xioffset* 6_pInt + 1_pInt:xioffset* 9_pInt,1:NofMyHomog)
         allocate(param(instance)%interfaceNormal(3*xioffset,NofMyHomog), source = 0.0_pReal)
         param(instance)%residual  => homogState(homog)%state(xioffset* 9_pInt + 1_pInt,1:NofMyHomog)
         param(instance)%stepLength=> homogState(homog)%state(xioffset* 9_pInt + 2_pInt,1:NofMyHomog)

       case(fullrankone_ID)


     end select

     phasefracMapping(homog)%p => mappingHomogenization(1,:,:)
     do grI = 1_pInt, homogenization_Ngrains(homog)
       allocate  (phasefrac(grI,homog)%p(NofMyHomog))
     enddo  

   endif myHomog2
 enddo initializeInstances
 
! state init
 elementLoop: do el = 1_pInt, mesh_NcpElems
   IpLoop: do ip = 1_pInt, FE_Nips(FE_geomtype(mesh_element(2,el)))
     homog = material_homog(ip,el)
     micro = mesh_element(4,el)
     instance = homogenization_typeInstance(homog)
     myHomog3: if (homogenization_type(homog) == HOMOGENIZATION_MULTIPHASE_ID) then
       do grI = 1_pInt, homogenization_Ngrains(homog)
         phasefrac(grI,homog)%p(phasefracMapping(homog)%p(ip,el)) = microstructure_fraction(grI,micro)
       enddo    
       instance = homogenization_typeInstance(homog)
       select case(param(instance)%mixtureID)
         case(isostrain_ID)
           
         case(isofPKstress_ID)
           offset = mappingHomogenization(1,ip,el)
           param(instance)%residual  (offset) = 0.0_pReal
           param(instance)%stepLength(offset) = 1.0_pReal
         
           do grI = 1_pInt, homogenization_Ngrains(homog)
             crystallite_F0(1:3,1:3,grI,ip,el) = &
               math_I3 + &
               math_mul33x33(math_transpose33(math_EulerToR(material_EulerAngles(1:3,grI,ip,el))), &
                             math_mul33x33(lattice_initialPlasticStrain(1:3,1:3,material_phase(grI,ip,el)), &
                                           math_EulerToR(material_EulerAngles(1:3,grI,ip,el))))
             crystallite_Fe(1:3,1:3,grI,ip,el)  = &
               math_mul33x33(crystallite_F0(1:3,1:3,grI,ip,el), &
                             math_inv33(math_mul33x33(crystallite_Fi0(1:3,1:3,grI,ip,el), & 
                                                      crystallite_Fp0(1:3,1:3,grI,ip,el))))
           enddo                                           
           do grI = 1_pInt, homogenization_Ngrains(homog)-1
             do grJ = grI+1_pInt, homogenization_Ngrains(homog)
               xioffset = 9_pInt*((grI-1_pInt)*(2_pInt*homogenization_Ngrains(homog)-grI)/2_pInt + grJ - grI - 1_pInt)
               param(instance)%newIter(xioffset+1:xioffset+9,offset) = &
                 reshape(crystallite_F0(1:3,1:3,grI,ip,el) - &
                         crystallite_F0(1:3,1:3,grJ,ip,el), shape = [9])
               param(instance)%oldIter(xioffset+1:xioffset+9,offset) = &
                 reshape(crystallite_F0(1:3,1:3,grI,ip,el) - &
                         crystallite_F0(1:3,1:3,grJ,ip,el), shape = [9])
             enddo
           enddo

           homogState(homog)%state0   (1:homogState(homog)%sizeState,offset) = &
             homogState(homog)%state(1:homogState(homog)%sizeState,offset)
           homogState(homog)%subState0(1:homogState(homog)%sizeState,offset) = &
             homogState(homog)%state(1:homogState(homog)%sizeState,offset)

                                           
         case(isoCauchystress_ID)
           offset = mappingHomogenization(1,ip,el)
           param(instance)%residual  (offset) = 0.0_pReal
           param(instance)%stepLength(offset) = 1.0_pReal

           do grI = 1_pInt, homogenization_Ngrains(homog)
             crystallite_F0(1:3,1:3,grI,ip,el) = &
               math_I3 + &
               math_mul33x33(math_transpose33(math_EulerToR(material_EulerAngles(1:3,grI,ip,el))), &
                             math_mul33x33(lattice_initialPlasticStrain(1:3,1:3,material_phase(grI,ip,el)), &
                                           math_EulerToR(material_EulerAngles(1:3,grI,ip,el))))
             crystallite_Fe(1:3,1:3,grI,ip,el)  = &
               math_mul33x33(crystallite_F0(1:3,1:3,grI,ip,el), &
                             math_inv33(math_mul33x33(crystallite_Fi0(1:3,1:3,grI,ip,el), & 
                                                      crystallite_Fp0(1:3,1:3,grI,ip,el))))
           enddo                                           

           homogState(homog)%state0   (1:homogState(homog)%sizeState,offset) = &
             homogState(homog)%state(1:homogState(homog)%sizeState,offset)
           homogState(homog)%subState0(1:homogState(homog)%sizeState,offset) = &
             homogState(homog)%state(1:homogState(homog)%sizeState,offset)


         case(partialrankone_ID) 
           offset = mappingHomogenization(1,ip,el)
           param(instance)%residual  (offset) = 0.0_pReal
           param(instance)%stepLength(offset) = 1.0_pReal
         
           do grI = 1_pInt, homogenization_Ngrains(homog)
             crystallite_F0(1:3,1:3,grI,ip,el) = &
               math_I3 + &
               math_mul33x33(math_transpose33(math_EulerToR(material_EulerAngles(1:3,grI,ip,el))), &
                             math_mul33x33(lattice_initialPlasticStrain(1:3,1:3,material_phase(grI,ip,el)), &
                                           math_EulerToR(material_EulerAngles(1:3,grI,ip,el))))
             crystallite_Fe(1:3,1:3,grI,ip,el)  = &
               math_mul33x33(crystallite_F0(1:3,1:3,grI,ip,el), &
                             math_inv33(math_mul33x33(crystallite_Fi0(1:3,1:3,grI,ip,el), & 
                                                      crystallite_Fp0(1:3,1:3,grI,ip,el))))
           enddo                                           

           do grI = 1_pInt, homogenization_Ngrains(homog)-1 
             do grJ = grI+1_pInt, homogenization_Ngrains(homog)
               xioffset = 3_pInt*((grI-1_pInt)*(2_pInt*homogenization_Ngrains(homog)-grI)/2_pInt + grJ - grI - 1_pInt) 
               param(instance)%newIter(xioffset+1:xioffset+3,offset) = &
                math_mul33x3(crystallite_F0(1:3,1:3,grI,ip,el) - &
                             crystallite_F0(1:3,1:3,grJ,ip,el), &
                             param(instance)%interfaceNormal(xioffset+1_pInt:xioffset+3_pInt,offset)) 
               param(instance)%oldIter(xioffset+1:xioffset+3,offset) = &
                math_mul33x3(crystallite_F0(1:3,1:3,grI,ip,el) - crystallite_F0(1:3,1:3,grJ,ip,el), &
                             param(instance)%interfaceNormal(xioffset+1_pInt:xioffset+3_pInt,offset)) 
             enddo
           enddo
           
           homogState(homog)%state0   (1:homogState(homog)%sizeState,offset) = &
             homogState(homog)%state(1:homogState(homog)%sizeState,offset)
           homogState(homog)%subState0(1:homogState(homog)%sizeState,offset) = &
             homogState(homog)%state(1:homogState(homog)%sizeState,offset)           
           
         case(fullrankone_ID)

       end select
     endif myHomog3
   enddo IPLoop
 enddo elementLoop    

end subroutine homogenization_multiphase_init


!--------------------------------------------------------------------------------------------------
!> @brief partitions the deformation gradient onto the constituents
!--------------------------------------------------------------------------------------------------
subroutine homogenization_multiphase_partitionDeformation(F,F0,avgF,avgF0,ip,el)
 use numerics, only: &
   err_phasefr_tolabs
 use math, only: &
   math_I3, &
   math_inv33, &
   math_mul33x33, &
   math_Mandel6to33, &
   math_tensorproduct33    
 use material, only: &
   material_homog, &
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
 
 homog = material_homog(ip,el)
 instance = homogenization_typeInstance(homog)
 offset = phasefracMapping(homog)%p(ip,el)
 select case(param(instance)%mixtureID)
   case(isostrain_ID)
     do grI = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(grI,homog)%p(offset) > err_phasefr_tolabs) &
         F(1:3,1:3,grI) = avgF
     enddo
   
   case(isofPKstress_ID)
     do grI = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(grI,homog)%p(offset) > err_phasefr_tolabs) &
         F(1:3,1:3,grI) = avgF
     enddo
     
     do grI = 1_pInt, homogenization_Ngrains(homog)-1_pInt
       do grJ = grI+1_pInt, homogenization_Ngrains(homog)
         if (      phasefrac(grI,homog)%p(offset) > err_phasefr_tolabs &
             .and. phasefrac(grJ,homog)%p(offset) > err_phasefr_tolabs) then
           xioffset = 9_pInt*((grI-1_pInt)*(2_pInt*homogenization_Ngrains(homog)-grI)/2_pInt + grJ - grI - 1_pInt)
           F(1:3,1:3,grI) = &
             F(1:3,1:3,grI) + &
             phasefrac(grJ,homog)%p(offset)* &
             reshape(param(instance)%newIter(xioffset+1:xioffset+9,offset), shape = [3,3])
           F(1:3,1:3,grJ) = &
             F(1:3,1:3,grJ) - &
             phasefrac(grI,homog)%p(offset)* &
             reshape(param(instance)%newIter(xioffset+1:xioffset+9,offset), shape = [3,3])
         endif
       enddo
     enddo

   case(isoCauchystress_ID)
     avgLdt = math_mul33x33(avgF - avgF0,math_inv33(avgF))
     do grI = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(grI,homog)%p(offset) > err_phasefr_tolabs) &
         Ldt(1:3,1:3,grI) = avgLdt
     enddo
     
     do grI = 1_pInt, homogenization_Ngrains(homog)-1_pInt
       do grJ = grI+1_pInt, homogenization_Ngrains(homog)
         if (      phasefrac(grI,homog)%p(offset) > err_phasefr_tolabs &
             .and. phasefrac(grJ,homog)%p(offset) > err_phasefr_tolabs) then
           xioffset = 9_pInt*((grI-1_pInt)*(2_pInt*homogenization_Ngrains(homog)-grI)/2_pInt + grJ - grI - 1_pInt)
           Ldt(1:3,1:3,grI) = &
             Ldt(1:3,1:3,grI) + &
             phasefrac(grJ,homog)%p(offset)* &
             math_Mandel6to33(param(instance)%newIter(xioffset+1:xioffset+6,offset))
           Ldt(1:3,1:3,grJ) = &
             Ldt(1:3,1:3,grJ) - &
             phasefrac(grI,homog)%p(offset)* &
             math_Mandel6to33(param(instance)%newIter(xioffset+1:xioffset+6,offset))
         endif
       enddo
     enddo

     do grI = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(grI,homog)%p(offset) > err_phasefr_tolabs) &
         F(1:3,1:3,grI) = math_mul33x33(math_inv33(math_I3 - Ldt(1:3,1:3,grI)),F0(1:3,1:3,grI))
     enddo                          

   case(partialrankone_ID)
     do grI = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(grI,homog)%p(offset) > err_phasefr_tolabs) &
         F(1:3,1:3,grI) = avgF
     enddo
     
     do grI = 1_pInt, homogenization_Ngrains(homog)-1_pInt
       do grJ = grI+1_pInt, homogenization_Ngrains(homog)
         if (      phasefrac(grI,homog)%p(offset) > err_phasefr_tolabs &
             .and. phasefrac(grJ,homog)%p(offset) > err_phasefr_tolabs) then
           xioffset = 3_pInt*((grI-1_pInt)*(2_pInt*homogenization_Ngrains(homog)-grI)/2_pInt + grJ - grI - 1_pInt)
           F(1:3,1:3,grI) = &
             F(1:3,1:3,grI) + &
             phasefrac(grJ,homog)%p(offset)* &
             math_tensorproduct33(param(instance)%newIter(xioffset+1:xioffset+3,offset), &
                                  param(instance)%interfaceNormal(xioffset+1_pInt:xioffset+3_pInt,offset))
           F(1:3,1:3,grJ) = &
             F(1:3,1:3,grJ) - &
             phasefrac(grI,homog)%p(offset)* &
             math_tensorproduct33(param(instance)%newIter(xioffset+1:xioffset+3,offset), &
                                  param(instance)%interfaceNormal(xioffset+1_pInt:xioffset+3_pInt,offset))
         endif
       enddo
     enddo

   case(fullrankone_ID)

 end select

end subroutine homogenization_multiphase_partitionDeformation


!--------------------------------------------------------------------------------------------------
!> @brief derive average stress and stiffness from constituent quantities 
!--------------------------------------------------------------------------------------------------
subroutine homogenization_multiphase_averageStressAndItsTangent(avgP,dAvgPdAvgF, &
                                                                avgF,avgF0,P,dPdF,F,F0,ip,el)
 use IO, only: &
   IO_error
 use numerics, only: &
   err_phasefr_tolabs
 use math, only: &
   math_I3, &
   math_tensorcomp3333, &
   math_tensorcomptransp3333, &
   math_tensorproduct3333, &
   math_mul3333xx3333, &
   math_mul33x33, &
   math_inv33, &
   math_det33, &
   math_transpose33, &
   math_invSym3333
 use material, only: &
   material_homog, &
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
 
 homog = material_homog(ip,el)
 offset = phasefracMapping(homog)%p(ip,el)
 instance = homogenization_typeInstance(homog)
 select case(param(instance)%mixtureID)
   case(isostrain_ID)
     avgP = 0.0_pReal
     dAvgPdAvgF = 0.0_pReal
     do gr = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(gr,homog)%p(offset) > err_phasefr_tolabs) then
         avgP = avgP + &
            phasefrac(gr,homog)%p(offset)*P(1:3,1:3,gr) 
         dAvgPdAvgF = dAvgPdAvgF + &
              phasefrac(gr,homog)%p(offset)*dPdF(1:3,1:3,1:3,1:3,gr) 
       endif
     enddo   

   case(isofPKstress_ID)                
     avgP = 0.0_pReal
     dAvgFdP = 0.0_pReal
     do gr = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(gr,homog)%p(offset) > err_phasefr_tolabs) then
         avgP = P(1:3,1:3,gr)
         dAvgFdP = dAvgFdP + phasefrac(gr,homog)%p(offset)* &
                             math_invSym3333(dPdF(1:3,1:3,1:3,1:3,gr))
       endif
     enddo 
     dAvgPdAvgF = math_invSym3333(dAvgFdP)

   case(isoCauchystress_ID)                        
     dAvgLdCauchy = 0.0_pReal
     do gr = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(gr,homog)%p(offset) > err_phasefr_tolabs) then
         Cauchy = &
           math_mul33x33(P(1:3,1:3,gr),math_transpose33(F(1:3,1:3,gr)))/math_det33(F(1:3,1:3,gr))              
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
       endif
     enddo   
     dCauchydAvgL = math_invSym3333(dAvgLdCauchy)
     avgP = math_mul33x33(Cauchy,math_inv33(math_transpose33(avgF)))*math_det33(avgF)              
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

   case(partialrankone_ID)
     avgP = 0.0_pReal
     dAvgPdAvgF = 0.0_pReal
     do gr = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(gr,homog)%p(offset) > err_phasefr_tolabs) then
         avgP = avgP + &
            phasefrac(gr,homog)%p(offset)*P(1:3,1:3,gr) 
         dAvgPdAvgF = dAvgPdAvgF + &
              phasefrac(gr,homog)%p(offset)*dPdF(1:3,1:3,1:3,1:3,gr) 
       endif
     enddo

   case(fullrankone_ID)

 end select

end subroutine homogenization_multiphase_averageStressAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief state update for different mixture rules 
!--------------------------------------------------------------------------------------------------
function homogenization_multiphase_updateState(P,dPdF,F,F0,iter,ip,el)
 use IO, only: &
   IO_error
 use numerics, only: &
   err_phasefr_tolabs
 use math, only: &
   math_I3, &
   math_Mandel33to6, &
   math_Mandel6to33, &
   math_Mandel3333to66, &
   math_tensorcomp3333, &
   math_tensorcomptransp3333, &
   math_tensorproduct33 , &  
   math_tensorproduct3333, &
   math_mul3333xx3333, &
   math_mul33x33, &
   math_mul33x3 , &      
   math_inv33, &
   math_det33, &
   math_transpose33
 use material, only: &
   material_homog, &
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
 real(pReal),   dimension(3,3,3,3,homogenization_maxNgrains), intent(in) :: &
   dPdF                                                                                             !< array of current grain stiffnesses
 integer(pInt),                                               intent(in) :: &
   iter, ip, el                                                                                     !< state loop iter, ip, element number
 real(pReal),   dimension(3,3)                                           :: &
   avgR                                                                                             !< average stress
 real(pReal),   dimension(:),                 allocatable                :: &
   residual, &
   detFiInv
 real(pReal),   dimension(:,:),               allocatable                :: &
   jacobian, &
   cauchy 
 real(pReal),   dimension(:,:,:),             allocatable                :: &
   dS_dL
 integer(pInt), dimension(:),                 allocatable                :: &
   ipiv, &
   active  
 integer(pInt) :: &
   homog, & 
   instance, &
   offset, &
   xioffsetI, &
   xioffsetJ, &
   xioffsetIN, &
   xioffsetJN, &
   xisize, &
   grI, grJ, grII, grIJ, grJI, grJJ, ii, jj, &
   Nactive, &
   ierr
 real(pReal) :: &
   stressTol

 external :: &
   dgesv

 homog = material_homog(ip,el)
 instance = homogenization_typeInstance(homog)
 offset = phasefracMapping(homog)%p(ip,el)
 myMixRule: select case(param(instance)%mixtureID)
   case(isostrain_ID) myMixRule
     homogenization_multiphase_updateState = [.true., .true.]

   case(isofPKstress_ID) myMixRule
     Nactive = 0_pInt
     do grI = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(grI,homog)%p(offset) > err_phasefr_tolabs) Nactive = Nactive + 1_pInt
     enddo
     allocate(active(Nactive), source = 0_pInt)

     avgR = 0.0_pReal
     Nactive = 0_pInt
     do grI = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(grI,homog)%p(offset) > err_phasefr_tolabs) then
         Nactive = Nactive + 1_pInt
         active(Nactive) = grI
         avgR = avgR + phasefrac(grI,homog)%p(offset)*P(1:3,1:3,grI)
       endif
     enddo     
     
     xisize = 9_pInt*(Nactive-1_pInt)
     allocate(residual(xisize       ), source=0.0_pReal)
     allocate(jacobian(xisize,xisize), source=0.0_pReal)
     do grI = 1_pInt, Nactive-1_pInt
       xioffsetI = 9_pInt*(grI-1_pInt)
       residual(xioffsetI+1_pInt:xioffsetI+9_pInt) = &
         reshape(P(1:3,1:3,active(grI)) - P(1:3,1:3,active(Nactive)),[9])
     enddo 
     
     stressTol = max(            param(instance)%absTol, &
                     norm2(avgR)*param(instance)%relTol)
     convergencefPK: if (     norm2(residual) < stressTol &
                         .or. Nactive         < 2_pInt) then
       homogenization_multiphase_updateState = [.true., .true.]
     elseif (     iter == 1_pInt &
             .or. norm2(residual) < param(instance)%residual(offset)) then convergencefPK           ! not converged, but improved norm of residuum (always proceed in first iteration)...
       homogenization_multiphase_updateState = [.false., .true.]
       param(instance)%residual  (offset) = norm2(residual)                                         ! ...remember old values and...
       param(instance)%stepLength(offset) = &
        min(1.0_pReal,2.0_pReal*param(instance)%stepLength(offset))                                 ! ...proceed with normal step length (calculate new search direction)
       param(instance)%oldIter(:,offset) = param(instance)%newIter(:,offset)

       do grI = 1_pInt,  Nactive-1
         xioffsetI = 9_pInt*(grI - 1_pInt)
         jacobian(xioffsetI+1_pInt:xioffsetI+9_pInt,xioffsetI+1_pInt:xioffsetI+9_pInt) = &
           reshape(dPdF(1:3,1:3,1:3,1:3,active(grI)),shape=[9,9])
         do grJ = 1_pInt,  Nactive-1
           xioffsetJ = 9_pInt*(grJ - 1_pInt)
           jacobian(xioffsetI+1_pInt:xioffsetI+9_pInt,xioffsetJ+1_pInt:xioffsetJ+9_pInt) = &
             jacobian(xioffsetI+1_pInt:xioffsetI+9_pInt,xioffsetJ+1_pInt:xioffsetJ+9_pInt) - &
             phasefrac(active(grJ),homog)%p(offset)* &
             (reshape(dPdF(1:3,1:3,1:3,1:3,active(grI))    ,shape=[9,9]) - &
              reshape(dPdF(1:3,1:3,1:3,1:3,active(Nactive)),shape=[9,9]))
         enddo    
       enddo    
       allocate(ipiv(xisize))
       call dgesv(xisize,1,jacobian,xisize,ipiv,residual,xisize,ierr)                               !< solve Jacobian * delta state = -residual for delta state
       if (ierr == 0_pInt) then
         param(instance)%searchDir(:,offset) = 0.0_pReal
         do grI = 1_pInt, Nactive-1_pInt
           xioffsetIN = 9_pInt*((active(grI)-1_pInt)*(2_pInt*homogenization_Ngrains(homog)-active(grI))/2_pInt + &
                               active(Nactive) - active(grI) - 1_pInt)
           xioffsetI = 9_pInt*(grI - 1_pInt)
           param(instance)%searchDir(xioffsetIN+1_pInt:xioffsetIN+9_pInt,offset) = &
             residual(xioffsetI+1_pInt:xioffsetI+9_pInt)
           do grJ = grI+1_pInt, Nactive-1_pInt
             xioffsetIN = 9_pInt*((active(grI)-1_pInt)*(2_pInt*homogenization_Ngrains(homog)-active(grI))/2_pInt + &
                                 active(grJ) - active(grI) - 1_pInt)
             xioffsetJ = 9_pInt*(grJ - 1_pInt)
             param(instance)%searchDir(xioffsetIN+1_pInt:xioffsetIN+9_pInt,offset) = &
               residual(xioffsetI+1_pInt:xioffsetI+9_pInt) - &
               residual(xioffsetJ+1_pInt:xioffsetJ+9_pInt)               
           enddo
         enddo
         param(instance)%newIter  (:,offset) = &
           param(instance)%oldIter  (:,offset) - &
           param(instance)%stepLength(offset)* &
           param(instance)%searchDir(:,offset)
       else
        call IO_error(400_pInt,el=el,ip=ip,ext_msg='homogenization multiphase')
       endif        
     else convergencefPK                                                                            ! not converged and residuum not improved...
       homogenization_multiphase_updateState = [.false., .true.]
       param(instance)%stepLength(offset) = &
         param(instance)%stepLength(offset)/2.0_pReal                                               ! ...try with smaller step length in same direction
       param(instance)%newIter  (:,offset) = &
         param(instance)%oldIter  (:,offset) - &
         param(instance)%stepLength(offset)* &
         param(instance)%searchDir(:,offset)
     endif convergencefPK    
     
   case(isoCauchystress_ID) myMixRule
     Nactive = 0_pInt
     do grI = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(grI,homog)%p(offset) > err_phasefr_tolabs) Nactive = Nactive + 1_pInt
     enddo
     allocate(active(Nactive), source = 0_pInt)
     allocate(cauchy(6,Nactive), source=0.0_pReal)
     allocate(detFiInv(Nactive), source=0.0_pReal)

     avgR = 0.0_pReal
     Nactive = 0_pInt
     do grI = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(grI,homog)%p(offset) > err_phasefr_tolabs) then
         Nactive = Nactive + 1_pInt
         active(Nactive) = grI
         detFiInv(Nactive) = 1.0_pReal/math_det33(F(1:3,1:3,grI))
         cauchy(1:6,Nactive) = &
           math_Mandel33to6(math_mul33x33(P(1:3,1:3,grI),math_transpose33(F(1:3,1:3,grI))))* &
           detFiInv(Nactive)
         avgR = avgR + phasefrac(grI,homog)%p(offset)*math_Mandel6to33(cauchy(1:6,Nactive))
       endif
     enddo     

     xisize = 6_pInt*(Nactive-1_pInt)
     allocate(residual(xisize       ), source=0.0_pReal)
     allocate(jacobian(xisize,xisize), source=0.0_pReal)
     do grI = 1_pInt, Nactive-1_pInt
       xioffsetI = 6_pInt*(grI-1_pInt)
       residual(xioffsetI+1_pInt:xioffsetI+6_pInt) = &
         cauchy(1:6,grI) - cauchy(1:6,Nactive)
     enddo 
     
     stressTol = max(            param(instance)%absTol, &
                     norm2(avgR)*param(instance)%relTol)
     convergenceCauchy: if (     norm2(residual) < stressTol &
                            .or. Nactive         < 2_pInt) then
       homogenization_multiphase_updateState = [.true., .true.]
     elseif (     iter == 1_pInt &
             .or. norm2(residual) < param(instance)%residual(offset)) then convergenceCauchy        ! not converged, but improved norm of residuum (always proceed in first iteration)...
       homogenization_multiphase_updateState = [.false., .true.]
       param(instance)%residual  (offset) = norm2(residual)                                         ! ...remember old values and...
       param(instance)%stepLength(offset) = &
        min(1.0_pReal,2.0_pReal*param(instance)%stepLength(offset))                                 ! ...proceed with normal step length (calculate new search direction)
       param(instance)%oldIter(:,offset) = param(instance)%newIter(:,offset)

       allocate(dS_dL(6,6,Nactive), source = 0.0_pReal)
       do grI = 1_pInt,  Nactive
         dS_dL(1:6,1:6,grI) = &
           math_Mandel3333to66( &
             math_mul3333xx3333( &
               math_mul3333xx3333( &
                 math_tensorcomp3333( &
                                     math_I3, &
                                     math_transpose33(F(1:3,1:3,active(grI)))* &
                                     detFiInv(grI) &
                                    ), &
                                  dPdF(1:3,1:3,1:3,1:3,active(grI)) &
                                 ) + &
               math_tensorcomptransp3333( &
                                         P(1:3,1:3,active(grI)), &
                                         math_I3*detFiInv(grI) &
                                        ) - &  
               math_mul3333xx3333( &
                 math_tensorproduct3333( &
                                        math_Mandel6to33(cauchy(1:6,grI)), &
                                        math_transpose33(F(1:3,1:3,active(grI))) &
                                       ), &
                 math_tensorcomp3333( &
                                     math_inv33(F(1:3,1:3,active(grI))), &
                                     math_inv33(F(1:3,1:3,active(grI))) &
                                    ) & 
                                 ), &
             math_mul3333xx3333( &
               math_tensorcomp3333(math_I3,F0(1:3,1:3,active(grI))), &
               math_tensorcomp3333( &
                                   math_mul33x33(F(1:3,1:3,active(grI)), &
                                                 math_inv33(F0(1:3,1:3,active(grI)))), &
                                   math_mul33x33(F(1:3,1:3,active(grI)), &
                                                 math_inv33(F0(1:3,1:3,active(grI))))  &
                                  ) &
                                ) &
                               ) &
                              )   
       enddo
       do grI = 1_pInt,  Nactive-1
         xioffsetI = 6_pInt*(grI - 1_pInt)
         jacobian(xioffsetI+1_pInt:xioffsetI+6_pInt,xioffsetI+1_pInt:xioffsetI+6_pInt) = &
           dS_dL(1:6,1:6,grI)
         do grJ = 1_pInt,  Nactive-1
           xioffsetJ = 6_pInt*(grJ - 1_pInt)
           jacobian(xioffsetI+1_pInt:xioffsetI+6_pInt,xioffsetJ+1_pInt:xioffsetJ+6_pInt) = &
             jacobian(xioffsetI+1_pInt:xioffsetI+6_pInt,xioffsetJ+1_pInt:xioffsetJ+6_pInt) - &
             phasefrac(active(grJ),homog)%p(offset)* (dS_dL(1:6,1:6,grI) - dS_dL(1:6,1:6,Nactive))
         enddo    
       enddo    
       allocate(ipiv(xisize))
       call dgesv(xisize,1,jacobian,xisize,ipiv,residual,xisize,ierr)                               !< solve Jacobian * delta state = -residual for delta state
       if (ierr == 0_pInt) then
         param(instance)%searchDir(:,offset) = 0.0_pReal
         do grI = 1_pInt, Nactive-1_pInt
           xioffsetIN = 6_pInt*((active(grI)-1_pInt)*(2_pInt*homogenization_Ngrains(homog)-active(grI))/2_pInt + &
                               active(Nactive) - active(grI) - 1_pInt)
           xioffsetI = 6_pInt*(grI - 1_pInt)
           param(instance)%searchDir(xioffsetIN+1_pInt:xioffsetIN+6_pInt,offset) = &
             residual(xioffsetI+1_pInt:xioffsetI+6_pInt)
           do grJ = grI+1_pInt, Nactive-1_pInt
             xioffsetIN = 6_pInt*((active(grI)-1_pInt)*(2_pInt*homogenization_Ngrains(homog)-active(grI))/2_pInt + &
                                 active(grJ) - active(grI) - 1_pInt)
             xioffsetJ = 6_pInt*(grJ - 1_pInt)
             param(instance)%searchDir(xioffsetIN+1_pInt:xioffsetIN+6_pInt,offset) = &
               residual(xioffsetI+1_pInt:xioffsetI+6_pInt) - &
               residual(xioffsetJ+1_pInt:xioffsetJ+6_pInt)               
           enddo
         enddo
         param(instance)%newIter  (:,offset) = &
           param(instance)%oldIter  (:,offset) - &
           param(instance)%stepLength(offset)* &
           param(instance)%searchDir(:,offset)
       else
        call IO_error(400_pInt,el=el,ip=ip,ext_msg='homogenization multiphase')
       endif        
     else convergenceCauchy                                                                         ! not converged and residuum not improved...
       homogenization_multiphase_updateState = [.false., .true.]
       param(instance)%stepLength(offset) = &
         param(instance)%stepLength(offset)/2.0_pReal                                               ! ...try with smaller step length in same direction
       param(instance)%newIter  (:,offset) = &
         param(instance)%oldIter  (:,offset) - &
         param(instance)%stepLength(offset)* &
         param(instance)%searchDir(:,offset)
     endif convergenceCauchy    
     

   case(partialrankone_ID) myMixRule    
     Nactive = 0_pInt
     do grI = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(grI,homog)%p(offset) > err_phasefr_tolabs) Nactive = Nactive + 1_pInt
     enddo
     allocate(active(Nactive), source = 0_pInt)
    
     avgR = 0.0_pReal
     Nactive = 0_pInt
     do grI = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(grI,homog)%p(offset) > err_phasefr_tolabs) then
         Nactive = Nactive + 1_pInt 
         active(Nactive) = grI  
         avgR = avgR + phasefrac(grI,homog)%p(offset)*P(1:3,1:3,grI)
       endif
     enddo     
     
     xisize = 3_pInt*Nactive*(Nactive-1_pInt)/2_pInt
     allocate(residual(xisize       ), source=0.0_pReal)
     allocate(jacobian(xisize,xisize), source=0.0_pReal)      
     do grI = 1_pInt, Nactive-1_pInt
       do grJ = grI+1_pInt, Nactive
         xioffsetIN = 3_pInt*((active(grI)-1_pInt)* &
                              (2_pInt*homogenization_Ngrains(homog)-active(grI))/2_pInt + &
                              active(grJ) - active(grI) - 1_pInt)
         xioffsetI  = 3_pInt*((grI-1_pInt)*(2_pInt*Nactive-grI)/2_pInt + grJ - grI - 1_pInt)                                      
         residual(xioffsetI+1_pInt:xioffsetI+3_pInt) = &                      
           phasefrac(active(grI),homog)%p(offset)* &
           phasefrac(active(grJ),homog)%p(offset)* &
           math_mul33x3(P(1:3,1:3,active(grI)) - P(1:3,1:3,active(grJ)), &
                        param(instance)%interfaceNormal(xioffsetIN+1_pInt:xioffsetIN+3_pInt,offset))    
       enddo
     enddo 
     
     stressTol = max(            param(instance)%absTol, &
                     norm2(avgR)*param(instance)%relTol)
     convergencefRankone: if (     norm2(residual) < stressTol &
                              .or. Nactive         < 2_pInt) then
       homogenization_multiphase_updateState = [.true., .true.]  
     elseif (     iter == 1_pInt &
             .or. norm2(residual) < param(instance)%residual(offset)) then convergencefRankone      ! not converged, but improved norm of residuum (always proceed in first iteration)...
       homogenization_multiphase_updateState = [.false., .true.]
       
       param(instance)%residual  (offset) = norm2(residual)                                         ! ...remember old values and...
       param(instance)%stepLength(offset) = &
        min(1.0_pReal,2.0_pReal*param(instance)%stepLength(offset))                                 ! ...proceed with normal step length (calculate new search direction)
       param(instance)%oldIter(:,offset) = param(instance)%newIter(:,offset)
       
       do grII = 1_pInt,  Nactive-1; do grIJ = grII+1_pInt, Nactive
         xioffsetIN = 3_pInt*((active(grII)-1_pInt)* &
                              (2_pInt*homogenization_Ngrains(homog)-active(grII))/2_pInt + &
                              active(grIJ) - active(grII) - 1_pInt)
         xioffsetI  = 3_pInt*((grII-1_pInt)*(2_pInt*Nactive-grII)/2_pInt + grIJ - grII - 1_pInt)                                      
         do grJI = 1_pInt,  Nactive-1; do grJJ = grJI+1_pInt, Nactive
           xioffsetJN = 3_pInt*((active(grJI)-1_pInt)* &
                                (2_pInt*homogenization_Ngrains(homog)-active(grJI))/2_pInt + &
                                active(grJJ) - active(grJI) - 1_pInt)
           xioffsetJ  = 3_pInt*((grJI-1_pInt)*(2_pInt*Nactive-grJI)/2_pInt + grJJ - grJI - 1_pInt)                                      
           if (grII == grJI) &
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                 
                 jacobian(xioffsetI+ii,xioffsetJ+jj) = &
                   jacobian(xioffsetI+ii,xioffsetJ+jj) +  &
                   phasefrac(active(grII),homog)%p(offset)* &
                   phasefrac(active(grIJ),homog)%p(offset)* &
                   phasefrac(active(grJJ),homog)%p(offset)* &
                   sum(dPdF(ii,1:3,jj,1:3,active(grII))* &
                       math_tensorproduct33(param(instance)%interfaceNormal(xioffsetIN+1_pInt:xioffsetIN+3_pInt,offset), &
                                            param(instance)%interfaceNormal(xioffsetJN+1_pInt:xioffsetJN+3_pInt,offset)))                            
           if (grII == grJJ) &
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                 
                 jacobian(xioffsetI+ii,xioffsetJ+jj) = &
                   jacobian(xioffsetI+ii,xioffsetJ+jj) -  &
                   phasefrac(active(grII),homog)%p(offset)* &
                   phasefrac(active(grIJ),homog)%p(offset)* &
                   phasefrac(active(grJI),homog)%p(offset)* &
                   sum(dPdF(ii,1:3,jj,1:3,active(grII))* &
                       math_tensorproduct33(param(instance)%interfaceNormal(xioffsetIN+1_pInt:xioffsetIN+3_pInt,offset), &
                                            param(instance)%interfaceNormal(xioffsetJN+1_pInt:xioffsetJN+3_pInt,offset)))                            
           if (grIJ == grJI) &
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                 
                 jacobian(xioffsetI+ii,xioffsetJ+jj) = &
                   jacobian(xioffsetI+ii,xioffsetJ+jj) -  &
                   phasefrac(active(grII),homog)%p(offset)* &
                   phasefrac(active(grIJ),homog)%p(offset)* &
                   phasefrac(active(grJJ),homog)%p(offset)* &
                   sum(dPdF(ii,1:3,jj,1:3,active(grIJ))* &
                       math_tensorproduct33(param(instance)%interfaceNormal(xioffsetIN+1_pInt:xioffsetIN+3_pInt,offset), &
                                            param(instance)%interfaceNormal(xioffsetJN+1_pInt:xioffsetJN+3_pInt,offset)))                            
           if (grIJ == grJJ) &
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                 
                 jacobian(xioffsetI+ii,xioffsetJ+jj) = &
                   jacobian(xioffsetI+ii,xioffsetJ+jj) +  &
                   phasefrac(active(grII),homog)%p(offset)* &
                   phasefrac(active(grIJ),homog)%p(offset)* &
                   phasefrac(active(grJI),homog)%p(offset)* &
                   sum(dPdF(ii,1:3,jj,1:3,active(grIJ))* &
                       math_tensorproduct33(param(instance)%interfaceNormal(xioffsetIN+1_pInt:xioffsetIN+3_pInt,offset), &
                                            param(instance)%interfaceNormal(xioffsetJN+1_pInt:xioffsetJN+3_pInt,offset)))                            
         enddo; enddo
       enddo; enddo  

       allocate(ipiv(xisize))
       call dgesv(xisize,1,jacobian,xisize,ipiv,residual,xisize,ierr)                               !< solve Jacobian * delta state = -residual for delta state
       if (ierr == 0_pInt) then
         param(instance)%searchDir(:,offset) = 0.0_pReal
         do grI = 1_pInt, Nactive-1_pInt
           do grJ = grI+1_pInt, Nactive
             xioffsetIN = 3_pInt*((active(grI)-1_pInt)*(2_pInt*homogenization_Ngrains(homog)-active(grI))/2_pInt + &
                                  active(grJ) - active(grI) - 1_pInt)
             xioffsetI = 3_pInt*((grI-1_pInt)*(2_pInt*Nactive-grI)/2_pInt + grJ - grI - 1_pInt) 
             param(instance)%searchDir(xioffsetIN+1_pInt:xioffsetIN+3_pInt,offset) = &
               residual(xioffsetI+1_pInt:xioffsetI+3_pInt)
           enddo
         enddo
         param(instance)%newIter  (:,offset) = &
           param(instance)%oldIter  (:,offset) - &
           param(instance)%stepLength(offset)* &
           param(instance)%searchDir(:,offset)
       else
        call IO_error(400_pInt,el=el,ip=ip,ext_msg='homogenization multiphase')
       endif        
     else convergencefRankone                                                                       ! not converged and residuum not improved...
       homogenization_multiphase_updateState = [.false., .true.]
       param(instance)%stepLength(offset) = &
         param(instance)%stepLength(offset)/2.0_pReal                                               ! ...try with smaller step length in same direction
       param(instance)%newIter  (:,offset) = &
         param(instance)%oldIter  (:,offset) - &
         param(instance)%stepLength(offset)* &
         param(instance)%searchDir(:,offset)
     endif convergencefRankone

   case(fullrankone_ID) myMixRule

 end select myMixRule

end function homogenization_multiphase_updateState


!--------------------------------------------------------------------------------------------------
!> @brief flux function for each phase 
!--------------------------------------------------------------------------------------------------
function homogenization_multiphase_getPhaseFlux(gradPhi,ip,el)
 use math, only: &
   PI
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   homogenization_maxNgrains, &
   homogenization_typeInstance
 
 implicit none
 real(pReal),   dimension(3,homogenization_maxNgrains) :: &
   homogenization_multiphase_getPhaseFlux, &
   selfFlux
 real(pReal),   dimension(3,homogenization_maxNgrains), intent(in) :: gradPhi
 integer(pInt),                                         intent(in) :: ip, el                !< element number
 integer(pInt) :: &
   grI, grJ, &
   homog, & 
   instance

 homog = material_homog(ip,el)
 instance = homogenization_typeInstance(homog)
 selfFlux = 0.0_pReal
 do grI = 1_pInt, homogenization_Ngrains(homog)
   do grJ = grI + 1_pInt, homogenization_Ngrains(homog)
     selfFlux(1:3,grI) = selfFlux(1:3,grI) - &
        param(instance)%InterfaceEnergy(grI,grJ)*gradPhi(1:3,grJ) 
     selfFlux(1:3,grJ) = selfFlux(1:3,grJ) - &
        param(instance)%InterfaceEnergy(grI,grJ)*gradPhi(1:3,grI) 
   enddo
 enddo    
 
 homogenization_multiphase_getPhaseFlux = 0.0_pReal
 do grI = 1_pInt, homogenization_Ngrains(homog)
   do grJ = 1_pInt, homogenization_Ngrains(homog)
     homogenization_multiphase_getPhaseFlux(1:3,grI) = &
       homogenization_multiphase_getPhaseFlux(1:3,grI) + &
       param(instance)%InterfaceMobility(grI,grJ)* &
       (selfFlux(1:3,grI) - selfFlux(1:3,grJ))
   enddo
 enddo    
 homogenization_multiphase_getPhaseFlux = &
   (8.0_pReal*param(instance)%InterfaceWidth/PI/PI)* &
   homogenization_multiphase_getPhaseFlux

end function homogenization_multiphase_getPhaseFlux


!--------------------------------------------------------------------------------------------------
!> @brief flux tangent function for each phase 
!--------------------------------------------------------------------------------------------------
function homogenization_multiphase_getPhaseFluxTangent(ip,el)
 use math, only: &
   PI, &
   math_I3
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   homogenization_maxNgrains, &
   homogenization_typeInstance
 
 implicit none
 real(pReal),   dimension(homogenization_maxNgrains,homogenization_maxNgrains) :: &
   homogenization_multiphase_getPhaseFluxTangent, &
   selfFluxTangent
 integer(pInt), intent(in) :: ip, el                                                                !< element number
 integer(pInt) :: &
   grI, grJ, grK, &
   homog, & 
   instance

 homog = material_homog(ip,el)
 instance = homogenization_typeInstance(homog)
 selfFluxTangent = 0.0_pReal
 do grI = 1_pInt, homogenization_Ngrains(homog)
   do grJ = grI + 1_pInt, homogenization_Ngrains(homog)
     selfFluxTangent(grI,grJ) = selfFluxTangent(grI,grJ) - &
        param(instance)%InterfaceEnergy(grI,grJ) 
     selfFluxTangent(grJ,grI) = selfFluxTangent(grJ,grI) - &
        param(instance)%InterfaceEnergy(grI,grJ) 
   enddo
 enddo    
 
 homogenization_multiphase_getPhaseFluxTangent = 0.0_pReal
 do grI = 1_pInt, homogenization_Ngrains(homog)
   do grJ = 1_pInt, homogenization_Ngrains(homog)
     do grK = 1_pInt, homogenization_Ngrains(homog)
       homogenization_multiphase_getPhaseFluxTangent(grI,grK) = &
         homogenization_multiphase_getPhaseFluxTangent(grI,grK) + &
         param(instance)%InterfaceMobility(grI,grJ)* &
         (selfFluxTangent(grI,grK) - selfFluxTangent(grJ,grK))
     enddo
   enddo
 enddo    
 homogenization_multiphase_getPhaseFluxTangent = &
   (8.0_pReal*param(instance)%InterfaceWidth/PI/PI)* &
   homogenization_multiphase_getPhaseFluxTangent

end function homogenization_multiphase_getPhaseFluxTangent


!--------------------------------------------------------------------------------------------------
!> @brief source function for each phase 
!--------------------------------------------------------------------------------------------------
function homogenization_multiphase_getPhaseSource(phi,ip,el)
 use material, only: &
   material_homog, &
   phaseAt, &
   phase_source, &
   phase_Nsources, &
   homogenization_Ngrains, &
   homogenization_maxNgrains, &
   homogenization_typeInstance, &
   SOURCE_elastic_energy_ID, &
   SOURCE_plastic_energy_ID, &
   SOURCE_chemical_energy_ID
 use crystallite, only: &
   crystallite_Tstar_v
 use constitutive, only: &
   constitutive_homogenizedC
 use source_elastic_energy, only: &
   source_elastic_energy_getRateAndItsTangent
 use source_plastic_energy, only: &
   source_plastic_energy_getRateAndItsTangent
 use source_chemical_energy, only: &
   source_chemical_energy_getRateAndItsTangent
 
 implicit none
 real(pReal),   dimension(homogenization_maxNgrains) :: &
   homogenization_multiphase_getPhaseSource, &
   selfSource
 real(pReal),   dimension(homogenization_maxNgrains), intent(in) :: phi
 integer(pInt),                                       intent(in) :: ip, el                          !< element number
 integer(pInt) :: &
   grI, grJ, grK, &
   homog, & 
   instance, &
   source, &
   phase
 real(pReal) :: &
   tripleJunctionEnergy, &
   localSource, localSourceTangent  

 homog = material_homog(ip,el)
 instance = homogenization_typeInstance(homog)

 selfSource = 0.0_pReal
 do grI = 1_pInt, homogenization_Ngrains(homog)
   phase = phaseAt(grI,ip,el)
   do source = 1_pInt, phase_Nsources(phase)
     select case(phase_source(source,phase)) 
       case (SOURCE_elastic_energy_ID)
         call source_elastic_energy_getRateAndItsTangent(localSource, localSourceTangent, &
                                                         crystallite_Tstar_v(1:6,grI,ip,el), &
                                                         constitutive_homogenizedC(grI,ip,el))
       
       case (SOURCE_plastic_energy_ID)
         call source_plastic_energy_getRateAndItsTangent(localSource, localSourceTangent, &
                                                         grI,ip,el)
       
       case (SOURCE_chemical_energy_ID)                                                  
         call source_chemical_energy_getRateAndItsTangent(localSource, localSourceTangent, &
                                                          grI,ip,el)
       
       case default
         localSource = 0.0_pReal
       
     end select
     selfSource(grI) = selfSource(grI) + localSource
   enddo  
 enddo
 
 do grI = 1_pInt, homogenization_Ngrains(homog)
   do grJ = grI+1_pInt, homogenization_Ngrains(homog)
     selfSource(grI) = selfSource(grI) - &
       param(instance)%InterfaceEnergy(grI,grJ)*phi(grJ)
     selfSource(grJ) = selfSource(grJ) - &
       param(instance)%InterfaceEnergy(grI,grJ)*phi(grI)
     do grK = grJ+1_pInt, homogenization_Ngrains(homog)
       tripleJunctionEnergy = &
         max(max(param(instance)%InterfaceEnergy(grJ,grK), &
                 param(instance)%InterfaceEnergy(grK,grI)), &
                 param(instance)%InterfaceEnergy(grI,grJ))
       selfSource(grI) = selfSource(grI) - &
         tripleJunctionEnergy*phi(grJ)*phi(grK)
       selfSource(grJ) = selfSource(grJ) - &
         tripleJunctionEnergy*phi(grK)*phi(grI)
       selfSource(grK) = selfSource(grK) - &
         tripleJunctionEnergy*phi(grI)*phi(grJ)
     enddo
   enddo
   if (phi(grI) < param(instance)%phasefracTol) &
     selfSource(grI) = selfSource(grI) + &
       ((phi(grI) - param(instance)%phasefracTol)/param(instance)%phasefracTol)**2.0_pReal
 enddo

 homogenization_multiphase_getPhaseSource = 0.0_pReal
 do grI = 1_pInt, homogenization_Ngrains(homog)
   do grJ = 1_pInt, homogenization_Ngrains(homog)
     homogenization_multiphase_getPhaseSource(grI) = &
       homogenization_multiphase_getPhaseSource(grI) + &
       param(instance)%InterfaceMobility(grI,grJ)* &
       (selfSource(grI) - selfSource(grJ))
   enddo
 enddo    
 homogenization_multiphase_getPhaseSource = &
   (8.0_pReal/param(instance)%InterfaceWidth)* &
   homogenization_multiphase_getPhaseSource
 
end function homogenization_multiphase_getPhaseSource


!--------------------------------------------------------------------------------------------------
!> @brief source tangent function for each phase 
!--------------------------------------------------------------------------------------------------
function homogenization_multiphase_getPhaseSourceTangent(phi,ip,el)
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   homogenization_maxNgrains, &
   homogenization_typeInstance
 
 implicit none
 real(pReal),   dimension(homogenization_maxNgrains,homogenization_maxNgrains) :: &
   homogenization_multiphase_getPhaseSourceTangent, &
   selfSourceTangent
 real(pReal),   dimension(homogenization_maxNgrains), intent(in) :: phi
 integer(pInt), intent(in) :: ip, el                                                                !< element number
 integer(pInt) :: &
   grI, grJ, grK, &
   homog, & 
   instance
 real(pReal) :: &
   tripleJunctionEnergy  

 homog = material_homog(ip,el)
 instance = homogenization_typeInstance(homog)
 selfSourceTangent = 0.0_pReal
 do grI = 1_pInt, homogenization_Ngrains(homog)
   do grJ = grI+1_pInt, homogenization_Ngrains(homog)
     selfSourceTangent(grI,grJ) = selfSourceTangent(grI,grJ) - &
       param(instance)%InterfaceEnergy(grI,grJ)
     selfSourceTangent(grJ,grI) = selfSourceTangent(grJ,grI) - &
       param(instance)%InterfaceEnergy(grI,grJ)
     do grK = grJ+1_pInt, homogenization_Ngrains(homog)
       tripleJunctionEnergy = &
         max(max(param(instance)%InterfaceEnergy(grJ,grK), &
                 param(instance)%InterfaceEnergy(grK,grI)), &
                 param(instance)%InterfaceEnergy(grI,grJ))
       selfSourceTangent(grI,grJ) = selfSourceTangent(grI,grJ) - &
         tripleJunctionEnergy*phi(grK)
       selfSourceTangent(grI,grK) = selfSourceTangent(grI,grK) - &
         tripleJunctionEnergy*phi(grJ)
       selfSourceTangent(grJ,grK) = selfSourceTangent(grJ,grK) - &
         tripleJunctionEnergy*phi(grI)
       selfSourceTangent(grJ,grI) = selfSourceTangent(grJ,grI) - &
         tripleJunctionEnergy*phi(grK)
       selfSourceTangent(grK,grI) = selfSourceTangent(grK,grI) - &
         tripleJunctionEnergy*phi(grJ)
       selfSourceTangent(grK,grJ) = selfSourceTangent(grK,grJ) - &
         tripleJunctionEnergy*phi(grI)
     enddo
   enddo  
   if (phi(grI) < param(instance)%phasefracTol) &
     selfSourceTangent(grI,grI) = selfSourceTangent(grI,grI) + &
       2.0_pReal*(phi(grI) - param(instance)%phasefracTol)/param(instance)%phasefracTol/param(instance)%phasefracTol
 enddo
 
 homogenization_multiphase_getPhaseSourceTangent = 0.0_pReal
 do grI = 1_pInt, homogenization_Ngrains(homog)
   do grJ = 1_pInt, homogenization_Ngrains(homog)
     do grK = 1_pInt, homogenization_Ngrains(homog)
       homogenization_multiphase_getPhaseSourceTangent(grI,grK) = &
         homogenization_multiphase_getPhaseSourceTangent(grI,grK) + &
         param(instance)%InterfaceMobility(grI,grJ)* &
         (selfSourceTangent(grI,grK) - selfSourceTangent(grJ,grK)) 
     enddo
   enddo
 enddo
 homogenization_multiphase_getPhaseSourceTangent = &
   (8.0_pReal/param(instance)%InterfaceWidth)* &
   homogenization_multiphase_getPhaseSourceTangent

end function homogenization_multiphase_getPhaseSourceTangent


!--------------------------------------------------------------------------------------------------
!> @brief set module wide phase fractions 
!--------------------------------------------------------------------------------------------------
subroutine homogenization_multiphase_putPhaseFrac(frac,ip,el)
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   homogenization_typeInstance, &
   phasefracMapping, &
   phasefrac
 
 implicit none
 integer(pInt),                                                intent(in)  :: ip, el                !< element number
 real(pReal),   dimension (homogenization_Ngrains(material_homog(ip,el))),     &
                                                               intent(in)  :: frac             !< array of current phase fractions
 integer(pInt) :: &
   gr, &
   homog, & 
   instance, &
   offset

 homog = material_homog(ip,el)
 instance = homogenization_typeInstance(homog)
 offset = phasefracMapping(homog)%p(ip,el)
 do gr = 1_pInt, homogenization_Ngrains(homog)
   phasefrac(gr,homog)%p(offset) =  frac(gr) 
 enddo

end subroutine homogenization_multiphase_putPhaseFrac


!--------------------------------------------------------------------------------------------------
!> @brief set interface normal 
!--------------------------------------------------------------------------------------------------
subroutine homogenization_multiphase_putInterfaceNormals(interfaceNormal,grI,grJ,ip,el)
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   homogenization_typeInstance, &
   phasefracMapping
 
 implicit none
 integer(pInt),                                                intent(in)  :: ip, el                !< element number
 integer(pInt),                                                intent(in)  :: grI, grJ              !< interface between grain I and J
 real(pReal),   dimension (3),                                 intent(in)  :: interfaceNormal       !< interface normal
 integer(pInt) :: &
   homog, & 
   instance, &
   offset, &
   xioffset

 homog = material_homog(ip,el)
 instance = homogenization_typeInstance(homog)
 offset = phasefracMapping(homog)%p(ip,el)
 xioffset = 3_pInt*((grI-1_pInt)*(2_pInt*homogenization_Ngrains(homog)-grI)/2_pInt + grJ - grI - 1_pInt)
 if (param(instance)%mixtureID == partialrankone_ID) &
   param(instance)%interfaceNormal(xioffset+1:xioffset+3,offset) = interfaceNormal

end subroutine homogenization_multiphase_putInterfaceNormals


!--------------------------------------------------------------------------------------------------
!> @brief return component concentration at material point 
!--------------------------------------------------------------------------------------------------
function homogenization_multiphase_getComponentConc(chempot,conc0,ip,el)
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   material_phase, &
   material_maxNcomponents, &
   phase_chemicalFE, &
   CHEMICALFE_quadenergy_ID, &
   CHEMICALFE_thermodynamic_ID, &
   phasefracMapping, &
   phasefrac
 use chemicalFE_quadenergy, only: &
   chemicalFE_quadenergy_getConc  
 use chemicalFE_thermodynamic, only: &
   chemicalFE_thermodynamic_getConc  
 
 implicit none
 integer(pInt),                                   intent(in) :: &
   ip, el                                                                                      !< element number
 real(pReal), dimension(material_maxNcomponents), intent(in) :: &
   chempot, conc0
 real(pReal), dimension(material_maxNcomponents) :: &
   homogenization_multiphase_getComponentConc
 integer(pInt) :: &
   gr, &
   homog, & 
   offset

 homog = material_homog(ip,el)
 offset = phasefracMapping(homog)%p(ip,el)
 homogenization_multiphase_getComponentConc = 0.0_pReal
 do gr = 1_pInt, homogenization_Ngrains(homog)
   chemicalFEType: select case (phase_chemicalFE(material_phase(gr,ip,el)))
     case (CHEMICALFE_quadenergy_ID) chemicalFEType
       homogenization_multiphase_getComponentConc = &
         homogenization_multiphase_getComponentConc + &
         phasefrac(gr,homog)%p(offset)* &
         chemicalFE_quadenergy_getConc(chempot,gr,ip,el)

     case (CHEMICALFE_thermodynamic_ID) chemicalFEType
       homogenization_multiphase_getComponentConc = &
         homogenization_multiphase_getComponentConc + &
         phasefrac(gr,homog)%p(offset)* &
         chemicalFE_thermodynamic_getConc(chempot,conc0,300.0_pReal,gr,ip,el)

   end select chemicalFEType
 enddo

end function homogenization_multiphase_getComponentConc


!--------------------------------------------------------------------------------------------------
!> @brief return component concentration tangent at material point 
!--------------------------------------------------------------------------------------------------
function homogenization_multiphase_getComponentConcTangent(chempot,conc0,ip,el)
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   material_phase, &
   material_maxNcomponents, &
   phase_chemicalFE, &
   CHEMICALFE_quadenergy_ID, &
   CHEMICALFE_thermodynamic_ID, &
   phasefracMapping, &
   phasefrac
 use chemicalFE_quadenergy, only: &
   chemicalFE_quadenergy_getConcTangent  
 use chemicalFE_thermodynamic, only: &
   chemicalFE_thermodynamic_getConcTangent  
 
 implicit none
 integer(pInt),                                   intent(in) :: &
   ip, el                                                                                      !< element number
 real(pReal), dimension(material_maxNcomponents), intent(in) :: &
   chempot, conc0
 real(pReal), dimension(material_maxNcomponents,material_maxNcomponents) :: &
   homogenization_multiphase_getComponentConcTangent
 integer(pInt) :: &
   gr, &
   homog, & 
   offset

 homog = material_homog(ip,el)
 offset = phasefracMapping(homog)%p(ip,el)
 homogenization_multiphase_getComponentConcTangent = 0.0_pReal
 do gr = 1_pInt, homogenization_Ngrains(homog)
   chemicalFEType: select case (phase_chemicalFE(material_phase(gr,ip,el)))
     case (CHEMICALFE_quadenergy_ID) chemicalFEType
       homogenization_multiphase_getComponentConcTangent = &
         homogenization_multiphase_getComponentConcTangent + &
         phasefrac(gr,homog)%p(offset)* &
         chemicalFE_quadenergy_getConcTangent(gr,ip,el)

     case (CHEMICALFE_thermodynamic_ID) chemicalFEType
       homogenization_multiphase_getComponentConcTangent = &
         homogenization_multiphase_getComponentConcTangent + &
         phasefrac(gr,homog)%p(offset)* &
         chemicalFE_thermodynamic_getConcTangent(chempot,conc0,300.0_pReal,gr,ip,el)

   end select chemicalFEType
 enddo

end function homogenization_multiphase_getComponentConcTangent


!--------------------------------------------------------------------------------------------------
!> @brief return component mobility at material point 
!--------------------------------------------------------------------------------------------------
function homogenization_multiphase_getComponentMobility(ip,el)
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   material_phase, &
   material_maxNcomponents, &
   phase_chemicalFE, &
   CHEMICALFE_quadenergy_ID, &
   CHEMICALFE_thermodynamic_ID, &
   phasefracMapping, &
   phasefrac
 use chemicalFE_quadenergy, only: &
   chemicalFE_quadenergy_getMobility  
 use chemicalFE_thermodynamic, only: &
   chemicalFE_thermodynamic_getMobility  
 
 implicit none
 integer(pInt),                                   intent(in) :: &
   ip, el                                                                                      !< element number
 real(pReal), dimension(material_maxNcomponents) :: &
   homogenization_multiphase_getComponentMobility
 integer(pInt) :: &
   gr, &
   homog, & 
   offset

 homog = material_homog(ip,el)
 offset = phasefracMapping(homog)%p(ip,el)
 homogenization_multiphase_getComponentMobility = 0.0_pReal
 do gr = 1_pInt, homogenization_Ngrains(homog)
   chemicalFEType: select case (phase_chemicalFE(material_phase(gr,ip,el)))
     case (CHEMICALFE_quadenergy_ID) chemicalFEType
       homogenization_multiphase_getComponentMobility = &
         homogenization_multiphase_getComponentMobility + &
         phasefrac(gr,homog)%p(offset)* &
         chemicalFE_quadenergy_getMobility(gr,ip,el)

     case (CHEMICALFE_thermodynamic_ID) chemicalFEType
       homogenization_multiphase_getComponentMobility = &
         homogenization_multiphase_getComponentMobility + &
         phasefrac(gr,homog)%p(offset)* &
         chemicalFE_thermodynamic_getMobility(gr,ip,el)

   end select chemicalFEType
 enddo

end function homogenization_multiphase_getComponentMobility


!--------------------------------------------------------------------------------------------------
!> @brief set module wide component concentrations 
!--------------------------------------------------------------------------------------------------
subroutine homogenization_multiphase_putComponentConc(chempot,conc0,ip,el)
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   material_phase, &
   material_maxNcomponents, &
   phase_chemicalFE, &
   CHEMICALFE_quadenergy_ID, &
   CHEMICALFE_thermodynamic_ID
 use chemicalFE_quadenergy, only: &
   chemicalFE_quadenergy_putConc  
 use chemicalFE_thermodynamic, only: &
   chemicalFE_thermodynamic_putConc  
 
 implicit none
 integer(pInt),                                   intent(in) :: &
   ip, el                                                                                      !< element number
 real(pReal), dimension(material_maxNcomponents), intent(in) :: &
   chempot, conc0
 integer(pInt) :: &
   gr, &
   homog

 homog = material_homog(ip,el)
 do gr = 1_pInt, homogenization_Ngrains(homog)
   chemicalFEType: select case (phase_chemicalFE(material_phase(gr,ip,el)))
     case (CHEMICALFE_quadenergy_ID) chemicalFEType
       call chemicalFE_quadenergy_putConc(chempot,gr,ip,el)

     case (CHEMICALFE_thermodynamic_ID) chemicalFEType
       call chemicalFE_thermodynamic_putConc(chempot,conc0,300.0_pReal,gr,ip,el)

   end select chemicalFEType
 enddo

end subroutine homogenization_multiphase_putComponentConc


!--------------------------------------------------------------------------------------------------
!> @brief return array of homogenization results for post file inclusion 
!--------------------------------------------------------------------------------------------------
pure function homogenization_multiphase_postResults(ip,el,avgP,avgF)
 use mesh, only: &
   mesh_ipCoordinates
 use material, only: &
   phasefrac, &
   phasefracMapping, &
   material_homog, &
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
                         (homogenization_typeInstance(material_homog(ip,el)))) :: &
   homogenization_multiphase_postResults
 
 integer(pInt) :: &
   homID, &
   o, c, gr
   
 c = 0_pInt
 homID = homogenization_typeInstance(material_homog(ip,el))
 homogenization_multiphase_postResults = 0.0_pReal
 
 do o = 1_pInt,homogenization_Noutput(material_homog(ip,el))
   select case(param(homID)%outputID(o))
     case (nconstituents_ID)
       homogenization_multiphase_postResults(c+1_pInt) = real(homogenization_Ngrains(material_homog(ip,el)),pReal)
       c = c + 1_pInt
     case (avgdefgrad_ID)
       homogenization_multiphase_postResults(c+1_pInt:c+9_pInt) = reshape(avgF,[9])
       c = c + 9_pInt
     case (avgfirstpiola_ID)
       homogenization_multiphase_postResults(c+1_pInt:c+9_pInt) = reshape(avgP,[9])
       c = c + 9_pInt
     case (phasefrac_ID)
       do gr = 1_pInt, homogenization_Ngrains(material_homog(ip,el))
         homogenization_multiphase_postResults(c+1_pInt) = &
         phasefrac(gr,material_homog(ip,el))%p(phasefracMapping(material_homog(ip,el))%p(ip,el))
         c = c + 1_pInt
       enddo  
     case (ipcoords_ID)
       homogenization_multiphase_postResults(c+1_pInt:c+3_pInt) = mesh_ipCoordinates(1:3,ip,el)                       ! current ip coordinates
       c = c + 3_pInt
    end select
 enddo

end function homogenization_multiphase_postResults

end module homogenization_multiphase