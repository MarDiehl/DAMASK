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
 integer(pInt),             dimension(:,:),   allocatable, target, public :: &
   homogenization_multiphase_sizePostResult
 
 character(len=64),         dimension(:,:),   allocatable, target, public :: &
  homogenization_multiphase_output                                                                  !< name of each post result output
 integer(pInt),             dimension(:),     allocatable, target, public :: &
   homogenization_multiphase_Noutput                                                                !< number of outputs per homog instance

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

 enum, bind(c) 
   enumerator :: linear_ID, &
                 moelans_ID
 end enum

 type, private :: tParameters                                                                       !< container type for internal constitutive parameters
   integer(kind(undefined_ID)), dimension(:), allocatable :: & 
     outputID
   integer(kind(isostrain_ID)) :: &
     mixtureID = isostrain_ID
   integer(kind(linear_ID)) :: &
     interpolationID = linear_ID
   real(pReal), dimension(:,:), allocatable :: &
     interfaceMobility, &
     interfaceEnergy, &
     interfaceNormalA, &                                                                            ! orientation of interface normal, alpha, in rad
     interfaceNormalB                                                                               ! orientation of interface normal, beta, in rad
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

 type(tParameters), dimension(:), allocatable, private :: param                                     !< containers of constitutive parameters (len Ninstance)

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
   homogenization_multiphase_calcPhaseFrac, &
   homogenization_multiphase_putInterfaceNormals, &
   homogenization_multiphase_postResults

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine homogenization_multiphase_init
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use debug, only: &
   debug_HOMOGENIZATION, &
   debug_level, &
   debug_levelBasic
 use IO, only: &
   IO_error, &
   IO_timeStamp
 use math, only: &
   math_I3, &
   math_inv33, &
   math_mul33x33, &
   math_mul33x3 , &   
   math_transpose33, &
   math_EulerToR, &
   math_VecMtoSymNN
 use lattice, only: &
   lattice_initialPlasticStrain
 use material
 use config
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
 integer(pInt) :: &
   o, el, ip, grI, grJ, outputSize
 integer :: &
   maxNinstance, &
   homog, &
   micro, &
   phase, &
   instance, &
   offset, &
   xioffset
 integer :: &
   NofMyHomog                                                                                       ! no pInt (stores a system dependen value from 'count'
 character(len=65536) :: &
   tag  = ''
 real(pReal) :: &
   interfaceNormal(3)  
 character(len=65536), dimension(:), allocatable :: outputs
 character(len=65536), dimension(0), parameter :: emptyStringArray = [character(len=65536)::]
 integer(kind(undefined_ID)) :: outputID                                                           !< ID of each post result output
 
 write(6,'(/,a)')   ' <<<+-  homogenization_'//HOMOGENIZATION_MULTIPHASE_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = count(homogenization_type == HOMOGENIZATION_MULTIPHASE_ID)
 if (maxNinstance == 0) return
 
 if (iand(debug_level(debug_HOMOGENIZATION),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 allocate(homogenization_multiphase_sizePostResult(maxval(homogenization_Noutput),maxNinstance), &
                                                                            source=0_pInt)
 allocate(homogenization_multiphase_Noutput(maxNinstance),                  source=0_pInt)
 allocate(homogenization_multiphase_output(maxval(homogenization_Noutput),maxNinstance))
          homogenization_multiphase_output = ''
 
 allocate(param(maxNinstance))

 parsingFile: do homog = 1_pInt, material_Nhomogenization
   if (homogenization_type(homog) /= HOMOGENIZATION_MULTIPHASE_ID) cycle
   instance = homogenization_typeInstance(homog)
   tag = config_homogenization(homog)%getString('mixture_rule', defaultVal='isostrain')
   select case(trim(tag))
     case('isostrain')
       param(instance)%mixtureID = isostrain_ID
     case('isofpkstress')
       param(instance)%mixtureID = isofPKstress_ID
     case('isocauchystress')
       param(instance)%mixtureID = isoCauchystress_ID
     case('partialrankone')
       param(instance)%mixtureID = partialrankone_ID
     case('fullrankone')
       param(instance)%mixtureID = fullrankone_ID
     case default
       call IO_error(211_pInt,el=instance,ext_msg='mixture_rule ('//HOMOGENIZATION_multiphase_label//')')
   end select
   tag = config_homogenization(homog)%getString('interpolation_type', defaultVal='linear')
   select case(trim(tag))
     case('linear')
       param(instance)%interpolationID = linear_ID
     case('moelans')
       param(instance)%interpolationID = moelans_ID
     case default
       call IO_error(211_pInt,el=instance,ext_msg='interpolation_type ('//HOMOGENIZATION_multiphase_label//')')
   end select
   param(instance)%interfaceMobility = &
     math_VecMtoSymNN(config_homogenization(homog)%getFloats('interface_mobility'), &
                      homogenization_Ngrains(homog))
   param(instance)%interfaceEnergy = &
     math_VecMtoSymNN(config_homogenization(homog)%getFloats('interface_energy'), &
                      homogenization_Ngrains(homog))
   param(instance)%interfaceNormalA = &
     math_VecMtoSymNN(config_homogenization(homog)%getFloats('interfacenormala'), &
                      homogenization_Ngrains(homog))
   param(instance)%interfaceNormalB = &
     math_VecMtoSymNN(config_homogenization(homog)%getFloats('interfacenormalb'), &
                      homogenization_Ngrains(homog))
   param(instance)%InterfaceWidth = config_homogenization(homog)%getFloat('interface_width',defaultVal=0.0_pReal)
   param(instance)%absTol = config_homogenization(homog)%getFloat('abs_tol',defaultVal=0.0_pReal)
   param(instance)%relTol = config_homogenization(homog)%getFloat('rel_tol',defaultVal=0.0_pReal)
   param(instance)%phaseFracTol = config_homogenization(homog)%getFloat('phasefrac_tol',defaultVal=1e-4_pReal)
   
   outputs = config_homogenization(homog)%getStrings('(output)',defaultVal=emptyStringArray)
   allocate(param(instance)%outputID(0))
   do o=1_pInt, size(outputs)
     outputID = undefined_ID
     select case(outputs(o))
       case('nconstituents','ngrains')
         outputID = nconstituents_ID
         outputSize = 1_pInt
       case('ipcoords')
         outputID = ipcoords_ID
         outputSize = 3_pInt
       case('avgdefgrad','avgf')
         outputID = avgdefgrad_ID
         outputSize = 9_pInt
       case('avgp','avgfirstpiola','avg1stpiola')
         outputID = avgfirstpiola_ID
         outputSize = 9_pInt
       case('phasefrac','phasefraction')
         outputID = phasefrac_ID
         outputSize = homogenization_Ngrains(homog)
      end select

      if (outputID /= undefined_ID) then
        homogenization_multiphase_output(o,instance) = outputs(o)
        homogenization_multiphase_sizePostResult(o,instance) = outputSize
        param(instance)%outputID = [param(instance)%outputID , outputID]
      endif
   enddo
 enddo parsingFile

 sanityChecks: do homog = 1_pInt, material_Nhomogenization
   if (homogenization_type(homog) /= HOMOGENIZATION_MULTIPHASE_ID) cycle
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
 enddo sanityChecks

 initializeInstances: do homog = 1_pInt, material_Nhomogenization
   if (homogenization_type(homog) /= HOMOGENIZATION_MULTIPHASE_ID) cycle
   NofMyHomog = count(material_homog == homog)
   instance = homogenization_typeInstance(homog)

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
       homogState(homog)%sizeState = &    
         15*homogenization_Ngrains(homog)*(homogenization_Ngrains(homog) - 1_pInt)/2_pInt + 2_pInt

   end select
   homogState(homog)%sizePostResults = sum(homogenization_multiphase_sizePostResult(:,instance))
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
       xioffset = homogenization_Ngrains(homog)*(homogenization_Ngrains(homog) - 1_pInt)/2_pInt
       param(instance)%newIter   => homogState(homog)% &
                                      state(xioffset* 0_pInt + 1_pInt:xioffset* 5_pInt,1:NofMyHomog)
       param(instance)%oldIter   => homogState(homog)% &
                                      state(xioffset* 5_pInt + 1_pInt:xioffset* 10_pInt,1:NofMyHomog)
       param(instance)%searchDir => homogState(homog)% &
                                      state(xioffset* 10_pInt + 1_pInt:xioffset* 15_pInt,1:NofMyHomog)
       param(instance)%residual  => homogState(homog)%state(xioffset* 15_pInt + 1_pInt,1:NofMyHomog)
       param(instance)%stepLength=> homogState(homog)%state(xioffset* 15_pInt + 2_pInt,1:NofMyHomog)

   end select

   phasefracMapping(homog)%p => mappingHomogenization(1,:,:)
   allocate  (phasefrac(homog)%p(homogenization_Ngrains(homog),NofMyHomog))

 enddo initializeInstances
 
! state init
 elementLoop: do el = 1_pInt, mesh_NcpElems
   IpLoop: do ip = 1_pInt, FE_Nips(FE_geomtype(mesh_element(2,el)))
     homog = material_homog(ip,el)
     micro = mesh_element(4,el)
     instance = homogenization_typeInstance(homog)
     myHomog3: if (homogenization_type(homog) == HOMOGENIZATION_MULTIPHASE_ID) then
       do grI = 1_pInt, homogenization_Ngrains(homog)
         phase = material_phase(grI,ip,el)
         phasefrac(homog)%p(grI,phasefracMapping(homog)%p(ip,el)) = microstructure_fraction(grI,micro)
       enddo
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
               xioffset = 5_pInt*((grI-1_pInt)*(2_pInt*homogenization_Ngrains(homog)-grI)/2_pInt + grJ - grI - 1_pInt) 
               param(instance)%newIter(xioffset+1_pInt,offset) = param(instance)%interfaceNormalA(grI,grJ)
               param(instance)%newIter(xioffset+2_pInt,offset) = param(instance)%interfaceNormalB(grI,grJ)
               interfaceNormal(1)= &
                   sin(param(instance)%newIter(xioffset+1_pInt,offset))* &
                   cos(param(instance)%newIter(xioffset+2_pInt,offset))
               interfaceNormal(2)= &
                   sin(param(instance)%newIter(xioffset+1_pInt,offset))* &
                   sin(param(instance)%newIter(xioffset+2_pInt,offset))
               interfaceNormal(3)= &
                   cos(param(instance)%newIter(xioffset+1_pInt,offset))                                            
               param(instance)%newIter(xioffset+3:xioffset+5,offset) = &
                math_mul33x3(crystallite_F0(1:3,1:3,grI,ip,el) - &
                             crystallite_F0(1:3,1:3,grJ,ip,el), &
                             interfaceNormal) 
               param(instance)%oldIter(xioffset+1:xioffset+5,offset) = &
                param(instance)%newIter(xioffset+1:xioffset+5,offset)
             enddo
           enddo

           homogState(homog)%state0   (1:homogState(homog)%sizeState,offset) = &
             homogState(homog)%state(1:homogState(homog)%sizeState,offset)
           homogState(homog)%subState0(1:homogState(homog)%sizeState,offset) = &
             homogState(homog)%state(1:homogState(homog)%sizeState,offset)


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
 real(pReal),   dimension (3)                             :: interfaceNormal 
 
 homog = material_homog(ip,el)
 instance = homogenization_typeInstance(homog)
 offset = phasefracMapping(homog)%p(ip,el)
 select case(param(instance)%mixtureID)
   case(isostrain_ID)
     do grI = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(homog)%p(grI,offset) > err_phasefr_tolabs) &
         F(1:3,1:3,grI) = avgF
     enddo
   
   case(isofPKstress_ID)
     do grI = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(homog)%p(grI,offset) > err_phasefr_tolabs) &
         F(1:3,1:3,grI) = avgF
     enddo
     
     do grI = 1_pInt, homogenization_Ngrains(homog)-1_pInt
       do grJ = grI+1_pInt, homogenization_Ngrains(homog)
         if (      phasefrac(homog)%p(grI,offset) > err_phasefr_tolabs &
             .and. phasefrac(homog)%p(grJ,offset) > err_phasefr_tolabs) then
           xioffset = 9_pInt*((grI-1_pInt)*(2_pInt*homogenization_Ngrains(homog)-grI)/2_pInt + grJ - grI - 1_pInt)
           F(1:3,1:3,grI) = &
             F(1:3,1:3,grI) + &
             phasefrac(homog)%p(grJ,offset)* &
             reshape(param(instance)%newIter(xioffset+1:xioffset+9,offset), shape = [3,3])
           F(1:3,1:3,grJ) = &
             F(1:3,1:3,grJ) - &
             phasefrac(homog)%p(grI,offset)* &
             reshape(param(instance)%newIter(xioffset+1:xioffset+9,offset), shape = [3,3])
         endif
       enddo
     enddo

   case(isoCauchystress_ID)
     avgLdt = math_mul33x33(avgF - avgF0,math_inv33(avgF))
     do grI = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(homog)%p(grI,offset) > err_phasefr_tolabs) &
         Ldt(1:3,1:3,grI) = avgLdt
     enddo
     
     do grI = 1_pInt, homogenization_Ngrains(homog)-1_pInt
       do grJ = grI+1_pInt, homogenization_Ngrains(homog)
         if (      phasefrac(homog)%p(grI,offset) > err_phasefr_tolabs &
             .and. phasefrac(homog)%p(grJ,offset) > err_phasefr_tolabs) then
           xioffset = 9_pInt*((grI-1_pInt)*(2_pInt*homogenization_Ngrains(homog)-grI)/2_pInt + grJ - grI - 1_pInt)
           Ldt(1:3,1:3,grI) = &
             Ldt(1:3,1:3,grI) + &
             phasefrac(homog)%p(grJ,offset)* &
             math_Mandel6to33(param(instance)%newIter(xioffset+1:xioffset+6,offset))
           Ldt(1:3,1:3,grJ) = &
             Ldt(1:3,1:3,grJ) - &
             phasefrac(homog)%p(grI,offset)* &
             math_Mandel6to33(param(instance)%newIter(xioffset+1:xioffset+6,offset))
         endif
       enddo
     enddo

     do grI = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(homog)%p(grI,offset) > err_phasefr_tolabs) &
         F(1:3,1:3,grI) = math_mul33x33(math_inv33(math_I3 - Ldt(1:3,1:3,grI)),F0(1:3,1:3,grI))
     enddo                          

   case(partialrankone_ID)
     do grI = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(homog)%p(grI,offset) > err_phasefr_tolabs) &
         F(1:3,1:3,grI) = avgF
     enddo
     
     do grI = 1_pInt, homogenization_Ngrains(homog)-1_pInt
       do grJ = grI+1_pInt, homogenization_Ngrains(homog)
         if (      phasefrac(homog)%p(grI,offset) > err_phasefr_tolabs &
             .and. phasefrac(homog)%p(grJ,offset) > err_phasefr_tolabs) then
           xioffset = 3_pInt*((grI-1_pInt)*(2_pInt*homogenization_Ngrains(homog)-grI)/2_pInt + grJ - grI - 1_pInt)
           F(1:3,1:3,grI) = &
             F(1:3,1:3,grI) + &
             phasefrac(homog)%p(grJ,offset)* &
             math_tensorproduct33(param(instance)%newIter(xioffset+1:xioffset+3,offset), &
                                  param(instance)%interfaceNormal(xioffset+1_pInt:xioffset+3_pInt,offset))
           F(1:3,1:3,grJ) = &
             F(1:3,1:3,grJ) - &
             phasefrac(homog)%p(grI,offset)* &
             math_tensorproduct33(param(instance)%newIter(xioffset+1:xioffset+3,offset), &
                                  param(instance)%interfaceNormal(xioffset+1_pInt:xioffset+3_pInt,offset))
         endif
       enddo
     enddo

   case(fullrankone_ID)
     do grI = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(homog)%p(grI,offset) > err_phasefr_tolabs) &
         F(1:3,1:3,grI) = avgF
     enddo

     do grI = 1_pInt, homogenization_Ngrains(homog)-1_pInt
       do grJ = grI+1_pInt, homogenization_Ngrains(homog)
         if (      phasefrac(homog)%p(grI,offset) > err_phasefr_tolabs &
             .and. phasefrac(homog)%p(grJ,offset) > err_phasefr_tolabs) then
           xioffset = 5_pInt*((grI-1_pInt)*(2_pInt*homogenization_Ngrains(homog)-grI)/2_pInt + grJ - grI - 1_pInt)
           interfaceNormal(1) = &
             sin(param(instance)%newIter(xioffset+1_pInt,offset))* &
             cos(param(instance)%newIter(xioffset+2_pInt,offset))
           interfaceNormal(2) = &
             sin(param(instance)%newIter(xioffset+1_pInt,offset))* &
             sin(param(instance)%newIter(xioffset+2_pInt,offset))
           interfaceNormal(3) = &
             cos(param(instance)%newIter(xioffset+1_pInt,offset))           
           F(1:3,1:3,grI) = &
             F(1:3,1:3,grI) + &
             phasefrac(homog)%p(grJ,offset)* &
             math_tensorproduct33(param(instance)%newIter(xioffset+3:xioffset+5,offset),interfaceNormal)
           F(1:3,1:3,grJ) = &
             F(1:3,1:3,grJ) - &
             phasefrac(homog)%p(grI,offset)* &
             math_tensorproduct33(param(instance)%newIter(xioffset+3:xioffset+5,offset),interfaceNormal)
         endif
       enddo
     enddo   

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
       if (phasefrac(homog)%p(gr,offset) > err_phasefr_tolabs) then
         avgP = avgP + &
            phasefrac(homog)%p(gr,offset)*P(1:3,1:3,gr) 
         dAvgPdAvgF = dAvgPdAvgF + &
              phasefrac(homog)%p(gr,offset)*dPdF(1:3,1:3,1:3,1:3,gr) 
       endif
     enddo   

   case(isofPKstress_ID)                
     avgP = 0.0_pReal
     dAvgFdP = 0.0_pReal
     do gr = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(homog)%p(gr,offset) > err_phasefr_tolabs) then
         avgP = P(1:3,1:3,gr)
         dAvgFdP = dAvgFdP + phasefrac(homog)%p(gr,offset)* &
                             math_invSym3333(dPdF(1:3,1:3,1:3,1:3,gr))
       endif
     enddo 
     dAvgPdAvgF = math_invSym3333(dAvgFdP)

   case(isoCauchystress_ID)                        
     dAvgLdCauchy = 0.0_pReal
     do gr = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(homog)%p(gr,offset) > err_phasefr_tolabs) then
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
            phasefrac(homog)%p(gr,offset)*dLdCauchy 
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
       if (phasefrac(homog)%p(gr,offset) > err_phasefr_tolabs) then
         avgP = avgP + &
            phasefrac(homog)%p(gr,offset)*P(1:3,1:3,gr) 
         dAvgPdAvgF = dAvgPdAvgF + &
              phasefrac(homog)%p(gr,offset)*dPdF(1:3,1:3,1:3,1:3,gr) 
       endif
     enddo

   case(fullrankone_ID)
     avgP = 0.0_pReal
     dAvgPdAvgF = 0.0_pReal
     do gr = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(homog)%p(gr,offset) > err_phasefr_tolabs) then
         avgP = avgP + &
            phasefrac(homog)%p(gr,offset)*P(1:3,1:3,gr) 
         dAvgPdAvgF = dAvgPdAvgF + &
              phasefrac(homog)%p(gr,offset)*dPdF(1:3,1:3,1:3,1:3,gr) 
       endif
     enddo
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
   math_identity2nd, &
   math_mul33xx33, &
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
   avgR
 real(pReal)                                                             :: &
   interfaceNormalAIJ, &
   interfaceNormalBIJ, &
   interfaceNormalAKL, &
   interfaceNormalBKL                                                                               !< normal angles for interface IJ and KL
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
       if (phasefrac(homog)%p(grI,offset) > err_phasefr_tolabs) Nactive = Nactive + 1_pInt
     enddo
     allocate(active(Nactive), source = 0_pInt)

     avgR = 0.0_pReal
     Nactive = 0_pInt
     do grI = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(homog)%p(grI,offset) > err_phasefr_tolabs) then
         Nactive = Nactive + 1_pInt
         active(Nactive) = grI
         avgR = avgR + phasefrac(homog)%p(grI,offset)*P(1:3,1:3,grI)
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
             phasefrac(homog)%p(active(grJ),offset)* &
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
       if (phasefrac(homog)%p(grI,offset) > err_phasefr_tolabs) Nactive = Nactive + 1_pInt
     enddo
     allocate(active(Nactive), source = 0_pInt)
     allocate(cauchy(6,Nactive), source=0.0_pReal)
     allocate(detFiInv(Nactive), source=0.0_pReal)

     avgR = 0.0_pReal
     Nactive = 0_pInt
     do grI = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(homog)%p(grI,offset) > err_phasefr_tolabs) then
         Nactive = Nactive + 1_pInt
         active(Nactive) = grI
         detFiInv(Nactive) = 1.0_pReal/math_det33(F(1:3,1:3,grI))
         cauchy(1:6,Nactive) = &
           math_Mandel33to6(math_mul33x33(P(1:3,1:3,grI),math_transpose33(F(1:3,1:3,grI))))* &
           detFiInv(Nactive)
         avgR = avgR + phasefrac(homog)%p(grI,offset)*math_Mandel6to33(cauchy(1:6,Nactive))
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
             phasefrac(homog)%p(active(grJ),offset)* (dS_dL(1:6,1:6,grI) - dS_dL(1:6,1:6,Nactive))
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
       if (phasefrac(homog)%p(grI,offset) > err_phasefr_tolabs) Nactive = Nactive + 1_pInt
     enddo
     allocate(active(Nactive), source = 0_pInt)
    
     avgR = 0.0_pReal
     Nactive = 0_pInt
     do grI = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(homog)%p(grI,offset) > err_phasefr_tolabs) then
         Nactive = Nactive + 1_pInt 
         active(Nactive) = grI  
         avgR = avgR + phasefrac(homog)%p(grI,offset)*P(1:3,1:3,grI)
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
           phasefrac(homog)%p(active(grI),offset)* &
           phasefrac(homog)%p(active(grJ),offset)* &
           math_mul33x3(P(1:3,1:3,active(grI)) - P(1:3,1:3,active(grJ)), &
                        param(instance)%interfaceNormal(xioffsetIN+1_pInt:xioffsetIN+3_pInt,offset))    
       enddo
     enddo 
     
     stressTol = max(            param(instance)%absTol, &
                     norm2(avgR)*param(instance)%relTol)
     convergencePartialRankOne: if (     norm2(residual) < stressTol &
                              .or. Nactive         < 2_pInt) then
       homogenization_multiphase_updateState = [.true., .true.]  
     elseif (     iter == 1_pInt &
             .or. norm2(residual) < param(instance)%residual(offset)) then convergencePartialRankOne! not converged, but improved norm of residuum (always proceed in first iteration)...
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
                   phasefrac(homog)%p(active(grII),offset)* &
                   phasefrac(homog)%p(active(grIJ),offset)* &
                   phasefrac(homog)%p(active(grJJ),offset)* &
                   sum(dPdF(ii,1:3,jj,1:3,active(grII))* &
                       math_tensorproduct33(param(instance)%interfaceNormal(xioffsetIN+1_pInt:xioffsetIN+3_pInt,offset), &
                                            param(instance)%interfaceNormal(xioffsetJN+1_pInt:xioffsetJN+3_pInt,offset)))                            
           if (grII == grJJ) &
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                 
                 jacobian(xioffsetI+ii,xioffsetJ+jj) = &
                   jacobian(xioffsetI+ii,xioffsetJ+jj) -  &
                   phasefrac(homog)%p(active(grII),offset)* &
                   phasefrac(homog)%p(active(grIJ),offset)* &
                   phasefrac(homog)%p(active(grJI),offset)* &
                   sum(dPdF(ii,1:3,jj,1:3,active(grII))* &
                       math_tensorproduct33(param(instance)%interfaceNormal(xioffsetIN+1_pInt:xioffsetIN+3_pInt,offset), &
                                            param(instance)%interfaceNormal(xioffsetJN+1_pInt:xioffsetJN+3_pInt,offset)))                            
           if (grIJ == grJI) &
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                 
                 jacobian(xioffsetI+ii,xioffsetJ+jj) = &
                   jacobian(xioffsetI+ii,xioffsetJ+jj) -  &
                   phasefrac(homog)%p(active(grII),offset)* &
                   phasefrac(homog)%p(active(grIJ),offset)* &
                   phasefrac(homog)%p(active(grJJ),offset)* &
                   sum(dPdF(ii,1:3,jj,1:3,active(grIJ))* &
                       math_tensorproduct33(param(instance)%interfaceNormal(xioffsetIN+1_pInt:xioffsetIN+3_pInt,offset), &
                                            param(instance)%interfaceNormal(xioffsetJN+1_pInt:xioffsetJN+3_pInt,offset)))                            
           if (grIJ == grJJ) &
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                 
                 jacobian(xioffsetI+ii,xioffsetJ+jj) = &
                   jacobian(xioffsetI+ii,xioffsetJ+jj) +  &
                   phasefrac(homog)%p(active(grII),offset)* &
                   phasefrac(homog)%p(active(grIJ),offset)* &
                   phasefrac(homog)%p(active(grJI),offset)* &
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
     else convergencePartialRankOne                                                                 ! not converged and residuum not improved...
       homogenization_multiphase_updateState = [.false., .true.]
       param(instance)%stepLength(offset) = &
         param(instance)%stepLength(offset)/2.0_pReal                                               ! ...try with smaller step length in same direction
       param(instance)%newIter  (:,offset) = &
         param(instance)%oldIter  (:,offset) - &
         param(instance)%stepLength(offset)* &
         param(instance)%searchDir(:,offset)
     endif convergencePartialRankOne

   case(fullrankone_ID) myMixRule
     Nactive = 0_pInt
     do grI = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(homog)%p(grI,offset) > err_phasefr_tolabs) Nactive = Nactive + 1_pInt
     enddo
     allocate(active(Nactive), source = 0_pInt)
    
     avgR = 0.0_pReal
     Nactive = 0_pInt
     do grI = 1_pInt, homogenization_Ngrains(homog)
       if (phasefrac(homog)%p(grI,offset) > err_phasefr_tolabs) then
         Nactive = Nactive + 1_pInt 
         active(Nactive) = grI  
         avgR = avgR + phasefrac(homog)%p(grI,offset)*P(1:3,1:3,grI)
       endif
     enddo     
     
     xisize = 5_pInt*Nactive*(Nactive-1_pInt)/2_pInt
     allocate(residual(xisize       ), source=0.0_pReal)
     allocate(jacobian(xisize,xisize), source=0.0_pReal)      
     do grI = 1_pInt, Nactive-1_pInt
       do grJ = grI+1_pInt, Nactive
         xioffsetIN = 5_pInt*((active(grI)-1_pInt)* &
                              (2_pInt*homogenization_Ngrains(homog)-active(grI))/2_pInt + &
                              active(grJ) - active(grI) - 1_pInt)
         xioffsetI  = 5_pInt*((grI-1_pInt)*(2_pInt*Nactive-grI)/2_pInt + grJ - grI - 1_pInt)
         interfaceNormalAIJ = param(instance)%newIter(xioffsetIN+1,offset)
         interfaceNormalBIJ = param(instance)%newIter(xioffsetIN+2,offset)
         residual(xioffsetI+1_pInt) = &                                                             ! RA= residual for interface normal alpha
           phasefrac(homog)%p(active(grI),offset)* &
           phasefrac(homog)%p(active(grJ),offset)* &
           math_mul33xx33( &
               transpose(P(1:3,1:3,active(grI)) - P(1:3,1:3,active(grJ))), &
               math_tensorproduct33( &
                                    [ cos(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                      cos(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                     -sin(interfaceNormalAIJ)], &
                                    param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset) &
                                   ) & 
                         )                
         residual(xioffsetI+2_pInt) = &                                                             ! RB= residual for interface normal beta
           phasefrac(homog)%p(active(grI),offset)* &
           phasefrac(homog)%p(active(grJ),offset)* &
           math_mul33xx33( &
               transpose(P(1:3,1:3,active(grI)) - P(1:3,1:3,active(grJ))), &
               math_tensorproduct33( &
                                    [-sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                      sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                      0.0_preal], &
                                    param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset) &
                                   ) &
                         )
         residual(xioffsetI+3_pInt:xioffsetI+5_pInt) = &                                            ! RH= residual for interface jump
           phasefrac(homog)%p(active(grI),offset)* &
           phasefrac(homog)%p(active(grJ),offset)* &
           math_mul33x3(P(1:3,1:3,active(grI)) - P(1:3,1:3,active(grJ)), &
                        [sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                         sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                         cos(interfaceNormalAIJ) ] )                      
       enddo
     enddo 
     
     stressTol = max(            param(instance)%absTol, &
                     norm2(avgR)*param(instance)%relTol)
     convergenceFullRankOne: if (     norm2(residual) < stressTol &
                                 .or. Nactive         < 2_pInt) then
       homogenization_multiphase_updateState = [.true., .true.]  
       
     elseif (     iter == 1_pInt &
             .or. norm2(residual) < param(instance)%residual(offset)) then convergenceFullRankOne   ! not converged, but improved norm of residuum (always proceed in first iteration)...
       homogenization_multiphase_updateState = [.false., .true.]
       
       param(instance)%residual  (offset) = norm2(residual)                                         ! ...remember old values and...
       param(instance)%stepLength(offset) = &
        min(1.0_pReal,2.0_pReal*param(instance)%stepLength(offset))                                 ! ...proceed with normal step length (calculate new search direction)
       param(instance)%oldIter(:,offset) = param(instance)%newIter(:,offset)
                
       do grII = 1_pInt,  Nactive-1; do grIJ = grII+1_pInt, Nactive
         xioffsetIN = 5_pInt*((active(grII)-1_pInt)* &
                              (2_pInt*homogenization_Ngrains(homog)-active(grII))/2_pInt + &
                              active(grIJ) - active(grII) - 1_pInt)
         xioffsetI  = 5_pInt*((grII-1_pInt)*(2_pInt*Nactive-grII)/2_pInt + grIJ - grII - 1_pInt)
         interfaceNormalAIJ = param(instance)%newIter(xioffsetIN+1,offset)
         interfaceNormalBIJ = param(instance)%newIter(xioffsetIN+2,offset)
         jacobian(xioffsetI+1,xioffsetI+1) = &                                                      ! jacobian (alpha,alpha) RA_A  
           jacobian(xioffsetI+1,xioffsetI+1) +  &
           phasefrac(homog)%p(active(grII),offset)* &
           phasefrac(homog)%p(active(grIJ),offset)* &
           math_mul33xx33(P(1:3,1:3,active(grII)) - P(1:3,1:3,active(grIJ)), &
                          math_tensorproduct33( &
                            param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset) , &
                            [-sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                             -sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                             -cos(interfaceNormalAIJ)] &
                                              ))
         jacobian(xioffsetI+1,xioffsetI+2) = &        !                                             ! jacobian (alpha,beta) RA_B
           jacobian(xioffsetI+1,xioffsetI+2) +  &
           phasefrac(homog)%p(active(grII),offset)* &
           phasefrac(homog)%p(active(grIJ),offset)* &
           math_mul33xx33(P(1:3,1:3,active(grII)) - P(1:3,1:3,active(grIJ)), &
                          math_tensorproduct33( &
                            param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset) , &
                            [-cos(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                              cos(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                              0.0_pReal]      ))
         forall (ii = 1_pInt:3_pInt) &                                                                                                
           jacobian(xioffsetI+1,xioffsetI+2+ii) = &                                                 ! jacobian (alpha,H) RA_H   
             jacobian(xioffsetI+1,xioffsetI+2+ii) +  &
             phasefrac(homog)%p(active(grII),offset)* &
             phasefrac(homog)%p(active(grIJ),offset)* &
             sum((P(ii,1:3,active(grII)) - P(ii,1:3,active(grIJ)))* &
                 [ cos(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                   cos(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                  -sin(interfaceNormalAIJ)] )
         jacobian(xioffsetI+2,xioffsetI+1) = &                                                      ! jacobian (beta,alpha) RB_A  
           jacobian(xioffsetI+2,xioffsetI+1) +  &
           phasefrac(homog)%p(active(grII),offset)* &
           phasefrac(homog)%p(active(grIJ),offset)* &
           math_mul33xx33(P(1:3,1:3,active(grII)) - P(1:3,1:3,active(grIJ)), &
                          math_tensorproduct33( &
                            param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                            [-cos(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                              cos(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                              0.0_pReal]      ))
         jacobian(xioffsetI+2,xioffsetI+2) = &                                                      ! jacobian (beta,beta) RB_B 
           jacobian(xioffsetI+2,xioffsetI+2) +  &
           phasefrac(homog)%p(active(grII),offset)* &
           phasefrac(homog)%p(active(grIJ),offset)* &
           math_mul33xx33(P(1:3,1:3,active(grII)) - P(1:3,1:3,active(grIJ)), &
                          math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset) , &
                                               [-sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                                -sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                                 0.0_pReal]) &
                          )
         forall (ii = 1_pInt:3_pInt) &                                                                                                
           jacobian(xioffsetI+2,xioffsetI+2+ii) = &                                                 ! jacobian (beta,H) RB_H   
             jacobian(xioffsetI+2,xioffsetI+2+ii) +  &
             phasefrac(homog)%p(active(grII),offset)* &
             phasefrac(homog)%p(active(grIJ),offset)* &
             sum((P(ii,1:3,active(grII)) - P(ii,1:3,active(grIJ)))* &
                 [-sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                  sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                  0.0_pReal])                                  
         forall (ii = 1_pInt:3_pInt) &                                                                                                
           jacobian(xioffsetI+2+ii,xioffsetI+1) = &                                                  ! jacobian (H,alpha) RH_A 
             jacobian(xioffsetI+2+ii,xioffsetI+1) +  &
             phasefrac(homog)%p(active(grII),offset)* &
             phasefrac(homog)%p(active(grIJ),offset)* &
             sum((P(ii,1:3,active(grII)) - P(ii,1:3,active(grIJ)))* &
                 [ cos(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                   cos(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                  -sin(interfaceNormalAIJ)] )                        
         forall (ii = 1_pInt:3_pInt) &                                                                                                
           jacobian(xioffsetI+2+ii,xioffsetI+2) = &                                                  ! jacobian (H,beta) RH_B  
             jacobian(xioffsetI+2+ii,xioffsetI+2) +  &
             phasefrac(homog)%p(active(grII),offset)* &
             phasefrac(homog)%p(active(grIJ),offset)* &
             sum(  (P(ii,1:3,active(grII)) - P(ii,1:3,active(grIJ)) )* &
                   [-sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                     sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                     0.0_pReal] )                                              
         do grJI = 1_pInt,  Nactive-1; do grJJ = grJI+1_pInt, Nactive
           xioffsetJN = 5_pInt*((active(grJI)-1_pInt)* &
                                (2_pInt*homogenization_Ngrains(homog)-active(grJI))/2_pInt + &
                                active(grJJ) - active(grJI) - 1_pInt)
           xioffsetJ  = 5_pInt*((grJI-1_pInt)*(2_pInt*Nactive-grJI)/2_pInt + grJJ - grJI - 1_pInt)
           interfaceNormalAKL = param(instance)%newIter(xioffsetJN+1,offset)
           interfaceNormalBKL = param(instance)%newIter(xioffsetJN+2,offset)                  
           if (grII == grJI) then                                                                   ! terms when I=K        
             jacobian(xioffsetI+1,xioffsetJ+1) = &                                                  ! jacobian (alpha,alpha) RA_A (1,1)
               jacobian(xioffsetI+1,xioffsetJ+1) +  &                        
               phasefrac(homog)%p(active(grII),offset)* &
               phasefrac(homog)%p(active(grIJ),offset)* &
               phasefrac(homog)%p(active(grJJ),offset)* &
               sum(dPdF(1:3,1:3,1:3,1:3,active(grII))* &
                   math_tensorproduct3333( &
                          math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                               [ cos(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                                 cos(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                                 -sin(interfaceNormalAIJ)] ), & 
                          math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset), &
                                               [ cos(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                                 cos(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                                 -sin(interfaceNormalAKL)] ) &
                                         )   &
                  )         
             jacobian(xioffsetI+1,xioffsetJ+2) = &                                                  ! jacobian (alpha,beta) RA_B  (1,2)
               jacobian(xioffsetI+1,xioffsetJ+2) +  &                        
               phasefrac(homog)%p(active(grII),offset)* &
               phasefrac(homog)%p(active(grIJ),offset)* &
               phasefrac(homog)%p(active(grJJ),offset)* &
               sum( dPdF(1:3,1:3,1:3,1:3,active(grII))* &
                    math_tensorproduct3333( &
                         math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                              [cos(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                               cos(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                            -sin(interfaceNormalAIJ)] ), & 
                          math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset), &
                                               [-sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                                sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                                0.0_pReal]              ) &
                                          )   &
                  )                                 
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                                       ! jacobian (alpha,H) RA_H  (1,3:5)         
               avgR(ii,jj) =sum(dPdF(1:3,1:3,ii,jj,active(grII))* &
                                math_tensorproduct33( &
                                   param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                   [cos(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                    cos(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                   -sin(interfaceNormalAIJ)] &
                                                    ))
             forall (ii = 1_pInt:3_pInt) &                                                              
                 jacobian(xioffsetI+1,xioffsetJ+2+ii) = &                                             
                   jacobian(xioffsetI+1,xioffsetJ+2+ii) +  &                        
                   phasefrac(homog)%p(active(grII),offset)* &
                   phasefrac(homog)%p(active(grIJ),offset)* &
                   phasefrac(homog)%p(active(grJJ),offset)* &
                   sum(avgR(ii,1:3)*[sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                     sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                     cos(interfaceNormalAKL)] ) 
             jacobian(xioffsetI+2,xioffsetJ+1) = &                                                  ! jacobian (beta,alpha) RB_A  (2,1)
               jacobian(xioffsetI+2,xioffsetJ+1) +  &                        
               phasefrac(homog)%p(active(grII),offset)* &
               phasefrac(homog)%p(active(grIJ),offset)* &
               phasefrac(homog)%p(active(grJJ),offset)* &
               sum(dPdF(1:3,1:3,1:3,1:3,active(grII))* &
                   math_tensorproduct3333( &
                         math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                              [-sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                                sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                                0.0_Preal] ), & 
                         math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset), &
                                              [cos(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                               cos(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                               -sin(interfaceNormalAKL)] ) &
                                         )   &
                  )
             jacobian(xioffsetI+2,xioffsetJ+2) = &                                                  ! jacobian (beta,beta) RB_B (2,2)
               jacobian(xioffsetI+2,xioffsetJ+2) +  &                        
               phasefrac(homog)%p(active(grII),offset)* &
               phasefrac(homog)%p(active(grIJ),offset)* &
               phasefrac(homog)%p(active(grJJ),offset)* &
               sum(dPdF(1:3,1:3,1:3,1:3,active(grII))* &
                   math_tensorproduct3333( &
                        math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                             [-sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                              sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                              0.0_Preal] ), & 
                        math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset), &
                                             [-sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                             sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                             0.0_pReal] ) &
                                         )   &
                  )                                                                         
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                                        ! jacobian (beta,H) RB_H (2,3:5)
               avgR(ii,jj) = sum(dPdF(1:3,1:3,ii,jj,active(grII))* &
                                 math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                                      [-sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                                        sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                                        0.0_Preal] ) )
             forall (ii = 1_pInt:3_pInt) &                                                              
               jacobian(xioffsetI+2,xioffsetJ+2+ii) = &          
                 jacobian(xioffsetI+2,xioffsetJ+2+ii) +  &                        
                 phasefrac(homog)%p(active(grII),offset)* &
                 phasefrac(homog)%p(active(grIJ),offset)* &
                 phasefrac(homog)%p(active(grJJ),offset)* &
                 sum(avgR(ii,1:3)*[sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                   sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                   cos(interfaceNormalAKL)])                                                                 
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                                       ! jacobian (H,alpha) RH1_A (3:5,1) 
               avgR(ii,jj) = sum(dPdF(ii,jj,1:3,1:3,active(grII))* &
                                 math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset), &
                                                      [ cos(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                                        cos(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                                       -sin(interfaceNormalAKL)] ))
             forall (ii = 1_pInt:3_pInt) &                                                              
               jacobian(xioffsetI+2+ii,xioffsetJ+1) = &                                          
                 jacobian(xioffsetI+2+ii,xioffsetJ+1) +  &                        
                 phasefrac(homog)%p(active(grII),offset)* &
                 phasefrac(homog)%p(active(grIJ),offset)* &
                 phasefrac(homog)%p(active(grJJ),offset)* &
                 sum(avgR(ii,1:3)*[sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                   sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                   cos(interfaceNormalAIJ)]) 
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                                                           
               avgR(ii,jj) = sum(dPdF(ii,jj,1:3,1:3,active(grII))* &
                                 math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset) , &
                                                      [-sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                                        sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                                        0.0_pReal] ) )
             forall (ii = 1_pInt:3_pInt) &                                                              
               jacobian(xioffsetI+2+ii,xioffsetJ+2) = &                                             ! jacobian (H,beta) RH2_B (3:5,2)
                 jacobian(xioffsetI+2+ii,xioffsetJ+2) +  &                        
                 phasefrac(homog)%p(active(grII),offset)* &
                 phasefrac(homog)%p(active(grIJ),offset)* &
                 phasefrac(homog)%p(active(grJJ),offset)* &
                 sum(avgR(ii,1:3)*[sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                   sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                   cos(interfaceNormalAIJ)]) 
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                                       ! jacobian (H,H) RH_H (3:5,3:5)
               jacobian(xioffsetI+2+ii,xioffsetJ+2+jj) = &
                 jacobian(xioffsetI+2+ii,xioffsetJ+2+jj) +  &
                 phasefrac(homog)%p(active(grII),offset)* &
                 phasefrac(homog)%p(active(grIJ),offset)* &
                 phasefrac(homog)%p(active(grJJ),offset)* &
                 sum(dPdF(ii,1:3,jj,1:3,active(grII))* &
                     math_tensorproduct33( &
                                [sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                cos(interfaceNormalAIJ)], &
                                [sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                cos(interfaceNormalAKL)] &
                                         ) &
                    )                                                                   
           endif
!-------------------------------------------------------------------------------------
           if (grII == grJJ) then                                                                   ! subtracting contribution from L 
             jacobian(xioffsetI+1,xioffsetJ+1) = &                                                  ! jacobian (A,A) RA_A (1,1)
               jacobian(xioffsetI+1,xioffsetJ+1) -  &                        
               phasefrac(homog)%p(active(grII),offset)* &
               phasefrac(homog)%p(active(grIJ),offset)* &
               phasefrac(homog)%p(active(grJI),offset)* &
               sum(dPdF(1:3,1:3,1:3,1:3,active(grII))* &
                   math_tensorproduct3333( &
                     math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                          [ cos(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                            cos(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                           -sin(interfaceNormalAIJ)]), & 
                     math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset), &
                                          [ cos(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                            cos(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                           -sin(interfaceNormalAKL)] ) & 
                                         )   &
                  )
             jacobian(xioffsetI+1,xioffsetJ+2) = &                                                  ! jacobian (A,B) RA_B (1,2)
               jacobian(xioffsetI+1,xioffsetJ+2)  -  &                        
               phasefrac(homog)%p(active(grII),offset)* &
               phasefrac(homog)%p(active(grIJ),offset)* &
               phasefrac(homog)%p(active(grJI),offset)* &
               sum(dPdF(1:3,1:3,1:3,1:3,active(grII))* &
                   math_tensorproduct3333( &
                     math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                          [ cos(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                            cos(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                           -sin(interfaceNormalAIJ)] ), & 
                     math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset), &
                                          [-sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                            sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                            0.0_pReal] ) &
                                         )   &
                  )
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                                       ! jacobian (A,H) RA_H (1,3:5)
               avgR(ii,jj) = sum(dPdF(1:3,1:3,ii,jj,active(grII))* &
                                 math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                                      [ cos(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                                        cos(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                                       -sin(interfaceNormalAIJ)]) )
             forall (ii = 1_pInt:3_pInt) &                                                              
                 jacobian(xioffsetI+1,xioffsetJ+2+ii) = &        
                   jacobian(xioffsetI+1,xioffsetJ+2+ii) -  &                        
                   phasefrac(homog)%p(active(grII),offset)* &
                   phasefrac(homog)%p(active(grIJ),offset)* &
                   phasefrac(homog)%p(active(grJI),offset)* &
                   sum(avgR(ii,1:3)*[sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                     sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                     cos(interfaceNormalAKL)])                                    
             jacobian(xioffsetI+2,xioffsetJ+1) = &                                                  ! jacobian (B,A) RB_A (2,1)
               jacobian(xioffsetI+2,xioffsetJ+1) -  &                        
               phasefrac(homog)%p(active(grII),offset)* &
               phasefrac(homog)%p(active(grIJ),offset)* &
               phasefrac(homog)%p(active(grJI),offset)* &
               sum(dPdF(1:3,1:3,1:3,1:3,active(grII))* &
                   math_tensorproduct3333( &
                        math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                             [-sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                               sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                               0.0_Preal]), & 
                        math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset), &
                                             [ cos(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                               cos(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                              -sin(interfaceNormalAKL)] ) &
                                         )   &
                  )
             jacobian(xioffsetI+2,xioffsetJ+2) = &                                                  ! jacobian (B,B) RB_B (2,2)
               jacobian(xioffsetI+2,xioffsetJ+2) -  &                        
               phasefrac(homog)%p(active(grII),offset)* &
               phasefrac(homog)%p(active(grIJ),offset)* &
               phasefrac(homog)%p(active(grJI),offset)* &
               sum(dPdF(1:3,1:3,1:3,1:3,active(grII))* &
                   math_tensorproduct3333( &
                       math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                            [-sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                              sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                              0.0_Preal] ), & 
                       math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset), &
                                            [-sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                              sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                              0.0_pReal] ) &
                                        )   &
                 )
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                                       ! jacobian (B,H) RB_H (2,3:5)
               avgR(ii,jj) = sum(dPdF(1:3,1:3,ii,jj,active(grII))* &
                                 math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                                      [-sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                                        sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                                        0.0_Preal]) )
             forall (ii = 1_pInt:3_pInt) &                                                              
               jacobian(xioffsetI+2,xioffsetJ+2+ii) = &                                       
                 jacobian(xioffsetI+2,xioffsetJ+2+ii) - &                        
                 phasefrac(homog)%p(active(grII),offset)* &
                 phasefrac(homog)%p(active(grIJ),offset)* &
                 phasefrac(homog)%p(active(grJI),offset)* &
                 sum(avgR(ii,1:3)*[sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                   sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                   cos(interfaceNormalAKL)])                                  
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                                       ! jacobian (H,A) RH_A (3:5,1)   
               avgR(ii,jj) = sum(dPdF(ii,jj,1:3,1:3,active(grII))* &
                                 math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset), &
                                                      [cos(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                                       cos(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                                       sin(interfaceNormalAKL)])  )
             forall (ii = 1_pInt:3_pInt) &                                                              
               jacobian(xioffsetI+2+ii,xioffsetJ+1) = &        
                 jacobian(xioffsetI+2+ii,xioffsetJ+1) -  &                        
                 phasefrac(homog)%p(active(grII),offset)* &
                 phasefrac(homog)%p(active(grIJ),offset)* &
                 phasefrac(homog)%p(active(grJI),offset)* &
                 sum(avgR(ii,1:3)*[sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                   sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                   cos(interfaceNormalAIJ)]) 
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                                       ! jacobian (H,B) RH_B (3:5,2) 
               avgR(ii,jj) = sum(dPdF(ii,jj,1:3,1:3,active(grII))* &
                                 math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset) , &
                                                      [-sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                                        sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                                        0.0_pReal] ) )
             forall (ii = 1_pInt:3_pInt) &                                                                         
               jacobian(xioffsetI+2+ii,xioffsetJ+2) = &                          
                 jacobian(xioffsetI+2+ii,xioffsetJ+2) -  &                        
                 phasefrac(homog)%p(active(grII),offset)* &
                 phasefrac(homog)%p(active(grIJ),offset)* &
                 phasefrac(homog)%p(active(grJI),offset)* &
                 sum(avgR(ii,1:3)*[sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                   sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                   cos(interfaceNormalAIJ)]) 
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                                       ! jacobian (H,H) RH_H (3:5,3:5)
               jacobian(xioffsetI+2+ii,xioffsetJ+2+jj) = &
                 jacobian(xioffsetI+2+ii,xioffsetJ+2+jj) -  &
                 phasefrac(homog)%p(active(grII),offset)* &
                 phasefrac(homog)%p(active(grIJ),offset)* &
                 phasefrac(homog)%p(active(grJI),offset)* &
                 sum(dPdF(ii,1:3,jj,1:3,active(grII))* &
                     math_tensorproduct33([sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                           sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                           cos(interfaceNormalAIJ)], &
                                          [sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                           sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                           cos(interfaceNormalAKL)])  )                                                                   
           endif
           if (grIJ == grJI) then                                                                   ! KL terms 
             jacobian(xioffsetI+1,xioffsetJ+1) = &                                                  ! jacobian (A,A) RA_A (1,1)
               jacobian(xioffsetI+1,xioffsetJ+1) -  &                        
               phasefrac(homog)%p(active(grII),offset)* &
               phasefrac(homog)%p(active(grIJ),offset)* &
               phasefrac(homog)%p(active(grJJ),offset)* &
               sum(dPdF(1:3,1:3,1:3,1:3,active(grIJ))* &
                   math_tensorproduct3333( &
                       math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                            [ cos(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                              cos(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                             -sin(interfaceNormalAIJ)]), & 
                       math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset), &
                                            [ cos(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                              cos(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                             -sin(interfaceNormalAKL)]) &
                                         )   &
                   )
             jacobian(xioffsetI+1,xioffsetJ+2) = &                                                  ! jacobian (A,B) RA_B (1,2)
               jacobian(xioffsetI+1,xioffsetJ+2) -  &                        
               phasefrac(homog)%p(active(grII),offset)* &
               phasefrac(homog)%p(active(grIJ),offset)* &
               phasefrac(homog)%p(active(grJJ),offset)* &
               sum(dPdF(1:3,1:3,1:3,1:3,active(grIJ))* &
                   math_tensorproduct3333( &
                      math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                           [ cos(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                             cos(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                            -sin(interfaceNormalAIJ)] ), & 
                      math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset), &
                                           [-sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                             sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                             0.0_pReal]) &
                                         )   &
                  )                                          
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                                       ! jacobian (A,H) RA_H (1,3:5)
               avgR(ii,jj) = sum(dPdF(1:3,1:3,ii,jj,active(grIJ))* &
                                 math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                                      [ cos(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                                        cos(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                                       -sin(interfaceNormalAIJ)]) )
             forall (ii = 1_pInt:3_pInt) &                                                              
               jacobian(xioffsetI+1,xioffsetJ+2+ii) = &        !  
                 jacobian(xioffsetI+1,xioffsetJ+2+ii) -  &                         
                 phasefrac(homog)%p(active(grII),offset)* &
                 phasefrac(homog)%p(active(grIJ),offset)* &
                 phasefrac(homog)%p(active(grJJ),offset)* &
                 sum(avgR(ii,1:3)*[sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                   sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                   cos(interfaceNormalAKL)])                   
             jacobian(xioffsetI+2,xioffsetJ+1) = &                                                  ! jacobian (B,A) RB_A (2,1)
               jacobian(xioffsetI+2,xioffsetJ+1) -  &                       
               phasefrac(homog)%p(active(grII),offset)* &
               phasefrac(homog)%p(active(grIJ),offset)* &
               phasefrac(homog)%p(active(grJJ),offset)* &
               sum(dPdF(1:3,1:3,1:3,1:3,active(grIJ))* &
                   math_tensorproduct3333( &
                       math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                           [-sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                             sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                             0.0_Preal] ), & 
                       math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset), &
                                           [ cos(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                             cos(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                            -sin(interfaceNormalAKL)] ) &
                                          )   &
                  )
             jacobian(xioffsetI+2,xioffsetJ+2) = &                                                  ! jacobian (B,B) RB_B (2,2)
               jacobian(xioffsetI+2,xioffsetJ+2) -  &                        
               phasefrac(homog)%p(active(grII),offset)* &
               phasefrac(homog)%p(active(grIJ),offset)* &
               phasefrac(homog)%p(active(grJJ),offset)* &
               sum(dPdF(1:3,1:3,1:3,1:3,active(grIJ))* &
                   math_tensorproduct3333( &
                             math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                                 [-sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                                   sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                                   0.0_Preal] ), & 
                             math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset), &
                                                 [-sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                                   sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                                   0.0_pReal] ) &
                                         )   &
                  )
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                                       ! jacobian (B,H) RB_H (2,3:5)
               avgR(ii,jj) = sum(dPdF(1:3,1:3,ii,jj,active(grIJ))* &
                                 math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                                      [-sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                                        sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                                        0.0_Preal] ) )
             forall (ii = 1_pInt:3_pInt) &                                                              
               jacobian(xioffsetI+2,xioffsetJ+2+ii) = &        
                              jacobian(xioffsetI+2,xioffsetJ+2+ii) -  &                       
                              phasefrac(homog)%p(active(grII),offset)* &
                              phasefrac(homog)%p(active(grIJ),offset)* &
                              phasefrac(homog)%p(active(grJJ),offset)* &
                              sum(avgR(ii,1:3)*[sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                                sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                                cos(interfaceNormalAKL)])                                  
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                                       ! jacobian (H,A) RH_A (3:5,1)
               avgR(ii,jj) = sum(dPdF(ii,jj,1:3,1:3,active(grIJ))* &
                                 math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset), &
                                                      [ cos(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                                        cos(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                                       -sin(interfaceNormalAKL)]) )
             forall (ii = 1_pInt:3_pInt) &                                                              
               jacobian(xioffsetI+2+ii,xioffsetJ+1) = &        !  
                 jacobian(xioffsetI+2+ii,xioffsetJ+1) -  &                         
                 phasefrac(homog)%p(active(grII),offset)* &
                 phasefrac(homog)%p(active(grIJ),offset)* &
                 phasefrac(homog)%p(active(grJJ),offset)* &
                 sum(avgR(ii,1:3)*[sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                   sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                   cos(interfaceNormalAIJ)]) 
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                                       ! jacobian (H,B) RH_B (3:5,2)     
               avgR(ii,jj) = sum(dPdF(ii,jj,1:3,1:3,active(grIJ))* &
                                 math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset) , &
                                                      [-sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                                        sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                                        0.0_pReal] ) )
             forall (ii = 1_pInt:3_pInt) &                                                              
               jacobian(xioffsetI+2+ii,xioffsetJ+2) = &        !  
                 jacobian(xioffsetI+2+ii,xioffsetJ+2) -  &                         
                 phasefrac(homog)%p(active(grII),offset)* &
                 phasefrac(homog)%p(active(grIJ),offset)* &
                 phasefrac(homog)%p(active(grJJ),offset)* &
                 sum(avgR(ii,1:3)*[ sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                    sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                    cos(interfaceNormalAIJ)]) 
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                 
               jacobian(xioffsetI+2+ii,xioffsetJ+2+jj) = &                                          ! jacobian (H,H) RH_H (3:5,3:5)
                 jacobian(xioffsetI+2+ii,xioffsetJ+2+jj) -  &
                 phasefrac(homog)%p(active(grII),offset)* &
                 phasefrac(homog)%p(active(grIJ),offset)* &
                 phasefrac(homog)%p(active(grJJ),offset)* &
                 sum(dPdF(ii,1:3,jj,1:3,active(grIJ))* &
                     math_tensorproduct33([ sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                            sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                            cos(interfaceNormalAIJ)], &
                                          [ sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                            sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                            cos(interfaceNormalAKL)] ) ) 
           endif
           if (grIJ == grJJ) then                                                                   ! adding contributions from JJ
             jacobian(xioffsetI+1,xioffsetJ+1) = &                                                  ! jacobian (A,A) RA_A (1,1)
               jacobian(xioffsetI+1,xioffsetJ+1) +  &                        
               phasefrac(homog)%p(active(grII),offset)* &
               phasefrac(homog)%p(active(grIJ),offset)* &
               phasefrac(homog)%p(active(grJI),offset)* &
               sum(dPdF(1:3,1:3,1:3,1:3,active(grIJ))* &
                   math_tensorproduct3333( &
                      math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                           [ cos(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                             cos(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                            -sin(interfaceNormalAIJ)] ), & 
                      math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset), &
                                           [ cos(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                             cos(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                            -sin(interfaceNormalAKL)] ) &
                                          )   &
                  )
             jacobian(xioffsetI+1,xioffsetJ+2) = &                                                  ! jacobian (A,B) RA_B (1,2)
               jacobian(xioffsetI+1,xioffsetJ+2) +  &                        
               phasefrac(homog)%p(active(grII),offset)* &
               phasefrac(homog)%p(active(grIJ),offset)* &
               phasefrac(homog)%p(active(grJI),offset)* &
               sum(dPdF(1:3,1:3,1:3,1:3,active(grIJ))* &
                   math_tensorproduct3333( &
                         math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                              [ cos(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                                cos(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                               -sin(interfaceNormalAIJ)] ), & 
                         math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset), &
                                              [-sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                                sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                                0.0_pReal] ) &
                                              )   &
                       )                                                 
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                                       ! jacobian (A,H) RA_H (1,3:5)
               avgR(ii,jj) = sum(dPdF(1:3,1:3,ii,jj,active(grIJ))* &
                                 math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                                      [ cos(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                                        cos(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                                       -sin(interfaceNormalAIJ)] ))
             forall (ii = 1_pInt:3_pInt) &                                                              
               jacobian(xioffsetI+1,xioffsetJ+2+ii) = &        !  
                 jacobian(xioffsetI+1,xioffsetJ+2+ii) +  &                         
                 phasefrac(homog)%p(active(grII),offset)* &
                 phasefrac(homog)%p(active(grIJ),offset)* &
                 phasefrac(homog)%p(active(grJI),offset)* &
                 sum(avgR(ii,1:3)*[sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                   sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                   cos(interfaceNormalAKL)])             
             jacobian(xioffsetI+2,xioffsetJ+1) = &                                                  ! jacobian (B,A) RB_A (2,1)
               jacobian(xioffsetI+2,xioffsetJ+1) +  &                       
               phasefrac(homog)%p(active(grII),offset)* &
               phasefrac(homog)%p(active(grIJ),offset)* &
               phasefrac(homog)%p(active(grJI),offset)* &
               sum(dPdF(1:3,1:3,1:3,1:3,active(grIJ))* &
                   math_tensorproduct3333( &
                       math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                            [-sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                              sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                              0.0_Preal] ), & 
                       math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset), &
                                            [ cos(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                              cos(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                             -sin(interfaceNormalAKL)] )  &
                                         )   &
                  )                                            
             jacobian(xioffsetI+2,xioffsetJ+2) = &                                                  ! jacobian (B,B) RB_B (2,2)
               jacobian(xioffsetI+2,xioffsetJ+2) +  &                        
               phasefrac(homog)%p(active(grII),offset)* &
               phasefrac(homog)%p(active(grIJ),offset)* &
               phasefrac(homog)%p(active(grJI),offset)* &
               sum(dPdF(1:3,1:3,1:3,1:3,active(grIJ))* &
                   math_tensorproduct3333( &
                       math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                            [-sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                              sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                              0.0_Preal] ), & 
                       math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset), &
                                           [-sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                             sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                             0.0_pReal] ) &
                                         )   &
                  )                
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                                       ! jacobian (B,H) RB_H (2,3:5)
               avgR(ii,jj) = sum(dPdF(1:3,1:3,ii,jj,active(grIJ))* &
                                 math_tensorproduct33(param(instance)%newIter(xioffsetIN+3:xioffsetIN+5,offset), &
                                                      [-sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                                        sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                                        0.0_Preal] ) )
             forall (ii = 1_pInt:3_pInt) &                                                              
               jacobian(xioffsetI+2,xioffsetJ+2+ii) = &         
                 jacobian(xioffsetI+2,xioffsetJ+2+ii) +  &                       
                 phasefrac(homog)%p(active(grII),offset)* &
                 phasefrac(homog)%p(active(grIJ),offset)* &
                 phasefrac(homog)%p(active(grJI),offset)* &
                 sum(avgR(ii,1:3)*[sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                   sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                   cos(interfaceNormalAKL)])                                  
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                                       ! jacobian (H,A) RH_A (3:5,1)  
               avgR(ii,jj) = sum(dPdF(ii,jj,1:3,1:3,active(grIJ))* &
                                 math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset), &
                                                      [ cos(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                                        cos(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                                       -sin(interfaceNormalAKL)] ))
             forall (ii = 1_pInt:3_pInt) &                                                              
               jacobian(xioffsetI+2+ii,xioffsetJ+1) = &          
                 jacobian(xioffsetI+2+ii,xioffsetJ+1) +  &                         
                 phasefrac(homog)%p(active(grII),offset)* &
                 phasefrac(homog)%p(active(grIJ),offset)* &
                 phasefrac(homog)%p(active(grJI),offset)* &
                 sum(avgR(ii,1:3)*[sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                   sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                   cos(interfaceNormalAIJ)]) 
             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                           
               avgR(ii,jj) = sum(dPdF(ii,jj,1:3,1:3,active(grIJ))* &
                             math_tensorproduct33(param(instance)%newIter(xioffsetJN+3:xioffsetJN+5,offset) , &
                                                  [-sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                                    sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                                    0.0_pReal] ) )
             forall (ii = 1_pInt:3_pInt) &                                                              
               jacobian(xioffsetI+2+ii,xioffsetJ+2) = &                                         ! jacobian (H,B) RH_B (3:5,2)  
                 jacobian(xioffsetI+2+ii,xioffsetJ+2) +  &                         
                 phasefrac(homog)%p(active(grII),offset)* &
                 phasefrac(homog)%p(active(grIJ),offset)* &
                 phasefrac(homog)%p(active(grJI),offset)* &
                 sum(avgR(ii,1:3)*[sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                   sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                   cos(interfaceNormalAIJ)]) 

             forall (ii = 1_pInt:3_pInt,jj = 1_pInt:3_pInt) &                 
               jacobian(xioffsetI+2+ii,xioffsetJ+2+jj) = &                                        ! jacobian (H,H) RH_H (3:5,3:5)
                 jacobian(xioffsetI+2+ii,xioffsetJ+2+jj) +  &
                 phasefrac(homog)%p(active(grII),offset)* &
                 phasefrac(homog)%p(active(grIJ),offset)* &
                 phasefrac(homog)%p(active(grJI),offset)* &
                 sum(dPdF(ii,1:3,jj,1:3,active(grIJ))* &
                     math_tensorproduct33([sin(interfaceNormalAIJ)*cos(interfaceNormalBIJ), &
                                           sin(interfaceNormalAIJ)*sin(interfaceNormalBIJ), &
                                           cos(interfaceNormalAIJ)], &
                                          [sin(interfaceNormalAKL)*cos(interfaceNormalBKL), &
                                           sin(interfaceNormalAKL)*sin(interfaceNormalBKL), &
                                           cos(interfaceNormalAKL)] ) )                                                                             
           endif
         enddo; enddo
       enddo; enddo  
       jacobian = jacobian + stressTol*math_identity2nd(xisize)                                     ! perturbation by identity so that jacobian is invertible
       allocate(ipiv(xisize))
       call dgesv(xisize,1,jacobian,xisize,ipiv,residual,xisize,ierr)                               !< solve Jacobian * delta state = -residual for delta state
       if (ierr == 0_pInt) then
         param(instance)%searchDir(:,offset) = 0.0_pReal
         do grI = 1_pInt, Nactive-1_pInt
           do grJ = grI+1_pInt, Nactive
             xioffsetIN = 5_pInt*((active(grI)-1_pInt)*(2_pInt*homogenization_Ngrains(homog)-active(grI))/2_pInt + &
                                  active(grJ) - active(grI) - 1_pInt)
             xioffsetI  = 5_pInt*((grI-1_pInt)*(2_pInt*Nactive-grI)/2_pInt + grJ - grI - 1_pInt) 
             param(instance)%searchDir(xioffsetIN+1_pInt:xioffsetIN+5_pInt,offset) = &
               residual(xioffsetI+1_pInt:xioffsetI+5_pInt)
           enddo
         enddo
         param(instance)%newIter  (:,offset) = &
           param(instance)%oldIter  (:,offset) - &
           param(instance)%stepLength(offset)* &
           param(instance)%searchDir(:,offset)
       else
        call IO_error(400_pInt,el=el,ip=ip,ext_msg='homogenization multiphase')
       endif        
     else convergenceFullRankOne                                                                    ! not converged and residuum not improved...
       homogenization_multiphase_updateState = [.false., .true.]
       param(instance)%stepLength(offset) = &
         param(instance)%stepLength(offset)/2.0_pReal                                               ! ...try with smaller step length in same direction
       param(instance)%newIter  (:,offset) = &
         param(instance)%oldIter  (:,offset) - &
         param(instance)%stepLength(offset)* &
         param(instance)%searchDir(:,offset)
     endif convergenceFullRankOne
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
   SOURCE_chemical_energy_ID, &
   SOURCE_stochastic_phase_nucleation_ID
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
 use source_stochastic_phase_nucleation, only: &
   source_stochastic_phase_nucleation_getRateAndItsTangent   
 
 implicit none
 real(pReal),   dimension(homogenization_maxNgrains) :: &
   homogenization_multiphase_getPhaseSource, &
   selfSource, &
   drivingSource, &
   totDrivingSource, &
   nucleationSource
 real(pReal),   dimension(homogenization_maxNgrains), intent(in) :: phi
 integer(pInt),                                       intent(in) :: ip, el                          !< element number
 integer(pInt) :: &
   grI, grJ, grK, &
   homog, & 
   instance, &
   source, &
   phase
 real(pReal),   dimension(homogenization_maxNgrains,homogenization_maxNgrains) :: &
   interpTangent
 real(pReal) :: &
   tripleJunctionEnergy, &
   localSource, localSourceTangent

 homog = material_homog(ip,el)
 instance = homogenization_typeInstance(homog)

 nucleationSource = 0.0_pReal
 do grI = 1_pInt, homogenization_Ngrains(homog)
   phase = phaseAt(grI,ip,el)
   do source = 1_pInt, phase_Nsources(phase)
     select case(phase_source(source,phase)) 
       case (SOURCE_stochastic_phase_nucleation_ID)                                                  
         call source_stochastic_phase_nucleation_getRateAndItsTangent(localSource, localSourceTangent, &
                                                                      phi(grI),grI,ip,el) 
       case default
         localSource = 0.0_pReal
       
     end select
     nucleationSource(grI) = nucleationSource(grI) - localSource
   enddo
 enddo
 
 drivingSource = 0.0_pReal
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
     drivingSource(grI) = drivingSource(grI) - localSource
   enddo
 enddo
 
 totDrivingSource = 0.0_pReal
 interpTangent = homogenization_multiphase_calcPhaseFracTangent(phi,ip,el)
 do grI = 1_pInt, homogenization_Ngrains(homog)
   do grJ = 1_pInt, homogenization_Ngrains(homog)
     totDrivingSource(grI) = &
       totDrivingSource(grI) + &
       interpTangent(grJ,grI)*drivingSource(grJ)
   enddo
 enddo    

 selfSource = 0.0_pReal
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

 selfSource = 8.0_pReal*selfSource/param(instance)%InterfaceWidth
 
 homogenization_multiphase_getPhaseSource = 0.0_pReal
 do grI = 1_pInt, homogenization_Ngrains(homog)
   do grJ = 1_pInt, homogenization_Ngrains(homog)
     homogenization_multiphase_getPhaseSource(grI) = &
       homogenization_multiphase_getPhaseSource(grI) + &
       param(instance)%InterfaceMobility(grI,grJ)* &
       ((selfSource      (grI) - selfSource      (grJ))  + &
        (totDrivingSource(grI) - totDrivingSource(grJ))) + &
        (nucleationSource(grI) - nucleationSource(grJ)) 
   enddo
 enddo    
 
end function homogenization_multiphase_getPhaseSource


!--------------------------------------------------------------------------------------------------
!> @brief source tangent function for each phase 
!--------------------------------------------------------------------------------------------------
function homogenization_multiphase_getPhaseSourceTangent(phi,ip,el)
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
   SOURCE_chemical_energy_ID, &
   SOURCE_stochastic_phase_nucleation_ID
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
 use source_stochastic_phase_nucleation, only: &
  source_stochastic_phase_nucleation_getRateAndItsTangent   
 
 implicit none
 real(pReal),   dimension(homogenization_maxNgrains) :: &
   drivingSource
 real(pReal),   dimension(homogenization_maxNgrains,homogenization_maxNgrains) :: &
   homogenization_multiphase_getPhaseSourceTangent, &
   selfSourceTangent, &
   drivingSourceTangent, &
   nucleationSourceTangent
 real(pReal),   dimension(homogenization_maxNgrains), intent(in) :: phi
 integer(pInt), intent(in) :: ip, el                                                                !< element number
 integer(pInt) :: &
   grI, grJ, grK, &
   homog, & 
   source, &
   phase, &
   instance
 real(pReal) :: &
   tripleJunctionEnergy, &
   localSource, localSourceTangent  
 real(pReal),   dimension(homogenization_maxNgrains, &
                          homogenization_maxNgrains, &
                          homogenization_maxNgrains) :: &
   interp2ndTangent

 homog = material_homog(ip,el)
 instance = homogenization_typeInstance(homog)

 nucleationSourceTangent = 0.0_pReal
 do grI = 1_pInt, homogenization_Ngrains(homog)
   phase = phaseAt(grI,ip,el)
   do source = 1_pInt, phase_Nsources(phase)
     select case(phase_source(source,phase)) 
        case (SOURCE_stochastic_phase_nucleation_ID)                                                  
         call source_stochastic_phase_nucleation_getRateAndItsTangent(localSource, localSourceTangent, &
                                                                      phi(grI),grI,ip,el) 
                                                          
       case default
         localSourceTangent = 0.0_pReal
       
     end select
     nucleationSourceTangent(grI,grI) = &
       nucleationSourceTangent(grI,grI) - localSourceTangent
   enddo  
 enddo
 
 drivingSource = 0.0_pReal
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
     drivingSource(grI) = drivingSource(grI) - localSource
   enddo 
 enddo

 interp2ndTangent = homogenization_multiphase_calcPhaseFrac2ndTangent(phi,ip,el)
 drivingSourceTangent = 0.0_pReal
 do grK = 1_pInt, homogenization_Ngrains(homog)
   do grJ = 1_pInt, homogenization_Ngrains(homog)
     do grI = 1_pInt, homogenization_Ngrains(homog)
       drivingSourceTangent(grK,grJ) =  drivingSourceTangent(grK,grJ) + &
                                        interp2ndTangent(grI,grK,grJ)* &
                                        drivingSource(grI) 
     enddo
   enddo
 enddo  
     
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
 
 selfSourceTangent = 8.0_pReal*selfSourceTangent/param(instance)%InterfaceWidth

 homogenization_multiphase_getPhaseSourceTangent = 0.0_pReal
 do grI = 1_pInt, homogenization_Ngrains(homog)
   do grJ = 1_pInt, homogenization_Ngrains(homog)
     do grK = 1_pInt, homogenization_Ngrains(homog)
       homogenization_multiphase_getPhaseSourceTangent(grI,grK) = &
         homogenization_multiphase_getPhaseSourceTangent(grI,grK) + &
         param(instance)%InterfaceMobility(grI,grJ)* &
         ((selfSourceTangent      (grI,grK) - selfSourceTangent      (grJ,grK))  + &
          (drivingSourceTangent   (grI,grK) - drivingSourceTangent   (grJ,grK))) + &
          (nucleationSourceTangent(grI,grK) - nucleationSourceTangent(grJ,grK)) 
     enddo
   enddo
 enddo

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
   phasefrac(homog)%p(gr,offset) =  frac(gr) 
 enddo

end subroutine homogenization_multiphase_putPhaseFrac

!--------------------------------------------------------------------------------------------------
!> @brief calculates phase fractions (i.e. interpolation functions) for given phase field values
!--------------------------------------------------------------------------------------------------
function homogenization_multiphase_calcPhaseFrac(phi,ip,el)
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   homogenization_typeInstance
 
 implicit none
 integer(pInt),                                                intent(in)  :: ip, el                !< element number
 real(pReal),   dimension (homogenization_Ngrains(material_homog(ip,el))),     &
                                                               intent(in)  :: phi             !< array of current phase fractions
 real(pReal),   dimension (homogenization_Ngrains(material_homog(ip,el)))  :: &
   homogenization_multiphase_calcPhaseFrac
 integer(pInt) :: &
   homog, & 
   instance

 homog = material_homog(ip,el)
 instance = homogenization_typeInstance(homog)
 
 select case(param(instance)%interpolationID)
   case(linear_ID)
     homogenization_multiphase_calcPhaseFrac = phi
   
   case(moelans_ID)
     homogenization_multiphase_calcPhaseFrac = phi*phi/sum(phi*phi)
     
 end select

end function homogenization_multiphase_calcPhaseFrac


!--------------------------------------------------------------------------------------------------
!> @brief calculates tangent of the phase fractions (i.e. interpolation functions) 
!--------------------------------------------------------------------------------------------------
function homogenization_multiphase_calcPhaseFracTangent(phi,ip,el)
 use math, only: &
   math_identity2nd
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   homogenization_typeInstance
 
 implicit none
 integer(pInt),                                                intent(in)  :: ip, el                !< element number
 real(pReal),   dimension (homogenization_Ngrains(material_homog(ip,el))),     &
                                                               intent(in)  :: phi             !< array of current phase fractions
 real(pReal),   dimension (homogenization_Ngrains(material_homog(ip,el)), &
                           homogenization_Ngrains(material_homog(ip,el)))  :: &
   homogenization_multiphase_calcPhaseFracTangent
 integer(pInt) :: &
   grI, grJ, &
   homog, & 
   instance

 homog = material_homog(ip,el)
 instance = homogenization_typeInstance(homog)
 
 select case(param(instance)%interpolationID)
   case(linear_ID)
     homogenization_multiphase_calcPhaseFracTangent = &
       math_identity2nd(homogenization_Ngrains(homog))
   
   case(moelans_ID)
     homogenization_multiphase_calcPhaseFracTangent = 0.0_pReal
     do grI = 1_pInt, homogenization_Ngrains(homog)
       homogenization_multiphase_calcPhaseFracTangent(grI,grI) = &
         2.0_pReal*phi(grI)/sum(phi*phi)
       do grJ = 1_pInt, homogenization_Ngrains(homog)
         homogenization_multiphase_calcPhaseFracTangent(grI,grJ) = &
           homogenization_multiphase_calcPhaseFracTangent(grI,grJ) - &
           2.0_pReal*phi(grI)*phi(grI)*phi(grJ)/sum(phi*phi)/sum(phi*phi)
       enddo
     enddo      
     
 end select

end function homogenization_multiphase_calcPhaseFracTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates tangent of the phase fractions (i.e. interpolation functions) 
!--------------------------------------------------------------------------------------------------
function homogenization_multiphase_calcPhaseFrac2ndTangent(phi,ip,el)
 use math, only: &
   math_identity2nd
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   homogenization_typeInstance
 
 implicit none
 integer(pInt),                                                intent(in)  :: ip, el                !< element number
 real(pReal),   dimension (homogenization_Ngrains(material_homog(ip,el))),     &
                                                               intent(in)  :: phi             !< array of current phase fractions
 real(pReal),   dimension (homogenization_Ngrains(material_homog(ip,el)), &
                           homogenization_Ngrains(material_homog(ip,el)), &
                           homogenization_Ngrains(material_homog(ip,el)))  :: &
   homogenization_multiphase_calcPhaseFrac2ndTangent
 real(pReal),   dimension (homogenization_Ngrains(material_homog(ip,el)), &
                           homogenization_Ngrains(material_homog(ip,el)))  :: &
   kronDelta
 integer(pInt) :: &
   grI, grJ, grK, &
   homog, & 
   instance

 homog = material_homog(ip,el)
 instance = homogenization_typeInstance(homog)
 
 select case(param(instance)%interpolationID)
   case(linear_ID)
     homogenization_multiphase_calcPhaseFrac2ndTangent = 0.0_pReal
   
   case(moelans_ID)
     kronDelta = math_identity2nd(homogenization_Ngrains(homog))
     homogenization_multiphase_calcPhaseFrac2ndTangent = 0.0_pReal
     do grI = 1_pInt, homogenization_Ngrains(homog)
       do grJ = 1_pInt, homogenization_Ngrains(homog)
         do grK = 1_pInt, homogenization_Ngrains(homog)
           homogenization_multiphase_calcPhaseFrac2ndTangent(grI,grJ,grK) = &
             2.0_pReal*KronDelta(grI,grJ)*kronDelta(grI,grK)/sum(phi*phi) - &
             4.0_pReal*KronDelta(grI,grJ)*phi(grI)*phi(grK)/sum(phi*phi)/sum(phi*phi) - &
             4.0_pReal*kronDelta(grI,grK)*phi(grI)*phi(grJ)/sum(phi*phi)/sum(phi*phi) - &
             2.0_pReal*kronDelta(grJ,grK)*phi(grI)*phi(grI)/sum(phi*phi)/sum(phi*phi) + &
             8.0_pReal*phi(grI)*phi(grI)*phi(grJ)*phi(grK)/sum(phi*phi)/sum(phi*phi)/sum(phi*phi)
         enddo                                           
       enddo                            
     enddo 
     
 end select

end function homogenization_multiphase_calcPhaseFrac2ndTangent


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
 real(pReal),  dimension(sum(homogenization_multiphase_sizePostResult(:, &
                         homogenization_typeInstance(material_homog(ip,el))))) :: &
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
         phasefrac(material_homog(ip,el))%p(gr,phasefracMapping(material_homog(ip,el))%p(ip,el))
         c = c + 1_pInt
       enddo  
     case (ipcoords_ID)
       homogenization_multiphase_postResults(c+1_pInt:c+3_pInt) = mesh_ipCoordinates(1:3,ip,el)                       ! current ip coordinates
       c = c + 3_pInt
    end select
 enddo

end function homogenization_multiphase_postResults

end module homogenization_multiphase
