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
 real(pReal),                 dimension(:),   allocatable,         private :: &
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
                 isostress_ID, &
                 rankone_ID, &
                 laminate_ID
 end enum

 integer(kind(isostrain_ID)), dimension(:),   allocatable,         private :: &
  homogenization_multiphase_mixtureID                                                                !< ID of mixture rule


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
   math_EulerToR
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
   instance
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
             case('isostress')
               homogenization_multiphase_mixtureID(i) = isostress_ID
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
         homogState(homog)%sizeState =                                          0_pInt
         
       case(isostress_ID)
         homogState(homog)%sizeState = (homogenization_Ngrains(homog) - 1_pInt)*9_pInt
         
       case(rankone_ID)
         homogState(homog)%sizeState = (homogenization_Ngrains(homog) - 1_pInt)*3_pInt

       case(laminate_ID)
         homogState(homog)%sizeState = (homogenization_Ngrains(homog) - 1_pInt)*6_pInt

     end select
     homogState(homog)%sizePostResults = homogenization_multiphase_sizePostResults(instance)
     allocate(homogState(homog)%state0   (homogState(homog)%sizeState,NofMyHomog), source=0.0_pReal)
     allocate(homogState(homog)%subState0(homogState(homog)%sizeState,NofMyHomog), source=0.0_pReal)
     allocate(homogState(homog)%state    (homogState(homog)%sizeState,NofMyHomog), source=0.0_pReal)

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
     myHomog2: if (homogenization_type(homog) == HOMOGENIZATION_MULTIPHASE_ID) then
       do gr = 1_pInt, homogenization_Ngrains(homog)
         phasefrac(gr,homog)%p(phasefracMapping(homog)%p(ip,el)) = microstructure_fraction(gr,micro)
       enddo    
       instance = homogenization_typeInstance(homog)
       select case(homogenization_multiphase_mixtureID(instance))
         case(isostrain_ID)
           
         case(isostress_ID)
           do gr = 1_pInt, homogenization_Ngrains(homog)-1_pInt
             homogState(homog)%state0   (9_pInt*gr-8_pInt:9_pInt*gr,mappingHomogenization(1,ip,el)) = &
             reshape(math_mul33x33(math_transpose33(math_EulerToR(material_EulerAngles(1:3,gr,ip,el))), &
                                   crystallite_Fp0(1:3,1:3,gr,ip,el)) - &
                     math_I3,[9])
             homogState(homog)%subState0(9_pInt*gr-8_pInt:9_pInt*gr,mappingHomogenization(1,ip,el)) = &
             reshape(math_mul33x33(math_transpose33(math_EulerToR(material_EulerAngles(1:3,gr,ip,el))), &
                                   crystallite_Fp0(1:3,1:3,gr,ip,el)) - &
                     math_I3,[9])
             homogState(homog)%state    (9_pInt*gr-8_pInt:9_pInt*gr,mappingHomogenization(1,ip,el)) = &
             reshape(math_mul33x33(math_transpose33(math_EulerToR(material_EulerAngles(1:3,gr,ip,el))), &
                                   crystallite_Fp0(1:3,1:3,gr,ip,el)) - &
                     math_I3,[9])
           enddo    
           do gr = 1_pInt, homogenization_Ngrains(homog)-1
             crystallite_F0(1:3,1:3,gr,ip,el) = &
               math_I3 + &
               reshape(homogState(homog)%state(9_pInt*gr-8_pInt:9_pInt*gr, &
                                               mappingHomogenization(1,ip,el)),[3,3])
             do o = 1_pInt,  homogenization_Ngrains(homog)-1
               crystallite_F0(1:3,1:3,gr,ip,el) = &
                 crystallite_F0(1:3,1:3,gr,ip,el) - &
                 phasefrac(o,homog)%p(phasefracMapping(homog)%p(ip,el))* &
                 reshape(homogState(homog)%state(9_pInt*o-8_pInt:9_pInt*o, &
                                                 mappingHomogenization(1,ip,el)),[3,3])
             enddo
             crystallite_Fe(1:3,1:3,gr,ip,el) = math_mul33x33(crystallite_F0(1:3,1:3,gr,ip,el), &
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
subroutine homogenization_multiphase_partitionDeformation(F,avgF,ip,el)
 use mesh, only: &
   mesh_element
 use material, only: &
   phasefrac, &
   phasefracMapping, &
   homogState, &
   mappingHomogenization, &
   homogenization_maxNgrains, &
   homogenization_Ngrains, &
   homogenization_typeInstance
 
 implicit none
 real(pReal),   dimension (3,3,homogenization_maxNgrains), intent(out) :: F                         !< partioned def grad per grain
 real(pReal),   dimension (3,3),                           intent(in)  :: avgF                      !< my average def grad
 integer(pInt),                                            intent(in)  :: &
   el, ip                                                                                           !< element number
 integer(pInt) :: homog, instance, grI, grJ  
 
 homog = mesh_element(3,el)
 instance = homogenization_typeInstance(homog)
 select case(homogenization_multiphase_mixtureID(instance))
   case(isostrain_ID)
     F(1:3,1:3,1:homogenization_Ngrains(mesh_element(3,el)))= &
       spread(avgF,3,homogenization_Ngrains(mesh_element(3,el)))
   
   case(isostress_ID)
     do grI = 1_pInt, homogenization_Ngrains(homog)-1
       F(1:3,1:3,grI) = &
         avgF + &
         reshape(homogState(homog)%state(9_pInt*grI-8_pInt:9_pInt*grI, &
                                         mappingHomogenization(1,ip,el)),[3,3])
       do grJ = 1_pInt,  homogenization_Ngrains(homog)-1
         F(1:3,1:3,grI) = &
           F(1:3,1:3,grI) - &
           phasefrac(grJ,homog)%p(phasefracMapping(homog)%p(ip,el))* &
           reshape(homogState(homog)%state(9_pInt*grJ-8_pInt:9_pInt*grJ, &
                                           mappingHomogenization(1,ip,el)),[3,3])
       enddo
     enddo    
     F(1:3,1:3,homogenization_Ngrains(homog)) = avgF
     do grJ = 1_pInt,  homogenization_Ngrains(homog)-1
       F(1:3,1:3,homogenization_Ngrains(homog)) = &
         F(1:3,1:3,homogenization_Ngrains(homog)) - &
         phasefrac(grJ,homog)%p(phasefracMapping(homog)%p(ip,el))* &
         reshape(homogState(homog)%state(9_pInt*grJ-8_pInt:9_pInt*grJ, &
                                         mappingHomogenization(1,ip,el)),[3,3])
     enddo

   case(rankone_ID)

   case(laminate_ID)

 end select

end subroutine homogenization_multiphase_partitionDeformation


!--------------------------------------------------------------------------------------------------
!> @brief derive average stress and stiffness from constituent quantities 
!--------------------------------------------------------------------------------------------------
subroutine homogenization_multiphase_averageStressAndItsTangent(avgP,dAvgPdAvgF,P,dPdF,ip,el)
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_maxNgrains, &
   homogenization_Ngrains, &
   phasefracMapping, &
   phasefrac
 
 implicit none
 real(pReal),   dimension (3,3),                               intent(out) :: avgP                  !< average stress at material point
 real(pReal),   dimension (3,3,3,3),                           intent(out) :: dAvgPdAvgF            !< average stiffness at material point
 real(pReal),   dimension (3,3,homogenization_maxNgrains),     intent(in)  :: P                     !< array of current grain stresses
 real(pReal),   dimension (3,3,3,3,homogenization_maxNgrains), intent(in)  :: dPdF                  !< array of current grain stiffnesses
 integer(pInt),                                                intent(in)  :: ip, el                !< element number
 integer(pInt) :: &
   gr, &
   homog, & 
   offset

 homog = mesh_element(3,el)
 offset = phasefracMapping(homog)%p(ip,el)
 avgP = 0.0_pReal
 dAvgPdAvgF = 0.0_pReal
 do gr = 1_pInt, homogenization_Ngrains(homog)
   avgP = avgP + &
          phasefrac(gr,homog)%p(offset)*P(1:3,1:3,gr) 
   dAvgPdAvgF = dAvgPdAvgF + &
                phasefrac(gr,homog)%p(offset)*dPdF(1:3,1:3,1:3,1:3,gr) 
 enddo

end subroutine homogenization_multiphase_averageStressAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief state update for different mixture rules 
!--------------------------------------------------------------------------------------------------
function homogenization_multiphase_updateState(ip,el)
 use IO, only: &
   IO_error
 use mesh, only: &
   mesh_element
 use material, only: &
   phasefrac, &
   phasefracMapping, &
   homogState, &
   homogenization_Ngrains, &
   homogenization_typeInstance
 use crystallite, only: &
   crystallite_P, &
   crystallite_dPdF
 
 implicit none
 logical,       dimension (2)               :: homogenization_multiphase_updateState                !< average stress at material point
 integer(pInt),                  intent(in) :: ip, el                                               !< element number
 real(pReal),   dimension(:),   allocatable :: residual  
 real(pReal),   dimension(:,:), allocatable :: jacobian
 integer(pInt), dimension(9)                :: ipiv                                                 !< needed for matrix inversion by LAPACK
 real(pReal),   dimension(3,3)              :: avgP                                                 !< average stress at material point
 integer(pInt) :: &
   homog, & 
   instance, &
   offset, grI, grJ, &
   ierr
 real(pReal) :: &
   stressTol

 external :: &
   dgesv

 homog = mesh_element(3,el)
 instance = homogenization_typeInstance(homog)
 offset = phasefracMapping(homog)%p(ip,el)
 myMixRule: select case(homogenization_multiphase_mixtureID(instance))
   case(isostrain_ID) myMixRule
     homogenization_multiphase_updateState = [.true., .true.]

   case(isostress_ID) myMixRule
     allocate(residual(homogState(homog)%sizeState), source=0.0_pReal)
     allocate(jacobian(homogState(homog)%sizeState, &
                       homogState(homog)%sizeState), source=0.0_pReal)
     avgP = crystallite_P(1:3,1:3,homogenization_Ngrains(homog),ip,el)
     do grI = 1_pInt, homogenization_Ngrains(homog)-1
       residual(9_pInt*grI-8_pInt:9_pInt*grI) = &
        reshape(crystallite_P(1:3,1:3,grI                          ,ip,el) - &
                crystallite_P(1:3,1:3,homogenization_Ngrains(homog),ip,el), [9_pInt])
       avgP = avgP + crystallite_P(1:3,1:3,grI,ip,el)
     enddo
     avgP = avgP/real(homogenization_Ngrains(homog),pReal)
     stressTol = max(norm2(avgP)*homogenization_multiphase_relTol(instance), &
                                 homogenization_multiphase_absTol(instance))
     convergence: if (norm2(residual) < stressTol) then
       homogenization_multiphase_updateState = [.true., .true.]
       exit convergence
     else convergence
       homogenization_multiphase_updateState = [.false., .true.]
       do grI = 1_pInt,  homogenization_Ngrains(homog)-1
         jacobian(9_pInt*grI-8_pInt:9_pInt*grI,9_pInt*grI-8_pInt:9_pInt*grI) = &
           reshape(crystallite_dPdF(1:3,1:3,1:3,1:3,grI,ip,el),[9,9])
         do grJ = 1_pInt,  homogenization_Ngrains(homog)-1
           jacobian(9_pInt*grI-8_pInt:9_pInt*grI,9_pInt*grJ-8_pInt:9_pInt*grJ) = &
             jacobian(9_pInt*grI-8_pInt:9_pInt*grI,9_pInt*grJ-8_pInt:9_pInt*grJ) - &
             phasefrac(grJ,homog)%p(offset)* &
             reshape(crystallite_dPdF(1:3,1:3,1:3,1:3,grI,ip,el),[9,9]) + &
             phasefrac(grJ,homog)%p(offset)* &
             reshape(crystallite_dPdF(1:3,1:3,1:3,1:3,homogenization_Ngrains(homog),ip,el),[9,9])
         enddo
       enddo    
       call dgesv(homogState(homog)%sizeState,1, &
                  jacobian, &
                  homogState(homog)%sizeState,ipiv, &
                  residual, &
                  homogState(homog)%sizeState,ierr)                                                   !< solve Jacobian * delta state = -residual for delta state
       if (ierr == 0_pInt) then
         homogState(homog)%state(1:homogState(homog)%sizeState,offset) = &
           homogState(homog)%state(1:homogState(homog)%sizeState,offset) - &
           residual
       else
        call IO_error(400_pInt,el=el,ip=ip,ext_msg='homogenization multiphase')
       endif        
     endif convergence 
     
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
