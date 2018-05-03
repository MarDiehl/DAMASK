!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating kinematics resulting from solute strain
!> @details to be done
!--------------------------------------------------------------------------------------------------
module kinematics_solute_strain
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   kinematics_solute_strain_offset, &                                                                         !< which kinematics is my current mechanism?
   kinematics_solute_strain_instance                                                                          !< instance of kinematics mechanism

 type, private :: tParameters                                                                                 !< container type for internal constitutive parameters
   real(pReal), pointer, dimension(:)   :: &                                                        
     StrainCoeff, &
     EqConc
 end type

 type(tParameters), dimension(:), allocatable, private :: param                                       !< containers of constitutive parameters (len Ninstance)

 public :: &
   kinematics_solute_strain_init, &
   kinematics_solute_strain_initialStrain, &
   kinematics_solute_strain_LiAndItsTangent, &
   kinematics_solute_strain_getMechanicalComponentPotential

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine kinematics_solute_strain_init(fileUnit)
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use debug, only: &
   debug_level,&
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
 use material, only: &
   phase_kinematics, &
   phase_Nkinematics, &
   phase_Ncomponents, &
   KINEMATICS_solute_strain_label, &
   KINEMATICS_solute_strain_ID, &
   material_Nphase, &
   MATERIAL_partPhase

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: maxNinstance,phase,instance,kinematics,j
 character(len=65536) :: &
   tag     = '', &
   line    = ''

 write(6,'(/,a)')   ' <<<+-  kinematics_'//KINEMATICS_solute_strain_LABEL//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = int(count(phase_kinematics == KINEMATICS_solute_strain_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(kinematics_solute_strain_offset(material_Nphase), source=0_pInt)
 allocate(kinematics_solute_strain_instance(material_Nphase), source=0_pInt)
 do phase = 1, material_Nphase
   kinematics_solute_strain_instance(phase) = count(phase_kinematics(:,1:phase) == kinematics_solute_strain_ID)
   do kinematics = 1, phase_Nkinematics(phase)
     if (phase_kinematics(kinematics,phase) == kinematics_solute_strain_ID) &
       kinematics_solute_strain_offset(phase) = kinematics
   enddo    
 enddo
   
 allocate(param(maxNinstance))

 rewind(fileUnit)
 phase = 0_pInt
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= MATERIAL_partPhase)         ! wind forward to <phase>
   line = IO_read(fileUnit)
 enddo
 
 parsingFile: do while (trim(line) /= IO_EOF)                                                       ! read through sections of phase part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit                                                                                           
   endif   
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next phase section
     phase = phase + 1_pInt                                                                         ! advance phase section counter
     cycle                                                                                          ! skip to next line
   endif
   if (phase > 0_pInt ) then; if (any(phase_kinematics(:,phase) == KINEMATICS_solute_strain_ID)) then         ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran
     instance = kinematics_solute_strain_instance(phase)                                                         ! which instance of my damage is present phase
     allocate(param(instance)%StrainCoeff(phase_Ncomponents(phase)), source=0.0_pReal)
     allocate(param(instance)%EqConc     (phase_Ncomponents(phase)), source=0.0_pReal)
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key...
     select case(tag)
       case ('component_straincoeff','component_straincoefficient')
         if (chunkPos(1) /= phase_Ncomponents(phase) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//KINEMATICS_solute_strain_LABEL//')')
         do j = 1_pInt, phase_Ncomponents(phase)
             param(instance)%StrainCoeff(j) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo

       case ('component_eqconc','component_equillibriumconc')
         if (chunkPos(1) /= phase_Ncomponents(phase) + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//KINEMATICS_solute_strain_LABEL//')')
         do j = 1_pInt, phase_Ncomponents(phase)
             param(instance)%EqConc(j) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo

       case default

     end select
   endif; endif
 enddo parsingFile

end subroutine kinematics_solute_strain_init

!--------------------------------------------------------------------------------------------------
!> @brief  report initial thermal strain based on current temperature deviation from reference
!--------------------------------------------------------------------------------------------------
pure function kinematics_solute_strain_initialStrain(ipc, ip, el)
 use material, only: &
   material_phase, &
   phase_Ncomponents, &
   chemicalConc, &
   chemConcMapping
 use math, only: &
   math_I3
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(3,3) :: &
   kinematics_solute_strain_initialStrain                                                           !< initial solute strain (should be small strain, though)
 integer(pInt) :: &
   phase, &
   instance, &
   cp
   
 phase = material_phase(ipc,ip,el)
 instance = kinematics_solute_strain_instance(phase)
 kinematics_solute_strain_initialStrain = 0.0_pReal
 do cp = 1_pInt, phase_Ncomponents(phase)
   kinematics_solute_strain_initialStrain = &
     kinematics_solute_strain_initialStrain + &
     param(instance)%StrainCoeff(cp)* &
     ( &
      chemicalConc(phase)%p(cp,chemConcMapping(phase)%p(ipc,ip,el)) - &
      param(instance)%EqConc(cp) &
     )*math_I3
 enddo
  
end function kinematics_solute_strain_initialStrain

!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient  
!--------------------------------------------------------------------------------------------------
subroutine kinematics_solute_strain_LiAndItsTangent(Li, dLi_dTstar3333, ipc, ip, el)
 use material, only: &
   material_phase, &
   phase_Ncomponents, &
   chemicalConcRate, &
   chemConcMapping
 use math, only: &
   math_I3
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(out), dimension(3,3) :: &
   Li                                                                                               !< thermal velocity gradient
 real(pReal),   intent(out), dimension(3,3,3,3) :: &
   dLi_dTstar3333                                                                                   !< derivative of Li with respect to Tstar (4th-order tensor defined to be zero)
 integer(pInt) :: &
   phase, &
   instance, &
   cp
   
 phase = material_phase(ipc,ip,el)
 instance = kinematics_solute_strain_instance(phase)
 Li = 0.0_pReal
 do cp = 1_pInt, phase_Ncomponents(phase)
   Li = Li + &
        param(instance)%StrainCoeff(cp)* &
        chemicalConcRate(phase)%p(cp,chemConcMapping(phase)%p(ipc,ip,el))
 enddo
 dLi_dTstar3333 = 0.0_pReal 
  
end subroutine kinematics_solute_strain_LiAndItsTangent

!--------------------------------------------------------------------------------------------------
!> @brief  return mechanical contribution to component chemical potential 
!--------------------------------------------------------------------------------------------------
function kinematics_solute_strain_getMechanicalComponentPotential(Tstar_v, ipc, ip, el)
 use material, only: &
   material_phase, &
   phase_Ncomponents
 use math, only: &
   math_Mandel6to33, &
   math_trace33
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in), dimension(6) :: &
   Tstar_v                                                                                               !< thermal velocity gradient
 real(pReal),   dimension(phase_Ncomponents(material_phase(ipc,ip,el))) :: &
   kinematics_solute_strain_getMechanicalComponentPotential                                                                                   !< derivative of Li with respect to Tstar (4th-order tensor defined to be zero)
 integer(pInt) :: &
   phase, &
   instance, &
   cp
   
 phase = material_phase(ipc,ip,el)
 instance = kinematics_solute_strain_instance(phase)
 kinematics_solute_strain_getMechanicalComponentPotential = 0.0_pReal
 do cp = 1_pInt, phase_Ncomponents(phase)
   kinematics_solute_strain_getMechanicalComponentPotential(cp) = &
        param(instance)%StrainCoeff(cp)* &
        math_trace33(math_Mandel6to33(Tstar_v))
 enddo
  
end function kinematics_solute_strain_getMechanicalComponentPotential

end module kinematics_solute_strain
