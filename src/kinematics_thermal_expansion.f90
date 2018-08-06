!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating kinematics resulting from thermal expansion
!> @details to be done
!--------------------------------------------------------------------------------------------------
module kinematics_thermal_expansion
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   kinematics_thermal_expansion_offset, &                                                                         !< which kinematics is my current mechanism?
   kinematics_thermal_expansion_instance                                                                          !< instance of kinematics mechanism

 type, private :: tParameters                                                                                 !< container type for internal constitutive parameters
   real(pReal) :: &                                                        
     TRef = 300.0_pReal
   real(pReal), dimension(3,3,3) :: &
     alpha = 0.0_pReal
                                                            
 end type

 type(tParameters), dimension(:), allocatable, private :: param                                       !< containers of constitutive parameters (len Ninstance)

 public :: &
   kinematics_thermal_expansion_init, &
   kinematics_thermal_expansion_initialStrain, &
   kinematics_thermal_expansion_LiAndItsTangent

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine kinematics_thermal_expansion_init(fileUnit)
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
   kinematics_thermal_expansion_label, &
   kinematics_thermal_expansion_ID
 use config, only: &
   material_Nphase, &
   MATERIAL_partPhase

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: maxNinstance,phase,instance,kinematics,i
 character(len=65536) :: &
   tag     = '', &
   line    = ''

 write(6,'(/,a)')   ' <<<+-  kinematics_'//kinematics_thermal_expansion_LABEL//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = int(count(phase_kinematics == kinematics_thermal_expansion_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(param(maxNinstance))
 allocate(kinematics_thermal_expansion_offset(material_Nphase), source=0_pInt)
 allocate(kinematics_thermal_expansion_instance(material_Nphase), source=0_pInt)

 do phase = 1, material_Nphase
   kinematics_thermal_expansion_instance(phase) = count(phase_kinematics(:,1:phase) == kinematics_thermal_expansion_ID)
   do kinematics = 1, phase_Nkinematics(phase)
     if (phase_kinematics(kinematics,phase) == kinematics_thermal_expansion_ID) &
       kinematics_thermal_expansion_offset(phase) = kinematics
   enddo
 enddo
 
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
   if (phase > 0_pInt ) then; if (any(phase_kinematics(:,phase) == kinematics_thermal_expansion_ID)) then         ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran
     instance = kinematics_thermal_expansion_instance(phase)                                                         ! which instance of my damage is present phase
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key...
     select case(tag)
       case ('thermal_expansion11')
         do i = 2_pInt, min(4,chunkPos(1))                                                        ! read up to three parameters (constant, linear, quadratic with T)
           param(instance)%alpha(1,1,i-1) = IO_floatValue(line,chunkPos,i)
         enddo

       case ('thermal_expansion22')
         do i = 2_pInt, min(4,chunkPos(1))                                                        ! read up to three parameters (constant, linear, quadratic with T)
           param(instance)%alpha(2,2,i-1) = IO_floatValue(line,chunkPos,i)
         enddo

       case ('thermal_expansion33')
         do i = 2_pInt, min(4,chunkPos(1))                                                        ! read up to three parameters (constant, linear, quadratic with T)
           param(instance)%alpha(3,3,i-1) = IO_floatValue(line,chunkPos,i)
         enddo

       case ('reference_temperature')
         param(instance)%TRef = IO_floatValue(line,chunkPos,2)

       case default

     end select
   endif; endif
 enddo parsingFile

end subroutine kinematics_thermal_expansion_init

!--------------------------------------------------------------------------------------------------
!> @brief  report initial thermal strain based on current temperature deviation from reference
!--------------------------------------------------------------------------------------------------
pure function kinematics_thermal_expansion_initialStrain(ipc, ip, el)
 use material, only: &
   material_phase, &
   temperature, &
   thermalMapping
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(3,3) :: &
   kinematics_thermal_expansion_initialStrain                                                       !< initial thermal strain (should be small strain, though)
 integer(pInt) :: &
   phase, &
   instance, &
   offset
   
 phase  = material_phase(ipc,ip,el)
 instance = kinematics_thermal_expansion_instance(phase)
 offset = thermalMapping(phase)%p(ipc,ip,el)
 
 kinematics_thermal_expansion_initialStrain = &
   (temperature(phase)%p(offset) - param(instance)%TRef)**1 / 1. * &
   param(instance)%alpha(1:3,1:3,1) + &                                                  ! constant  coefficient
   (temperature(phase)%p(offset) - param(instance)%TRef)**2 / 2. * &
   param(instance)%alpha(1:3,1:3,2) + &                                                  ! linear    coefficient
   (temperature(phase)%p(offset) - param(instance)%TRef)**3 / 3. * &
   param(instance)%alpha(1:3,1:3,3)                                                      ! quadratic coefficient
  
end function kinematics_thermal_expansion_initialStrain

!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient  
!--------------------------------------------------------------------------------------------------
subroutine kinematics_thermal_expansion_LiAndItsTangent(Li, dLi_dTstar3333, ipc, ip, el)
 use material, only: &
   material_phase, &
   temperature, &
   temperatureRate, &
   thermalMapping
 
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
   offset
 real(pReal) :: &
   T, TRef, TDot  
   
 phase = material_phase(ipc,ip,el)
 instance = kinematics_thermal_expansion_instance(phase)
 offset = thermalMapping(phase)%p(ipc,ip,el)
 T = temperature(phase)%p(offset)
 TDot = temperatureRate(phase)%p(offset)
 TRef = param(instance)%TRef
 
 Li = TDot * ( &
               param(instance)%alpha(1:3,1:3,1)*(T - TRef)**0 &                           ! constant  coefficient
             + param(instance)%alpha(1:3,1:3,2)*(T - TRef)**1 &                           ! linear    coefficient
             + param(instance)%alpha(1:3,1:3,3)*(T - TRef)**2 &                           ! quadratic coefficient
             ) / &
      (1.0_pReal \
            + param(instance)%alpha(1:3,1:3,1)*(T - TRef)**1 / 1. &
            + param(instance)%alpha(1:3,1:3,2)*(T - TRef)**2 / 2. &
            + param(instance)%alpha(1:3,1:3,3)*(T - TRef)**3 / 3. &
      )
 dLi_dTstar3333 = 0.0_pReal 
  
end subroutine kinematics_thermal_expansion_LiAndItsTangent

end module kinematics_thermal_expansion
