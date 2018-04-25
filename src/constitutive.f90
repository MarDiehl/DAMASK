!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief elasticity, plasticity, internal microstructure state
!--------------------------------------------------------------------------------------------------
module constitutive
 use prec, only: &
   pInt

 implicit none
 private
 integer(pInt), public, protected :: &
   constitutive_plasticity_maxSizePostResults, &
   constitutive_plasticity_maxSizeDotState, &
   constitutive_chemicalFE_maxSizePostResults, &
   constitutive_chemicalFE_maxSizeDotState, &
   constitutive_source_maxSizePostResults, &
   constitutive_source_maxSizeDotState

 public :: &
   constitutive_init, &
   constitutive_homogenizedC, &
   constitutive_microstructure, &
   constitutive_LpAndItsTangent, &
   constitutive_LiAndItsTangent, &
   constitutive_initialFi, &
   constitutive_TandItsTangent, &
   constitutive_collectDotState, &
   constitutive_collectDeltaState, &
   constitutive_postResults

 private :: &
   constitutive_hooke_TandItsTangent

contains


!--------------------------------------------------------------------------------------------------
!> @brief allocates arrays pointing to array of the various constitutive modules
!--------------------------------------------------------------------------------------------------
subroutine constitutive_init()
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use prec, only: &
   pReal
 use debug, only: &
   debug_constitutive, &
   debug_levelBasic
 use numerics, only: &
   worldrank
 use IO, only: &
   IO_error, &
   IO_open_file, &
   IO_checkAndRewind, &
   IO_open_jobFile_stat, &
   IO_write_jobFile, &
   IO_write_jobIntFile, &
   IO_timeStamp
 use mesh, only: &
   FE_geomtype
 use material, only: &
   material_phase, &
   material_Nphase, &
   material_localFileExt, &
   material_configFile, &
   phase_name, &
   phase_plasticity, &
   phase_plasticityInstance, &
   phase_chemicalFE, &
   phase_chemicalFEInstance, &
   phase_Nsources, &
   phase_source, &
   phase_kinematics, &
   ELASTICITY_hooke_ID, &
   PLASTICITY_none_ID, &
   PLASTICITY_isotropic_ID, &
   PLASTICITY_phenopowerlaw_ID, &
   PLASTICITY_dislotwin_ID, &
   PLASTICITY_disloucla_ID, &
   PLASTICITY_nonlocal_ID ,&
   CHEMICALFE_none_ID, &
   CHEMICALFE_quadenergy_ID, &
   CHEMICALFE_thermodynamic_ID, &
   SOURCE_thermal_dissipation_ID, &
   SOURCE_thermal_externalheat_ID, &
   SOURCE_elastic_energy_ID, &
   SOURCE_plastic_energy_ID, &
   SOURCE_chemical_energy_ID, &
   KINEMATICS_cleavage_opening_ID, &
   KINEMATICS_slipplane_opening_ID, &
   KINEMATICS_thermal_expansion_ID, &
   ELASTICITY_HOOKE_label, &
   PLASTICITY_NONE_label, &
   PLASTICITY_ISOTROPIC_label, &
   PLASTICITY_PHENOPOWERLAW_label, &
   PLASTICITY_DISLOTWIN_label, &
   PLASTICITY_DISLOUCLA_label, &
   PLASTICITY_NONLOCAL_label, &
   CHEMICALFE_none_label, &
   CHEMICALFE_quadenergy_label, &
   CHEMICALFE_thermodynamic_label, &
   SOURCE_thermal_dissipation_label, &
   SOURCE_thermal_externalheat_label, &
   SOURCE_elastic_energy_label, &
   SOURCE_plastic_energy_label, &
   SOURCE_chemical_energy_label, &
   plasticState, &
   chemicalState, &
   sourceState

 use plastic_none
 use plastic_isotropic
 use plastic_phenopowerlaw
 use plastic_dislotwin
 use plastic_disloucla
 use plastic_nonlocal
 use chemicalFE_none
 use chemicalFE_quadenergy
 use chemicalFE_thermodynamic
 use source_thermal_dissipation
 use source_thermal_externalheat
 use source_elastic_energy
 use source_plastic_energy
 use source_chemical_energy
 use kinematics_cleavage_opening
 use kinematics_slipplane_opening
 use kinematics_thermal_expansion

 implicit none
 integer(pInt), parameter :: FILEUNIT = 200_pInt
 integer(pInt) :: &
   o, &                                                                                             !< counter in output loop
   p, &                                                                                             !< counter in phase loop
   s, &                                                                                             !< counter in source loop
   ins                                                                                              !< instance of plasticity/source

 integer(pInt), dimension(:,:), pointer :: thisSize
 integer(pInt), dimension(:)  , pointer :: thisNoutput
 character(len=64), dimension(:,:), pointer :: thisOutput
 character(len=32) :: outputName                                                                    !< name of output, intermediate fix until HDF5 output is ready
 logical :: knownPlasticity, knownchemicalFE, knownSource, nonlocalConstitutionPresent
 nonlocalConstitutionPresent = .false.

!--------------------------------------------------------------------------------------------------
! open material.config
 if (.not. IO_open_jobFile_stat(FILEUNIT,material_localFileExt)) &                                  ! no local material configuration present...
   call IO_open_file(FILEUNIT,material_configFile)                                                  ! ... open material.config file

!--------------------------------------------------------------------------------------------------
! parse plasticities from config file
 if (any(phase_plasticity == PLASTICITY_NONE_ID))          call plastic_none_init
 if (any(phase_plasticity == PLASTICITY_ISOTROPIC_ID))     call plastic_isotropic_init(FILEUNIT)
 if (any(phase_plasticity == PLASTICITY_PHENOPOWERLAW_ID)) call plastic_phenopowerlaw_init(FILEUNIT)
 if (any(phase_plasticity == PLASTICITY_DISLOTWIN_ID))     call plastic_dislotwin_init(FILEUNIT)
 if (any(phase_plasticity == PLASTICITY_DISLOUCLA_ID))     call plastic_disloucla_init(FILEUNIT)
 if (any(phase_plasticity == PLASTICITY_NONLOCAL_ID)) then
  call plastic_nonlocal_init(FILEUNIT)
  call plastic_nonlocal_stateInit()
 endif

!--------------------------------------------------------------------------------------------------
! parse chemical FE models from config file
 if (any(phase_chemicalFE == CHEMICALFE_none_ID))          call chemicalFE_none_init
 if (any(phase_chemicalFE == CHEMICALFE_quadenergy_ID))    call chemicalFE_quadenergy_init(FILEUNIT)
 if (any(phase_chemicalFE == CHEMICALFE_thermodynamic_ID)) call chemicalFE_thermodynamic_init(FILEUNIT)

!--------------------------------------------------------------------------------------------------
! parse source mechanisms from config file
 call IO_checkAndRewind(FILEUNIT)
 if (any(phase_source == SOURCE_thermal_dissipation_ID))     call source_thermal_dissipation_init(FILEUNIT)
 if (any(phase_source == SOURCE_thermal_externalheat_ID))    call source_thermal_externalheat_init(FILEUNIT)
 if (any(phase_source == SOURCE_elastic_energy_ID))          call source_elastic_energy_init
 if (any(phase_source == SOURCE_plastic_energy_ID))          call source_plastic_energy_init
 if (any(phase_source == SOURCE_chemical_energy_ID))         call source_chemical_energy_init

!--------------------------------------------------------------------------------------------------
! parse kinematic mechanisms from config file
 call IO_checkAndRewind(FILEUNIT)
 if (any(phase_kinematics == KINEMATICS_cleavage_opening_ID))  call kinematics_cleavage_opening_init(FILEUNIT)
 if (any(phase_kinematics == KINEMATICS_slipplane_opening_ID)) call kinematics_slipplane_opening_init(FILEUNIT)
 if (any(phase_kinematics == KINEMATICS_thermal_expansion_ID)) call kinematics_thermal_expansion_init(FILEUNIT)
 close(FILEUNIT)

 write(6,'(/,a)')   ' <<<+-  constitutive init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 mainProcess: if (worldrank == 0) then
!--------------------------------------------------------------------------------------------------
! write description file for constitutive output
   call IO_write_jobFile(FILEUNIT,'outputConstitutive')
   PhaseLoop: do p = 1_pInt,material_Nphase
     activePhase: if (any(material_phase == p)) then
       ins = phase_plasticityInstance(p)
       knownPlasticity = .true.                                                                     ! assume valid
       plasticityType: select case(phase_plasticity(p))
         case (PLASTICITY_NONE_ID) plasticityType
           outputName = PLASTICITY_NONE_label
           thisNoutput => null()
           thisOutput => null()
           thisSize   => null()
         case (PLASTICITY_ISOTROPIC_ID) plasticityType
           outputName = PLASTICITY_ISOTROPIC_label
           thisNoutput => plastic_isotropic_Noutput
           thisOutput => plastic_isotropic_output
           thisSize   => plastic_isotropic_sizePostResult
         case (PLASTICITY_PHENOPOWERLAW_ID) plasticityType
           outputName = PLASTICITY_PHENOPOWERLAW_label
           thisNoutput => plastic_phenopowerlaw_Noutput
           thisOutput => plastic_phenopowerlaw_output
           thisSize   => plastic_phenopowerlaw_sizePostResult
         case (PLASTICITY_DISLOTWIN_ID) plasticityType
           outputName = PLASTICITY_DISLOTWIN_label
           thisNoutput => plastic_dislotwin_Noutput
           thisOutput => plastic_dislotwin_output
           thisSize   => plastic_dislotwin_sizePostResult
         case (PLASTICITY_DISLOUCLA_ID) plasticityType
           outputName = PLASTICITY_DISLOUCLA_label
           thisNoutput => plastic_disloucla_Noutput
           thisOutput => plastic_disloucla_output
           thisSize   => plastic_disloucla_sizePostResult
         case (PLASTICITY_NONLOCAL_ID) plasticityType
           outputName = PLASTICITY_NONLOCAL_label
           thisNoutput => plastic_nonlocal_Noutput
           thisOutput => plastic_nonlocal_output
           thisSize   => plastic_nonlocal_sizePostResult
         case default plasticityType
           knownPlasticity = .false.
       end select plasticityType
       write(FILEUNIT,'(/,a,/)') '['//trim(phase_name(p))//']'
       if (knownPlasticity) then
         write(FILEUNIT,'(a)') '(plasticity)'//char(9)//trim(outputName)
         if (phase_plasticity(p) /= PLASTICITY_NONE_ID) then
           OutputPlasticityLoop: do o = 1_pInt,thisNoutput(ins)
             write(FILEUNIT,'(a,i4)') trim(thisOutput(o,ins))//char(9),thisSize(o,ins)
           enddo OutputPlasticityLoop
         endif
       endif
       ins = phase_chemicalFEInstance(p)
       knownchemicalFE = .true.                                                                     ! assume valid
       chemicalFEType: select case(phase_chemicalFE(p))
         case (CHEMICALFE_none_ID) chemicalFEType
           outputName = CHEMICALFE_none_label
           thisNoutput => null()
           thisOutput => null()
           thisSize   => null()
         case (CHEMICALFE_quadenergy_ID) chemicalFEType
           outputName = CHEMICALFE_quadenergy_label
           thisNoutput => chemicalFE_quadenergy_Noutput
           thisOutput => chemicalFE_quadenergy_output
           thisSize   => chemicalFE_quadenergy_sizePostResult
         case (CHEMICALFE_thermodynamic_ID) chemicalFEType
           outputName = CHEMICALFE_thermodynamic_label
           thisNoutput => chemicalFE_thermodynamic_Noutput
           thisOutput => chemicalFE_thermodynamic_output
           thisSize   => chemicalFE_thermodynamic_sizePostResult
         case default chemicalFEType
           knownchemicalFE = .false.
       end select chemicalFEType
       if (knownchemicalFE) then
         write(FILEUNIT,'(a)') '(chemicalfe)'//char(9)//trim(outputName)
         if (phase_chemicalFE(p) /= CHEMICALFE_none_ID) then
           OutputChemicalFELoop: do o = 1_pInt,thisNoutput(ins)
             write(FILEUNIT,'(a,i4)') trim(thisOutput(o,ins))//char(9),thisSize(o,ins)
           enddo OutputChemicalFELoop
         endif
       endif
       SourceLoop: do s = 1_pInt, phase_Nsources(p)
         knownSource = .true.                                                                       ! assume valid
         sourceType: select case (phase_source(s,p))
           case (SOURCE_thermal_dissipation_ID) sourceType
             ins = source_thermal_dissipation_instance(p)
             outputName = SOURCE_thermal_dissipation_label
             thisNoutput => source_thermal_dissipation_Noutput
             thisOutput => source_thermal_dissipation_output
             thisSize   => source_thermal_dissipation_sizePostResult
           case (SOURCE_thermal_externalheat_ID) sourceType
             ins = source_thermal_externalheat_instance(p)
             outputName = SOURCE_thermal_externalheat_label
             thisNoutput => source_thermal_externalheat_Noutput
             thisOutput => source_thermal_externalheat_output
             thisSize   => source_thermal_externalheat_sizePostResult
           case (SOURCE_elastic_energy_ID) sourceType
             ins = source_elastic_energy_instance(p)
             outputName = SOURCE_elastic_energy_label
             thisNoutput => source_elastic_energy_Noutput
             thisOutput => source_elastic_energy_output
             thisSize   => source_elastic_energy_sizePostResult
           case (SOURCE_plastic_energy_ID) sourceType
             ins = source_plastic_energy_instance(p)
             outputName = SOURCE_plastic_energy_label
             thisNoutput => source_plastic_energy_Noutput
             thisOutput => source_plastic_energy_output
             thisSize   => source_plastic_energy_sizePostResult
           case (SOURCE_chemical_energy_ID) sourceType
             ins = source_chemical_energy_instance(p)
             outputName = SOURCE_chemical_energy_label
             thisNoutput => source_chemical_energy_Noutput
             thisOutput => source_chemical_energy_output
             thisSize   => source_chemical_energy_sizePostResult
           case default sourceType
             knownSource = .false.
         end select sourceType
         if (knownSource) then
           write(FILEUNIT,'(a)') '(source)'//char(9)//trim(outputName)
           OutputSourceLoop: do o = 1_pInt,thisNoutput(ins)
             write(FILEUNIT,'(a,i4)') trim(thisOutput(o,ins))//char(9),thisSize(o,ins)
           enddo OutputSourceLoop
         endif
       enddo SourceLoop
     endif activePhase
   enddo PhaseLoop
   close(FILEUNIT)
 endif mainProcess

 constitutive_plasticity_maxSizeDotState = 0_pInt
 constitutive_plasticity_maxSizePostResults = 0_pInt
 constitutive_chemicalFE_maxSizeDotState = 0_pInt
 constitutive_chemicalFE_maxSizePostResults = 0_pInt
 constitutive_source_maxSizeDotState = 0_pInt
 constitutive_source_maxSizePostResults = 0_pInt

 PhaseLoop2:do p = 1_pInt,material_Nphase
!--------------------------------------------------------------------------------------------------
! partition and inititalize state
   plasticState(p)%partionedState0 = plasticState(p)%State0
   plasticState(p)%State           = plasticState(p)%State0
   chemicalState(p)%partionedState0= chemicalState(p)%State0
   chemicalState(p)%State          = chemicalState(p)%State0
   forall(s = 1_pInt:phase_Nsources(p))
     sourceState(p)%p(s)%partionedState0 = sourceState(p)%p(s)%State0
     sourceState(p)%p(s)%State           = sourceState(p)%p(s)%State0
   end forall
!--------------------------------------------------------------------------------------------------
! determine max size of state and output
   constitutive_plasticity_maxSizeDotState    = max(constitutive_plasticity_maxSizeDotState,    &
                                                    plasticState(p)%sizeDotState)
   constitutive_plasticity_maxSizePostResults = max(constitutive_plasticity_maxSizePostResults, &
                                                    plasticState(p)%sizePostResults)
   constitutive_chemicalFE_maxSizeDotState    = max(constitutive_chemicalFE_maxSizeDotState,    &
                                                    plasticState(p)%sizeDotState)
   constitutive_chemicalFE_maxSizePostResults = max(constitutive_chemicalFE_maxSizePostResults, &
                                                    plasticState(p)%sizePostResults)
   constitutive_source_maxSizeDotState        = max(constitutive_source_maxSizeDotState, &
                                                    maxval(sourceState(p)%p(:)%sizeDotState))
   constitutive_source_maxSizePostResults     = max(constitutive_source_maxSizePostResults, &
                                                    maxval(sourceState(p)%p(:)%sizePostResults))
 enddo PhaseLoop2


#ifdef TODO
!--------------------------------------------------------------------------------------------------
! report
 constitutive_maxSizeState       = maxval(constitutive_sizeState)
 constitutive_plasticity_maxSizeDotState    = maxval(constitutive_sizeDotState)

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) then
   write(6,'(a32,1x,7(i8,1x))')   'constitutive_state0:          ', shape(constitutive_state0)
   write(6,'(a32,1x,7(i8,1x))')   'constitutive_partionedState0: ', shape(constitutive_partionedState0)
   write(6,'(a32,1x,7(i8,1x))')   'constitutive_subState0:       ', shape(constitutive_subState0)
   write(6,'(a32,1x,7(i8,1x))')   'constitutive_state:           ', shape(constitutive_state)
   write(6,'(a32,1x,7(i8,1x))')   'constitutive_aTolState:       ', shape(constitutive_aTolState)
   write(6,'(a32,1x,7(i8,1x))')   'constitutive_dotState:        ', shape(constitutive_dotState)
   write(6,'(a32,1x,7(i8,1x))')   'constitutive_deltaState:      ', shape(constitutive_deltaState)
   write(6,'(a32,1x,7(i8,1x))')   'constitutive_sizeState:       ', shape(constitutive_sizeState)
   write(6,'(a32,1x,7(i8,1x))')   'constitutive_sizeDotState:    ', shape(constitutive_sizeDotState)
   write(6,'(a32,1x,7(i8,1x),/)') 'constitutive_sizePostResults: ', shape(constitutive_sizePostResults)
   write(6,'(a32,1x,7(i8,1x))')   'maxSizeState:       ', constitutive_maxSizeState
   write(6,'(a32,1x,7(i8,1x))')   'maxSizeDotState:    ', constitutive_plasticity_maxSizeDotState
   write(6,'(a32,1x,7(i8,1x))')   'maxSizePostResults: ', constitutive_plasticity_maxSizePostResults
 endif
 flush(6)
#endif


end subroutine constitutive_init


!--------------------------------------------------------------------------------------------------
!> @brief returns the homogenize elasticity matrix
!--------------------------------------------------------------------------------------------------
function constitutive_homogenizedC(ipc,ip,el)
 use prec, only: &
   pReal
 use material, only: &
   phase_plasticity, &
   material_phase, &
   PLASTICITY_DISLOTWIN_ID, &
   PLASTICITY_DISLOUCLA_ID
 use plastic_dislotwin, only: &
   plastic_dislotwin_homogenizedC
 use lattice, only: &
   lattice_C66
 use material, only: &
   material_phase, &
   material_homog, &
   phase_NstiffnessModifiers, &
   phase_stiffnessModifier

 implicit none
 real(pReal), dimension(6,6) :: constitutive_homogenizedC
 integer(pInt), intent(in) :: &
   ipc, &                                                                                            !< component-ID of integration point
   ip, &                                                                                             !< integration point
   el 
 integer(pInt) :: &
   ho, d                                                                                                 !< element

 plasticityType: select case (phase_plasticity(material_phase(ipc,ip,el)))
   case (PLASTICITY_DISLOTWIN_ID) plasticityType
     constitutive_homogenizedC = plastic_dislotwin_homogenizedC(ipc,ip,el)
   case default plasticityType
     constitutive_homogenizedC = lattice_C66(1:6,1:6,material_phase (ipc,ip,el))
 end select plasticityType

 ho = material_homog(ip,el)
 ModifierLoop: do d = 1_pInt, phase_NstiffnessModifiers(material_phase(ipc,ip,el))
   ModifierType: select case(phase_stiffnessModifier(d,material_phase(ipc,ip,el)))

   end select ModifierType
 enddo ModifierLoop

end function constitutive_homogenizedC

!--------------------------------------------------------------------------------------------------
!> @brief calls microstructure function of the different constitutive models
!--------------------------------------------------------------------------------------------------
subroutine constitutive_microstructure(orientations, Fe, Fp, ipc, ip, el)
 use prec, only: &
   pReal
 use material, only: &
   phase_plasticity, &
   material_phase, &
   material_homog, &
   temperature, &
   thermalMapping, &
   PLASTICITY_dislotwin_ID, &
   PLASTICITY_disloucla_ID, &
   PLASTICITY_nonlocal_ID
 use plastic_nonlocal, only: &
   plastic_nonlocal_microstructure
 use plastic_dislotwin, only: &
   plastic_dislotwin_microstructure
 use plastic_disloucla, only: &
   plastic_disloucla_microstructure

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),   intent(in), dimension(3,3) :: &
   Fe, &                                                                                            !< elastic deformation gradient
   Fp                                                                                               !< plastic deformation gradient
 integer(pInt) :: &
   ho, &                                                                                            !< homogenization
   tme                                                                                              !< thermal member position
 real(pReal),   intent(in), dimension(:,:,:,:) :: &
   orientations                                                                                     !< crystal orientations as quaternions

 ho = material_homog(ip,el)
 tme = thermalMapping(ho)%p(ip,el)

 plasticityType: select case (phase_plasticity(material_phase(ipc,ip,el)))
   case (PLASTICITY_DISLOTWIN_ID) plasticityType
     call plastic_dislotwin_microstructure(temperature(ho)%p(tme),ipc,ip,el)
   case (PLASTICITY_DISLOUCLA_ID) plasticityType
     call plastic_disloucla_microstructure(temperature(ho)%p(tme),ipc,ip,el)
   case (PLASTICITY_NONLOCAL_ID) plasticityType
     call plastic_nonlocal_microstructure (Fe,Fp,ip,el)
 end select plasticityType

end subroutine constitutive_microstructure


!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient
!--------------------------------------------------------------------------------------------------
subroutine constitutive_LpAndItsTangent(Lp, dLp_dTstar3333, dLp_dFi3333, Tstar_v, Fi, ipc, ip, el)
 use prec, only: &
   pReal
 use math, only: &
   math_mul33x33, &
   math_Mandel6to33, &
   math_Mandel33to6, &
   math_Plain99to3333
 use material, only: &
   phase_plasticity, &
   material_phase, &
   material_homog, &
   temperature, &
   thermalMapping, &
   PLASTICITY_NONE_ID, &
   PLASTICITY_ISOTROPIC_ID, &
   PLASTICITY_PHENOPOWERLAW_ID, &
   PLASTICITY_DISLOTWIN_ID, &
   PLASTICITY_DISLOUCLA_ID, &
   PLASTICITY_NONLOCAL_ID
 use plastic_isotropic, only: &
   plastic_isotropic_LpAndItsTangent
 use plastic_phenopowerlaw, only: &
   plastic_phenopowerlaw_LpAndItsTangent
 use plastic_dislotwin, only: &
   plastic_dislotwin_LpAndItsTangent
 use plastic_disloucla, only: &
   plastic_disloucla_LpAndItsTangent
 use plastic_nonlocal, only: &
   plastic_nonlocal_LpAndItsTangent

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),   intent(in),  dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola-Kirchhoff stress
 real(pReal),   intent(in),  dimension(3,3) :: &
   Fi                                                                                               !< intermediate deformation gradient
 real(pReal),   intent(out), dimension(3,3) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal),   intent(out), dimension(3,3,3,3) :: &
   dLp_dTstar3333, &                                                                                !< derivative of Lp with respect to Tstar (4th-order tensor)
   dLp_dFi3333                                                                                      !< derivative of Lp with respect to Fi (4th-order tensor)
 real(pReal), dimension(6) :: &
   Mstar_v                                                                                          !< Mandel stress work conjugate with Lp
 real(pReal), dimension(9,9) :: &
   dLp_dMstar                                                                                       !< derivative of Lp with respect to Mstar (4th-order tensor)
 real(pReal), dimension(3,3) :: &
   temp_33
 integer(pInt) :: &
   ho, &                                                                                            !< homogenization
   tme                                                                                              !< thermal member position
 integer(pInt) :: &
   i, j

 ho = material_homog(ip,el)
 tme = thermalMapping(ho)%p(ip,el)

 Mstar_v = math_Mandel33to6(math_mul33x33(math_mul33x33(transpose(Fi),Fi),math_Mandel6to33(Tstar_v)))

 plasticityType: select case (phase_plasticity(material_phase(ipc,ip,el)))
   case (PLASTICITY_NONE_ID) plasticityType
     Lp = 0.0_pReal
     dLp_dMstar = 0.0_pReal
   case (PLASTICITY_ISOTROPIC_ID) plasticityType
     call plastic_isotropic_LpAndItsTangent(Lp,dLp_dMstar,Mstar_v,ipc,ip,el)
   case (PLASTICITY_PHENOPOWERLAW_ID) plasticityType
     call plastic_phenopowerlaw_LpAndItsTangent(Lp,dLp_dMstar,Mstar_v,ipc,ip,el)
   case (PLASTICITY_NONLOCAL_ID) plasticityType
     call plastic_nonlocal_LpAndItsTangent(Lp,dLp_dMstar,Mstar_v, &
                                           temperature(ho)%p(tme),ip,el)
   case (PLASTICITY_DISLOTWIN_ID) plasticityType
     call plastic_dislotwin_LpAndItsTangent(Lp,dLp_dMstar,Mstar_v, &
                                            temperature(ho)%p(tme),ipc,ip,el)
   case (PLASTICITY_DISLOUCLA_ID) plasticityType
     call plastic_disloucla_LpAndItsTangent(Lp,dLp_dMstar,Mstar_v, &
                                            temperature(ho)%p(tme), ipc,ip,el)
 end select plasticityType

 dLp_dTstar3333 = math_Plain99to3333(dLp_dMstar)
 temp_33 = math_mul33x33(Fi,math_Mandel6to33(Tstar_v))
 forall(i = 1_pInt:3_pInt, j = 1_pInt:3_pInt) &
   dLp_dFi3333(i,j,1:3,1:3) = math_mul33x33(temp_33,transpose(dLp_dTstar3333(i,j,1:3,1:3))) + &
                              math_mul33x33(math_mul33x33(Fi,dLp_dTstar3333(i,j,1:3,1:3)),math_Mandel6to33(Tstar_v))

 temp_33 = math_mul33x33(transpose(Fi),Fi)

 forall(i = 1_pInt:3_pInt, j = 1_pInt:3_pInt) &
   dLp_dTstar3333(i,j,1:3,1:3) = math_mul33x33(temp_33,dLp_dTstar3333(i,j,1:3,1:3))

end subroutine constitutive_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient
!--------------------------------------------------------------------------------------------------
subroutine constitutive_LiAndItsTangent(Li, dLi_dTstar3333, dLi_dFi3333, Tstar_v, Fi, ipc, ip, el)
 use prec, only: &
   pReal
 use math, only: &
   math_I3, &
   math_inv33, &
   math_det33, &
   math_mul33x33
 use material, only: &
   phase_plasticity, &
   material_phase, &
   phase_kinematics, &
   phase_Nkinematics, &
   PLASTICITY_isotropic_ID, &
   KINEMATICS_cleavage_opening_ID, &
   KINEMATICS_slipplane_opening_ID, &
   KINEMATICS_thermal_expansion_ID
 use plastic_isotropic, only: &
   plastic_isotropic_LiAndItsTangent
 use kinematics_cleavage_opening, only: &
   kinematics_cleavage_opening_LiAndItsTangent
 use kinematics_slipplane_opening, only: &
   kinematics_slipplane_opening_LiAndItsTangent
 use kinematics_thermal_expansion, only: &
   kinematics_thermal_expansion_LiAndItsTangent

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),   intent(in),  dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola-Kirchhoff stress
 real(pReal),   intent(in),  dimension(3,3) :: &
   Fi                                                                                               !< intermediate deformation gradient
 real(pReal),   intent(out), dimension(3,3) :: &
   Li                                                                                               !< intermediate velocity gradient
 real(pReal),   intent(out), dimension(3,3,3,3) :: &
   dLi_dTstar3333, &                                                                                !< derivative of Li with respect to Tstar (4th-order tensor)
   dLi_dFi3333
 real(pReal), dimension(3,3) :: &
   my_Li                                                                                            !< intermediate velocity gradient
 real(pReal), dimension(3,3,3,3) :: &
   my_dLi_dTstar
 real(pReal), dimension(3,3) :: &
   FiInv, &
   temp_33
 real(pReal) :: &
   detFi
 integer(pInt) :: &
   k                                                                                                !< counter in kinematics loop
 integer(pInt) :: &
   i, j

 Li = 0.0_pReal
 dLi_dTstar3333  = 0.0_pReal
 dLi_dFi3333     = 0.0_pReal

 plasticityType: select case (phase_plasticity(material_phase(ipc,ip,el)))
   case (PLASTICITY_isotropic_ID) plasticityType
     call plastic_isotropic_LiAndItsTangent(my_Li, my_dLi_dTstar, Tstar_v, ipc, ip, el)
   case default plasticityType
     my_Li = 0.0_pReal
     my_dLi_dTstar = 0.0_pReal
 end select plasticityType

 Li = Li + my_Li
 dLi_dTstar3333 = dLi_dTstar3333 + my_dLi_dTstar

 KinematicsLoop: do k = 1_pInt, phase_Nkinematics(material_phase(ipc,ip,el))
   kinematicsType: select case (phase_kinematics(k,material_phase(ipc,ip,el)))
     case (KINEMATICS_cleavage_opening_ID) kinematicsType
       call kinematics_cleavage_opening_LiAndItsTangent(my_Li, my_dLi_dTstar, Tstar_v, ipc, ip, el)
     case (KINEMATICS_slipplane_opening_ID) kinematicsType
       call kinematics_slipplane_opening_LiAndItsTangent(my_Li, my_dLi_dTstar, Tstar_v, ipc, ip, el)
     case (KINEMATICS_thermal_expansion_ID) kinematicsType
       call kinematics_thermal_expansion_LiAndItsTangent(my_Li, my_dLi_dTstar, ipc, ip, el)
     case default kinematicsType
       my_Li = 0.0_pReal
       my_dLi_dTstar = 0.0_pReal
   end select kinematicsType
   Li = Li + my_Li
   dLi_dTstar3333 = dLi_dTstar3333 + my_dLi_dTstar
 enddo KinematicsLoop

 FiInv = math_inv33(Fi)
 detFi = math_det33(Fi)
 Li = math_mul33x33(math_mul33x33(Fi,Li),FiInv)*detFi                                               !< push forward to intermediate configuration
 temp_33 = math_mul33x33(FiInv,Li)
 forall(i = 1_pInt:3_pInt, j = 1_pInt:3_pInt)
   dLi_dTstar3333(1:3,1:3,i,j) = math_mul33x33(math_mul33x33(Fi,dLi_dTstar3333(1:3,1:3,i,j)),FiInv)*detFi
   dLi_dFi3333   (1:3,1:3,i,j) = dLi_dFi3333(1:3,1:3,i,j) + Li*FiInv(j,i)
   dLi_dFi3333   (1:3,i,1:3,j) = dLi_dFi3333(1:3,i,1:3,j) + math_I3*temp_33(j,i) + Li*FiInv(j,i)
 end forall

end subroutine constitutive_LiAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief  collects initial intermediate deformation gradient
!--------------------------------------------------------------------------------------------------
pure function constitutive_initialFi(ipc, ip, el)
 use prec, only: &
   pReal
 use math, only: &
   math_I3, &
   math_inv33, &
   math_mul33x33
 use material, only: &
   phase_kinematics, &
   phase_Nkinematics, &
   material_phase, &
   KINEMATICS_thermal_expansion_ID
 use kinematics_thermal_expansion, only: &
   kinematics_thermal_expansion_initialStrain

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(3,3) :: &
   constitutive_initialFi                                                                           !< composite initial intermediate deformation gradient
 integer(pInt) :: &
   k                                                                                                !< counter in kinematics loop

 constitutive_initialFi = math_I3

 KinematicsLoop: do k = 1_pInt, phase_Nkinematics(material_phase(ipc,ip,el))                        !< Warning: small initial strain assumption
   kinematicsType: select case (phase_kinematics(k,material_phase(ipc,ip,el)))
     case (KINEMATICS_thermal_expansion_ID) kinematicsType
       constitutive_initialFi = &
         constitutive_initialFi + kinematics_thermal_expansion_initialStrain(ipc, ip, el)
   end select kinematicsType
 enddo KinematicsLoop

end function constitutive_initialFi


!--------------------------------------------------------------------------------------------------
!> @brief returns the 2nd Piola-Kirchhoff stress tensor and its tangent with respect to
!> the elastic deformation gradient depending on the selected elastic law (so far no case switch
!! because only hooke is implemented
!--------------------------------------------------------------------------------------------------
subroutine constitutive_TandItsTangent(T, dT_dFe, dT_dFi, Fe, Fi, ipc, ip, el)
 use prec, only: &
   pReal

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),   intent(in),  dimension(3,3) :: &
   Fe, &                                                                                            !< elastic deformation gradient
   Fi                                                                                               !< intermediate deformation gradient
 real(pReal),   intent(out), dimension(3,3) :: &
   T                                                                                                !< 2nd Piola-Kirchhoff stress tensor
 real(pReal),   intent(out), dimension(3,3,3,3) :: &
   dT_dFe, &                                                                                        !< derivative of 2nd P-K stress with respect to elastic deformation gradient
   dT_dFi                                                                                           !< derivative of 2nd P-K stress with respect to intermediate deformation gradient

 call constitutive_hooke_TandItsTangent(T, dT_dFe, dT_dFi, Fe, Fi, ipc, ip, el)


end subroutine constitutive_TandItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief returns the 2nd Piola-Kirchhoff stress tensor and its tangent with respect to
!> the elastic deformation gradient using hookes law
!--------------------------------------------------------------------------------------------------
subroutine constitutive_hooke_TandItsTangent(T, dT_dFe, dT_dFi, Fe, Fi, ipc, ip, el)
 use prec, only: &
   pReal
 use math, only : &
   math_mul33x33, &
   math_transpose33, &
   math_mul3333xx33, &
   math_mul3333xx3333, &
   math_tensorcomp3333, &
   math_tensorcomptransp3333, &
   math_Mandel66to3333, &
   math_I3

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),   intent(in),  dimension(3,3) :: &
   Fe, &                                                                                            !< elastic deformation gradient
   Fi                                                                                               !< intermediate deformation gradient
 real(pReal),   intent(out), dimension(3,3) :: &
   T                                                                                                !< 2nd Piola-Kirchhoff stress tensor in lattice configuration
 real(pReal),   intent(out), dimension(3,3,3,3) :: &
   dT_dFe, &                                                                                        !< derivative of 2nd P-K stress with respect to elastic deformation gradient
   dT_dFi
 real(pReal), dimension(3,3) :: E
 real(pReal), dimension(3,3,3,3) :: C, dEi_dFe, dEi_dFi

 C = math_Mandel66to3333(constitutive_homogenizedC(ipc,ip,el))
 E = 0.5_pReal*(math_mul33x33(math_transpose33(Fe),Fe)-math_I3)                                            !< Green-Lagrange strain in unloaded configuration
 T = math_mul3333xx33(C,math_mul33x33(math_mul33x33(transpose(Fi),E),Fi))                           !< 2PK stress in lattice configuration in work conjugate with GL strain pulled back to lattice configuration

 dEi_dFe = (math_tensorcomp3333(math_transpose33(math_mul33x33(Fe,Fi)),Fi) + &
            math_tensorcomptransp3333(math_transpose33(Fi),math_mul33x33(Fe,Fi)))/2.0_pReal
 dEi_dFi = (math_tensorcomp3333(math_transpose33(math_mul33x33(E,Fi)),math_I3) + &
            math_tensorcomptransp3333(math_I3,math_mul33x33(E,Fi)))/2.0_pReal 
 dT_dFe = math_mul3333xx3333(C,dEi_dFe)
 dT_dFi = math_mul3333xx3333(C,dEi_dFi)

end subroutine constitutive_hooke_TandItsTangent

!--------------------------------------------------------------------------------------------------
!> @brief contains the constitutive equation for calculating the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine constitutive_collectDotState(Tstar_v, FeArray, FpArray, subdt, subfracArray,ipc, ip, el)
 use prec, only: &
   pReal, &
   pLongInt
 use debug, only: &
   debug_cumDotStateCalls, &
   debug_cumDotStateTicks, &
   debug_level, &
   debug_constitutive, &
   debug_levelBasic
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   phase_plasticity, &
   phase_chemicalFE, &
   phase_source, &
   phase_Nsources, &
   material_phase, &
   material_homog, &
   temperature, &
   thermalMapping, &
   homogenization_maxNgrains, &
   PLASTICITY_none_ID, &
   PLASTICITY_isotropic_ID, &
   PLASTICITY_phenopowerlaw_ID, &
   PLASTICITY_dislotwin_ID, &
   PLASTICITY_disloucla_ID, &
   PLASTICITY_nonlocal_ID, &
   CHEMICALFE_quadenergy_ID, &
   CHEMICALFE_thermodynamic_ID, &
   SOURCE_thermal_externalheat_ID
 use plastic_isotropic, only:  &
   plastic_isotropic_dotState
 use plastic_phenopowerlaw, only: &
   plastic_phenopowerlaw_dotState
 use plastic_dislotwin, only: &
   plastic_dislotwin_dotState
 use plastic_disloucla, only: &
   plastic_disloucla_dotState
 use plastic_nonlocal, only: &
   plastic_nonlocal_dotState
 use chemicalFE_quadenergy, only: &
   chemicalFE_quadenergy_dotState  
 use chemicalFE_thermodynamic, only: &
   chemicalFE_thermodynamic_dotState  
 use source_thermal_externalheat, only: &
   source_thermal_externalheat_dotState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),  intent(in) :: &
   subdt                                                                                            !< timestep
 real(pReal),  intent(in), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   subfracArray                                                                                     !< subfraction of timestep
 real(pReal),  intent(in), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   FeArray, &                                                                                       !< elastic deformation gradient
   FpArray                                                                                          !< plastic deformation gradient
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)
 integer(pLongInt) :: &
   tick = 0_pLongInt, &
   tock = 0_pLongInt, &
   tickrate, &
   maxticks
 integer(pInt) :: &
   ho, &                                                                                            !< homogenization
   tme, &                                                                                           !< thermal member position
   s                                                                                                !< counter in source loop

 if (iand(debug_level(debug_constitutive), debug_levelBasic) /= 0_pInt) &
   call system_clock(count=tick,count_rate=tickrate,count_max=maxticks)

 ho  = material_homog(    ip,el)
 tme = thermalMapping(ho)%p(ip,el)

 plasticityType: select case (phase_plasticity(material_phase(ipc,ip,el)))
   case (PLASTICITY_ISOTROPIC_ID) plasticityType
     call plastic_isotropic_dotState    (Tstar_v,ipc,ip,el)
   case (PLASTICITY_PHENOPOWERLAW_ID) plasticityType
     call plastic_phenopowerlaw_dotState(Tstar_v,ipc,ip,el)
   case (PLASTICITY_DISLOTWIN_ID) plasticityType
     call plastic_dislotwin_dotState    (Tstar_v,temperature(ho)%p(tme), &
                                         ipc,ip,el)
   case (PLASTICITY_DISLOUCLA_ID) plasticityType
     call plastic_disloucla_dotState    (Tstar_v,temperature(ho)%p(tme), &
                                         ipc,ip,el)
   case (PLASTICITY_NONLOCAL_ID) plasticityType
     call plastic_nonlocal_dotState     (Tstar_v,FeArray,FpArray,temperature(ho)%p(tme), &
                                         subdt,subfracArray,ip,el)
 end select plasticityType

 chemicalType: select case (phase_chemicalFE(material_phase(ipc,ip,el)))
   case (CHEMICALFE_quadenergy_ID) chemicalType
     call chemicalFE_quadenergy_dotState   (ipc,ip,el)
   case (CHEMICALFE_thermodynamic_ID) chemicalType
     call chemicalFE_thermodynamic_dotState(ipc,ip,el)
 end select chemicalType

 SourceLoop: do s = 1_pInt, phase_Nsources(material_phase(ipc,ip,el))
    sourceType: select case (phase_source(s,material_phase(ipc,ip,el)))
     case (SOURCE_thermal_externalheat_ID) sourceType
       call source_thermal_externalheat_dotState(         ipc, ip, el)
   end select sourceType
 enddo SourceLoop

 if (iand(debug_level(debug_constitutive), debug_levelBasic) /= 0_pInt) then
   call system_clock(count=tock,count_rate=tickrate,count_max=maxticks)
   !$OMP CRITICAL (debugTimingDotState)
   debug_cumDotStateCalls = debug_cumDotStateCalls + 1_pInt
   debug_cumDotStateTicks = debug_cumDotStateTicks + tock-tick
   !$OMP FLUSH (debug_cumDotStateTicks)
   if (tock < tick) debug_cumDotStateTicks  = debug_cumDotStateTicks + maxticks
   !$OMP END CRITICAL (debugTimingDotState)
 endif
end subroutine constitutive_collectDotState

!--------------------------------------------------------------------------------------------------
!> @brief for constitutive models having an instantaneous change of state
!> will return false if delta state is not needed/supported by the constitutive model
!--------------------------------------------------------------------------------------------------
subroutine constitutive_collectDeltaState(Tstar_v, Fe, ipc, ip, el)
 use prec, only: &
   pReal, &
   pLongInt
 use debug, only: &
   debug_cumDeltaStateCalls, &
   debug_cumDeltaStateTicks, &
   debug_level, &
   debug_constitutive, &
   debug_levelBasic
 use material, only: &
   phase_plasticity, &
   material_phase, &
   PLASTICITY_NONLOCAL_ID
 use plastic_nonlocal, only: &
   plastic_nonlocal_deltaState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),   intent(in),  dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola-Kirchhoff stress
 real(pReal),   intent(in), dimension(3,3) :: &
   Fe                                                                                               !< elastic deformation gradient
 integer(pLongInt) :: &
   tick, tock, &
   tickrate, &
   maxticks

 if (iand(debug_level(debug_constitutive), debug_levelBasic) /= 0_pInt) &
   call system_clock(count=tick,count_rate=tickrate,count_max=maxticks)

 if(phase_plasticity(material_phase(ipc,ip,el)) == PLASTICITY_NONLOCAL_ID) &
   call plastic_nonlocal_deltaState(Tstar_v,ip,el)

 if (iand(debug_level(debug_constitutive), debug_levelBasic) /= 0_pInt) then
   call system_clock(count=tock,count_rate=tickrate,count_max=maxticks)
   !$OMP CRITICAL (debugTimingDeltaState)
     debug_cumDeltaStateCalls = debug_cumDeltaStateCalls + 1_pInt
     debug_cumDeltaStateTicks = debug_cumDeltaStateTicks + tock-tick
     !$OMP FLUSH (debug_cumDeltaStateTicks)
     if (tock < tick) debug_cumDeltaStateTicks  = debug_cumDeltaStateTicks + maxticks
   !$OMP END CRITICAL (debugTimingDeltaState)
 endif

end subroutine constitutive_collectDeltaState


!--------------------------------------------------------------------------------------------------
!> @brief returns array of constitutive results
!--------------------------------------------------------------------------------------------------
function constitutive_postResults(Tstar_v, FeArray, ipc, ip, el)
 use prec, only: &
   pReal
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   plasticState, &
   chemicalState, &
   sourceState, &
   phase_plasticity, &
   phase_chemicalFE, &
   phase_source, &
   phase_Nsources, &
   material_phase, &
   material_homog, &
   temperature, &
   thermalMapping, &
   homogenization_maxNgrains, &
   PLASTICITY_NONE_ID, &
   PLASTICITY_ISOTROPIC_ID, &
   PLASTICITY_PHENOPOWERLAW_ID, &
   PLASTICITY_DISLOTWIN_ID, &
   PLASTICITY_DISLOUCLA_ID, &
   PLASTICITY_NONLOCAL_ID, &
   CHEMICALFE_quadenergy_ID, &
   CHEMICALFE_thermodynamic_ID
 use plastic_isotropic, only: &
   plastic_isotropic_postResults
 use plastic_phenopowerlaw, only: &
   plastic_phenopowerlaw_postResults
 use plastic_dislotwin, only: &
   plastic_dislotwin_postResults
 use plastic_disloucla, only: &
   plastic_disloucla_postResults
 use plastic_nonlocal, only: &
   plastic_nonlocal_postResults
 use chemicalFE_quadenergy, only: &
   chemicalFE_quadenergy_postResults  
 use chemicalFE_thermodynamic, only: &
   chemicalFE_thermodynamic_postResults  

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(plasticState (material_phase(ipc,ip,el))%sizePostResults + &
                        chemicalState(material_phase(ipc,ip,el))%sizePostResults + &
                        sum(sourceState(material_phase(ipc,ip,el))%p(:)%sizePostResults)) :: &
   constitutive_postResults
 real(pReal),  intent(in), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   FeArray                                                                                          !< elastic deformation gradient
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)
 integer(pInt) :: &
   startPos, endPos
 integer(pInt) :: &
   ho, &                                                                                            !< homogenization
   tme, &                                                                                           !< thermal member position
   s                                                                                                !< counter in source loop

 constitutive_postResults = 0.0_pReal

 ho = material_homog(    ip,el)
 tme = thermalMapping(ho)%p(ip,el)

 startPos = 1_pInt
 endPos = plasticState(material_phase(ipc,ip,el))%sizePostResults

 plasticityType: select case (phase_plasticity(material_phase(ipc,ip,el)))
   case (PLASTICITY_ISOTROPIC_ID) plasticityType
     constitutive_postResults(startPos:endPos) = plastic_isotropic_postResults(Tstar_v,ipc,ip,el)
   case (PLASTICITY_PHENOPOWERLAW_ID) plasticityType
     constitutive_postResults(startPos:endPos) = &
       plastic_phenopowerlaw_postResults(Tstar_v,ipc,ip,el)
   case (PLASTICITY_DISLOTWIN_ID) plasticityType
     constitutive_postResults(startPos:endPos) = &
       plastic_dislotwin_postResults(Tstar_v,temperature(ho)%p(tme),ipc,ip,el)
   case (PLASTICITY_DISLOUCLA_ID) plasticityType
     constitutive_postResults(startPos:endPos) = &
       plastic_disloucla_postResults(Tstar_v,temperature(ho)%p(tme),ipc,ip,el)
   case (PLASTICITY_NONLOCAL_ID) plasticityType
     constitutive_postResults(startPos:endPos) = &
       plastic_nonlocal_postResults (Tstar_v,FeArray,ip,el)
 end select plasticityType

 startPos = endPos + 1_pInt
 endPos = endPos + chemicalState(material_phase(ipc,ip,el))%sizePostResults
 
 chemicalFEType: select case (phase_chemicalFE(material_phase(ipc,ip,el)))
   case (CHEMICALFE_quadenergy_ID) chemicalFEType
     constitutive_postResults(startPos:endPos) = &
       chemicalFE_quadenergy_postResults(ipc,ip,el)
   case (CHEMICALFE_thermodynamic_ID) chemicalFEType
     constitutive_postResults(startPos:endPos) = &
       chemicalFE_thermodynamic_postResults(ipc,ip,el)
 end select chemicalFEType

 SourceLoop: do s = 1_pInt, phase_Nsources(material_phase(ipc,ip,el))
   startPos = endPos + 1_pInt
   endPos = endPos + sourceState(material_phase(ipc,ip,el))%p(s)%sizePostResults
   sourceType: select case (phase_source(s,material_phase(ipc,ip,el)))
   end select sourceType
 enddo SourceLoop

end function constitutive_postResults

end module constitutive
