!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Mixture homogenization models
!--------------------------------------------------------------------------------------------------
module solute_flux
 use prec, only: &
   pInt, &
   pReal, &
   p_2Dvec
 
 implicit none
 private
 integer(pInt),             dimension(:,:),   allocatable, target, public :: &
   solute_flux_sizePostResult
 
 character(len=64),         dimension(:,:),   allocatable, target, public :: &
  solute_flux_output                                                                                !< name of each post result output
 integer(pInt),             dimension(:),     allocatable, target, public :: &
   solute_flux_Noutput                                                                              !< number of outputs per homog instance

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 conc_ID
 end enum

 type, private :: tParameters                                                                       !< container type for internal constitutive parameters
   integer(kind(undefined_ID)), dimension(:), allocatable :: & 
     outputID
   real(pReal),   dimension(:), allocatable :: &
     initialChemPot
   type(p_2Dvec), dimension(:), allocatable :: &
     interfaceChemPot
 end type

 type(tParameters), dimension(:), allocatable, private :: param                                     !< containers of constitutive parameters (len Ninstance)

 public :: &
   solute_flux_init, &
   solute_flux_getInitialComponentPotential, &  
   solute_flux_phase_calComponentConcandTangent, &
   solute_flux_calComponentConcandTangent, &
   solute_flux_getComponentConc, &
   solute_flux_getComponentMobility, &
   solute_flux_getEffChargeNumber, &
   solute_flux_calAndPutComponentConcRate, &
   solute_flux_postResults

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine solute_flux_init
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
 use mesh, only: &
   FE_Nips, &
   FE_geomtype, &
   mesh_NcpElems, &
   mesh_element
 use math, only: &
   math_VecMtoSymNN
 
 implicit none
 integer(pInt) :: &
   o, el, ip, gr, cp, outputsize
 integer(pInt) :: sizeState
 integer :: &
   maxNinstance, &
   homog, &
   phase, &
   instance
 integer :: &
   NofmySolute                                                                                       ! no pInt (stores a system dependen value from 'count'
 character(len=65536) :: cpStr
 character(len=65536), dimension(:), allocatable :: outputs
 character(len=65536), dimension(0), parameter :: emptyStringArray = [character(len=65536)::]
 integer(kind(undefined_ID)) :: outputID                                                           !< ID of each post result output
 
 write(6,'(/,a)')   ' <<<+-  solute_'//SOLUTE_flux_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = count(solute_type == SOLUTE_flux_ID)
 if (maxNinstance == 0) return
 
 if (iand(debug_level(debug_HOMOGENIZATION),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 allocate(solute_flux_sizePostResult(maxval(homogenization_Noutput),maxNinstance),source=0_pInt)
 allocate(solute_flux_Noutput        (maxNinstance),                              source=0_pInt)
 allocate(solute_flux_output        (maxval(homogenization_Noutput),maxNinstance))
          solute_flux_output = ''
 
 allocate(param(maxNinstance))

 parsingFile: do homog = 1_pInt, material_Nhomogenization
  if (solute_type(homog) /= SOLUTE_flux_ID) cycle
  instance = solute_typeInstance(homog)
  param(instance)%initialChemPot = &
                     config_homogenization(homog)%getFloats('initial_chemicalpotential')
  allocate(param(instance)%interfaceChemPot(homogenization_Ncomponents(homog)))
  do cp = 1_pInt, homogenization_Ncomponents(homog)
    write(cpStr,'(i0)') cp
    if (config_homogenization(homog)%keyExists('interface_chemicalpotential_'//trim(cpStr))) then
      param(instance)%interfaceChemPot(cp)%p = &
        math_VecMtoSymNN(config_homogenization(homog)%getFloats('interface_chemicalpotential_'//trim(cpStr)), &
                         homogenization_Ngrains(homog))
    else
      allocate(param(instance)%interfaceChemPot(cp)%p(0,0))
    endif
  enddo  
    
  outputs = config_homogenization(homog)%getStrings('(output)',defaultVal=emptyStringArray)
  allocate(param(instance)%outputID(0))
  do o=1_pInt, size(outputs)
    outputID = undefined_ID
    select case(outputs(o))
      case('conc','concentration')
        outputID = conc_ID
        outputSize = homogenization_Ncomponents(homog)
     end select

     if (outputID /= undefined_ID) then
       solute_flux_output(o,instance) = outputs(o)
       solute_flux_sizePostResult(o,instance) = outputSize
       param(instance)%outputID = [param(instance)%outputID , outputID]
     endif
   enddo
 enddo parsingFile

 initializeInstances: do homog = 1_pInt, material_Nhomogenization
   mySolute: if (solute_type(homog) == SOLUTE_flux_ID) then
     NofmySolute = count(material_homog == homog)
     instance = solute_typeInstance(homog)

! allocate state arrays
     sizeState = 0_pInt
     soluteState(homog)%sizeState = sizeState
     soluteState(homog)%sizePostResults = sum(solute_flux_sizePostResult(:,instance))
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
     mySolute2: if (solute_type(homog) == SOLUTE_flux_ID) then
       do gr = 1_pInt, homogenization_Ngrains(homog)
         phase = material_phase(gr,ip,el)
         if (homogenization_Ncomponents(homog) /= phase_Ncomponents(phase)) &
           call IO_error(211_pInt,el=el,ip=ip,g=gr,ext_msg='mismatch in # of homog and phase components')
       enddo
     endif mySolute2
   enddo IPLoop
 enddo elementLoop    

end subroutine solute_flux_init


!--------------------------------------------------------------------------------------------------
!> @brief return initial component chemical potential at material point 
!--------------------------------------------------------------------------------------------------
function solute_flux_getInitialComponentPotential(ip,el)
 use material, only: &
   material_homog, &
   solute_typeInstance, &
   homogenization_maxNcomponents
 
 implicit none
 integer(pInt),                                         intent(in) :: &
   ip, el                                                                                      !< element number
 real(pReal), dimension(homogenization_maxNcomponents) :: &
   solute_flux_getInitialComponentPotential

 solute_flux_getInitialComponentPotential = &
   param(solute_typeInstance(material_homog(ip,el)))%initialChemPot

end function solute_flux_getInitialComponentPotential


!--------------------------------------------------------------------------------------------------
!> @brief calculates and returns component concentrations and tangents at ipc,ip,el
!--------------------------------------------------------------------------------------------------
subroutine solute_flux_phase_calComponentConcandTangent(Conc_local,dConcdChemPot_local,dConcdGradC_local, & 
                                                  ChemPot,GradC,ipc,ip,el)
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   material_phase, &
   homogenization_Ncomponents, &
   homogenization_maxNcomponents, &
   solute_typeInstance, &
   phase_chemicalFE, &
   CHEMICALFE_quadenergy_ID, &
   CHEMICALFE_thermodynamic_ID, &
   phasefracMapping, &
   phasefrac, &
   phase_Nkinematics, &
   phase_kinematics, &
   KINEMATICS_solute_strain_ID
 use crystallite, only: &
   crystallite_Tstar_v
 use kinematics_solute_strain, only: &
   kinematics_solute_strain_getMechanicalComponentPotential
 use chemicalFE_quadenergy, only: &
   chemicalFE_quadenergy_calConcandTangent  
 use chemicalFE_thermodynamic, only: &
   chemicalFE_thermodynamic_calConcandTangent  
 
 implicit none
 integer(pInt),                                         intent(in)  :: &
   ipc, ip, el                                                                                      !< element number
 real(pReal), dimension(homogenization_maxNcomponents), intent(in)  :: &
   ChemPot, &
   GradC
 real(pReal), dimension(homogenization_maxNcomponents), intent(out) :: &
   Conc_local
 real(pReal), dimension(homogenization_maxNcomponents, &
                        homogenization_maxNcomponents), intent(out) :: &
   dConcdChemPot_local, &
   dConcdGradC_local
 real(pReal), dimension(homogenization_maxNcomponents) :: &
   MechChemPot, &
   InftChemPot
 integer(pInt) :: &
   gr, grI, grJ, k, &
   homog, &
   instance, & 
   offset

 homog = material_homog(ip,el)
 instance = solute_typeInstance(homog)
 offset = phasefracMapping(homog)%p(ip,el)
 
 InftChemPot = 0.0_pReal
 do gr = 1_pInt, homogenization_Ncomponents(homog)
   do grI = 1_pInt, size(param(instance)%interfaceChemPot(gr)%p,1)
     do grJ = grI+1_pInt, size(param(instance)%interfaceChemPot(gr)%p,2)
       InftChemPot(gr) = &
         InftChemPot(gr) + &
         param(instance)%interfaceChemPot(gr)%p(grI,grJ)* &
         phasefrac(homog)%p(grI,offset)* &
         phasefrac(homog)%p(grJ,offset)
     enddo
   enddo
 enddo        
  
 MechChemPot = 0.0_pReal
  KinematicsLoop: do k = 1_pInt, phase_Nkinematics(material_phase(ipc,ip,el))
   kinematicsType: select case (phase_kinematics(k,material_phase(ipc,ip,el)))
    case (KINEMATICS_solute_strain_ID) kinematicsType
      MechChemPot = MechChemPot + &
        kinematics_solute_strain_getMechanicalComponentPotential(crystallite_Tstar_v(1:6,ipc,ip,el), &
                                                                    ipc,ip,el)

    end select kinematicsType
  enddo KinematicsLoop

  chemicalFEType: select case (phase_chemicalFE(material_phase(ipc,ip,el)))
    case (CHEMICALFE_quadenergy_ID) chemicalFEType
      call chemicalFE_quadenergy_calConcandTangent(Conc_local,dConcdChemPot_local,dConcdGradC_local, & 
                                                    ChemPot,GradC,MechChemPot,InftChemPot,ipc,ip,el)

     case (CHEMICALFE_thermodynamic_ID) chemicalFEType
       call chemicalFE_thermodynamic_calConcandTangent(Conc_local,dConcdChemPot_local,dConcdGradC_local, & 
                                                       ChemPot,GradC,MechChemPot,InftChemPot,ipc,ip,el)

     case default
       Conc_local = 0.0_pReal
       dConcdChemPot_local = 0.0_pReal
       dConcdGradC_local = 0.0_pReal
       
  end select chemicalFEType

end subroutine solute_flux_phase_calComponentConcandTangent



!--------------------------------------------------------------------------------------------------
!> @brief calculates and returns component concentrations and tangents at material point 
!--------------------------------------------------------------------------------------------------
subroutine solute_flux_calComponentConcandTangent(Conc,dConcdChemPot,dConcdGradC, & 
                                                  ChemPot,GradC,ip,el)
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   material_phase, &
   homogenization_Ncomponents, &
   homogenization_maxNcomponents, &
   solute_typeInstance, &
   phase_chemicalFE, &
   CHEMICALFE_quadenergy_ID, &
   CHEMICALFE_thermodynamic_ID, &
   phasefracMapping, &
   phasefrac, &
   phase_Nkinematics, &
   phase_kinematics, &
   KINEMATICS_solute_strain_ID
  
 implicit none
 integer(pInt),                                         intent(in)  :: &
   ip, el                                                                                      !< element number
 real(pReal), dimension(homogenization_maxNcomponents), intent(in)  :: &
   ChemPot, &
   GradC
 real(pReal), dimension(homogenization_maxNcomponents), intent(out) :: &
   Conc
 real(pReal), dimension(homogenization_maxNcomponents, &
                        homogenization_maxNcomponents), intent(out) :: &
   dConcdChemPot, &
   dConcdGradC
 real(pReal), dimension(homogenization_maxNcomponents)  :: &
   Conc_local
 real(pReal), dimension(homogenization_maxNcomponents, &
                        homogenization_maxNcomponents)  :: &
   dConcdChemPot_local, &
   dConcdGradC_local  
 integer(pInt) :: &
   gr, grI, grJ, k, &
   homog, &
   instance, & 
   offset

 homog = material_homog(ip,el)
 instance = solute_typeInstance(homog)
 offset = phasefracMapping(homog)%p(ip,el)
 Conc = 0.0_pReal
 dConcdChemPot = 0.0_pReal
 dConcdGradC = 0.0_pReal

 do gr = 1_pInt, homogenization_Ngrains(homog)
   call solute_flux_phase_calComponentConcandTangent(Conc_local,dConcdChemPot_local,dConcdGradC_local, & 
                                                  ChemPot,GradC,gr,ip,el)

   Conc = Conc + phasefrac(homog)%p(gr,offset)*Conc_local
   dConcdChemPot = dConcdChemPot + phasefrac(homog)%p(gr,offset)*dConcdChemPot_local
   dConcdGradC = dConcdGradC + phasefrac(homog)%p(gr,offset)*dConcdGradC_local
 enddo

end subroutine solute_flux_calComponentConcandTangent


!--------------------------------------------------------------------------------------------------
!> @brief return previously calculated component concentration at material point 
!--------------------------------------------------------------------------------------------------
function solute_flux_getComponentConc(ip,el)
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   material_phase, &
   homogenization_maxNcomponents, &
   phasefracMapping, &
   phasefrac, &
   chemConcMapping, &
   chemicalConc
 
 implicit none
 integer(pInt),                                         intent(in) :: &
   ip, el                                                                                      !< element number
 real(pReal), dimension(homogenization_maxNcomponents) :: &
   solute_flux_getComponentConc
 integer(pInt) :: &
   gr, cp, &
   homog, &
   phase

 homog = material_homog(ip,el)
 solute_flux_getComponentConc = 0.0_pReal
 do gr = 1_pInt, homogenization_Ngrains(homog)
   phase = material_phase(gr,ip,el)
   do cp = 1_pInt, homogenization_maxNcomponents
     solute_flux_getComponentConc(cp) = &
       solute_flux_getComponentConc(cp) + &
       phasefrac(homog)%p(gr,phasefracMapping(homog)%p(ip,el))* &
       chemicalConc(phase)%p(cp,chemConcMapping(phase)%p(gr,ip,el))
   enddo
 enddo

end function solute_flux_getComponentConc


!--------------------------------------------------------------------------------------------------
!> @brief return component mobility at material point 
!--------------------------------------------------------------------------------------------------
function solute_flux_getComponentMobility(ip,el)
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   material_phase, &
   homogenization_maxNcomponents, &
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
 integer(pInt),                             intent(in) :: &
   ip, el                                                                                      !< element number
 real(pReal), dimension(homogenization_maxNcomponents) :: &
   solute_flux_getComponentMobility
 integer(pInt) :: &
   gr, &
   homog, & 
   offset

 homog = material_homog(ip,el)
 offset = phasefracMapping(homog)%p(ip,el)
 solute_flux_getComponentMobility = 0.0_pReal
 do gr = 1_pInt, homogenization_Ngrains(homog)
   chemicalFEType: select case (phase_chemicalFE(material_phase(gr,ip,el)))
     case (CHEMICALFE_quadenergy_ID) chemicalFEType
       solute_flux_getComponentMobility = &
         solute_flux_getComponentMobility + &
         phasefrac(homog)%p(gr,offset)* &
         chemicalFE_quadenergy_getMobility(gr,ip,el)

     case (CHEMICALFE_thermodynamic_ID) chemicalFEType
       solute_flux_getComponentMobility = &
         solute_flux_getComponentMobility + &
         phasefrac(homog)%p(gr,offset)* &
         chemicalFE_thermodynamic_getMobility(gr,ip,el)

   end select chemicalFEType
 enddo

end function solute_flux_getComponentMobility


!--------------------------------------------------------------------------------------------------
!> @brief return component mobility at material point 
!--------------------------------------------------------------------------------------------------
function solute_flux_getEffChargeNumber(ip,el)
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   material_phase, &
   homogenization_maxNcomponents, &
   phase_chemicalFE, &
   CHEMICALFE_quadenergy_ID, &
   CHEMICALFE_thermodynamic_ID, &
   phasefracMapping, &
   phasefrac
 use chemicalFE_quadenergy, only: &
   chemicalFE_quadenergy_getEffChargeNumber  
 use chemicalFE_thermodynamic, only: &
   chemicalFE_thermodynamic_getEffChargeNumber  
 
 implicit none
 integer(pInt),                             intent(in) :: &
   ip, el                                                                                      !< element number
 real(pReal), dimension(homogenization_maxNcomponents) :: &
   solute_flux_getEffChargeNumber
 integer(pInt) :: &
   gr, &
   homog, & 
   offset

 homog = material_homog(ip,el)
 offset = phasefracMapping(homog)%p(ip,el)
 solute_flux_getEffChargeNumber = 0.0_pReal
 do gr = 1_pInt, homogenization_Ngrains(homog)
   chemicalFEType: select case (phase_chemicalFE(material_phase(gr,ip,el)))
     case (CHEMICALFE_quadenergy_ID) chemicalFEType
       solute_flux_getEffChargeNumber = &
         solute_flux_getEffChargeNumber + &
         phasefrac(homog)%p(gr,offset)* &
         chemicalFE_quadenergy_getEffChargeNumber(gr,ip,el)

     case (CHEMICALFE_thermodynamic_ID) chemicalFEType
       solute_flux_getEffChargeNumber = &
         solute_flux_getEffChargeNumber + &
         phasefrac(homog)%p(gr,offset)* &
         chemicalFE_thermodynamic_getEffChargeNumber(gr,ip,el)

   end select chemicalFEType
 enddo

end function solute_flux_getEffChargeNumber


!--------------------------------------------------------------------------------------------------
!> @brief set module wide component concentrations 
!--------------------------------------------------------------------------------------------------
subroutine solute_flux_calAndPutComponentConcRate(ChemPot,GradC,dt,ip,el)
 use material, only: &
   material_homog, &
   homogenization_Ngrains, &
   material_phase, &
   homogenization_maxNcomponents, &
   chemicalConcRate, &
   chemicalConc0, &
   chemConcMapping
 implicit none
 integer(pInt),                                         intent(in) :: &
   ip, el                                                                                      !< element number
 real(pReal),                                           intent(in) :: &
   dt
 real(pReal), dimension(homogenization_maxNcomponents), intent(in) :: &
   ChemPot, &
   GradC
 real(pReal), dimension(homogenization_maxNcomponents) :: &
   Conc_local
 real(pReal), dimension(homogenization_maxNcomponents, &
                        homogenization_maxNcomponents) :: &
   dConcdChemPot_local, &
   dConcdGradC_local
 integer(pInt) :: &
   gr, &
   cp, &
   phase, &
   homog

 homog = material_homog(ip,el)
 do gr = 1_pInt, homogenization_Ngrains(homog)
   phase = material_phase(gr,ip,el)
   
   call solute_flux_phase_calComponentConcandTangent(Conc_local,dConcdChemPot_local,dConcdGradC_local, & 
                                                  ChemPot,GradC,gr,ip,el)
   
   do cp = 1_pInt, homogenization_maxNcomponents
     chemicalConcRate(phase)%p(cp,chemConcMapping(phase)%p(gr,ip,el)) = &
       (Conc_local(cp) - chemicalConc0(phase)%p(cp,chemConcMapping(phase)%p(gr,ip,el)))/dt
   enddo
 enddo

end subroutine solute_flux_calAndPutComponentConcRate


!--------------------------------------------------------------------------------------------------
!> @brief return array of homogenization results for post file inclusion 
!--------------------------------------------------------------------------------------------------
function solute_flux_postResults(ip,el)
 use material, only: &
   material_homog, &
   homogenization_Ncomponents, &
   solute_typeInstance, &
   homogenization_Noutput
 
 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),  dimension(sum(solute_flux_sizePostResult(:,solute_typeInstance(material_homog(ip,el))))) :: &
   solute_flux_postResults
 
 integer(pInt) :: &
   homID, &
   o, c
   
 c = 0_pInt
 homID = solute_typeInstance(material_homog(ip,el))
 solute_flux_postResults = 0.0_pReal
 
 do o = 1_pInt,homogenization_Noutput(material_homog(ip,el))
   select case(param(homID)%outputID(o))
     case (conc_ID)
       solute_flux_postResults(c+1_pInt:c+homogenization_Ncomponents(material_homog(ip,el))) = &
         solute_flux_getComponentConc(ip,el)
       c = c + homogenization_Ncomponents(material_homog(ip,el))

   end select
 enddo

end function solute_flux_postResults

end module solute_flux
