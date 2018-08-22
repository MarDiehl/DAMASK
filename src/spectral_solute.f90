!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Multiphase PETSc solver
!--------------------------------------------------------------------------------------------------
module spectral_solute
#include <petsc/finclude/petscsnes.h>
#include <petsc/finclude/petscdmda.h>
 use PETScdmda
 use PETScsnes
 use prec, only: & 
   pInt, &
   pReal
 use spectral_utilities, only: &
   tSolutionState, &
   tSolutionParams

 implicit none
 private

!--------------------------------------------------------------------------------------------------
real(pReal),                                          parameter,           private :: &
   electronic_charge            = 1.60217662e-19_pReal                                                                          !< electronic charge in coulombs

! derived types
 type(tSolutionParams), private :: params

!--------------------------------------------------------------------------------------------------
! PETSc data
 SNES, private :: solute_snes
 Vec,  private :: solution_current, solution_lastInc

!--------------------------------------------------------------------------------------------------
! common pointwise data
 integer(pInt)                     , private              :: Ncomponents
 real   (pReal)                    , private              :: delta(3), BMat(3,4)
 real   (pReal), dimension(:,:,:,:), private, allocatable :: conc_current, conc_stagInc, conc_lastInc 
 PetscInt                          , private              :: xstart,ystart,zstart,xend,yend,zend

 public :: &
   spectral_solute_init, &
   spectral_solute_solution, &
   spectral_solute_forward
 external :: &
   PETScErrorF

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields and fills them with data, potentially from restart info
!--------------------------------------------------------------------------------------------------
subroutine spectral_solute_init
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use IO, only: &
   IO_error, &
   IO_intOut, &
   IO_read_realFile, &
   IO_timeStamp
 use numerics, only: &
   worldrank, &
   worldsize
 use mesh, only: &
   geomSize, &
   grid, &
   grid3
 use material, only: &
   homogenization_Ncomponents, &
   homogenization_maxNcomponents, &
   homogenization_active, &
   material_homog   
 use config, only: &  
   material_Nhomogenization
 use solute_flux, only: &
   solute_flux_getInitialComponentPotential, &
   solute_flux_getComponentConc
   
 implicit none
 DM :: solute_grid
 PetscScalar, pointer :: x_scal(:,:,:,:)
 PetscErrorCode       :: ierr
 integer(pInt), dimension(:), allocatable :: localK  
 integer(pInt) :: proc
 integer(pInt) :: i, j, k, cell, homog

 write(6,'(/,a)') ' <<<+-  spectral_solute init  -+>>>'
 write(6,'(a15,a)')   ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
 do homog = 1_pInt, material_Nhomogenization
   if (homogenization_active(homog)) then
     if (homogenization_Ncomponents(homog) /= homogenization_maxNcomponents) &
       call IO_error(211_pInt,ext_msg='# of components mismatch in spectral_solute')
   endif
 enddo 
 Ncomponents = homogenization_maxNcomponents

!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
 call SNESCreate(PETSC_COMM_WORLD,solute_snes,ierr);    CHKERRQ(ierr)
 call SNESSetOptionsPrefix(solute_snes,'solute_',ierr); CHKERRQ(ierr) 
 allocate(localK(worldsize), source = 0); localK(worldrank+1) = grid3
 do proc = 1, worldsize
   call MPI_Bcast(localK(proc),1,MPI_INTEGER,proc-1,PETSC_COMM_WORLD,ierr)
 enddo  
 call DMDACreate3d(PETSC_COMM_WORLD, &
        DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, &                         ! cut off stencil at boundary
        DMDA_STENCIL_BOX, &                                                                         ! Moore (26) neighborhood around central point
        grid(1),grid(2),grid(3), &                                                                  ! global grid
        1, 1, worldsize, &
        2*Ncomponents, 1, &                                                                         ! #dof, ghost boundary width
        [grid(1)],[grid(2)],localK, &                                                                  ! local grid
        solute_grid,ierr)                                                                           ! handle, error
 CHKERRQ(ierr)
 call SNESSetDM(solute_snes,solute_grid,ierr); CHKERRQ(ierr)
 call DMsetFromOptions(solute_grid,ierr); CHKERRQ(ierr)
 call DMsetUp(solute_grid,ierr); CHKERRQ(ierr)
 call DMCreateGlobalVector(solute_grid,solution_current,ierr); CHKERRQ(ierr)                        ! current phase field vector 
 call DMCreateGlobalVector(solute_grid,solution_lastInc,ierr); CHKERRQ(ierr)                        ! last increment phase field vector 
 call DMSNESSetFunctionLocal(solute_grid,spectral_solute_formResidual,PETSC_NULL_SNES,ierr)         ! residual vector of same shape as solution vector
 CHKERRQ(ierr) 
 call DMSNESSetJacobianLocal(solute_grid,spectral_solute_formJacobian,PETSC_NULL_SNES,ierr)         ! function to evaluate stiffness matrix
 CHKERRQ(ierr)
 call SNESSetMaxLinearSolveFailures(solute_snes, huge(1), ierr); CHKERRQ(ierr)                      ! ignore linear solve failures 
 call SNESSetFromOptions(solute_snes,ierr); CHKERRQ(ierr)                                           ! pull it all together with additional cli arguments

!--------------------------------------------------------------------------------------------------
! init fields                 
 call VecSet(solution_current,0.0,ierr); CHKERRQ(ierr)
 call VecSet(solution_lastInc,0.0,ierr); CHKERRQ(ierr)
 call DMDAGetCorners(solute_grid,xstart,ystart,zstart,xend,yend,zend,ierr)                      ! local grid extent (corner coordinates and size)
 CHKERRQ(ierr) 
 xend = xstart+xend-1                                                                               ! upper bound = lower + size
 yend = ystart+yend-1
 zend = zstart+zend-1
 allocate(conc_current(1:Ncomponents,xstart:xend,ystart:yend,zstart:zend))
 allocate(conc_stagInc(1:Ncomponents,xstart:xend,ystart:yend,zstart:zend))
 allocate(conc_lastInc(1:Ncomponents,xstart:xend,ystart:yend,zstart:zend))

 cell = 0_pInt
 call DMDAVecGetArrayF90(solute_grid,solution_current,x_scal,ierr)
 CHKERRQ(ierr)                  
 do k = zstart, zend;  do j = ystart, yend;  do i = xstart, xend
   cell = cell + 1_pInt
   homog = material_homog(1,cell)
   x_scal(          0:  Ncomponents-1,i,j,k) = solute_flux_getInitialComponentPotential(1,cell)
   x_scal(Ncomponents:2*Ncomponents-1,i,j,k) = solute_flux_getComponentConc(1,cell)
   conc_current(    1:  Ncomponents  ,i,j,k) = solute_flux_getComponentConc(1,cell)
 enddo; enddo; enddo
 call DMDAVecRestoreArrayF90(solute_grid,solution_current,x_scal,ierr)
 CHKERRQ(ierr)
  
!--------------------------------------------------------------------------------------------------
! init FEM                 
 delta = geomSize/real(grid,pReal)                                                                  ! grid spacing
 BMat = 0.0_pReal
 BMat(1,1) = -1.0_pReal/delta(1)
 BMat(1,2) =  1.0_pReal/delta(1)
 BMat(2,1) = -1.0_pReal/delta(2)
 BMat(2,3) =  1.0_pReal/delta(2)
 BMat(3,1) = -1.0_pReal/delta(3)
 BMat(3,4) =  1.0_pReal/delta(3)
   
end subroutine spectral_solute_init
  
!--------------------------------------------------------------------------------------------------
!> @brief solution for the FEM scheme with internal iterations
!--------------------------------------------------------------------------------------------------
type(tSolutionState) function spectral_solute_solution(timeinc,timeinc_old)
 use numerics, only: &
   err_conc_tolAbs
 use spectral_utilities, only: &
   wgt
 use solute_flux, only: &
   solute_flux_calAndPutComponentConcRate

 implicit none

!--------------------------------------------------------------------------------------------------
! input data for solution
 real(pReal), intent(in) :: &
   timeinc, &                                                                                       !< increment in time for current solution
   timeinc_old                                                                                      !< increment in time of last increment
 DM :: da_local
 PetscScalar, pointer :: x_scal_current(:,:,:,:), x_scal_lastInc(:,:,:,:)
 real(pReal), dimension(Ncomponents) :: stagNorm, minConc, maxConc, avgConc
 integer(pInt) :: i, j, k, cell, comp

!--------------------------------------------------------------------------------------------------
! PETSc Data
 PetscErrorCode :: ierr   
 SNESConvergedReason :: reason

!--------------------------------------------------------------------------------------------------
! set module wide available data 
 params%timeinc = timeinc
 params%timeincOld = timeinc_old

!--------------------------------------------------------------------------------------------------
! solve BVP 
 call SNESSolve(solute_snes,PETSC_NULL_VEC,solution_current,ierr)
 CHKERRQ(ierr)
 call SNESGetConvergedReason(solute_snes,reason,ierr); CHKERRQ(ierr)
 spectral_solute_solution%converged = reason > 0

!--------------------------------------------------------------------------------------------------
! get local fields 
 call SNESGetDM(solute_snes,da_local,ierr); CHKERRQ(ierr)
 call DMDAVecGetArrayF90(da_local,solution_current,x_scal_current,ierr)
 CHKERRQ(ierr)
 call DMDAVecGetArrayF90(da_local,solution_lastInc,x_scal_lastInc,ierr)
 CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! updating phase state 
 avgConc  =  0.0_pReal
 minConc  =  huge(1.0_pReal)
 maxConc  = -huge(1.0_pReal)
 stagNorm = -huge(1.0_pReal)
 cell = 0
 do k = zstart, zend; do j = ystart, yend; do i = xstart, xend
   cell = cell + 1
   do comp = 1, Ncomponents
     minConc (comp) = min(minConc(comp),conc_current(comp,i,j,k))
     maxConc (comp) = max(maxConc(comp),conc_current(comp,i,j,k))
     avgConc (comp) = avgConc(comp) + conc_current(comp,i,j,k)
     stagNorm(comp) = max(stagNorm(comp),abs(conc_current(comp,i,j,k) - conc_stagInc(comp,i,j,k)))
   enddo    
   call solute_flux_calAndPutComponentConcRate(x_scal_current(          0:  Ncomponents-1,i,j,k), &
                                               x_scal_current(Ncomponents:2*Ncomponents-1,i,j,k), &
                                               params%timeinc,1,cell)
 enddo; enddo; enddo
 call MPI_Allreduce(MPI_IN_PLACE,minConc ,Ncomponents  ,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD,ierr)
 call MPI_Allreduce(MPI_IN_PLACE,maxConc ,Ncomponents  ,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD,ierr)
 call MPI_Allreduce(MPI_IN_PLACE,avgConc ,Ncomponents  ,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
 call MPI_Allreduce(MPI_IN_PLACE,stagNorm,Ncomponents  ,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD,ierr)

!--------------------------------------------------------------------------------------------------
! restore local fields 
 call DMDAVecRestoreArrayF90(da_local,solution_lastInc,x_scal_lastInc,ierr)
 CHKERRQ(ierr)
 call DMDAVecRestoreArrayF90(da_local,solution_current,x_scal_current,ierr)
 CHKERRQ(ierr)
 
!--------------------------------------------------------------------------------------------------
! reporting 
 conc_stagInc = conc_current
 spectral_solute_solution%stagConverged = all(stagNorm < err_conc_tolAbs)

 if (spectral_solute_solution%converged) then 
   write(6,'(/,a)') ' ... solute converged .........................................'
   do comp = 1, Ncomponents
     write(6,'(/,a,e12.5,2x,e12.5,2x,e11.4,/)',advance='no') ' Minimum|Maximum|Avg Conc          = ',&
       minConc(comp), maxConc(comp), avgConc(comp)*wgt
   enddo
   write(6,'(/,a)') ' ==========================================================================='
   flush(6) 
 endif 
end function spectral_solute_solution


!--------------------------------------------------------------------------------------------------
!> @brief forms the residual vector
!--------------------------------------------------------------------------------------------------
subroutine spectral_solute_formResidual(da_local,solution_current_local,residual_current_local,dummy,ierr)
 use numerics, only: &
   charLength
 use solute_flux, only: &
   solute_flux_calComponentConcandTangent, &
   solute_flux_getComponentMobility, &
   solute_flux_getEffChargeNumber
 use electrical_conduction, only: &
   electrical_conduction_getavgElectricalField_from_currentDensity

 implicit none

 DM                   :: da_local
 Vec                  :: solution_current_local, &
                         residual_current_local
 PetscScalar, pointer :: solution_current_scal(:,:,:,:), &
                         residual_current_scal(:,:,:,:)
 PetscScalar          :: solution_current_elem(4,2*Ncomponents), &
                         residual_current_elem(4,2*Ncomponents)
 real(pReal)          :: Conc         (Ncomponents), &
                         dConcdChemPot(Ncomponents,Ncomponents), &
                         dConcdGradC  (Ncomponents,Ncomponents), &
                         gradChempot      (3,          Ncomponents), &
                         elecMigrateForce (3,          Ncomponents), &
                         totdrivForce     (3,          Ncomponents), &
                         totFlux          (3,          Ncomponents), &
                         avgElfield       (3                      )
 PetscInt             :: i, j, k, cell
 PetscObject          :: dummy
 PetscErrorCode       :: ierr
  integer(pInt)       :: grad_dim


!--------------------------------------------------------------------------------------------------
! constructing residual
 call VecSet(residual_current_local,0.0,ierr);CHKERRQ(ierr)                                         ! set residual to zero
 call DMDAVecGetArrayF90(da_local,residual_current_local,residual_current_scal,ierr)
 CHKERRQ(ierr)                                                                                      ! residual_current_local(PETSc) -> residual_current_scal(our) with layout from da_local
 call DMDAVecGetArrayF90(da_local,solution_current_local,solution_current_scal,ierr)                ! sol_cur_local(PETSc, phase frac) -> sol_cut_scal(our) ...
 CHKERRQ(ierr)
 cell = 0
 do k = zstart, zend; do j = ystart, yend; do i = xstart, xend
   cell = cell + 1                                                                                  ! chemical potential at three neighboring points for differetiation
   solution_current_elem(1,            1:  Ncomponents) = solution_current_scal(          0:  Ncomponents-1,i  ,j  ,k  )      
   solution_current_elem(1,Ncomponents+1:2*Ncomponents) = solution_current_scal(Ncomponents:2*Ncomponents-1,i  ,j  ,k  )
   solution_current_elem(2,            1:  Ncomponents) = solution_current_scal(          0:  Ncomponents-1,i+1,j  ,k  )      
   solution_current_elem(2,Ncomponents+1:2*Ncomponents) = solution_current_scal(Ncomponents:2*Ncomponents-1,i+1,j  ,k  )
   solution_current_elem(3,            1:  Ncomponents) = solution_current_scal(          0:  Ncomponents-1,i  ,j+1,k  )      
   solution_current_elem(3,Ncomponents+1:2*Ncomponents) = solution_current_scal(Ncomponents:2*Ncomponents-1,i  ,j+1,k  )
   solution_current_elem(4,            1:  Ncomponents) = solution_current_scal(          0:  Ncomponents-1,i  ,j  ,k+1)      
   solution_current_elem(4,Ncomponents+1:2*Ncomponents) = solution_current_scal(Ncomponents:2*Ncomponents-1,i  ,j  ,k+1)
   
   call solute_flux_calComponentConcandTangent(Conc,dConcdChemPot,dConcdGradC, & 
                                               solution_current_scal(          0:  Ncomponents-1,i,j,k), &
                                               solution_current_scal(Ncomponents:2*Ncomponents-1,i,j,k), &
                                               1,cell)
   conc_current(1:Ncomponents,i,j,k) = Conc
   
   call electrical_conduction_getavgElectricalField_from_currentDensity(avgElfield, 1, cell)
   
   gradChempot      = matmul(BMat,solution_current_elem(1:4,1:Ncomponents))
   
   do grad_dim = 1,3
    elecMigrateForce(grad_dim,1:Ncomponents) = solute_flux_getEffChargeNumber(1,cell) * electronic_charge * avgElfield(grad_dim) 
   enddo
   
   totdrivForce =  gradChempot + elecMigrateForce
   totFlux = transpose(spread(solute_flux_getComponentMobility(1,cell),dim=2,ncopies=3))* totdrivForce
   
   residual_current_elem = 0.0_pReal
   residual_current_elem(1  ,            1:  Ncomponents) = &
   residual_current_elem(1  ,            1:  Ncomponents) + &
     conc_current(1:Ncomponents,i,j,k) - &
     conc_lastInc(1:Ncomponents,i,j,k)
   
   residual_current_elem(1:4,            1:  Ncomponents) = &
   residual_current_elem(1:4,            1:  Ncomponents) + &
     params%timeinc*transpose(spread(solute_flux_getComponentMobility(1,cell),dim=2,ncopies=4))* &
     matmul(transpose(BMat),matmul(BMat,solution_current_elem(1:4,            1:  Ncomponents)))
   
   residual_current_elem(1  ,Ncomponents+1:2*Ncomponents) = &
   residual_current_elem(1  ,Ncomponents+1:2*Ncomponents) + &
     solution_current_elem(1  ,Ncomponents+1:2*Ncomponents) - conc_current(1:Ncomponents,i,j,k)
   
   residual_current_elem(1:4,Ncomponents+1:2*Ncomponents) = &
   residual_current_elem(1:4,Ncomponents+1:2*Ncomponents) + &
     charLength*charLength* &
     matmul(transpose(BMat),matmul(BMat,solution_current_elem(1:4,Ncomponents+1:2*Ncomponents)))
     
   residual_current_scal(0:2*Ncomponents-1,i  ,j  ,k  ) = &
   residual_current_scal(0:2*Ncomponents-1,i  ,j  ,k  ) + residual_current_elem(1,1:2*Ncomponents)
   residual_current_scal(0:2*Ncomponents-1,i+1,j  ,k  ) = &
   residual_current_scal(0:2*Ncomponents-1,i+1,j  ,k  ) + residual_current_elem(2,1:2*Ncomponents)
   residual_current_scal(0:2*Ncomponents-1,i  ,j+1,k  ) = &
   residual_current_scal(0:2*Ncomponents-1,i  ,j+1,k  ) + residual_current_elem(3,1:2*Ncomponents)
   residual_current_scal(0:2*Ncomponents-1,i  ,j  ,k+1) = &
   residual_current_scal(0:2*Ncomponents-1,i  ,j  ,k+1) + residual_current_elem(4,1:2*Ncomponents)
 enddo; enddo; enddo
 call DMDAVecRestoreArrayF90(da_local,solution_current_local,solution_current_scal,ierr)            ! put back to PETSc data structure
 CHKERRQ(ierr)
 call DMDAVecRestoreArrayF90(da_local,residual_current_local,residual_current_scal,ierr)
 CHKERRQ(ierr)
   
end subroutine spectral_solute_formResidual


!--------------------------------------------------------------------------------------------------
!> @brief forms the FEM stiffness matrix
!--------------------------------------------------------------------------------------------------
subroutine spectral_solute_formJacobian(da_local,solution_current_local,Jac_pre,Jac,dummy,ierr)
 use math, only: &
   math_identity2nd
 use numerics, only: &
   charLength
 use solute_flux, only: &
   solute_flux_calComponentConcandTangent, &
   solute_flux_getComponentMobility

 implicit none

 DM                   :: da_local
 Vec                  :: solution_current_local
 Mat                  :: Jac, Jac_pre
 MatStencil           :: col(4,8*Ncomponents)
 PetscScalar, pointer :: solution_current_scal(:,:,:,:)
 PetscScalar          :: k_ele(8*Ncomponents,8*Ncomponents)
 real(pReal)          :: Conc         (Ncomponents), &
                         mobility     (Ncomponents), &
                         dConcdChemPot(Ncomponents,Ncomponents), &
                         dConcdGradC  (Ncomponents,Ncomponents)
 PetscInt             :: i, j, k, cell, comp
 PetscObject          :: dummy
 PetscErrorCode       :: ierr

!--------------------------------------------------------------------------------------------------
! assemble matrix
 call MatSetOption(Jac,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr); CHKERRQ(ierr)
 call MatSetOption(Jac,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr); CHKERRQ(ierr)
 call MatZeroEntries(Jac,ierr); CHKERRQ(ierr)
 call DMDAVecGetArrayF90(da_local,solution_current_local,solution_current_scal,ierr)
 CHKERRQ(ierr)
 cell = 0
 do k = zstart, zend; do j = ystart, yend; do i = xstart, xend
   cell = cell + 1
   do comp = 0, 2*Ncomponents-1 
     col(MatStencil_i,4*comp+1) = i
     col(MatStencil_j,4*comp+1) = j
     col(MatStencil_k,4*comp+1) = k
     col(MatStencil_c,4*comp+1) = comp
     col(MatStencil_i,4*comp+2) = i+1
     col(MatStencil_j,4*comp+2) = j
     col(MatStencil_k,4*comp+2) = k
     col(MatStencil_c,4*comp+2) = comp
     col(MatStencil_i,4*comp+3) = i
     col(MatStencil_j,4*comp+3) = j+1
     col(MatStencil_k,4*comp+3) = k
     col(MatStencil_c,4*comp+3) = comp
     col(MatStencil_i,4*comp+4) = i
     col(MatStencil_j,4*comp+4) = j
     col(MatStencil_k,4*comp+4) = k+1
     col(MatStencil_c,4*comp+4) = comp
   enddo

   call solute_flux_calComponentConcandTangent(Conc,dConcdChemPot,dConcdGradC, & 
                                               solution_current_scal(          0:  Ncomponents-1,i,j,k), &
                                               solution_current_scal(Ncomponents:2*Ncomponents-1,i,j,k), &
                                               1,cell)
   
   k_ele = 0.0
   k_ele(              1:4*Ncomponents:4,              1:4*Ncomponents:4) = &
     k_ele(              1:4*Ncomponents:4,              1:4*Ncomponents:4) + &
     dConcdChemPot
   k_ele(              1:4*Ncomponents:4,4*Ncomponents+1:8*Ncomponents:4) = &
     k_ele(              1:4*Ncomponents:4,4*Ncomponents+1:8*Ncomponents:4) + &
     dConcdGradC
   k_ele(4*Ncomponents+1:8*Ncomponents:4,              1:4*Ncomponents:4) = &
     k_ele(4*Ncomponents+1:8*Ncomponents:4,              1:4*Ncomponents:4) - &
     dConcdChemPot
   k_ele(4*Ncomponents+1:8*Ncomponents:4,4*Ncomponents+1:8*Ncomponents:4) = &
     k_ele(4*Ncomponents+1:8*Ncomponents:4,4*Ncomponents+1:8*Ncomponents:4) + &
     math_identity2nd(Ncomponents) - dConcdGradC
   mobility = solute_flux_getComponentMobility(1,cell)
   do comp = 0, Ncomponents-1
     k_ele(4*(            comp)+1:4*(            comp+1),4*(            comp)+1:4*(            comp+1)) = &
       k_ele(4*(            comp)+1:4*(            comp+1),4*(            comp)+1:4*(            comp+1)) + &
       params%timeinc*mobility(comp+1)*matmul(transpose(BMat),BMat)
     k_ele(4*(Ncomponents+comp)+1:4*(Ncomponents+comp+1),4*(Ncomponents+comp)+1:4*(Ncomponents+comp+1)) = &
       k_ele(4*(Ncomponents+comp)+1:4*(Ncomponents+comp+1),4*(Ncomponents+comp)+1:4*(Ncomponents+comp+1)) + &
       charLength*charLength*matmul(transpose(BMat),BMat)
   enddo  
   k_ele = transpose(k_ele)
   call MatSetValuesStencil(Jac,8*Ncomponents,col,8*Ncomponents,col,k_ele,ADD_VALUES,ierr)
   CHKERRQ(ierr)
 enddo; enddo; enddo
 call MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 call MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 call MatAssemblyBegin(Jac_pre,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 call MatAssemblyEnd(Jac_pre,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 call DMDAVecRestoreArrayF90(da_local,solution_current_local,solution_current_scal,ierr)
 CHKERRQ(ierr)
   
end subroutine spectral_solute_formJacobian


!--------------------------------------------------------------------------------------------------
!> @brief forwarding routine
!--------------------------------------------------------------------------------------------------
subroutine spectral_solute_forward()
 use spectral_utilities, only: &
   cutBack

 implicit none
 PetscErrorCode :: ierr

 if (cutBack) then 
   call VecCopy(solution_lastInc,solution_current,ierr); CHKERRQ(ierr)
 else
   conc_lastInc = conc_current
   conc_stagInc = conc_current
   call VecCopy(solution_current,solution_lastInc,ierr); CHKERRQ(ierr)                              ! also store current rate in rate_lastInc to recover upon cutback!!
 endif

end subroutine spectral_solute_forward

end module spectral_solute
