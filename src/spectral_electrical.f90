!--------------------------------------------------------------------------------------------------
!> @author Arka Lahiri, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Multiphase PETSc solver
!--------------------------------------------------------------------------------------------------
module spectral_electrical
#include <petsc/finclude/petscsnes.h>
#include <petsc/finclude/petscdmda.h>
use PETScdmda
use PETScsnes
 use prec, only: & 
   pInt, &
   pReal
 use math, only: &
   math_I3
 use spectral_utilities, only: &
   tSolutionState, &
   tSolutionParams
 use numerics, only: &
   worldrank, &
   worldsize

 implicit none
 private


!--------------------------------------------------------------------------------------------------
! derived types
 type(tSolutionParams), private :: params

!--------------------------------------------------------------------------------------------------
! PETSc data
 SNES, private :: electrical_snes
 Vec,  private :: solution_current, solution_lastInc, solution_stagInc

!--------------------------------------------------------------------------------------------------
! common pointwise data
 real   (pReal)                    , private              :: delta(3), doubledelta(3), BMat(3,4)
 PetscInt                          , private              :: xstart,ystart,zstart,xend,yend,zend
 
 real(pReal), private, dimension(3) :: &
   avgCurrentDensity=0_pReal                                                                        !the average of current density over the entire domain   
 
 integer(pInt), private :: &
   totalIter = 0_pInt                                                                               !< total iteration in current increment
   
 real(pReal), private :: &  
   minPotential, maxPotential, totPotential
 
    
 public :: &
   spectral_electrical_init, &
   spectral_electrical_solution, &
   spectral_electrical_forward
 external :: &
   PETScErrorF

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields and fills them with data, potentially from restart info
!--------------------------------------------------------------------------------------------------
subroutine spectral_electrical_init
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
 use mesh, only: &
   geomSize, &
   grid, &
   grid3
 
 implicit none
 DM :: electrical_grid
 PetscErrorCode       :: ierr
 integer(pInt), dimension(:), allocatable :: localK  
 integer(pInt) :: proc
 integer(pInt) :: i, j, k, cell, phase


 write(6,'(/,a)') ' <<<+-  spectral_electrical init  -+>>>'
 write(6,'(a15,a)')   ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
 call SNESCreate(PETSC_COMM_WORLD,electrical_snes,ierr);    CHKERRQ(ierr)
 call SNESSetOptionsPrefix(electrical_snes,'electrical_',ierr); CHKERRQ(ierr) 
 allocate(localK(worldsize), source = 0); localK(worldrank+1) = grid3
 do proc = 1, worldsize
   call MPI_Bcast(localK(proc),1,MPI_INTEGER,proc-1,PETSC_COMM_WORLD,ierr)
 enddo  
 call DMDACreate3d(PETSC_COMM_WORLD, &
        DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, &                         ! cut off stencil at boundary
        DMDA_STENCIL_BOX, &                                                                         ! Moore (26) neighborhood around central point
        grid(1),grid(2),grid(3), &                                                                  ! global grid
        1, 1, worldsize, &
        1, 1, &                                                                                     ! #dof, ghost boundary width
        [grid (1)],[grid(2)],localK, &                                                                  ! local grid
        electrical_grid,ierr)                                                                       ! handle, error
 CHKERRQ(ierr)
 call SNESSetDM(electrical_snes,electrical_grid,ierr); CHKERRQ(ierr)
 call DMsetFromOptions(electrical_grid,ierr); CHKERRQ(ierr)
 call DMsetUp(electrical_grid,ierr); CHKERRQ(ierr)
 call DMCreateGlobalVector(electrical_grid,solution_current,ierr); CHKERRQ(ierr)                    ! current electric potential vector 
 call DMCreateGlobalVector(electrical_grid,solution_lastInc,ierr); CHKERRQ(ierr)                    ! last increment electric potential vector 
 call DMCreateGlobalVector(electrical_grid,solution_stagInc,ierr); CHKERRQ(ierr)                    ! staggered increment electric potential vector 
 call DMSNESSetFunctionLocal(electrical_grid,spectral_electrical_formResidual,PETSC_NULL_SNES,ierr)! residual vector of same shape as solution vector
 CHKERRQ(ierr) 
 call DMSNESSetJacobianLocal(electrical_grid,spectral_electrical_formJacobian,PETSC_NULL_SNES,ierr)! function to evaluate stiffness matrix
 CHKERRQ(ierr)
 call SNESSetConvergenceTest(electrical_snes,spectral_electrical_converged,PETSC_NULL_SNES,PETSC_NULL_FUNCTION,ierr)  ! specify custom convergence check function "_converged"
 CHKERRQ(ierr)
 call SNESSetMaxLinearSolveFailures(electrical_snes, huge(1), ierr); CHKERRQ(ierr)                  ! ignore linear solve failures 
 call SNESSetFromOptions(electrical_snes,ierr); CHKERRQ(ierr)                                       ! pull it all together with additional cli arguments

!--------------------------------------------------------------------------------------------------
! init fields                 
 call VecSet(solution_current,0.0,ierr); CHKERRQ(ierr)
 call VecSet(solution_lastInc,0.0,ierr); CHKERRQ(ierr)
 call VecSet(solution_stagInc,0.0,ierr); CHKERRQ(ierr)
 call DMDAGetCorners(electrical_grid,xstart,ystart,zstart,xend,yend,zend,ierr)                      ! local grid extent (corner coordinates and size)
 CHKERRQ(ierr) 
 xend = xstart+xend-1                                                                               ! upper bound = lower + size
 yend = ystart+yend-1
 zend = zstart+zend-1

!--------------------------------------------------------------------------------------------------
! init FEM                 
 delta = geomSize/real(grid,pReal)                                                                  ! grid spacing
 doubledelta = 2.0_pReal*delta                                                                        ! double grid spacing            
 BMat = 0.0_pReal
 BMat(1,1) = -1.0_pReal/delta(1)
 BMat(1,2) =  1.0_pReal/delta(1)
 BMat(2,1) = -1.0_pReal/delta(2)
 BMat(2,3) =  1.0_pReal/delta(2)
 BMat(3,1) = -1.0_pReal/delta(3)
 BMat(3,4) =  1.0_pReal/delta(3)
  
end subroutine spectral_electrical_init
  
!--------------------------------------------------------------------------------------------------
!> @brief solution for the FEM scheme with internal iterations
!--------------------------------------------------------------------------------------------------
type(tSolutionState) function spectral_electrical_solution(avgElectricalField)
 use IO, only: &
   IO_error
 use numerics, only: &
   err_elecPot_tolRel, &
   err_elecPot_tolAbs
 use spectral_utilities, only: &
   tVectorBoundaryCondition
 use electrical_conduction, only: &
   electrical_conduction_calAndPutCurrentDensity, &
   electrical_conduction_PutElectricPotential
 use FEsolving, only: &
   terminallyIll
   

 implicit none

!--------------------------------------------------------------------------------------------------
! input data for solution
 type(tVectorBoundaryCondition),      intent(in) :: &
   avgElectricalField
 DM :: da_local
 Vec :: solution_current_local
 PetscScalar, pointer :: x_scal_current(:,:,:), x_scal_lastInc(:,:,:), x_scal_stagInc(:,:,:)
 real(pReal) :: Potential_current,Potential_lastInc, Potential_stagInc
 real(pReal) :: stagNorm, solnNorm
 integer(pInt) :: i, j, k, cell
 real(pReal) :: ElectricalField(3)

!--------------------------------------------------------------------------------------------------
! PETSc Data
 PetscErrorCode :: ierr   
 SNESConvergedReason :: reason

!--------------------------------------------------------------------------------------------------
! set module wide available data 
 params%electricalField_BC = avgElectricalField%values

!--------------------------------------------------------------------------------------------------
! solve BVP 
 call SNESSolve(electrical_snes,PETSC_NULL_VEC,solution_current,ierr)
 CHKERRQ(ierr)
 call SNESGetConvergedReason(electrical_snes,reason,ierr); CHKERRQ(ierr)
 spectral_electrical_solution%converged = reason > 0
 spectral_electrical_solution%iterationsNeeded = totalIter
 spectral_electrical_solution%termIll = terminallyIll
 terminallyIll = .false.
 if (reason == -4) call IO_error(893_pInt)  

!--------------------------------------------------------------------------------------------------
! get local fields 
 call SNESGetDM(electrical_snes,da_local,ierr); CHKERRQ(ierr)
 call DMGetLocalVector(da_local,solution_current_local,ierr); CHKERRQ(ierr)
 call DMGlobalToLocalBegin(da_local,solution_current,INSERT_VALUES,solution_current_local,ierr)                       !< retrieve my partition of global solution vector
 CHKERRQ(ierr)
 call DMGlobalToLocalEnd(da_local,solution_current,INSERT_VALUES,solution_current_local,ierr)
 CHKERRQ(ierr)
 call DMDAVecGetArrayF90(da_local,solution_current_local,x_scal_current,ierr)
 CHKERRQ(ierr)
 call DMDAVecGetArrayF90(da_local,solution_lastInc,x_scal_lastInc,ierr)
 CHKERRQ(ierr)
 call DMDAVecGetArrayF90(da_local,solution_stagInc,x_scal_stagInc,ierr)
 CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! updating current density state 
 totpotential   =  0.0_pReal
 minpotential   =  huge(1.0_pReal)
 maxpotential   = -huge(1.0_pReal)
 stagNorm       = -huge(1.0_pReal)
 solnNorm       = -huge(1.0_pReal)
  cell = 0
  do k = zstart, zend; do j = ystart, yend; do i = xstart, xend
   cell = cell + 1
   potential_stagInc = x_scal_stagInc(i,j,k)
   potential_current = x_scal_current(i,j,k)
   potential_lastInc = x_scal_lastInc(i,j,k)
   
   minPotential = min(minPotential,Potential_current)
   maxPotential = max(maxPotential,Potential_current)
   totPotential = totPotential + Potential_current
   stagNorm = max(stagNorm,abs(Potential_current - Potential_stagInc))
   solnNorm = max(solnNorm,abs(Potential_current                    ))
     
   ElectricalField(1)= ((x_scal_current(i+1,  j   ,k  ) - &
                         x_scal_current(i-1,  j   ,k  ))/(doubledelta(1))) + params%electricalField_BC(1)
   ElectricalField(2)= ((x_scal_current(i  , j+1  ,k  ) - &
                         x_scal_current(i  , j-1  ,k  ))/(doubledelta(2))) + params%electricalField_BC(2)      
   ElectricalField(3)= ((x_scal_current(i ,  j   ,k+1 ) - &
                         x_scal_current(i ,  j   ,k-1 ))/(doubledelta(3))) + params%electricalField_BC(3)                                               
   
   call electrical_conduction_calAndPutCurrentDensity(ElectricalField,1,cell)
   call electrical_conduction_PutElectricPotential(potential_current,1,cell)

enddo; enddo; enddo

 call MPI_Allreduce(MPI_IN_PLACE,minPotential    ,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD,ierr)
 call MPI_Allreduce(MPI_IN_PLACE,maxPotential    ,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD,ierr)
 call MPI_Allreduce(MPI_IN_PLACE,totPotential    ,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
 call MPI_Allreduce(MPI_IN_PLACE,stagNorm        ,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD,ierr)
 call MPI_Allreduce(MPI_IN_PLACE,solnNorm        ,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD,ierr)

!--------------------------------------------------------------------------------------------------
! restore local fields 
 call DMDAVecRestoreArrayF90(da_local,solution_stagInc,x_scal_stagInc,ierr)
 CHKERRQ(ierr)
 call DMDAVecRestoreArrayF90(da_local,solution_lastInc,x_scal_lastInc,ierr)
 CHKERRQ(ierr)
 call DMDAVecRestoreArrayF90(da_local,solution_current_local,x_scal_current,ierr)
 CHKERRQ(ierr)
 call DMRestoreLocalVector(da_local,solution_current_local,ierr); CHKERRQ(ierr)
 call VecCopy(solution_current,solution_stagInc,ierr); CHKERRQ(ierr)
 
end function spectral_electrical_solution


!--------------------------------------------------------------------------------------------------
!> @brief forms the residual vector
!--------------------------------------------------------------------------------------------------
subroutine spectral_electrical_formResidual(da_local,solution_current_local,residual_current_local,dummy,ierr)
use electrical_conduction, only: &
   electrical_conduction_getFlux
use spectral_utilities, only: &
   wgt
 implicit none

 DM                   :: da_local
 Vec                  :: solution_current_local, &
                         residual_current_local
 PetscScalar, pointer :: solution_current_scal(:,:,:), &
                         residual_current_scal(:,:,:)
 PetscScalar          :: solution_current_elem(4), &
                         residual_current_elem(4)
 PetscScalar          :: Potential_current, &
                         Potential_grad(3)
 real(pReal)          :: Potential_flux(3)
 PetscInt             :: i, j, k, cell
 PetscObject          :: dummy
 PetscErrorCode       :: ierr
   
 
!--------------------------------------------------------------------------------------------------
! constructing residual
 call VecSet(residual_current_local,0.0,ierr);CHKERRQ(ierr)                                         ! set residual to zero
 call DMDAVecGetArrayF90(da_local,residual_current_local,residual_current_scal,ierr)
 CHKERRQ(ierr)                                                                                      ! residual_current_local(PETSc) -> residual_current_scal(our) with layout from da_local
 call DMDAVecGetArrayF90(da_local,solution_current_local,solution_current_scal,ierr)                ! sol_cur_local(PETSc, phase frac) -> sol_cut_scal(our) ...
 CHKERRQ(ierr)
 
 cell = 0
 avgCurrentDensity = 0_pReal
 do k = zstart, zend; do j = ystart, yend; do i = xstart, xend
   cell = cell + 1
   solution_current_elem(1) = solution_current_scal(i  ,j  ,k  )                                    ! phase fracs at three neighboring points for differetiation
   solution_current_elem(2) = solution_current_scal(i+1,j  ,k  )
   solution_current_elem(3) = solution_current_scal(i  ,j+1,k  )
   solution_current_elem(4) = solution_current_scal(i  ,j  ,k+1)
   
   Potential_current = solution_current_scal(i,j,k)                 
   Potential_grad = matmul(BMat,solution_current_elem) + params%electricalField_BC
   
   call electrical_conduction_getFlux(Potential_flux,Potential_grad,1,cell) 
   
   avgCurrentDensity(1:3) = avgCurrentDensity(1:3) + Potential_flux(1:3)
      
   residual_current_elem = matmul(transpose(BMat), Potential_flux)
   residual_current_scal(i  ,j  ,k  ) = &
   residual_current_scal(i  ,j  ,k  ) + residual_current_elem(1)
   residual_current_scal(i+1,j  ,k  ) = &
   residual_current_scal(i+1,j  ,k  ) + residual_current_elem(2)
   residual_current_scal(i  ,j+1,k  ) = &
   residual_current_scal(i  ,j+1,k  ) + residual_current_elem(3)
   residual_current_scal(i  ,j  ,k+1) = &
   residual_current_scal(i  ,j  ,k+1) + residual_current_elem(4)
   
enddo; enddo; enddo
 
 call MPI_Allreduce(MPI_IN_PLACE,avgCurrentDensity    ,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
   
 avgCurrentDensity=avgCurrentDensity*wgt
  
 call DMDAVecRestoreArrayF90(da_local,solution_current_local,solution_current_scal,ierr)            ! put back to PETSc data structure
 CHKERRQ(ierr)
 call DMDAVecRestoreArrayF90(da_local,residual_current_local,residual_current_scal,ierr)
 CHKERRQ(ierr)
 
end subroutine spectral_electrical_formResidual


!--------------------------------------------------------------------------------------------------
!> @brief forms the FEM stiffness matrix
!--------------------------------------------------------------------------------------------------
subroutine spectral_electrical_formJacobian(da_local,solution_current_local,Jac_pre,Jac,dummy,ierr)
 use math, only: &
   math_identity2nd
use electrical_conduction, only: &
   electrical_conduction_getFluxTangent

 implicit none

 DM                   :: da_local
 Vec                  :: solution_current_local
 Mat                  :: Jac, Jac_pre
 MatStencil           :: col(4,4)
 PetscScalar          :: k_ele(4,4)
 real(pReal)          :: Potential_fluxTangent(3,3)
 PetscInt             :: i, j, k, cell
 PetscObject          :: dummy
 PetscErrorCode       :: ierr

!--------------------------------------------------------------------------------------------------
! assemble matrix
 call MatSetOption(Jac,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr); CHKERRQ(ierr)
 call MatSetOption(Jac,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr); CHKERRQ(ierr)
 call MatZeroEntries(Jac,ierr); CHKERRQ(ierr)
 cell = 0
 do k = zstart, zend; do j = ystart, yend; do i = xstart, xend
   cell = cell + 1
   col(MatStencil_i,1) = i
   col(MatStencil_j,1) = j
   col(MatStencil_k,1) = k
   col(MatStencil_c,1) = 0
   col(MatStencil_i,2) = i+1
   col(MatStencil_j,2) = j
   col(MatStencil_k,2) = k
   col(MatStencil_c,2) = 0
   col(MatStencil_i,3) = i
   col(MatStencil_j,3) = j+1
   col(MatStencil_k,3) = k
   col(MatStencil_c,3) = 0
   col(MatStencil_i,4) = i
   col(MatStencil_j,4) = j
   col(MatStencil_k,4) = k+1
   col(MatStencil_c,4) = 0

   call electrical_conduction_getFluxTangent(Potential_fluxTangent,1,cell)

   k_ele = matmul(transpose(BMat),matmul(Potential_fluxTangent,BMat))
   k_ele = transpose(k_ele)
   call MatSetValuesStencil(Jac,4,col,4,col,k_ele,ADD_VALUES,ierr)
   CHKERRQ(ierr)
 enddo; enddo; enddo
 call MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 call MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 call MatAssemblyBegin(Jac_pre,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 call MatAssemblyEnd(Jac_pre,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 
end subroutine spectral_electrical_formJacobian


!--------------------------------------------------------------------------------------------------
!> @brief convergence check
!--------------------------------------------------------------------------------------------------
subroutine spectral_electrical_converged(snes_local,PETScIter,xnorm,snorm,fnorm,reason,dummy,ierr)
 use numerics, only: &
   itmax, &
   itmin, &
   err_elecPot_tolRel, &
   err_elecPot_tolAbs 
 use FEsolving, only: &
   terminallyIll

 implicit none
 SNES :: snes_local
 PetscInt :: PETScIter
 PetscReal :: &
   xnorm, &
   snorm, &
   fnorm
 SNESConvergedReason :: reason
 PetscObject :: dummy
 PetscErrorCode :: ierr
 real(pReal) :: &
   divTol

 divTol = max(maxval(abs(avgCurrentDensity))*err_elecPot_tolRel   ,err_elecPot_tolAbs)

 converged: if (fnorm/divTol < 1.0_pReal) then         
   reason = 1
 else converged
   reason = 0
 endif converged

!--------------------------------------------------------------------------------------------------
! report
write(6,'(/,a)') ' ... electrical converged .....................................'
write(6,'(/,a,e12.5,2x,e12.5,2x,e11.4,/)',advance='no') ' Minimum|Maximum|Tot Potential           = ',&
   minPotential, maxPotential, totPotential
write(6,'(/,a)') ' ==========================================================================='
flush(6) 
!--------------------------------------------------------------------------------------------------- 

end subroutine spectral_electrical_converged



!--------------------------------------------------------------------------------------------------
!> @brief forwarding routine
!--------------------------------------------------------------------------------------------------
subroutine spectral_electrical_forward(avgElectricalField)
 use spectral_utilities, only: &
   cutBack
 use electrical_conduction, only: &
    electrical_conduction_calAndPutCurrentDensity, &
    electrical_conduction_PutElectricPotential
 use spectral_utilities, only: &
   tVectorBoundaryCondition    
    

 implicit none
 type(tVectorBoundaryCondition),      intent(in) :: &
   avgElectricalField
 DM :: da_local
 Vec :: solution_current_local
 PetscScalar, pointer :: solution_current_scal(:,:,:)
 integer(pInt) :: i, j, k, cell
 real(pReal) :: ElectricalField(3), &
                ElectricPotential
 PetscErrorCode :: ierr
 

 call SNESGetDM(electrical_snes,da_local,ierr); CHKERRQ(ierr)
 if (cutBack) then 
   call VecCopy(solution_lastInc,solution_current,ierr); CHKERRQ(ierr)
   call VecCopy(solution_lastInc,solution_stagInc,ierr); CHKERRQ(ierr)
   call DMGetLocalVector(da_local,solution_current_local,ierr); CHKERRQ(ierr)
   call DMGlobalToLocalBegin(da_local,solution_current,INSERT_VALUES,solution_current_local,ierr)                       !< retrieve my partition of global solution vector
   CHKERRQ(ierr)
   call DMGlobalToLocalEnd(da_local,solution_current,INSERT_VALUES,solution_current_local,ierr)
   CHKERRQ(ierr)
   call DMDAVecGetArrayF90(da_local,solution_current_local,solution_current_scal,ierr)
   CHKERRQ(ierr)
   cell = 0
   do k = zstart, zend; do j = ystart, yend; do i = xstart, xend                                    ! walk through my (sub)domain
     cell = cell + 1      
     ElectricalField(1)= ((solution_current_scal(i+1,  j   ,k  ) - &
                           solution_current_scal(i-1,  j   ,k  ))/(doubledelta(1))) + avgElectricalField%values(1)
     ElectricalField(2)= ((solution_current_scal(i , j+1   ,k  ) - &
                           solution_current_scal(i,  j-1   ,k  ))/(doubledelta(2))) + avgElectricalField%values(2)        
     ElectricalField(3)= ((solution_current_scal(i ,  j    ,k+1) - &
                           solution_current_scal(i,   j    ,k-1))/(doubledelta(3))) + avgElectricalField%values(3)                                                  
                           
     ElectricPotential = solution_current_scal(i,  j   ,k)                      
     
     call electrical_conduction_calAndPutCurrentDensity(ElectricalField,1,cell)
     call electrical_conduction_PutElectricPotential(ElectricPotential,1,cell)
   enddo; enddo; enddo
   call DMDAVecRestoreArrayF90(da_local,solution_current_local,solution_current_scal,ierr)
   CHKERRQ(ierr)
   call DMRestoreLocalVector(da_local,solution_current_local,ierr); CHKERRQ(ierr)
 else
   call VecCopy(solution_current,solution_lastInc,ierr); CHKERRQ(ierr)                              ! also store current rate in rate_lastInc to recover upon cutback!!
   call VecCopy(solution_current,solution_stagInc,ierr); CHKERRQ(ierr)                              ! also store current rate in rate_lastInc to recover upon cutback!!
 endif
 
end subroutine spectral_electrical_forward

end module spectral_electrical