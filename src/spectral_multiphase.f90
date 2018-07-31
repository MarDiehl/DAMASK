!--------------------------------------------------------------------------------------------------
!> @author Chuanlai Liu, Shanghai Jiao Tong University
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Multiphase PETSc solver
!--------------------------------------------------------------------------------------------------
module spectral_multiphase
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
 SNES, private :: multiphase_snes
 Vec,  private :: solution_current, solution_lastInc

!--------------------------------------------------------------------------------------------------
! common pointwise data
 integer(pInt)                     , private              :: NActivePhases
 real   (pReal)                    , private              :: delta(3), BMat(3,4)
 PetscInt                          , private              :: xstart,ystart,zstart,xend,yend,zend

 public :: &
   spectral_multiphase_init, &
   spectral_multiphase_solution, &
   spectral_multiphase_forward
 external :: &
   PETScErrorF

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields and fills them with data, potentially from restart info
!--------------------------------------------------------------------------------------------------
subroutine spectral_multiphase_init
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran >4.6 at the moment)
 use IO, only: &
   IO_error, &
   IO_intOut, &
   IO_read_realFile, &
   IO_timeStamp
 use mesh, only: &
   geomSize, &
   grid, &
   grid3
 use material, only: &
   homogenization_Ngrains, &
   homogenization_maxNgrains, &
   homogenization_active, &
   material_homog, &   
   phasefracMapping, &
   phasefrac
 use config, only: &  
   material_Nhomogenization
   
 implicit none
 DM :: multiphase_grid
 Vec :: uBound, lBound
 PetscScalar, pointer :: x_scal(:,:,:,:)
 PetscErrorCode       :: ierr
 integer(pInt), dimension(:), allocatable :: localK  
 integer(pInt) :: proc
 integer(pInt) :: i, j, k, cell, homog, phase
 character(len=100) :: snes_type


 write(6,'(/,a)') ' <<<+-  spectral_multiphase init  -+>>>'
 write(6,'(a15,a)')   ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
 do homog = 1_pInt, material_Nhomogenization
   if (homogenization_active(homog)) then
     if (homogenization_Ngrains(homog) /= homogenization_maxNgrains) &
       call IO_error(211_pInt,ext_msg='# of phases mismatch in spectral_multiphase')
   endif
 enddo 
 NActivePhases = homogenization_maxNgrains
 if (NActivePhases < 2_pInt) call IO_error(211_pInt,ext_msg='# of phases < 2 in spectral_multiphase')

!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
 call SNESCreate(PETSC_COMM_WORLD,multiphase_snes,ierr);    CHKERRQ(ierr)
 call SNESSetOptionsPrefix(multiphase_snes,'multiphase_',ierr); CHKERRQ(ierr) 
 allocate(localK(worldsize), source = 0); localK(worldrank+1) = grid3
 do proc = 1, worldsize
   call MPI_Bcast(localK(proc),1,MPI_INTEGER,proc-1,PETSC_COMM_WORLD,ierr)
 enddo  
 call DMDACreate3d(PETSC_COMM_WORLD, &
        DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, &                         ! cut off stencil at boundary
        DMDA_STENCIL_BOX, &                                                                         ! Moore (26) neighborhood around central point
        grid(1),grid(2),grid(3), &                                                                  ! global grid
        1, 1, worldsize, &
        NActivePhases, 1, &                                                                         ! #dof, ghost boundary width
        [grid(1)],[grid(2)],localK, &                                                                  ! local grid
        multiphase_grid,ierr)                                                                       ! handle, error
 CHKERRQ(ierr)
 call SNESSetDM(multiphase_snes,multiphase_grid,ierr); CHKERRQ(ierr)
 call DMsetFromOptions(multiphase_grid,ierr); CHKERRQ(ierr)
 call DMsetUp(multiphase_grid,ierr); CHKERRQ(ierr)
 call DMCreateGlobalVector(multiphase_grid,solution_current,ierr); CHKERRQ(ierr)                    ! current phase field vector 
 call DMCreateGlobalVector(multiphase_grid,solution_lastInc,ierr); CHKERRQ(ierr)                    ! last increment phase field vector 
 call DMSNESSetFunctionLocal(multiphase_grid,spectral_multiphase_formResidual,PETSC_NULL_SNES,ierr) ! residual vector of same shape as solution vector
 CHKERRQ(ierr) 
 call DMSNESSetJacobianLocal(multiphase_grid,spectral_multiphase_formJacobian,PETSC_NULL_SNES,ierr) ! function to evaluate stiffness matrix
 CHKERRQ(ierr)
 call SNESSetMaxLinearSolveFailures(multiphase_snes, huge(1), ierr); CHKERRQ(ierr)                  ! ignore linear solve failures 
 call SNESSetFromOptions(multiphase_snes,ierr); CHKERRQ(ierr)                                       ! pull it all together with additional cli arguments
 call SNESGetType(multiphase_snes,snes_type,ierr); CHKERRQ(ierr)
 if (trim(snes_type) == 'vinewtonrsls' .or. &
     trim(snes_type) == 'vinewtonssls') then
   call DMGetGlobalVector(multiphase_grid,lBound,ierr); CHKERRQ(ierr)
   call DMGetGlobalVector(multiphase_grid,uBound,ierr); CHKERRQ(ierr)
   call VecSet(lBound,0.0,ierr); CHKERRQ(ierr)
   call VecSet(uBound,1.0,ierr); CHKERRQ(ierr)
   call SNESVISetVariableBounds(multiphase_snes,lBound,uBound,ierr)                                 !< variable bounds for variational inequalities like contact mechanics, damage etc.
   call DMRestoreGlobalVector(multiphase_grid,lBound,ierr); CHKERRQ(ierr)
   call DMRestoreGlobalVector(multiphase_grid,uBound,ierr); CHKERRQ(ierr)
 endif

!--------------------------------------------------------------------------------------------------
! init fields                 
 call VecSet(solution_current,0.0,ierr); CHKERRQ(ierr)
 call VecSet(solution_lastInc,0.0,ierr); CHKERRQ(ierr)
 call DMDAGetCorners(multiphase_grid,xstart,ystart,zstart,xend,yend,zend,ierr)                      ! local grid extent (corner coordinates and size)
 CHKERRQ(ierr) 
 xend = xstart+xend-1                                                                               ! upper bound = lower + size
 yend = ystart+yend-1
 zend = zstart+zend-1

 cell = 0_pInt
 call DMDAVecGetArrayF90(multiphase_grid,solution_current,x_scal,ierr)
 CHKERRQ(ierr)                  
 do k = zstart, zend;  do j = ystart, yend;  do i = xstart, xend
   cell = cell + 1_pInt
   homog = material_homog(1,cell)
   do phase = 1, NActivePhases
     x_scal(phase-1,i,j,k) =  phasefrac(homog)%p(phase,phasefracMapping(homog)%p(1,cell))
   enddo
 enddo; enddo; enddo
 call DMDAVecRestoreArrayF90(multiphase_grid,solution_current,x_scal,ierr)
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
   
end subroutine spectral_multiphase_init
  
!--------------------------------------------------------------------------------------------------
!> @brief solution for the FEM scheme with internal iterations
!--------------------------------------------------------------------------------------------------
type(tSolutionState) function spectral_multiphase_solution(timeinc,timeinc_old)
 use numerics, only: &
   err_phasefr_tolAbs
 use material, only: &
   material_homog, &
   phasefracMapping, &
   phasefrac
 use spectral_utilities, only: &
   tBoundaryCondition, &
   wgt
 use homogenization_multiphase, only: &
   homogenization_multiphase_putPhaseFrac, &
   homogenization_multiphase_putInterfaceNormals

 implicit none

!--------------------------------------------------------------------------------------------------
! input data for solution
 real(pReal), intent(in) :: &
   timeinc, &                                                                                       !< increment in time for current solution
   timeinc_old                                                                                      !< increment in time of last increment
 DM :: da_local
 Vec :: solution_current_local
 PetscScalar, pointer :: x_scal_current(:,:,:,:), x_scal_lastInc(:,:,:,:)
 real(pReal), dimension(NActivePhases) :: phi_current, phi_lastInc, phi_stagInc
 real(pReal), dimension(NActivePhases) :: stagNorm, minPhi, maxPhi, avgPhi
 integer(pInt) :: i, j, k, cell, phase, phaseI, phaseJ, homog
 real(pReal) :: interfaceNormal(3)

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
 call SNESSolve(multiphase_snes,PETSC_NULL_VEC,solution_current,ierr)
 CHKERRQ(ierr)
 call SNESGetConvergedReason(multiphase_snes,reason,ierr); CHKERRQ(ierr)
 spectral_multiphase_solution%converged = reason > 0

!--------------------------------------------------------------------------------------------------
! get local fields 
 call SNESGetDM(multiphase_snes,da_local,ierr); CHKERRQ(ierr)
 call DMGetLocalVector(da_local,solution_current_local,ierr); CHKERRQ(ierr)
 call DMGlobalToLocalBegin(da_local,solution_current,INSERT_VALUES,solution_current_local,ierr)                       !< retrieve my partition of global solution vector
 CHKERRQ(ierr)
 call DMGlobalToLocalEnd(da_local,solution_current,INSERT_VALUES,solution_current_local,ierr)
 CHKERRQ(ierr)
 call DMDAVecGetArrayF90(da_local,solution_current_local,x_scal_current,ierr)
 CHKERRQ(ierr)
 call DMDAVecGetArrayF90(da_local,solution_lastInc,x_scal_lastInc,ierr)
 CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! updating phase state 
 avgPhi   =  0.0_pReal
 minPhi   =  huge(1.0_pReal)
 maxPhi   = -huge(1.0_pReal)
 stagNorm = -huge(1.0_pReal)
 cell = 0
 do k = zstart, zend; do j = ystart, yend; do i = xstart, xend
   cell = cell + 1
   homog = material_homog(1,cell)
   do phase = 1, NActivePhases
     phi_stagInc(phase) = phasefrac(homog)%p(phase,phasefracMapping(homog)%p(1,cell))
     phi_current(phase) = x_scal_current(phase-1,i,j,k)
     phi_lastInc(phase) = x_scal_lastInc(phase-1,i,j,k)
   enddo
   do phase = 1, NActivePhases
     minPhi(phase) = min(minPhi(phase),phi_current(phase))
     maxPhi(phase) = max(maxPhi(phase),phi_current(phase))
     avgPhi(phase) = avgPhi(phase) + phi_current(phase)
     stagNorm(phase) = max(stagNorm(phase),abs(phi_current(phase) - phi_stagInc(phase)))
   enddo
   do phaseI = 1, NActivePhases-1
     do phaseJ = phaseI+1, NActivePhases
       interfaceNormal(1) = (x_scal_current(phaseI-1,i+1,j  ,k  ) - &
                             x_scal_current(phaseI-1,i-1,j  ,k  ))/delta(1) - &
                            (x_scal_current(phaseJ-1,i+1,j  ,k  ) - &
                             x_scal_current(phaseJ-1,i-1,j  ,k  ))/delta(1)
       interfaceNormal(2) = (x_scal_current(phaseI-1,i  ,j+1,k  ) - &
                             x_scal_current(phaseI-1,i  ,j-1,k  ))/delta(2) - &
                            (x_scal_current(phaseJ-1,i  ,j+1,k  ) - &
                             x_scal_current(phaseJ-1,i  ,j-1,k  ))/delta(2)
       interfaceNormal(3) = (x_scal_current(phaseI-1,i  ,j  ,k+1) - &
                             x_scal_current(phaseI-1,i  ,j  ,k-1))/delta(3) - &
                            (x_scal_current(phaseJ-1,i  ,j  ,k+1) - &
                             x_scal_current(phaseJ-1,i  ,j  ,k-1))/delta(3)
       if (norm2(interfaceNormal) > 0.0_pReal) then
         interfaceNormal = interfaceNormal/norm2(interfaceNormal)
       else
         interfaceNormal = 0.0_pReal
       endif
       call homogenization_multiphase_putInterfaceNormals(interfaceNormal,phaseI,phaseJ,1,cell)
     enddo
   enddo    
   call homogenization_multiphase_putPhaseFrac(phi_current,1,cell)
 enddo; enddo; enddo
 call MPI_Allreduce(MPI_IN_PLACE,minPhi    ,NActivePhases,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD,ierr)
 call MPI_Allreduce(MPI_IN_PLACE,maxPhi    ,NActivePhases,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD,ierr)
 call MPI_Allreduce(MPI_IN_PLACE,avgPhi    ,NActivePhases,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
 call MPI_Allreduce(MPI_IN_PLACE,stagNorm,  NActivePhases,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD,ierr)

!--------------------------------------------------------------------------------------------------
! restore local fields 
 call DMDAVecRestoreArrayF90(da_local,solution_lastInc,x_scal_lastInc,ierr)
 CHKERRQ(ierr)
 call DMDAVecRestoreArrayF90(da_local,solution_current_local,x_scal_current,ierr)
 CHKERRQ(ierr)
 call DMRestoreLocalVector(da_local,solution_current_local,ierr); CHKERRQ(ierr)
 
!--------------------------------------------------------------------------------------------------
! reporting 
 spectral_multiphase_solution%stagConverged = all(stagNorm < err_phasefr_tolAbs)

 if (spectral_multiphase_solution%converged) then 
   write(6,'(/,a)') ' ... multiphase converged .....................................'
   do phase = 1, NActivePhases
     write(6,'(/,a,e12.5,2x,e12.5,2x,e11.4,/)',advance='no') ' Minimum|Maximum|Avg Phi           = ',&
       minPhi(phase), maxPhi(phase), avgPhi(phase)*wgt
   enddo
   write(6,'(/,a)') ' ==========================================================================='
   flush(6) 
 endif 
end function spectral_multiphase_solution


!--------------------------------------------------------------------------------------------------
!> @brief forms the residual vector
!--------------------------------------------------------------------------------------------------
subroutine spectral_multiphase_formResidual(da_local,solution_current_local,residual_current_local,dummy,ierr)
 use homogenization_multiphase, only: &
   homogenization_multiphase_getPhaseFlux, &
   homogenization_multiphase_getPhaseSource

 implicit none

 DM                   :: da_local
 Vec                  :: solution_current_local, &
                         residual_current_local
 PetscScalar, pointer :: solution_current_scal(:,:,:,:), &
                         solution_lastInc_scal(:,:,:,:), &
                         residual_current_scal(:,:,:,:)
 PetscScalar          :: solution_current_elem(4,NActivePhases), &
                         residual_current_elem(4,NActivePhases)
 PetscScalar          :: phi_current(NActivePhases), &
                         phi_lastInc(NActivePhases), &
                         phi_grad(3,NActivePhases)
 real(pReal)          :: phi_source(NActivePhases), &
                         phi_flux(3,NActivePhases)
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
 call DMDAVecGetArrayF90(da_local,solution_lastInc,      solution_lastInc_scal,ierr)                ! solution_lastInc(PETSc) -> sol_last_scal(our, but phase fracs)
 CHKERRQ(ierr)
 cell = 0
 do k = zstart, zend; do j = ystart, yend; do i = xstart, xend
   cell = cell + 1
   solution_current_elem(1,1:NActivePhases) = solution_current_scal(0:NActivePhases-1,i  ,j  ,k  )  ! phase fracs at three neighboring points for differetiation
   solution_current_elem(2,1:NActivePhases) = solution_current_scal(0:NActivePhases-1,i+1,j  ,k  )
   solution_current_elem(3,1:NActivePhases) = solution_current_scal(0:NActivePhases-1,i  ,j+1,k  )
   solution_current_elem(4,1:NActivePhases) = solution_current_scal(0:NActivePhases-1,i  ,j  ,k+1)
   
   phi_current(1:NActivePhases) = solution_current_scal(0:NActivePhases-1,i,j,k)                 
   phi_lastInc(1:NActivePhases) = solution_lastInc_scal(0:NActivePhases-1,i,j,k)                 
   phi_grad = matmul(BMat,solution_current_elem)
   
   phi_flux   = homogenization_multiphase_getPhaseFlux  (phi_grad   ,1,cell)
   phi_source = homogenization_multiphase_getPhaseSource(phi_current,1,cell)
   
   residual_current_elem = 0.0_pReal
   residual_current_elem(1  ,1:NActivePhases) = residual_current_elem(1,1:NActivePhases) + &
                                               (phi_current(1:NActivePhases) - &
                                                phi_lastInc(1:NActivePhases) - &
                                                params%timeinc*phi_source(1:NActivePhases))
   residual_current_elem(1:4,1:NActivePhases) = residual_current_elem(1:4,1:NActivePhases) + &
                                                params%timeinc* &
                                                matmul(transpose(BMat),phi_flux(1:3,1:NActivePhases))
   residual_current_scal(0:NActivePhases-1,i  ,j  ,k  ) = &
   residual_current_scal(0:NActivePhases-1,i  ,j  ,k  ) + residual_current_elem(1,1:NActivePhases)
   residual_current_scal(0:NActivePhases-1,i+1,j  ,k  ) = &
   residual_current_scal(0:NActivePhases-1,i+1,j  ,k  ) + residual_current_elem(2,1:NActivePhases)
   residual_current_scal(0:NActivePhases-1,i  ,j+1,k  ) = &
   residual_current_scal(0:NActivePhases-1,i  ,j+1,k  ) + residual_current_elem(3,1:NActivePhases)
   residual_current_scal(0:NActivePhases-1,i  ,j  ,k+1) = &
   residual_current_scal(0:NActivePhases-1,i  ,j  ,k+1) + residual_current_elem(4,1:NActivePhases)
 enddo; enddo; enddo
 call DMDAVecRestoreArrayF90(da_local,solution_lastInc      ,solution_lastInc_scal,ierr)
 CHKERRQ(ierr)
 call DMDAVecRestoreArrayF90(da_local,solution_current_local,solution_current_scal,ierr)            ! put back to PETSc data structure
 CHKERRQ(ierr)
 call DMDAVecRestoreArrayF90(da_local,residual_current_local,residual_current_scal,ierr)
 CHKERRQ(ierr)
   
end subroutine spectral_multiphase_formResidual


!--------------------------------------------------------------------------------------------------
!> @brief forms the FEM stiffness matrix
!--------------------------------------------------------------------------------------------------
subroutine spectral_multiphase_formJacobian(da_local,solution_current_local,Jac_pre,Jac,dummy,ierr)
 use math, only: &
   math_identity2nd
 use homogenization_multiphase, only: &
   homogenization_multiphase_getPhaseFluxTangent, &
   homogenization_multiphase_getPhaseSourceTangent

 implicit none

 DM                   :: da_local
 Vec                  :: solution_current_local
 Mat                  :: Jac, Jac_pre
 MatStencil           :: col(4,4*NActivePhases)
 PetscScalar, pointer :: solution_current_scal(:,:,:,:)
 PetscScalar          :: k_ele(4*NActivePhases,4*NActivePhases)
 PetscScalar          :: phi_current(NActivePhases)
 real(pReal)          :: phi_sourceTangent(NActivePhases,NActivePhases), &
                         phi_fluxTangent  (NActivePhases,NActivePhases)
 PetscInt             :: i, j, k, cell, phaseI, phaseJ
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
   do phaseI = 0, NActivePhases-1 
     col(MatStencil_i,4*phaseI+1) = i
     col(MatStencil_j,4*phaseI+1) = j
     col(MatStencil_k,4*phaseI+1) = k
     col(MatStencil_c,4*phaseI+1) = phaseI
     col(MatStencil_i,4*phaseI+2) = i+1
     col(MatStencil_j,4*phaseI+2) = j
     col(MatStencil_k,4*phaseI+2) = k
     col(MatStencil_c,4*phaseI+2) = phaseI
     col(MatStencil_i,4*phaseI+3) = i
     col(MatStencil_j,4*phaseI+3) = j+1
     col(MatStencil_k,4*phaseI+3) = k
     col(MatStencil_c,4*phaseI+3) = phaseI
     col(MatStencil_i,4*phaseI+4) = i
     col(MatStencil_j,4*phaseI+4) = j
     col(MatStencil_k,4*phaseI+4) = k+1
     col(MatStencil_c,4*phaseI+4) = phaseI
   enddo

   phi_current = solution_current_scal(0:NActivePhases-1,i,j,k)                 
   phi_fluxTangent   = homogenization_multiphase_getPhaseFluxTangent  (           1,cell)
   phi_sourceTangent = homogenization_multiphase_getPhaseSourceTangent(phi_current,1,cell)

   k_ele = 0.0
   k_ele(1:4*NActivePhases:4,1:4*NActivePhases:4) = &
     k_ele(1:4*NActivePhases:4,1:4*NActivePhases:4) + &
     math_identity2nd(NActivePhases) - params%timeinc*phi_sourceTangent
   do phaseI = 0, NActivePhases-1; do phaseJ = 0, NActivePhases-1
     k_ele(4*phaseI+1:4*(phaseI+1),4*phaseJ+1:4*(phaseJ+1)) = &
       k_ele(4*phaseI+1:4*(phaseI+1),4*phaseJ+1:4*(phaseJ+1)) + &
       params%timeinc*phi_fluxTangent(phaseI+1,phaseJ+1)* &
       matmul(transpose(BMat),BMat)
   enddo; enddo  
   k_ele = transpose(k_ele)
   call MatSetValuesStencil(Jac,4*NActivePhases,col,4*NActivePhases,col,k_ele,ADD_VALUES,ierr)
   CHKERRQ(ierr)
 enddo; enddo; enddo
 call MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 call MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 call MatAssemblyBegin(Jac_pre,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 call MatAssemblyEnd(Jac_pre,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 call DMDAVecRestoreArrayF90(da_local,solution_current_local,solution_current_scal,ierr)
 CHKERRQ(ierr)
   
end subroutine spectral_multiphase_formJacobian


!--------------------------------------------------------------------------------------------------
!> @brief forwarding routine
!--------------------------------------------------------------------------------------------------
subroutine spectral_multiphase_forward()
 use spectral_utilities, only: &
   cutBack
 use homogenization_multiphase, only: &
   homogenization_multiphase_putPhaseFrac, &
   homogenization_multiphase_putInterfaceNormals

 implicit none
 DM :: da_local
 Vec :: solution_current_local
 PetscScalar, pointer :: solution_current_scal(:,:,:,:)
 integer(pInt) :: i, j, k, cell, phaseI, phaseJ
 real(pReal) :: interfaceNormal(3)
 PetscErrorCode :: ierr

 call SNESGetDM(multiphase_snes,da_local,ierr); CHKERRQ(ierr)
 if (cutBack) then 
   call VecCopy(solution_lastInc,solution_current,ierr); CHKERRQ(ierr)
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
     do phaseI = 1, NActivePhases-1
       do phaseJ = phaseI+1, NActivePhases
         interfaceNormal(1) = (solution_current_scal(phaseI-1,i+1,j  ,k  ) - &
                               solution_current_scal(phaseI-1,i-1,j  ,k  ))/delta(1) - &
                              (solution_current_scal(phaseJ-1,i+1,j  ,k  ) - &
                               solution_current_scal(phaseJ-1,i-1,j  ,k  ))/delta(1)
         interfaceNormal(2) = (solution_current_scal(phaseI-1,i  ,j+1,k  ) - &
                               solution_current_scal(phaseI-1,i  ,j-1,k  ))/delta(2) - &
                              (solution_current_scal(phaseJ-1,i  ,j+1,k  ) - &
                               solution_current_scal(phaseJ-1,i  ,j-1,k  ))/delta(2)
         interfaceNormal(3) = (solution_current_scal(phaseI-1,i  ,j  ,k+1) - &
                               solution_current_scal(phaseI-1,i  ,j  ,k-1))/delta(3) - &
                              (solution_current_scal(phaseJ-1,i  ,j  ,k+1) - &
                               solution_current_scal(phaseJ-1,i  ,j  ,k-1))/delta(3)
         if (norm2(interfaceNormal) > 0.0_pReal) then
           interfaceNormal = interfaceNormal/norm2(interfaceNormal)
         else
           interfaceNormal = 0.0_pReal  
         endif  
         call homogenization_multiphase_putInterfaceNormals(interfaceNormal,phaseI,phaseJ,1,cell)
       enddo
     enddo    
     call homogenization_multiphase_putPhaseFrac(solution_current_scal(0:NActivePhases-1,i,j,k),1,cell)
   enddo; enddo; enddo
   call DMDAVecRestoreArrayF90(da_local,solution_current_local,solution_current_scal,ierr)
   CHKERRQ(ierr)
   call DMRestoreLocalVector(da_local,solution_current_local,ierr); CHKERRQ(ierr)
 else
   call VecCopy(solution_current,solution_lastInc,ierr); CHKERRQ(ierr)                              ! also store current rate in rate_lastInc to recover upon cutback!!
 endif

end subroutine spectral_multiphase_forward

end module spectral_multiphase