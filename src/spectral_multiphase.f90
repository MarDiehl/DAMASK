!--------------------------------------------------------------------------------------------------
!> @author Chuanlai Liu, Shanghai Jiao Tong University
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Multiphase PETSc solver
!--------------------------------------------------------------------------------------------------
module spectral_multiphase
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
#include <petsc/finclude/petsc.h90>

!--------------------------------------------------------------------------------------------------
! derived types
 type(tSolutionParams), private :: params

!--------------------------------------------------------------------------------------------------
! PETSc data
 SNES, private :: multiphase_snes
 Vec,  private :: solution_current, solution_lastInc

!--------------------------------------------------------------------------------------------------
! common pointwise data
 integer(pInt) , private :: NPhases
 real   (pReal), private :: delta(3),  nablaMat(7)
 PetscInt      , private :: xstart,ystart,zstart,xend,yend,zend

 public :: &
   spectral_multiphase_init, &
   spectral_multiphase_solution, &
   spectral_multiphase_forward, &
   spectral_multiphase_destroy
 external :: &
   VecDestroy, &
   DMDestroy, &
   DMDACreate3D, &
   DMCreateGlobalVector, &
   DMDASNESSetFunctionLocal, &
   PETScFinalize, &
   SNESDestroy, &
   SNESGetNumberFunctionEvals, &
   SNESGetIterationNumber, &
   SNESSolve, &
   SNESSetDM, &
   SNESGetConvergedReason, &
   SNESSetConvergenceTest, &
   SNESSetFromOptions, &
   SNESCreate, &
   MPI_Abort, &
   MPI_Bcast, &
   MPI_Allreduce

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
   material_Nhomogenization, &
   phasefracMapping, &
   phasefrac
   
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
 NPhases = homogenization_maxNgrains
 if (NPhases < 2_pInt) call IO_error(211_pInt,ext_msg='# of phases < 2 in spectral_multiphase')

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
        NPhases, 1, &                                                                               ! #dof, ghost boundary width
        grid (1),grid(2),localK, &                                                                  ! local grid
        multiphase_grid,ierr)                                                                       ! handle, error
 CHKERRQ(ierr)
 call SNESSetDM(multiphase_snes,multiphase_grid,ierr); CHKERRQ(ierr)
 call DMCreateGlobalVector(multiphase_grid,solution_current,ierr); CHKERRQ(ierr)                    ! current phase field vector 
 call DMCreateGlobalVector(multiphase_grid,solution_lastInc,ierr); CHKERRQ(ierr)                    ! last increment phase field vector 
 call DMSNESSetFunctionLocal(multiphase_grid,spectral_multiphase_formResidual,PETSC_NULL_OBJECT,ierr)! residual vector of same shape as solution vector
 CHKERRQ(ierr) 
 call DMSNESSetJacobianLocal(multiphase_grid,spectral_multiphase_formJacobian,PETSC_NULL_OBJECT,ierr)! function to evaluate stiffness matrix
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
 CHKERRQ(ierr)                   !< get the data out of PETSc to work with
 do k = zstart, zend;  do j = ystart, yend;  do i = xstart, xend
   cell = cell + 1_pInt
   homog = material_homog(1,cell)
   do phase = 1, NPhases
     x_scal(phase-1,i,j,k) =  phasefrac(phase,homog)%p(phasefracMapping(homog)%p(1,cell))
   enddo
 enddo; enddo; enddo
 call DMDAVecRestoreArrayF90(multiphase_grid,solution_current,x_scal,ierr)
 CHKERRQ(ierr)
  
!--------------------------------------------------------------------------------------------------
! init FEM                 
 delta = geomSize/real(grid,pReal)                                                                  ! grid spacing
 nablaMat = [ 1.0/delta(1)/delta(1), 1.0/delta(1)/delta(1), &
              1.0/delta(2)/delta(2), 1.0/delta(2)/delta(2), &
              1.0/delta(3)/delta(3), 1.0/delta(3)/delta(3), &
             -2.0/delta(1)/delta(1) -2.0/delta(2)/delta(2) -2.0/delta(3)/delta(3)]
   
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
   homogenization_multiphase_putPhaseFrac

 implicit none

!--------------------------------------------------------------------------------------------------
! input data for solution
 real(pReal), intent(in) :: &
   timeinc, &                                                                                       !< increment in time for current solution
   timeinc_old                                                                                      !< increment in time of last increment
 DM :: da_local
 PetscScalar, pointer :: x_scal_current(:,:,:,:), x_scal_lastInc(:,:,:,:)
 real(pReal), dimension(NPhases) :: phi_current, phi_lastInc, phi_stagInc
 real(pReal), dimension(NPhases) :: stagNorm, minPhi, maxPhi, avgPhi
 integer(pInt) :: i, j, k, cell, phase, homog

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
 call SNESSolve(multiphase_snes,PETSC_NULL_OBJECT,solution_current,ierr)
 CHKERRQ(ierr)
 call SNESGetConvergedReason(multiphase_snes,reason,ierr); CHKERRQ(ierr)
 spectral_multiphase_solution%converged = reason > 0

!--------------------------------------------------------------------------------------------------
! get local fields 
 call SNESGetDM(multiphase_snes,da_local,ierr); CHKERRQ(ierr)
 call DMDAVecGetArrayF90(da_local,solution_current,x_scal_current,ierr)
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
   do phase = 1, NPhases
     phi_stagInc(phase) = phasefrac(phase,homog)%p(phasefracMapping(homog)%p(1,cell))
     phi_current(phase) = x_scal_current(phase-1,i,j,k)
     phi_lastInc(phase) = x_scal_lastInc(phase-1,i,j,k)
   enddo
   do phase = 1, NPhases
     minPhi(phase) = min(minPhi(phase),phi_current(phase))
     maxPhi(phase) = max(maxPhi(phase),phi_current(phase))
     avgPhi(phase) = avgPhi(phase) + phi_current(phase)
     stagNorm(phase) = max(stagNorm(phase),abs(phi_current(phase) - phi_stagInc(phase)))
   enddo
   call homogenization_multiphase_putPhaseFrac(phi_current,1,cell)
 enddo; enddo; enddo
 call MPI_Allreduce(MPI_IN_PLACE,minPhi    ,NPhases,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD,ierr)
 call MPI_Allreduce(MPI_IN_PLACE,maxPhi    ,NPhases,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD,ierr)
 call MPI_Allreduce(MPI_IN_PLACE,avgPhi    ,NPhases,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
 call MPI_Allreduce(MPI_IN_PLACE,stagNorm,  NPhases,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD,ierr)

!--------------------------------------------------------------------------------------------------
! restore local fields 
 call DMDAVecRestoreArrayF90(da_local,solution_lastInc,x_scal_lastInc,ierr)
 CHKERRQ(ierr)
 call DMDAVecRestoreArrayF90(da_local,solution_current,x_scal_current,ierr)
 CHKERRQ(ierr)
 
!--------------------------------------------------------------------------------------------------
! reporting 
 spectral_multiphase_solution%stagConverged = all(stagNorm < err_phasefr_tolAbs)

 if (spectral_multiphase_solution%converged) then 
   write(6,'(/,a)') ' ... multiphase converged .....................................'
   do phase = 1, NPhases
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
 use numerics, only: &
   err_phasefr_tolAbs
 use homogenization_multiphase, only: &
   homogenization_multiphase_getPhaseFlux, &
   homogenization_multiphase_getPhaseSource, &
   homogenization_multiphase_getPhaseMobility

 implicit none

 DM                   :: da_local
 Vec                  :: solution_current_local, &
                         residual_current_local
 PetscScalar, pointer :: solution_current_scal(:,:,:,:), &
                         solution_lastInc_scal(:,:,:,:), &
                         residual_current_scal(:,:,:,:)
 PetscScalar          :: phi_current(NPhases), &
                         phi_lastInc(NPhases), &
                         phi_nabla  (NPhases), &
                         phi_fluxRate(NPhases), &
                         phi_bulkRate(NPhases), &
                         phi_intfRate_prediction(NPhases), &
                         phi_intfRate_correction(NPhases), &
                         phi_mobility(NPhases,NPhases)
 logical              :: phi_active(NPhases)
 PetscInt             :: i, j, k, cell, phaseI, phaseJ
 PetscObject          :: dummy
 PetscErrorCode       :: ierr

!--------------------------------------------------------------------------------------------------
! constructing residual
 call VecSet(residual_current_local,0.0,ierr);CHKERRQ(ierr)                                         ! set residual to zero
 call DMDAVecGetArrayF90(da_local,residual_current_local,residual_current_scal,ierr)
 CHKERRQ(ierr)                                                                                      ! residual_current_local(PETSc) -> residual_current_scal(our) with layout from da_local
 call DMDAVecGetArrayF90(da_local,solution_current_local,solution_current_scal,ierr)                ! sol_cur_local(PETSc, phase fracs) -> sol_cut_scal(our) ...
 CHKERRQ(ierr)
 call DMDAVecGetArrayF90(da_local,solution_lastInc,      solution_lastInc_scal,ierr)                
 CHKERRQ(ierr)
 cell = 0
 do k = zstart, zend; do j = ystart, yend; do i = xstart, xend
   cell = cell + 1
   phi_current(1:NPhases) = &
     solution_current_scal(0:NPhases-1,i,j,k)                 
   phi_lastInc(1:NPhases) = &
     solution_lastInc_scal(0:NPhases-1,i,j,k)                 
   phi_nabla  = &
     (solution_current_scal(0:NPhases-1,i+1,j  ,k  ) - &
      solution_current_scal(0:NPhases-1,i  ,j  ,k  ))/delta(1)/delta(1) - &
     (solution_current_scal(0:NPhases-1,i  ,j  ,k  ) - &
      solution_current_scal(0:NPhases-1,i-1,j  ,k  ))/delta(1)/delta(1) + &
     (solution_current_scal(0:NPhases-1,i  ,j+1,k  ) - &
      solution_current_scal(0:NPhases-1,i  ,j  ,k  ))/delta(2)/delta(2) - &
     (solution_current_scal(0:NPhases-1,i  ,j  ,k  ) - &
      solution_current_scal(0:NPhases-1,i  ,j-1,k  ))/delta(2)/delta(2) + &
     (solution_current_scal(0:NPhases-1,i  ,j  ,k+1) - &
      solution_current_scal(0:NPhases-1,i  ,j  ,k  ))/delta(3)/delta(3) - &
     (solution_current_scal(0:NPhases-1,i  ,j  ,k  ) - &
      solution_current_scal(0:NPhases-1,i  ,j  ,k-1))/delta(3)/delta(3)

   phi_fluxRate = homogenization_multiphase_getPhaseFlux    (phi_nabla  ,1,cell)
   phi_bulkRate = homogenization_multiphase_getPhaseSource  (phi_current,1,cell)
   phi_mobility = homogenization_multiphase_getPhaseMobility(            1,cell) 
   
   phi_intfRate_prediction = phi_current - phi_lastInc
   do phaseI = 1, NPhases; do phaseJ = 1, NPhases
     phi_intfRate_prediction(phaseI) = &
       phi_intfRate_prediction(phaseI) - &
       params%timeinc* &
       phi_mobility(phaseI,phaseJ)* &
       ((phi_fluxRate(phaseI) + phi_bulkRate(phaseI)) - &
        (phi_fluxRate(phaseJ) + phi_bulkRate(phaseJ)))
   enddo; enddo
   phi_active =      (phi_current <= 0.0 + err_phasefr_tolAbs .and. phi_intfRate_prediction >= 0.0) &
                .or. (phi_current >= 1.0 - err_phasefr_tolAbs .and. phi_intfRate_prediction <= 0.0)
   
   phi_intfRate_correction = 0.0_pReal
   do phaseI = 1, NPhases
     if (.not. phi_active(phaseI)) then
       phi_intfRate_correction(phaseI) = phi_current(phaseI) - phi_lastInc(phaseI)
       do phaseJ = 1, NPhases
         if (.not. phi_active(phaseJ)) then
           phi_intfRate_correction(phaseI) = &
             phi_intfRate_correction(phaseI) - &
             params%timeinc* &
             phi_mobility(phaseI,phaseJ)* &
             ((phi_fluxRate(phaseI) + phi_bulkRate(phaseI)) - &
              (phi_fluxRate(phaseJ) + phi_bulkRate(phaseJ)))
         endif
       enddo 
     else
       phi_intfRate_correction(phaseI) = phi_intfRate_prediction(phaseI)
     endif
   enddo  
   residual_current_scal(0:NPhases-1,i,j,k) = phi_intfRate_correction
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
 use numerics, only: &
   err_phasefr_tolAbs
 use homogenization_multiphase, only: &
   homogenization_multiphase_getPhaseFlux, &
   homogenization_multiphase_getPhaseFluxTangent, &
   homogenization_multiphase_getPhaseSource, &
   homogenization_multiphase_getPhaseSourceTangent, &
   homogenization_multiphase_getPhaseMobility

 implicit none

 DM                   :: da_local
 Vec                  :: solution_current_local
 Mat                  :: Jac, Jac_pre
 MatStencil           :: col(4,7*NPhases), &
                         row(4,  NPhases)   
 PetscScalar          :: k_correction(NPhases,7*NPhases)
 PetscScalar, pointer :: solution_current_scal(:,:,:,:), &
                         solution_lastInc_scal(:,:,:,:)
 PetscScalar          :: phi_current(NPhases), &
                         phi_lastInc(NPhases), &
                         phi_nabla  (NPhases), &
                         phi_fluxRate(NPhases), &
                         phi_bulkRate(NPhases), &
                         phi_intfRate_prediction(NPhases), &
                         phi_fluxRateTangent(NPhases,NPhases), &
                         phi_bulkRateTangent(NPhases,NPhases), &
                         phi_mobility(NPhases,NPhases)
 logical              :: phi_active(NPhases)
 PetscInt             :: i, j, k, cell, phaseI, phaseJ, phaseK
 PetscObject          :: dummy
 PetscErrorCode       :: ierr

!--------------------------------------------------------------------------------------------------
! constructing residual
 call MatSetOption(Jac,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr); CHKERRQ(ierr)
 call MatSetOption(Jac,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr); CHKERRQ(ierr)
 call MatZeroEntries(Jac,ierr); CHKERRQ(ierr)
 call DMDAVecGetArrayF90(da_local,solution_current_local,solution_current_scal,ierr)                ! sol_cur_local(PETSc, phase fracs) -> sol_cut_scal(our) ...
 CHKERRQ(ierr)
 call DMDAVecGetArrayF90(da_local,solution_lastInc,      solution_lastInc_scal,ierr)                
 CHKERRQ(ierr)
 cell = 0
 do k = zstart, zend; do j = ystart, yend; do i = xstart, xend
   cell = cell + 1
   do phaseI = 0, NPhases-1 
     col(MatStencil_i,7*phaseI+1) = i+1
     col(MatStencil_j,7*phaseI+1) = j
     col(MatStencil_k,7*phaseI+1) = k
     col(MatStencil_c,7*phaseI+1) = phaseI
     col(MatStencil_i,7*phaseI+2) = i-1
     col(MatStencil_j,7*phaseI+2) = j
     col(MatStencil_k,7*phaseI+2) = k
     col(MatStencil_c,7*phaseI+2) = phaseI
     col(MatStencil_i,7*phaseI+3) = i
     col(MatStencil_j,7*phaseI+3) = j+1
     col(MatStencil_k,7*phaseI+3) = k
     col(MatStencil_c,7*phaseI+3) = phaseI
     col(MatStencil_i,7*phaseI+4) = i
     col(MatStencil_j,7*phaseI+4) = j-1
     col(MatStencil_k,7*phaseI+4) = k
     col(MatStencil_c,7*phaseI+4) = phaseI
     col(MatStencil_i,7*phaseI+5) = i
     col(MatStencil_j,7*phaseI+5) = j
     col(MatStencil_k,7*phaseI+5) = k+1
     col(MatStencil_c,7*phaseI+5) = phaseI
     col(MatStencil_i,7*phaseI+6) = i
     col(MatStencil_j,7*phaseI+6) = j
     col(MatStencil_k,7*phaseI+6) = k-1
     col(MatStencil_c,7*phaseI+6) = phaseI
     col(MatStencil_i,7*phaseI+7) = i
     col(MatStencil_j,7*phaseI+7) = j
     col(MatStencil_k,7*phaseI+7) = k
     col(MatStencil_c,7*phaseI+7) = phaseI
     row(MatStencil_i,  phaseI+1) = i
     row(MatStencil_j,  phaseI+1) = j
     row(MatStencil_k,  phaseI+1) = k
     row(MatStencil_c,  phaseI+1) = phaseI
   enddo

   phi_current(1:NPhases) = &
     solution_current_scal(0:NPhases-1,i,j,k)                 
   phi_lastInc(1:NPhases) = &
     solution_lastInc_scal(0:NPhases-1,i,j,k)                 
   phi_nabla  = &
     (solution_current_scal(0:NPhases-1,i+1,j  ,k  ) - &
      solution_current_scal(0:NPhases-1,i  ,j  ,k  ))/delta(1)/delta(1) - &
     (solution_current_scal(0:NPhases-1,i  ,j  ,k  ) - &
      solution_current_scal(0:NPhases-1,i-1,j  ,k  ))/delta(1)/delta(1) + &
     (solution_current_scal(0:NPhases-1,i  ,j+1,k  ) - &
      solution_current_scal(0:NPhases-1,i  ,j  ,k  ))/delta(2)/delta(2) - &
     (solution_current_scal(0:NPhases-1,i  ,j  ,k  ) - &
      solution_current_scal(0:NPhases-1,i  ,j-1,k  ))/delta(2)/delta(2) + &
     (solution_current_scal(0:NPhases-1,i  ,j  ,k+1) - &
      solution_current_scal(0:NPhases-1,i  ,j  ,k  ))/delta(3)/delta(3) - &
     (solution_current_scal(0:NPhases-1,i  ,j  ,k  ) - &
      solution_current_scal(0:NPhases-1,i  ,j  ,k-1))/delta(3)/delta(3)

   phi_fluxRate = homogenization_multiphase_getPhaseFlux  (phi_nabla  ,1,cell)
   phi_bulkRate = homogenization_multiphase_getPhaseSource(phi_current,1,cell)
   phi_mobility = homogenization_multiphase_getPhaseMobility(          1,cell) 
   
   phi_intfRate_prediction = phi_current - phi_lastInc
   do phaseI = 1, NPhases; do phaseJ = 1, NPhases
     phi_intfRate_prediction(phaseI) = &
       phi_intfRate_prediction(phaseI) - &
       params%timeinc* &
       phi_mobility(phaseI,phaseJ)* &
       ((phi_fluxRate(phaseI) + phi_bulkRate(phaseI)) - &
        (phi_fluxRate(phaseJ) + phi_bulkRate(phaseJ)))
   enddo; enddo
   phi_active =      (phi_current <= 0.0 + err_phasefr_tolAbs .and. phi_intfRate_prediction >= 0.0) &
                .or. (phi_current >= 1.0 - err_phasefr_tolAbs .and. phi_intfRate_prediction <= 0.0)
   
   phi_fluxRateTangent = homogenization_multiphase_getPhaseFluxTangent  (1,cell)
   phi_bulkRateTangent = homogenization_multiphase_getPhaseSourceTangent(1,cell)
   
   k_correction = 0.0_pReal
   do phaseI = 0, NPhases-1
     if (.not. phi_active(phaseI+1)) then
       k_correction(phaseI+1,7*phaseI+7) = 1.0
       do phaseJ = 0, NPhases-1
         do phaseK = 0, NPhases-1
           if ((.not. phi_active(phaseJ+1)) .and. &
               (.not. phi_active(phaseK+1))) then
             k_correction(phaseI+1,7*phaseK+7) = &
               k_correction(phaseI+1,7*phaseK+7) - &
               params%timeinc* &
               phi_mobility(phaseI+1,phaseJ+1)* &
               (phi_bulkRateTangent(phaseI+1,phaseK+1) - &
                phi_bulkRateTangent(phaseJ+1,phaseK+1))
             k_correction(phaseI+1,7*phaseK+1:7*phaseK+7) = &
               k_correction(phaseI+1,7*phaseK+1:7*phaseK+7) - &
               params%timeinc* &
               phi_mobility(phaseI+1,phaseJ+1)* &
               (phi_fluxRateTangent(phaseI+1,phaseK+1) - &
                phi_fluxRateTangent(phaseJ+1,phaseK+1))*nablaMat
           endif
         enddo
       enddo
     endif
   enddo      
   call MatSetValuesStencil(Jac,NPhases,row,7*NPhases,col,transpose(k_correction),ADD_VALUES,ierr)
   CHKERRQ(ierr)
 enddo; enddo; enddo
 call MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 call MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 call MatAssemblyBegin(Jac_pre,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 call MatAssemblyEnd(Jac_pre,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 call DMDAVecRestoreArrayF90(da_local,solution_lastInc      ,solution_lastInc_scal,ierr)
 CHKERRQ(ierr)
 call DMDAVecRestoreArrayF90(da_local,solution_current_local,solution_current_scal,ierr)            ! put back to PETSc data structure
 CHKERRQ(ierr)

end subroutine spectral_multiphase_formJacobian


!--------------------------------------------------------------------------------------------------
!> @brief forwarding routine
!--------------------------------------------------------------------------------------------------
subroutine spectral_multiphase_forward()
 use spectral_utilities, only: &
   cutBack
 use homogenization_multiphase, only: &
   homogenization_multiphase_putPhaseFrac

 implicit none
 DM :: da_local
 PetscScalar, pointer :: solution_current_scal(:,:,:,:)
 integer(pInt) :: i, j, k, cell
 PetscErrorCode :: ierr

 call SNESGetDM(multiphase_snes,da_local,ierr); CHKERRQ(ierr)
 if (cutBack) then 
   call VecCopy(solution_lastInc,solution_current,ierr); CHKERRQ(ierr)
   call DMDAVecGetArrayF90(da_local,solution_current,solution_current_scal,ierr)
   CHKERRQ(ierr)
   cell = 0
   do k = zstart, zend; do j = ystart, yend; do i = xstart, xend                                    ! walk through my (sub)domain
     cell = cell + 1
     call homogenization_multiphase_putPhaseFrac(solution_current_scal(0:NPhases-1,i,j,k),1,cell)
   enddo; enddo; enddo
   call DMDAVecRestoreArrayF90(da_local,solution_current,solution_current_scal,ierr)
   CHKERRQ(ierr)
 else
   call VecCopy(solution_current,solution_lastInc,ierr); CHKERRQ(ierr)                              ! also store current rate in rate_lastInc to recover upon cutback!!
 endif

end subroutine spectral_multiphase_forward

!--------------------------------------------------------------------------------------------------
!> @brief destroy routine
!--------------------------------------------------------------------------------------------------
subroutine spectral_multiphase_destroy()

 implicit none
 PetscErrorCode :: ierr

 call VecDestroy(solution_current,ierr); CHKERRQ(ierr)
 call VecDestroy(solution_lastInc,ierr); CHKERRQ(ierr)
 call SNESDestroy(multiphase_snes,ierr); CHKERRQ(ierr)

end subroutine spectral_multiphase_destroy

end module spectral_multiphase
