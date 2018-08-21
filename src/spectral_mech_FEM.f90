!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief FEM PETSc solver
!--------------------------------------------------------------------------------------------------
module spectral_mech_FEM
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

 implicit none
 private

 character (len=*), parameter, public :: &
   DAMASK_spectral_SolverFEM_label = 'fem'
   
!--------------------------------------------------------------------------------------------------
! derived types
 type(tSolutionParams), private :: params

!--------------------------------------------------------------------------------------------------
! PETSc data
 DM,   private :: mech_grid
 SNES, private :: mech_snes
 Vec,  private :: solution_current, solution_lastInc, solution_rate

!--------------------------------------------------------------------------------------------------
! common pointwise data
 real(pReal), private, dimension(:,:,:,:,:), allocatable ::  F_current, P_current, F_lastInc
 real(pReal), private :: detJ, delta(3), BMat(3,8), HGMat(8,8)
 PetscInt,    private :: xstart,ystart,zstart,xend,yend,zend

!--------------------------------------------------------------------------------------------------
! stress, stiffness and compliance average etc.
 real(pReal), private, dimension(3,3) :: &
   F_aim = math_I3, &
   F_aim_lastIter = math_I3, &
   F_aim_lastInc = math_I3, &
   P_av = 0.0_pReal, &
   F_aimDot=0.0_pReal
 character(len=1024), private :: incInfo   
 real(pReal), private, dimension(3,3,3,3) :: &
   C_volAvg = 0.0_pReal, &                                                                          !< current volume average stiffness 
   C_volAvgLastInc = 0.0_pReal, &                                                                   !< previous volume average stiffness
   C_minMaxAvg = 0.0_pReal, &                                                                       !< current (min+max)/2 stiffness
   S = 0.0_pReal                                                                                    !< current compliance (filled up with zeros)
 real(pReal), private :: err_stress, err_div
 logical, private :: ForwardData
 integer(pInt), private :: &
   totalIter = 0_pInt                                                                               !< total iteration in current increment
 real(pReal), private, dimension(3,3) :: mask_stress = 0.0_pReal

 public :: &
   FEM_init, &
   FEM_solution, &
   FEM_forward
 external :: &
   PETScErrorF

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields and fills them with data, potentially from restart info
!--------------------------------------------------------------------------------------------------
subroutine FEM_init
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use IO, only: &
   IO_intOut, &
   IO_read_realFile, &
   IO_timeStamp
 use debug, only: &
  debug_level, &
  debug_spectral, &
  debug_spectralRestart
 use FEsolving, only: &
   restartInc
 use numerics, only: &
   worldrank, &
   worldsize
 use DAMASK_interface, only: &
   getSolverJobName
 use homogenization, only: &
   materialpoint_F0
 use spectral_utilities, only: &
   utilities_constitutiveResponse, &
   utilities_updateIPcoords, &
   wgt
 use mesh, only: &
   geomSize, &
   grid, &
   grid3
 use math, only: &
   math_invSym3333
   
 implicit none
 PetscErrorCode       :: ierr
 real(pReal)          :: &
   temp33_Real(3,3) = 0.0_pReal, HGCoeff = 0e-2_pReal
 integer(pInt), dimension(:), allocatable :: localK  
 integer(pInt) :: proc
 character(len=1024) :: rankStr

 mainProcess: if (worldrank == 0_pInt) then
   write(6,'(/,a)') ' <<<+-  spectral_mech_FEM init  -+>>>'
   write(6,'(a15,a)')   ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

!--------------------------------------------------------------------------------------------------
! allocate global fields
 allocate (F_current (3,3,grid(1),grid(2),grid3),source = 0.0_pReal)
 allocate (P_current (3,3,grid(1),grid(2),grid3),source = 0.0_pReal)
 allocate (F_lastInc (3,3,grid(1),grid(2),grid3),source = 0.0_pReal)
    
!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
 call SNESCreate(PETSC_COMM_WORLD,mech_snes,ierr); CHKERRQ(ierr)
 call SNESSetOptionsPrefix(mech_snes,'mech_',ierr);CHKERRQ(ierr) 
 allocate(localK(worldsize), source = 0); localK(worldrank+1) = grid3
 do proc = 1, worldsize
   call MPI_Bcast(localK(proc),1,MPI_INTEGER,proc-1,PETSC_COMM_WORLD,ierr)
 enddo  
 call DMDACreate3d(PETSC_COMM_WORLD, &
        DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, &                         ! cut off stencil at boundary
        DMDA_STENCIL_BOX, &                                                                         ! Moore (26) neighborhood around central point
        grid(1),grid(2),grid(3), &                                                                  ! global grid
        1, 1, worldsize, &
        3, 1, &                                                                                     ! #dof (F tensor), ghost boundary width (domain overlap)
        [grid(1)],[grid(2)],localK, &                                                                 ! local grid
        mech_grid,ierr)                                                                             ! handle, error
 CHKERRQ(ierr)
 call DMDASetUniformCoordinates(mech_grid,0.0,geomSize(1),0.0,geomSize(2),0.0,geomSize(3),ierr)     ! set dimensions of grid
 CHKERRQ(ierr)
 call SNESSetDM(mech_snes,mech_grid,ierr); CHKERRQ(ierr)
 call DMsetFromOptions(mech_grid,ierr); CHKERRQ(ierr)
 call DMsetUp(mech_grid,ierr); CHKERRQ(ierr)
 call DMCreateGlobalVector(mech_grid,solution_current,ierr); CHKERRQ(ierr)                          ! current global displacement vector 
 call DMCreateGlobalVector(mech_grid,solution_lastInc,ierr); CHKERRQ(ierr)                          ! last increment global displacement vector 
 call DMCreateGlobalVector(mech_grid,solution_rate   ,ierr); CHKERRQ(ierr)                          ! current global velocity vector 
 call DMSNESSetFunctionLocal(mech_grid,FEM_formResidual,PETSC_NULL_SNES,ierr)                       ! residual vector of same shape as solution vector
 CHKERRQ(ierr) 
 call DMSNESSetJacobianLocal(mech_grid,FEM_formJacobian,PETSC_NULL_SNES,ierr)                       ! function to evaluate stiffness matrix
 CHKERRQ(ierr)
 call SNESSetConvergenceTest(mech_snes,FEM_converged,PETSC_NULL_SNES,PETSC_NULL_FUNCTION,ierr) 
 CHKERRQ(ierr)                                                                                      ! specify custom convergence check function "_converged"
 call SNESSetMaxLinearSolveFailures(mech_snes, huge(1), ierr); CHKERRQ(ierr)                        ! ignore linear solve failures 
 call SNESSetFromOptions(mech_snes,ierr); CHKERRQ(ierr)                                             ! pull it all together with additional cli arguments

!--------------------------------------------------------------------------------------------------
! init fields                 
 call VecSet(solution_current,0.0,ierr);CHKERRQ(ierr)
 call VecSet(solution_lastInc,0.0,ierr);CHKERRQ(ierr)
 call VecSet(solution_rate   ,0.0,ierr);CHKERRQ(ierr)
 call DMDAGetCorners(mech_grid,xstart,ystart,zstart,xend,yend,zend,ierr)                            ! local grid extent
 CHKERRQ(ierr) 
 xend = xstart+xend-1
 yend = ystart+yend-1
 zend = zstart+zend-1
 delta = geomSize/real(grid,pReal)                                                                  ! grid spacing
 detJ = product(delta)                                                                              ! cell volume 
 BMat = reshape(real([-1.0_pReal/delta(1),-1.0_pReal/delta(2),-1.0_pReal/delta(3), &
                       1.0_pReal/delta(1),-1.0_pReal/delta(2),-1.0_pReal/delta(3), &
                      -1.0_pReal/delta(1), 1.0_pReal/delta(2),-1.0_pReal/delta(3), &
                       1.0_pReal/delta(1), 1.0_pReal/delta(2),-1.0_pReal/delta(3), &
                      -1.0_pReal/delta(1),-1.0_pReal/delta(2), 1.0_pReal/delta(3), &
                       1.0_pReal/delta(1),-1.0_pReal/delta(2), 1.0_pReal/delta(3), &
                      -1.0_pReal/delta(1), 1.0_pReal/delta(2), 1.0_pReal/delta(3), &
                       1.0_pReal/delta(1), 1.0_pReal/delta(2), 1.0_pReal/delta(3)],pReal), [3,8])/4.0_pReal                          ! shape function derivative matrix 
 HGMat = matmul(transpose(reshape([ 1.0_pReal, 1.0_pReal, 1.0_pReal,-1.0_pReal, &
                                    1.0_pReal,-1.0_pReal,-1.0_pReal, 1.0_pReal, &
                                   -1.0_pReal, 1.0_pReal,-1.0_pReal, 1.0_pReal, &
                                   -1.0_pReal,-1.0_pReal, 1.0_pReal,-1.0_pReal, &
                                   -1.0_pReal,-1.0_pReal, 1.0_pReal, 1.0_pReal, &
                                   -1.0_pReal, 1.0_pReal,-1.0_pReal,-1.0_pReal, &
                                    1.0_pReal,-1.0_pReal,-1.0_pReal,-1.0_pReal, &
                                    1.0_pReal, 1.0_pReal, 1.0_pReal, 1.0_pReal], [4,8])), &
                reshape([ 1.0_pReal, 1.0_pReal, 1.0_pReal,-1.0_pReal, &
                          1.0_pReal,-1.0_pReal,-1.0_pReal, 1.0_pReal, &
                         -1.0_pReal, 1.0_pReal,-1.0_pReal, 1.0_pReal, &
                         -1.0_pReal,-1.0_pReal, 1.0_pReal,-1.0_pReal, &
                         -1.0_pReal,-1.0_pReal, 1.0_pReal, 1.0_pReal, &
                         -1.0_pReal, 1.0_pReal,-1.0_pReal,-1.0_pReal, &
                          1.0_pReal,-1.0_pReal,-1.0_pReal,-1.0_pReal, &
                          1.0_pReal, 1.0_pReal, 1.0_pReal, 1.0_pReal], [4,8]))* &
         HGCoeff*(delta(1)*delta(2) + delta(2)*delta(3) + delta(3)*delta(1))/16.0_pReal             ! hourglass stabilisation matrix
 
 restart: if (restartInc > 0_pInt) then                                                     
   if (iand(debug_level(debug_spectral),debug_spectralRestart) /= 0) then
     write(6,'(/,a,'//IO_intOut(restartInc)//',a)') &
     'reading values of increment ', restartInc, ' from file'
     flush(6)
   endif
   write(rankStr,'(a1,i0)')'_',worldrank
   call IO_read_realFile(777,'F'//trim(rankStr),trim(getSolverJobName()),size(F_current))
   read (777,rec=1) F_current; close (777)
   call IO_read_realFile(777,'F_lastInc'//trim(rankStr),trim(getSolverJobName()),size(F_lastInc))
   read (777,rec=1) F_lastInc; close (777)
   call IO_read_realFile(777,'F_aimDot',trim(getSolverJobName()),size(F_aimDot))
   read (777,rec=1) F_aimDot; close (777)
   F_aim         = sum(sum(sum(F_current,dim=5),dim=4),dim=3) * wgt                                 ! average of F
   F_aim_lastInc = sum(sum(sum(F_lastInc,dim=5),dim=4),dim=3) * wgt                                 ! average of F_lastInc 
 elseif (restartInc == 0_pInt) then restart
   F_lastInc = spread(spread(spread(math_I3,3,grid(1)),4,grid(2)),5,grid3)                          ! initialize to identity
   F_current = spread(spread(spread(math_I3,3,grid(1)),4,grid(2)),5,grid3)
 endif restart

 materialpoint_F0 = reshape(F_lastInc, [3,3,1,product(grid(1:2))*grid3])                            ! set starting condition for materialpoint_stressAndItsTangent
 call Utilities_updateIPcoords(F_current)
 call Utilities_constitutiveResponse(P_current,temp33_Real,C_volAvg,C_minMaxAvg, &                  ! stress field, stress avg, global average of stiffness and (min+max)/2
                                     F_current, &                                                   ! target F
                                     0.0_pReal, &                                                   ! time increment
                                     math_I3)                                                       ! no rotation of boundary condition

 restartRead: if (restartInc > 1_pInt) then                                                    
   if (iand(debug_level(debug_spectral),debug_spectralRestart)/= 0 .and. worldrank == 0_pInt) &
     write(6,'(/,a,'//IO_intOut(restartInc-1_pInt)//',a)') &
     'reading more values of increment', restartInc - 1_pInt, 'from file'
   flush(6)
   call IO_read_realFile(777,'C_volAvg',trim(getSolverJobName()),size(C_volAvg))
   read (777,rec=1) C_volAvg
   close (777)
   call IO_read_realFile(777,'C_volAvgLastInc',trim(getSolverJobName()),size(C_volAvgLastInc))
   read (777,rec=1) C_volAvgLastInc
   close (777)
   call IO_read_realFile(777,'C_ref',trim(getSolverJobName()),size(C_minMaxAvg))
   read (777,rec=1) C_minMaxAvg
   close (777)
 endif restartRead
   
end subroutine FEM_init
  
!--------------------------------------------------------------------------------------------------
!> @brief solution for the FEM scheme with internal iterations
!--------------------------------------------------------------------------------------------------
type(tSolutionState) function &
  FEM_solution(incInfoIn,timeinc,timeinc_old,loadCaseTime,stress_BC,rotation)
 use IO, only: &
   IO_error
 use spectral_utilities, only: &
   tTensorBoundaryCondition, &
   utilities_maskedCompliance
 use FEsolving, only: &
   terminallyIll

 implicit none

!--------------------------------------------------------------------------------------------------
! input data for solution
 real(pReal), intent(in) :: &
   timeinc, &                                                                                       !< increment in time for current solution
   timeinc_old, &                                                                                   !< increment in time of last increment
   loadCaseTime                                                                                     !< remaining time of current load case
 type(tTensorBoundaryCondition),      intent(in) :: &
   stress_BC
 character(len=*), intent(in) :: &
   incInfoIn
 real(pReal), dimension(3,3), intent(in) :: rotation
 
!--------------------------------------------------------------------------------------------------
! PETSc Data
 PetscErrorCode :: ierr   
 SNESConvergedReason :: reason
 incInfo = incInfoIn

!--------------------------------------------------------------------------------------------------
! update stiffness (and gamma operator)
 S = Utilities_maskedCompliance(rotation,stress_BC%maskLogical,C_volAvg) 
 
!--------------------------------------------------------------------------------------------------
! set module wide availabe data 
 mask_stress       = stress_BC%maskFloat
 params%stress_BC  = stress_BC%values
 params%rotation_BC= rotation
 params%timeinc    = timeinc
 params%timeincOld = timeinc_old

!--------------------------------------------------------------------------------------------------
! solve BVP 
 call SNESSolve(mech_snes,PETSC_NULL_VEC,solution_current,ierr)
 CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! check convergence
 FEM_solution%iterationsNeeded = totalIter
 FEM_solution%termIll = terminallyIll
 terminallyIll = .false.
 call SNESGetConvergedReason(mech_snes,reason,ierr);CHKERRQ(ierr)
 FEM_solution%converged = reason > 0
 if (reason == SNES_DIVERGED_FNORM_NAN) call IO_error(893_pInt)

end function FEM_solution


!--------------------------------------------------------------------------------------------------
!> @brief forms the residual vector
!--------------------------------------------------------------------------------------------------
subroutine FEM_formResidual(da_local,x_local,f_local,dummy,ierr)
 use numerics, only: &
   itmax, &
   itmin
 use numerics, only: &
   worldrank
 use mesh, only: &
   grid
 use math, only: &
   math_rotate_backward33, &
   math_transpose33, &
   math_mul3333xx33
 use debug, only: &
   debug_level, &
   debug_spectral, &
   debug_spectralRotation
 use spectral_utilities, only: &
   utilities_constitutiveResponse
 use IO, only: &
   IO_intOut 
 use FEsolving, only: &
   terminallyIll
 use homogenization, only: &
   materialpoint_dPdF

 implicit none
 DM                   :: da_local
 Vec                  :: x_local, f_local
 PetscScalar, pointer :: x_scal(:,:,:,:), f_scal(:,:,:,:)
 PetscScalar          :: x_elem(8,3),  f_elem(8,3)
 PetscInt             :: i, ii, j, jj, k, kk, ctr, ele
 PetscInt             :: PETScIter, nfuncs
 PetscObject          :: dummy
 PetscErrorCode       :: ierr

 call SNESGetNumberFunctionEvals(mech_snes,nfuncs,ierr); CHKERRQ(ierr)
 call SNESGetIterationNumber(mech_snes,PETScIter,ierr); CHKERRQ(ierr)

 if(nfuncs== 0 .and. PETScIter == 0) totalIter = -1_pInt                                            ! new increment
 newIteration: if (totalIter <= PETScIter) then
!--------------------------------------------------------------------------------------------------
! report begin of new iteration
   totalIter = totalIter + 1_pInt
   if (worldrank == 0_pInt) then
     write(6,'(1x,a,3(a,'//IO_intOut(itmax)//'))') trim(incInfo), &
                    ' @ Iteration ', itmin, '≤',totalIter, '≤', itmax
     if (iand(debug_level(debug_spectral),debug_spectralRotation) /= 0) &
       write(6,'(/,a,/,3(3(f12.7,1x)/))',advance='no') ' deformation gradient aim (lab) =', &
                                   math_transpose33(math_rotate_backward33(F_aim,params%rotation_BC))
     write(6,'(/,a,/,3(3(f12.7,1x)/))',advance='no') ' deformation gradient aim =', &
                                 math_transpose33(F_aim)
     flush(6)
   endif
 endif newIteration

!--------------------------------------------------------------------------------------------------
! get deformation gradient
 call DMDAVecGetArrayF90(da_local,x_local,x_scal,ierr);CHKERRQ(ierr)
 do k = zstart, zend; do j = ystart, yend; do i = xstart, xend
   ctr = 0
   do kk = 0, 1; do jj = 0, 1; do ii = 0, 1
     ctr = ctr + 1
     x_elem(ctr,1:3) = x_scal(0:2,i+ii,j+jj,k+kk)
   enddo; enddo; enddo
   ii = i-xstart+1; jj = j-ystart+1; kk = k-zstart+1
   F_current(1:3,1:3,ii,jj,kk) = F_aim + transpose(matmul(BMat,x_elem)) 
 enddo; enddo; enddo
 call DMDAVecRestoreArrayF90(da_local,x_local,x_scal,ierr);CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! evaluate constitutive response
 call Utilities_constitutiveResponse(P_current,P_av,C_volAvg,C_minmaxAvg, &
                                     F_current,params%timeinc,params%rotation_BC)
 call MPI_Allreduce(MPI_IN_PLACE,terminallyIll,1,MPI_LOGICAL,MPI_LOR,PETSC_COMM_WORLD,ierr)
 ForwardData = .false.
  
!--------------------------------------------------------------------------------------------------
! stress BC handling
 F_aim_lastIter = F_aim
 F_aim = F_aim - math_mul3333xx33(S, ((P_av - params%stress_BC)))                                        ! S = 0.0 for no bc
 err_stress = maxval(abs(mask_stress * (P_av - params%stress_BC)))                                       ! mask = 0.0 for no bc

!--------------------------------------------------------------------------------------------------
! constructing residual
 call VecSet(f_local,0.0,ierr);CHKERRQ(ierr)
 call DMDAVecGetArrayF90(da_local,f_local,f_scal,ierr);CHKERRQ(ierr)
 call DMDAVecGetArrayF90(da_local,x_local,x_scal,ierr);CHKERRQ(ierr)
 ele = 0
 do k = zstart, zend; do j = ystart, yend; do i = xstart, xend
   ctr = 0
   do kk = 0, 1; do jj = 0, 1; do ii = 0, 1
     ctr = ctr + 1
     x_elem(ctr,1:3) = x_scal(0:2,i+ii,j+jj,k+kk)
   enddo; enddo; enddo
   ii = i-xstart+1; jj = j-ystart+1; kk = k-zstart+1
   ele = ele + 1
   f_elem = matmul(transpose(BMat),transpose(P_current(1:3,1:3,ii,jj,kk)))*detJ + &
            matmul(HGMat,x_elem)*(materialpoint_dPdF(1,1,1,1,1,ele) + &
                                  materialpoint_dPdF(2,2,2,2,1,ele) + &
                                  materialpoint_dPdF(3,3,3,3,1,ele))/3.0_pReal
   ctr = 0
   do kk = 0, 1; do jj = 0, 1; do ii = 0, 1
     ctr = ctr + 1
     f_scal(0:2,i+ii,j+jj,k+kk) = f_scal(0:2,i+ii,j+jj,k+kk) + f_elem(ctr,1:3)
   enddo; enddo; enddo
 enddo; enddo; enddo
 call DMDAVecRestoreArrayF90(da_local,x_local,x_scal,ierr);CHKERRQ(ierr)
 call DMDAVecRestoreArrayF90(da_local,f_local,f_scal,ierr);CHKERRQ(ierr)
 
!--------------------------------------------------------------------------------------------------
! applying boundary conditions
 call DMDAVecGetArrayF90(da_local,f_local,f_scal,ierr);CHKERRQ(ierr)
 if (zstart == 0) then         
   f_scal(0:2,xstart,ystart,zstart) = 0.0 
   f_scal(0:2,xend+1,ystart,zstart) = 0.0 
   f_scal(0:2,xstart,yend+1,zstart) = 0.0 
   f_scal(0:2,xend+1,yend+1,zstart) = 0.0 
 endif
 if (zend + 1 == grid(3)) then
   f_scal(0:2,xstart,ystart,zend+1) = 0.0 
   f_scal(0:2,xend+1,ystart,zend+1) = 0.0 
   f_scal(0:2,xstart,yend+1,zend+1) = 0.0 
   f_scal(0:2,xend+1,yend+1,zend+1) = 0.0 
 endif
 call DMDAVecRestoreArrayF90(da_local,f_local,f_scal,ierr);CHKERRQ(ierr) 
 
end subroutine FEM_formResidual


!--------------------------------------------------------------------------------------------------
!> @brief forms the FEM stiffness matrix
!--------------------------------------------------------------------------------------------------
subroutine FEM_formJacobian(da_local,x_local,Jac_pre,Jac,dummy,ierr)
 use mesh, only: &
   mesh_ipCoordinates
 use homogenization, only: &
   materialpoint_dPdF

 implicit none

 DM                                   :: da_local
 Vec                                  :: x_local, coordinates
 Mat                                  :: Jac_pre, Jac
 MatStencil                           :: row(4,24), col(4,24)
 PetscScalar,                 pointer :: x_scal(:,:,:,:)
 PetscScalar                          :: K_ele(24,24), BMatFull(9,24)
 PetscInt                             :: i, ii, j, jj, k, kk, ctr, ele
 PetscInt                             :: nrows, rows(3)
 PetscScalar                          :: diag
 PetscObject                          :: dummy
 MatNullSpace                         :: matnull
 PetscErrorCode                       :: ierr
 
 BMatFull = 0.0
 BMatFull(1:3,1 :8 ) = BMat
 BMatFull(4:6,9 :16) = BMat
 BMatFull(7:9,17:24) = BMat
 call MatSetOption(Jac,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr); CHKERRQ(ierr)
 call MatSetOption(Jac,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr); CHKERRQ(ierr)
 call MatZeroEntries(Jac,ierr); CHKERRQ(ierr)
 ele = 0
 do k = zstart, zend; do j = ystart, yend; do i = xstart, xend
   ctr = 0
   do kk = 0, 1; do jj = 0, 1; do ii = 0, 1
     ctr = ctr + 1
     col(MatStencil_i,ctr   ) = i+ii
     col(MatStencil_j,ctr   ) = j+jj
     col(MatStencil_k,ctr   ) = k+kk
     col(MatStencil_c,ctr   ) = 0
     col(MatStencil_i,ctr+8 ) = i+ii
     col(MatStencil_j,ctr+8 ) = j+jj
     col(MatStencil_k,ctr+8 ) = k+kk
     col(MatStencil_c,ctr+8 ) = 1
     col(MatStencil_i,ctr+16) = i+ii
     col(MatStencil_j,ctr+16) = j+jj
     col(MatStencil_k,ctr+16) = k+kk
     col(MatStencil_c,ctr+16) = 2
   enddo; enddo; enddo
   row = col
   ele = ele + 1
   K_ele = 0.0
   K_ele(1 :8 ,1 :8 ) = HGMat*(materialpoint_dPdF(1,1,1,1,1,ele) + &
                               materialpoint_dPdF(2,2,2,2,1,ele) + &
                               materialpoint_dPdF(3,3,3,3,1,ele))/3.0_pReal
   K_ele(9 :16,9 :16) = HGMat*(materialpoint_dPdF(1,1,1,1,1,ele) + &
                               materialpoint_dPdF(2,2,2,2,1,ele) + &
                               materialpoint_dPdF(3,3,3,3,1,ele))/3.0_pReal
   K_ele(17:24,17:24) = HGMat*(materialpoint_dPdF(1,1,1,1,1,ele) + &
                               materialpoint_dPdF(2,2,2,2,1,ele) + &
                               materialpoint_dPdF(3,3,3,3,1,ele))/3.0_pReal
   K_ele = K_ele + &
           matmul(transpose(BMatFull), &
                  matmul(reshape(reshape(materialpoint_dPdF(1:3,1:3,1:3,1:3,1,ele), &
                                         shape=[3,3,3,3], order=[2,1,4,3]),shape=[9,9]),BMatFull))*detJ
   call MatSetValuesStencil(Jac,24,row,24,col,K_ele,ADD_VALUES,ierr)
   CHKERRQ(ierr)
 enddo; enddo; enddo
 call MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 call MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 call MatAssemblyBegin(Jac_pre,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 call MatAssemblyEnd(Jac_pre,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 
!--------------------------------------------------------------------------------------------------
! applying boundary conditions
 nrows = 3
 rows = [0, 1, 2]
 diag = (C_volAvg(1,1,1,1)/delta(1)/delta(1) + &
         C_volAvg(2,2,2,2)/delta(2)/delta(2) + &
         C_volAvg(3,3,3,3)/delta(3)/delta(3))*detJ
 call MatZeroRowsColumns(Jac,nrows,rows,diag,PETSC_NULL_VEC,PETSC_NULL_VEC,ierr)
 CHKERRQ(ierr)
 call DMGetGlobalVector(da_local,coordinates,ierr);CHKERRQ(ierr)
 call DMDAVecGetArrayF90(da_local,coordinates,x_scal,ierr);CHKERRQ(ierr)
 ele = 0
 do k = zstart, zend; do j = ystart, yend; do i = xstart, xend
   ele = ele + 1
   x_scal(0:2,i,j,k) = mesh_ipCoordinates(1:3,1,ele)
 enddo; enddo; enddo
 call DMDAVecRestoreArrayF90(da_local,coordinates,x_scal,ierr);CHKERRQ(ierr)                       ! initialise to undeformed coordinates (ToDo: use ip coordinates)
 call MatNullSpaceCreateRigidBody(coordinates,matnull,ierr);CHKERRQ(ierr)                           ! get rigid body deformation modes
 call DMRestoreGlobalVector(da_local,coordinates,ierr);CHKERRQ(ierr)
 call MatSetNullSpace(Jac,matnull,ierr); CHKERRQ(ierr)
 call MatSetNearNullSpace(Jac,matnull,ierr); CHKERRQ(ierr)
 call MatNullSpaceDestroy(matnull,ierr); CHKERRQ(ierr)

end subroutine FEM_formJacobian


!--------------------------------------------------------------------------------------------------
!> @brief convergence check
!--------------------------------------------------------------------------------------------------
subroutine FEM_converged(snes_local,PETScIter,xnorm,snorm,fnorm,reason,dummy,ierr)
 use numerics, only: &
   itmax, &
   itmin, &
   err_div_tolRel, &
   err_div_tolAbs, &
   err_stress_tolRel, &
   err_stress_tolAbs, &
   worldrank
 use mesh, only: &
   geomSize
 use spectral_utilities, only: &
   scaledGeomSize, &
   wgt
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
   divTol, &
   stressTol 
 
 err_div = fnorm*sqrt(wgt)*geomSize(1)/scaledGeomSize(1)/detJ
 divTol    = max(maxval(abs(P_av))*err_div_tolRel,err_div_tolAbs)
 stressTol = max(maxval(abs(P_av))*err_stress_tolrel,err_stress_tolabs)
 
 converged: if ((totalIter >= itmin .and. &
                           all([ err_div/divTol, &
                                 err_stress/stressTol       ] < 1.0_pReal)) &
             .or.    terminallyIll) then  
   reason = 1
 elseif (totalIter >= itmax) then converged
   reason = -1
 else converged
   reason = 0
 endif converged

!--------------------------------------------------------------------------------------------------
! report
 if (worldrank == 0_pInt) then
   write(6,'(1/,a)') ' ... reporting .............................................................'
   write(6,'(1/,a,f12.2,a,es8.2,a,es9.2,a)') ' error divergence = ', &
            err_div/divTol,  ' (',err_div,' / m, tol =',divTol,')'
   write(6,'(a,f12.2,a,es8.2,a,es9.2,a)')   ' error stress BC =  ', &
            err_stress/stressTol, ' (',err_stress, ' Pa,  tol =',stressTol,')' 
   write(6,'(/,a)') ' ==========================================================================='
   flush(6) 
 endif
 
end subroutine FEM_converged

!--------------------------------------------------------------------------------------------------
!> @brief forwarding routine
!--------------------------------------------------------------------------------------------------
subroutine FEM_forward(guess,timeinc,timeinc_old,loadCaseTime,deformation_BC,stress_BC,rotation)
 use math, only: &
   math_mul33x33
 use spectral_utilities, only: &
   utilities_updateIPcoords, &
   tTensorBoundaryCondition, &
   cutBack
  use homogenization, only: &
    materialpoint_F0
  use mesh, only: &
    grid, &
    grid3
  use CPFEM2, only: &
    CPFEM_age

 implicit none
 real(pReal), intent(in) :: &
   timeinc_old, &
   timeinc, &
   loadCaseTime                                                                                     !< remaining time of current load case
 type(tTensorBoundaryCondition),      intent(in) :: &
   stress_BC, &
   deformation_BC
 real(pReal), dimension(3,3), intent(in) :: rotation
 logical, intent(in) :: &
   guess
 PetscErrorCode :: ierr

  if (cutBack) then
    C_volAvg    = C_volAvgLastInc                                                                  ! QUESTION: where is this required?
  else
    call CPFEM_age()                                                                                 ! age state and kinematics
    call utilities_updateIPcoords(F_current)

    C_volAvgLastInc    = C_volAvg
 
    F_aimDot = merge(stress_BC%maskFloat*(F_aim-F_aim_lastInc)/timeinc_old, 0.0_pReal, guess)
    F_aim_lastInc = F_aim

    !--------------------------------------------------------------------------------------------------
    ! calculate rate for aim
    if     (deformation_BC%myType=='l') then                                                          ! calculate F_aimDot from given L and current F
      F_aimDot = &
      F_aimDot + deformation_BC%maskFloat * math_mul33x33(deformation_BC%values, F_aim_lastInc)
    elseif(deformation_BC%myType=='fdot') then                                                        ! F_aimDot is prescribed
      F_aimDot = &
      F_aimDot + deformation_BC%maskFloat * deformation_BC%values
    elseif (deformation_BC%myType=='f') then                                                          ! aim at end of load case is prescribed
      F_aimDot = &
      F_aimDot + deformation_BC%maskFloat * (deformation_BC%values - F_aim_lastInc)/loadCaseTime
    endif


    if (guess) then
      call VecWAXPY(solution_rate,-1.0,solution_lastInc,solution_current,ierr)
      CHKERRQ(ierr)
      call VecScale(solution_rate,1.0/timeinc_old,ierr); CHKERRQ(ierr)
    else
      call VecSet(solution_rate,0.0,ierr); CHKERRQ(ierr)
    endif
    call VecCopy(solution_current,solution_lastInc,ierr); CHKERRQ(ierr)
    F_lastInc        = F_current                                                                      ! winding F forward
    materialpoint_F0 = reshape(F_lastInc, [3,3,1,product(grid(1:2))*grid3])                           ! set starting condition for materialpoint_stressAndItsTangent
  endif

!--------------------------------------------------------------------------------------------------
! update average and local deformation gradients
  F_aim = F_aim_lastInc + F_aimDot * timeinc
  call VecAXPY(solution_current,timeinc,solution_rate,ierr); CHKERRQ(ierr)
  
end subroutine FEM_forward

end module spectral_mech_FEM
