!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Basic PETSc solver
!--------------------------------------------------------------------------------------------------
module spectral_mech_basic
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
#include <petsc/finclude/petsc.h90>

 character (len=*), parameter, public :: &
   DAMASK_spectral_SolverbasicPETSc_label = 'basicpetsc'
   
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
   basicPETSc_init, &
   basicPETSc_solution, &
   basicPETSc_forward, &
   basicPETSc_destroy
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
subroutine basicPETSc_init
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran >4.6 at the moment)
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
   utilities_updateGamma, &
   wgt
 use mesh, only: &
   geomSize, &
   grid, &
   grid3
 use math, only: &
   math_invSym3333
   
 implicit none
 PetscErrorCode       :: ierr
 real(pReal)          :: P_avg(3,3)
 integer(pInt), dimension(:), allocatable :: localK  
 integer(pInt) :: proc
 character(len=1024) :: rankStr
 
 external :: &
   SNESCreate, &
   SNESSetOptionsPrefix, &
   DMDACreate3D, &
   SNESSetDM, &
   DMCreateGlobalVector, &
   DMDASNESSetFunctionLocal, &
   SNESGetConvergedReason, &
   SNESSetConvergenceTest, &
   SNESSetFromOptions
   
 write(6,'(/,a)') ' <<<+-  DAMASK_spectral_solverBasicPETSc init  -+>>>'
 write(6,'(/,a)') ' Shanthraj et al., International Journal of Plasticity, 66:31–45, 2015'
 write(6,'(/,a)') ' https://doi.org/10.1016/j.ijplas.2014.02.006'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

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
        3, 0, &                                                                                     ! #dof (F tensor), ghost boundary width (domain overlap)
        grid (1),grid (2),localK, &                                                                 ! local grid
        mech_grid,ierr)                                                                             ! handle, error
 CHKERRQ(ierr)
 call DMDASetUniformCoordinates(mech_grid,0.0,geomSize(1),0.0,geomSize(2),0.0,geomSize(3),ierr)     ! set dimensions of grid
 CHKERRQ(ierr)
 call SNESSetDM(mech_snes,mech_grid,ierr); CHKERRQ(ierr)
 call DMCreateGlobalVector(mech_grid,solution_current,ierr); CHKERRQ(ierr)                          ! current global displacement vector 
 call DMCreateGlobalVector(mech_grid,solution_lastInc,ierr); CHKERRQ(ierr)                          ! last increment global displacement vector 
 call DMCreateGlobalVector(mech_grid,solution_rate   ,ierr); CHKERRQ(ierr)                          ! current global velocity vector 
 call DMSNESSetFunctionLocal(mech_grid,basicPETSc_formResidual,PETSC_NULL_OBJECT,ierr)              ! residual vector of same shape as solution vector
 CHKERRQ(ierr) 
 call SNESSetConvergenceTest(mech_snes,basicPETSc_converged,PETSC_NULL_OBJECT,PETSC_NULL_FUNCTION,ierr) 
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
 call Utilities_constitutiveResponse(P_current,P_avg,C_volAvg,C_minMaxAvg, &                  ! stress field, stress avg, global average of stiffness and (min+max)/2
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
 call Utilities_updateGamma(C_minMaxAvg,.true.)
   
end subroutine basicPETSc_init
  
!--------------------------------------------------------------------------------------------------
!> @brief solution for the FEM scheme with internal iterations
!--------------------------------------------------------------------------------------------------
type(tSolutionState) function &
  basicPETSc_solution(incInfoIn,timeinc,timeinc_old,loadCaseTime,stress_BC,rotation)
 use IO, only: &
   IO_error
 use numerics, only: &
   update_gamma
 use spectral_utilities, only: &
   tBoundaryCondition, &
   utilities_maskedCompliance, &
   utilities_updateGamma
 use FEsolving, only: &
   restartWrite, &
   terminallyIll

 implicit none

!--------------------------------------------------------------------------------------------------
! input data for solution
 real(pReal), intent(in) :: &
   timeinc, &                                                                                       !< increment in time for current solution
   timeinc_old, &                                                                                   !< increment in time of last increment
   loadCaseTime                                                                                     !< remaining time of current load case
 type(tBoundaryCondition),      intent(in) :: &
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
 if (update_gamma) call utilities_updateGamma(C_minMaxAvg,restartWrite)
 
!--------------------------------------------------------------------------------------------------
! set module wide availabe data 
 mask_stress       = stress_BC%maskFloat
 params%stress_BC  = stress_BC%values
 params%rotation_BC= rotation
 params%timeinc    = timeinc
 params%timeincOld = timeinc_old

!--------------------------------------------------------------------------------------------------
! solve BVP 
 call SNESSolve(mech_snes,PETSC_NULL_OBJECT,solution_current,ierr)
 CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! check convergence
 basicPETSc_solution%iterationsNeeded = totalIter
 basicPETSc_solution%termIll = terminallyIll
 terminallyIll = .false.
 call SNESGetConvergedReason(mech_snes,reason,ierr);CHKERRQ(ierr)
 basicPETSc_solution%converged = reason > 0
 if (reason == SNES_DIVERGED_FNORM_NAN) call IO_error(893_pInt)

end function basicPETSc_solution


!--------------------------------------------------------------------------------------------------
!> @brief forms the residual vector
!--------------------------------------------------------------------------------------------------
subroutine basicPETSc_formResidual(da_local,x_local,f_local,dummy,ierr)
 use numerics, only: &
   itmax, &
   itmin
 use numerics, only: &
   worldrank
 use mesh, only: &
   grid, grid3
 use math, only: &
   math_rotate_backward33, &
   math_transpose33, &
   math_mul3333xx33
 use debug, only: &
   debug_level, &
   debug_spectral, &
   debug_spectralRotation
 use spectral_utilities, only: &
   vectorField_real, &
   tensorField_real, &
   utilities_FFTvectorForward, &
   utilities_FFTvectorBackward, &
   utilities_FFTtensorForward, &
   utilities_FFTtensorBackward, &
   utilities_fourierVectorGradient, &
   utilities_fourierTensorDivergence, &
   utilities_fourierVectorGreenConvolution, &
   utilities_divergenceRMS, &
   utilities_constitutiveResponse
 use IO, only: &
   IO_intOut 
 use FEsolving, only: &
   terminallyIll

 implicit none
 DM                   :: da_local
 Vec                  :: x_local, f_local
 PetscScalar, pointer :: x_scal(:,:,:,:), f_scal(:,:,:,:)
 PetscInt             :: i, ii, j, jj, k, kk, ele
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
 vectorField_real = 0.0_pReal
 tensorField_real = 0.0_pReal
 vectorField_real(1:3,1:grid(1),1:grid(2),1:grid3) = x_scal(0:2,xstart:xend,ystart:yend,zstart:zend)
 call DMDAVecRestoreArrayF90(da_local,x_local,x_scal,ierr);CHKERRQ(ierr)
 call utilities_FFTvectorForward
 call utilities_fourierVectorGradient
 call utilities_FFTtensorBackward
 do k = 1_pInt, grid3; do j = 1_pInt, grid(2); do i = 1_pInt, grid(1)
   F_current(1:3,1:3,i,j,k) = math_rotate_backward33(F_aim,params%rotation_BC) + tensorField_real(1:3,1:3,i,j,k) 
 enddo; enddo; enddo

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
 vectorField_real = 0.0_pReal
 tensorField_real = 0.0_pReal
 tensorField_real(1:3,1:3,1:grid(1),1:grid(2),1:grid3) = P_current
 call utilities_FFTtensorForward
 err_div = utilities_divergenceRMS()
 call utilities_fourierTensorDivergence
 call utilities_fourierVectorGreenConvolution
 call utilities_FFTvectorBackward
 call VecSet(f_local,0.0,ierr);CHKERRQ(ierr)
 call DMDAVecGetArrayF90(da_local,f_local,f_scal,ierr);CHKERRQ(ierr)
 ele = 0
 do k = zstart, zend; do j = ystart, yend; do i = xstart, xend
   ii = i-xstart+1; jj = j-ystart+1; kk = k-zstart+1
   f_scal(0:2,i,j,k) = vectorField_real(1:3,ii,jj,kk)
 enddo; enddo; enddo
 call DMDAVecRestoreArrayF90(da_local,f_local,f_scal,ierr);CHKERRQ(ierr)
 
end subroutine basicPETSc_formResidual


!--------------------------------------------------------------------------------------------------
!> @brief convergence check
!--------------------------------------------------------------------------------------------------
subroutine basicPETSc_converged(snes_local,PETScIter,xnorm,snorm,fnorm,reason,dummy,ierr)
 use numerics, only: &
   itmax, &
   itmin, &
   err_div_tolRel, &
   err_div_tolAbs, &
   err_stress_tolRel, &
   err_stress_tolAbs, &
   worldrank
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
 
end subroutine basicPETSc_converged

!--------------------------------------------------------------------------------------------------
!> @brief forwarding routine
!--------------------------------------------------------------------------------------------------
subroutine basicPETSc_forward(guess,timeinc,timeinc_old,loadCaseTime,deformation_BC,stress_BC,rotation)
 use math, only: &
   math_mul33x33
  use numerics, only: &
    worldrank 
 use spectral_utilities, only: &
   utilities_updateIPcoords, &
   tBoundaryCondition, &
   cutBack
  use IO, only: &
    IO_write_JobRealFile
  use homogenization, only: &
    materialpoint_F0
  use mesh, only: &
    grid, &
    grid3
  use CPFEM2, only: &
    CPFEM_age
  use FEsolving, only: &
    restartWrite

 implicit none
 real(pReal), intent(in) :: &
   timeinc_old, &
   timeinc, &
   loadCaseTime                                                                                     !< remaining time of current load case
 type(tBoundaryCondition),      intent(in) :: &
   stress_BC, &
   deformation_BC
 real(pReal), dimension(3,3), intent(in) :: rotation
 logical, intent(in) :: &
   guess
 PetscErrorCode :: ierr

 character(len=32) :: rankStr

  if (cutBack) then
    C_volAvg    = C_volAvgLastInc                                                                  ! QUESTION: where is this required?
  else
!--------------------------------------------------------------------------------------------------
! restart information for spectral solver
    if (restartWrite) then                                                                           ! QUESTION: where is this logical properly set?
      write(6,'(/,a)') ' writing converged results for restart'
      flush(6)

      if (worldrank == 0_pInt) then
        call IO_write_jobRealFile(777,'C_volAvg',size(C_volAvg))
        write (777,rec=1) C_volAvg; close(777)
        call IO_write_jobRealFile(777,'C_volAvgLastInc',size(C_volAvgLastInc))
        write (777,rec=1) C_volAvgLastInc; close(777)
        call IO_write_jobRealFile(777,'F_aimDot',size(F_aimDot))
        write (777,rec=1) F_aimDot; close(777)
      endif

      write(rankStr,'(a1,i0)')'_',worldrank
      call IO_write_jobRealFile(777,'F'//trim(rankStr),size(F_current))                                      ! writing deformation gradient field to file
      write (777,rec=1) F_current; close (777)
      call IO_write_jobRealFile(777,'F_lastInc'//trim(rankStr),size(F_lastInc))                      ! writing F_lastInc field to file
      write (777,rec=1) F_lastInc; close (777)
    endif

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
  
end subroutine basicPETSc_forward

!--------------------------------------------------------------------------------------------------
!> @brief destroy routine
!--------------------------------------------------------------------------------------------------
subroutine basicPETSc_destroy()

 implicit none
 PetscErrorCode :: ierr

 call VecDestroy(solution_current,ierr); CHKERRQ(ierr)
 call VecDestroy(solution_lastInc,ierr); CHKERRQ(ierr)
 call VecDestroy(solution_rate   ,ierr); CHKERRQ(ierr)
 call SNESDestroy(mech_snes,ierr); CHKERRQ(ierr)
 call DMDestroy(mech_grid,ierr); CHKERRQ(ierr)

end subroutine basicPETSc_destroy

end module spectral_mech_basic
