!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief quit subroutine
!> @details exits the program and reports current time and duration. Exit code 0 signals
!> everything is fine. Exit code 1 signals an error, message according to IO_error. Exit code
!> 2 signals no severe problems, but some increments did not converge
!--------------------------------------------------------------------------------------------------
subroutine quit(stop_id)
#include <petsc/finclude/petscsys.h>
#ifdef _OPENMP
 use MPI, only: &
   MPI_finalize
#endif
 use PetscSys
 use hdf5

 implicit none
 integer, intent(in) :: stop_id
 integer, dimension(8) :: dateAndTime                                                               ! type default integer
 integer :: error
 PetscErrorCode :: ierr = 0

 call h5open_f(error)
 if (error /= 0) write(6,'(a,i5)') ' Error in h5open_f ',error                                      ! prevents error if not opened yet
 call h5close_f(error)
 if (error /= 0) write(6,'(a,i5)') ' Error in h5close_f ',error

 call PETScFinalize(ierr)
 CHKERRQ(ierr)

#ifdef _OPENMP
 call MPI_finalize(error)
 if (error /= 0) write(6,'(a,i5)') ' Error in MPI_finalize',error
#endif
 
 call date_and_time(values = dateAndTime)
 write(6,'(/,a)') 'DAMASK terminated on:'
 write(6,'(a,2(i2.2,a),i4.4)') 'Date:               ',dateAndTime(3),'/',&
                                                      dateAndTime(2),'/',&
                                                      dateAndTime(1)
 write(6,'(a,2(i2.2,a),i2.2)') 'Time:               ',dateAndTime(5),':',&
                                                      dateAndTime(6),':',&
                                                      dateAndTime(7)

 if (stop_id == 0 .and. ierr == 0 .and. error == 0) stop 0                                          ! normal termination
 if (stop_id == 2 .and. ierr == 0 .and. error == 0) stop 2                                          ! not all incs converged
 stop 1                                                                                             ! error (message from IO_error)

end subroutine quit
