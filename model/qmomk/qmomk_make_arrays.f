!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: QMOMK_MAKE_ARRAYS                                      C
!  Purpose: DES - Initialize/Restart QMOMK arrays                      C
!                                                                      C
!                                                                      C
!  Author: Alberto Passalacqua                        Date:            C
!  Reviewer:                                          Date:            C
!  Comments:                                                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

SUBROUTINE qmomk_make_arrays

  USE compar
  USE exit, only: mfix_exit
  USE funits
  USE geometry
  USE param1
  USE qmom_kinetic_equation
  USE run

  IMPLICIT NONE

  INTEGER CHECK_MPI

  IF(COORDINATES == 'CYLINDRICAL') THEN
     WRITE (UNIT_LOG, *) ' '
     WRITE (UNIT_LOG, *) 'Cylindrical coordinates are being used. STOP'
     WRITE (UNIT_LOG, *) 'QMOMK should only be run using cartesian coordinates.'
     WRITE (*, *) ' '
     WRITE (*, *) 'Cylindrical coordinates are being used. STOP'
     WRITE (*, *) 'QMOMK should only be run using cartesian coordinates.'
     CALL MFIX_EXIT(myPE)
  END IF

  CHECK_MPI = NODESI * NODESJ * NODESK
  IF((CHECK_MPI.NE.1).AND.(QMOMK)) THEN
     WRITE (UNIT_LOG, *) ' '
     WRITE (UNIT_LOG, *) 'QMOMK being run on multiple processors. STOP'
     WRITE (UNIT_LOG, *) 'QMOMK should only be run serially on one processor.'
     WRITE (*, *) ' '
     WRITE (*, *) 'QMOMK being run on multiple processors. STOP'
     WRITE (*, *) 'QMOMK should only be run serially on one processor.'
     CALL MFIX_EXIT(myPE)
  END IF

  IF(RUN_TYPE == 'RESTART_1') THEN !  Read Restart
     CALL QMOMK_READ_RESTART
     WRITE(*,*) 'QMOMK_RES file read at Time= ', TIME
     WRITE(UNIT_LOG,*) 'QMOMK_RES file read at Time= ', TIME
  ELSE IF (RUN_TYPE == 'RESTART_2') THEN
     WRITE(UNIT_LOG,*) 'Restart 2 is not implemented with QMOMK'
     WRITE(*,*) 'Restart 2 is not implemented with QMOMK'
     CALL MFIX_EXIT(myPE)
  END IF

  RETURN
END SUBROUTINE qmomk_make_arrays
