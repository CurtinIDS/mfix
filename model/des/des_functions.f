
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Function: DES_GETINDEXFROMPOS
!  Purpose: Knowing the current particle x, y, or z position determine
!  the associated i, j or k index by searching a specific region defined
!  by a lower and upper limit on the index that is known to contain
!  the particles position
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      INTEGER FUNCTION DES_GETINDEXFROMPOS(LIM1,LIM2,PART_POS,&
         GRID_POS,AXIS,AXIS_INDEX)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE funits
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! given i, j, or k index values defining the upper and lower limits of
! the region in space known to contain the particle
      INTEGER, INTENT (IN) :: LIM1, LIM2
! given current x, y, or z position of the particle
      DOUBLE PRECISION, INTENT (IN) :: PART_POS
! given position for the cell faces on the fluid grid in the x, y or z
! direction
      DOUBLE PRECISION, DIMENSION(:), INTENT (IN) :: GRID_POS
! given axis (x, y, or z) and associated index (i, j, or k)
      CHARACTER(LEN=1), INTENT (IN) :: AXIS, AXIS_INDEX
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! index & loop counter
      INTEGER IND
!-----------------------------------------------

! error condition
      IND = -1

      IF (LIM1 <= LIM2) THEN
         DO IND = LIM1, LIM2
            IF (PART_POS >= GRID_POS(IND-1) .AND. &
                PART_POS <  GRID_POS(IND)) EXIT
         ENDDO
      ELSEIF (LIM1 > LIM2) THEN
         DO IND = LIM1, LIM2, -1
            IF (PART_POS >= GRID_POS(IND-1) .AND. &
                PART_POS <  GRID_POS(IND)) EXIT
         ENDDO
      ENDIF

      IF (IND == -1) THEN
         WRITE (UNIT_LOG, 1001) AXIS_INDEX, AXIS, PART_POS, &
            AXIS_INDEX, LIM1, AXIS_INDEX, LIM2
         WRITE (*,1001) AXIS_INDEX, AXIS, PART_POS, &
            AXIS_INDEX, LIM1, AXIS_INDEX, LIM2
!         CALL MFIX_EXIT(myPE)
      ENDIF

      DES_GETINDEXFROMPOS = IND

      RETURN

 1001 FORMAT(/1X,70('*')//,' From: DES_GETINDEXFROMPOS',/,' Message: ',&
         'Could not identify the ', A, ' index associated with the ',&
         'particles ',/10X,A, '-position= ',ES15.5,' between the ',&
         'limits ',A,'=',I5,/10X,' and ',A,'=',I5,'.',/1X,70('*')/)

      END FUNCTION DES_GETINDEXFROMPOS




