!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: WRITE_USR1 (L)                                         !
!  Purpose: Write user-defined output                                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE WRITE_USR1(L)

      use constant
      use discretelement
      use output
      use run

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: L

! Loop index
      INTEGER :: NP

      CHARACTER*11, PARAMETER :: FNAME = 'POST_L1.dat'
      INTEGER, PARAMETER :: FUNIT = 2030
      LOGICAL :: FEXISTS

      DOUBLE PRECISION :: X_ERR, Y_ERR, lPHI
      DOUBLE PRECISION :: L1, MAX_L1, SUM_L1

! This file is only designed to output if L=1
      IF(L/=1) RETURN

      INQUIRE(FILE=FNAME, EXIST=FEXISTS)
      IF(.NOT.FEXISTS) THEN
         OPEN(UNIT=FUNIT, FILE=FNAME, STATUS='NEW')
         WRITE(FUNIT,"(5x,A,5X,A,4x,A,3x,A)")'TIME','CYCLE', &
            'MAX L1-NORM','SUM L1-NORM'
      ELSE
         OPEN(UNIT=FUNIT,FILE=FNAME,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF

      MAX_L1 = 0.0d0
      SUM_L1 = 0.0d0
      DO NP = 1, 64

         lPHI = dble(NP-1) * 2.0d0*PI/64
         X_ERR = (0.15d0*cos(lPHI) + 0.50d0) - DES_POS_NEW(NP,1)
         Y_ERR = (0.15d0*sin(lPHI) + 0.75d0) - DES_POS_NEW(NP,2)

         L1 = sqrt(X_ERR**2 + Y_ERR**2)
         SUM_L1 = SUM_L1 + L1
         MAX_L1 = max(L1, MAX_L1)
      END DO

      write(FUNIT,"(2(3X,f7.4),2(3x,g11.5))") &
         TIME, TIME/USR_DT(L), MAX_L1, SUM_L1

      CLOSE(FUNIT)

      RETURN
      END SUBROUTINE WRITE_USR1
