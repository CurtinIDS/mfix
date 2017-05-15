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

      INTEGER :: NP

      CHARACTER*11, PARAMETER :: FNAME = 'POST_L1.dat'
      INTEGER, PARAMETER :: FUNIT = 2030
      LOGICAL :: FEXISTS

      INTEGER :: I, MAX_I
      INTEGER :: J, MAX_J

      DOUBLE PRECISION :: X_ERR, Y_ERR, Z_ERR
      DOUBLE PRECISION :: lTHETA, lPHI

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

      NP = 0

      MAX_I = 32
      DO I=1, MAX_I

         lTHETA = dble(I-1) * PI / dble(MAX_I -1)

         MAX_J = merge(2*I-1, 2*(MAX_I-I)+1, I < MAX_I/2+1)
         DO J=1, MAX_J

            NP = NP + 1

            lPHI = J*2.0d0*PI/(MAX_J)

            X_ERR = (0.15d0*sin(lTHETA)*cos(lPHI) + 0.35d0) - DES_POS_NEW(NP,1)
            Y_ERR = (0.15d0*sin(lTHETA)*sin(lPHI) + 0.35d0) - DES_POS_NEW(NP,2)
            Z_ERR = (0.15d0*cos(lTHETA) + 0.35d0)  - DES_POS_NEW(NP,3)

            L1 = sqrt(X_ERR**2 + Y_ERR**2 + Z_ERR**2)
            SUM_L1 = SUM_L1 + L1
            MAX_L1 = max(L1, MAX_L1)
         ENDDO
      ENDDO

      write(FUNIT,"(2(3X,f7.4),2(3x,g11.5))") &
         TIME, TIME/USR_DT(L), MAX_L1, SUM_L1

      CLOSE(FUNIT)

      RETURN
      END SUBROUTINE WRITE_USR1
