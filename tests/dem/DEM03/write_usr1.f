!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_USR1 (L)                                         C
!  Purpose: Write user-defined output                                  C
!                                                                      C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_USR1(L)

      use run, only: TIME

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: L

      SELECT CASE(L)
      CASE(1); CALL WRITE_DES_OUT(TIME)
      END SELECT

      RETURN
      END SUBROUTINE WRITE_USR1


!......................................................................!
!  Subroutine: WRITE_DES_Out                                           !
!                                                                      !
!  Purpose: Calculate the position and velocity (1D Y-axis) of a free  !
!  falling particle. Compare the results to the MFIX-DEM solultion.    !
!                                                                      !
!  Author: J.Musser                                   Date:  Jan-13    !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!......................................................................!
      SUBROUTINE WRITE_DES_Out(lTime)

      Use discretelement
      Use run
      Use usr

      IMPLICIT NONE

! Passed variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION, INTENT(IN) :: lTime


! Local variables
!---------------------------------------------------------------------//
! file name
      CHARACTER*64 :: FNAME1, FNAME2
! logical used for testing is the data file already exists
      LOGICAL :: F_EXISTS1, F_EXISTS2
! file unit for heat transfer data
      INTEGER, PARAMETER :: uPos1 = 2030
      INTEGER, PARAMETER :: uPos2 = 2031

      DOUBLE PRECISION, SAVE :: RK4_TIME = 0.0d0
      DOUBLE PRECISION :: RK4_DT, RK4_DT_LAST
      DOUBLE PRECISION TIME_INTERVAL
      DOUBLE PRECISION, PARAMETER :: RK4_DT_DEFAULT = 1.0d-6
      INTEGER :: I, RK4_STEPS


! Open the files.
      OPEN(UNIT=uPos1,FILE='POST_POS1.dat', &
         POSITION="APPEND",STATUS='OLD')

      OPEN(UNIT=uPos2,FILE='POST_POS2.dat', &
         POSITION="APPEND",STATUS='OLD')

! Calculate the value for the RK4 solutions.
      TIME_INTERVAL = lTime - RK4_TIME
      IF(TIME_INTERVAL .LE. RK4_DT) THEN
         RK4_STEPS = 1
         RK4_DT = TIME_INTERVAL
         RK4_DT_LAST = UNDEFINED
      ELSE
         RK4_STEPS = floor(real(TIME_INTERVAL/RK4_DT_DEFAULT))
         RK4_DT = RK4_DT_DEFAULT
         RK4_DT_LAST = lTime - (RK4_TIME + RK4_STEPS*RK4_DT)
      ENDIF

      DO I=1, RK4_STEPS
         CALL RK4_V4(RK4_DT, gY1, gX1, gY2, gX2)
         RK4_TIME = RK4_TIME + RK4_DT
      ENDDO


      IF(RK4_DT_LAST .NE. UNDEFINED) THEN
         CALL RK4_V4(RK4_DT_LAST, gY1, gX1, gY2, gX2)
         RK4_TIME = RK4_TIME + RK4_DT_LAST
      ENDIF


      if(gY1/=gY1 .OR. gY2/=gY2 .OR. gX1/=gX1 .OR. gX2/=gX2) then
         write(*,*)' NAN found... exiting'
         stop
      endif


! Write the results to a file.
      WRITE(uPos1,"(3x,F15.8,5X,F15.8,2(3x,F15.8))") lTime, gY1,   &
         DES_POS_new(1,2), (ABS(gY1 - DES_POS_new(1,2))/ABS(gY1))*100

      WRITE(uPos2,"(3x,F15.8,5X,F15.8,2(3x,F15.8))") lTime, gY2,   &
         DES_POS_new(2,2), (ABS(gY2 - DES_POS_new(2,2))/ABS(gY1))*100

      CLOSE(uPos1)
      CLOSE(uPos2)


      END SUBROUTINE WRITE_DES_Out
