!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR0                                                   C
!  Purpose: This routine is called before the time loop starts and is  C
!           user-definable.  The user may insert code in this routine  C
!           or call appropriate user defined subroutines.  This        C
!           can be used for setting constants and checking errors in   C
!           data.  This routine is not called from an IJK loop, hence  C
!           all indices are undefined.                                 C
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
!
      SUBROUTINE USR0

      USE constant
      USE funits
      USE param
      USE param1
      USE physprop
      USE toleranc
      USE usr
      use error_manager

      IMPLICIT NONE

      INCLUDE 'usrnlst.inc'

!-----------------------------------------------

      DOUBLE PRECISION :: SUM

      CALL INIT_ERR_MSG('USR0')

      Allocate(  N_Sh (DIMENSION_3, DIMENSION_M) )

      IF(PAFC .EQ. UNDEFINED ) THEN
         WRITE(ERR_MSG,1100)
          CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
 1100 FORMAT('*** PAFC not specified')

      IF(PAA .NE. UNDEFINED)THEN
         SUM = PAFC + PAA
         IF(.NOT.COMPARE(ONE,SUM) )THEN
            WRITE(ERR_MSG, 1200) SUM
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 1200 FORMAT('*** PAFC + PAA = ',F10.5,/'It should equal 1.0')

       ELSE
          PAA = 1.0 - PAFC
       ENDIF
!
!  Function of the ash-layer void fraction
      f_EP_A = (0.25 + 0.75 * ( 1.0 - PAA )) ** 2.5

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE USR0
