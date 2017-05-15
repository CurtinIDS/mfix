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
      SUBROUTINE USR0

      USE compar
      USE constant
      USE exit, only: mfix_exit
      USE funits
      USE param
      USE param1
      USE physprop
      USE toleranc
      USE usr

      IMPLICIT NONE

      INCLUDE 'usrnlst.inc'

! Allocate the variable for Sherwood number
      Allocate( N_Sh (DIMENSION_3, DIMENSION_M) )

! Verify that PAFC and PAA were specified in the dat file and that
! the values are physical.
      IF(C(1) .EQ. UNDEFINED ) THEN
         IF(DMP_LOG) WRITE(*,1100)
         CALL MFIX_EXIT(myPE)
      ELSEIF(C(1) < 0.0d0  .OR. C(1) > 1.0d0) THEN
         IF(DMP_LOG) WRITE(*,1101)
         CALL MFIX_EXIT(myPE)
      ELSE
         PAFC = C(1)
      ENDIF

 1100 FORMAT(2/,1x,70('*'),'From: USR0',/' Error 1100: PAFC not ',     &
         'specified in C(1).'/1x,70('*'))

 1101 FORMAT(2/,1x,70('*'),'From: USR0',/' Error 1101: PAFC specified',&
         ' in C(1) is unphysical.'/1x,70('*'))

      IF(C(2) .EQ. UNDEFINED) THEN
         PAA = 1.0d0 - PAFC
      ELSE
         PAA = C(2)
      ENDIF

      IF(.NOT.COMPARE(ONE,(PAFC+PAA)) )THEN
         IF(DMP_LOG) WRITE(*,1102)
         CALL MFIX_EXIT(myPE)
      ENDIF

 1102 FORMAT(2/,1x,70('*'),'From: USR0',/' Error 1102: Sum of PAFC ',  &
         'and PAA does not equal 1.0d0.'/1x,70('*'))

!  Function of the ash-layer void fraction
      f_EP_A = (0.25 + 0.75 * ( 1.0 - PAA )) ** 2.5

      RETURN
      END SUBROUTINE USR0
