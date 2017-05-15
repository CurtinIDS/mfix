!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_LOCATION                                           C
!  Purpose: Process the LOCATION sub-menu for post-processing          C
!                                                                      C
!  Author: P. Nicoletti                               Date: 13-MAR-92  C
!  Reviewer:                                                           C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IMIN1,IMAX2, JMIN1, JMAX2, KMIN1, KMAX2       C
!  Variables modified: LOC_X, LOC_Y, LOC_Z                             C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_LOCATION
!
!
      Use param
      Use param1
      Use geometry
      Use post3d
      IMPLICIT NONE
!
      WRITE (*,1000) IMIN1,IMAX2,JMIN1,JMAX2,KMIN1,KMAX2
      WRITE (*,*) ' ENTER LOC_X,LOC_Y,LOC_Z'
      READ  (*,*) LOC_X,LOC_Y,LOC_Z
      RETURN
1000  FORMAT(1X,'IMIN1 = ' , I5,10X,'IMAX2 = ', I5,/,&
             1X,'JMIN1 = ' , I5,10X,'JMAX2 = ', I5,/,&
             1X,'KMIN1 = ' , I5,10X,'KMAX2 = ', I5,/,/,/)
      END
