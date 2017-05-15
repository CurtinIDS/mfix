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
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
      Use usr, only     : x_g_1_ex, x_g_2_ex
      Use usr, only     : lnorms_x_g_1, lnorms_x_g_2
      Use usr, only     : xtr, ytr, ztr
      Use usr, only     : de_x_g_1, de_x_g_2
      Use param, only   : dimension_3

      use ps, only: PS_MASSFLOW_G

      IMPLICIT NONE
!-----------------------------------------------
!
!  Include files defining common blocks here
!
!
!  Define local variables here
!
!
!  Include files defining statement functions here
!
!
!  Insert user-defined code here
!

      write(*,*) 'Forcing a negative mass flow'
      PS_MASSFLOW_G(2) = -PS_MASSFLOW_G(1)


! allocate variables defined in usr_mod.f
        allocate(x_g_1_ex(dimension_3))
        allocate(x_g_2_ex(dimension_3))

        allocate(lnorms_x_g_1(3))
        allocate(lnorms_x_g_2(3))

        allocate(xtr(dimension_3))
        allocate(ytr(dimension_3))
        allocate(ztr(dimension_3))

        allocate(de_x_g_1(dimension_3))
        allocate(de_x_g_2(dimension_3))


      RETURN
      END SUBROUTINE USR0
