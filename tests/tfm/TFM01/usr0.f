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
      Use usr, only     : p_g_ex, u_g_ex, v_g_ex
      Use usr, only     : lnorms_p_g, lnorms_u_g, lnorms_v_g
      Use usr, only     : xtr, ytr, ztr
      Use usr, only     : de_p_g, de_u_g, de_v_g
      Use param, only   : dimension_3
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

! allocate variables defined in usr_mod.f
        allocate(p_g_ex(dimension_3))
        allocate(u_g_ex(dimension_3))
        allocate(v_g_ex(dimension_3))

        allocate(lnorms_p_g(3))
        allocate(lnorms_u_g(3))
        allocate(lnorms_v_g(3))

        allocate(xtr(dimension_3))
        allocate(ytr(dimension_3))
        allocate(ztr(dimension_3))

        allocate(de_p_g(dimension_3))
        allocate(de_u_g(dimension_3))
        allocate(de_v_g(dimension_3))


      RETURN
      END SUBROUTINE USR0
