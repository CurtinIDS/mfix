!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_DPoDY (TAVG,NTAVG,NPOINTS)                        C
!  Purpose: Update the time-averaged DPoDY calculation                 C
!           Should FLAG be referenced ?????                            C
!                                                                      C
!  Author: P. Nicoletti, M. Syamlal                   Date: 16-FEB-92  C
!  Reviewer:                                                           C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IMIN1, IMAX1, JMIN1, JMAX1, KMIN1, KMAX1      C
!                        P_g, DY                                       C
!  Variables modified: I,J,K,IJK,IJPK                                  C
!                                                                      C
!  Local variables: DPoDY, NCOUNT                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_DPoDY(TAVG,NTAVG,NPOINTS)
!
!
      Use param
      Use param1
      Use fldvar
      Use indices
      Use physprop
      Use geometry
      Use compar
      Use functions

      IMPLICIT NONE
!
!     Passed arguments
!
      DOUBLE PRECISION TAVG(DIMENSION_3,*)
      INTEGER          NPOINTS , NTAVG
      INTEGER        I, J, K, IJK, IJPK
!
!     local variables
!
      INTEGER          NCOUNT
      DOUBLE PRECISION DPoDY
!
      DO J = JMIN1,JMAX1
         DPoDY = 0.0
         NCOUNT = 0
         DO K = KMIN1,KMAX1
            DO I = IMIN1,IMAX1
               IJK = FUNIJK(I,J,K)
               IJPK = FUNIJK(I,J+1,K)
               NCOUNT = NCOUNT + 1
               DPoDY = DPoDY &
                 -( P_g(IJPK) - P_g(IJK) ) * 2.0 / (DY(J+1)+DY(J))
            END DO
         END DO
         TAVG(J,NTAVG) = TAVG(J,NTAVG) + DPoDY / REAL(NCOUNT)
      END DO
      NPOINTS = NPOINTS + 1
!
      RETURN
      END
