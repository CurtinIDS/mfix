!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_VTC_SS(VxTC_ss, IER)                              C
!  Purpose: multiply solid-solid granular energy transfer coefficient  C
!           with cell volume
!                                                                      C
!  Author: J. Galvin                                 Date:             C
!  Reviewer: S. benyahia                             Date:             C
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
      SUBROUTINE CALC_VTC_SS(VXTC_SS)
!
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!     Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE indices
      USE physprop
      USE compar
      USE kintheory
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!     Local variables
!-----------------------------------------------

!                      Indices
      INTEGER          IJK
!
!                      Index of both solids phases
      INTEGER          L, M, LM
!
!                      Volume x interphase transfer coefficient
      DOUBLE PRECISION VxTC_ss(DIMENSION_3, DIMENSION_LM)
!-----------------------------------------------
!
      DO M = 1, MMAX
          DO L = 1, MMAX
               LM = FUNLM(L,M)
               IF (L .NE. M) THEN
!!!$omp  parallel do private(IJK)
                    DO IJK = ijkstart3, ijkend3
                         IF (FLUID_AT(IJK)) THEN
                              VXTC_SS(IJK,LM) = ED_ss_ip(IJK,LM)*VOL(IJK)
                         ELSE
                              VXTC_SS(IJK,LM) = ZERO
                         ENDIF
                    ENDDO
               ENDIF
          ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE CALC_VTC_SS
!-----------------------------------------------

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
