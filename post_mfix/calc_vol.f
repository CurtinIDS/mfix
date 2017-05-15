!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_VOL                                               C
!  Purpose: Calculate the volume of the cells                          C
!                                                                      C
!  Author: P. Nicoletti                               Date: 01-AUG-92  C
!  Reviewer:                                                           C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: DX,DY,DZ,IMIN1,IMAX1,IMAX2,JMIN1,JMAX1,JMAX2  C
!                        KMIN1,KMAX1,KMAX2,XDIST_SC,XDIST_VEC          C
!                        YDIST_SC,YDIST_VEC,ZDIST_SC,ZDIST_VEC         C
!  Variables modified: I,J,K,VOL,VOL_U,VOL_V,VOL_W                  C
!                                                                      C
!  Local variables: LEN_K,LEN_K_KP,KEN_J,LEN_J_JP,LEN_I,LEN_I_IP       C
!                   XI,XIH                                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
SUBROUTINE CALC_VOL
  !
  !
  USE param
  USE param1
  USE geometry
  USE indices
  USE fldvar
  USE physprop
  USE post3d
  USE compar
  USE functions

  IMPLICIT NONE
  INTEGER  I, J, K, IJK
  INTEGER, EXTERNAL :: FUNIJK_LOC

  DO K = 1,KMAX2
     DO J = 1,JMAX2
        DO I = 1,IMAX2
           IJK = FUNIJK_LOC(I,J,K)
           VOL(IJK)    = DX(I)   * DY(J)   * X(I)*DZ(K)
           VOL_U(IJK)  = DX_E(I) * DY(J)   * X_E(I)*DZ(K)
           VOL_V(IJK)  = DX(I)   * DY_N(J) * X(I)*DZ(K)
           VOL_W(IJK)  = DX(I)   * DY(J)   * X(I)*DZ_T(K)
        END DO
     END DO
  END DO
  !
  RETURN

CONTAINS

  REAL FUNCTION DZ_T(K)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: K
    DZ_T = HALF * (DZ(K) + DZ(KP1(K)))
  END FUNCTION DZ_T

  REAL FUNCTION DY_N(J)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: J
    DY_N = HALF * (DY(J) + DY(JP1(J)))
  END FUNCTION DY_N

  REAL FUNCTION DX_E(I)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: I
    DX_E = HALF * (DX(I) + DX(IP1(I)))
  END FUNCTION DX_E

END SUBROUTINE CALC_VOL
