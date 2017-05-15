! This file is included by fun_avg_mod.f

! These function are are defined in this file so it can be included by
! individual subprograms in order to inline the (short) functions
! defined in this file.

! All functions that may benefits from inlining should be defined in this file.

! For code that is not performance-critical, do:

! SUBROUTINE subname
! USE fun_avg
! ...
! END SUBROUTINE subname

! For code that is performance-critical, do:

! SUBROUTINE subname
! ...
! CONTAINS
! INCLUDE 'fun_avg.inc'
! END SUBROUTINE subname

! Arithmetic average at i+1/2 location (i+1/2, j, k)
! USE IJKE for scalar.
! F(IJK) F(IPJK)  I
      DOUBLE PRECISION FUNCTION AVG_X(XXXm, XXXp, xL)
      USE geometry, only: fx_e, fx_e_bar
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) ::XXXp, XXXm
      INTEGER, INTENT(IN) :: xL
      AVG_X = FX_E(xL) * XXXm + FX_E_bar(xL) * XXXp
      END FUNCTION AVG_X


! average (of U) at i+1 location (i+1, j, k)
! F(IJK) F(IPJK)  IP
      DOUBLE PRECISION FUNCTION AVG_X_E(XXXm, XXXp, xL)
      USE geometry, only: fx, fx_bar
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: XXXp, XXXm
      INTEGER, INTENT(IN) :: xL
      AVG_X_E = FX(xL) * XXXm + FX_bar(xL) * XXXp
      END FUNCTION AVG_X_E


! average at j+1/2 location (i, j+1/2, k)
! Use IJKN for scalar
! F(IJK) F(IJPK) J
      DOUBLE PRECISION FUNCTION AVG_Y(XXXm, XXXp, xL)
      USE geometry, only: fy_n, fy_n_bar
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: XXXp, XXXm
      INTEGER, INTENT(IN) :: xL
      AVG_Y = FY_N(xL) * XXXm + FY_N_bar(xL) * XXXp
      END FUNCTION AVG_Y


! average (of V) at j+1 location (i, j+1, k)
! F(IJK) F(IJPK)
      DOUBLE PRECISION FUNCTION AVG_Y_N(XXXm,  XXXp)
      USE param1, only: half
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: XXXp, XXXm
      AVG_Y_N = HALF *( XXXm + XXXp )
      END FUNCTION AVG_Y_N


! average at k+1/2 location (i, j, k+1/2)
! Use IJKT for scalar
! F(IJK) F(IJKP) K
      DOUBLE PRECISION FUNCTION AVG_Z(XXXm, XXXp, xL)
      USE geometry, only: fz_t, fz_t_bar
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: XXXp, XXXm
      INTEGER, INTENT(IN) :: xL
      AVG_Z = FZ_T(xL) * XXXm + FZ_T_bar(xL) * XXXp
      END FUNCTION AVG_Z


! average (of W) at k+1 location (i, j, k+1)
! F(IJK) F(IJKP)
      DOUBLE PRECISION FUNCTION AVG_Z_T(XXXm,  XXXp)
      USE param1, only: half
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: XXXp, XXXm
      AVG_Z_T = HALF *( XXXm + XXXp )
      END FUNCTION AVG_Z_T

! Harmonic average at i+1/2 location (i+1/2, j, k)
! Use IJKE for scalar
! F(IJK) F(IPJK)  I
      DOUBLE PRECISION FUNCTION AVG_X_h(XXXm, XXXp, xL)
      USE geometry, only: fx_e, fx_e_bar
      USE param1, only: SMALL_NUMBER
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: XXXp, XXXm
      INTEGER, INTENT(IN) :: xL
      AVG_X_h = XXXm * XXXp / &
         MAX(SMALL_NUMBER, (FX_E(xL)*XXXm + FX_E_bar(xL)*XXXp) )
      END FUNCTION AVG_X_h


! Harmonic averageat j+1/2 location (i, j+1/2, k)
! Use IJKN for scalar
! F(IJK) F(IJPK) J
      DOUBLE PRECISION FUNCTION AVG_Y_h(XXXm,  XXXp,   xL)
      USE geometry, only: fy_n, fy_n_bar
      USE param1, only: small_number
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: XXXp, XXXm
      INTEGER, INTENT(IN) :: xL
      AVG_Y_h = XXXm * XXXp / &
         MAX(SMALL_NUMBER, (FY_N(xL)*XXXm + FY_N_bar(xL)*XXXp) )
      END FUNCTION AVG_Y_h


! Harmonic average at k+1/2 location (i, j, k+1/2)
! Use IJKT for scalar
! F(IJK) F(IJKP) K
      DOUBLE PRECISION FUNCTION AVG_Z_h(XXXm,  XXXp,   xL)
      USE geometry, only: fz_t, fz_t_bar
      USE param1, only: small_number
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: XXXp, XXXm
      INTEGER, INTENT(IN) :: xL
      AVG_Z_h = XXXm * XXXp / &
         MAX(SMALL_NUMBER, (FZ_T(xL)*XXXm + FZ_T_bar(xL)*XXXp) )
      END FUNCTION AVG_Z_h


! Harmonic averages for possibly negative scalars (sof, Aug 22 2006)
! Harmonic average at i+1/2 location (i+1/2, j, k)
! Use IJKE for scalar
! F(IJK) F(IPJK)  I
      DOUBLE PRECISION FUNCTION AVG_X_S(XXXm,  XXXp, xL)
      USE geometry, only: fx_e, fx_e_bar
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: XXXp, XXXm
      INTEGER, INTENT(IN) :: xL
      AVG_X_S = XXXm * XXXp / &
         (1D-30 + ( FX_E(xL) * XXXm + FX_E_bar(xL) * XXXp) )
      END FUNCTION AVG_X_S


! For possibly negative scalars:
! Harmonic average at j+1/2 location (i, j+1/2, k)
! Use IJKN for scalar
! F(IJK) F(IJPK) J
      DOUBLE PRECISION FUNCTION AVG_Y_S(XXXm, XXXp, xL)
      USE geometry, only: fy_n, fy_n_bar
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: XXXp, XXXm
      INTEGER, INTENT(IN) :: xL
      AVG_Y_S = XXXm * XXXp / &
         (1D-30 + ( FY_N(xL) * XXXm + FY_N_bar(xL) * XXXp) )
      END FUNCTION AVG_Y_S


! For possibly negative scalars:
! Harmonic average at k+1/2 location (i, j, k+1/2)
! Use IJKT for scalar
! F(IJK) F(IJKP) K
      DOUBLE PRECISION FUNCTION AVG_Z_S(XXXm,  XXXp,   xL)
      USE geometry, only: fz_t, fz_t_bar
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: XXXp, XXXm
      INTEGER, INTENT(IN) :: xL
      AVG_Z_S = XXXm * XXXp / &
         (1D-30 + ( FZ_T(xL) * XXXm + FZ_T_bar(xL) * XXXp) )
      END FUNCTION AVG_Z_S
