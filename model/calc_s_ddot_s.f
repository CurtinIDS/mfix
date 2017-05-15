!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_S_DDOT_S(IJK1, IJK2, FCELL, COM, M, DEL_DOT_U,    C
!                             S_DDOT_S, S_dd)                          C
!                                                                      C
!  Purpose: Calculate del.u, S:S and S_xx, S_yy or S_zz at the         C
!           boundary for use in frictional boundary condition          C
!                                                                      C
!                                                                      C
!  Author: Anuj Srivastava, Princeton University      Date: 4-APR-98   C
!  Reviewer:                                          Date:            C
!                                                                      C
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
      SUBROUTINE CALC_S_DDOT_S(IJK1,IJK2,FCELL,COM,M,DEL_DOT_U,S_DDOT_S,S_DD)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE constant
      USE fldvar
      USE geometry
      USE indices
      USE compar
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      IJK indices for wall cell and fluid cell
      INTEGER          IJK1, IJK2
!
!                      Other indices
      INTEGER          IPJK2, IPJMK2, IPJPK2, IPJKM2, IPJKP2
      INTEGER          IJPK2, IJPKM2, IJPKP2, IMJPK2
      INTEGER          IJKP2, IJMKP2, IMJKP2

!
!                      The location (e,w,n...) of fluid cell
      CHARACTER        FCELL
!
!                      Velocity component (U, V, W)
      CHARACTER        COM
!
!                      Solids phase index
      INTEGER          M
!
!                      del.u
      DOUBLE PRECISION DEL_DOT_U
!
!                      S:S
      DOUBLE PRECISION S_DDOT_S

!                      S_dd (dd is the relevant direction x,y or z)
      DOUBLE PRECISION S_dd
!
!  The location where D (rate of strain tensor) is calculated in
!  located at the center of 4
!  cells - 2 fluid cells and 2 wall cells. Coordinates of this are i,j,k

!                      U_s at the north (i, j+1/2, k)
      DOUBLE PRECISION U_s_N
!
!                      U_s at the south (i, j-1/2, k)
      DOUBLE PRECISION U_s_S

!                      U_s at the east (i+1/2, j, k)
      DOUBLE PRECISION U_s_E
!
!                      U_s at the west (i-1/2, j, k)
      DOUBLE PRECISION U_s_W
!
!                      U_s at the top (i, j, k+1/2)
      DOUBLE PRECISION U_s_T
!
!                      U_s at the bottom (i, j, k-1/2)
      DOUBLE PRECISION U_s_B
!
!                      U_s at the center (i, j, k)
!                      Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION U_s_C
!
!                      V_s at the north (i, j+1/2, k)
      DOUBLE PRECISION V_s_N
!
!                      V_s at the south (i, j-1/2, k)
      DOUBLE PRECISION V_s_S
!
!                      V_s at the east (i+1/2, j, k)
      DOUBLE PRECISION V_s_E
!
!                      V_s at the west (i-1/2, j, k)
      DOUBLE PRECISION V_s_W

!                      V_s at the top (i, j, k+1/2)
      DOUBLE PRECISION V_s_T
!
!                      V_s at the bottom (i, j, k-1/2)
      DOUBLE PRECISION V_s_B

!                      W_s at the north (i, j+1/2, k)
      DOUBLE PRECISION W_s_N
!
!                      W_s at the south (i, j-1/2, k)
      DOUBLE PRECISION W_s_S
!
!                      W_s at the east (i+1/2, j, k)
      DOUBLE PRECISION W_s_E
!
!                      W_s at the west (1-1/2, j, k)
      DOUBLE PRECISION W_s_W
!
!                      W_s at the top (i, j, k+1/2)
      DOUBLE PRECISION W_s_T
!
!                      W_s at the bottom (i, j, k-1/2)
      DOUBLE PRECISION W_s_B
!
!                      W_s at the center (i, j, k).
!                      Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION W_s_C
!-----------------------------------------------

      SELECT CASE (TRIM(COM))
      CASE ('U')
         SELECT CASE (TRIM(FCELL))
         CASE ('N')
            IPJK2 = IP_OF(IJK2)
            IPJMK2 = JM_OF(IPJK2)
!
            U_S_N = U_S(IJK2,M)
!
            U_S_S = U_S(IJK1,M)
!
            U_S_E = AVG_Y(AVG_X_E(U_S(IJK1,M),U_S(IPJMK2,M),I_OF(IPJMK2)),&
               AVG_X_E(U_S(IJK2,M),U_S(IPJK2,M),I_OF(IPJK2)),J_OF(IPJMK2))
!
            U_S_W = AVG_Y(AVG_X_E(U_S(IM_OF(IJK1),M),U_S(IJK1,M),I_OF(IJK1)),&
               AVG_X_E(U_S(IM_OF(IJK2),M),U_S(IJK2,M),I_OF(IJK2)),J_OF(IJK1))
!
            U_S_T = AVG_Y(AVG_Z(U_S(IJK1,M),U_S(KP_OF(IJK1),M),K_OF(IJK1)),&
               AVG_Z(U_S(IJK2,M),U_S(KP_OF(IJK2),M),K_OF(IJK2)),J_OF(IJK1))
!
            U_S_B = AVG_Y(AVG_Z(U_S(KM_OF(IJK1),M),U_S(IJK1,M),K_OF(KM_OF(&
               IJK1))),AVG_Z(U_S(KM_OF(IJK2),M),U_S(IJK2,M),K_OF(KM_OF(IJK2))&
               ),J_OF(IJK1))
!
            V_S_N = AVG_X(AVG_Y_N(V_S(IJK1,M),V_S(IJK2,M)),AVG_Y_N(V_S(&
               IPJMK2,M),V_S(IPJK2,M)),I_OF(IJK2))
!
            V_S_S = AVG_X(AVG_Y_N(V_S(JM_OF(IJK1),M),V_S(IJK1,M)),AVG_Y_N(&
               V_S(JM_OF(IPJMK2),M),V_S(IPJMK2,M)),I_OF(IJK1))
!
            V_S_E = ZERO
!
            V_S_W = ZERO
!
            V_S_T = ZERO
!
            V_S_B = ZERO
!
            W_S_N = AVG_X(AVG_Z_T(W_S(KM_OF(IJK2),M),W_S(IJK2,M)),AVG_Z_T(&
               W_S(KM_OF(IPJK2),M),W_S(IPJK2,M)),I_OF(IJK2))
!
            W_S_S = AVG_X(AVG_Z_T(W_S(KM_OF(IJK1),M),W_S(IJK1,M)),AVG_Z_T(&
               W_S(KM_OF(IPJMK2),M),W_S(IPJMK2,M)),I_OF(IJK1))
!
            W_S_E = AVG_Y(AVG_Z_T(W_S(KM_OF(IPJMK2),M),W_S(IPJMK2,M)),AVG_Z_T&
               (W_S(KM_OF(IPJK2),M),W_S(IPJK2,M)),J_OF(IPJMK2))
!
            W_S_W = AVG_Y(AVG_Z_T(W_S(KM_OF(IJK1),M),W_S(IJK1,M)),AVG_Z_T(&
               W_S(KM_OF(IJK2),M),W_S(IJK2,M)),J_OF(IJK1))
!
            W_S_T = AVG_X(AVG_Y(W_S(IJK1,M),W_S(IJK2,M),J_OF(IJK1)),AVG_Y(&
               W_S(IPJMK2,M),W_S(IPJK2,M),J_OF(IPJMK2)),I_OF(IJK1))
!
            W_S_B = AVG_X(AVG_Y(W_S(KM_OF(IJK1),M),W_S(KM_OF(IJK2),M),J_OF(&
               KM_OF(IJK1))),AVG_Y(W_S(KM_OF(IPJMK2),M),W_S(KM_OF(IPJK2),M),&
               J_OF(KM_OF(IPJMK2))),I_OF(IJK1))
!
            IF (CYLINDRICAL) THEN
               U_S_C = AVG_Y(U_S(IJK1,M),U_S(IJK2,M),J_OF(IJK1))
               W_S_C = AVG_X(W_S_W,W_S_E,I_OF(IJK1))
            ELSE
               U_S_C = ZERO
               W_S_C = ZERO
            ENDIF
!
            CALL SDDOTS (IJK1, FCELL, 'Y', U_S_N, U_S_S, U_S_E, U_S_W, U_S_T, &
               U_S_B, U_S_C, V_S_N, V_S_S, V_S_E, V_S_W, V_S_T, V_S_B, W_S_N, &
               W_S_S, W_S_E, W_S_W, W_S_T, W_S_B, W_S_C, DEL_DOT_U, S_DDOT_S, &
               S_DD)
!
!
!
!
         CASE ('S')
            IPJK2 = IP_OF(IJK2)
            IPJPK2 = JP_OF(IPJK2)
!
            U_S_N = U_S(IJK1,M)
!
            U_S_S = U_S(IJK2,M)
!
            U_S_E = AVG_Y(AVG_X_E(U_S(IJK2,M),U_S(IPJK2,M),I_OF(IPJK2)),&
               AVG_X_E(U_S(IJK1,M),U_S(IPJPK2,M),I_OF(IPJPK2)),J_OF(IPJK2))
!
            U_S_W = AVG_Y(AVG_X_E(U_S(IM_OF(IJK2),M),U_S(IJK2,M),I_OF(IJK2)),&
               AVG_X_E(U_S(IM_OF(IJK1),M),U_S(IJK1,M),I_OF(IJK1)),J_OF(IJK2))
!
            U_S_T = AVG_Y(AVG_Z(U_S(IJK2,M),U_S(KP_OF(IJK2),M),K_OF(IJK2)),&
               AVG_Z(U_S(IJK1,M),U_S(KP_OF(IJK1),M),K_OF(IJK1)),J_OF(IJK2))
!
            U_S_B = AVG_Y(AVG_Z(U_S(KM_OF(IJK2),M),U_S(IJK2,M),K_OF(KM_OF(&
               IJK2))),AVG_Z(U_S(KM_OF(IJK1),M),U_S(IJK1,M),K_OF(KM_OF(IJK1))&
               ),J_OF(IJK2))
!
            V_S_N = AVG_X(AVG_Y_N(V_S(IJK2,M),V_S(IJK1,M)),AVG_Y_N(V_S(IPJK2&
               ,M),V_S(IPJPK2,M)),I_OF(IJK1))
!
            V_S_S = AVG_X(AVG_Y_N(V_S(JM_OF(IJK2),M),V_S(IJK2,M)),AVG_Y_N(&
               V_S(JM_OF(IPJK2),M),V_S(IPJK2,M)),I_OF(IJK2))
!
            V_S_E = ZERO
!
            V_S_W = ZERO
!
            V_S_T = ZERO
!
            V_S_B = ZERO
!
            W_S_N = AVG_X(AVG_Z_T(W_S(KM_OF(IJK1),M),W_S(IJK1,M)),AVG_Z_T(&
               W_S(KM_OF(IPJPK2),M),W_S(IPJPK2,M)),I_OF(IJK1))
!
            W_S_S = AVG_X(AVG_Z_T(W_S(KM_OF(IJK2),M),W_S(IJK2,M)),AVG_Z_T(&
               W_S(KM_OF(IPJK2),M),W_S(IPJK2,M)),I_OF(IJK2))
!
            W_S_E = AVG_Y(AVG_Z_T(W_S(KM_OF(IPJK2),M),W_S(IPJK2,M)),AVG_Z_T(&
               W_S(KM_OF(IPJPK2),M),W_S(IPJPK2,M)),J_OF(IPJK2))
!
            W_S_W = AVG_Y(AVG_Z_T(W_S(KM_OF(IJK2),M),W_S(IJK2,M)),AVG_Z_T(&
               W_S(KM_OF(IJK1),M),W_S(IJK1,M)),J_OF(IJK2))
!
            W_S_T = AVG_X(AVG_Y(W_S(IJK2,M),W_S(IJK1,M),J_OF(IJK2)),AVG_Y(&
               W_S(IPJK2,M),W_S(IPJPK2,M),J_OF(IPJK2)),I_OF(IJK2))
!
            W_S_B = AVG_X(AVG_Y(W_S(KM_OF(IJK2),M),W_S(KM_OF(IJK1),M),J_OF(&
               KM_OF(IJK2))),AVG_Y(W_S(KM_OF(IPJK2),M),W_S(KM_OF(IPJPK2),M),&
               J_OF(KM_OF(IPJK2))),I_OF(IJK2))
!
            IF (CYLINDRICAL) THEN
               U_S_C = AVG_Y(U_S(IJK2,M),U_S(IJK1,M),J_OF(IJK2))
               W_S_C = AVG_X(W_S_W,W_S_E,I_OF(IJK2))
            ELSE
               U_S_C = ZERO
               W_S_C = ZERO
            ENDIF
!
            CALL SDDOTS (IJK2, FCELL, 'Y', U_S_N, U_S_S, U_S_E, U_S_W, U_S_T, &
               U_S_B, U_S_C, V_S_N, V_S_S, V_S_E, V_S_W, V_S_T, V_S_B, W_S_N, &
               W_S_S, W_S_E, W_S_W, W_S_T, W_S_B, W_S_C, DEL_DOT_U, S_DDOT_S, &
               S_DD)
!
!
         CASE ('T')
            IPJK2 = IP_OF(IJK2)
            IPJKM2 = KM_OF(IPJK2)
!
!
            U_S_N = AVG_Z(AVG_Y(U_S(IJK1,M),U_S(JP_OF(IJK1),M),J_OF(IJK1)),&
               AVG_Y(U_S(IJK2,M),U_S(JP_OF(IJK2),M),J_OF(IJK2)),K_OF(IJK1))
!
            U_S_S = AVG_Z(AVG_Y(U_S(JM_OF(IJK1),M),U_S(IJK1,M),J_OF(JM_OF(&
               IJK1))),AVG_Y(U_S(JM_OF(IJK2),M),U_S(IJK2,M),J_OF(JM_OF(IJK2))&
               ),K_OF(IJK1))
!
            U_S_E = AVG_Z(AVG_X_E(U_S(IJK1,M),U_S(IPJKM2,M),I_OF(IPJKM2)),&
               AVG_X_E(U_S(IJK2,M),U_S(IPJK2,M),I_OF(IPJK2)),K_OF(IPJKM2))
!
            U_S_W = AVG_Z(AVG_X_E(U_S(IM_OF(IJK1),M),U_S(IJK1,M),I_OF(IJK1)),&
               AVG_X_E(U_S(IM_OF(IJK2),M),U_S(IJK2,M),I_OF(IJK2)),K_OF(IJK1))
!
            U_S_T = U_S(IJK2,M)
!
            U_S_B = U_S(IJK1,M)
!
            V_S_N = AVG_X(AVG_Z(V_S(IJK1,M),V_S(IJK2,M),K_OF(IJK1)),AVG_Z(&
               V_S(IPJKM2,M),V_S(IPJK2,M),K_OF(IPJKM2)),I_OF(IJK1))
!
            V_S_S = AVG_X(AVG_Z(V_S(JM_OF(IJK1),M),V_S(JM_OF(IJK2),M),K_OF(&
               JM_OF(IJK1))),AVG_Z(V_S(JM_OF(IPJKM2),M),V_S(JM_OF(IPJK2),M),&
               K_OF(JM_OF(IPJKM2))),I_OF(IJK1))
!
            V_S_E = AVG_Z(AVG_Y_N(V_S(JM_OF(IPJKM2),M),V_S(IPJKM2,M)),AVG_Y_N&
               (V_S(JM_OF(IPJK2),M),V_S(IPJK2,M)),K_OF(IPJKM2))
!
            V_S_W = AVG_Z(AVG_Y_N(V_S(JM_OF(IJK1),M),V_S(IJK1,M)),AVG_Y_N(&
               V_S(JM_OF(IJK2),M),V_S(IJK2,M)),K_OF(IJK1))
!
            V_S_T = AVG_X(AVG_Y_N(V_S(JM_OF(IJK2),M),V_S(IJK2,M)),AVG_Y_N(&
               V_S(JM_OF(IPJK2),M),V_S(IPJK2,M)),I_OF(IJK2))
!
            V_S_B = AVG_X(AVG_Y_N(V_S(JM_OF(IJK1),M),V_S(IJK1,M)),AVG_Y_N(&
               V_S(JM_OF(IPJKM2),M),V_S(IPJKM2,M)),I_OF(IJK1))
!
            W_S_N = ZERO
!
            W_S_S = ZERO
!
            W_S_E = ZERO
!
            W_S_W = ZERO
!
            W_S_T = AVG_X(AVG_Z_T(W_S(IJK1,M),W_S(IJK2,M)),AVG_Z_T(W_S(&
               IPJKM2,M),W_S(IPJK2,M)),I_OF(IJK2))
!
            W_S_B = AVG_X(AVG_Z_T(W_S(KM_OF(IJK1),M),W_S(IJK1,M)),AVG_Z_T(&
               W_S(KM_OF(IPJKM2),M),W_S(IPJKM2,M)),I_OF(IJK1))
!
            IF (CYLINDRICAL) THEN
               U_S_C = AVG_Z(U_S(IJK1,M),U_S(IJK2,M),K_OF(IJK1))
               W_S_C = AVG_X(W_S(IJK1,M),W_S(IPJKM2,M),I_OF(IJK1))
            ELSE
               U_S_C = ZERO
               W_S_C = ZERO
            ENDIF
!
            CALL SDDOTS (IJK1, FCELL, 'Y', U_S_N, U_S_S, U_S_E, U_S_W, U_S_T, &
               U_S_B, U_S_C, V_S_N, V_S_S, V_S_E, V_S_W, V_S_T, V_S_B, W_S_N, &
               W_S_S, W_S_E, W_S_W, W_S_T, W_S_B, W_S_C, DEL_DOT_U, S_DDOT_S, &
               S_DD)
!
!
         CASE ('B')
            IPJK2 = IP_OF(IJK2)
            IPJKP2 = KP_OF(IPJK2)
!
            U_S_N = AVG_Z(AVG_Y(U_S(IJK2,M),U_S(JP_OF(IJK2),M),J_OF(IJK2)),&
               AVG_Y(U_S(IJK1,M),U_S(JP_OF(IJK1),M),J_OF(IJK1)),K_OF(IJK2))
!
            U_S_S = AVG_Z(AVG_Y(U_S(JM_OF(IJK2),M),U_S(IJK2,M),J_OF(JM_OF(&
               IJK2))),AVG_Y(U_S(JM_OF(IJK1),M),U_S(IJK1,M),J_OF(JM_OF(IJK1))&
               ),K_OF(IJK2))
!
            U_S_E = AVG_Z(AVG_X_E(U_S(IJK2,M),U_S(IPJK2,M),I_OF(IPJK2)),&
               AVG_X_E(U_S(IJK1,M),U_S(IPJKP2,M),I_OF(IPJKP2)),K_OF(IPJK2))
!
            U_S_W = AVG_Z(AVG_X_E(U_S(IM_OF(IJK2),M),U_S(IJK2,M),I_OF(IJK2)),&
               AVG_X_E(U_S(IM_OF(IJK1),M),U_S(IJK1,M),I_OF(IJK1)),K_OF(IJK2))
!
            U_S_T = U_S(IJK1,M)
!
            U_S_B = U_S(IJK2,M)
!
            V_S_N = AVG_X(AVG_Z(V_S(IJK2,M),V_S(IJK1,M),K_OF(IJK2)),AVG_Z(&
               V_S(IPJK2,M),V_S(IPJKP2,M),K_OF(IPJK2)),I_OF(IJK2))
!
            V_S_S = AVG_X(AVG_Z(V_S(JM_OF(IJK2),M),V_S(JM_OF(IJK1),M),K_OF(&
               JM_OF(IJK2))),AVG_Z(V_S(JM_OF(IPJK2),M),V_S(JM_OF(IPJKP2),M),&
               K_OF(JM_OF(IPJK2))),I_OF(IJK2))
!
            V_S_E = AVG_Z(AVG_Y_N(V_S(JM_OF(IPJK2),M),V_S(IPJK2,M)),AVG_Y_N(&
               V_S(JM_OF(IPJKP2),M),V_S(IPJKP2,M)),K_OF(IPJK2))
!
            V_S_W = AVG_Z(AVG_Y_N(V_S(JM_OF(IJK2),M),V_S(IJK2,M)),AVG_Y_N(&
               V_S(JM_OF(IJK1),M),V_S(IJK1,M)),K_OF(IJK2))
!
            V_S_T = AVG_X(AVG_Y_N(V_S(JM_OF(IJK1),M),V_S(IJK1,M)),AVG_Y_N(&
               V_S(JM_OF(IPJKP2),M),V_S(IPJKP2,M)),I_OF(IJK1))
!
            V_S_B = AVG_X(AVG_Y_N(V_S(JM_OF(IJK2),M),V_S(IJK2,M)),AVG_Y_N(&
               V_S(JM_OF(IPJK2),M),V_S(IPJK2,M)),I_OF(IJK2))
!
            W_S_N = ZERO
!
            W_S_S = ZERO
!
            W_S_E = ZERO
!
            W_S_W = ZERO
!
            W_S_T = AVG_X(AVG_Z_T(W_S(IJK2,M),W_S(IJK1,M)),AVG_Z_T(W_S(IPJK2&
               ,M),W_S(IPJKP2,M)),I_OF(IJK1))
!
            W_S_B = AVG_X(AVG_Z_T(W_S(KM_OF(IJK2),M),W_S(IJK2,M)),AVG_Z_T(&
               W_S(KM_OF(IPJK2),M),W_S(IPJK2,M)),I_OF(IJK2))
!
            IF (CYLINDRICAL) THEN
               U_S_C = AVG_Z(U_S(IJK2,M),U_S(IJK1,M),K_OF(IJK2))
               W_S_C = AVG_X(W_S(IJK2,M),W_S(IPJK2,M),I_OF(IJK2))
            ELSE
               U_S_C = ZERO
               W_S_C = ZERO
            ENDIF
!
            CALL SDDOTS (IJK2, FCELL, 'Y', U_S_N, U_S_S, U_S_E, U_S_W, U_S_T, &
               U_S_B, U_S_C, V_S_N, V_S_S, V_S_E, V_S_W, V_S_T, V_S_B, W_S_N, &
               W_S_S, W_S_E, W_S_W, W_S_T, W_S_B, W_S_C, DEL_DOT_U, S_DDOT_S, &
               S_DD)
!
         END SELECT
!
      CASE ('V')
         SELECT CASE (TRIM(FCELL))
         CASE ('T')
            IJPK2 = JP_OF(IJK2)
            IJPKM2 = KM_OF(IJPK2)
!
            U_S_N = AVG_Z(AVG_X_E(U_S(IM_OF(IJPKM2),M),U_S(IJPKM2,M),I_OF(&
               IJPKM2)),AVG_X_E(U_S(IM_OF(IJPK2),M),U_S(IJPK2,M),I_OF(IJPK2))&
               ,K_OF(IJPKM2))
!
            U_S_S = AVG_Z(AVG_X_E(U_S(IM_OF(IJK1),M),U_S(IJK1,M),I_OF(IJK1)),&
               AVG_X_E(U_S(IM_OF(IJK2),M),U_S(IJK2,M),I_OF(IJK2)),K_OF(IJK1))
!
            U_S_E = AVG_Z(AVG_Y(U_S(IJK1,M),U_S(IJPKM2,M),J_OF(IJK1)),AVG_Y(&
               U_S(IJK2,M),U_S(IJPK2,M),J_OF(IJK2)),K_OF(IJK1))
!
            U_S_W = AVG_Z(AVG_Y(U_S(IM_OF(IJK1),M),U_S(IM_OF(IJPKM2),M),J_OF(&
               IM_OF(IJK1))),AVG_Y(U_S(IM_OF(IJK2),M),U_S(IM_OF(IJPK2),M),&
               J_OF(IM_OF(IJK2))),K_OF(IJK1))
!
            U_S_T = AVG_Y(AVG_X_E(U_S(IM_OF(IJK2),M),U_S(IJK2,M),I_OF(IJK2)),&
               AVG_X_E(U_S(IM_OF(IJPK2),M),U_S(IJPK2,M),I_OF(IJPK2)),J_OF(&
               IJK2))
!
            U_S_B = AVG_Y(AVG_X_E(U_S(IM_OF(IJK1),M),U_S(IJK1,M),I_OF(IJK1)),&
               AVG_X_E(U_S(IM_OF(IJPKM2),M),U_S(IJPKM2,M),I_OF(IJPKM2)),J_OF(&
               IJK1))
!
            V_S_N = AVG_Z(AVG_Y_N(V_S(IJK1,M),V_S(IJPKM2,M)),AVG_Y_N(V_S(&
               IJK2,M),V_S(IJPK2,M)),K_OF(IJPKM2))
!
            V_S_S = AVG_Z(AVG_Y_N(V_S(JM_OF(IJK1),M),V_S(IJK1,M)),AVG_Y_N(&
               V_S(JM_OF(IJK2),M),V_S(IJK2,M)),K_OF(IJK1))
!
            V_S_E = AVG_Z(AVG_X(V_S(IJK1,M),V_S(IP_OF(IJK1),M),I_OF(IJK1)),&
               AVG_X(V_S(IJK2,M),V_S(IP_OF(IJK2),M),I_OF(IJK2)),K_OF(IJK1))
!
            V_S_W = AVG_Z(AVG_X(V_S(IM_OF(IJK1),M),V_S(IJK1,M),I_OF(IM_OF(&
               IJK1))),AVG_X(V_S(IM_OF(IJK2),M),V_S(IJK2,M),I_OF(IM_OF(IJK2))&
               ),K_OF(IJK1))
!
            V_S_T = V_S(IJK2,M)
!
            V_S_B = V_S(IJK1,M)
!
            W_S_N = ZERO
!
            W_S_S = ZERO
!
            W_S_E = ZERO
!
            W_S_W = ZERO
!
            W_S_T = AVG_Y(AVG_Z_T(W_S(IJK1,M),W_S(IJK2,M)),AVG_Z_T(W_S(&
               IJPKM2,M),W_S(IJPK2,M)),J_OF(IJK2))
!
            W_S_B = AVG_Y(AVG_Z_T(W_S(KM_OF(IJK1),M),W_S(IJK1,M)),AVG_Z_T(&
               W_S(KM_OF(IJPKM2),M),W_S(IJPKM2,M)),J_OF(IJK1))
!
            IF (CYLINDRICAL) THEN
               U_S_C = AVG_Z(U_S_B,U_S_T,K_OF(IJK1))
               W_S_C = AVG_Y(W_S(IJK1,M),W_S(IJPKM2,M),J_OF(IJK1))
            ELSE
               U_S_C = ZERO
               W_S_C = ZERO
            ENDIF
!
            CALL SDDOTS (IJK1, FCELL, 'N', U_S_N, U_S_S, U_S_E, U_S_W, U_S_T, &
               U_S_B, U_S_C, V_S_N, V_S_S, V_S_E, V_S_W, V_S_T, V_S_B, W_S_N, &
               W_S_S, W_S_E, W_S_W, W_S_T, W_S_B, W_S_C, DEL_DOT_U, S_DDOT_S, &
               S_DD)
!
!
         CASE ('B')
            IJPK2 = JP_OF(IJK2)
            IJPKP2 = KP_OF(IJPK2)
!
            U_S_N = AVG_Z(AVG_X_E(U_S(IM_OF(IJPK2),M),U_S(IJPK2,M),I_OF(IJPK2&
               )),AVG_X_E(U_S(IM_OF(IJPKP2),M),U_S(IJPKP2,M),I_OF(IJPKP2)),&
               K_OF(IJPK2))
!
            U_S_S = AVG_Z(AVG_X_E(U_S(IM_OF(IJK2),M),U_S(IJK2,M),I_OF(IJK2)),&
               AVG_X_E(U_S(IM_OF(IJK1),M),U_S(IJK1,M),I_OF(IJK1)),K_OF(IJK2))
!
            U_S_E = AVG_Z(AVG_Y(U_S(IJK2,M),U_S(IJPK2,M),J_OF(IJK2)),AVG_Y(&
               U_S(IJK1,M),U_S(IJPKP2,M),J_OF(IJK1)),K_OF(IJK2))
!
            U_S_W = AVG_Z(AVG_Y(U_S(IM_OF(IJK2),M),U_S(IM_OF(IJPK2),M),J_OF(&
               IM_OF(IJK2))),AVG_Y(U_S(IM_OF(IJK1),M),U_S(IM_OF(IJPKP2),M),&
               J_OF(IM_OF(IJK1))),K_OF(IJK2))
!
            U_S_T = AVG_Y(AVG_X_E(U_S(IM_OF(IJK1),M),U_S(IJK1,M),I_OF(IJK1)),&
               AVG_X_E(U_S(IM_OF(IJPKP2),M),U_S(IJPKP2,M),I_OF(IJPKP2)),J_OF(&
               IJK1))
!
            U_S_B = AVG_Y(AVG_X_E(U_S(IM_OF(IJK2),M),U_S(IJK2,M),I_OF(IJK2)),&
               AVG_X_E(U_S(IM_OF(IJPK2),M),U_S(IJPK2,M),I_OF(IJPK2)),J_OF(&
               IJK2))
!
            V_S_N = AVG_Z(AVG_Y_N(V_S(IJK2,M),V_S(IJPK2,M)),AVG_Y_N(V_S(IJK1&
               ,M),V_S(IJPKP2,M)),K_OF(IJPK2))
!
            V_S_S = AVG_Z(AVG_Y_N(V_S(JM_OF(IJK2),M),V_S(IJK2,M)),AVG_Y_N(&
               V_S(JM_OF(IJK1),M),V_S(IJK1,M)),K_OF(IJK2))
!
            V_S_E = AVG_Z(AVG_X(V_S(IJK2,M),V_S(IP_OF(IJK2),M),I_OF(IJK2)),&
               AVG_X(V_S(IJK1,M),V_S(IP_OF(IJK1),M),I_OF(IJK1)),K_OF(IJK2))
!
            V_S_W = AVG_Z(AVG_X(V_S(IM_OF(IJK2),M),V_S(IJK2,M),I_OF(IM_OF(&
               IJK2))),AVG_X(V_S(IM_OF(IJK1),M),V_S(IJK1,M),I_OF(IM_OF(IJK1))&
               ),K_OF(IJK2))
!
            V_S_T = V_S(IJK1,M)
!
            V_S_B = V_S(IJK2,M)
!
            W_S_N = ZERO
!
            W_S_S = ZERO
!
            W_S_E = ZERO
!
            W_S_W = ZERO
!
            W_S_T = AVG_Y(AVG_Z_T(W_S(IJK2,M),W_S(IJK1,M)),AVG_Z_T(W_S(IJPK2&
               ,M),W_S(IJPKP2,M)),J_OF(IJK1))
!
            W_S_B = AVG_Y(AVG_Z_T(W_S(KM_OF(IJK2),M),W_S(IJK2,M)),AVG_Z_T(&
               W_S(KM_OF(IJPK2),M),W_S(IJPK2,M)),J_OF(IJK2))
!
            IF (CYLINDRICAL) THEN
               U_S_C = AVG_Z(U_S_B,U_S_T,K_OF(IJK2))
               W_S_C = AVG_Y(W_S(IJK2,M),W_S(IJPK2,M),J_OF(IJK2))
            ELSE
               U_S_C = ZERO
               W_S_C = ZERO
            ENDIF
!
            CALL SDDOTS (IJK2, FCELL, 'N', U_S_N, U_S_S, U_S_E, U_S_W, U_S_T, &
               U_S_B, U_S_C, V_S_N, V_S_S, V_S_E, V_S_W, V_S_T, V_S_B, W_S_N, &
               W_S_S, W_S_E, W_S_W, W_S_T, W_S_B, W_S_C, DEL_DOT_U, S_DDOT_S, &
               S_DD)
!
!
         CASE ('E')
            IJPK2 = JP_OF(IJK2)
            IMJPK2 = IM_OF(IJPK2)
!
            U_S_N = ZERO
!
            U_S_S = ZERO
!
            U_S_E = AVG_Y(AVG_X_E(U_S(IJK1,M),U_S(IJK2,M),I_OF(IJK2)),AVG_X_E&
               (U_S(IMJPK2,M),U_S(IJPK2,M),I_OF(IJPK2)),J_OF(IJK2))
!
            U_S_W = AVG_Y(AVG_X_E(U_S(IM_OF(IJK1),M),U_S(IJK1,M),I_OF(IJK1)),&
               AVG_X_E(U_S(IM_OF(IMJPK2),M),U_S(IMJPK2,M),I_OF(IMJPK2)),J_OF(&
               IJK1))
!
            U_S_T = ZERO
!
            U_S_B = ZERO
!
            V_S_N = AVG_X(AVG_Y_N(V_S(IJK1,M),V_S(IMJPK2,M)),AVG_Y_N(V_S(&
               IJK2,M),V_S(IJPK2,M)),I_OF(IMJPK2))
!
            V_S_S = AVG_X(AVG_Y_N(V_S(JM_OF(IJK1),M),V_S(IJK1,M)),AVG_Y_N(&
               V_S(JM_OF(IJK2),M),V_S(IJK2,M)),I_OF(IJK1))
!
            V_S_E = V_S(IJK2,M)
!
            V_S_W = V_S(IJK1,M)
!
            V_S_T = AVG_X(AVG_Z(V_S(IJK1,M),V_S(KP_OF(IJK1),M),K_OF(IJK1)),&
               AVG_Z(V_S(IJK2,M),V_S(KP_OF(IJK2),M),K_OF(IJK2)),I_OF(IJK1))
!
            V_S_B = AVG_X(AVG_Z(V_S(KM_OF(IJK1),M),V_S(IJK1,M),K_OF(KM_OF(&
               IJK1))),AVG_Z(V_S(KM_OF(IJK2),M),V_S(IJK2,M),K_OF(KM_OF(IJK2))&
               ),I_OF(IJK1))
!
            W_S_N = AVG_X(AVG_Z_T(W_S(KM_OF(IMJPK2),M),W_S(IMJPK2,M)),AVG_Z_T&
               (W_S(KM_OF(IJPK2),M),W_S(IJPK2,M)),I_OF(IMJPK2))
!
            W_S_S = AVG_X(AVG_Z_T(W_S(KM_OF(IJK1),M),W_S(IJK1,M)),AVG_Z_T(&
               W_S(KM_OF(IJK2),M),W_S(IJK2,M)),I_OF(IJK1))
!
            W_S_E = AVG_Y(AVG_Z_T(W_S(KM_OF(IJK2),M),W_S(IJK2,M)),AVG_Z_T(&
               W_S(KM_OF(IJPK2),M),W_S(IJPK2,M)),J_OF(IJK2))
!
            W_S_W = AVG_Y(AVG_Z_T(W_S(KM_OF(IJK1),M),W_S(IJK1,M)),AVG_Z_T(&
               W_S(KM_OF(IMJPK2),M),W_S(IMJPK2,M)),J_OF(IJK1))
!
            W_S_T = AVG_X(AVG_Y(W_S(IJK1,M),W_S(IMJPK2,M),J_OF(IJK1)),AVG_Y(&
               W_S(IJK2,M),W_S(IJPK2,M),J_OF(IJK2)),I_OF(IJK1))
!
            W_S_B = AVG_X(AVG_Y(W_S(KM_OF(IJK1),M),W_S(KM_OF(IMJPK2),M),J_OF(&
               KM_OF(IJK1))),AVG_Y(W_S(KM_OF(IJK2),M),W_S(KM_OF(IJPK2),M),&
               J_OF(KM_OF(IJK2))),I_OF(IJK1))
!
            IF (CYLINDRICAL) THEN
               U_S_C = AVG_Y(U_S(IJK1,M),U_S(IMJPK2,M),J_OF(IJK1))
               W_S_C = AVG_X(W_S_W,W_S_E,I_OF(IJK1))
            ELSE
               U_S_C = ZERO
               W_S_C = ZERO
            ENDIF
!
            CALL SDDOTS (IJK1, FCELL, 'Y', U_S_N, U_S_S, U_S_E, U_S_W, U_S_T, &
               U_S_B, U_S_C, V_S_N, V_S_S, V_S_E, V_S_W, V_S_T, V_S_B, W_S_N, &
               W_S_S, W_S_E, W_S_W, W_S_T, W_S_B, W_S_C, DEL_DOT_U, S_DDOT_S, &
               S_DD)
!
!
         CASE ('W')
            IJPK2 = JP_OF(IJK2)
            IPJPK2 = IP_OF(IJPK2)
!
            U_S_N = ZERO
!
            U_S_S = ZERO
!
            U_S_E = AVG_Y(AVG_X_E(U_S(IJK2,M),U_S(IJK1,M),I_OF(IJK1)),AVG_X_E&
               (U_S(IJPK2,M),U_S(IPJPK2,M),I_OF(IPJPK2)),J_OF(IJK1))
!
            U_S_W = AVG_Y(AVG_X_E(U_S(IM_OF(IJK2),M),U_S(IJK2,M),I_OF(IJK2)),&
               AVG_X_E(U_S(IM_OF(IJPK2),M),U_S(IJPK2,M),I_OF(IJPK2)),J_OF(&
               IJK2))
!
            U_S_T = ZERO
!
            U_S_B = ZERO
!
            V_S_N = AVG_X(AVG_Y_N(V_S(IJK2,M),V_S(IJPK2,M)),AVG_Y_N(V_S(IJK1&
               ,M),V_S(IPJPK2,M)),I_OF(IJPK2))
!
            V_S_S = AVG_X(AVG_Y_N(V_S(JM_OF(IJK2),M),V_S(IJK2,M)),AVG_Y_N(&
               V_S(JM_OF(IJK1),M),V_S(IJK1,M)),I_OF(IJK2))
!
            V_S_E = V_S(IJK1,M)
!
            V_S_W = V_S(IJK2,M)
!
            V_S_T = AVG_X(AVG_Z(V_S(IJK2,M),V_S(KP_OF(IJK2),M),K_OF(IJK2)),&
               AVG_Z(V_S(IJK1,M),V_S(KP_OF(IJK1),M),K_OF(IJK1)),I_OF(IJK2))
!
            V_S_B = AVG_X(AVG_Z(V_S(KM_OF(IJK2),M),V_S(IJK2,M),K_OF(KM_OF(&
               IJK2))),AVG_Z(V_S(KM_OF(IJK1),M),V_S(IJK1,M),K_OF(KM_OF(IJK1))&
               ),I_OF(IJK2))
!
            W_S_N = AVG_X(AVG_Z_T(W_S(KM_OF(IJPK2),M),W_S(IJPK2,M)),AVG_Z_T(&
               W_S(KM_OF(IPJPK2),M),W_S(IPJPK2,M)),I_OF(IJPK2))
!
            W_S_S = AVG_X(AVG_Z_T(W_S(KM_OF(IJK2),M),W_S(IJK2,M)),AVG_Z_T(&
               W_S(KM_OF(IJK1),M),W_S(IJK1,M)),I_OF(IJK2))
!
            W_S_E = AVG_Y(AVG_Z_T(W_S(KM_OF(IJK1),M),W_S(IJK1,M)),AVG_Z_T(&
               W_S(KM_OF(IPJPK2),M),W_S(IPJPK2,M)),J_OF(IJK1))
!
            W_S_W = AVG_Y(AVG_Z_T(W_S(KM_OF(IJK2),M),W_S(IJK2,M)),AVG_Z_T(&
               W_S(KM_OF(IJPK2),M),W_S(IJPK2,M)),J_OF(IJK2))
!
            W_S_T = AVG_X(AVG_Y(W_S(IJK2,M),W_S(IJPK2,M),J_OF(IJK2)),AVG_Y(&
               W_S(IJK1,M),W_S(IPJPK2,M),J_OF(IJK1)),I_OF(IJK2))
!
            W_S_B = AVG_X(AVG_Y(W_S(KM_OF(IJK2),M),W_S(KM_OF(IJPK2),M),J_OF(&
               KM_OF(IJK2))),AVG_Y(W_S(KM_OF(IJK1),M),W_S(KM_OF(IPJPK2),M),&
               J_OF(KM_OF(IJK1))),I_OF(IJK2))
!
            IF (CYLINDRICAL) THEN
               U_S_C = AVG_Y(U_S(IJK2,M),U_S(IJPK2,M),J_OF(IJK2))
               W_S_C = AVG_X(W_S_W,W_S_E,I_OF(IJK2))
            ELSE
               U_S_C = ZERO
               W_S_C = ZERO
            ENDIF
!
            CALL SDDOTS (IJK2, FCELL, 'Y', U_S_N, U_S_S, U_S_E, U_S_W, U_S_T, &
               U_S_B, U_S_C, V_S_N, V_S_S, V_S_E, V_S_W, V_S_T, V_S_B, W_S_N, &
               W_S_S, W_S_E, W_S_W, W_S_T, W_S_B, W_S_C, DEL_DOT_U, S_DDOT_S, &
               S_DD)
!
         END SELECT
!
      CASE ('W')
         SELECT CASE (TRIM(FCELL))
         CASE ('N')
            IJKP2 = KP_OF(IJK2)
            IJMKP2 = JM_OF(IJKP2)
!
            U_S_N = AVG_Z(AVG_X_E(U_S(IM_OF(IJK2),M),U_S(IJK2,M),I_OF(IJK2)),&
               AVG_X_E(U_S(IM_OF(IJKP2),M),U_S(IJKP2,M),I_OF(IJKP2)),K_OF(&
               IJK2))
!
            U_S_S = AVG_Z(AVG_X_E(U_S(IM_OF(IJK1),M),U_S(IJK1,M),I_OF(IJK1)),&
               AVG_X_E(U_S(IM_OF(IJMKP2),M),U_S(IJMKP2,M),I_OF(IJMKP2)),K_OF(&
               IJK1))
!
            U_S_E = AVG_Z(AVG_Y(U_S(IJK1,M),U_S(IJK2,M),J_OF(IJK1)),AVG_Y(&
               U_S(IJMKP2,M),U_S(IJKP2,M),J_OF(IJMKP2)),K_OF(IJK1))
!
            U_S_W = AVG_Z(AVG_Y(U_S(IM_OF(IJK1),M),U_S(IM_OF(IJK2),M),J_OF(&
               IM_OF(IJK1))),AVG_Y(U_S(IM_OF(IJMKP2),M),U_S(IM_OF(IJKP2),M),&
               J_OF(IM_OF(IJMKP2))),K_OF(IJK1))
!
            U_S_T = AVG_Y(AVG_X_E(U_S(IM_OF(IJMKP2),M),U_S(IJMKP2,M),I_OF(&
               IJMKP2)),AVG_X_E(U_S(IM_OF(IJKP2),M),U_S(IJKP2,M),I_OF(IJKP2))&
               ,J_OF(IJMKP2))
!
            U_S_B = AVG_Y(AVG_X_E(U_S(IM_OF(IJK1),M),U_S(IJK1,M),I_OF(IJK1)),&
               AVG_X_E(U_S(IM_OF(IJK2),M),U_S(IJK2,M),I_OF(IJK2)),J_OF(IJK1))
!
            V_S_N = AVG_Z(AVG_Y_N(V_S(IJK1,M),V_S(IJK2,M)),AVG_Y_N(V_S(&
               IJMKP2,M),V_S(IJKP2,M)),K_OF(IJK2))
!
            V_S_S = AVG_Z(AVG_Y_N(V_S(JM_OF(IJK1),M),V_S(IJK1,M)),AVG_Y_N(&
               V_S(JM_OF(IJMKP2),M),V_S(IJMKP2,M)),K_OF(IJK1))
!
            V_S_E = ZERO
!
            V_S_W = ZERO
!
            V_S_T = ZERO
!
            V_S_B = ZERO
!
            W_S_N = W_S(IJK2,M)
!
            W_S_S = W_S(IJK1,M)
!
            W_S_E = AVG_Y(AVG_X(W_S(IJK1,M),W_S(IP_OF(IJK1),M),I_OF(IJK1)),&
               AVG_X(W_S(IJK2,M),W_S(IP_OF(IJK2),M),I_OF(IJK2)),J_OF(IJK1))
!
            W_S_W = AVG_Y(AVG_X(W_S(IM_OF(IJK1),M),W_S(IJK1,M),I_OF(IM_OF(&
               IJK1))),AVG_X(W_S(IM_OF(IJK2),M),W_S(IJK2,M),I_OF(IM_OF(IJK2))&
               ),J_OF(IJK1))
!
            W_S_T = AVG_Y(AVG_Z_T(W_S(IJK1,M),W_S(IJMKP2,M)),AVG_Z_T(W_S(&
               IJK2,M),W_S(IJKP2,M)),J_OF(IJMKP2))
!
            W_S_B = AVG_Y(AVG_Z_T(W_S(KM_OF(IJK1),M),W_S(IJK1,M)),AVG_Z_T(&
               W_S(KM_OF(IJK2),M),W_S(IJK2,M)),J_OF(IJK1))
!
            IF (CYLINDRICAL) THEN
               U_S_C = AVG_Z(U_S_B,U_S_T,K_OF(IJK1))
               W_S_C = AVG_Y(W_S(IJK1,M),W_S(IJK2,M),J_OF(IJK1))
            ELSE
               U_S_C = ZERO
               W_S_C = ZERO
            ENDIF
!
            CALL SDDOTS (IJK1, FCELL, 'N', U_S_N, U_S_S, U_S_E, U_S_W, U_S_T, &
               U_S_B, U_S_C, V_S_N, V_S_S, V_S_E, V_S_W, V_S_T, V_S_B, W_S_N, &
               W_S_S, W_S_E, W_S_W, W_S_T, W_S_B, W_S_C, DEL_DOT_U, S_DDOT_S, &
               S_DD)
!
!
         CASE ('S')
            IJKP2 = KP_OF(IJK2)
            IJPKP2 = JP_OF(IJKP2)
!
            U_S_N = AVG_Z(AVG_X_E(U_S(IM_OF(IJK1),M),U_S(IJK1,M),I_OF(IJK1)),&
               AVG_X_E(U_S(IM_OF(IJPKP2),M),U_S(IJPKP2,M),I_OF(IJPKP2)),K_OF(&
               IJK1))
!
            U_S_S = AVG_Z(AVG_X_E(U_S(IM_OF(IJK2),M),U_S(IJK2,M),I_OF(IJK2)),&
               AVG_X_E(U_S(IM_OF(IJKP2),M),U_S(IJKP2,M),I_OF(IJKP2)),K_OF(&
               IJK2))
!
            U_S_E = AVG_Z(AVG_Y(U_S(IJK2,M),U_S(IJK1,M),J_OF(IJK2)),AVG_Y(&
               U_S(IJKP2,M),U_S(IJPKP2,M),J_OF(IJKP2)),K_OF(IJK2))
!
            U_S_W = AVG_Z(AVG_Y(U_S(IM_OF(IJK2),M),U_S(IM_OF(IJK1),M),J_OF(&
               IM_OF(IJK2))),AVG_Y(U_S(IM_OF(IJKP2),M),U_S(IM_OF(IJPKP2),M),&
               J_OF(IM_OF(IJKP2))),K_OF(IJK2))
!
            U_S_T = AVG_Y(AVG_X_E(U_S(IM_OF(IJKP2),M),U_S(IJKP2,M),I_OF(IJKP2&
               )),AVG_X_E(U_S(IM_OF(IJPKP2),M),U_S(IJPKP2,M),I_OF(IJPKP2)),&
               J_OF(IJKP2))
!
            U_S_B = AVG_Y(AVG_X_E(U_S(IM_OF(IJK2),M),U_S(IJK2,M),I_OF(IJK2)),&
               AVG_X_E(U_S(IM_OF(IJK1),M),U_S(IJK1,M),I_OF(IJK1)),J_OF(IJK2))
!
            V_S_N = AVG_Z(AVG_Y_N(V_S(IJK2,M),V_S(IJK1,M)),AVG_Y_N(V_S(IJKP2&
               ,M),V_S(IJPKP2,M)),K_OF(IJK1))
!
            V_S_S = AVG_Z(AVG_Y_N(V_S(JM_OF(IJK2),M),V_S(IJK2,M)),AVG_Y_N(&
               V_S(JM_OF(IJKP2),M),V_S(IJKP2,M)),K_OF(IJK2))
!
            V_S_E = ZERO
!
            V_S_W = ZERO
!
            V_S_T = ZERO
!
            V_S_B = ZERO
!
            W_S_N = W_S(IJK1,M)
!
            W_S_S = W_S(IJK2,M)
!
            W_S_E = AVG_Y(AVG_X(W_S(IJK2,M),W_S(IP_OF(IJK2),M),I_OF(IJK2)),&
               AVG_X(W_S(IJK1,M),W_S(IP_OF(IJK1),M),I_OF(IJK1)),J_OF(IJK2))
!
            W_S_W = AVG_Y(AVG_X(W_S(IM_OF(IJK2),M),W_S(IJK2,M),I_OF(IM_OF(&
               IJK2))),AVG_X(W_S(IM_OF(IJK1),M),W_S(IJK1,M),I_OF(IM_OF(IJK1))&
               ),J_OF(IJK2))
!
            W_S_T = AVG_Y(AVG_Z_T(W_S(IJK2,M),W_S(IJKP2,M)),AVG_Z_T(W_S(IJK1&
               ,M),W_S(IJPKP2,M)),J_OF(IJKP2))
!
            W_S_B = AVG_Y(AVG_Z_T(W_S(KM_OF(IJK2),M),W_S(IJK2,M)),AVG_Z_T(&
               W_S(KM_OF(IJK1),M),W_S(IJK1,M)),J_OF(IJK2))
!
            IF (CYLINDRICAL) THEN
               U_S_C = AVG_Z(U_S_B,U_S_T,K_OF(IJK2))
               W_S_C = AVG_Y(W_S(IJK2,M),W_S(IJK1,M),J_OF(IJK2))
            ELSE
               U_S_C = ZERO
               W_S_C = ZERO
            ENDIF
!
            CALL SDDOTS (IJK2, FCELL, 'N', U_S_N, U_S_S, U_S_E, U_S_W, U_S_T, &
               U_S_B, U_S_C, V_S_N, V_S_S, V_S_E, V_S_W, V_S_T, V_S_B, W_S_N, &
               W_S_S, W_S_E, W_S_W, W_S_T, W_S_B, W_S_C, DEL_DOT_U, S_DDOT_S, &
               S_DD)
!
!
         CASE ('E')
            IJKP2 = KP_OF(IJK2)
            IMJKP2 = IM_OF(IJKP2)
!
            U_S_N = ZERO
!
            U_S_S = ZERO
!
            U_S_E = AVG_Z(AVG_X_E(U_S(IJK1,M),U_S(IJK2,M),I_OF(IJK2)),AVG_X_E&
               (U_S(IMJKP2,M),U_S(IJKP2,M),I_OF(IJKP2)),K_OF(IJK2))
!
            U_S_W = AVG_Z(AVG_X_E(U_S(IM_OF(IJK1),M),U_S(IJK1,M),I_OF(IJK1)),&
               AVG_X_E(U_S(IM_OF(IMJKP2),M),U_S(IMJKP2,M),I_OF(IMJKP2)),K_OF(&
               IJK1))
!
            U_S_T = ZERO
!
            U_S_B = ZERO
!
            V_S_N = AVG_X(AVG_Z(V_S(IJK1,M),V_S(IMJKP2,M),K_OF(IJK1)),AVG_Z(&
               V_S(IJK2,M),V_S(IJKP2,M),K_OF(IJK2)),I_OF(IJK1))
!
            V_S_S = AVG_X(AVG_Z(V_S(JM_OF(IJK1),M),V_S(JM_OF(IMJKP2),M),K_OF(&
               JM_OF(IJK1))),AVG_Z(V_S(JM_OF(IJK2),M),V_S(JM_OF(IJKP2),M),&
               K_OF(JM_OF(IJK2))),I_OF(IJK1))
!
            V_S_E = AVG_Z(AVG_Y_N(V_S(JM_OF(IJK2),M),V_S(IJK2,M)),AVG_Y_N(&
               V_S(JM_OF(IJKP2),M),V_S(IJKP2,M)),K_OF(IJK2))
!
            V_S_W = AVG_Z(AVG_Y_N(V_S(JM_OF(IJK1),M),V_S(IJK1,M)),AVG_Y_N(&
               V_S(JM_OF(IMJKP2),M),V_S(IMJKP2,M)),K_OF(IJK1))
!
            V_S_T = AVG_X(AVG_Y_N(V_S(JM_OF(IMJKP2),M),V_S(IMJKP2,M)),AVG_Y_N&
               (V_S(JM_OF(IJKP2),M),V_S(IJKP2,M)),I_OF(IMJKP2))
!
            V_S_B = AVG_X(AVG_Y_N(V_S(JM_OF(IJK1),M),V_S(IJK1,M)),AVG_Y_N(&
               V_S(JM_OF(IJK2),M),V_S(IJK2,M)),I_OF(IJK1))
!
            W_S_N = AVG_X(AVG_Y(W_S(IJK1,M),W_S(JP_OF(IJK1),M),J_OF(IJK1)),&
               AVG_Y(W_S(IJK2,M),W_S(JP_OF(IJK2),M),J_OF(IJK2)),I_OF(IJK1))
!
            W_S_S = AVG_X(AVG_Y(W_S(JM_OF(IJK1),M),W_S(IJK1,M),J_OF(JM_OF(&
               IJK1))),AVG_Y(W_S(JM_OF(IJK2),M),W_S(IJK2,M),J_OF(JM_OF(IJK2))&
               ),I_OF(IJK1))
!
            W_S_E = W_S(IJK2,M)
!
            W_S_W = W_S(IJK1,M)
!
            W_S_T = AVG_X(AVG_Z_T(W_S(IJK1,M),W_S(IMJKP2,M)),AVG_Z_T(W_S(&
               IJK2,M),W_S(IJKP2,M)),I_OF(IMJKP2))
!
            W_S_B = AVG_X(AVG_Z_T(W_S(KM_OF(IJK1),M),W_S(IJK1,M)),AVG_Z_T(&
               W_S(KM_OF(IJK2),M),W_S(IJK2,M)),I_OF(IJK1))
!
            IF (CYLINDRICAL) THEN
               U_S_C = AVG_Z(U_S(IJK1,M),U_S(IMJKP2,M),K_OF(IJK1))
               W_S_C = AVG_X(W_S(IJK1,M),W_S(IJK2,M),I_OF(IJK1))
            ELSE
               U_S_C = ZERO
               W_S_C = ZERO
            ENDIF
!
            CALL SDDOTS (IJK1, FCELL, 'Y', U_S_N, U_S_S, U_S_E, U_S_W, U_S_T, &
               U_S_B, U_S_C, V_S_N, V_S_S, V_S_E, V_S_W, V_S_T, V_S_B, W_S_N, &
               W_S_S, W_S_E, W_S_W, W_S_T, W_S_B, W_S_C, DEL_DOT_U, S_DDOT_S, &
               S_DD)
!
!
         CASE ('W')
            IJKP2 = KP_OF(IJK2)
            IPJKP2 = IP_OF(IJKP2)
!
            U_S_N = ZERO
!
            U_S_S = ZERO
!
            U_S_E = AVG_Z(AVG_X_E(U_S(IJK2,M),U_S(IJK1,M),I_OF(IJK1)),AVG_X_E&
               (U_S(IJKP2,M),U_S(IPJKP2,M),I_OF(IPJKP2)),K_OF(IJK1))
!
            U_S_W = AVG_Z(AVG_X_E(U_S(IM_OF(IJK2),M),U_S(IJK2,M),I_OF(IJK2)),&
               AVG_X_E(U_S(IM_OF(IJKP2),M),U_S(IJKP2,M),I_OF(IJKP2)),K_OF(&
               IJK2))
!
            U_S_T = ZERO
!
            U_S_B = ZERO
!
            V_S_N = AVG_X(AVG_Z(V_S(IJK2,M),V_S(IJKP2,M),K_OF(IJK2)),AVG_Z(&
               V_S(IJK1,M),V_S(IPJKP2,M),K_OF(IJK1)),I_OF(IJK2))
!
            V_S_S = AVG_X(AVG_Z(V_S(JM_OF(IJK2),M),V_S(JM_OF(IJKP2),M),K_OF(&
               JM_OF(IJK2))),AVG_Z(V_S(JM_OF(IJK1),M),V_S(JM_OF(IPJKP2),M),&
               K_OF(JM_OF(IJK1))),I_OF(IJK2))
!
            V_S_E = AVG_Z(AVG_Y_N(V_S(JM_OF(IJK1),M),V_S(IJK1,M)),AVG_Y_N(&
               V_S(JM_OF(IPJKP2),M),V_S(IPJKP2,M)),K_OF(IJK1))
!
            V_S_W = AVG_Z(AVG_Y_N(V_S(JM_OF(IJK2),M),V_S(IJK2,M)),AVG_Y_N(&
               V_S(JM_OF(IJKP2),M),V_S(IJKP2,M)),K_OF(IJK2))
!
            V_S_T = AVG_X(AVG_Y_N(V_S(JM_OF(IJKP2),M),V_S(IJKP2,M)),AVG_Y_N(&
               V_S(JM_OF(IPJKP2),M),V_S(IPJKP2,M)),I_OF(IJKP2))
!
            V_S_B = AVG_X(AVG_Y_N(V_S(JM_OF(IJK2),M),V_S(IJK2,M)),AVG_Y_N(&
               V_S(JM_OF(IJK1),M),V_S(IJK1,M)),I_OF(IJK2))
!
            W_S_N = AVG_X(AVG_Y(W_S(IJK2,M),W_S(JP_OF(IJK2),M),J_OF(IJK2)),&
               AVG_Y(W_S(IJK1,M),W_S(JP_OF(IJK1),M),J_OF(IJK1)),I_OF(IJK2))
!
            W_S_S = AVG_X(AVG_Y(W_S(JM_OF(IJK2),M),W_S(IJK2,M),J_OF(JM_OF(&
               IJK2))),AVG_Y(W_S(JM_OF(IJK1),M),W_S(IJK1,M),J_OF(JM_OF(IJK1))&
               ),I_OF(IJK2))
!
            W_S_E = W_S(IJK1,M)
!
            W_S_W = W_S(IJK2,M)
!
            W_S_T = AVG_X(AVG_Z_T(W_S(IJK2,M),W_S(IJKP2,M)),AVG_Z_T(W_S(IJK1&
               ,M),W_S(IPJKP2,M)),I_OF(IJKP2))
!
            W_S_B = AVG_X(AVG_Z_T(W_S(KM_OF(IJK2),M),W_S(IJK2,M)),AVG_Z_T(&
               W_S(KM_OF(IJK1),M),W_S(IJK1,M)),I_OF(IJK2))
!
            IF (CYLINDRICAL) THEN
               U_S_C = AVG_Z(U_S(IJK2,M),U_S(IJKP2,M),K_OF(IJK2))
               W_S_C = AVG_X(W_S(IJK2,M),W_S(IJK1,M),I_OF(IJK2))
            ELSE
               U_S_C = ZERO
               W_S_C = ZERO
            ENDIF
!
            CALL SDDOTS (IJK2, FCELL, 'Y', U_S_N, U_S_S, U_S_E, U_S_W, U_S_T, &
               U_S_B, U_S_C, V_S_N, V_S_S, V_S_E, V_S_W, V_S_T, V_S_B, W_S_N, &
               W_S_S, W_S_E, W_S_W, W_S_T, W_S_B, W_S_C, DEL_DOT_U, S_DDOT_S, &
               S_DD)
!
         END SELECT
      END SELECT
      RETURN
      END SUBROUTINE CALC_S_DDOT_S
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SUBROUTINE SDDOTS(IJK, NORMAL, DZVALUE,                C
!                 U_s_N, U_s_S, U_s_E, U_s_W, U_s_T, U_s_B, U_s_C,     C
!                 V_s_N, V_s_S, V_s_E,  V_s_W, V_s_T, V_s_B,           C
!                 W_s_N, W_s_S, W_s_E, W_s_W, W_s_T, W_s_B, W_s_C,     C
!                 DEL_DOT_U, S_DDOT_S)                                 C
!                                                                      C
!  Purpose: Calculate del.U (trace of D_s), S:S and S_xx, S_yy or S_zz C
!           at the boundary                                            C
!                                                                      C
!                                                                      C
!  Author: Anuj Srivastava, Princeton University      Date: 4-APR-98   C
!  Reviewer:                                          Date:            C
!                                                                      C
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
      SUBROUTINE SDDOTS(IJK, NORMAL, DZVALUE, U_S_N, U_S_S, U_S_E, U_S_W, U_S_T&
         , U_S_B, U_S_C, V_S_N, V_S_S, V_S_E, V_S_W, V_S_T, V_S_B, W_S_N, W_S_S&
         , W_S_E, W_S_W, W_S_T, W_S_B, W_S_C, DEL_DOT_U, S_DDOT_S, S_DD)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE constant
      USE fldvar
      USE geometry
      USE indices
      USE compar
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!                      Strain rate tensor components for mth solids phase
       DOUBLE PRECISION D_s(3,3)

!                      direction of normal to wall
       CHARACTER        NORMAL

!                      tells what value should odelta_Z should take
       CHARACTER        DZVALUE

!  The location where D is calculated in located at the center of 4
!  cells - 2 fluid cells and 2 wall cells. Coordinates of this are i,j,k

!                      U_s at the north (i, j+1/2, k)
       DOUBLE PRECISION U_s_N
!
!                      U_s at the south (i, j-1/2, k)
       DOUBLE PRECISION U_s_S

!                      U_s at the east (i+1/2, j, k)
       DOUBLE PRECISION U_s_E
!
!                      U_s at the west (i-1/2, j, k)
       DOUBLE PRECISION U_s_W
!
!                      U_s at the top (i, j, k+1/2)
       DOUBLE PRECISION U_s_T
!
!                      U_s at the bottom (i, j, k-1/2)
       DOUBLE PRECISION U_s_B
!
!                      U_s at the center (i, j, k)
!                      Calculated for Cylindrical coordinates only.
       DOUBLE PRECISION U_s_C
!
!                      V_s at the north (i, j+1/2, k)
       DOUBLE PRECISION V_s_N
!
!                      V_s at the south (i, j-1/2, k)
       DOUBLE PRECISION V_s_S
!
!                      V_s at the east (i+1/2, j, k)
       DOUBLE PRECISION V_s_E
!
!                      V_s at the west (i-1/2, j, k)
       DOUBLE PRECISION V_s_W

!                      V_s at the top (i, j, k+1/2)
       DOUBLE PRECISION V_s_T
!
!                      V_s at the bottom (i, j, k-1/2)
       DOUBLE PRECISION V_s_B

!                      W_s at the north (i, j+1/2, k)
       DOUBLE PRECISION W_s_N
!
!                      W_s at the south (i, j-1/2, k)
       DOUBLE PRECISION W_s_S
!
!                      W_s at the east (i+1/2, j, k)
       DOUBLE PRECISION W_s_E
!
!                      W_s at the west (1-1/2, j, k)
       DOUBLE PRECISION W_s_W
!
!                      W_s at the top (i, j, k+1/2)
       DOUBLE PRECISION W_s_T
!
!                      W_s at the bottom (i, j, k-1/2)
       DOUBLE PRECISION W_s_B
!
!                      W_s at the center (i, j, k).
!                      Calculated for Cylindrical coordinates only.
       DOUBLE PRECISION W_s_C

!                      del.u
       DOUBLE PRECISION DEL_DOT_U

!                      S:S
       DOUBLE PRECISION S_DDOT_S

!                      S_dd (d is x,y, or z)
       DOUBLE PRECISION S_dd

!                      trace of D
       DOUBLE PRECISION TRACE_D

!                      trace of the square of D
       DOUBLE PRECISION TRACE_sD

!                      Local indices
       INTEGER          IJK, I, J, K, I1, I2

       DOUBLE PRECISION odelta_Z

!-----------------------------------------------
!
!         Define I, J, K
!     IJK is the cell whose coordinates are i-1/2, j-1/2 and k-1/2
!     (this cell actually has two of the above coordinates. The other
!     coordinate is i, j or k)
!     This is necessary because we are using oDX_E, oDY_N, and oDZ_T to
!     calculate D_s
!
!
      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)
!
      IF (DZVALUE == 'Y') THEN
         ODELTA_Z = (ODZ(K)+ODZ(K_OF(IP_OF(IJK))))/2D0
      ELSE
         ODELTA_Z = ODZ_T(K)
      ENDIF
!
!
!
!         Find components of Mth solids phase continuum strain rate
!         tensor, D_s, at i,j,k
!
      D_S(1,1) = (U_S_E - U_S_W)*ODX_E(I)
      D_S(1,2) = HALF*((U_S_N - U_S_S)*ODY_N(J)+(V_S_E-V_S_W)*ODX_E(I))
      D_S(1,3) = HALF*((W_S_E - W_S_W)*ODX_E(I)+(U_S_T-U_S_B)*(OX_E(I)*ODELTA_Z&
         )-W_S_C*OX_E(I))
      D_S(2,1) = D_S(1,2)
      D_S(2,2) = (V_S_N - V_S_S)*ODY_N(J)
      D_S(2,3)=HALF*((V_S_T-V_S_B)*(OX_E(I)*ODELTA_Z)+(W_S_N-W_S_S)*ODY_N(J))
      D_S(3,1) = D_S(1,3)
      D_S(3,2) = D_S(2,3)
      D_S(3,3) = (W_S_T - W_S_B)*(OX_E(I)*ODELTA_Z) + U_S_C*OX_E(I)
!
!         Calculate the trace of D_s
      TRACE_D = D_S(1,1) + D_S(2,2) + D_S(3,3)
!
!         Calculate trace of the square of D_s
      TRACE_SD = 0.D0                            !Initialize the totalizer
      DO I1 = 1, 3
         TRACE_SD = TRACE_SD + SUM(D_S(I1,:)*D_S(I1,:))
         I2 = 4
      END DO
      DEL_DOT_U = TRACE_D
!
      S_DDOT_S = TRACE_SD - (TRACE_D*TRACE_D)/3.D0
!
      S_DDOT_S = DMAX1(1D-10,S_DDOT_S)
!
      IF (NORMAL=='E' .OR. NORMAL=='W') THEN
         S_DD = D_S(1,1) - TRACE_D/3D0
!
      ELSE IF (NORMAL=='N' .OR. NORMAL=='S') THEN
         S_DD = D_S(2,2) - TRACE_D/3D0
!
      ELSE
         S_DD = D_S(3,3) - TRACE_D/3D0
      ENDIF
!
      RETURN
      END SUBROUTINE SDDOTS
