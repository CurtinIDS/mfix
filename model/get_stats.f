!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_STATS(IER)                                         C
!  Purpose: Get statistics for stalled or diverged iterations          C
!                                                                      C
!  Author: M. Syamlal                                 Date: 11-FEB-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!                                                                      C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:  None                                              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE GET_STATS()
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE physprop
      USE geometry
      USE fldvar
      USE indices
      USE funits
      USE residual
      USE run
      USE compar
      USE functions

      use machine

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                      Indices
      INTEGER          I, J, K, IJK, M

      DOUBLE PRECISION tmp
!
!                      direction (1 - x, 2 - y, 3- z)
      INTEGER          Dir_g, Dir_s(DIMENSION_M)
!
!                      Courant number
      DOUBLE PRECISION NC_g, NC_s(DIMENSION_M)
!
!                      Maximum P_star
      DOUBLE PRECISION Ps
!
!                      Locations of maxima
      INTEGER          IJK_NC_g, IJK_Ps, IJK_NC_s(DIMENSION_M)


!
!-----------------------------------------------

      DIR_G = 1
      NC_G = LARGE_NUMBER
      IJK_NC_G = 0
!
      PS = 0.
      IJK_PS = 0
!
      DO IJK = IJKMIN1, IJKMAX1
         IF (FLUID_AT(IJK)) THEN
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
!
            IF (U_G(IJK) == ZERO) THEN
               TMP = LARGE_NUMBER
            ELSE IF (U_G(IJK) > ZERO) THEN
               TMP = DX(IP1(I))/U_G(IJK)
            ELSE
               TMP = -DX(I)/U_G(IJK)
            ENDIF
!
            IF (TMP < NC_G) THEN
               NC_G = TMP
               IJK_NC_G = IJK
               DIR_G = 1
            ENDIF
!
            IF (V_G(IJK) == ZERO) THEN
               TMP = LARGE_NUMBER
            ELSE IF (V_G(IJK) > ZERO) THEN
               TMP = DY(JP1(J))/V_G(IJK)
            ELSE
               TMP = -DY(J)/V_G(IJK)
            ENDIF
!
            IF (TMP < NC_G) THEN
               NC_G = TMP
               IJK_NC_G = IJK
               DIR_G = 2
            ENDIF
!
            IF (W_G(IJK) == ZERO) THEN
               TMP = LARGE_NUMBER
            ELSE IF (W_G(IJK) > ZERO) THEN
               TMP = DZ(KP1(K))*X(I)/W_G(IJK)
            ELSE
               TMP = -DZ(K)*X(I)/W_G(IJK)
            ENDIF
!
            IF (TMP < NC_G) THEN
               NC_G = TMP
               IJK_NC_G = IJK
               DIR_G = 3
            ENDIF
!
            IF (P_STAR(IJK) > PS) THEN
               PS = P_STAR(IJK)
               IJK_PS = IJK
            ENDIF
!
         ENDIF
      END DO
      M = 1
      IF (MMAX > 0) THEN
         DIR_S(:MMAX) = 1
         NC_S(:MMAX) = LARGE_NUMBER
         IJK_NC_S(:MMAX) = 0
         M = MMAX + 1
      ENDIF
      DO M = 1, MMAX
         DO IJK = IJKMIN1, IJKMAX1
            IF (FLUID_AT(IJK)) THEN
               I = I_OF(IJK)
               J = J_OF(IJK)
               K = K_OF(IJK)
               IF (U_S(IJK,M) == ZERO) THEN
                  TMP = LARGE_NUMBER
               ELSE IF (U_S(IJK,M) > ZERO) THEN
                  TMP = DX(IP1(I))/U_S(IJK,M)
               ELSE
                  TMP = -DX(I)/U_S(IJK,M)
               ENDIF
!
               IF (TMP < NC_S(M)) THEN
                  NC_S(M) = TMP
                  IJK_NC_S(M) = IJK
                  DIR_S(M) = 1
               ENDIF
!
               IF (V_S(IJK,M) == ZERO) THEN
                  TMP = LARGE_NUMBER
               ELSE IF (V_S(IJK,M) > ZERO) THEN
                  TMP = DY(JP1(J))/V_S(IJK,M)
               ELSE
                  TMP = -DY(J)/V_S(IJK,M)
               ENDIF
!
               IF (TMP < NC_S(M)) THEN
                  NC_S(M) = TMP
                  IJK_NC_S(M) = IJK
                  DIR_S(M) = 2
               ENDIF
!
               IF (W_S(IJK,M) == ZERO) THEN
                  TMP = LARGE_NUMBER
               ELSE IF (W_S(IJK,M) > ZERO) THEN
                  TMP = DZ(KP1(K))*X(I)/W_S(IJK,M)
               ELSE
                  TMP = -DZ(K)*X(I)/W_S(IJK,M)
               ENDIF
!
               IF (TMP < NC_S(M)) THEN
                  NC_S(M) = TMP
                  IJK_NC_S(M) = IJK
                  DIR_S(M) = 3
               ENDIF
!
            ENDIF
         END DO
      END DO
      CALL START_LOG
      IF(DMP_LOG)WRITE (UNIT_LOG, *) 'Gas phase:'
      IF(DMP_LOG)WRITE (UNIT_LOG, '(A, G12.3, A, I6, A, I1)') ' Minimum Courant No = ', &
         NC_G*ODT, '  Location = ', IJK_NC_G, '  Direction = ', DIR_G
!
      DO M = 1, MMAX
         IF(DMP_LOG)WRITE (UNIT_LOG, '(A, I2, A)') ' Solids phase (', M, '):'
         IF(DMP_LOG)WRITE (UNIT_LOG, '(A, G12.3, A, I6, A, I1)') ' Minimum Courant No = '&
            , NC_S(M)*ODT, '  Location = ', IJK_NC_S(M), '  Direction = ', &
            DIR_S(M)
      END DO
      IF(DMP_LOG)WRITE (UNIT_LOG, '(A, G12.3, A, I6)') ' Maximum P_star = ', PS, &
         '  Location = ', IJK_PS
!
!      IF(DMP_LOG)WRITE(UNIT_LOG,'(A, G12.3, A, I6)')
!     & " Maximum P_g residual = ", MAX_RESID(RESID_P, 0),
!     & "  Location = ", IJK_RESID(RESID_P, 0)
!      IF(DMP_LOG)WRITE(UNIT_LOG,'(A, G12.3, A, I6)')
!     & " Maximum P_s residual = ", MAX_RESID(RESID_p, 1),
!     & "  Location = ", IJK_RESID(RESID_p, 1)
!
      IF(DMP_LOG)WRITE (UNIT_LOG, *)
!
      CALL END_LOG
!
      RETURN
      END SUBROUTINE GET_STATS

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!//? DO loop limits still running between IJKMIN1, IJKMAX1
