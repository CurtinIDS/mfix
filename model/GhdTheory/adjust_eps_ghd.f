!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ADJUST_EPS_GHD                                         C
!  Purpose: Eliminate the solids phases that occupy only very small    C
!           fractions of the computational cell volume for GHD theory  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 28-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 29-JAN-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Include CONSTANTS.INC                                      C
!  Author: M. Syamlal                                 Date: 7-FEB-92   C
!  Reviewer: P. Nicoletti                             Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References: None                                C
!                                                                      C
!  Variables referenced: FLAG, RO_s, ROP_s, U_s, V_s, W_s, EP_g        C
!                        IMAX1, JMAX1, KMAX1, MMAX, IMIN1, JMIN1, KMIN1C
!  Variables modified: ROP_s, U_s, V_s, W_s, EP_g, ROP_g, I, J, K, IJK,C
!                      RO_s                                            C
!                                                                      C
!  Local variables: EPSUM                                              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE ADJUST_EPS_GHD
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE toleranc
      USE constant
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE compar
      USE sendrecv
      USE ghdtheory
      USE fun_avg
      USE functions
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
      INTEGER          I, J, K, IJK, IJKE, IJKN, IJKT
!
!                      Solids phase
      INTEGER          M
!                      Sum of (very small) solids volume fractions that
!                      are set to zero.
      DOUBLE PRECISION epsMixE, epsMixN, epsMixT, epSolid, epSolidE(smax), epSolidN(smax), epSolidT(smax)
      LOGICAL          DiluteCellE, DiluteCellN, DiluteCellT
!
! First set solids cell-center density in very dilute cells to zero
!
      DO IJK = ijkstart3, ijkend3
               IF (FLUID_AT(IJK)) THEN
                  DO M = 1, SMAX
                     epSolid = ROP_S(IJK,M)/RO_S(IJK,M)
                     IF (epSolid < ZERO_EP_S) THEN

!  Remove solids in very small quantities and set solids velocity to zero
!  if there is outflow from the present cell.

                        ROP_S(IJK,MMAX) = ROP_S(IJK,MMAX) - ROP_S(IJK,M) ! mmax = mixture phase
                        EP_G(IJK) = EP_G(IJK) + epSolid
                        ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
                        ROP_S(IJK,M) = ZERO
                     ENDIF
                  END DO
               ENDIF
      END DO
!
! now set velocities at cell faces to zero
!
      DO IJK = ijkstart3, ijkend3
               I = I_OF(IJK)
               J = J_OF(IJK)
               K = K_OF(IJK)
               IJKE = EAST_OF(IJK)
               IJKN = NORTH_OF(IJK)
               IJKT = TOP_OF(IJK)
               IF (FLUID_AT(IJK)) THEN
                  epsMixE = ZERO
                  epsMixN = ZERO
                  epsMixT = ZERO
                  DiluteCellE = .FALSE.
                  DiluteCellN = .FALSE.
                  DiluteCellT = .FALSE.
                  DO M = 1, SMAX
                     epSolidE(M) = AVG_X(EP_S(IJK,M),EP_S(IJKE,M),I)
                     epSolidN(M) = AVG_Y(EP_S(IJK,M),EP_S(IJKN,M),J)
                     epSolidT(M) = AVG_Z(EP_S(IJK,M),EP_S(IJKT,M),K)

                     epsMixE = epsMixE + epSolidE(M)*RO_S(IJK,M)
                     epsMixN = epsMixN + epSolidN(M)*RO_S(IJK,M)
                     epsMixT = epsMixT + epSolidT(M)*RO_S(IJK,M)

                     IF (epSolidE(M) < ZERO_EP_S) THEN
                        epsMixE = epsMixE - epSolidE(M)*RO_S(IJK,M)
                        epSolidE(M) = ZERO
                        U_S(IJK,M) = ZERO
                        JoiX(IJK,M) = ZERO
                        DiluteCellE = .TRUE. ! to make sure at least one species is very dilute before
                                             ! re-calculating mixture velocity
                     ENDIF

                     IF (epSolidN(M) < ZERO_EP_S) THEN
                        epsMixN = epsMixN - epSolidN(M)*RO_S(IJK,M)
                        epSolidN(M) = ZERO
                        V_S(IJK,M) = ZERO
                        JoiY(IJK,M) = ZERO
                        DiluteCellN = .TRUE.
                     ENDIF

                     IF (epSolidT(M) < ZERO_EP_S) THEN
                        epsMixT = epsMixT - epSolidT(M)*RO_S(IJK,M)
                        epSolidT(M) = ZERO
                        W_S(IJK,M) = ZERO
                        JoiZ(IJK,M) = ZERO
                        DiluteCellT = .TRUE.
                     ENDIF
                  END DO
!
! compute corrected mixture velcoity and species mass fluxes based on their definition.
!
                     IF (epsMixE > ZERO .AND. DiluteCellE) THEN
                       U_S(IJK,MMAX) = ZERO
                       DO M = 1, SMAX
                         U_S(IJK,MMAX) = U_S(IJK,MMAX) + U_S(IJK,M)*epSolidE(M)*RO_S(IJK,M)
                       ENDDO
                       U_S(IJK,MMAX) = U_S(IJK,MMAX)/epsMixE
                       DO M = 1, SMAX
                         JoiX(IJK,M) = epSolidE(M)*RO_S(IJK,M)*(U_S(IJK,M)-U_S(IJK,MMAX))
                       ENDDO
                     ELSEIF(epsMixE == ZERO) THEN
                       U_S(IJK,MMAX) = ZERO
                     ENDIF

                     IF (epsMixN > ZERO .AND. DiluteCellN) THEN
                       V_S(IJK,MMAX) = ZERO
                       DO M = 1, SMAX
                         V_S(IJK,MMAX) = V_S(IJK,MMAX) + V_S(IJK,M)*epSolidN(M)*RO_S(IJK,M)
                       ENDDO
                       V_S(IJK,MMAX) = V_S(IJK,MMAX)/epsMixN
                       DO M = 1, SMAX
                         JoiY(IJK,M) = epSolidN(M)*RO_S(IJK,M)*(V_S(IJK,M)-V_S(IJK,MMAX))
                       ENDDO
                     ELSEIF(epsMixN == ZERO) THEN
                       V_S(IJK,MMAX) = ZERO
                     ENDIF

                     IF(.NOT.NO_K) THEN  ! for 3D cases only
                       IF (epsMixT > ZERO .AND. DiluteCellT) THEN
                         W_S(IJK,MMAX) = ZERO
                         DO M = 1, SMAX
                           W_S(IJK,MMAX) = W_S(IJK,MMAX) + W_S(IJK,M)*epSolidT(M)*RO_S(IJK,M)
                         ENDDO
                         W_S(IJK,MMAX) = W_S(IJK,MMAX)/epsMixT
                         DO M = 1, SMAX
                           JoiZ(IJK,M) = epSolidT(M)*RO_S(IJK,M)*(W_S(IJK,M)-W_S(IJK,MMAX))
                         ENDDO
                       ELSEIF(epsMixT == ZERO) THEN
                         W_S(IJK,MMAX) = ZERO
                       ENDIF
                     ENDIF
               ENDIF
      END DO

      RETURN
      END SUBROUTINE ADJUST_EPS_GHD

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 200 1119 Changed the limits for the triple loop, do k=kmin1,kmax1=>kstart1,kend1
!// 400 Added sendrecv module and send_recv calls for COMMunication
