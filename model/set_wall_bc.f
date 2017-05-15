!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_WALL_BC                                             C
!  Purpose: Set wall boundary conditions                               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 29-JAN-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Add calculations for mass outflow boundary condition       C
!  Author: M. Syamlal                                 Date: 23-OCT-92  C
!  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
!  Revision Number: 2                                                  C
!  Purpose: Revised for MFIX 2.0. This subroutine is different from    C
!           old set_wall_bc.                                           C
!  Author: M. Syamlal                                 Date: 18-JUL-96  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: BC_DEFINED, BC_TYPE, BC_JJ_PS, BC_I_w,        C
!                        BC_I_e, BC_J_s, BC_J_n, BC_K_b, BC_K_t,       C
!                        ISTART3, IEND3, JSTART3, JEND3, KSTART3,      C
!                        KEND3, ISTART2, IEND2, JSTART2, JEND2,        C
!                        KSTART2, KEND3, IMAX2, JMAX2, KMAX2, MMAX,    C
!                        W_g, W_S in fluid cell adjacent to wall cell  C
!  Variables modified: W_g, W_S in wall cell                           C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_WALL_BC()

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE bc
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE funits
      USE compar
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Local index for boundary condition
      INTEGER :: L
! indices
      INTEGER :: IJK, IPJK
! Starting & ending I index
      INTEGER :: I1, I2
! Starting & ending J index
      INTEGER :: J1, J2
! Starting and ending K index
      INTEGER :: K1, K2
!-----------------------------------------------


! Set the boundary conditions
      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L)) THEN

! The range of boundary cells
            I1 = BC_I_W(L)
            I2 = BC_I_E(L)
            J1 = BC_J_S(L)
            J2 = BC_J_N(L)
            K1 = BC_K_B(L)
            K2 = BC_K_T(L)

            SELECT CASE (BC_TYPE_ENUM(L))
               CASE (FREE_SLIP_WALL)
                  CALL SET_WALL_BC1 (I1, I2, J1, J2, K1, K2, &
                                     BC_JJ_PS(L))

               CASE (NO_SLIP_WALL)
                  CALL SET_WALL_BC1 (I1, I2, J1, J2, K1, K2, &
                                     BC_JJ_PS(L))

               CASE (PAR_SLIP_WALL)
! updating the boundary velocity may improve convergence
            END SELECT
         ENDIF
      ENDDO


! The above section did not address bc_type=undefined (which by default
! is either a ns wall, or if i=1 and cylindrical, a fs wall) or
! bc_type='dummy' conditions. The section below will handle both events
! since default_wall_at will register as true
      K1 = 1
      DO J1 = JSTART3, JEND3
         DO I1 = ISTART3, IEND3
            IF(K1.NE.KSTART2) EXIT
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) CALL SET_WALL_BC1 (I1, I1,&
                J1, J1, K1, K1, 0)
         ENDDO
      ENDDO

! top xy-plane
      K1 = KMAX2
      DO J1 = JSTART3, JEND3
         DO I1 = ISTART3, IEND3
            IF(K1.NE.KEND2) EXIT
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) CALL SET_WALL_BC1 (I1, I1, &
               J1, J1, K1, K1, 0)
         ENDDO
      ENDDO

! south xz-plane
      J1 = 1
      DO K1 = KSTART3, KEND3
         DO I1 = ISTART3, IEND3
            IF(J1.NE.JSTART2) EXIT
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) CALL SET_WALL_BC1 (I1, I1, &
                J1, J1, K1, K1, 0)
         ENDDO
      ENDDO

! north xz-plane
      J1 = JMAX2
      DO K1 = KSTART3, KEND3
         DO I1 = ISTART3, IEND3
            IF(J1.NE.JEND2) EXIT
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) CALL SET_WALL_BC1 (I1, I1, &
               J1, J1, K1, K1, 0)
         ENDDO
      ENDDO

! west zy-plane
      I1 = 1
      DO K1 = KSTART3, KEND3
         DO J1 = JSTART3, JEND3
            IF(I1.NE.ISTART2) EXIT
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) CALL SET_WALL_BC1 (I1, I1, &
                J1, J1, K1, K1, 0)

! For cylindrical coordinates the azimuthal component should be zero
! at center (this forces no-slip for the azimuthal component at what
! is a free-slip wall)
            IF (CYLINDRICAL .AND. XMIN==ZERO) THEN
               IPJK = IP_OF(IJK)
               W_G(IJK) = -W_G(IPJK)
               IF (MMAX > 0) THEN
                  W_S(IJK,:MMAX) = -W_S(IPJK,:MMAX)
               ENDIF
            ENDIF
         ENDDO
      ENDDO

! east zy-plane
      I1 = IMAX2
      DO K1 = KSTART3, KEND3
         DO J1 = JSTART3, JEND3
            IF(I1.NE.IEND2) EXIT
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) CALL SET_WALL_BC1 (I1, I1, &
                J1, J1, K1, K1, 0)
         ENDDO
      ENDDO
      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

      END SUBROUTINE SET_WALL_BC

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_WALL_BC1                                            C
!                                                                      C
!  Purpose: Set U, V, and W components for the specified cells by      C
!           copying the same or negative values from near by fluid     C
!           cell                                                       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JAN-92  C
!  Reviewer: M. Syamlal, S. Venkatesan, P. Nicoletti  Date: 29-JAN-92  C
!            W. Rogers                                                 C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: V_g, W_g, U_g, V_s, W_s, U_s in fluid cell    C
!                        adjacent to wall cell                         C
!  Variables modified: V_g, W_g, U_g, V_s, W_s, U_s in wall cell       C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_WALL_BC1(II1, II2, JJ1, JJ2, KK1, KK2, &
                              BC_JJ_PSL)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE bc
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE funits
      USE compar
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Starting and ending I index
      INTEGER, INTENT(IN) :: II1, II2
! Starting and ending J index
      INTEGER, INTENT(IN) :: JJ1, JJ2
! Starting and ending K index
      INTEGER, INTENT(IN) :: KK1, KK2
! Johnson-Jackson boundary condition: 0= no, 1=yes
      INTEGER, INTENT(IN) :: BC_JJ_PSL
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Sign with legal values +1 or -1
      DOUBLE PRECISION :: SIGN0
! Local indices near wall cell
      INTEGER :: I, J, K
      INTEGER :: IJK, IMJK, IJMK, IJKM, IPJK, IJPK, IJKP
      INTEGER :: I1, I2, J1, J2, K1, K2
! Local index for a fluid cell near the wall cell
      INTEGER :: LFLUID
!-----------------------------------------------

! Limit I1, I2 and all to local processor first ghost layer
      I1 = II1
      I2 = II2
      J1 = JJ1
      J2 = JJ2
      K1 = KK1
      K2 = KK2

      IF(I1.LE.IEND2)   I1 = MAX(I1, ISTART2)
      IF(J1.LE.JEND2)   J1 = MAX(J1, JSTART2)
      IF(K1.LE.KEND2)   K1 = MAX(K1, KSTART2)
      IF(I2.GE.ISTART2) I2 = MIN(I2, IEND2)
      IF(J2.GE.JSTART2) J2 = MIN(J2, JEND2)
      IF(K2.GE.KSTART2) K2 = MIN(K2, KEND2)

      DO K = K1, K2
         DO J = J1, J2
            DO I = I1, I2
               IJK = FUNIJK(I,J,K)

               IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells


               IF(NS_WALL_AT(IJK))THEN
                  SIGN0 = -ONE
               ELSE
! note fs_wall occurs for default_wall_at at i==1 if cylindrical and
! xmin=zero
                  SIGN0 = ONE
               ENDIF

               IF (WALL_AT(IJK)) THEN
                  IMJK = IM_OF(IJK)
                  IJMK = JM_OF(IJK)
                  IJKM = KM_OF(IJK)
                  IPJK = IP_OF(IJK)
                  IJPK = JP_OF(IJK)
                  IJKP = KP_OF(IJK)

! Fluid cell at West
                  IF (.NOT.WALL_AT(IMJK)) THEN
                     LFLUID = IMJK
! Wall cell at North
                     IF (WALL_AT(IJPK)) THEN
                        V_G(IJK) = SIGN0*V_G(LFLUID)
                        IF(BC_JJ_PSL==0) CALL EQUAL(V_S,IJK,SIGN0,V_S,LFLUID)
                     ENDIF
! Wall cell at Top
                     IF (WALL_AT(IJKP)) THEN
                        W_G(IJK) = SIGN0*W_G(LFLUID)
                        IF(BC_JJ_PSL==0) CALL EQUAL(W_S,IJK,SIGN0,W_S,LFLUID)
                     ENDIF
                  ENDIF

! Fluid cell at East
                  IF (.NOT.WALL_AT(IPJK)) THEN
                     LFLUID = IPJK
! Wall cell at North
                     IF (WALL_AT(IJPK)) THEN
                        V_G(IJK) = SIGN0*V_G(LFLUID)
                        IF(BC_JJ_PSL==0) CALL EQUAL(V_S,IJK,SIGN0,V_S,LFLUID)
                     ENDIF
! Wall cell at Top
                     IF (WALL_AT(IJKP)) THEN
                        W_G(IJK) = SIGN0*W_G(LFLUID)
                        IF(BC_JJ_PSL==0) CALL EQUAL(W_S,IJK,SIGN0,W_S,LFLUID)
                     ENDIF
                  ENDIF


! Fluid cell at South
                  IF (.NOT.WALL_AT(IJMK)) THEN
                     LFLUID = IJMK
! Wall cell at East
                     IF (WALL_AT(IPJK)) THEN
                        U_G(IJK) = SIGN0*U_G(LFLUID)
                        IF(BC_JJ_PSL==0) CALL EQUAL(U_S,IJK,SIGN0,U_S,LFLUID)
                     ENDIF
! Wall cell at Top
                     IF (WALL_AT(IJKP)) THEN
                        W_G(IJK) = SIGN0*W_G(LFLUID)
                        IF(BC_JJ_PSL==0) CALL EQUAL(W_S,IJK,SIGN0,W_S,LFLUID)
                     ENDIF
                  ENDIF

! Fluid cell at North
                  IF (.NOT.WALL_AT(IJPK)) THEN
                     LFLUID = IJPK
! Wall cell at East
                     IF (WALL_AT(IPJK)) THEN
                        U_G(IJK) = SIGN0*U_G(LFLUID)
                        IF(BC_JJ_PSL==0) CALL EQUAL(U_S,IJK,SIGN0,U_S,LFLUID)
                     ENDIF
! Wall cell at Top
                     IF (WALL_AT(IJKP)) THEN
                        W_G(IJK) = SIGN0*W_G(LFLUID)
                        IF(BC_JJ_PSL==0) CALL EQUAL(W_S,IJK,SIGN0,W_S,LFLUID)
                     ENDIF
                  ENDIF


                  IF (DO_K) THEN
! Fluid cell at Bottom
                     IF (.NOT.WALL_AT(IJKM)) THEN
                        LFLUID = IJKM
! Wall cell at East
                        IF (WALL_AT(IPJK)) THEN
                           U_G(IJK) = SIGN0*U_G(LFLUID)
                           IF (BC_JJ_PSL == 0) CALL EQUAL(U_S, IJK, &
                              SIGN0, U_S, LFLUID)
                        ENDIF
! Wall cell at North
                        IF (WALL_AT(IJPK)) THEN
                           V_G(IJK) = SIGN0*V_G(LFLUID)
                           IF (BC_JJ_PSL == 0) CALL EQUAL(V_S, IJK, &
                              SIGN0, V_S, LFLUID)
                        ENDIF
                     ENDIF

! Fluid cell at Top
                     IF (.NOT.WALL_AT(IJKP)) THEN
                        LFLUID = IJKP
! Wall cell at East
                        IF (WALL_AT(IPJK)) THEN
                           U_G(IJK) = SIGN0*U_G(LFLUID)
                           IF (BC_JJ_PSL == 0) CALL EQUAL(U_S, IJK, &
                              SIGN0, U_S, LFLUID)
                        ENDIF
! Wall cell at North
                        IF (WALL_AT(IJPK)) THEN
                           V_G(IJK) = SIGN0*V_G(LFLUID)
                           IF (BC_JJ_PSL == 0) CALL EQUAL(V_S, IJK, &
                              SIGN0, V_S, LFLUID)
                        ENDIF
                     ENDIF
                  ENDIF   ! end if (do_k)

               ENDIF   ! end if (wall_at(ijk))
            ENDDO   ! end do loop (i = i1, i2)
         ENDDO   ! end do loop (j = j1, j2)
      ENDDO   ! end do loop (k = k1, k2)

      RETURN

      CONTAINS

        INCLUDE 'functions.inc'

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: EQUAL                                                   C
!  Purpose: Loop on the number of solids phases to set a variable      C
!           equal to the value or negative value of another variable   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 29-JAN-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE EQUAL(ARRAY1, IJK1, SIGN0, ARRAY2, IJK2)

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! First array
      DOUBLE PRECISION, INTENT(OUT) :: ARRAY1 (DIMENSION_3, *)
! Second array
      DOUBLE PRECISION, INTENT(IN) :: ARRAY2 (DIMENSION_3, *)
! IJK index for the first array
      INTEGER, INTENT(IN) :: IJK1
! IJK index for the second array
      INTEGER, INTENT(IN) :: IJK2
! Sign to be used when setting ARRAY1.  Legal values
! are + or - 1.0.
      DOUBLE PRECISION, INTENT(IN) :: SIGN0
!-----------------------------------------------
! Local variables
!-----------------------------------------------
!-----------------------------------------------

      IF (MMAX > 0) THEN
         ARRAY1(IJK1,:MMAX) = SIGN0*ARRAY2(IJK2,:MMAX)
      ENDIF

      RETURN
      END SUBROUTINE EQUAL

      END SUBROUTINE SET_WALL_BC1
