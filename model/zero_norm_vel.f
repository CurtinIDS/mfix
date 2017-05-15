!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: ZERO_NORM_VEL                                           C
!  Purpose: Set the velocity component normal to a wall to zero        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 14-MAY-92  C
!  Reviewer: S. Venkatesan                            Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE ZERO_NORM_VEL

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE geometry
      USE physprop
      USE fldvar
      USE indices
      USE is
      USE compar
      USE discretelement
      USE mfix_pic
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
!
      INTEGER :: ISV
! Indicies
      INTEGER :: IJK, IMJK, IJMK, IJKM
      INTEGER :: M

!!$omp  parallel do private( IMJK, IJMK, IJKM)

      DO IJK = ijkstart3, ijkend3

         IF (.NOT.WALL_AT(IJK)) THEN
            IF (IP_AT_E(IJK)) U_G(IJK) = ZERO
            IF (IP_AT_N(IJK)) V_G(IJK) = ZERO
            IF (IP_AT_T(IJK)) W_G(IJK) = ZERO
         ELSE
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)
            U_G(IJK) = ZERO
            V_G(IJK) = ZERO
            W_G(IJK) = ZERO
            IF (.NOT.(CYCLIC_AT(IJK) .AND. (I_OF(IJK)==IMAX2 .OR. &
                I_OF(IJK)==IMAX3))) U_G(IMJK) = ZERO
            IF (.NOT.(CYCLIC_AT(IJK) .AND. (J_OF(IJK)==JMAX2 .OR. &
                J_OF(IJK)==JMAX3))) V_G(IJMK) = ZERO
            IF (.NOT.(CYCLIC_AT(IJK) .AND. (K_OF(IJK)==KMAX2 .OR. &
                K_OF(IJK)==KMAX3))) W_G(IJKM) = ZERO
         ENDIF
      ENDDO   ! end do (ijk=ijkstart3,ijkend3)


      IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID .OR.MPPIC) THEN
         DO M = 1, MMAX
!!$omp  parallel do private( ISV,  IMJK, IJMK, IJKM)
            DO IJK = ijkstart3, ijkend3

               IF (.NOT.WALL_AT(IJK)) THEN
                  IF (IP_AT_E(IJK)) THEN
                     U_S(IJK,M) = ZERO
                  ELSEIF (SIP_AT_E(IJK)) THEN
                     ISV = FLAG_E(IJK) - 1000
                     U_S(IJK,M) = IS_VEL_S(ISV,M)
                  ENDIF
                  IF (IP_AT_N(IJK)) THEN
                     V_S(IJK,M) = ZERO
                  ELSEIF (SIP_AT_N(IJK)) THEN
                     ISV = FLAG_N(IJK) - 1000
                     V_S(IJK,M) = IS_VEL_S(ISV,M)
                  ENDIF
                  IF (IP_AT_T(IJK)) THEN
                     W_S(IJK,M) = ZERO
                  ELSEIF (SIP_AT_T(IJK)) THEN
                     ISV = FLAG_T(IJK) - 1000
                     W_S(IJK,M) = IS_VEL_S(ISV,M)
                  ENDIF
               ELSE
                  IMJK = IM_OF(IJK)
                  IJMK = JM_OF(IJK)
                  IJKM = KM_OF(IJK)
                  U_S(IJK,M) = ZERO
                  V_S(IJK,M) = ZERO
                  W_S(IJK,M) = ZERO
                  IF (.NOT.(CYCLIC_AT(IJK) .AND. (I_OF(IJK)==IMAX2 .OR. &
                      I_OF(IJK)==IMAX3))) U_S(IMJK,M) = ZERO
                  IF (.NOT.(CYCLIC_AT(IJK) .AND. (J_OF(IJK)==JMAX2 .OR. &
                      J_OF(IJK)==JMAX3))) V_S(IJMK,M) = ZERO
                  IF (.NOT.(CYCLIC_AT(IJK) .AND. (K_OF(IJK)==KMAX2 .OR. &
                      K_OF(IJK)==KMAX3))) W_S(IJKM,M) = ZERO
               ENDIF
            ENDDO   ! end do (ijk=ijkstart3,ijkend3)
         ENDDO   ! end do (m=1,mmax)
      ENDIF   ! endif (.not.discrete_element)

      RETURN

      END SUBROUTINE ZERO_NORM_VEL



