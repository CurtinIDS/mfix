!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: ADJUST_EPS                                              C
!  Purpose: Eliminate the solids phases that occupy only very small    C
!           fractions of the computational cell volume                 C
!                                                                      C
!  Author: M. Syamlal                                 Date: 28-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 29-JAN-92  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE ADJUST_EPS

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1, only: zero
! minimum value of solids volume fraction tracked
      USE toleranc, only: zero_ep_s
! x,y,z-components of solids velocity
      USE fldvar, only: u_s, v_s, w_s
! gas void fracition, bulk density and density
      USE fldvar, only: ep_g, rop_g, ro_g
! solids phase particle bulk density and material density
      USE fldvar, only: rop_s, ro_s
! kinetic theories
      USE run, only: kt_type_enum
      USE run, only: ghd_2007
! number of solids phases
      USE physprop, only: mmax, smax

! needed for function.inc and other quantities
      USE geometry
      USE indices
      USE compar
! for sendrecv calls
      USE sendrecv
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IJK
! Solids phase
      INTEGER :: M
! Sum of (very small) solids volume fractions that are set to zero
      DOUBLE PRECISION :: EPSUM
! Sum of solids volume fractions
      DOUBLE PRECISION :: epsMix, epSolid
!-----------------------------------------------


      DO K = Kstart1, Kend1
         DO J = Jstart1, Jend1
            DO I = Istart1, Iend1

               IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells

               IJK = FUNIJK(I,J,K)
               IF (FLUID_AT(IJK)) THEN
! initialize
                  EPSUM = ZERO
                  epsMix = ZERO
                  DO M = 1, SMAX
! why not use function ep_s?
                     epSolid = ROP_S(IJK,M)/RO_S(IJK,M)
                     epsMix = epsMix +  epSolid

                     IF (epSolid < ZERO_EP_S) THEN
! Remove solids in very small quantities and set solids velocity to zero
! if there is outflow from the present cell.
                        EPSUM = EPSUM + epSolid

                        IF(KT_TYPE_ENUM == GHD_2007) &
                           ROP_S(IJK,MMAX) = ROP_S(IJK,MMAX) - ROP_S(IJK,M)

                        ROP_S(IJK,M) = ZERO
                        U_S(IJK,M) = MIN(U_S(IJK,M),ZERO)
                        V_S(IJK,M) = MIN(V_S(IJK,M),ZERO)
                        W_S(IJK,M) = MIN(W_S(IJK,M),ZERO)
                        U_S(IM_OF(IJK),M) = MAX(U_S(IM_OF(IJK),M),ZERO)
                        V_S(JM_OF(IJK),M) = MAX(V_S(JM_OF(IJK),M),ZERO)
                        W_S(KM_OF(IJK),M) = MAX(W_S(KM_OF(IJK),M),ZERO)
                     ENDIF
                  ENDDO

                  epsMix = epsMix - EPSUM
                  IF(KT_TYPE_ENUM == GHD_2007 .AND. epsMix < ZERO_EP_S) THEN
                    U_S(IJK,MMAX) = MIN(U_S(IJK,MMAX),ZERO)
                    V_S(IJK,MMAX) = MIN(V_S(IJK,MMAX),ZERO)
                    W_S(IJK,MMAX) = MIN(W_S(IJK,MMAX),ZERO)
                    U_S(IM_OF(IJK),MMAX) = MAX(U_S(IM_OF(IJK),MMAX),ZERO)
                    V_S(JM_OF(IJK),MMAX) = MAX(V_S(JM_OF(IJK),MMAX),ZERO)
                    W_S(KM_OF(IJK),MMAX) = MAX(W_S(KM_OF(IJK),MMAX),ZERO)
                  ENDIF

                  EP_G(IJK) = EP_G(IJK) + EPSUM
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

! Communicate field variables calculated in the do i,j,k loop
      call send_recv(ROP_S,2)
      call send_recv(U_S,2)
      call send_recv(V_S,2)
      call send_recv(W_S,2)
      call send_recv(EP_G,2)
      call send_recv(ROP_G,2)

      RETURN
      END SUBROUTINE ADJUST_EPS
