!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_SOR                                                 C
!  Purpose: Solve system of linear system using SOR method             C
!           Successive over-relaxation                                 C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-AUG-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE LEQ_SOR(VNAME, VNO, VAR, A_M, B_M, ITMAX, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE leqsol
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! variable name
      CHARACTER(LEN=*), INTENT(IN) :: Vname
! variable number (not really used here; see calling subroutine)
      INTEGER, INTENT(IN) :: VNO
! variable
!     e.g., pp_g, epp, rop_g, rop_s, u_g, u_s, v_g, v_s, w_g,
!     w_s, T_g, T_s, x_g, x_s, Theta_m, scalar, K_Turb_G,
!     e_Turb_G
      DOUBLE PRECISION, INTENT(INOUT) :: Var(DIMENSION_3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: &
                        A_m(DIMENSION_3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: &
                        B_m(DIMENSION_3)
! maximum number of iterations (generally leq_it)
      INTEGER, INTENT(IN) :: ITMAX
! error indicator
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
! OVERRELAXATION FACTOR
      DOUBLE PRECISION, PARAMETER :: OMEGA = 1.0  !1.2
      integer :: iidebug
      parameter( iidebug = 0 )
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Variable
      DOUBLE PRECISION :: Var_tmp(DIMENSION_3)
! Indices
      INTEGER :: IJK
      INTEGER :: ITER

      DOUBLE PRECISION oAm
!-----------------------------------------------

!!$omp parallel do private(IJK,OAM)
      DO IJK = ijkstart3, ijkend3
         IF(.NOT.IS_ON_myPE_owns(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE

         OAM = ONE/A_M(IJK,0)
         A_M(IJK,0) = ONE
         A_M(IJK,-2) = A_M(IJK,-2)*OAM
         A_M(IJK,-1) = A_M(IJK,-1)*OAM
         A_M(IJK,1) = A_M(IJK,1)*OAM
         A_M(IJK,2) = A_M(IJK,2)*OAM
         A_M(IJK,-3) = A_M(IJK,-3)*OAM
         A_M(IJK,3) = A_M(IJK,3)*OAM
         B_M(IJK) = B_M(IJK)*OAM
      ENDDO

      DO ITER = 1, ITMAX
         IF (DO_K) THEN

!!$omp parallel do private(IJK)
            DO IJK = ijkstart3, ijkend3
              IF(.NOT.IS_ON_myPE_owns(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE
              VAR_tmp(IJK) = VAR(IJK) + OMEGA*(B_M(IJK)-&
                 A_M(IJK,-1)*VAR(IM_OF(IJK))-A_M(IJK,1)*VAR(IP_OF(IJK))-&
                 A_M(IJK,-2)*VAR(JM_OF(IJK))-A_M(IJK,2)*VAR(JP_OF(IJK))-&
                 A_M(IJK,-3)*VAR(KM_OF(IJK))-A_M(IJK,3)*VAR(KP_OF(IJK))-&
                 VAR(IJK))
            ENDDO
         ELSE

!!$omp parallel do private(IJK)
           DO IJK = ijkstart3, ijkend3
             IF(.NOT.IS_ON_myPE_owns(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE
             VAR_tmp(IJK) = VAR(IJK) + OMEGA*(B_M(IJK)-&
                  A_M(IJK,-2)*VAR(JM_OF(IJK))-A_M(IJK,2)*VAR(JP_OF(IJK))-&
                  A_M(IJK,-1)*VAR(IM_OF(IJK))-A_M(IJK,1)*VAR(IP_OF(IJK))-&
                  VAR(IJK))
           ENDDO
         ENDIF

      call send_recv(var,2)
      ENDDO

!!$omp parallel do private(IJK)
      DO IJK = ijkstart3, ijkend3
        VAR(IJK) = VAR_tmp(IJK)
      ENDDO

      ITER_TOT(VNO) = ITER


      RETURN
      END SUBROUTINE LEQ_SOR
