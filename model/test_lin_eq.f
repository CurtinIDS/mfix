!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: TEST_LIN_EQ(A_m, LEQIT, LEQMETHOD, LEQSWEEP, LEQTOL, LEQPC, TEST, IER)
!  Purpose: Routine for testing the accuracy of linear equation solver C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 4-JUN-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!  **** See Solve_energy_eq.f for an example  ****                     C
!
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
!
      SUBROUTINE TEST_LIN_EQ(A_M, LEQIT, LEQMETHOD, LEQSWEEP, LEQTOL, LEQPC,  TEST, IER)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE indices
      USE compar
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER TEST  !=0 use the passed A_m; =1 construct a random A_m
      INTEGER IER
      DOUBLE PRECISION, DIMENSION(DIMENSION_3,-3:3) :: A_M
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IJKERR
      INTEGER IJK, IpJK, ImJK, IJpK, IJmK, IJKp, IJKm
      DOUBLE PRECISION, DIMENSION(DIMENSION_3,-3:3) :: Am
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: Bm, X_ACT, X_SOL
      DOUBLE PRECISION :: ERR, ERRMAX, ERRSUM, XSUM
      CHARACTER(LEN=80), DIMENSION(8) :: LINE
!
!                      linear equation solver method and iterations
      INTEGER          LEQMETHOD, LEQIT
      CHARACTER(LEN=4) ::   LEQSWEEP
      DOUBLE PRECISION LEQTOL
      CHARACTER(LEN=4) ::   LEQPC
      REAL  :: Harvest

!  Initialize the random number generator
!
      CALL RANDOM_SEED
!
!  Fill the A and x arrays with random numbers, but ensuring that
!  the matrix is diagonally dominant
!
      DO IJK = IJKSTART3, IJKEND3
         CALL RANDOM_NUMBER(HARVEST)
         X_ACT(IJK) = DBLE(HARVEST) + 1.E-5
         X_SOL(IJK) = 0.0
         IF (TEST == 0) THEN
            Am(IJK,-3) = A_M(IJK,-3)
            Am(IJK,-2) = A_M(IJK,-2)
            Am(IJK,-1) = A_M(IJK,-1)
            Am(IJK,0) = A_M(IJK,0)
            Am(IJK,1) = A_M(IJK,1)
            Am(IJK,2) = A_M(IJK,2)
            Am(IJK,3) = A_M(IJK,3)
         ELSE
            CALL RANDOM_NUMBER(HARVEST)
            Am(IJK,-3) = DBLE(HARVEST)
            CALL RANDOM_NUMBER(HARVEST)
            Am(IJK,-2) = DBLE(HARVEST)
            CALL RANDOM_NUMBER(HARVEST)
            Am(IJK,-1) = DBLE(HARVEST)
            CALL RANDOM_NUMBER(HARVEST)
            Am(IJK,0) = -DBLE(MAX(HARVEST,0.1))*70.
            CALL RANDOM_NUMBER(HARVEST)
            Am(IJK,1) = DBLE(HARVEST)
            CALL RANDOM_NUMBER(HARVEST)
            Am(IJK,2) = DBLE(HARVEST)
            CALL RANDOM_NUMBER(HARVEST)
            Am(IJK,3) = DBLE(HARVEST)
         ENDIF
      END DO


!
!!!$omp  parallel do private( IJK, IpJK, ImJK, IJpK, IJmK, IJKp, IJKm)
      DO IJK = ijkstart3, ijkend3

            ImJK = IM_OF(IJK)
            IJmK = JM_OF(IJK)
            IpJK = IP_OF(IJK)
            IJpK = JP_OF(IJK)
            Bm(IJK) = Am(IJK,0)*X_ACT(IJK)
            IF(I_OF(IJK) > 1) Bm(IJK) = Bm(IJK) + Am(IJK,west)*X_ACT(ImJK)
            IF(I_OF(IJK) < IMAX2) Bm(IJK) = Bm(IJK) +Am(IJK,east)*X_ACT(IpJK)
            IF(J_OF(IJK) > 1) Bm(IJK) = Bm(IJK) + Am(IJK,south)*X_ACT(IJmK)
            IF(J_OF(IJK) < JMAX2) Bm(IJK) = Bm(IJK) + Am(IJK,north)*X_ACT(IJpK)
            IF (DO_K) THEN
               IJKm = KM_OF(IJK)
               IJKp = KP_OF(IJK)
               IF(K_OF(IJK) > 1) Bm(IJK) = Bm(IJK) + Am(IJK,bottom)*X_ACT(IJKm)
               IF(K_OF(IJK) < KMAX2) Bm(IJK) = Bm(IJK) + Am(IJK,top)*X_ACT(IJKp)
            ENDIF
      END DO

!
!  Solve the linear equation
!
      CALL SOLVE_LIN_EQ ('Test', 1, X_SOL, Am, Bm, 0, LEQIT, LEQMETHOD, LEQSWEEP, LEQTOL, LEQPC, IER)
!
!  Check the solution
!
      ERRSUM = 0.0
      XSUM = 0.0
      ERRMAX = 0.0
      IJKERR = 0
      DO IJK = ijkstart3, ijkend3
         IF (X_ACT(IJK) /= 0.0) THEN
            ERR = ABS(X_SOL(IJK)-X_ACT(IJK))/X_ACT(IJK)
         ELSE IF (X_SOL(IJK) == 0.0) THEN
            ERR = 0.0
         ELSE
            ERR = 1.0D32
         ENDIF
         IF (ERR > ERRMAX) THEN
            ERRMAX = ERR
            IJKERR = IJK
         ENDIF
         ERRSUM = ERRSUM + ABS(X_SOL(IJK)-X_ACT(IJK))
         XSUM = XSUM + ABS(X_ACT(IJK))
!        print *, ijk, i_of(ijk), j_of(ijk), nint(err * 100.)
      END DO
      IF (XSUM /= 0.0) THEN
         ERR = ERRSUM/XSUM
      ELSE IF (ERRSUM == 0.0) THEN
         ERR = 0.0
      ELSE
         ERR = 1.0D32
      ENDIF
!
      IF (ERR < LEQTOL) THEN
         IER = 0
         LINE(1) = 'Message: Lin equation solution satisfies tolerance.'
      ELSE
         IER = 1
         LINE(1) = 'Error: Lin equation solution does not satisfy tolerance!'
      ENDIF
!
      WRITE (LINE(2), *) 'Average normalized error = ', ERR
      WRITE (LINE(3), *) 'Max normalized error = ', ERRMAX
      WRITE (LINE(4), *) 'Location of max error = ', IJKERR
      WRITE (LINE(5), *) 'Xa and Xs @ max error = ', X_ACT(IJKERR), X_SOL(IJKERR)
      WRITE (LINE(6), *) 'Sample values of actual (Xa) and solution (Xs):'
      WRITE (LINE(7), '(A,G12.5, A, I6, A, G12.5, A, I6, A, G12.5)') 'Xa(1)=', &
         X_ACT(1), '  Xa(', IJKMAX2/2, ')=', X_ACT(IJKMAX2/2), '  Xa(', IJKMAX2, ')=', &
         X_ACT(IJKMAX2)
      WRITE (LINE(8), '(A,G12.5, A, I6, A, G12.5, A, I6, A, G12.5)') 'Xs(1)=', &
         X_SOL(1), '  Xs(', IJKMAX2/2, ')=', X_SOL(IJKMAX2/2), '  Xs(', IJKMAX2&
         , ')=', X_SOL(IJKMAX2)
!
      CALL WRITE_ERROR ('TEST_LIN_EQ', LINE, 8)
!
      RETURN
      END SUBROUTINE TEST_LIN_EQ




