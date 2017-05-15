!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CHECK_CONVERGENCE                                       C
!  Purpose: Monitor convergence                                        C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 8-JUL-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CHECK_CONVERGENCE(NIT, errorpercent, MUSTIT)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE geometry
      USE indices
      USE mpi_utility
      USE param
      USE param1
      USE physprop
      USE residual
      USE run
      USE scalars, only :NScalar
      USE toleranc
      USE utilities, ONLY: check_vel_bound
      IMPLICIT NONE
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
! Maximum % error allowed in fluid continuity
!      DOUBLE PRECISION, PARAMETER :: MaxErrorPercent = 1.0E-6
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Iteration number
      INTEGER, INTENT(IN) :: NIT
! %error in fluid mass balance
      DOUBLE PRECISION, INTENT(IN) :: errorpercent
! value tells whether to iterate (1) or not (0).
      INTEGER, INTENT(INOUT) :: MUSTIT
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! sum of residuals
      DOUBLE PRECISION :: SUM, SUM_T, SUM_X, SUM_Th
! max of residuals
      DOUBLE PRECISION :: maxres
! index
      INTEGER :: L, M, maxL, maxM, NN, maxN
! to indicate undefined residual in species eq at the
! beginning of iterations
      LOGICAL :: NO_RESID
!-----------------------------------------------

! sum the residuals from correction equation (pressure and/or
! solids), continuity equations (gas and/or solids) and momentum
! equations (gas and solids)

! add pressure correction residual
!      if(abs(errorpercent) > MaxErrorPercent)then
        SUM = RESID(RESID_P,0)
!      else
!        SUM = zero
!      endif


! add solids correction residual
      IF(MMAX > 0) SUM = SUM + RESID(RESID_P,1)

! add continuity equation residuals
      DO M = 0, MMAX
         SUM = SUM + RESID(RESID_RO,M)
      ENDDO
! add momentum equation residuals
      DO M = 0, MMAX
         SUM = SUM + RESID(RESID_U,M)
      ENDDO
      DO M = 0, MMAX
         SUM = SUM + RESID(RESID_V,M)
      ENDDO
      IF (DO_K) THEN
         DO M = 0, MMAX
            SUM = SUM + RESID(RESID_W,M)
         ENDDO
      ENDIF
!      call global_all_sum(SUM)


! sum the granular energy equation residuals
      SUM_Th = zero
      IF (GRANULAR_ENERGY) THEN
         DO M = 1, MMAX
            SUM_Th = SUM_Th + RESID(RESID_TH,M)
         ENDDO
      ENDIF


! sum the energy equation residuals
      SUM_T = ZERO
      IF (ENERGY_EQ) THEN
         DO M = 0, MMAX
            SUM_T = SUM_T + RESID(RESID_T,M)
         END DO
      ENDIF
!      call global_all_sum(SUM_T)

! sum the species equation residuals
      SUM_X = ZERO
      NO_RESID = .FALSE.
      DO M = 0, MMAX
         IF (SPECIES_EQ(M)) THEN
            DO NN = 1, NMAX(M)
               IF (RESID(RESID_X+(NN-1),M) == UNDEFINED) NO_RESID = .TRUE.
               SUM_X = SUM_X + RESID(RESID_X+(NN-1),M)
            ENDDO
         ENDIF
      ENDDO
!      call global_all_sum(SUM_X)
      IF (NO_RESID) SUM_X = TOL_RESID_X + ONE


! find the variable with maximum residual
      IF (RESID_INDEX(MAX_RESID_INDEX,1) == UNDEFINED_I) THEN
         MAXRES = ZERO
         DO L = 1, NRESID
            DO M = 0, MMAX
               IF (RESID(L,M) >= MAXRES) THEN
                  MAXRES = RESID(L,M)
                  MAXL = L
                  MAXM = M
                  IF (L >= RESID_X) THEN
                     MAXN = L - RESID_X + 1
                  ELSE
                     MAXN = UNDEFINED_I
                  ENDIF
               ENDIF
            END DO
         END DO
         IF (MAXN == UNDEFINED_I) THEN
            WRITE (RESID_STRING(MAX_RESID_INDEX), '(A1,I1)') RESID_PREFIX(MAXL)&
               , MAXM
         ELSE
            WRITE (RESID_STRING(MAX_RESID_INDEX), '(A1,I1,I2.0)') 'X', MAXM, &
               MAXN
         ENDIF
      ENDIF

      IF (GROUP_RESID) THEN
         RESID_GRP(HYDRO_GRP) = SUM
         IF(GRANULAR_ENERGY) RESID_GRP(THETA_GRP) = SUM_TH
         IF(ENERGY_EQ) RESID_GRP(ENERGY_GRP) = SUM_T
         IF(ANY_SPECIES_EQ) RESID_GRP(SPECIES_GRP) = SUM_X
         IF(NScalar > 0) RESID_GRP(SCALAR_GRP) = RESID(RESID_sc,0)
         IF(K_EPSILON) RESID_GRP(KE_GRP) = RESID(RESID_ke,0)
      ENDIF

! Every 5 iterations detect whether the run is stalled by checking
! that the total residual has decreased.
      IF(DETECT_STALL .AND. MOD(NIT,5) == 0) THEN
         IF(NIT > 10) THEN
            IF(SUM5_RESID <= SUM) THEN
! The run is stalled. Reduce the time step.
               IF(.NOT.PERSISTENT_MODE) THEN
                  MUSTIT = 2
                  RETURN
! Forces the max number of iterations for DT=DT_MIN
               ELSEIF(DT > DT_MIN) THEN
                  MUSTIT = 1
                  RETURN
               ENDIF
            ENDIF
         ENDIF
         SUM5_RESID = SUM
      ENDIF

! Require at least two iterations.
      IF(NIT == 1) THEN
         MUSTIT = 1
         RETURN
      ENDIF

! total residual
      IF(SUM<=TOL_RESID .AND. SUM_T<=TOL_RESID_T .AND. &
         RESID(RESID_sc,0)<=TOL_RESID_Scalar .AND. SUM_X<=TOL_RESID_X &
        .AND. RESID(RESID_ke,0)<=TOL_RESID_K_Epsilon &
        .AND. SUM_Th <=TOL_RESID_Th)THEN
         MUSTIT = 0                              !converged
      ELSEIF (SUM>=TOL_DIVERGE .OR. SUM_T>=TOL_DIVERGE .OR.&
              RESID(RESID_sc,0)>= TOL_DIVERGE .OR. SUM_X>=TOL_DIVERGE&
              .OR. RESID(RESID_ke,0)>= TOL_DIVERGE &
              .OR. SUM_Th >= TOL_DIVERGE ) THEN
         IF (NIT /= 1) THEN
            MUSTIT = 2                           !diverged
         ELSE
            MUSTIT = 1                           !not converged
         ENDIF
      ELSE
         MUSTIT = 1                              !not converged
      ENDIF


! to check upper bound (speed of sound) limit for gas and
! solids velocity components
! only check velocity if any of the momentum equations are solved
      DO M = 0,MMAX
         IF (MOMENTUM_X_EQ(M) .OR. MOMENTUM_Y_EQ(M) .OR. &
             MOMENTUM_Z_EQ(M)) THEN
            IF(CHECK_VEL_BOUND()) MUSTIT = 2     !divergence
            EXIT  ! only need to call check_vel_bound once
         ENDIF
      ENDDO


      RETURN
      END SUBROUTINE CHECK_CONVERGENCE


