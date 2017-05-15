!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: CHECK_DATA_30                                          !
!  Author: M. Syamlal                                 Date: 27-OCT-92  !
!                                                                      !
!  Purpose: Check whether the sum of reaction rates is zero and the    !
!  sum of mass fractions is 1.0 and EP_g >= EP_Star.                   !
!           and EP_g >= EP_star. Set miscellaneous constants           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_DATA_30

! Global variables: (common to sub-functions)
!---------------------------------------------------------------------//
      use run, only: ENERGY_EQ, SPECIES_EQ
      use compar, only: ISTART2, IEND2
      use compar, only: JSTART2, JEND2
      use compar, only: KSTART2, KEND2
      use physprop, only: SMAX, NMAX
      use param1, only: ZERO, ONE
      use rxns, only: USE_RRATES
      use mms, only: USE_MMS
      use param, only: DIMENSION_M

      use error_manager

      IMPLICIT NONE

! Local variables:
!---------------------------------------------------------------------//
! Flag for error message headers.
      INTEGER :: ERR_TAG
! Flags to report stats on species equations.
      LOGICAL :: REPORT_SPECIES(0:DIMENSION_M)
!......................................................................!


! Check physical properties in inflow/outflow cells.
      IF(.NOT.USE_MMS) CALL CHECK_FLOW_CELL_PROPS
! Verify physical values for field variables.
      CALL CHECK_PHYSICAL_BOUNDS
! Generate species report if needed.
      IF(ANY(REPORT_SPECIES)) CALL REPORT_SPECIES_STATS
! Check interphase mass transfer.
      IF(USE_RRATES) CALL CHECK_RXN_MASS_BALANCE

      RETURN

      CONTAINS

!----------------------------------------------------------------------!
!  Module name: CHECK_FLOW_PROPS                                       !
!  Author: M. Syamlal                                 Date: 27-OCT-92  !
!                                                                      !
!  Purpose: Verify that inflow/outflow cells do not contain physical   !
!  properties for specified variables.                                 !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_FLOW_CELL_PROPS

! Global variables:
!---------------------------------------------------------------------//
      use visc_s, only: MU_S, LAMBDA_S
      use visc_g, only: MU_GT, LAMBDA_GT
      use physprop, only: K_G, K_S
      use physprop, only: DIF_s, DIF_g

      use functions, only: FUNIJK
      use functions, only: FLOW_AT
      use compar, only: DEAD_CELL_AT

      IMPLICIT NONE

! Local variables:
!---------------------------------------------------------------------//
! Loop counters
      INTEGER :: I, J, K, IJK
      INTEGER :: M, N
! Integer error flag.
      INTEGER :: IER
!......................................................................!

      CALL INIT_ERR_MSG("CHECK_FLOW_CELL_PROPS")

      IER = 0
      ERR_TAG=2000

      DO K = KSTART2, KEND2
      DO J = JSTART2, JEND2
      DO I = ISTART2, IEND2

         IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells

         IJK = FUNIJK(I,J,K)

         IF(FLOW_AT(IJK)) THEN

! Turbulent viscosity of fluid phase.
            IF(MU_gt(IJK) /= ZERO) CALL REPORT_ERROR                   &
               (IER, I, J, K, MU_GT(IJK), '/=', ZERO, 'MU_GT')

! Granular second coefficient of viscosity.
            IF(LAMBDA_gt(IJK) /= ZERO) CALL REPORT_ERROR               &
               (IER, I, J, K, LAMBDA_gt(IJK), '/=', ZERO, 'LAMBDA_GT')
! Gas conductivity.
            IF(K_g(IJK) /= ZERO) CALL REPORT_ERROR                     &
               (IER, I, J, K, K_g(IJK),'/=',ZERO,'K_G')
! Diffusivity of gas species N.
            DO N = 1, NMAX(0)
               IF(DIF_g(IJK, N) /= ZERO) CALL REPORT_ERROR             &
                  (IER, I, J, K, DIF_g(IJK,N),'/=',ZERO,'DIF_G',N)
            ENDDO

            DO M = 1, SMAX
! Granular first coefficient of (shear) viscosity.
               IF(MU_s(IJK, M) /= ZERO) CALL REPORT_ERROR              &
                  (IER, I, J, K, MU_s(IJK,M),'/=',ZERO,'MU_S',M)
! Granular second coefficient of viscosity.
               IF(LAMBDA_s(IJK, M) /= ZERO) CALL REPORT_ERROR          &
                  (IER, I, J, K, LAMBDA_s(IJK,M),'/=',ZERO,'LAMBDA_S',M)
! Solids thermal conductivity.
               IF(K_s(IJK, M) /= ZERO) CALL REPORT_ERROR               &
                  (IER, I, J, K, K_s(IJK,M),'/=',ZERO,'K_S', M)
! Diffusivity of solids phase M, species N.
               DO N = 1, NMAX(M)
                  IF(DIF_s(IJK,M,N) /= ZERO) CALL REPORT_ERROR         &
                     (IER, I,J,K, DIF_s(IJK,M,N),'/=',ZERO,'DIF_S',M,N)
               ENDDO
            ENDDO
         ENDIF

      ENDDO
      ENDDO
      ENDDO

      IF(IER /= 0) THEN
         WRITE(ERR_MSG,"('End of Report.')")
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., ABORT=.TRUE.)
! Close DMP logs when running in interactive mode.
         CALL CLOSE_PE_LOG
      ENDIF

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_FLOW_CELL_PROPS



!----------------------------------------------------------------------!
!                                                                      !
!  Module name: CHECK_PHYSICAL_BOUNDS                                  !
!  Author: M. Syamlal                                 Date: 27-OCT-92  !
!                                                                      !
!  Purpose: Verify that fluid cells have physical values for the       !
!  specified variables.                                                !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_PHYSICAL_BOUNDS

! Global variables:
!---------------------------------------------------------------------//
      use toleranc, only: TMIN, TMAX, TOL_COM
      use fldvar, only: X_G, X_S
      use fldvar, only: T_G, T_S
      use fldvar, only: ROP_s
      use rxns, only: RoX_GC, RoX_SC
      use rxns, only: R_GP, R_SP
      use visc_s, only: MU_S
      use physprop, only: K_G, K_S
      use physprop, only: C_PG, C_PS
      use physprop, only: DIF_s, DIF_g
      use physprop, only: MU_G
      use physprop, only: MW_MIX_G

      use mpi_utility, only: GLOBAL_ALL_SUM
      use mpi_utility, only: GLOBAL_ALL_OR

      use functions, only: FUNIJK
      use functions, only: WALL_AT
      use compar, only: DEAD_CELL_AT

      IMPLICIT NONE

! Local variables:
!---------------------------------------------------------------------//
      INTEGER :: I, J, K, IJK
      INTEGER :: M, N
! Integer error flag.
      INTEGER :: IER

      CALL INIT_ERR_MSG("CHECK_PHYSICAL_BOUNDS")

      IER = 0
      ERR_TAG = 3000
      REPORT_SPECIES = .FALSE.

      DO K = KSTART2, KEND2
      DO J = JSTART2, JEND2
      DO I = ISTART2, IEND2
         IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
         IJK = FUNIJK(I,J,K)
         IF (.NOT.WALL_AT(IJK)) THEN

! Gas viscosity must be non-negative.
            IF(MU_G(IJK) < ZERO) CALL REPORT_ERROR                     &
               (IER, I, J, K, MU_G(IJK), '<', ZERO, 'MU_G')

! Mixture molecular weight must be positive.
            IF(MW_MIX_G(IJK) <= ZERO) CALL REPORT_ERROR                &
               (IER, I, J, K,MW_MIX_G(IJK), '<=', ZERO, 'MW_MIX_G')

! Verify thermodynamic properties when solving energy equations:
            IF(ENERGY_EQ) THEN
! Gas conductivity must be non-negative.
               IF(K_G(IJK) < ZERO) CALL REPORT_ERROR                   &
                  (IER, I, J, K, K_G(IJK), '<', ZERO, 'K_G')
! Gas phase specific heat must be positive.
               IF(C_PG(IJK) <= ZERO) CALL REPORT_ERROR                 &
                  (IER, I, J, K, C_PG(IJK), '<=', ZERO, 'C_PG')
! Verify that the gas phase temperature is within the bounds.
               IF(T_G(IJK) <= TMIN) CALL REPORT_ERROR                  &
                  (IER, I, J, K, T_G(IJK), '<=', TMIN, 'T_G')
               IF(T_G(IJK) >= TMAX) CALL REPORT_ERROR                  &
                  (IER, I, J, K, T_G(IJK), '>=', TMAX, 'T_G')

! Diffusivity of gas species N must be non-negative.
               DO N = 1, NMAX(0)
                  IF(DIF_g(IJK, N) < ZERO) CALL REPORT_ERROR           &
                     (IER, I,J,K, DIF_g(IJK,N),'<',ZERO, 'DIF_G',N)
                  IF(X_g(IJK,N) < ZERO) CALL REPORT_ERROR              &
                     (IER, I, J, K, X_g(IJK,N), '<', ZERO, 'X_G',N)
                  IF(X_g(IJK,N) > ONE) CALL REPORT_ERROR               &
                     (IER, I, J, K, X_g(IJK,N), '>', ONE, 'X_G',N)
               ENDDO
            ENDIF

! Verify that the gas phase mass fractons sum to one.
            IF(SPECIES_EQ(0)) THEN
               IF(ABS(ONE - sum(X_G(IJK,1:NMAX(0)))) > TOL_COM) &
                  REPORT_SPECIES(0) = .TRUE.

! Verify that the rates of formation and consumption adhear to the
! expected coding restraints. (non-negative)
               IF(USE_RRATES) THEN
                  DO N = 1, NMAX(0)
                     IF(R_GP(IJK,N) < ZERO) CALL REPORT_ERROR          &
                       (IER,I,J,K, R_GP(IJK,N),'<',ZERO,'R_GP',N)
                     IF(ROX_GC(IJK,N) < ZERO) CALL REPORT_ERROR        &
                       (IER,I,J,K,RoX_GC(IJK,N),'<',ZERO,'RoX_GC',N)
                  ENDDO
               ENDIF

            ENDIF

            DO M = 1, SMAX

! Solids viscosity should be non-negativel.
               IF(MU_S(IJK,M) < ZERO) CALL REPORT_ERROR                &
                  (IER, I, J, K, MU_S(IJK,M), '<', ZERO, 'MU_S',M)

               IF(ENERGY_EQ) THEN

! Thermal conductivity must be non-negative.
                  IF(K_S(IJK,M) < ZERO) CALL REPORT_ERROR              &
                     (IER, I, J, K, K_S(IJK,M), '<', ZERO, 'K_S',M)

! Solids specific heat must be positive.
                  IF(C_PS(IJK,M) <= ZERO) CALL REPORT_ERROR            &
                     (IER, I, J, K, C_PS(IJK,M), '<=', ZERO, 'C_PS',M)

! Verify that the solids temperature is within the required bounds.
                  IF(T_S(IJK,M) <= TMIN) CALL REPORT_ERROR             &
                     (IER, I, J, K, T_S(IJK,M), '<=', TMIN, 'T_S',M)
                  IF(T_S(IJK,M) >= TMAX) CALL REPORT_ERROR             &
                     (IER, I, J, K, T_S(IJK,M), '>=', TMAX, 'T_S',M)

               ENDIF

               DO N = 1, NMAX(M)
                  IF(DIF_s(IJK, M, N) < ZERO) CALL REPORT_ERROR        &
                     (IER, I, J, K, DIF_s(IJK,M,N),'<',ZERO,'DIF_S',M,N)
               ENDDO

! Sum of solids mass fractions should be one
               IF(SPECIES_EQ(M)) THEN
                  IF(ROP_S(IJK,M) /= ZERO) THEN
                     IF(ABS(ONE-sum(X_S(IJK,M,1:NMAX(M)))) > TOL_COM) &
                        REPORT_SPECIES(M) = .TRUE.
                  ENDIF

                  DO N = 1, NMAX(M)
                     IF(X_s(IJK, M, N) < ZERO) CALL REPORT_ERROR       &
                        (IER, I,J,K, X_s(IJK,M,N),'<',ZERO, 'X_S',M,N)
                     IF(X_s(IJK, M, N) > ONE) CALL REPORT_ERROR        &
                        (IER, I,J,K, X_s(IJK,M,N),'>',ONE, 'X_S',M,N)
                  ENDDO
               ENDIF

! Verify that the rates of formation and consumption adhear to the
! expected coding restraints. (non-negative)
               IF(SPECIES_EQ(M) .AND. USE_RRATES) THEN
                  DO N = 1, NMAX(0)
                     IF(R_SP(IJK,M,N) < ZERO) CALL REPORT_ERROR        &
                       (IER,I,J,K,R_SP(IJK,M,N),'<',ZERO,'R_SP',M,N)
                     IF(ROX_GC(IJK,N) < ZERO) CALL REPORT_ERROR        &
                       (IER,I,J,K,RoX_SC(IJK,M,N),'<',ZERO,'RoX_SC',M,N)
                  ENDDO
               ENDIF
            ENDDO

         ENDIF ! IF(.NOT.WALL_AT(IJK))
      ENDDO ! DO I = ISTART2, IEND2
      ENDDO ! DO J = JSTART2, JEND2
      ENDDO ! DO K = KSTART2, KEND2

      CALL GLOBAL_ALL_SUM(IER)

      IF(IER /= 0) THEN
         WRITE(ERR_MSG,"('End of Report.')")
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., ABORT=.TRUE.)
         CALL CLOSE_PE_LOG
      ENDIF

      CALL GLOBAL_ALL_OR(REPORT_SPECIES)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_PHYSICAL_BOUNDS



!----------------------------------------------------------------------!
!                                                                      !
!  Module name: REPORT_SPECIES_STATS                                   !
!  Author: M. Syamlal                                 Date: 27-OCT-92  !
!                                                                      !
!  Purpose: Collect information on the sums of species mass fractions. !
!  This routine is call as-needed.                                     !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE REPORT_SPECIES_STATS

! Global variables:
!---------------------------------------------------------------------//
      use fldvar, only: X_G, X_S
      use fldvar, only: ROP_S
      use run, only: TIME
      use mpi_utility, only: GLOBAL_ALL_SUM
      use mpi_utility, only: GLOBAL_ALL_MAX
      use mpi_utility, only: GLOBAL_ALL_MIN

      use functions, only: FUNIJK
      use functions, only: WALL_AT

      use param1, only: UNDEFINED
      use compar, only: DEAD_CELL_AT

      IMPLICIT NONE

! Local variables:
!---------------------------------------------------------------------//
      INTEGER :: I, J, K, IJK, L
      INTEGER :: M

      INTEGER :: COUNT(0:DIMENSION_M, 9)
      INTEGER :: MIN_LOC(0:DIMENSION_M)
      INTEGER :: MAX_LOC(0:DIMENSION_M)

      DOUBLE PRECISION :: MIN_VAL(0:DIMENSION_M)
      DOUBLE PRECISION :: MAX_VAL(0:DIMENSION_M)

      DOUBLE PRECISION :: lSUM
!......................................................................!

      CALL INIT_ERR_MSG('REPORT_SPECIES_STATS')

      WRITE(ERR_MSG, 4000) TIME
      CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)

 4000 FORMAT('Time = ',G12.5,/'Warning: The sum of mass fractions ',   &
      'is not equal to one.')

! Initialize the counters
      COUNT = 0
      MAX_VAL = -UNDEFINED
      MIN_VAL =  UNDEFINED

! Collect stats on species across the domain.
      DO K = KSTART2, KEND2
      DO J = JSTART2, JEND2
      DO I = ISTART2, IEND2
         IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
         IJK = FUNIJK(I,J,K)
         IF (.NOT.WALL_AT(IJK)) THEN

            IF(SPECIES_EQ(0)) THEN
               lSUM = sum(X_G(IJK,1:NMAX(0)))
               IF (lSUM < MIN_VAL(0)) THEN
                  MIN_VAL(0) = lSUM
                  MIN_LOC(0) = IJK
               ENDIF
               IF (lSUM > MAX_VAL(0)) THEN
                  MAX_VAL(0) = lSUM
                  MAX_LOC(0) = IJK
               ENDIF
               IF (lSUM < 0.9) THEN
                  COUNT(0,1) = COUNT(0,1) + 1 ! < 0.9pP
               ELSE IF (lSUM < 0.99) THEN
                  COUNT(0,2) = COUNT(0,2) + 1 ! 0.9    - 0.99
               ELSE IF (lSUM < 0.999) THEN
                  COUNT(0,3) = COUNT(0,3) + 1 ! 0.99   - 0.999
               ELSE IF (lSUM < 0.9999) THEN
                  COUNT(0,4) = COUNT(0,4) + 1 ! 0.999  - 0.9999
               ELSE IF (lSUM < 1.0001) THEN
                  COUNT(0,5) = COUNT(0,5) + 1 ! 0.9999 - 1.0001
               ELSE IF (lSUM < 1.001) THEN
                  COUNT(0,6) = COUNT(0,6) + 1 ! 1.0001 - 1.001
               ELSE IF (lSUM < 1.01) THEN
                  COUNT(0,7) = COUNT(0,7) + 1 ! 1.001  - 1.01
               ELSE IF (lSUM < 1.1) THEN
                  COUNT(0,8) = COUNT(0,8) + 1 ! 1.01   - 1.1
               ELSE
                  COUNT(0,9) = COUNT(0,9) + 1 ! > 1.1
               ENDIF
            ENDIF


            DO M = 1, SMAX
!  Sum of solids mass fractions should be one
               IF(SPECIES_EQ(M)) THEN
                  lSUM = sum(X_S(IJK,M,1:NMAX(M)) )
                  IF(ROP_S(IJK,M) /= ZERO) THEN
                     IF(lSUM < MIN_VAL(M)) THEN
                        MIN_VAL(M) = lSUM
                        MIN_LOC(M) = IJK
                     ENDIF
                     IF (lSUM > MAX_VAL(M)) THEN
                        MAX_VAL(M) = lSUM
                        MAX_LOC(M) = IJK
                     ENDIF
                     IF (lSUM < 0.9) THEN
                        COUNT(M,1) = COUNT(M,1) + 1 ! < 0.9
                     ELSE IF (lSUM < 0.99) THEN
                        COUNT(M,2) = COUNT(M,2) + 1 ! 0.9    - 0.99
                     ELSE IF (lSUM < 0.999) THEN
                        COUNT(M,3) = COUNT(M,3) + 1 ! 0.99   - 0.999
                     ELSE IF (lSUM < 0.9999) THEN
                        COUNT(M,4) = COUNT(M,4) + 1 ! 0.999  - 0.9999
                     ELSE IF (lSUM < 1.0001) THEN
                        COUNT(M,5) = COUNT(M,5) + 1 ! 0.9999 - 1.0001
                     ELSE IF (lSUM < 1.001) THEN
                        COUNT(M,6) = COUNT(M,6) + 1 ! 1.0001 - 1.001
                     ELSE IF (lSUM < 1.01) THEN
                        COUNT(M,7) = COUNT(M,7) + 1 ! 1.001  - 1.01
                     ELSE IF (lSUM < 1.1) THEN
                        COUNT(M,8) = COUNT(M,8) + 1 ! 1.01   - 1.1
                     ELSE
                        COUNT(M,9) = COUNT(M,9) + 1 ! > 1.1
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO


         ENDIF ! IF(.NOT.WALL_AT(IJK))
      ENDDO ! DO I = ISTART2, IEND2
      ENDDO ! DO J = JSTART2, JEND2
      ENDDO ! DO K = KSTART2, KEND2


      CALL GLOBAL_ALL_SUM(COUNT)
      CALL GLOBAL_ALL_MAX(MAX_VAL)
      CALL GLOBAL_ALL_MIN(MIN_VAL)


 4100 FORMAT(//'Statistics of sum of gas species mass fraction',/      &
         1X,'Minimum sum of X_g=',G12.5,/1X,'Maximum sum of X_g=',     &
         G12.5,2/,3x,'Sum of X_g',T20,'No of Cells  Distribution')

      IF(REPORT_SPECIES(0)) THEN
         lSUM = SUM(COUNT(0,:))
         WRITE(ERR_MSG,4100) MIN_VAL(0), MAX_VAL(0)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         WRITE(ERR_MSG,4999) (COUNT(0,L),DBLE(COUNT(0,L))/lSUM,L=1,9)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDIF


 4200 FORMAT(//'Statistics of sum of solids phase ',I2,' species ',    &
         'mass fractions',/1X,'Minimum sum of X_s=',G12.5,/1X,         &
         'Maximum sum of X_s=',G12.5,2/,3x,'Sum of ',A,T20,'No of ',   &
         'Cells',2x,'Distribution')

      DO M=1,SMAX
         IF(REPORT_SPECIES(M)) THEN
            lSUM = SUM(COUNT(M,:))
            WRITE(ERR_MSG,4200) M,  MIN_VAL(M), MAX_VAL(M), &
               trim(iVAR('X_s',M))
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
            WRITE(ERR_MSG,4999) (COUNT(M,L),DBLE(COUNT(M,L))/lSUM,L=1,9)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         ENDIF
      ENDDO

      WRITE(ERR_MSG,"(/'End of report.')")
      CALL FLUSH_ERR_MSG(HEADER=.FALSE.)

 4999 FORMAT(                                    &
         1X,'<0.9            ',T20,I9,T33,G12.5,/&
         1X,' 0.9    - 0.99  ',T20,I9,T33,G12.5,/&
         1X,' 0.99   - 0.999 ',T20,I9,T33,G12.5,/&
         1X,' 0.999  - 0.9999',T20,I9,T33,G12.5,/&
         1X,' 0.9999 - 1.0001',T20,I9,T33,G12.5,/&
         1X,' 1.0001 - 1.001 ',T20,I9,T33,G12.5,/&
         1X,' 1.001  - 1.01  ',T20,I9,T33,G12.5,/&
         1X,' 1.01   - 1.1   ',T20,I9,T33,G12.5,/&
         1X,'>1.1            ',T20,I9,T33,G12.5)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE REPORT_SPECIES_STATS



!----------------------------------------------------------------------!
!                                                                      !
!  Module name: CHECK_RXN_MASS_BALANCE                                 !
!  Author: M. Syamlal                                 Date: 27-OCT-92  !
!                                                                      !
!  Purpose: Verify that the net interphase mass transfer rates sum to  !
!  zero. This check is not needed with updated rates implementation.   !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_RXN_MASS_BALANCE
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE toleranc
      USE fldvar
      USE rxns
      USE visc_s
      USE visc_g
      USE geometry
      USE run
      USE constant
      USE physprop
      USE indices
      USE funits
      USE mpi_utility
      USE discretelement
      USE mms
      USE functions
      use compar, only: DEAD_CELL_AT

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
! Loop counters:
      INTEGER :: I, J, K, IJK
      INTEGER :: M, L, LM
! Temp variable for summations.
      DOUBLE PRECISION :: lSUM

      INTEGER :: COUNT(DIMENSION_M+2, 2)
      INTEGER :: MAX_LOC(DIMENSION_M+2)
      DOUBLE PRECISION :: MAX_VAL(DIMENSION_M+2)
!......................................................................!

      CALL INIT_ERR_MSG('CHECK_RXN_MASS_BALANCE')


      COUNT = 0
      MAX_LOC = 0
      MAX_VAL = ZERO


      DO K = KSTART2, KEND2
      DO J = JSTART2, JEND2
      DO I = ISTART2, IEND2
         IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
         IJK = FUNIJK(I,J,K)
         IF (.NOT.WALL_AT(IJK)) THEN

!---------------------------------------------------------------------//

! The rate of interphase mass transfer must sum to zero over all phases.
            lSUM = SUM_R_G(IJK)
            IF(SMAX > 0) lSUM = lSUM + sum(SUM_R_S(IJK,1:SMAX))

            IF (ABS(lSUM) > SMALL_NUMBER) THEN
               IF(abs(lSUM) > abs(MAX_VAL(1))) THEN
                  MAX_VAL(1) = lSUM
                  MAX_LOC(1) = IJK
               ENDIF
! Bin the fluid cells: Good/Bad
               IF (ABS(lSUM) > TOL_COM) THEN
                  COUNT(1,2) = COUNT(1,2) + 1
               ELSE
                  COUNT(1,1) = COUNT(1,1) + 1
               ENDIF
            ENDIF

! Verify that the net rate of production (SUM_R_x) matches the the total
! amount of mass transferred from other phases.
            DO L = 0, SMAX
               IF (L == 0) THEN
                  lSUM = SUM_R_G(IJK)
               ELSE
                  lSUM = SUM_R_S(IJK,L)
               ENDIF
               DO M = 0, SMAX
                  IF (M > L) THEN
                     LM = L + 1 + (M - 1)*M/2
                     lSUM = lSUM - R_PHASE(IJK,LM)
                  ELSE IF (L > M) THEN
                     LM = M + 1 + (L - 1)*L/2
                     lSUM = lSUM + R_PHASE(IJK,LM)
                  ENDIF
               ENDDO

               IF(ABS(lSUM) > SMALL_NUMBER) THEN
                  IF(abs(lSUM) > abs(MAX_VAL(L+2))) THEN
                     MAX_VAL(L+2) = lSUM
                     MAX_LOC(L+2) = IJK
                  ENDIF

! Force an exit if an imbalance exceeds the tollerance.
                  IF(ABS(lSUM) > TOL_COM) then
                     COUNT(L+2,2) = COUNT(L+2,2) + 1
                  ELSE
! Count the number of cells that do not sum to one, but are below the
! tollerance for force and exit.
                     COUNT(L+2,1) = COUNT(L+2,1) + 1
                  ENDIF
               ENDIF
            END DO

         ENDIF ! IF(.NOT.WALL_AT(IJK))
      ENDDO ! DO I = ISTART2, IEND2
      ENDDO ! DO J = JSTART2, JEND2
      ENDDO ! DO K = KSTART2, KEND2

      CALL GLOBAL_ALL_SUM(COUNT)

      IF(sum(COUNT(:,2)) > 0) THEN

         CALL GLOBAL_ALL_MAX(MAX_VAL)

         WRITE(ERR_MSG,5000)
         CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)

         IF(COUNT(1,2) > 0) THEN
            WRITE(ERR_MSG, 5100) COUNT(1,1), COUNT(1,2), MAX_VAL(1)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         ENDIF

         DO L=0,SMAX
            IF(COUNT(L+2,2) > 0) THEN
               WRITE(ERR_MSG, 5200) L, COUNT(L+2,1), COUNT(L+2,2),&
                  MAX_VAL(L+2)
               CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
            ENDIF
         ENDDO

         WRITE(ERR_MSG,"(/'End of report.')")
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., ABORT=.TRUE.)
      ENDIF

      CALL FINL_ERR_MSG

      RETURN
 5000 FORMAT(5X,'Time = ',G12.5,/&
         'Message: One or more of the following errors detected:',/    &
         '  1. Discrepancies in the reaction rates.',/                 &
         '  2. The rate of production of phases (SUM_R_g or SUM_R_s)',/&
         '     and the interphase mass transfer rates (R_Phase) are',/ &
         '     inconsistent (in subroutine RRATES).',/4X,'I',T14,'J',  &
         T24,'K',T34,'M',T45,'Value')

 5100 FORMAT(//'Sum of all the reaction rates is not zero!',/,         &
         'Number of cells with discrepancy < error tolerance = ',I5,/, &
         'Number of cells with discrepancy > error tolerance = ',I5,/, &
         'Maximum discrepancy = ',G12.5)

 5200 FORMAT(//'Mesage: Production of phase ',I2,' not equal to ',     &
         'total mass transfer',/'from other phases!',/                 &
         'Number of cells with discrepancy < error tolerance = ',I9,/  &
         'Number of cells with discrepancy > error tolerance = ',I9,/  &
         'Maximum discrepancy = ',G12.4)


      END SUBROUTINE CHECK_RXN_MASS_BALANCE



!----------------------------------------------------------------------!
!  Subroutine: REPORT_ERROR                                            !
!                                                                      !
!  Purpose: Manage error messages for CHECK_DATA_30.                   !
!----------------------------------------------------------------------!
      SUBROUTINE REPORT_ERROR(pIER, pI, pJ, pK, VAL,  RELATION, BND, &
         VAR, LC1, LC2)

      INTEGER, INTENT(INOUT) :: pIER
      INTEGER, INTENT(IN) :: pI, pJ, pK
      DOUBLE PRECISION, INTENT(IN) :: BND, VAL
      CHARACTER(LEN=*), INTENT(IN) :: RELATION
      CHARACTER(LEN=*), INTENT(IN) :: VAR
      INTEGER, INTENT(IN), OPTIONAL :: LC1, LC2
      CHARACTER(LEN=32) :: VAR_FULL

      INTEGER :: lIER

      VAR_FULL=''
      IF(PRESENT(LC2)) THEN
         VAR_FULL = iVAR(VAR,LC1,LC2)
      ELSEIF(PRESENT(LC1)) THEN
         VAR_FULL = iVAR(VAR,LC1)
      ELSE
         VAR_FULL = VAR
      ENDIF

      IF(pIER == 0) THEN
         CALL OPEN_PE_LOG(lIER)
         SELECT CASE(ERR_TAG)
         CASE(2000); WRITE(ERR_MSG,2000)
         CASE(3000); WRITE(ERR_MSG,3000)
         END SELECT
         CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)
         pIER = 1
      ENDIF

 2000 FORMAT('Error 2000: Physical properties detected in flow cells.',&
         2/3X,'I',6x,'J',6x,'K',5x,'Value',8x,A,2x,'Bound',5x,'Variable')

 3000 FORMAT('Error 3000: Unphysical field variables detected.',&
         2/3X,'I',6x,'J',6x,'K',5x,'Value',8x,A,2x,'Bound',5x,'Variable')

      WRITE(ERR_MSG,9000) pI, pJ, pK, VAL, RELATION, BND, trim(VAR_FULL)
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

 9000 FORMAT(3(I6,1X),G12.4,1X,A,G12.4,1X,A)

      RETURN
      END SUBROUTINE REPORT_ERROR
      END SUBROUTINE CHECK_DATA_30
