!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_RESID_C                                            C
!  Purpose: Calculate residuals for continuity equations               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_RESID_C(VAR, A_M, B_M, M, NUM, DEN, &
         RESID, MAX_RESID, IJK_RESID)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1, ONLY: ZERO, ONE, UNDEFINED
      USE parallel
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE run
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Variable being evaluated
      DOUBLE PRECISION, INTENT(IN) :: Var(DIMENSION_3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Phase index
      INTEGER, INTENT(IN) :: M
! Numerator and denominator
      DOUBLE PRECISION, INTENT(OUT) :: NUM, DEN
! Average value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: RESID
! Maximum value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: MAX_RESID
! IJK of Maximum value of Residual
      INTEGER, INTENT(OUT) :: IJK_RESID
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, IJKW, IJKS, IJKB, IJKE, IJKN, IJKT
! Numerators and denominators
      DOUBLE PRECISION :: NUM1, DEN1
! Number of fluid cells
      INTEGER :: NCELLS
! New local variables for DMP version
      DOUBLE PRECISION, DIMENSION(ijksize3_all(myPE)) :: RESID_IJK
      DOUBLE PRECISION :: MAX_RESID_GL(0:numPEs-1), MAX_RESID_L(0:numPEs-1)
      INTEGER :: IJK_RESID_GL(0:numPEs-1), IJK_RESID_L(0:numPEs-1)
      INTEGER :: nproc
!-----------------------------------------------

! initializing values
      NUM = ZERO
      DEN = ZERO
      MAX_RESID = -ONE
      NCELLS = 0

!!$omp parallel do private( IJK )
      DO IJK = ijkstart3, ijkend3
          RESID_IJK(IJK) = ZERO
      ENDDO

!!$omp  parallel do private( IJK, IJKW, IJKS, IJKB, IJKE, IJKN, IJKT,  &
!!$omp&  NUM1, DEN1) &
!!$omp&  REDUCTION(+:NUM,DEN,NCELLS)
      DO IJK = ijkstart3, ijkend3
         IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE
         IF (FLUID_AT(IJK)) THEN
            IJKW = WEST_OF(IJK)
            IJKS = SOUTH_OF(IJK)
            IJKE = EAST_OF(IJK)
            IJKN = NORTH_OF(IJK)

! evaluating the residual at cell ijk:
!   RESp = B-sum(Anb*VARnb)-Ap*VARp
!   (where nb = neighbor cells and p = center/0 cell)
            NUM1 = B_M(IJK,M) - (A_M(IJK,0,M)*VAR(IJK)+A_M(IJK,east,M)*VAR(IJKE)+&
               A_M(IJK,west,M)*VAR(IJKW)+A_M(IJK,north,M)*VAR(IJKN)+A_M(IJK,south,M)*VAR(&
               IJKS))

            IF (DO_K) THEN
               IJKB = BOTTOM_OF(IJK)
               IJKT = TOP_OF(IJK)
               NUM1 = NUM1 - (A_M(IJK,top,M)*VAR(IJKT)+A_M(IJK,bottom,M)*VAR(IJKB))
            ENDIF

            NUM1 = ABS(NUM1)
            DEN1 = ABS(A_M(IJK,0,M)*VAR(IJK))
! storing value of residual at each ijk location
            RESID_IJK(IJK) = NUM1

! adding to terms that are accumulated
            NCELLS = NCELLS + 1
            NUM = NUM + NUM1
            DEN = DEN + DEN1
         ENDIF
      ENDDO

      IF(.not.debug_resid) RETURN

! Collecting all the information among all the procesors -
! determining the global sum
      call global_all_sum(NUM)
      call global_all_sum(DEN)
      call global_all_sum(NCELLS)

      IJK_RESID = 1
      MAX_RESID = RESID_IJK( IJK_RESID )
      DO IJK = ijkstart3, ijkend3
         IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE
         IF (RESID_IJK(IJK) > MAX_RESID) THEN
            IJK_RESID = IJK
            MAX_RESID = RESID_IJK( IJK_RESID )
         ENDIF
      ENDDO

! Determining the max residual
      do nproc=0,NumPEs-1
         if(nproc.eq.myPE) then
            MAX_RESID_L(nproc) = MAX_RESID
            IJK_RESID_L(nproc) = FUNIJK_GL(I_OF(IJK_RESID), J_OF(IJK_RESID), K_OF(IJK_RESID))
         else
            MAX_RESID_L(nproc) = 0.0
            IJK_RESID_L(nproc) = 0
         endif
      enddo

! Determining the maximum among all the procesors
      call global_all_max(MAX_RESID)

! Collecting all the information among all the procesors
      call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
      call global_all_sum(IJK_RESID_L, IJK_RESID_GL)

! Calling to determine the global IJK location w.r.t. serial version
      IJK_RESID = IJKMAX2
      do nproc=0,NumPEs-1
         if(MAX_RESID_GL(nproc).eq.MAX_RESID.and.IJK_RESID_GL(nproc).lt.IJK_RESID) then
            IJK_RESID = IJK_RESID_GL(nproc)
         endif
      enddo

! Normalizing the residual
      IF (DEN > ZERO) THEN
         RESID = NUM/DEN
         MAX_RESID = NCELLS*MAX_RESID/DEN
      ELSEIF (NUM == ZERO) THEN
         RESID = ZERO
         MAX_RESID = ZERO
         IJK_RESID = 0
      ELSE
         RESID = UNDEFINED
         MAX_RESID = UNDEFINED
!         WRITE (LINE, *) 'Warning: All center coefficients are zero.'
!         CALL WRITE_ERROR ('CALC_RESID_C', LINE, 1)
      ENDIF

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

      END SUBROUTINE CALC_RESID_C

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_RESID_S                                            C
!  Purpose: Calculate residuals for scalar equations                   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_RESID_S(VAR, A_M, B_M, M, NUM, DEN, &
         RESID, MAX_RESID, IJK_RESID, TOL)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE run

      USE fldvar
      USE physprop
      USE toleranc

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Variable being evaluated
      DOUBLE PRECISION, INTENT(IN) :: Var(DIMENSION_3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Phase index
      INTEGER, INTENT(IN) :: M
! Numerator and denominator
      DOUBLE PRECISION, INTENT(OUT) :: NUM, DEN
! Average value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: RESID
! Maximum value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: MAX_RESID
! IJK of Maximum value of Residual
      INTEGER, INTENT(OUT) :: IJK_RESID
! Ignore residual calculation for scalar values below this
      DOUBLE PRECISION, INTENT(IN) :: TOL
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP
! Numerators and denominators
      DOUBLE PRECISION :: NUM1, DEN1
! Number of fluid cells
      INTEGER :: NCELLS
! New local variables for DMP version
      DOUBLE PRECISION, DIMENSION(ijksize3_all(myPE)) :: RESID_IJK
      DOUBLE PRECISION :: MAX_RESID_GL(0:numPEs-1), MAX_RESID_L(0:numPEs-1)
      INTEGER :: IJK_RESID_GL(0:numPEs-1), IJK_RESID_L(0:numPEs-1)
      INTEGER :: nproc
!-----------------------------------------------

! initializing
      NUM = ZERO
      DEN = ZERO
      MAX_RESID = -ONE
      NCELLS = 0

!!$omp parallel do private( IJK )
      DO IJK = ijkstart3, ijkend3
         RESID_IJK(IJK) = ZERO
      ENDDO

!!$omp    parallel do &
!!$omp&   private(   IJK, IMJK, IJMK, IPJK, IJPK, IJKM, IJKP, &
!!$omp&   NUM1, DEN1) &
!!$omp&   REDUCTION(+:NUM, DEN,NCELLS)

      DO IJK = ijkstart3, ijkend3
         IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE

         IF (FLUID_AT(IJK) .AND. ABS(VAR(IJK)) > TOL) THEN
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)

            if(M/=0) then
               if(EP_S(IJK,M) <= DIL_EP_s) CYCLE
            endif


! evaluating the residual at cell ijk:
!   RESp = B-sum(Anb*VARnb)-Ap*VARp
!   (where nb = neighbor cells and p = center/0 cell)
            NUM1 = B_M(IJK,M) - (A_M(IJK,0,M)*VAR(IJK)+A_M(IJK,east,M)*VAR(IPJK)+&
               A_M(IJK,west,M)*VAR(IMJK)+A_M(IJK,north,M)*VAR(IJPK)+A_M(IJK,south,M)*VAR(&
               IJMK))
            IF (DO_K) THEN
               IJKM = KM_OF(IJK)
               IJKP = KP_OF(IJK)
               NUM1 = NUM1 - (A_M(IJK,top,M)*VAR(IJKP)+A_M(IJK,bottom,M)*VAR(IJKM))
            ENDIF

            NUM1 = ABS(NUM1)
            DEN1 = ABS(A_M(IJK,0,M)*VAR(IJK))
! storing value of residual at each ijk location
            RESID_IJK(IJK) = NUM1

! adding to terms that are accumulated
            NCELLS = NCELLS + 1
            NUM = NUM + NUM1
            DEN = DEN + DEN1
         ENDIF
      ENDDO

      IF(.not.debug_resid) RETURN

! Collecting all the information among all the procesors -
! determining the global sum
      call global_all_sum(NUM)
      call global_all_sum(DEN)
      call global_all_sum(NCELLS)

      IJK_RESID = 1
      MAX_RESID = RESID_IJK( IJK_RESID )
      DO IJK = ijkstart3, ijkend3
      IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE
         IF (RESID_IJK(IJK) > MAX_RESID) THEN
               IJK_RESID = IJK
               MAX_RESID = RESID_IJK( IJK_RESID )
         ENDIF
      ENDDO

! Determining the max residual
      do nproc=0,NumPEs-1
         if(nproc.eq.myPE) then
            MAX_RESID_L(nproc) = MAX_RESID
            IJK_RESID_L(nproc) = FUNIJK_GL(I_OF(IJK_RESID), J_OF(IJK_RESID), K_OF(IJK_RESID))
         else
            MAX_RESID_L(nproc) = 0.0
            IJK_RESID_L(nproc) = 0
         endif
      enddo

! Determining the maximum among all the procesors
      call global_all_max(MAX_RESID)

! Collecting all the information among all the procesors
      call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
      call global_all_sum(IJK_RESID_L, IJK_RESID_GL)

! Determining the global IJK location w.r.t. serial version
      IJK_RESID = IJKMAX2
      do nproc=0,NumPEs-1
        if(MAX_RESID_GL(nproc).eq.MAX_RESID.and.IJK_RESID_GL(nproc).lt.IJK_RESID) then
           IJK_RESID = IJK_RESID_GL(nproc)
        endif
      enddo

! Normalizing the residual
      IF (DEN > ZERO) THEN
         RESID = NUM/DEN
         MAX_RESID = NCELLS*MAX_RESID/DEN
      ELSEIF (NUM == ZERO) THEN
         RESID = ZERO
         MAX_RESID = ZERO
         IJK_RESID = 0
      ELSE
         RESID = UNDEFINED
         MAX_RESID = UNDEFINED
!         WRITE(LINE,*)'Message: All center coefficients are zero.'
!         CALL WRITE_ERROR('CALC_RESID_S', LINE, 1)
      ENDIF

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

      END SUBROUTINE CALC_RESID_S

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_RESID_pp                                           C
!  Purpose: Calculate residuals for pressure correction equation       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Comments:                                                           C
!  for correction equations the convergence for the corrections must   C
!  go to zero, therefore the vector b must go to zero. this value      C
!  cannot be normalized as the other equations are since the           C
!  denominator here will vanish.  thus the residual is normalized      C
!  based on its value in the first iteration                           C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_RESID_PP(B_M, NORM, NUM, DEN, RESID, MAX_RESID, &
         IJK_RESID)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE run
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Normalization factor
      DOUBLE PRECISION, INTENT(IN) :: NORM
! Numerator and denominator
      DOUBLE PRECISION, INTENT(OUT) :: NUM, DEN
! Average value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: RESID
! Maximum value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: MAX_RESID
! IJK of Maximum value of Residual
      INTEGER, INTENT(OUT) :: IJK_RESID
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK
! Number of fluid cells
      INTEGER :: NCELLS
! Numerators and denominators
      DOUBLE PRECISION :: NUM1, DEN1
! New local variables for DMP version
      DOUBLE PRECISION :: MAX_RESID_GL(0:numPEs-1), MAX_RESID_L(0:numPEs-1)
      INTEGER :: IJK_RESID_GL(0:numPEs-1), IJK_RESID_L(0:numPEs-1)
      INTEGER :: nproc
!-----------------------------------------------

! initializing values
      NUM = ZERO
      DEN = ZERO
      MAX_RESID = -ONE
      NCELLS = 0
      DEN1 = ONE
      IJK_RESID = 1

      DO IJK = ijkstart3, ijkend3
         IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE
         IF (FLUID_AT(IJK)) THEN

! evaluating the residual at cell ijk:
            NUM1 = ABS(B_M(IJK,0))
            IF (NUM1 > MAX_RESID) THEN
               MAX_RESID = NUM1
               IJK_RESID = IJK
            ENDIF

! adding to terms that are accumulated
            NCELLS = NCELLS + 1
            NUM = NUM + NUM1
            DEN = DEN + DEN1
         ENDIF
      ENDDO


      IF(.not.debug_resid) THEN
! Collecting all the information among all the procesors
         call global_all_sum(NUM)
         call global_all_sum(DEN)

! Normalizing the residual
         IF (DEN*NORM > ZERO) THEN
! if norm=1 then this simply becomes an unscaled 'average' residual
            RESID = NUM/(DEN*NORM)
         ELSEIF (NUM == ZERO) THEN
            RESID = ZERO
         ELSE
            RESID = LARGE_NUMBER
         ENDIF
      ELSE   ! if(debug_resid) branch

! Collecting all the information among all the procesors -
! determining the global sum
         call global_all_sum(NUM)
         call global_all_sum(DEN)
         call global_all_sum(NCELLS)

! Determining the max residual
         do nproc=0,NumPEs-1
            if(nproc.eq.myPE) then
               MAX_RESID_L(nproc) = MAX_RESID
               IJK_RESID_L(nproc) = FUNIJK_GL(I_OF(IJK_RESID), J_OF(IJK_RESID), K_OF(IJK_RESID))
            else
               MAX_RESID_L(nproc) = 0.0
               IJK_RESID_L(nproc) = 0
            endif
         enddo

! Determining the maximum among all the procesors
         call global_all_max(MAX_RESID)

! Collecting all the information among all the procesors
         call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
         call global_all_sum(IJK_RESID_L, IJK_RESID_GL)

! Determining the global IJK location w.r.t. serial version
         IJK_RESID = IJKMAX2
         do nproc=0,NumPEs-1
            if(MAX_RESID_GL(nproc).eq.MAX_RESID.and.IJK_RESID_GL(nproc).lt.IJK_RESID) then
               IJK_RESID = IJK_RESID_GL(nproc)
            endif
         enddo

! Normalizing the residual
         IF (DEN*NORM > ZERO) THEN
            RESID = NUM/(DEN*NORM)
            MAX_RESID = NCELLS*MAX_RESID/(DEN*NORM)
         ELSEIF (NUM == ZERO) THEN
            RESID = ZERO
            MAX_RESID = ZERO
            IJK_RESID = 0
         ELSE
            RESID = LARGE_NUMBER
            MAX_RESID = LARGE_NUMBER
!            WRITE (LINE, *) 'Warning: All center coefficients are zero.'
!             CALL WRITE_ERROR ('CALC_RESID_pp', LINE, 1)
         ENDIF

      ENDIF   ! end if/else debug_resid branch

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

      END SUBROUTINE CALC_RESID_PP

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_RESID_mb                                           C
!  Purpose: Calculate overall mass balance error                       C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 9-DEC-02   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_RESID_MB(INIT, ErrorPercent)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE bc
      USE check, only: calc_mass_fluxhr, accumulation
      USE compar
      USE constant
      USE fldvar
      USE geometry
      USE indices
      USE mflux
      USE mpi_utility
      USE parallel
      USE param
      USE param1
      USE physprop
      USE residual
      USE run, only: Added_Mass, dt, m_am, steady_state
      USE rxns
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Flag to check whether this is an initialization call
! 0 -> initialize old accumulation; 1 -> calc residual
      INTEGER, INTENT(IN) :: init
! Total mass balance error as a % of inflow
      DOUBLE PRECISION, INTENT(OUT) :: ErrorPercent(0:MMAX)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! phase index
      INTEGER :: M
! boundary index
      INTEGER :: L
! locally define dt, so that this routine works when dt is not defined
      DOUBLE PRECISION :: dt_local
! error index
      INTEGER :: IER
      DOUBLE PRECISION :: flux_in, flux_out, fin, fout
      DOUBLE PRECISION :: err, accum_new, denom
!-----------------------------------------------

      if(STEADY_STATE)then
        dt_local = ONE
      else
        dt_local = dt
      endif

      IF(init == 0) THEN
! Initialize this routine
! ---------------------------------------------------------------->>>

! Accumulation
        if(STEADY_STATE)then
          Accum_resid_g = ZERO
        else
          Accum_resid_g = Accumulation(ROP_g)
        endif
        DO M=1, MMAX
          if(STEADY_STATE)then
            Accum_resid_s(M) = ZERO
          else
            Accum_resid_s(M) = Accumulation(ROP_s(1,M))
          endif
        ENDDO
        RETURN
! end initialization
! ----------------------------------------------------------------<<<

      ELSE

! Calculate residual
! ---------------------------------------------------------------->>>
        if(STEADY_STATE)then
          Accum_new = - Accumulation(SUM_R_g) * dt_local
        else
          Accum_new = Accumulation(ROP_g) - Accumulation(SUM_R_g) * dt_local
        endif

        flux_out = zero
        flux_in = zero
        DO L = 1, DIMENSION_BC
          IF (BC_DEFINED(L)) THEN
!            call Calc_mass_flux(BC_I_W(L), BC_I_E(L), BC_J_S(L), &
!            BC_J_N(L), BC_K_B(L), BC_K_T(L), BC_PLANE(L), U_g, V_g, W_g, &
!            ROP_g, fin, fout, IER)
            IF(.NOT.Added_Mass) THEN
              call Calc_mass_fluxHR(BC_I_W(L), BC_I_E(L), BC_J_S(L), &
                 BC_J_N(L), BC_K_B(L), BC_K_T(L), BC_PLANE(L), &
                 Flux_gE, Flux_gN, Flux_gT, fin, fout, IER)
            ELSE
              call Calc_mass_fluxHR(BC_I_W(L), BC_I_E(L), BC_J_S(L), &
                 BC_J_N(L), BC_K_B(L), BC_K_T(L), BC_PLANE(L), Flux_gSE,&
                 Flux_gSN, Flux_gST, fin, fout, IER)
            ENDIF
            flux_out = flux_out + fout  * dt_local
            flux_in = flux_in + fin * dt_local
          ENDIF
        END DO

        Err = (accum_new - Accum_resid_g) - (flux_in - flux_out)
        denom = max(abs(accum_new), abs(Accum_resid_g), abs(flux_in), abs(flux_out))
        IF (denom /= ZERO) THEN
           ErrorPercent(0) = err*100./denom
        ELSE
           ErrorPercent(0) = err*100./SMALL_NUMBER
        ENDIF

        DO M =1, MMAX
          if(STEADY_STATE)then
            Accum_new =  - Accumulation(SUM_R_s(1,M)) * dt_local
          else
            Accum_new = Accumulation(ROP_s(1,M)) - Accumulation(SUM_R_s(1,M)) * dt_local
          endif

          flux_out = zero
          flux_in = zero
          DO L = 1, DIMENSION_BC
            IF (BC_DEFINED(L)) THEN
!              call Calc_mass_flux(BC_I_W(L), BC_I_E(L), BC_J_S(L), BC_J_N(L), &
!              BC_K_B(L), BC_K_T(L), BC_PLANE(L), U_s(1,M), V_s(1,M), W_s(1,M), &
!              ROP_s(1,M), fin, fout, IER)
              IF(.NOT.Added_Mass .OR. M /= M_AM) THEN
                call Calc_mass_fluxHR(BC_I_W(L), BC_I_E(L), BC_J_S(L), &
                   BC_J_N(L), BC_K_B(L), BC_K_T(L), BC_PLANE(L), &
                   Flux_sE(1,M), Flux_sN(1,M), Flux_sT(1,M), fin, fout, IER)
              ELSE
                call Calc_mass_fluxHR(BC_I_W(L), BC_I_E(L), BC_J_S(L), &
                   BC_J_N(L), BC_K_B(L), BC_K_T(L), BC_PLANE(L), &
                   Flux_sSE, Flux_sSN, Flux_sST, fin, fout, IER)
              ENDIF
              flux_out = flux_out + fout  * dt_local
              flux_in = flux_in + fin * dt_local
            ENDIF
          ENDDO

          Err = (accum_new - Accum_resid_s(M)) - (flux_in - flux_out)
          denom = max(abs(accum_new), abs(Accum_resid_s(M)), abs(flux_in), abs(flux_out))
          if(denom /= ZERO) THEN
            ErrorPercent(M) = err*100./denom
          else
            ErrorPercent(M) = err*100./SMALL_NUMBER
          endif
        ENDDO  ! end do m=1,mmax
! end calculate residual
! ----------------------------------------------------------------<<<

      ENDIF   ! end if/else (init==0)

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

      END SUBROUTINE CALC_RESID_MB

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_RESID_U                                            C
!  Purpose: Calculate residuals for u-momentum equations               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_RESID_U(U_M, V_M, W_M, A_M, B_M, M, NUM, DEN, &
         RESID, MAX_RESID, IJK_RESID)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE run
      USE fldvar
      USE physprop
      USE toleranc
      USE fun_avg

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! U velocity (x-dir)
      DOUBLE PRECISION, INTENT(IN) :: U_m(DIMENSION_3)
! V velocity (y-dir), used here for scaling
      DOUBLE PRECISION, INTENT(IN) :: V_m(DIMENSION_3)
! W velocity (z-dir), used here for scaling
      DOUBLE PRECISION, INTENT(IN) :: W_m(DIMENSION_3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Phase index
      INTEGER, INTENT(IN) :: M
! Numerator and denominator
      DOUBLE PRECISION, INTENT(OUT) :: NUM, DEN
! Average value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: RESID
! Maximum value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: MAX_RESID
! IJK of Maximum value of Residual
      INTEGER, INTENT(OUT) :: IJK_RESID
!-----------------------------------------------
!     Local variables
!-----------------------------------------------
! Velocity magnitude
      DOUBLE PRECISION :: VEL
! Indices
      INTEGER :: IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP
! Numerators and denominators
      DOUBLE PRECISION :: NUM1, DEN1
! Number of fluid cells
      INTEGER :: NCELLS
! New local variables for DMP version
      DOUBLE PRECISION, DIMENSION(ijksize3_all(myPE)) :: RESID_IJK
      DOUBLE PRECISION :: MAX_RESID_GL(0:numPEs-1), MAX_RESID_L(0:numPEs-1)
      INTEGER :: IJK_RESID_GL(0:numPEs-1), IJK_RESID_L(0:numPEs-1)
      INTEGER :: nproc

! Solids volume fraction at face
      DOUBLE PRECISION :: EPSA

!-----------------------------------------------

! initializing
      NUM = ZERO
      DEN = ZERO
      MAX_RESID = -ONE
      NCELLS = 0

!$omp  parallel default(none) &
!$omp  private( IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, &
!$omp           VEL,  NUM1, DEN1,EPSA) &
!$omp  shared (ijkstart3, ijkend3, resid_ijk, i_of, j_of, k_of,m,a_m,b_m,w_m,do_k,u_m,v_m,num,den,ncells)
!$omp do reduction(+:num, DEN, NCELLS )
      DO IJK = ijkstart3, ijkend3
        RESID_IJK(IJK) = ZERO

        IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE

! Skip walls where some values are undefined.
        IF(WALL_AT(IJK)) cycle

         if(m/=0) then
            EPSA = AVG_X(EP_S(IJK,M),EP_S(EAST_OF(IJK),M),I_OF(IJK))
            if(EPSA <= DIL_EP_s) CYCLE
         endif

         IF (.NOT.IP_AT_E(IJK)) THEN
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)

! evaluating the residual at cell ijk:
!   RESp = B-sum(Anb*VARnb)-Ap*VARp
!   (where nb = neighbor cells and p = center/0 cell)
            NUM1 = B_M(IJK,M) - (A_M(IJK,0,M)*U_M(IJK)+&
               A_M(IJK,east,M)*U_M(IPJK)+A_M(IJK,west,M)*U_M(IMJK)+&
               A_M(IJK,north,M)*U_M(IJPK)+A_M(IJK,south,M)*U_M(IJMK))
            IF (DO_K) THEN
               IJKM = KM_OF(IJK)
               IJKP = KP_OF(IJK)
               NUM1 = NUM1 - (A_M(IJK,top,M)*U_M(IJKP)+A_M(IJK,bottom,M)*U_M(IJKM))
            ENDIF

! Ignore momentum residual in stagnant regions.  Need an alternative
! criteria for residual scaling for such cases.
            VEL = SQRT(U_M(IJK)**2+V_M(IJK)**2+W_M(IJK)**2)
            IF (VEL > SMALL_NUMBER) THEN
               NUM1 = ABS(NUM1)
               DEN1 = ABS(A_M(IJK,0,M)*VEL)
! storing value of residual at each ijk location
               RESID_IJK(IJK) = NUM1
! adding to terms that are accumulated
               NCELLS = NCELLS + 1
               NUM = NUM + NUM1
               DEN = DEN + DEN1
            ENDIF
         ENDIF
      ENDDO
!$omp end parallel

      IF(.not.debug_resid) RETURN

! Collecting all the information among all the procesors -
! determining the global sum
      call global_all_sum(NUM)
      call global_all_sum(DEN)
      call global_all_sum(NCELLS)

      IJK_RESID = 1
      MAX_RESID = RESID_IJK( IJK_RESID )
      DO IJK = ijkstart3, ijkend3
         IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE
         IF (RESID_IJK( IJK ) > MAX_RESID) THEN
            IJK_RESID = IJK
            MAX_RESID = RESID_IJK( IJK_RESID )
         ENDIF
      ENDDO

! Determining the max residual
      do nproc=0,NumPEs-1
         if(nproc.eq.myPE) then
            MAX_RESID_L(nproc) = MAX_RESID
            IJK_RESID_L(nproc) = FUNIJK_GL(I_OF(IJK_RESID), J_OF(IJK_RESID), K_OF(IJK_RESID))
         else
            MAX_RESID_L(nproc) = 0.0
            IJK_RESID_L(nproc) = 0
         endif
      enddo

! Determining the maximum among all the procesors
      call global_all_max(MAX_RESID)

! Collecting all the information among all the procesors
      call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
      call global_all_sum(IJK_RESID_L, IJK_RESID_GL)

! Determining the global IJK location w.r.t. serial version
      IJK_RESID = IJKMAX2
      do nproc=0,NumPEs-1
        IF(MAX_RESID_GL(nproc).eq.MAX_RESID.and.&
           IJK_RESID_GL(nproc).lt.IJK_RESID) THEN
          IJK_RESID = IJK_RESID_GL(nproc)
        ENDIF
      ENDDO

! Normalizing the residual
      IF (DEN > ZERO) THEN
         RESID = NUM/DEN
         MAX_RESID = NCELLS*MAX_RESID/DEN
      ELSEIF (NUM == ZERO) THEN
         RESID = ZERO
         MAX_RESID = ZERO
         IJK_RESID = 0
      ELSE
         RESID = LARGE_NUMBER
         MAX_RESID = LARGE_NUMBER
!         WRITE (LINE, *) 'Warning: All center coefficients are zero.'
!         CALL WRITE_ERROR ('CALC_RESID_U', LINE, 1)
      ENDIF

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

      END SUBROUTINE CALC_RESID_U

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_RESID_V                                            C
!  Purpose: Calculate residuals for v-momentum equations               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_RESID_V(U_M, V_M, W_M, A_M, B_M, M, NUM, DEN, &
         RESID, MAX_RESID, IJK_RESID)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE run
      USE fldvar
      USE physprop
      USE toleranc
      USE fun_avg

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! U velocity (x-dir)
      DOUBLE PRECISION, INTENT(IN) :: U_m(DIMENSION_3)
! V velocity (y-dir), used here for scaling
      DOUBLE PRECISION, INTENT(IN) :: V_m(DIMENSION_3)
! W velocity (z-dir), used here for scaling
      DOUBLE PRECISION, INTENT(IN) :: W_m(DIMENSION_3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Phase index
      INTEGER, INTENT(IN) :: M
! Numerator and denominator
      DOUBLE PRECISION, INTENT(OUT) :: NUM, DEN
! Average value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: RESID
! Maximum value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: MAX_RESID
! IJK of Maximum value of Residual
      INTEGER, INTENT(OUT) :: IJK_RESID
!-----------------------------------------------
!     Local variables
!-----------------------------------------------
! Velocity magnitude
      DOUBLE PRECISION :: VEL
! Indices
      INTEGER :: IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP
! Numerators and denominators
      DOUBLE PRECISION :: NUM1, DEN1
! Number of fluid cells
      INTEGER :: NCELLS
! New local variables for DMP version
      DOUBLE PRECISION, DIMENSION(ijksize3_all(myPE)) :: RESID_IJK
      DOUBLE PRECISION :: MAX_RESID_GL(0:numPEs-1), MAX_RESID_L(0:numPEs-1)
      INTEGER :: IJK_RESID_GL(0:numPEs-1), IJK_RESID_L(0:numPEs-1)
      INTEGER :: nproc
! Solids volume fraction at face
      DOUBLE PRECISION :: EPSA

!-----------------------------------------------

! initializing
      NUM = ZERO
      DEN = ZERO
      MAX_RESID = -ONE
      NCELLS = 0

!$omp  parallel default(none) &
!$omp  private( IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, &
!$omp           VEL,  NUM1, DEN1,EPSA) &
!$omp  shared (ijkstart3, ijkend3, resid_ijk, i_of, j_of, k_of,m,a_m,b_m,w_m,do_k,u_m,v_m,num,den,ncells)
!$omp do reduction(+:num, DEN, NCELLS )
      DO IJK = ijkstart3, ijkend3
        RESID_IJK(IJK) = ZERO

        IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE

! Skip walls where some values are undefined.
        IF(WALL_AT(IJK)) cycle

         if(m/=0) then
            EPSA = AVG_Y(EP_S(IJK,M),EP_S(NORTH_OF(IJK),M),J_OF(IJK))
            if(EPSA <= DIL_EP_s) CYCLE
         endif

         IF (.NOT.IP_AT_N(IJK)) THEN
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)

! evaluating the residual at cell ijk:
!   RESp = B-sum(Anb*VARnb)-Ap*VARp
!   (where nb = neighbor cells and p = center/0 cell)
            NUM1 = B_M(IJK,M) - (A_M(IJK,0,M)*V_M(IJK)+&
               A_M(IJK,east,M)*V_M(IPJK)+A_M(IJK,west,M)*V_M(IMJK)+&
               A_M(IJK,north,M)*V_M(IJPK)+A_M(IJK,south,M)*V_M(IJMK))
            IF (DO_K) THEN
               IJKM = KM_OF(IJK)
               IJKP = KP_OF(IJK)
               NUM1 = NUM1 - (A_M(IJK,top,M)*V_M(IJKP)+A_M(IJK,bottom,M)*V_M(IJKM))
            ENDIF

! Ignore momentum residual in stagnant regions.  Need an alternative
! criteria for residual scaling for such cases.
            VEL = SQRT(U_M(IJK)**2+V_M(IJK)**2+W_M(IJK)**2)
            IF (VEL > SMALL_NUMBER) THEN
               NUM1 = ABS(NUM1)
               DEN1 = ABS(A_M(IJK,0,M)*VEL)
! storing value of residual at each ijk location
               RESID_IJK(IJK) = NUM1
! adding to terms that are accumulated
               NCELLS = NCELLS + 1
               NUM = NUM + NUM1
               DEN = DEN + DEN1
            ENDIF
         ENDIF
      ENDDO
!$omp end parallel

      if(.not.debug_resid) return

! Collecting all the information among all the procesors -
! determining the global sum
      call global_all_sum(NUM)
      call global_all_sum(DEN)
      call global_all_sum(NCELLS)

      IJK_RESID = 1
      MAX_RESID = RESID_IJK( IJK_RESID )
      DO IJK = ijkstart3, ijkend3
      IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE
          IF (RESID_IJK( IJK ) > MAX_RESID) THEN
              IJK_RESID = IJK
              MAX_RESID = RESID_IJK( IJK_RESID )
          ENDIF
      ENDDO

! Determining the max residual
      do nproc=0,NumPEs-1
         if(nproc.eq.myPE) then
            MAX_RESID_L(nproc) = MAX_RESID
            IJK_RESID_L(nproc) = FUNIJK_GL(I_OF(IJK_RESID), J_OF(IJK_RESID), K_OF(IJK_RESID))
         else
            MAX_RESID_L(nproc) = 0.0
            IJK_RESID_L(nproc) = 0
         endif
      enddo

! Determine the maximum among all the procesors
      call global_all_max(MAX_RESID)

! Collecting all the information among all the procesors
      call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
      call global_all_sum(IJK_RESID_L, IJK_RESID_GL)

! Determining the global IJK location w.r.t. serial version
      IJK_RESID = IJKMAX2
      do nproc=0,NumPEs-1
        if(MAX_RESID_GL(nproc).eq.MAX_RESID.and.IJK_RESID_GL(nproc).lt.IJK_RESID) then
          IJK_RESID = IJK_RESID_GL(nproc)
        endif
      enddo

! Normalizing the residual
      IF (DEN > ZERO) THEN
         RESID = NUM/DEN
         MAX_RESID = NCELLS*MAX_RESID/DEN
      ELSEIF (NUM == ZERO) THEN
         RESID = ZERO
         MAX_RESID = ZERO
         IJK_RESID = 0
      ELSE
         RESID = LARGE_NUMBER
         MAX_RESID = LARGE_NUMBER
!         WRITE (LINE, *) 'Warning: All center coefficients are zero.'
!         CALL WRITE_ERROR ('CALC_RESID_V', LINE, 1)
      ENDIF

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

      END SUBROUTINE CALC_RESID_V

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_RESID_W                                            C
!  Purpose: Calculate residuals for w-momentum equations               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_RESID_W(U_M, V_M, W_M, A_M, B_M, M, NUM, DEN, &
         RESID, MAX_RESID, IJK_RESID)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE run
      USE fldvar
      USE physprop
      USE toleranc
      USE fun_avg

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! U velocity (x-dir)
      DOUBLE PRECISION, INTENT(IN) :: U_m(DIMENSION_3)
! V velocity (y-dir), used here for scaling
      DOUBLE PRECISION, INTENT(IN) :: V_m(DIMENSION_3)
! W velocity (z-dir), used here for scaling
      DOUBLE PRECISION, INTENT(IN) :: W_m(DIMENSION_3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Phase index
      INTEGER, INTENT(IN) :: M
! Numerator and denominator
      DOUBLE PRECISION, INTENT(OUT) :: NUM, DEN
! Average value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: RESID
! Maximum value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: MAX_RESID
! IJK of Maximum value of Residual
      INTEGER, INTENT(OUT) :: IJK_RESID
!-----------------------------------------------
!     Local variables
!-----------------------------------------------
! Velocity magnitude
      DOUBLE PRECISION :: VEL
! Indices
      INTEGER :: IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP
! Numerators and denominators
      DOUBLE PRECISION :: NUM1, DEN1
! Number of fluid cells
      INTEGER :: NCELLS
! New local variables for DMP version
      DOUBLE PRECISION, DIMENSION(ijksize3_all(myPE)) :: RESID_IJK
      DOUBLE PRECISION :: MAX_RESID_GL(0:numPEs-1), MAX_RESID_L(0:numPEs-1)
      INTEGER :: IJK_RESID_GL(0:numPEs-1), IJK_RESID_L(0:numPEs-1)
      INTEGER :: nproc
! Solids volume fraction at face
      DOUBLE PRECISION :: EPSA

!-----------------------------------------------

! initializing
      NUM = ZERO
      DEN = ZERO
      MAX_RESID = -ONE
      NCELLS = 0

!$omp  parallel default(none) &
!$omp  private( IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, &
!$omp           VEL,  NUM1, DEN1,EPSA) &
!$omp  shared (ijkstart3, ijkend3, resid_ijk, i_of, j_of, k_of,m,a_m,b_m,w_m,do_k,u_m,v_m,num,den,ncells)
!$omp do reduction(+:num, DEN, NCELLS )
      DO IJK = ijkstart3, ijkend3
        RESID_IJK(IJK) = ZERO

        IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE

! Skip walls where some values are undefined.
        IF(WALL_AT(IJK)) cycle

         if(m/=0) then
            EPSA = AVG_Z(EP_S(IJK,M),EP_S(TOP_OF(IJK),M),K_OF(IJK))
            if(EPSA <= DIL_EP_s) CYCLE
         endif

         IF (.NOT.IP_AT_T(IJK)) THEN
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)

! evaluating the residual at cell ijk:
!   RESp = B-sum(Anb*VARnb)-Ap*VARp
!   (where nb = neighbor cells and p = center/0 cell)
            NUM1 = B_M(IJK,M) - (A_M(IJK,0,M)*W_M(IJK)+&
               A_M(IJK,east,M)*W_M(IPJK)+A_M(IJK,west,M)*W_M(IMJK)+&
               A_M(IJK,north,M)*W_M(IJPK)+A_M(IJK,south,M)*W_M(IJMK))
            IF (DO_K) THEN
               IJKM = KM_OF(IJK)
               IJKP = KP_OF(IJK)
               NUM1 = NUM1 - (A_M(IJK,top,M)*W_M(IJKP)+A_M(IJK,bottom,M)*W_M(IJKM))
            ENDIF

! Ignore momentum residual in stagnant regions.  Need an alternative
! criteria for residual scaling for such cases.
            VEL = SQRT(U_M(IJK)**2+V_M(IJK)**2+W_M(IJK)**2)
            IF (VEL > SMALL_NUMBER) THEN
               NUM1 = ABS(NUM1)
               DEN1 = ABS(A_M(IJK,0,M)*VEL)
! storing value of residual at each ijk location
               RESID_IJK(IJK) = NUM1
! adding to terms that are accumulated
               NCELLS = NCELLS + 1
               NUM = NUM + NUM1
               DEN = DEN + DEN1
            ENDIF
         ENDIF
      ENDDO
!$omp end parallel

      if(.not.debug_resid) return

! Collecting all the information among all the procesors -
! determining the global sum
      call global_all_sum(NUM)
      call global_all_sum(DEN)
      call global_all_sum(NCELLS)

      IJK_RESID = 1
      MAX_RESID = RESID_IJK( IJK_RESID )
      DO IJK = ijkstart3, ijkend3
      IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE

          IF (RESID_IJK( IJK ) > MAX_RESID) THEN
              IJK_RESID = IJK
              MAX_RESID = RESID_IJK( IJK_RESID )
          ENDIF
      ENDDO

! Determining the max residual
      do nproc=0,NumPEs-1
         if(nproc.eq.myPE) then
            MAX_RESID_L(nproc) = MAX_RESID
            IJK_RESID_L(nproc) = FUNIJK_GL(I_OF(IJK_RESID), J_OF(IJK_RESID), K_OF(IJK_RESID))
         else
            MAX_RESID_L(nproc) = 0.0
            IJK_RESID_L(nproc) = 0
         endif
      enddo

! Determining the maximum among all the procesors
      call global_all_max(MAX_RESID)

! Collecting all the information among all the procesors
      call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
      call global_all_sum(IJK_RESID_L, IJK_RESID_GL)

! Determining the global IJK location w.r.t. serial version
      IJK_RESID = IJKMAX2

      do nproc=0,NumPEs-1
        if(MAX_RESID_GL(nproc).eq.MAX_RESID.and.IJK_RESID_GL(nproc).lt.IJK_RESID) then
          IJK_RESID = IJK_RESID_GL(nproc)
        endif
      enddo

! Normalizing the residual
      IF (DEN > ZERO) THEN
         RESID = NUM/DEN
         MAX_RESID = NCELLS*MAX_RESID/DEN
      ELSE IF (NUM == ZERO) THEN
         RESID = ZERO
         MAX_RESID = ZERO
         IJK_RESID = 0
      ELSE
         RESID = LARGE_NUMBER
         MAX_RESID = LARGE_NUMBER
!         WRITE (LINE, *) 'Warning: All center coefficients are zero.'
!         CALL WRITE_ERROR ('CALC_RESID_W', LINE, 1)
      ENDIF

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

      END SUBROUTINE CALC_RESID_W
