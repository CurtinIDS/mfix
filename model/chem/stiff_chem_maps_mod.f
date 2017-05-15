      MODULE STIFF_CHEM_MAPS

      PRIVATE

! Variable Access:
!---------------------------------------------------------------------//

! Subroutine Access:
!---------------------------------------------------------------------//
      PUBLIC :: mapMFIXtoODE,   &
                mapODEtoMFIX

      contains



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: mapMFIXtoODE                                           !
!                                                                      !
!  Purpose: This routine maps MFIX variables into the ODE array.       !
!                                                                      !
!  Author: J.Musser                                   Date: 07-Feb-13  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE mapMFIXtoODE(lnD, lNEQ, loD, lVars)

! Global Variables:
!---------------------------------------------------------------------//
      use fldvar, only: ROP_g
      use fldvar, only: ROP_s
      use fldvar, only: T_g
      use fldvar, only: T_s
      use fldvar, only: X_g
      use fldvar, only: X_s

      use physprop, only: NMAX
      use physprop, only: MMAX


      implicit none


! Passed Variables:
!---------------------------------------------------------------------//
! Passed array dimensions
      INTEGER, intent(in) :: lnD  ! lNEQ
      INTEGER, intent(in) :: loD  ! lVars

! (1) Number of ODEs to be solve
! (2) Fluid cell index
      INTEGER, intent(in) :: lNEQ(lnD)
! Array of dependent variable initial values.
      DOUBLE PRECISION, intent(out)  :: lVars(loD)


! Local Variables:
!---------------------------------------------------------------------//
! Fluid cell index
      INTEGER :: IJK
! Loop indicies:
      INTEGER :: M  ! phase
      INTEGER :: N  ! species
! ODE Equation Counter
      INTEGER :: Node


! Initialize.
      IJK = lNEQ(2)
      lVars = 0.0d0
      Node = 1

! Gas phase density.
      lVars(Node) = ROP_G(IJK);             Node = Node + 1
! Gas phase temperature.
      lVars(Node) = T_G(IJK);               Node = Node + 1
! Gas phase species mass fractions.
      DO N=1,NMAX(0)
         lVars(Node) = X_G(IJK,N);          Node = Node + 1
      ENDDO

! Solids temperature.
      DO M = 1, MMAX
         lVars(Node) = T_S(IJK,M);          Node = Node + 1
      ENDDO

      DO M = 1, MMAX
         IF(lNEQ(2+M) == 1) THEN
! Solids volume fraction.
            lVars(Node) = ROP_S(IJK,M);     Node = Node + 1
! Solids phase species mass fractions.
            DO N=1,NMAX(M)
               lVars(Node) = X_S(IJK,M,N);  Node = Node + 1
            ENDDO
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE mapMFIXtoODE



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: mapODEtoMFIX                                           !
!                                                                      !
!  Purpose: This is a driver routine for mapping variables stored in   !
!  the ODE array back to MFIX field variables.                         !
!                                                                      !
!  Author: J.Musser                                   Date: 07-Feb-13  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE mapODEtoMFIX(lnD, lNEQ, loD, lVars)

! Global Variables:
!---------------------------------------------------------------------//
! Material density.
      use fldvar, only : RO_g, RO_s
! Bulk density.
      use fldvar, only : ROP_g, ROP_s
! Gas temperature. (K)
      use fldvar, only : T_g, T_s
! Species Mass frcations.
      use fldvar, only : X_g, X_s
! Number of solids phases.
      use physprop, only : MMAX
! Number of species comprising each phase.
      use physprop, only : NMAX
! Gas constant (cal/g.K)
      use constant, only : GAS_CONST
! Gas pressure
      use fldvar, only : P_g
! Gas phase volume fraction.
      use fldvar,   only : EP_g
! Species molecular weights
      use physprop, only : MW_g
! Gas phase mixture molecular weight.
      use physprop, only : MW_MIX_g
! Baseline/Unreaced solids density.
      use physprop, only: RO_s0
! Initial mass fraction of inert solids species.
      use physprop, only: X_S0
! Index of inert solids phase species.
      use physprop, only: INERT_SPECIES
! Rank ID and Rank of IO process.
      use compar,   only : myPE
! Runtime flag for variable solids density.
      use run, only: SOLVE_ROs

! Variable solids density calculation.
      USE eos, ONLY: EOSS

! Global Parameters:
!---------------------------------------------------------------------//
      use param1,   only : ONE
      use param1,   only : LARGE_NUMBER
      use param1,   only : SMALL_NUMBER


      implicit none


! Passed Variables:
!---------------------------------------------------------------------//
! (1) Number of ODEs to be solve
! (2) Fluid cell index
      INTEGER, intent(in) :: lnD
      INTEGER, intent(in) :: lNEQ(lnD)

! Array of dependent variable initial values.
      INTEGER, intent(in) :: loD
      DOUBLE PRECISION, intent(in)  :: lVars(loD)


! Local Variables:
!---------------------------------------------------------------------//
! Fluid Cell index.
      INTEGER :: IJK
      INTEGER :: L

! Loop indicies:
      INTEGER :: M    ! phase
      INTEGER :: N    ! species
      INTEGER :: Node ! ODE Equation Counter

! Error flags/counters.
      INTEGER :: countNaN
      LOGICAL :: writeMsg

      IJK = lNEQ(2)

! Check if any NaNs were returned from ODEPACK.  This shouldn't happend
! but errors in the usr_rates file can cause this to occur.
!-----------------------------------------------------------------------
      countNaN = 0
      writeMsg = .FALSE.
      NaN_lp: do l=1, loD
         if(lVars(l).NE.lVars(l)) then
            countNaN = countNan + 1
            writeMsg = .TRUE.
         endif
      enddo NaN_lp

      if(writeMsg) then
         write(*,"(3x,'From MapODEtoMFIX: NaNs Found! :: ',3(3x,I4))") &
            myPE, IJK, countNaN

         if(countNaN < loD) then
            do l=1, loD
               if(lVars(l).NE.lVars(l))                                &
                  write(*,"(5x,' NaN in Var ',I2)") l
            enddo
         endif
      endif

! Directly map the ODE values into the field variable names.
!-----------------------------------------------------------------------
! Initialize the loop counter for ODEs.
      Node = 1

! Gas phase density.
      ROP_G(IJK) = lVars(Node);                    Node = Node + 1
! Gas phase temperature.
      T_G(IJK) = lVars(Node);                      Node = Node + 1
! Gas phase species mass fractions.
      DO N=1,NMAX(0)
         X_G(IJK,N) = lVars(Node);                 Node = Node + 1
      ENDDO

! Solids temperature.
      DO M = 1, MMAX
         IF(ROP_s(IJK,M) > 1.0d-8) &
            T_S(IJK,M) = lVars(Node);              Node = Node + 1
      ENDDO

! Only map back what was calculated.
      DO M = 1, MMAX
         IF(lNEQ(2+M) == 1) THEN
! Solids bulk density.
            ROP_S(IJK,M) = lVars(Node);            Node = Node + 1
! Solids phase species mass fractions.
            DO N=1,NMAX(M)
               X_S(IJK,M,N) = lVars(Node);         Node = Node + 1
            ENDDO
! Update the solids density from speices composition. This update must
! come after the species data is mapped back into the MFIX variables.
            IF(SOLVE_ROs(M)) RO_S(IJK,M) = EOSS(RO_s0(M),              &
               X_s0(M,INERT_SPECIES(M)), X_s(IJK,M,INERT_SPECIES(M)))
         ENDIF
      ENDDO

! Calculate the gas volume fraction from solids volume fractions. Only
! update it's value if the solids equations are being solved.
      IF(sum(lNEQ(3:)) > 0) EP_G(IJK) = &
         ONE - sum(ROP_S(IJK,1:MMAX)/RO_S(IJK,1:MMAX))

! Gas phase bulk density is updated within the stiff solver (lVar(1)).
! Now that the gas phase volume fraction is updated, the gas phase
! density can be backed out. RO_g * EP_g = ROP_g
      IF(EP_g(IJK) > small_number) THEN
         RO_g(IJK) = ROP_g(IJK) / EP_g(IJK)
      ELSE
! This case shouldn't happen, however 'LARGE_NUMBER' is used to aid
! in tracking errors should this somehow become and issue.
         RO_g(IJK) = LARGE_NUMBER
      ENDIF

! Calculate the mixture molecular weight.
      MW_MIX_G(IJK) = sum(X_G(IJK,1:NMAX(0))/MW_g(1:NMAX(0)))
      MW_MIX_G(IJK) = ONE/MW_MIX_G(IJK)

! Calculate the gas phase pressure.
      P_G(IJK) = (RO_G(IJK)*GAS_CONST*T_G(IJK))/MW_MIX_G(IJK)

      RETURN
      END SUBROUTINE mapODEtoMFIX

      END MODULE STIFF_CHEM_MAPS
