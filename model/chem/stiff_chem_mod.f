      MODULE STIFF_CHEM

! External Routines.
!---------------------------------------------------------------------//
! Routine used to calculate the reaction rates and populate the
! fluid variable ODEs for ODEPACK.
      external STIFF_CHEM_RRATES
! Routine used to compare to values.
      LOGICAL, external :: COMPARE
! Routine used to calculate species enthalpies.
      DOUBLE PRECISION, external :: CALC_H0


! Runtime Flags:
!---------------------------------------------------------------------//
! Flag to invoke stiff chemistry solver.
      LOGICAL :: STIFF_CHEMISTRY
! Flag to invoke variable solids density.
!      LOGICAL :: VARIABLE_DENSITY
! Flag indicating if cell IJK is own by myPE.
      LOGICAL, dimension(:), allocatable :: notOwner

! ODEPACK Controlling parameters:
!---------------------------------------------------------------------//
! Dimension of ODEs solved in stiff solver.
      INTEGER :: ODE_DIMN_all
! Dimension of ODEs solved in stiff solver for gas phase only.
      INTEGER :: ODE_DIMN_g

! Dimension of ODEs solved in stiff solver for gas phase only.
      INTEGER :: NEQ_DIMN

! Indicates type of Error control.
      INTEGER :: ODE_ITOL
! Relative error tolerance paramter.
      DOUBLE PRECISION, DIMENSION(1) :: ODE_RTOL
! Absolute error tolerance parameter. (Dimension (ODE_DIMN))
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ODE_ATOL
! Declared length of RWORK.
      INTEGER :: ODE_LRW
! Declared length of IWORK.
      INTEGER :: ODE_LIW
! Jacobian type indicator.
      INTEGER :: ODE_JT
! The maximum number of steps ODEPACK may use to integrate.
      INTEGER :: STIFF_CHEM_MAX_STEPS
! Flag indicating that the max number of steps is unlimited.
      LOGICAL :: UNLIMITED_STEPS

! Explicit interface for ODEPACK
!---------------------------------------------------------------------//
      INTERFACE
         SUBROUTINE DLSODA (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, &
            ITASK,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, JT)
            external F
            INTEGER :: ITOL, ITASK, ISTATE, IOPT, LRW, LIW, JT
            INTEGER, dimension(2) :: NEQ
            INTEGER, dimension(LIW) :: IWORK
            DOUBLE PRECISION :: T, TOUT
            DOUBLE PRECISION :: JAC
            DOUBLE PRECISION, dimension(1) :: RTOL
            DOUBLE PRECISION, dimension(LRW) :: RWORK
            DOUBLE PRECISION, dimension(NEQ(1)) :: Y, ATOL
         END SUBROUTINE DLSODA
      END INTERFACE

      contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C
!     Module name: MCHEM_TIME_MARCH                                       C
!     Purpose: Called in time_march.f to do rxns calcs                    C
!                                                                         C
!     Author: Nan Xie                                   Date: 02-Aug-04   C
!     Reviewer:                                         Date:             C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE STIFF_CHEM_SOLVER(ODE_DT, iErr)

! External Module Procedures:
!---------------------------------------------------------------------//
      use stiff_chem_maps, only : mapODEtoMFIX
      use stiff_chem_maps, only : mapMFIXtoODE

! Global Variables:
!---------------------------------------------------------------------//
      use output,   only : FULL_LOG
      use param1,   only : zero
      use run,      only : TIME

      use mpi_utility

      use stiff_chem_dbg
      use stiff_chem_stats
      use functions

      implicit none

! Passed Variables:
!----------------------------------------------------------------------!
! Time integral length.
      DOUBLE PRECISION, intent(IN) :: ODE_DT
! Error Flag
      INTEGER, intent(OUT) :: iErr

! Local Variables:
!----------------------------------------------------------------------!
! Error flag -> Integration failed in one or more fluid cells.
      LOGICAL :: lErr_l   ! local

! Fluid Cell Index
      INTEGER :: IJK
! The maximum number of ODEs to solve.
      DOUBLE PRECISION, dimension(ODE_DIMN_all) :: ODE_VARS

! (1) :: Number of ODEs
! (2) :: Fluid cell index (IJK) passed into ODEPACK
! (:) :: Flag for solving solid phase M density and species equations.
      INTEGER, dimension(NEQ_DIMN) :: lNEQ
! Start time for integration
      DOUBLE PRECISION :: lT
! Stop time for integration
      DOUBLE PRECISION :: lTOUT
! Indicates type of Error control.
      INTEGER :: lITOL
! Relative error tolerance paramter.
      DOUBLE PRECISION :: lRTOL(1)
! Absolute error tolerance parameter. (Dimension (ODE_DIMN))
      DOUBLE PRECISION :: lATOL(ODE_DIMN_all)
! Index specifying the ODEPACK task.
      INTEGER :: lITASK
! Specifies the state of ODEPACK
      INTEGER :: lISTATE
! Flag indicating optional inputs are used.
      INTEGER :: lIOPT
! Array for REAL* work
      DOUBLE PRECISION :: RWORK(ODE_LRW)
! Declared length of RWORK.
      INTEGER :: lLRW
! Array for Integer work
      INTEGER :: IWORK(ODE_LIW)
! Declared length of IWORK.
      INTEGER :: lLIW
! Jacobain Routine (not used)
      DOUBLE PRECISION :: lJAC
! Jacobian type indicator.
      INTEGER :: lJT

! The number of attempts of a specific fluid cell.
      INTEGER :: lAtps

      LOGICAL :: lReset
      LOGICAL :: lIncpt

      lErr_l = .FALSE.

      CALL INIT_STIFF_CHEM_STATS

      IJK_LP: DO IJK = IJKSTART3, IJKEND3
         IF(notOwner(IJK)) cycle IJK_LP
         IF(FLUID_AT(IJK)) THEN

            lAtps = 0
            lReset = .FALSE.
            lIncpt = .FALSE.

! Forced restset of tolerance values.
            lRTOL = ODE_RTOL
            lATOL = ODE_ATOL

! Increment the attempts counter.
  50        lAtps = lAtps + 1

! Forced restset to initial values.
            lT    = 0.0d0
            lTOUT = ODE_DT
            lITOL = ODE_ITOL
            lLRW  = ODE_LRW
            lLIW  = ODE_LIW
            lJT   = ODE_JT

! Fixed parameters
            lITASK  = 1
            lISTATE = 1
            lIOPT   = 1

! Calculate the number of ODEs to solve.
            CALL CALC_ODE_COEFF(lNEQ, IJK)

! Clear the work arrays.
            IWORK = 0
            RWORK = ZERO

! The maximum number of internal steps ODEPACK may use to integrate over
! the time interval. The default value is 500.
            IWORK(6) = STIFF_CHEM_MAX_STEPS

            IF(CALC_REACTIONS(IJK)) THEN

! Map MFIX variables to ODE variables.
               CALL mapMFIXtoODE(NEQ_DIMN, lNEQ, ODE_DIMN_all, ODE_VARS)

! Store a copy of the original field variables. This allows for these
! values to be 'reset' in the event that the stiff solver fails.
               CALL UPDATE_ODE_OLD(ODE_DIMN_all, ODE_VARS)

! Clear the error flag.
 100           iErr = 0

! Integrate flow field variables to incorporate reactions.
               CALL DLSODA(STIFF_CHEM_RRATES, lNEQ, ODE_VARS, lT,      &
                  lTOUT, lITOL, lRTOL, lATOL, lITASK, lISTATE, lIOPT,  &
                  RWORK, lLRW, IWORK, lLIW, lJAC, lJT)

! Verify that the results are well defined.
               CALL CHECK_ODE_DATA(NEQ_DIMN, lNEQ, ODE_DIMN_all,       &
                  ODE_VARS, UNLIMITED_STEPS, lISTATE, iErr)

! Successfully Integrated ODEs.
               IF(iErr == 0) THEN
                  lReset = .FALSE.
! Additional integration steps are needed (lT < lTOUT).
               ELSEIF(iErr == -1) THEN
! Reste the state flag and keep integrating.
                  IF(UNLIMITED_STEPS) THEN
                     lISTATE = 2
                     goto 100
                  ELSE
                     lReset = .FALSE.
                     lIncpt = .TRUE.
                  ENDIF

! Too much accuracy was requested.
               ELSEIF(iErr == -2) THEN
                  IF(lAtps < 3) THEN
! Write to the error log file.
                     IF(ODE_DEBUG_LEVEL >= 2) CALL WRITE_ODE_LOG(iErr, &
                        NEQ_DIMN, lNEQ, ODE_DIMN_all, ODE_VARS)
! Reset the ODE variable array.
                     CALL RESET_ODE(ODE_DIMN_all, ODE_VARS, lAtps)
! Reset the field variables.
                     CALL mapODEtoMFIX(NEQ_DIMN, lNEQ,               &
                        ODE_DIMN_all, ODE_VARS)
! Loosen the convergence criteria and try again.
                     lRTOL = ODE_RTOL*10.0d0
                     lATOL = ODE_ATOL*10.0d0
                     goto 50
                  ELSE
! Write to the error log file.
                     IF(ODE_DEBUG_LEVEL >= 1) CALL WRITE_ODE_LOG(iErr, &
                        NEQ_DIMN, lNEQ, ODE_DIMN_all, ODE_VARS)
! Set the flag to reset the field variables to the initial values.
                     lReset = .TRUE.
                  ENDIF
! All other errors.
               ELSE
! Tighten the convergence criteria and try again.
                  IF(lAtps < 3) THEN
! Write to the error log file.
                     IF(ODE_DEBUG_LEVEL >= 2) CALL WRITE_ODE_LOG(iErr, &
                        NEQ_DIMN, lNEQ, ODE_DIMN_all, ODE_VARS)
! Reset the ODE variable array.
                     CALL RESET_ODE(ODE_DIMN_all, ODE_VARS, lAtps)
! Rest the filed variables to their original values.
                     CALL mapODEtoMFIX(NEQ_DIMN, lNEQ,                &
                        ODE_DIMN_all, ODE_VARS)
! Reduce the tolerances and try again.
                     lRTOL = ODE_RTOL/(10.0d0**lAtps)
                     lATOL = ODE_ATOL/(10.0d0**lAtps)
                     goto 50
                  ELSE
! Write to the error log file.
                     IF(ODE_DEBUG_LEVEL >= 1) CALL WRITE_ODE_LOG(iErr, &
                        NEQ_DIMN, lNEQ, ODE_DIMN_all, ODE_VARS)
! Set the flag to reset the field variables to the initial values.
                     lReset = .TRUE.
                  ENDIF

               ENDIF ! IF(iErr == 0)

! Reset the field variables.
               if(lReset) CALL RESET_ODE(ODE_DIMN_all, ODE_VARS, lAtps)
! Store the results in the field variables.
               CALL mapODEtoMFIX(NEQ_DIMN, lNEQ, ODE_DIMN_all, ODE_VARS)
! Collect solver stats.
               IF(FULL_LOG) CALL UPDATE_STIFF_CHEM_STATS(lNEQ, &
                  NEQ_DIMN, IWORK(11), ODE_DIMN_all, lAtps, lIncpt)


            ENDIF  ! EndIF CALC_REACTIONS
         ENDIF  ! IF(CALC_REACTIONS(IJK))
      END DO IJK_LP ! End Loop over fluiod Cells, IJK

!      gErr_l = .FALSE.
!      CALL GLOBAL_ALL_OR(lErr_l, gErr_l)
!      IF(gErr_l) CALL WRITE_VTU_FILE


      CALL FINALIZE_STIFF_SOLVER()

      IF(FULL_LOG) CALL WRITE_STIFF_CHEM_STATS()

      iErr = 0

      RETURN
      END SUBROUTINE STIFF_CHEM_SOLVER


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CALC_REACTIONS                                         !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 07-Feb-13  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      LOGICAL FUNCTION CALC_REACTIONS(IJK)

      use rxns,   only : NO_OF_RXNS

      implicit none

      INTEGER, intent(in) :: IJK

      DOUBLE PRECISION :: RATES(NO_OF_RXNS)

      DOUBLE PRECISION, parameter :: rLimit = 1.0d-8

! Initialize
      RATES = 0.0d0

! Calculate user defined reaction rates.
      CALL USR_RATES(IJK, RATES)

! If there is little to no reaction in the cell, then set the ODE
! Time to zero to avoid calling the stiff solver.
      CALC_REACTIONS = .TRUE.
      if(maxval(RATES) < rLimit) CALC_REACTIONS = .FALSE.


      RETURN
      END FUNCTION CALC_REACTIONS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CALC_DIMN_ODE                                          !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 07-Feb-13  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_ODE_COEFF(lNEQ, IJK)

      use fldvar, only : EP_S
      use physprop, only : MMAX
      use run, only : SPECIES_EQ

      implicit none

      INTEGER, intent(in)  :: IJK
      INTEGER, intent(out) :: lNEQ(NEQ_DIMN)

      INTEGER :: M

      LOGICAL :: USE_SOLIDS_ODEs

! Initialize.
      USE_SOLIDS_ODEs = .FALSE.
      lNEQ(2) = IJK
      lNEQ(3:) = 0

! If there is little to no solids in the cell, then set the ODE
! dimension to gas phase only.
      DO M=1, MMAX
         IF(SPECIES_EQ(M)) THEN
            IF(EP_s(IJK,M) > 1.0d-6) THEN
               USE_SOLIDS_ODEs = .TRUE.
               lNEQ(2+M) = 1
            ENDIF
         ENDIF
      ENDDO

      IF(USE_SOLIDS_ODEs)THEN
         lNEQ(1) = ODE_DIMN_all
      ELSE
         lNEQ(1) = ODE_DIMN_g
         lNEQ(3:) = 0
      ENDIF


      RETURN
      END SUBROUTINE CALC_ODE_COEFF


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: FINALIZE_STIFF_SOLVER                                  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 07-Feb-13  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE FINALIZE_STIFF_SOLVER

! Gas phase volume fraction
      use fldvar, only : EP_g
! Gas phase denisty
      use fldvar, only: RO_g, ROP_g
! Gas phase pressure
      use fldvar, only: P_g
! Gas phase temperature
      use fldvar, only: T_g
! Gas phase species mass fraction
      use fldvar, only: X_g
! Molecular weight of each gas phase species
      use physprop, only : MW_g
! Gas phase mixture molecular weight.
      use physprop, only : MW_MIX_g

! Solids phase bulk density
      use fldvar, only: ROP_S
! Solids phase temperature
      use fldvar, only: T_s
! Solids phase species mass fractions
      use fldvar, only: X_s

! Number of solids phases
      use physprop, only: MMAX
! Number of species in each phase.
      use physprop, only: NMAX

! Double precision:  ONE = 1.0d0
      use param1, only: ONE
! Universal gas constant
      use constant, only : GAS_CONST

      use compar
      use mpi_utility
      use sendrecv
      use functions
      use utilities

      implicit none

! Local loop indicies.
      INTEGER :: IJK  ! Fluid Cell index.
      INTEGER :: M    ! Solids phase index
      INTEGER :: NN    ! Species index

      CALL send_recv(EP_G,2)
      CALL send_recv(RO_G,2)
      CALL send_recv(ROP_G,2)
      CALL send_recv(T_G,2)

      DO NN=1,NMAX(0)
         CALL send_recv(X_G(:,NN),2)
         CALL BOUND_X (X_G(1,NN), IJKMAX2)
      ENDDO

      DO M = 1, MMAX
! Solids temperature.
         CALL send_recv(T_S(:,M),2)
! Solids volume fraction. (Constant Solids Density)
         CALL send_recv(ROP_S(:,M),2)
! Solids phase species mass fractions.
         DO NN=1,NMAX(M)
            CALL send_recv(X_S(:,M,NN),2)
            CALL BOUND_X (X_S(1,M,NN), IJKMAX2)
         ENDDO
      ENDDO

      DO IJK = ijkStart3, ijkEnd3
         IF(.NOT.FLUID_AT(IJK)) CYCLE
! Calculate the mixture molecular weight.
         MW_MIX_G(IJK) = sum(X_G(IJK,1:NMAX(0))/MW_g(1:NMAX(0)))
         MW_MIX_G(IJK) = ONE/MW_MIX_G(IJK)
! Calculate the gas phase pressure.
         P_G(IJK) = (RO_G(IJK)*GAS_CONST*T_G(IJK))/MW_MIX_G(IJK)
      ENDDO

      RETURN
      END SUBROUTINE FINALIZE_STIFF_SOLVER

      END MODULE STIFF_CHEM
