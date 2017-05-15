!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_DATA_CHEM                                        !
!  Author: J.Musser                                   Date: 02-Aug-13  !
!                                                                      !
!  Purpose: Check the chemical rxns namelist variables.                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_ODEPACK_STIFF_CHEM

! Global Variables:
!---------------------------------------------------------------------//
! Dimension of IJK arrays.
      use param, only: DIMENSION_3
! Double precision zero.
      use param1, only: ZERO
! Double precision value for undefined variables.
      use param1, only: UNDEFINED
! Constant gas phase density
      use physprop, only : RO_G0
! Runtime logical for solving energy equations.
      use run, only: ENERGY_EQ
! Run time logical for solving species equations.
      use run, only: SPECIES_EQ
! Net rate of gas phase production/consumption
      use rxns, only: SUM_R_g
! Net rate of solids phase production/consumption
      use rxns, only: SUM_R_s
! Run time logical for using stiff chemistry solver
      use stiff_chem, only: STIFF_CHEMISTRY
! Run time logicals for identifying cells owned by myPE
      use stiff_chem, only: notOwner

! Full access to the following modules:
!---------------------------------------------------------------------//
      use compar
      use geometry
      use indices

      use error_manager
      use functions

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
      INTEGER :: I, J, K, IJK

      CALL INIT_ERR_MSG('CHECK_ODEPACK_STIFF_CHEM')

! Verify that there is sufficient run complexity to use the stiff solver
      IF(STIFF_CHEMISTRY) THEN

! Energy equations must be solved.
         IF(.NOT.ENERGY_EQ) THEN
            WRITE(ERR_MSG,1004)'ENERGY_EQ = .FALSE.'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

! Must be compressible
         IF(RO_G0 /= UNDEFINED) THEN
            WRITE(ERR_MSG,1004)'RO_G0 /= UNDEFINED'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF(.NOT.SPECIES_EQ(0)) THEN
            WRITE(ERR_MSG,1004)'SPECIES_EQ(0) = .FALSE.'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

! The stiff chemistry solver only needs to loop over physical cells
! owned by a process (e.g., not ghost cells). To avoid having a
! triple do loop, this array is populated to identify the cells that
! are not owned.
         ALLOCATE( notOwner(DIMENSION_3) ); notOwner = .TRUE.
         do k=kstart, kend
         do j=jstart, jend
         do i=istart, iend
            ijk = funijk(i,j,k)
            notOwner(IJK) = .FALSE.
         enddo
         enddo
         enddo

! Initialize ODEPACK operating parameters.
         CALL ODEPACK_INIT

! Clear the interphase mass transfer terms as the stiff solver
! does no use them.
         IF(allocated(SUM_R_g)) SUM_R_g = ZERO
         IF(allocated(SUM_R_s)) SUM_R_s = ZERO

      ENDIF

      CALL FINL_ERR_MSG



 1004 FORMAT('Error 1004: ',                                           &
         'Invalid parameters for stiff chemistry solver!',//           &
         ' The following criteria must be satisfied:',/                &
         '   > Solving the energy equations.',/                        &
         '   > Compressible gas phase.',/                              &
         '   > Solving gas phase species equations.',//                &
         ' >>> Invalid Parameter: ',A,//                               &
         'Check the user documentation and correct the mfix.dat file.')

      RETURN
      END SUBROUTINE CHECK_ODEPACK_STIFF_CHEM



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C
!     Module name: MCHEM_ODEPACK_INIT                                     C
!     Purpose: controlling values for ODEAPCK(reference to ODEPACK manual)C
!                                                                         C
!     Author: Nan Xie                                   Date: 02-Aug-04   C
!     Reviewer:                                         Date:             C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE ODEPACK_INIT

! Global Variables:
!---------------------------------------------------------------------//
      use physprop, only: MMAX
      use physprop, only: NMAX

      use run, only: SPECIES_EQ

      use stiff_chem, only: NEQ_DIMN
      use stiff_chem, only: ODE_DIMN_g
      use stiff_chem, only: ODE_DIMN_all

      use stiff_chem, only: ODE_ITOL
      use stiff_chem, only: ODE_RTOL
      use stiff_chem, only: ODE_ATOL

      use stiff_chem, only: ODE_LRW
      use stiff_chem, only: ODE_LIW
      use stiff_chem, only: ODE_JT

      use stiff_chem, only: STIFF_CHEM_MAX_STEPS

      use stiff_chem_dbg, only: ALLOCATE_STIFF_CHEM_DBG
      use stiff_chem_stats, only: ALLOCATE_STIFF_CHEM_STATS

! Double precision value for undefined variables.
      use param1, only: UNDEFINED_I

      implicit none


! Local Variables:
!---------------------------------------------------------------------//
      INTEGER :: LRN, LRS
      INTEGER :: M

      LOGICAL :: LIMIT_MAX_STEPS

! Number of ODEs (maximum)
      NEQ_DIMN = 2 + MMAX

! Dimension of ODEs solved if only the gas phase is present:
! Gas density, temperature, and species
      ODE_DIMN_g = 2 + NMAX(0)
! Solids temperature excluding for all phases.
      ODE_DIMN_g = ODE_DIMN_g + MMAX

! Calculate the total number of ODEs that are solve.
      ODE_DIMN_all = ODE_DIMN_g
      DO M=1, MMAX
! Solids bulk density and species.
         IF(SPECIES_EQ(M)) ODE_DIMN_all = ODE_DIMN_all + (1 + NMAX(M))
      ENDDO

! Indicates type of Error control.
      ODE_ITOL = 2 ! :: EWT(i) = RTOL * ABS(Y(i)) * ATOL(i)

! Relative error tolerance parameter.
      ODE_RTOL(1) = 1.0D-5

! Absolute error tolerance parameter.
      IF(.NOT.(allocated(ODE_ATOL))) allocate(ODE_ATOL(ODE_DIMN_all))
      ODE_ATOL(:) = 1.0D-6  ! All Equations

! Declared length of RWORK.
      LRN = 20 + 16*ODE_DIMN_all
      LRS = 22 + 9*ODE_DIMN_all + (ODE_DIMN_all**2)
      ODE_LRW = max(LRN, LRS)

! Declared length of IWORK.
      ODE_LIW = 20 + ODE_DIMN_all

! Jacobian type indicator.
      ODE_JT = 2 ! Internally generated.


      IF(STIFF_CHEM_MAX_STEPS == UNDEFINED_I) THEN
         STIFF_CHEM_MAX_STEPS = 500000
         LIMIT_MAX_STEPS = .FALSE.
      ELSE
         LIMIT_MAX_STEPS = .TRUE.
      ENDIF


      CALL ALLOCATE_STIFF_CHEM_DBG(ODE_DIMN_all)
      CALL ALLOCATE_STIFF_CHEM_STATS

      return
      END SUBROUTINE ODEPACK_INIT
