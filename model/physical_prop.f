!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PHYSICAL_PROP                                           C
!  Purpose: Calculate the indicated physical properties that vary      C
!           with time if directed to do so by the corresponding flag   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-JUL-92  C
!  Reviewer: P. Nicoletti                             Date: 11-DEC-92  C
!                                                                      C
!  Revision: Mods for MFIX 2.0 (old name CALC_PHYSPROP)                C
!  Author: M. Syamlal                                 Date: 23-APR-96  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!    Perry, R.H., and C.H. Chilton, Chemical Engineer's Handbook, 5th  C
!      edition, McGraw-Hill Kogakusha, Tokyo, 1973.                    C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PHYSICAL_PROP(IER, LEVEL)

! Modules
!---------------------------------------------------------------------//
      use compar, only: numPEs, myPe, pe_io
      use funits, only: unit_log
      use exit, only: mfix_exit
      use mpi_utility, only: global_all_sum
      use physprop, only: smax
      use coeff, only: DENSITY
      use coeff, only: SP_HEAT
      use coeff, only: PSIZE
      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! Global error Flag.
      INTEGER, intent(inout) :: IER
      INTEGER, intent(in) :: LEVEL

! Local variables
!---------------------------------------------------------------------//
! Arrays for storing errors:
! 100 - Negative gas phase density
! 101 - Negative solids phase density
! 10x - Unclassified
! 900 - Invalid temperature in calc_CpoR
      INTEGER :: Err_l(0:numPEs-1)  ! local
      INTEGER :: Err_g(0:numPEs-1)  ! global
! loop index
      INTEGER :: M
!......................................................................!

! Initialize error flags.
      Err_l = 0

! Calculate density only. This is invoked several times within iterate,
! making it the most frequently called.
      if(LEVEL == 0) then
         if(DENSITY(0)) CALL PHYSICAL_PROP_ROg
           DO M=1,SMAX
              if(DENSITY(M)) CALL PHYSICAL_PROP_ROs(M)
           ENDDO

! Calculate everything except density. This is called at the start of
! each iteration.
      elseif(LEVEL == 1) then
         if(SP_HEAT(0)) CALL PHYSICAL_PROP_CPg
         DO M=1,SMAX
            if(SP_HEAT(M)) CALL PHYSICAL_PROP_CPs(M)
            if(PSIZE(M)) CALL PHYSICAL_PROP_Dp(M)
         ENDDO


! Calculate everything. This is invoked via calc_coeff_all as part of
! the initialization (before starting the time march) and at the start
! of each step step thereafter.
      elseif(LEVEL == 2) then
         if(DENSITY(0)) CALL PHYSICAL_PROP_ROg
         if(SP_HEAT(0)) CALL PHYSICAL_PROP_CPg
         DO M=1,SMAX
            if(DENSITY(M)) CALL PHYSICAL_PROP_ROs(M)
            if(SP_HEAT(M)) CALL PHYSICAL_PROP_CPs(M)
            if(PSIZE(M)) CALL PHYSICAL_PROP_Dp(M)
         ENDDO
      endif


! In case of negative density force exit from the physical property
! calculation routine and reduce the time step
      CALL global_all_sum(Err_l, Err_g)
      IER = maxval(Err_g)
      if(IER == 0) return


! Error handeling. - Local.
! An invalid temperature was found by calc_CpoR. This is a fatal run-
! time error and forces a call to MFIX_EXIT.
      IF(IER == 901 .OR. IER == 902) then
! Temperature is now bounded so this error will not occur.
         if(myPE == PE_IO) then
            write(*,2000) IER
            write(UNIT_LOG,2000) IER
         endif
         CALL MFIX_EXIT(myPE)
      ENDIF

      return

 2000 FORMAT(/1X,70('*')/' From: PHYSICAL_PROP',/' Fatal Error 2000:', &
         ' calc_CpoR reported an invalid temperature: 0x0', I3/,       &
         'See Cp.log for details. Calling MFIX_EXIT.',/1X,70('*')/)

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PHYSICAL_PROP_ROg                                       C
!  Purpose: Calculate the gas phase density.                           C
!                                                                      C
!  Author: J. Musser                                  Date: 28-JUN-13  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PHYSICAL_PROP_ROg

! Global variables:
!---------------------------------------------------------------------//
      use compar, only: ijkstart3, ijkend3
! Equation of State - GAS
      use eos, only: EOSG
! Gas phase species mass fractions.
      use fldvar, only: X_g
! Gas phase temperature.
      use fldvar, only: T_g
! Gas phase density (compressible).
      use fldvar, only: RO_g
! Gas phase pressure.
      use fldvar, only: P_g
! Gas phase volume fraction.
      use fldvar, only: EP_g
! Gas phase material density.
      use fldvar, only: ROP_g

      use functions, only: wall_at
! Run time flag for generating negative gas density log files
      use output, only: REPORT_NEG_DENSITY
      use param1, only: undefined, zero, one
      use physprop, only: database_read, mw_g, mw_mix_g, mw_avg, nmax
! Maximum value for molecular weight (divided by one)
      use toleranc, only: OMW_MAX
! Detect NaN in density for compressible flows.
      use utilities, only: mfix_isnan
! invoke user defined quantity
      USE usr_prop, only: usr_rog, calc_usr_prop
      USE usr_prop, only: gas_density
      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Average molecular weight
      DOUBLE PRECISION :: MW
! Loop indicies
      INTEGER :: IJK   ! Computational cell
! Flag to write log header
      LOGICAL :: wHeader
!......................................................................!

! Ensure that the database was read. This *should* have been caught by
! check_gas_phase but this call remains to prevent an accident.
      IF(.NOT.database_read) call read_database0(IER)

! User-defined function
      IF(USR_ROg) THEN
         CALL CALC_USR_PROP(Gas_Density,lm=0,lerr=err_l)
         RETURN
      ENDIF

! Initialize:
      wHeader = .TRUE.

      IJK_LP: DO IJK = IJKSTART3, IJKEND3
         IF(WALL_AT(IJK)) cycle IJK_LP
         IF (MW_AVG == UNDEFINED) THEN
! Calculating the average molecular weight of the fluid.
            MW = SUM(X_G(IJK,:NMAX(0))/MW_G(:NMAX(0)))
            MW = ONE/MAX(MW,OMW_MAX)
            MW_MIX_G(IJK) = MW
! Calculate the fluid density and bulk density
            RO_G(IJK) = EOSG(MW,P_G(IJK),T_G(IJK))
            ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
         ELSE
            RO_G(IJK) = EOSG(MW_AVG,P_G(IJK),T_G(IJK))
            ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
         ENDIF

         IF(RO_G(IJK) < ZERO .or. mfix_isnan(RO_G(IJK))) THEN
            Err_l(myPE) = 100
            IF(REPORT_NEG_DENSITY) CALL ROgErr_LOG(IJK, wHeader)
         ENDIF
      ENDDO IJK_LP

      RETURN
      END SUBROUTINE PHYSICAL_PROP_ROg


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PHYSICAL_PROP_ROs                                       C
!  Purpose: Calculate solids phase (variable) density.                 C
!                                                                      C
!  Author: J. Musser                                  Date: 28-JUN-13  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PHYSICAL_PROP_ROs(M)

! Global variables:
!---------------------------------------------------------------------//
      use compar, only: ijkstart3, ijkend3
! Equation of State - Solid
      use eos, only: EOSS
! Solid phase species mass fractions.
      use fldvar, only: X_s, ROP_s, RO_s
! Solid phase density (variable).
      use fldvar, only: ROP_s, RO_s

      use functions, only: wall_at
! Run time flag for generating negative density log files.
      use output, only: REPORT_NEG_DENSITY
      use param1, only: undefined, zero
! Baseline/Unreaced solids density
      use physprop, only: RO_s0
! Initial mass fraction of inert species
      use physprop, only: X_S0
! Index of inert solids phase species.
      use physprop, only: INERT_SPECIES
! Inert solids phase species mass fraction in dilute region.
      use physprop, only: DIL_INERT_X_VSD
! Factor to define dilute region where DIL_INERT_X_VSD is used
      use physprop, only: DIL_FACTOR_VSD
! Flag for variable solids density.
      use run, only: SOLVE_ROs
! Minimum solids volume fraction
      use toleranc, only: DIL_EP_s
! invoke user defined quantity
      USE usr_prop, only: usr_ros, calc_usr_prop
      USE usr_prop, only: solids_density
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local Variables:
!---------------------------------------------------------------------//
! Loop indicies
      INTEGER :: IJK
! Index of inert species
      INTEGER :: IIS
! Flag to write log header
      LOGICAL :: wHeader
! Minimum bulk density
      DOUBLE PRECISION :: minROPs
!......................................................................!

! User defined function
      IF(USR_ROs(m)) THEN
         CALL CALC_USR_PROP(Solids_Density,lm=M,lerr=err_l)
         RETURN
      ENDIF

      IF(SOLVE_ROs(M)) THEN
! Initialize header flag.
         wHeader = .TRUE.
! Set the index of the inert species
         IIS = INERT_SPECIES(M)
! Calculate the minimum solids denisty.
         minROPs = RO_s0(M)*(DIL_FACTOR_VSD*DIL_EP_s)

! Calculate the solids denisty over all cells.
         IJK_LP: DO IJK = IJKSTART3, IJKEND3
            IF(WALL_AT(IJK)) cycle IJK_LP
            IF(ROP_s(IJK,M) > minROPs) THEN
               RO_S(IJK,M) = EOSS(RO_s0(M), X_s0(M,IIS), &
                  X_s(IJK,M,IIS))
            ELSE
!               RO_s(IJK,M) = RO_s0(M)
               RO_S(IJK,M) = EOSS(RO_s0(M), X_s0(M,IIS), &
                  DIL_INERT_X_VSD(M))
            ENDIF

! Report errors.
            IF(RO_S(IJK,M) <= ZERO) THEN
               Err_l(myPE) = 101
               IF(REPORT_NEG_DENSITY) CALL ROsErr_LOG(IJK, M, wHeader)
            ENDIF
         ENDDO IJK_LP
      ENDIF

      RETURN
      END SUBROUTINE PHYSICAL_PROP_ROs


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PHYSICAL_PROP_CPg                                       C
!  Purpose: Calculate the gas phase constant pressure specific heat.   C
!                                                                      C
!  Author: J. Musser                                  Date: 28-JUN-13  C
!                                                                      C
! Notes:                                                               C
!  > Unit conversion: 1 cal = 4.183925 J                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PHYSICAL_PROP_CPg

! Global Variables:
!---------------------------------------------------------------------//
      use compar, only: ijkstart3, ijkend3
! Universal gas constant in cal/mol.K
      use constant, only: RGAS => GAS_CONST_cal
! Gas phase species mass fractions.
      use fldvar, only: X_g
! Gas phase temperature.
      use fldvar, only: T_g

      use functions, only: wall_at
      use output, only: REPORT_NEG_SPECIFICHEAT
      use param1, only: undefined, zero
      use physprop, only: database_read
      use physprop, only: MW_g, c_pg, nmax
! Function to calculate Cp over gas constant R
      use read_thermochemical, only: calc_CpoR
! Units: CGS/SI
      use run, only: UNITS
! invoke user defined quantity
      USE usr_prop, only: usr_cpg, calc_usr_prop
      USE usr_prop, only: gas_specificheat
      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Species specific heat.
      DOUBLE PRECISION :: lCp
! Loop indicies
      INTEGER :: IJK, NN
! Flag to write log header
      LOGICAL :: wHeader
! Error flag returned from calc_CpoR
      INTEGER :: lCP_Err
      INTEGER :: gCP_Err
!......................................................................!

! Ensure that the database was read. This *should* have been caught by
! check_gas_phase but this call remains to prevent an accident.
      IF(.NOT.database_read) CALL read_database0(IER)

! User defined function
      IF(USR_CPg) THEN
         CALL CALC_USR_PROP(Gas_SpecificHeat,lm=0,lerr=err_l)
         RETURN
      ENDIF

! Initialize
      wHeader = .TRUE.
      gCP_Err = 0
      lCP_Err = 0
      IJK_LP: DO IJK = IJKSTART3, IJKEND3
         IF(WALL_AT(IJK)) CYCLE IJK_LP
! Calculating an average specific heat for the fluid.
         C_PG(IJK) = ZERO
         DO NN = 1, NMAX(0)
            lCp = calc_CpoR(T_G(IJK), 0, NN)
! calc_cpor used to return error flag (101/102) on out of bound
! temperature. however, temperature is now bounded in routine.
            gCP_Err = max(gCP_Err, lCP_Err)
            C_PG(IJK) = C_PG(IJK) + X_g(IJK,NN) * lCp * RGAS / MW_g(NN)
         ENDDO

! report errors
         IF(C_PG(IJK) <= ZERO) THEN
!            gCP_err = 103
            IF(REPORT_NEG_SPECIFICHEAT) CALL CPgErr_LOG(IJK, wHeader)
         ENDIF

      ENDDO IJK_LP

! The database calculation always returns cal/g.K thus the following
! conversion is needed if using SI units.
! Only convert useful part of the array. When the array is re-indexed,
! elements past IJKEND3 are UNDEFINED and would keep growing if
! multiplied, leading to overflow.
      IF(UNITS == 'SI') THEN
         DO IJK = IJKSTART3, IJKEND3
           C_PG(IJK) = 4.183925d3 * C_PG(IJK)
         ENDDO
      ENDIF

! Increment the error to 900+ to invoke fatal exit.
      IF(gCP_Err /= 0) Err_l(myPE) = 800 + gCP_Err

      RETURN
      END SUBROUTINE PHYSICAL_PROP_CPg


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PHYSICAL_PROP_CPs                                       C
!  Purpose: Calculate solids phase constant pressure specific heat.    C
!                                                                      C
!  Author: J. Musser                                  Date: 28-JUN-13  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PHYSICAL_PROP_CPs(M)

! Global Variables:
!---------------------------------------------------------------------//
      use compar, only: ijkstart3, ijkend3
! Universal gas constant in cal/mol.K
      use constant, only: RGAS => GAS_CONST_cal
! Solids temperature
      use fldvar, only: T_s
! Solids species mass fractions
      use fldvar, only: X_s

      use functions, only: wall_at
      use output, only: REPORT_NEG_SPECIFICHEAT
      use param1, only: undefined, zero
      use physprop, only: database_read
      use physprop, only: MW_s, c_ps, nmax
! Function to calculate Cp over gas constant R
      use read_thermochemical, only: calc_CpoR
! Units: CGS/SI
      use run, only: UNITS
! invoke user defined quantity
      USE usr_prop, only: usr_cps, calc_usr_prop
      USE usr_prop, only: solids_specificheat
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables:
!---------------------------------------------------------------------//
! Local value for species specific heat (Cp)
      DOUBLE PRECISION :: lCp
! Loop indicies
      INTEGER :: IJK, NN
! Flag to write log header
      LOGICAL :: wHeader
! Error flag returned from calc_CpoR
      INTEGER :: lCP_Err
      INTEGER :: gCP_Err
!......................................................................!

! Ensure that the database was read. This *should* have been caught by
! check_solids_common_all but this call remains to prevent an accident.
      IF(.NOT.database_read) CALL read_database0(IER)

! User defined function
      IF(USR_CPs(M)) THEN
         CALL CALC_USR_PROP(Solids_SpecificHeat,lm=M,lerr=err_l)
         RETURN
      ENDIF

! Initialize
      wHeader = .TRUE.
      gCP_Err = 0
      lCP_Err = 0
      IJK_LP: DO IJK = IJKSTART3, IJKEND3
         IF(WALL_AT(IJK)) CYCLE IJK_LP
! Calculating an average specific heat for the fluid.
         C_PS(IJK, M) = ZERO
         DO NN = 1, NMAX(M)
            lCp = calc_CpoR(T_S(IJK,M), M, NN)
! calc_cpor used to return error flag (101/102) based on out of bound
! temperature. however, temperature is now bounded in routine
            gCP_Err = max(gCP_Err, lCP_Err)
            C_PS(IJK,M) = C_PS(IJK,M) + X_s(IJK,M,NN) * &
               (lCp * RGAS / MW_s(M,NN))
         ENDDO

! report errors
         IF(C_PS(IJK,M) <= ZERO) THEN
!            gCP_err = 104
            IF(REPORT_NEG_SPECIFICHEAT) CALL CPsErr_LOG(IJK, M, wHeader)
         ENDIF

      ENDDO IJK_LP

! The database calculation always returns cal/g.K thus the following
! conversion is needed if using SI units.
! Only convert useful part of the array. When the array is re-indexed,
! elements past IJKEND3 are UNDEFINED and would keep growing if
! multiplied, leading to overflow.
      IF(UNITS == 'SI') THEN
         DO IJK = IJKSTART3, IJKEND3
           C_PS(IJK,M) = 4.183925d3 * C_PS(IJK,M)
         ENDDO
      ENDIF

! Increment the error to 900+ to invoke fatal exit.
      IF(gCP_Err /= 0) Err_l(myPE) = 800 + gCP_Err

      RETURN
      END SUBROUTINE PHYSICAL_PROP_CPs


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PHYSICAL_PROP_DP                                        C
!  Purpose: Calculate solids phase diameter.                           C
!                                                                      C
!  Author: J. Musser                                  Date: 28-JUN-13  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PHYSICAL_PROP_Dp(M)

! Modules
!---------------------------------------------------------------------//
      use compar, only: ijkstart3, ijkend3
      use scalars, only: phase4scalar
      use fldvar, only: scalar
      use fldvar, only: D_p, EP_S
      use functions, only: wall_at
      use param1, only: small_number
      use run, only: CALL_DQMOM
      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! Loop indicies
      INTEGER :: IJK   ! Computational cell
! Map from true index to map.
      INTEGER :: lM
!......................................................................!

      IF(.NOT.CALL_DQMOM) return

      lM = phase4scalar(M) ! Map from scalar eq to solids phase

      IJK_LP: DO IJK = IJKSTART3, IJKEND3
         IF(WALL_AT(IJK)) CYCLE IJK_LP
         IF(EP_s(IJK,lM) > small_number) D_p(IJK,M)= Scalar(IJK,lM)
      ENDDO IJK_LP

      RETURN
      END SUBROUTINE PHYSICAL_PROP_Dp


      END SUBROUTINE PHYSICAL_PROP

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: ROgErr_LOG                                              C
!  Purpose: Record information about the location and conditions that  C
!           resulted in a negative gas phase density.                  C
!                                                                      C
!  Author: J. Musser                                  Date: 28-JUN-13  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE ROgErr_LOG(IJK, tHeader)

! Global Variables:
!---------------------------------------------------------------------//
      use compar, only: mype, numPEs
! Cutcell data.
      use cutcell, only: cartesian_grid
      use cutcell, only: cut_cell_at, small_cell_at
      use cutcell, only: xg_e, yg_n, zg_t
! Simulation time
      use run, only: TIME
! Gas phase temperature.
      use fldvar, only: T_g
! Gas phase density (compressible).
      use fldvar, only: RO_g
!
      use indices, only: i_of, j_of, k_of
! Gas phase pressure.
      use fldvar, only: P_g
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
      INTEGER, intent(in) :: IJK
      LOGICAL, intent(inout) :: tHeader

! Local variables:
!---------------------------------------------------------------------//
      LOGICAL :: lExists
      CHARACTER(LEN=255) :: lFile
      INTEGER, parameter :: lUnit = 4868
      LOGICAL, save :: fHeader = .TRUE.
!......................................................................!

      lFile = '';
      if(numPEs > 1) then
         write(lFile,"('ROgErr_',I4.4,'.log')") myPE
      else
         write(lFile,"('ROgErr.log')")
      endif
      inquire(file=trim(lFile),exist=lExists)
      if(lExists) then
         open(lUnit,file=trim(adjustl(lFile)),  &
            status='old', position='append', convert='big_endian')
      else
         open(lUnit,file=trim(adjustl(lFile)), status='new',&
              convert='big_endian')
      endif

      if(fHeader) then
         write(lUnit,1000)
         fHeader = .FALSE.
      endif

      if(tHeader) then
         write(lUnit,"(/2x,'Simulation time: ',g12.5)") TIME
         tHeader = .FALSE.
      endif

      write(lUnit,1001) IJK, I_OF(IJK), J_OF(IJK), K_OF(IJK)
      write(lUnit,"(6x,A,1X,g12.5)",ADVANCE='NO') 'RO_g:', RO_g(IJK)
      write(lUnit,"(2x,A,1X,g12.5)",ADVANCE='NO') 'P_g:', P_g(IJK)
      write(lUnit,"(2x,A,1X,g12.5)") 'T_g:', T_g(IJK)
      if(CARTESIAN_GRID) then
         write(lUnit,"(6x,A,1X,L1)",ADVANCE='NO') 'Cut Cell:', CUT_CELL_AT(IJK)
         write(lUnit,"(6x,A,1X,L1)") 'Small Cell:', SMALL_CELL_AT(IJK)
         write(lUnit,"(6x,'Coordinates (E/N/T): ',1X,3(2x, g17.8))") &
            xg_e(I_OF(IJK)), yg_n(J_of(ijk)), zg_t(k_of(ijk))
      endif

      close(lUnit)

      RETURN

 1000 FORMAT(2X,'One or more cells have reported a negative/NaN gas',  &
         ' density (RO_g(IJK)). If',/2x,'this is a persistent issue,', &
         ' lower UR_FAC(1) in mfix.dat.')

 1001 FORMAT(/4X,'IJK: ',I8,7X,'I: ',I4,'  J: ',I4,'  K: ',I4)

      END SUBROUTINE ROgErr_LOG


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ROsErr_LOG                                              !
!  Purpose: Record information about the location and conditions that  !
!           resulted in a negative solids phase density.               !
!                                                                      !
!  Author: J. Musser                                  Date: 09-Oct-13  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ROsErr_LOG(IJK, M, tHeader)

! Global variables:
!---------------------------------------------------------------------//
      use compar, only: mype, numPEs
! Cutcell data.
      use cutcell, only: cartesian_grid
      use cutcell, only: cut_cell_at, small_cell_at
      use cutcell, only: xg_e, yg_n, zg_t
! Solid phase species mass fractions.
      use fldvar, only: X_s
! Solid phase density (variable).
      use fldvar, only: RO_s
!
      use indices, only: i_of, j_of, k_of
! Baseline/Unreaced solids density
      use physprop, only: RO_s0
! Initial mass fraction of inert species
      use physprop, only: X_S0
! Index of inert solids phase species.
      use physprop, only: INERT_SPECIES
! Simulation time
      use run, only: TIME
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Fluid cell and solids phase indices
      INTEGER, intent(in) :: IJK, M
! Flag to output header
      LOGICAL, intent(inout) :: tHeader

! Local Variables:
!---------------------------------------------------------------------//
! Local aliase for inert species index
      INTEGER :: NN
! Local file values.
      LOGICAL :: lExists
      CHARACTER(LEN=255) :: lFile
      INTEGER, parameter :: lUnit = 4868
      LOGICAL, save :: fHeader = .TRUE.
!......................................................................!

      lFile = '';
      if(numPEs > 1) then
         write(lFile,"('ROsErr_',I4.4,'.log')") myPE
      else
         write(lFile,"('ROsErr.log')")
      endif
      inquire(file=trim(lFile),exist=lExists)
      if(lExists) then
         open(lUnit,file=trim(adjustl(lFile)), &
            status='old', position='append',convert='big_endian')
      else
         open(lUnit,file=trim(adjustl(lFile)), status='new',&
              convert='big_endian')
      endif

      if(fHeader) then
         write(lUnit,1000)
         fHeader = .FALSE.
      endif

      if(tHeader) then
         write(lUnit,"(/2x,'Simulation time: ',g12.5,'  Phase: ',I2)")&
             TIME, M
         tHeader = .FALSE.
      endif

      NN = INERT_SPECIES(M)
      write(lUnit,1001) IJK, I_OF(IJK), J_OF(IJK), K_OF(IJK)
      write(lUnit,"(6x,A,1X,g12.5)",advance='no') 'RO_s:', RO_s(IJK,M)
      write(lUnit,"(2x,A,1X,g12.5)",advance='no') 'Base:', RO_s0(M)
      write(lUnit,"(2x,A,1X,g12.5)",advance='no') 'X_s0:', X_s0(M,NN)
      write(lUnit,"(2x,A,1X,g12.5)") 'X_s:', X_s(IJK,M,NN)

      if(CARTESIAN_GRID) then
         write(lUnit,"(6x,A,1X,L1)",ADVANCE='NO') 'Cut Cell:', CUT_CELL_AT(IJK)
         write(lUnit,"(6x,A,1X,L1)") 'Small Cell:', SMALL_CELL_AT(IJK)
         write(lUnit,"(6x,'Coordinates (E/N/T): ',1X,3(2x, g17.8))") &
            xg_e(I_OF(IJK)), yg_n(J_of(ijk)), zg_t(k_of(ijk))
      endif

      close(lUnit)

      RETURN

 1000 FORMAT(2X,'One or more cells have reported a negative solids',   &
         ' density (RO_s(IJK)).')

 1001 FORMAT( 4X,'IJK: ',I8,7X,'I: ',I4,'  J: ',I4,'  K: ',I4)

      END SUBROUTINE ROsErr_LOG

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CPgErr_LOG                                              !
!  Purpose: Record information about the location and conditions that  !
!           resulted in a zero or negative gas specific heat.          !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CPgErr_LOG(IJK, tHeader)

! Global variables:
!---------------------------------------------------------------------//
      use compar, only: mype, numPEs
! Cutcell data.
      use cutcell, only: cartesian_grid
      use cutcell, only: cut_cell_at, small_cell_at
      use cutcell, only: xg_e, yg_n, zg_t
! Gas phase species mass fractions.
      use fldvar, only: T_g, X_g, EP_g
! Gas phase density (variable).
      use fldvar, only: RO_g
!
      use indices, only: i_of, j_of, k_of

      use physprop, only: MW_g, C_pg, nmax
! Simulation time
      use run, only: TIME
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Fluid cell and solids phase indices
      INTEGER, intent(in) :: IJK
! Flag to output header
      LOGICAL, intent(inout) :: tHeader

! Local Variables:
!---------------------------------------------------------------------//
! Local aliase for inert species index
      INTEGER :: NN
! Local file values.
      LOGICAL :: lExists
      CHARACTER(LEN=255) :: lFile
      INTEGER, parameter :: lUnit = 4868
      LOGICAL, save :: fHeader = .TRUE.
      CHARACTER(LEN=7) :: X_gN
!......................................................................!

      lFile = '';
      if(numPEs > 1) then
         write(lFile,"('CPgErr_',I4.4,'.log')") myPE
      else
         write(lFile,"('CPgErr.log')")
      endif
      inquire(file=trim(lFile),exist=lExists)
      if(lExists) then
         open(lUnit,file=trim(adjustl(lFile)), &
            status='old', position='append',convert='big_endian')
      else
         open(lUnit,file=trim(adjustl(lFile)), status='new',&
              convert='big_endian')
      endif

      if(fHeader) then
         write(lUnit,1000)
         fHeader = .FALSE.
      endif

      if(tHeader) then
         write(lUnit,"(/2x,'Simulation time: ',g12.5)") TIME
         tHeader = .FALSE.
      endif

      write(lUnit,1001) IJK, I_OF(IJK), J_OF(IJK), K_OF(IJK)
      write(lUnit,"(6x,A,1X,g12.5)",advance='no') 'C_PG:', C_PG(IJK)
      write(lUnit,"(2x,A,1X,g12.5)",advance='no') 'EP_G:', EP_g(IJK)
      write(lUnit,"(2x,A,1X,g12.5)") 'T_G:', T_g(IJK)
      write(lUnit,"(6x,A,1X)",advance='no') 'X_gN:'
      DO NN = 1,NMAX(0)
         write(lUnit,"(g12.5,2X)",advance='no') X_g(IJK,NN)
      ENDDO
      write(lUnit,fmt='(/)')

      if(CARTESIAN_GRID) then
         write(lUnit,"(6x,A,1X,L1)",ADVANCE='NO') 'Cut Cell:', CUT_CELL_AT(IJK)
         write(lUnit,"(6x,A,1X,L1)") 'Small Cell:', SMALL_CELL_AT(IJK)
         write(lUnit,"(6x,'Coordinates (E/N/T): ',1X,3(2x, g17.8))") &
            xg_e(I_OF(IJK)), yg_n(J_of(ijk)), zg_t(k_of(ijk))
      endif

      close(lUnit)

      RETURN

 1000 FORMAT(2X,'One or more cells have reported a negative or zero',&
         ' gas specific',/2X,'heat (C_pg(IJK)).')

 1001 FORMAT( 4X,'IJK: ',I8,7X,'I: ',I4,'  J: ',I4,'  K: ',I4)

      END SUBROUTINE CPgErr_LOG

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CPsErr_LOG                                              !
!  Purpose: Record information about the location and conditions that  !
!           resulted in a zero or negative solids specific heat.       !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CPsErr_LOG(IJK, M, tHeader)

! Global variables:
!---------------------------------------------------------------------//
      use compar, only: mype, numPEs
! Cutcell data.
      use cutcell, only: cartesian_grid
      use cutcell, only: cut_cell_at, small_cell_at
      use cutcell, only: xg_e, yg_n, zg_t
! Solid phase species mass fractions.
      use fldvar, only: T_s, X_s, EP_S
! Solid phase density (variable).
      use fldvar, only: RO_s
!
      use indices, only: i_of, j_of, k_of

      use physprop, only: MW_s, C_ps, nmax
! Simulation time
      use run, only: TIME
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Fluid cell and solids phase indices
      INTEGER, intent(in) :: IJK, M
! Flag to output header
      LOGICAL, intent(inout) :: tHeader

! Local Variables:
!---------------------------------------------------------------------//
! Local aliase for inert species index
      INTEGER :: NN
! Local file values.
      LOGICAL :: lExists
      CHARACTER(LEN=255) :: lFile
      INTEGER, parameter :: lUnit = 4868
      LOGICAL, save :: fHeader = .TRUE.
      CHARACTER(LEN=7) :: X_sN
!......................................................................!

      lFile = '';
      if(numPEs > 1) then
         write(lFile,"('CPsErr_',I4.4,'.log')") myPE
      else
         write(lFile,"('CPsErr.log')")
      endif
      inquire(file=trim(lFile),exist=lExists)
      if(lExists) then
         open(lUnit,file=trim(adjustl(lFile)), &
            status='old', position='append',convert='big_endian')
      else
         open(lUnit,file=trim(adjustl(lFile)), status='new',&
              convert='big_endian')
      endif

      if(fHeader) then
         write(lUnit,1000)
         fHeader = .FALSE.
      endif

      if(tHeader) then
         write(lUnit,"(/2x,'Simulation time: ',g12.5,'  Phase: ',I2)")&
             TIME, M
         tHeader = .FALSE.
      endif

      write(lUnit,1001) IJK, I_OF(IJK), J_OF(IJK), K_OF(IJK)
      write(lUnit,"(6x,A,1X,g12.5)",advance='no') 'C_PS:', C_PS(IJK,M)
      write(lUnit,"(2x,A,1X,g12.5)",advance='no') 'EP_S:', EP_s(IJK,M)
      write(lUnit,"(2x,A,1X,g12.5)") 'T_S:', T_s(IJK,M)
      write(lUnit,"(6x,A,1X)",advance='no') 'X_sN:'
      DO NN = 1,NMAX(M)
         write(lUnit,"(g12.5,2X)",advance='no') X_s(IJK,M,NN)
      ENDDO
      write(lUnit,fmt='(/)')

      if(CARTESIAN_GRID) then
         write(lUnit,"(6x,A,1X,L1)",ADVANCE='NO') 'Cut Cell:', CUT_CELL_AT(IJK)
         write(lUnit,"(6x,A,1X,L1)") 'Small Cell:', SMALL_CELL_AT(IJK)
         write(lUnit,"(6x,'Coordinates (E/N/T): ',1X,3(2x, g17.8))") &
            xg_e(I_OF(IJK)), yg_n(J_of(ijk)), zg_t(k_of(ijk))
      endif

      close(lUnit)

      RETURN

 1000 FORMAT(2X,'One or more cells have reported a negative or zero',&
         ' solids specific',/2X,'heat (C_ps(IJK)).')

 1001 FORMAT( 4X,'IJK: ',I8,7X,'I: ',I4,'  J: ',I4,'  K: ',I4)

      END SUBROUTINE CPsErr_LOG

