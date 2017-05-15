!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_RO_s                                                !
!                                                                      !
!  Author: J.Musser                                   Date: 09-Oct-13  !
!  Reviewer:                                                           !
!                                                                      !
!  Purpose: Initialize solids densities.                               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_RO_S

! Global Variables:
!---------------------------------------------------------------------/
! Number of solids phases.
      use physprop, only: MMAX
! Solids density field variable.
      use fldvar, only: RO_s, ROP_s
! Solid phase species mass fractions.
      use fldvar, only: X_s
! Initial mass fraction of inert species
      use physprop, only: X_S0
! Index of inert solids phase species.
      use physprop, only: INERT_SPECIES
! Inert solids phase species mass fraction in dilute region.
      use physprop, only: DIL_INERT_X_VSD
! Factor to define dilute region where DIL_INERT_X_VSD is used
      use physprop, only: DIL_FACTOR_VSD
! Run-time flag for variable soilds density
      use run, only: SOLVE_ROs
! Constant solids density.
      use physprop, only: RO_s0
! Minimum solids volume fraction
      use toleranc, only: DIL_EP_s

! Function for evaluating solids density.
      use eos, only: EOSS

! Modules needed to support function.inc
      use compar
      use geometry
      use indices
      use functions

      implicit none

! Local Variables:
!---------------------------------------------------------------------/
! Solids phase index
      INTEGER :: M
! Fluid cell index
      INTEGER :: IJK
! Index of the inert solids species.
      INTEGER :: IIS
! Flag for debugging.
      LOGICAL, parameter :: dbgMode = .FALSE.

      DOUBLE PRECISION :: minROPs

! Loop over all solids
      DO M=1,MMAX

! Variable solids density.
         IF (SOLVE_ROs(M)) THEN
! Set the index of the intert phase.
            IIS = INERT_SPECIES(M)
! Calculate the minimum solids denisty.
!            minROPs = RO_s0(M)*DIL_EP_s
            minROPs = RO_s0(M)*(DIL_FACTOR_VSD*DIL_EP_s)
! Debug/Development option.
            IF(dbgMode) CALL CHECK_SET_ROs()

! Calculate Ro_s in all fluid and flow boundary cells.
            DO IJK = ijkStart3, ijkEnd3
               IF(WALL_AT(IJK)) CYCLE
               IF(ROP_s(IJK,M) > minROPs) THEN
                  RO_S(IJK,M) = EOSS(RO_s0(M), X_s0(M,IIS),         &
                     X_s(IJK,M,IIS))
               ELSE
!                  RO_s(IJK,M) = RO_s0(M)
                  RO_S(IJK,M) = EOSS(RO_s0(M), X_s0(M,IIS),            &
                     DIL_INERT_X_VSD(M))
               ENDIF
            ENDDO
         ELSE
! Constant solids density.
            DO IJK = ijkstart3, ijkend3
               IF (WALL_AT(IJK)) CYCLE
               RO_S(IJK,M) = RO_S0(M)
            ENDDO
         ENDIF
      ENDDO

      RETURN

      CONTAINS

!``````````````````````````````````````````````````````````````````````!
!  Subroutine: CHECK_SET_ROs                                           !
!  Author: J.Musser                                   Date: 21-JAN-92  !
!                                                                      !
!  Purpose: Verify that all the variable solids density information is !
!           present for solids phase M.                                !
!                                                                      !
!           Note: The check_data routines should have caught any       !
!           problematic IC/BC specifications. This is included mainly  !
!           for development efforts.                                   !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE CHECK_SET_ROs()

      use exit, only: mfix_exit
! Flag for who writes
      use funits, only: DMP_LOG
! Solids species mass fractions.
      use fldvar, only: X_s

      use param1, only: zero
! Number of phase species.
      use physprop, only: NMAX
! Index of inert species.
      use physprop, only: INERT_SPECIES

      use toleranc

      implicit none

! Sum of solids phase mass fractions.
      DOUBLE PRECISION :: SUM_Xs
! Index of inert solids phase species.
      INTEGER :: INERT
! Integer Error Flag.
      INTEGER :: IER(2)

! Error file log.
      INTEGER, parameter :: lUnit = 8454
      LOGICAL :: lExists
      CHARACTER(LEN=64) :: lFName

! Initialize error flags.
      IER = 0

! Set the inert species index.
      INERT = INERT_SPECIES(M)

! Check all computational cells.
      DO IJK = ijkStart3, ijkEnd3
! Skip walls.
         IF (WALL_AT(IJK)) CYCLE
! Calculate the solids species mass fraction sum.
         SUM_Xs = sum(X_s(IJK,M,:NMAX(M)))
! Verify that the species mass fractions are specified and valid.
         IF(.NOT.compare(ONE,SUM_Xs)) IER(1) = IER(1)+1
! Verify that the inert species mass fraction is greater than zero.
         IF(X_s(IJK,M,INERT) <= ZERO) IER(2) = IER(2)+1

      ENDDO

! An error was detected. Open a log file.
      IF(sum(IER) /= 0) THEN
         lFName=''
         IF(numPEs == 1) THEN
            WRITE(lFName,"('setROs.log')")
         ELSE
            WRITE(lFName,"('setROs_',I6.6,'.log')") myPE
         ENDIF
         inquire(file=trim(lFName),exist=lExists)
         IF(lExists) THEN
            OPEN(CONVERT='BIG_ENDIAN',unit=lUnit,file=trim(lFName),status='replace')
         ELSE
            OPEN(CONVERT='BIG_ENDIAN',unit=lUnit,file=trim(lFName),status='new')
         ENDIF
      ENDIF


! An error was detected in species mass fraction sum.
      IF(IER(1) /= 0)THEN
         WRITE(lUnit,1100) myPE
! Skip walls.
         DO IJK = ijkStart3, ijkEnd3
            IF (WALL_AT(IJK)) CYCLE
! Calculate the solids species mass fraction sum.
            SUM_Xs = sum(X_s(IJK,M,:NMAX(M)))
! Verify that the species mass fractions are specified and valid.
            IF(.NOT.compare(ONE,SUM_Xs)) WRITE(lUnit,1101) IJK, SUM_Xs
         ENDDO
         WRITE(lUnit,9999)
      ENDIF

! An error was detected in inert species mass fraction.
      IF(IER(2) /= 0)THEN
         WRITE(lUnit,1200) myPE
         WRITE(lUnit,1201) M
         WRITE(lUnit,1202) INERT
! Skip walls.
         DO IJK = ijkStart3, ijkEnd3
            IF (WALL_AT(IJK)) CYCLE
! Calculate the solids species mass fraction sum.
! Verify that the species mass fractions are specified and valid.
            IF(X_s(IJK,M,INERT) <= ZERO) WRITE(lUnit,1203)             &
               IJK, X_s(IJK,M,INERT)
         ENDDO
         WRITE(lUnit,9999)
      ENDIF

! Close the file, cleanup, and exit.
      IF(sum(IER) /= 0) THEN
         CLOSE(lUnit)
         IF(DMP_LOG) THEN
         ENDIF
         CALL MFIX_EXIT(myPE)
      ENDIF


      RETURN

 1100 FORMAT(//1X,70('*')/' From: CHECK_SET_ROs',/,' Error 1100:',     &
         ' One or more fluid cells contain invalid species mass',/     &
         ' fractions which do NOT sum to one.'/,'   > myPE = ',I6)

 1101 FORMAT('   > sum(X_s(',I6,')) = ',g12.5)

 1200 FORMAT(//1X,70('*')/' From: CHECK_SET_ROs',/,' Error 1200:',     &
         ' One or more fluid cells contain an invalid species mass',/  &
         ' fraction for the inert material.'/,'   > myPE = ',I6)

 1201 FORMAT('   > Solid Phase: ',I2)

 1202 FORMAT('   > Inert species index: ',I4)

 1203 FORMAT('   > X_s(',I6,',INERT) = ',g12.5)

 9999 FORMAT(1x,70('*')/)

      END SUBROUTINE CHECK_SET_ROs

      END SUBROUTINE SET_RO_S
