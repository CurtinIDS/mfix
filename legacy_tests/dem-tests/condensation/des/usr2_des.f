!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: URS2_DES                                               !
!                                                                      !
!  Purpose: This routine is called within the discrete phase time loop !
!  after the source terms are applied and the time step updated. The   !
!  The user may insert code in this routine or call user defined       !
!  subroutines.                                                        !
!                                                                      !
!  This routien is called from the time loop, but no indicies (fluid   !
!  cell or particle) are defined.                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 06-Nov-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR2_DES

      Use des_rxns
      Use des_thermo
      Use discretelement
      Use run
      Use usr

      IMPLICIT NONE


! Local variables
!---------------------------------------------------------------------//
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.
      DOUBLE PRECISION, SAVE :: OUT_TIME
      DOUBLE PRECISION, PARAMETER :: OUT_dT = 0.1d0


      IF(FIRST_PASS) THEN
         CALL WRITE_DES_Tp(S_TIME)
         OUT_TIME = OUT_dT
         FIRST_PASS = .FALSE.
      ELSE
         IF(S_TIME >= OUT_TIME) THEN
            CALL WRITE_DES_Tp(S_TIME)
            OUT_TIME = OUT_TIME + OUT_dT
         ENDIF
      ENDIF

      RETURN

      contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_URS2                                               !
!                                                                      !
!  Purpose: This routine is called within the discrete phase time loop !
!  after the source terms are applied and the time step updated. The   !
!  The user may insert code in this routine or call user defined       !
!  subroutines.                                                        !
!                                                                      !
!  This routien is called from the time loop, but no indicies (fluid   !
!  cell or particle) are defined.                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 06-Nov-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE WRITE_DES_Tp(lTime)

      Use des_rxns
      Use des_thermo
      Use discretelement
      Use fldvar
      Use run
      Use usr

      IMPLICIT NONE

! Passed variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION lTime

! Local variables
!---------------------------------------------------------------------//
! file name
      CHARACTER*64 :: FNAME
! logical used for testing is the data file already exists
      LOGICAL :: F_EXISTS
! file unit for heat transfer data
      INTEGER, PARAMETER :: TP_UNIT = 2030

      FNAME = 'POST_Tp.dat'
      INQUIRE(FILE=FNAME,EXIST=F_EXISTS)
      IF (.NOT.F_EXISTS) THEN
         OPEN(UNIT=TP_UNIT,FILE=FNAME,STATUS='NEW')
         WRITE(TP_UNIT,"(7X,'Time',12X,'Tp')")
      ELSE
         OPEN(UNIT=TP_UNIT,FILE=FNAME,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF

! Write the data to a file.
      WRITE(TP_UNIT,"(2(3X,F12.8))")lTime, (DES_T_s(1)-273.15d0)
      CLOSE(TP_UNIT)

      RETURN
      END SUBROUTINE WRITE_DES_Tp

      END SUBROUTINE USR2_DES
