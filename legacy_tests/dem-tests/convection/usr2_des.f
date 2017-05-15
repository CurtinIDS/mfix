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


!-----------------------------------------------
! Local variables
!-----------------------------------------------
! file name
      CHARACTER*64 :: FNAME
! logical used for testing is the data file already exists
      LOGICAL :: F_EXISTS
! file unit for heat transfer data
      INTEGER, PARAMETER :: TP_UNIT = 2030

      DOUBLE PRECISION Tg, Tp0, Coeff, Tp,REL_ERR


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

      DOUBLE PRECISION Tg, Tp0, Coeff, Tp,REL_ERR


      FNAME = 'POST_Tp.dat'
      INQUIRE(FILE=FNAME,EXIST=F_EXISTS)
      IF (.NOT.F_EXISTS) THEN
         OPEN(UNIT=TP_UNIT,FILE=FNAME,STATUS='NEW')
         WRITE(TP_UNIT,"(7X,'Time',12X,'Tp',10X,'Tp_MFIX',8X,'REL ERR')")
      ELSE
         OPEN(UNIT=TP_UNIT,FILE=FNAME,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF

! Calculate the value for the analytic solution.
      Tg = 298.15d0
      Tp0 = 373.14994357d0
      Coeff =  -0.878629815377096

      Tp = Tg - (Tg - Tp0)*exp(coeff*lTime)

! Calculate the relative error.
      REL_ERR = (ABS(Tp - DES_T_S(1))/ABS(Tp))*100

! Write the data to a file.
      WRITE(TP_UNIT,"(4(3X,F12.8))")lTime,Tp,DES_T_S(1),REL_ERR
      CLOSE(TP_UNIT)

      END SUBROUTINE WRITE_DES_Tp

      END SUBROUTINE USR2_DES
