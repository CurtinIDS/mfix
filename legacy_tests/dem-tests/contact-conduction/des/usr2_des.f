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

! Passed variables
!---------------------------------------------------------------------//
! None

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
      Use run
      Use usr

      IMPLICIT NONE

! Passed variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION lTime


! Local variables
!---------------------------------------------------------------------//
! file name
      CHARACTER*64 :: FNAME1, FNAME2
! logical used for testing is the data file already exists
      LOGICAL :: F_EXISTS1, F_EXISTS2
! file unit for heat transfer data
      INTEGER, PARAMETER :: TP_UNIT1 = 2030
      INTEGER, PARAMETER :: TP_UNIT2 = 2031

      DOUBLE PRECISION Tp10, Tp20, Tp1, Tp2, REL_ERR1, REL_ERR2
      DOUBLE PRECISION A_C, B_C

! Open the files.
      FNAME1 = 'POST_TP1.dat'
      INQUIRE(FILE=FNAME1,EXIST=F_EXISTS1)
      IF (.NOT.F_EXISTS1) THEN
         OPEN(UNIT=TP_UNIT1,FILE=FNAME1,STATUS='NEW')
         WRITE(TP_UNIT1,"(7X,'Time',11X,'Tp1',10X,'Tp1_MFIX',6X,'REL ERR 1')")
      ELSE
         OPEN(UNIT=TP_UNIT1,FILE=FNAME1,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF

      FNAME2 = 'POST_TP2.dat'
      INQUIRE(FILE=FNAME2,EXIST=F_EXISTS2)
      IF (.NOT.F_EXISTS2) THEN
         OPEN(UNIT=TP_UNIT2,FILE=FNAME2,STATUS='NEW')
         WRITE(TP_UNIT2,"(7X,'Time',11X,'Tp2',10X,'Tp2_MFIX',6X,'REL ERR 2')")
      ELSE
         OPEN(UNIT=TP_UNIT2,FILE=FNAME2,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF

! Calculate the value for the analytic solutions.
      Tp10 = 298.15d0
      Tp20 = 453.15d0

      A_C = 0.01027998273952850
      B_C = 0.00410668176759131

      Tp1 = 1.0d0/(A_C + B_C)*(B_C*Tp10 + A_C*Tp20 + &
         A_C*(Tp10-Tp20)*exp(-(A_C + B_C)*lTime))

      Tp2 = 1.0d0/(A_C + B_C)*(B_C*Tp10 + A_C*Tp20 - &
         B_C*(Tp10-Tp20)*exp(-(A_C + B_C)*lTime))

! Calculate the relative error.
       REL_ERR1 = (ABS(Tp1 - DES_T_s(1))/ABS(Tp1))*100
       REL_ERR2 = (ABS(Tp2 - DES_T_s(2))/ABS(Tp2))*100

! Write the data to a file.
      WRITE(TP_UNIT1,"(4(3X,F12.8))")lTime,Tp1,DES_T_s(1),REL_ERR1
      CLOSE(TP_UNIT1)

      WRITE(TP_UNIT2,"(4(3X,F12.8))")lTime,Tp2,DES_T_s(2),REL_ERR2
      CLOSE(TP_UNIT2)


      END SUBROUTINE WRITE_DES_Tp

      END SUBROUTINE USR2_DES
